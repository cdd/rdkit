//
//  Copyright (C) 20126 Tad Hurst and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "FileParsers.h"
#include "FileParserUtils.h"
#include "MolSGroupParsing.h"
#include <RDGeneral/StreamOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/SubstanceGroup.h>

#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/FileParsers/SCSRUtils.h>
#include <GraphMol/FileParsers/FileWriters.h>
#include <GraphMol/MacroMol.h>
#include <GraphMol/FileParsers/MolFromMacroMol.h>

namespace RDKit {

enum AttachPointType : int {
  None=0,
  Al,
  Br,
  Cx,
  Dx,
  Ex,
  Fx,
  Gx,
  Hx,
  Ix,
  Jx,
  Hb=100,
  Hb1,
  Hb2,
  Hb3
};

namespace {
std::map<std::string, std::string> ScsrClassToMonomerMolClass {
  {"AA", "AminoAcid"},
  {"dAA", "AminoAcid"},
  {"MODAA", "AminoAcid"},
  {"MODdAA", "AminoAcid"},
  {"DNA", "NucleicAcid"},
  {"RNA", "NucleicAcid"},
  {"SUGAR", "NucleicAcid"},
  {"BASE", "NucleicAcid"},
  {"MODDNA", "NucleicAcid"},
  {"MODRNA", "NucleicAcid"},
  {"PHOSPHATE", "NucleicAcid"},
  {"XLINKAA", "Chemical"},
  {"XLINKdAA", "Chemical"},
  {"XLINKDNA", "Chemical"},
  {"XLINKRNA", "Chemical"},
  {"LINKER", "Chemical"},
  {"CHEM", "Chemical"},
  {"LGRP", "Chemical"},
  {"XLINKRNA", "Chemical"},
  {"XLINKRNA", "Chemical"},
};

std::map<std::string, AttachPointType> ScsrTypeToAttachPointType = {
        {"", AttachPointType::None},
        {"Al", AttachPointType::Al},
        {"Br", AttachPointType::Br},
        {"Cx", AttachPointType::Cx},
        {"Dx", AttachPointType::Dx},
        {"Ex", AttachPointType::Ex},
        {"Fx", AttachPointType::Fx},
        {"Gx", AttachPointType::Gx},
        {"Hx", AttachPointType::Hx},
        {"Ix", AttachPointType::Ix},
        {"Hb", AttachPointType::Hb},
        {"Hb1", AttachPointType::Hb1},
        {"Hb2", AttachPointType::Hb2},
        {"Hb3", AttachPointType::Hb3},
  };

std::map<const AttachPointType,std::string> AttachPointTypeToSCSRType = {
      {AttachPointType::Al, "Al"},
      {AttachPointType::Br, "Br"},
      {AttachPointType::Cx, "Cx"},
      {AttachPointType::Dx, "Dx"},
      {AttachPointType::Ex, "Ex"},
      {AttachPointType::Fx, "Fx"},
      {AttachPointType::Gx, "Gx"},
      {AttachPointType::Hx, "Hx"},
      {AttachPointType::Ix, "Ix"},
      {AttachPointType::Jx, "Jx"},
      {AttachPointType::Hb, "Hb"},
      {AttachPointType::Hb1, "Hb1"},
      {AttachPointType::Hb2, "Hb2"},
      {AttachPointType::Hb3, "Hb3"}
};

std::string getScsrClass(const MacroMolTemplate &macroMolTemplate) {
  switch( macroMolTemplate.getMonomerClass()) {
    case MonomerClass::NucleicAcid:
      {
        bool AlFound=false;
        bool BrFound=false;
        bool CxFound=false;

        for (auto leavingGroup : macroMolTemplate.getLeavingGroups()) {
          switch (leavingGroup.attachPoint) {
            case 1:
              AlFound = true;
              break;
            case 2:
              BrFound = true;
              break;
            case 3:
              CxFound = true;
          }
        }

        if (AlFound && BrFound && CxFound) {
          return "SUGAR";
        }
        if (AlFound && BrFound && !CxFound) {
          return "PHOSPHATE";
        }
        if (AlFound && !BrFound && !CxFound) {
          return "BASE";
        }
      }
      return "";
      
    case  MonomerClass::AminoAcid:
      return "AA";
    

    case MonomerClass::Chemical:
      return "CHEM";

    case MonomerClass::Other:
      throw FileParseException("Unknown macro monomer type");
  }
  throw FileParseException("Should Not happen");
}

void makeHydrogenBonds(const std::vector<bool> &donorFlags1,
                    const std::vector<bool> &donorFlags2,
                    MacroBondInfo &macroBondInfo) {
  if (donorFlags1.size() == 3 && donorFlags2.size() == 3) {
    // see if the two donor flag lists are
    // complementary

    bool complimentary = true;
    for (auto index = 0; index < 3; ++index) {
      if (donorFlags1[index] == donorFlags2[index]) {
        complimentary = false;  // both are donors or both are
                                // acceptors
        break;
      }
    }

  
    if (complimentary) {
      macroBondInfo.addBond(AttachPointType::Hb1,AttachPointType::Hb1,Bond::BondType::HYDROGEN );
      macroBondInfo.addBond(AttachPointType::Hb2,AttachPointType::Hb2,Bond::BondType::HYDROGEN );
      macroBondInfo.addBond(AttachPointType::Hb3,AttachPointType::Hb3,Bond::BondType::HYDROGEN );
      return;
    }
  }

  // not complimentary, but could be a wobble pair}

  std::string parings="";
  for (auto &donorFlags : {donorFlags1, donorFlags2}){
    for (auto donorFlag : donorFlags) {
      if (donorFlag) {
        parings += "D";
      } else {
        parings += "A";
      }
    }
    parings += ":";
  }

  // Check for wobble pairs.  Theses are:
  // G-U or U-G
  // I-A or A-I
  // I-C or C-I
  // I-U of U-I

  if (parings == "DDA:ADA:"  /* like G-U*/ ||
      parings == "DAD:DA:"    /* like A-I*/ ||
      parings == "AAD:DA:"    /* like C-I*/) {
    // use connection 2 and 3 from DDA,  and 1 and 2 from ADA
    macroBondInfo.addBond(AttachPointType::Hb2,AttachPointType::Hb1,Bond::BondType::HYDROGEN );
    macroBondInfo.addBond(AttachPointType::Hb3,AttachPointType::Hb2,Bond::BondType::HYDROGEN );
    return;
  } 

  if (parings == "ADA:DDA:"  /* like U-G*/ ||
      parings == "DA:DAD:"  /* like I-A*/ ||
      parings == "DA:AAD:"  /* like I-C*/ ) {
    // use connections 1 and 2 from the first, and 2 and 3 from the second
    macroBondInfo.addBond(AttachPointType::Hb1,AttachPointType::Hb2,Bond::BondType::HYDROGEN );
    macroBondInfo.addBond(AttachPointType::Hb2,AttachPointType::Hb3,Bond::BondType::HYDROGEN );
    return;
  } 

  if (parings == "DA:ADA:"  /* like I-G*/ ||
      parings == "ADA:DA:"  /* like G-I*/) {
    macroBondInfo.addBond(AttachPointType::Hb1,AttachPointType::Hb1,Bond::BondType::HYDROGEN );
    macroBondInfo.addBond(AttachPointType::Hb2,AttachPointType::Hb2,Bond::BondType::HYDROGEN );
      return;
  }

  // if here, we have no wobble bond, so just add
  // the bond between the first atoms on both sides
  // this is not a real pairing, but keeps the two bases associated with each other

  macroBondInfo.addBond(AttachPointType::Hb1,AttachPointType::Hb1,Bond::BondType::HYDROGEN );
  return;
}

void skipSpaces(const char *&linePtr) {
  while (*linePtr == ' ') {
    ++linePtr;  // skip spaces
  }
}

std::string getToken(const char *&linePtr, char delim) {
  skipSpaces(linePtr);
  unsigned int charCount = 0;
  while (linePtr[charCount] && linePtr[charCount] != delim) {
    ++charCount;
  }
  std::string res(linePtr, charCount);
  linePtr += charCount;
  if (*linePtr == delim) {
    ++linePtr;  // skip delim
  }
  return res;
}

std::string getQuotedToken(const char *&linePtr) {
  skipSpaces(linePtr);
  std::string res;
  if (*linePtr != '"') {
    return res;
  }
  ++linePtr;  // skip opening quote
  while (*linePtr && *linePtr != '"') {
    if (*linePtr == '\\') {
      // escaped char
      ++linePtr;
      if (*linePtr == '\0') {
        break;
      }
    }

    res += *linePtr++;
  }
  if (*linePtr != '"') {
    res = "";  // error: no closing quote
  } else {
    ++linePtr;  // skip closing quote
  }

  return res;
}

void parseTemplateLine(std::string lineStr,
                       unsigned int &line,  std::string &templateClass, std::string &templateSymbol, std::string &templateName, std::vector<std::pair<std::string,std::string>> &otherTokens) {

  // TEMPLATE 1 AA/Cya/Cya/ NATREPLACE=AA/A COMMENT=comment
  // FULLNAME=fullname CATEGORY=cat UNIQUEID=uniqueid CASNUMBER=xxxx,
  // COLLABORATOR=col, PROTECTION=prot

  // other attributes are allowed.  We capture them and ignore them, except
  // for writing them back out

  if (lineStr.substr(0, 9) != "TEMPLATE ") {
    std::ostringstream errout;
    errout << "Expected \"TEMPLATE\" at line  " << line;
    throw FileParseException(errout.str());
  }

  const char *linePtr = lineStr.c_str() + 9;
  std::string token = getToken(linePtr, ' ');  // Template ID
  if (token.empty()) {
    std::ostringstream errout;
    errout << "Expected a Template ID at line  " << line;
    throw FileParseException(errout.str());
  }
  // get the class and template symbol and name

  token = getToken(linePtr, ' ');  // Template ID
  if (token.empty()) {
    std::ostringstream errout;
    errout
        << "Type/Name(s) string of the form \"AA/Gly/G/\" was not found at line  "
        << line;
    throw FileParseException(errout.str());
  }

  // get the class and template names from the token
  std::vector<std::string> subTokens;
  boost::algorithm::split(subTokens, token, boost::algorithm::is_any_of("/"));
  if (subTokens.size() < 3 ){
    std::ostringstream errout;
    errout << "Type/Name(s) string is not of the form \"AA/Gly/G/\" at line  "
           << line;
    throw FileParseException(errout.str());
  }

  templateClass = subTokens[0];
  templateName = subTokens[1];
  templateSymbol = subTokens[2];

  // now parse attrs of the form ATTRNAME=VALUE or ATTRNAME="VALUE WITH
  // SPACES"

  otherTokens.clear();
  while (true) {
    std::string attrName = getToken(linePtr, '=');
    if (attrName.empty()) {
      break;
    }
    std::string attrValue;
    if (*linePtr == '"') {
      attrValue = getQuotedToken(linePtr);
    } else {
      attrValue = getToken(linePtr, ' ');
    }

    otherTokens.push_back(std::pair(attrName, attrValue));
    //templateMol->setProp(attrName, attrValue);
  }

  if (*linePtr != '\0') {
    std::ostringstream errout;
    errout
        << "extra characters at the end of a TEMPLATE definition line at line  "
        << line;
    throw FileParseException(errout.str());
  }

  return;
}

 SubstanceGroup &findMainSgroupForTemplate(
    RWMol &mol,
    std::string className, std::string templateName) {

  auto &sgroups = getSubstanceGroups(mol);

  for (unsigned int sgroupIdx = 0; sgroupIdx < sgroups.size(); ++sgroupIdx) {
    auto &sgroupToTest = sgroups[sgroupIdx];
    std::string sup;
    std::string sgroupAtomClass;
    if (sgroupToTest.getPropIfPresent("TYPE", sup) && sup == "SUP" &&
        sgroupToTest.getPropIfPresent("CLASS", sgroupAtomClass) &&
        sgroupAtomClass == className) {
      return sgroupToTest;
      break;
    }
  }

  throw FileParseException("No main sgroup found for template: " + className + "/" + templateName);
}

SubstanceGroup &findSgroupForSap(RWMol &mol, SubstanceGroup::AttachPoint &sap) {
  for (auto &sgroup : getSubstanceGroups(mol)) {
    std::string sgroupClass="";
    sgroup.getPropIfPresent("CLASS", sgroupClass);
    if (sgroupClass != "LGRP") {
      continue;
    } 

    auto &atoms = sgroup.getAtoms();
    if (std::find(atoms.begin(), atoms.end(), sap.lvIdx ) != atoms.end()) {
      return sgroup;
    }
  }

  throw FileParseException("No LGRP sgroup found for SAP");

}

void fixPhosphteConnectionPoints(MacroMol &macroMol) {
  //If there are two or more on the same attachAtomIdx, 
  // if will be random how there are assigned.  And this could lead to undesireable pairings.
  // (like R1-R1 or Al-Al)
  // The example that causes this is the Phospate group that has two leaving groups on the same P
  // atom.
  //
  // The macroMol reoutines are not supposed to know about the meanings of the attachment types,
  // 
  // this routine detects any such anomalies with P gropups and fixes them

  for (auto atom : macroMol.atoms()) {
    auto macroAtomInfo = atom->getMacroAtomInfo();
    if (!macroAtomInfo || macroAtomInfo->getMonomerClass() != MonomerClass::NucleicAcid ||
        macroAtomInfo->getSymbol() != "P") {
      continue;
    }

    for (auto bndI : boost::make_iterator_range(macroMol.getAtomBonds(atom))) {
      auto bond = (macroMol)[bndI];
      auto macroBondInfo = bond->getMacroBondInfo();
      if (!macroBondInfo) {
        continue;
      }

      for (auto &bondInfo : macroBondInfo->getBonds()) {
        if (bondInfo.beginAttachPt == 0 && bondInfo.endAttachPt == 0) {
          if (bond->getBeginAtomIdx() == atom->getIdx()) {  // first atom of bond is the P 
            bondInfo.beginAttachPt = 1;  
          } else {
            bondInfo.endAttachPt = 1;
          }
        } else if (bondInfo.beginAttachPt == 0 && bondInfo.endAttachPt == 0) {
          if (bond->getBeginAtomIdx() == atom->getIdx()) {  // first atom of bond is the P 
            bondInfo.beginAttachPt = 0;  
          } else {
            bondInfo.endAttachPt = 0;
          }
        }
      }
    }
  }
}


}  // namespace

using SCSRTemplateKey = std::pair<std::string, std::string>;

std::unique_ptr<MacroMol> SCSRUtils::MacroMolFromSCSRDataStream(
    std::istream &inStream, unsigned int &line,
    const v2::FileParsers::MolFileParserParams &params) {
  bool chiralityPossible = false;
  if (inStream.eof()) {
    return nullptr;
  }

  auto localParams = params;
  localParams.parsingSCSRMol = true;


  auto tempMol =
      RDKit::v2::FileParsers::MolFromMolDataStream(inStream, line, localParams);
  auto res = std::unique_ptr<RDKit::MacroMol>(new RDKit::MacroMol(tempMol));
  if (res->needsUpdatePropertyCache()) {
      res->updatePropertyCache(false);
  }


  // now get all of the templates
  // first these are placed in a vector of MacroMolTemplateBuilder items
  // so we can fix them up
  //then we convert them to MacroMolTemplates, which are const

  std::vector<std::unique_ptr<MacroMolTemplateBuilder>> templateBuilderLib;
  std::map<SCSRTemplateKey,  MacroMolTemplateBuilder *> bySymbolOrName;

  auto tempStr = RDKit::FileParserUtils::getV3000Line(&inStream, line);

  if (tempStr != "BEGIN TEMPLATE") {
    std::ostringstream errout;
    errout << "BEGIN TEMPLATE not found at line  " << line;
    throw RDKit::FileParseException(errout.str());
  }
  tempStr = RDKit::FileParserUtils::getV3000Line(&inStream, line);

  // TEMPLATE 1 AA/Cya/Cya/ NATREPLACE=AA/A COMMENT=comment
  // FULLNAME=fullname CATEGORY=cat UNIQUEID=uniqueid CASNUMBER=xxxx,
  // COLLABORATOR=col, PROTECTION=prot

  // other attributes are allowed.  We capture them and ignore them, except
  // for writing them back out

  while (tempStr.substr(0, 8) == "TEMPLATE") {
    auto templateMol = std::unique_ptr<RDKit::RWMol>(new RDKit::RWMol());
    std::string templateClass;
    std::string templateName;
    std::string templateSymbol;
    std::vector<std::pair<std::string,std::string>> otherTokens;

    parseTemplateLine( tempStr.c_str(), line, templateClass, templateSymbol, templateName, otherTokens);
    if (templateSymbol == "") {
      templateSymbol = templateName;
    }

    if (!ScsrClassToMonomerMolClass.contains(templateClass)) {
      std::ostringstream errout;
      errout << "Template class " << templateClass << " not recognized at line  " << line;
      throw RDKit::FileParseException(errout.str());
    }

    auto molComplete = false;
    RDKit::Conformer *conf = nullptr;
    try {
      unsigned int nAtoms = 0, nBonds = 0;
      bool expectMEND = false;
      molComplete = RDKit::FileParserUtils::ParseV3000CTAB(
          &inStream, line, templateMol.get(), conf, chiralityPossible, nAtoms,
          nBonds, params.strictParsing, expectMEND);
    } catch (RDKit::MolFileUnhandledFeatureException &e) {
      // unhandled mol file feature, show an error
      res.reset();
      delete conf;
      conf = nullptr;
      BOOST_LOG(rdErrorLog) << " Unhandled CTAB feature: '" << e.what()
                            << "'. Molecule skipped." << std::endl;

      if (!inStream.eof()) {
        tempStr = RDKit::getLine(inStream);
      }
      ++line;
      while (!inStream.eof() && !inStream.fail() &&
             tempStr.substr(0, 6) != "M  END" &&
             tempStr.substr(0, 4) != "$$$$") {
        tempStr = RDKit::getLine(inStream);
        ++line;
      }
      molComplete = !inStream.eof() || tempStr.substr(0, 6) == "M  END" ||
                    tempStr.substr(0, 4) == "$$$$";
    } catch (RDKit::FileParseException &e) {
      // catch our exceptions and throw them back after cleanup
      delete conf;
      conf = nullptr;
      throw e;
    }

    if (!molComplete) {
      delete conf;
      conf = nullptr;
      std::ostringstream errout;
      errout
          << "Problems encountered parsing Mol data, M  END or BEGIN TEMPLATE missing around line "
          << line;
      throw RDKit::FileParseException(errout.str());
    }

    if (templateMol) {
      auto tempParams = params;
      tempParams.sanitize = false;
      tempParams.removeHs = false;
      RDKit::FileParserUtils::finishMolProcessing(
          templateMol.get(), chiralityPossible, tempParams);

      std::unique_ptr<MacroMolTemplateBuilder> newTemplateBuilder(new MacroMolTemplateBuilder(
          *(templateMol.get()), monomerClassFromString(ScsrClassToMonomerMolClass[templateClass]),  templateName, templateSymbol, ""));
      newTemplateBuilder->setSubclass(templateClass);

      auto &molToUse = newTemplateBuilder->getMol();

      // find the main sgroup

      auto &sgroup = findMainSgroupForTemplate(molToUse, templateClass, templateSymbol);
      newTemplateBuilder->setMainGroup(sgroup.getAtoms());
     
      const auto &saps = sgroup.getAttachPoints();
      for (auto sap : saps) {

        // find the SUP for this SAP
        
        std::vector<unsigned int> atomIds;
        if (sap.lvIdx != -1) {
          auto &lvSgroup = findSgroupForSap(molToUse, sap);
          atomIds = lvSgroup.getAtoms();
        }
        MacroMolLeavingGroup leavingGroup(atomIds, sap.aIdx, sap.lvIdx, ScsrTypeToAttachPointType[sap.id] );
        newTemplateBuilder->addLeavingGroup(std::move(leavingGroup));
      }

      for (auto otherToken : otherTokens) {
        newTemplateBuilder->getMol().setProp(otherToken.first, otherToken.second);
      }

      if (bySymbolOrName.contains({templateClass, templateSymbol})) {
        std::ostringstream errout;
        errout
          << "Definition of Template " << templateClass << ":" << templateSymbol << " is previously defined: "
          << line;
        throw RDKit::FileParseException(errout.str());
      }


      if (bySymbolOrName.contains({templateClass, templateName}) && templateName != templateSymbol) {  // symbol and name could be the same
        std::ostringstream errout;
        errout
          << "Definition of Template " << templateClass << ":" << templateName << "is previously defined: "
          << line;
        throw RDKit::FileParseException(errout.str());
      }

      bySymbolOrName[SCSRTemplateKey(templateClass, templateSymbol)] = newTemplateBuilder.get();
      bySymbolOrName[SCSRTemplateKey(templateClass, templateName)] = newTemplateBuilder.get();
      templateBuilderLib.push_back(std::move(newTemplateBuilder));
    }

    tempStr = RDKit::FileParserUtils::getV3000Line(&inStream, line);
  }

  if (tempStr != "END TEMPLATE") {
    std::ostringstream errout;
    errout << "END TEMPLATE not found at line  " << line;
    throw RDKit::FileParseException(errout.str());
  }

  tempStr = RDKit::getLine(inStream);
  ++line;
  if (tempStr.substr(0, 6) != "M  END") {
    std::ostringstream errout;
    errout << "M  END not found at line  " << line;
    throw RDKit::FileParseException(errout.str());
  }

  // now validate the templates to the main mol
  // each atom in the mol that requires a template should match one,
  // and the ATTCHORD entries of each main atom must be consistent
  // with the template mol it matches.
  //
  // Also add the MacroAtomInfo, so that each macro atom is "proper"

  for (auto atom : res->atoms()) {
    std::string templateClass;
    std::string templateName;
    if (!atom->getPropIfPresent(common_properties::molAtomClass, templateClass) ||
        !atom->getPropIfPresent(common_properties::dummyLabel ,templateName)) {
          continue;  // not a macro atom
    }

    if(!bySymbolOrName.contains(SCSRTemplateKey(templateClass, templateName))) {
      std::ostringstream errout;
      errout
          << "Template for Macro atom " << templateClass << "/" <<templateName << "not found";

      throw RDKit::FileParseException(errout.str());
    }
    auto templateBuilder = bySymbolOrName[SCSRTemplateKey(templateClass, templateName)];

    auto templateSymbol =  templateBuilder->getSymbol();
    
    // add the MacroAtomInfo
    atom->setMacroAtomInfo(new MacroAtomInfo(templateSymbol, monomerClassFromString(ScsrClassToMonomerMolClass[templateClass])));

    auto mainSup =  findMainSgroupForTemplate(templateBuilder->getMol(), templateClass, templateName);

    std::vector<std::pair<unsigned int, std::string>> attchOrds;
    if (atom->hasProp(RDKit::common_properties::molAttachOrderTemplate)) {
      atom->getProp(RDKit::common_properties::molAttachOrderTemplate,
                    attchOrds);

      // check the ATTCHORD entries
      auto supAttachPoints = mainSup.getAttachPoints();
      for (const auto &attachOrd : attchOrds) {
        if (find_if(supAttachPoints.begin(), supAttachPoints.end(),
                    [&attachOrd](auto a) {
                      return a.id == attachOrd.second;
                    }) == supAttachPoints.end()) {
          std::ostringstream errout;
          errout << "No attachment point found in the template for atom "
                 << atom->getIdx() << " and attachment label "
                 << attachOrd.second;
          throw RDKit::FileParseException(errout.str());
        }
      }
    }
  }

  // now fix the bonds.   SCSR format has the connection points for macro
  // atoms defined in the main mol atoms, while MacroMol has that
  // information in the bonds.
  // It is possible for two macro atoms to have two or more bonds between them.
  // There are two cases currently known:  hydrogen bonds beetwen two BASE macro atoms,
  // and gem-disulfides (two cysteines adjacent to each other and also bonded by a 
  // disiulfide bond)
  // IN SCSR, for Hbonds, the attachment points can be either NOT used, and it is presumed
  // that the Hbond in the bond block indicates that here are two or more hbonds
  // between the BASES, OR, a single Hb attachment point is listed for each of the
  // BASE macro atoms, and two or three HB SAPs are listed in the correcponding 
  // Template definition.
  // 
  // For gem-disulfides, the two sets of attachPoints (Al or Br and Cx) are listed in 
  // each atom, and it is assumed that the two Cx attachPts are connected to each other.
  //
  // In MacroMols, the Hbond can either be a bond between the two atoms with NO 
  // bond info record (unspecific hbond for the atoms) or the bond can have 2 or more
  // entries in the bondInfo record for the specific H bonds
  // 
  // For gen-disulfides, the two connections are listed in the bondUInfo record
  //

  for (auto bond : res->bonds()) {
    auto atom1 = bond->getBeginAtom();
    auto atom2 = bond->getEndAtom();
    
    std::string atomClass1 = "";
    std::string dummyLabel1 = "";
    atom1->getPropIfPresent(common_properties::dummyLabel, dummyLabel1);
    atom1->getPropIfPresent(common_properties::molAtomClass, atomClass1);
  
    std::string atomClass2 = "";
    std::string dummyLabel2 = "";
    atom2->getPropIfPresent(common_properties::dummyLabel, dummyLabel2);
    atom2->getPropIfPresent(common_properties::molAtomClass, atomClass2);

    if (dummyLabel1 == "" && dummyLabel2 == "") {
      continue;  // neither atom is a macro atom - this is just a regular bond
    }

    std::vector<std::pair<unsigned int, std::string>> attchOrds1;
    std::vector<std::pair<unsigned int, std::string>> attchOrds2;
    atom1->getPropIfPresent(RDKit::common_properties::molAttachOrderTemplate,
                              attchOrds1);
    atom2->getPropIfPresent(RDKit::common_properties::molAttachOrderTemplate,
                                attchOrds2);

    
    std::vector<std::string> attchPoints1;
    std::vector<std::string> attchPoints2;

    for (auto attchOrd : attchOrds1) {
      if (attchOrd.first == atom2->getIdx()) {
        attchPoints1.push_back(attchOrd.second);
      }
    }
    for (auto attchOrd : attchOrds2) {
      if (attchOrd.first == atom1->getIdx()) {
        attchPoints2.push_back(attchOrd.second);
      }
    }

    if (attchPoints1.size() ==1 && attchPoints2.size() == 1) {  // most cases - just one bond between two atoms
      auto macroBondInfo = new MacroBondInfo();
      macroBondInfo->addBond(ScsrTypeToAttachPointType[attchPoints1[0]], ScsrTypeToAttachPointType[attchPoints2[0]], bond->getBondType());
      bond->setMacroBondInfo(macroBondInfo);
    } else if (attchPoints1.size() ==0 && attchPoints2.size() == 1) {  // regular atom bonded to macro atom
      auto macroBondInfo = new MacroBondInfo();
      macroBondInfo->addBond(-1, ScsrTypeToAttachPointType[attchPoints2[0]], bond->getBondType());
      bond->setMacroBondInfo(macroBondInfo);
    } else if (attchPoints1.size() ==1 && attchPoints2.size() == 0) {  // macro atom bonded to regular atom
      auto macroBondInfo = new MacroBondInfo();
      macroBondInfo->addBond(ScsrTypeToAttachPointType[attchPoints1[0]], -1, bond->getBondType());
      bond->setMacroBondInfo(macroBondInfo);
    } else if (attchPoints1.size() ==2 && attchPoints2.size() == 2) { 
      // two connections for both - must be like a gem-disullfide
      // there should be one Cx point on both sides
      
      std::string otherAttchPoint1 = "";
      std::string otherAttchPoint2 = "";
      if (attchPoints1[0] == "Cx" && attchPoints1[1] != "Cx") {
        otherAttchPoint1 = attchPoints1[1];
      } else if (attchPoints1[0] != "Cx" && attchPoints1[1] == "Cx") {
        otherAttchPoint1 = attchPoints1[0];
      } else {
        throw RDKit::FileParseException("Bond " + std::to_string(bond->getIdx()) + " has inconsistent attachOrds");
      }

      auto macroBondInfo = new MacroBondInfo();
      macroBondInfo->addBond(ScsrTypeToAttachPointType[otherAttchPoint1], ScsrTypeToAttachPointType[otherAttchPoint1], bond->getBondType());
      macroBondInfo->addBond(ScsrTypeToAttachPointType["Cx"], ScsrTypeToAttachPointType["Cx"], bond->getBondType());
      bond->setMacroBondInfo(macroBondInfo);

    } else if (attchPoints1.size() ==0 && attchPoints2.size() == 0 && bond->getBondType() == Bond::BondType::HYDROGEN) { 
      // no nothing for hbonds with no connection info
    } else {
      throw RDKit::FileParseException("Bond " + std::to_string(bond->getIdx()) + " has inconsistent attachOrds");
    }
  }

  // remove the atom's attchpts

  for (auto &atom : res->atoms()) {
    std::vector<std::pair<unsigned int, std::string>> attchOrds;
    if (atom->hasProp(RDKit::common_properties::molAttachOrderTemplate)) {
      atom->clearProp(RDKit::common_properties::molAttachOrderTemplate);
    }
  }

  for (auto &templateBuilder : templateBuilderLib) {
    auto &sGroups = RDKit::getSubstanceGroups(templateBuilder->getMol());

    // remove the SUPs from the template
    sGroups.clear();

    // make the actual templates from the Template builders

    res->addLocalTemplate(templateBuilder->build());
  }

  fixPhosphteConnectionPoints(*res);

  return res;
}

struct HbondQueryData {
      std::string d_smarts;
      std::string d_name;
      std::vector<SCSRUtils::HydrogenBondConnection> d_hBondConnections;
      std::shared_ptr<ROMol> dp_mol;
      HbondQueryData(const std::string &smarts, const std::string &name,
                      std::vector<SCSRUtils::HydrogenBondConnection> hBondConnections)
          : d_smarts(smarts),
              d_name(name),
              d_hBondConnections(hBondConnections),
              dp_mol{SmartsToMol(smarts)} {}
  };
  static const std::vector<HbondQueryData> hbondQueries = {
      {"[NH]1C=NC2=C1N=C([NH2])[NH]C2=O",
        "Guanine",
        {SCSRUtils::HydrogenBondConnection(7, true), SCSRUtils::HydrogenBondConnection(8, true),
        SCSRUtils::HydrogenBondConnection(10, false)}},
      {"[NH]1C=NC2=C1N=[CH]N=C2[NH2]",
        "Adenine",
        {SCSRUtils::HydrogenBondConnection(6, true), SCSRUtils::HydrogenBondConnection(7, false),
        SCSRUtils::HydrogenBondConnection(9, true)}},
      {"[NH]1C=CC([NH2])=NC1=O",
        "Cytosine",
        {SCSRUtils::HydrogenBondConnection(7, false), SCSRUtils::HydrogenBondConnection(5, false),
        SCSRUtils::HydrogenBondConnection(4, true)}},
      {"[NH]1C=CC(=O)[NH]C1=O",
        "Thyamine-Uracil",
        {SCSRUtils::HydrogenBondConnection(7, false), SCSRUtils::HydrogenBondConnection(5, true),
        SCSRUtils::HydrogenBondConnection(4, false)}},
      {"[NH]1C=NC2=C1N=C[NH]C2=O",
        "Inosine",
        {SCSRUtils::HydrogenBondConnection(7, true),
        SCSRUtils::HydrogenBondConnection(9, false)}}};

namespace {
void getBaseTemplates(MacroMol &macroMol, std::vector<const MacroMolTemplate *> &baseTemplates, std::vector<const MacroMolTemplate *> &nonBaseTemplates) {
 
  for (auto templateMol : macroMol.getLocalTemplateLibrary().getTemplates()) {
    if (getScsrClass(*templateMol) == "BASE") {
      baseTemplates.push_back(templateMol);
    } else {
      nonBaseTemplates.push_back(templateMol);
    }
  }

  return;
}

bool isAtomABaseReference(const Atom &atom,  std::vector<const MacroMolTemplate *> &baseTemplates) {
  auto macroAtomInfo = atom.getMacroAtomInfo();
  if (!macroAtomInfo) {
    return false;
  }
  auto atomClass = macroAtomInfo->getMonomerClass();
  auto atomSymbol = macroAtomInfo->getSymbol();
  auto baseTemplateItr = std::find_if(baseTemplates.begin(), baseTemplates.end(),
    [atomClass, atomSymbol](auto &bt) {
        return bt->getMonomerClass() == atomClass && bt->getSymbol() == atomSymbol; }
    );

  if (baseTemplateItr != baseTemplates.end()) {
      return true;
  }

  return false;
}

void RemoveBaseHBonds(std::unique_ptr<MacroMol> &macroMol,  bool removeHBondsFlag ){

  // see if we have anything to do

  std::vector<const MacroMolTemplate *> baseTemplates;
  std::vector<const MacroMolTemplate *> nonBaseTemplates;
  getBaseTemplates(*(macroMol.get()), baseTemplates, nonBaseTemplates);

  if (baseTemplates.empty()) {
    return;  // nothing to do
  }
 
  std::unique_ptr<MacroMolTemplateLibrary> newTemplateLibrary(new MacroMolTemplateLibrary());

  // preserve all templates that are NOT BASE class


  for (auto localTemplatePtr : nonBaseTemplates) {
    newTemplateLibrary->addTemplate(std::unique_ptr<MacroMolTemplate>(new MacroMolTemplate(*localTemplatePtr)));
  }

  bool newTemplateLibraryNeeded = false;
  for (auto localTemplate : baseTemplates) {
    auto leavingGroups = localTemplate->getLeavingGroups();

    // see if it has any  Hb Hb1, Hb2 or Hb3.  If not, nothing needs to be done to this one


    if (std::find_if(leavingGroups.begin(), leavingGroups.end(), 
          [](auto lg) {
                  return AttachPointTypeToSCSRType[(AttachPointType)lg.attachPoint].substr(0,2) == "Hb";}
              ) == leavingGroups.end()) {
      newTemplateLibrary->addTemplate(std::unique_ptr<MacroMolTemplate>(new MacroMolTemplate(*localTemplate)));
      continue;
    }
  
    std::unique_ptr<MacroMolTemplateBuilder> newTemplateBuilder(new MacroMolTemplateBuilder(localTemplate->getMol(), 
                          MonomerClass::NucleicAcid,
                          localTemplate->getName(),localTemplate->getSymbol(),
                          localTemplate->getOriginalData()));
    newTemplateBuilder->setMainGroup(localTemplate->getMainAtomIdxs());

    // copy the leaving groups for this template - but not Hb if it is found

    for (auto &leavingGroup : leavingGroups) {
      if (AttachPointTypeToSCSRType[(AttachPointType)leavingGroup.attachPoint].substr(0,2) != "Hb") {
        MacroMolLeavingGroup newLeavingGroup(leavingGroup);
        newTemplateBuilder->addLeavingGroup(std::move(newLeavingGroup));
      }
    }

    newTemplateLibrary->addTemplate(newTemplateBuilder->build());
    newTemplateLibraryNeeded = true;
  }

  if (newTemplateLibraryNeeded) {
    macroMol->setLocalTemplateLibrary(std::move(newTemplateLibrary));
  }


  if (removeHBondsFlag) {
      // remove any hbond connection info from the main mol bonds

    macroMol->beginBatchEdit();
    for (auto &bond : macroMol->bonds()) {

      if (bond->getBondType() != Bond::BondType::HYDROGEN) {
        continue;
      }

      if (isAtomABaseReference(*bond->getBeginAtom(), baseTemplates) || isAtomABaseReference(*bond->getEndAtom(), baseTemplates)) {
        macroMol->removeBond(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
      }
    }
    macroMol->commitBatchEdit();
  }

  else {

    for (auto &bond : macroMol->bonds()) {

      auto macroBondInfo = bond->getMacroBondInfo();

      if (!macroBondInfo) {
        continue;
      }

      if (isAtomABaseReference(*bond->getBeginAtom(), baseTemplates) || isAtomABaseReference(*bond->getEndAtom(), baseTemplates)) {
        // if any of the bond items is an "Hb*" type, they should all be 
        if (!macroBondInfo->getBonds().empty()) {
          auto bondInfo = macroBondInfo->getBonds()[0];
          if (AttachPointTypeToSCSRType[(AttachPointType)bondInfo.beginAttachPt].substr(0,2) == "Hb" ||
              AttachPointTypeToSCSRType[(AttachPointType)bondInfo.endAttachPt].substr(0,2) == "Hb" ) {
            bond->setBondType(Bond::BondType::HYDROGEN);
            bond->setMacroBondInfo(nullptr);
          }
        }
      }
    }
  }
}

void AddBaseHBonds(std::unique_ptr<MacroMol> &macroMol) {

  
  RemoveBaseHBonds(macroMol, false);

    // see of we have anything to do

  std::vector<const MacroMolTemplate *> baseTemplates;
  std::vector<const MacroMolTemplate *> nonBaseTemplates;
  getBaseTemplates(*(macroMol.get()), baseTemplates, nonBaseTemplates);

  if (baseTemplates.empty()) {
    return;  // nothing to do
  }
 
  
  std::unique_ptr<MacroMolTemplateLibrary> newTemplateLibrary(new MacroMolTemplateLibrary());

  
  // make sure any BASE templates have the HBond definitions

  std::map<std::string, std::vector<bool>> donorFlagsMap;

  for (auto localTemplate : nonBaseTemplates) {
    newTemplateLibrary->addTemplate(std::unique_ptr<MacroMolTemplate>(new MacroMolTemplate(*localTemplate)));
  }

  for (auto localTemplate : baseTemplates) {
    auto leavingGroups = localTemplate->getLeavingGroups();
    auto templateSymbol = localTemplate->getSymbol();

    std::unique_ptr<MacroMolTemplateBuilder> newTemplateBuilder(new MacroMolTemplateBuilder(localTemplate->getMol(), 
                          MonomerClass::NucleicAcid,
                          localTemplate->getName(),templateSymbol,
                          localTemplate->getOriginalData()));
    newTemplateBuilder->setMainGroup(localTemplate->getMainAtomIdxs());

    // copy the leaving groups for this template - but not Hb if it is found

    for (auto &leavingGroup : leavingGroups) {
      if (leavingGroup.attachPoint != AttachPointType::Hb) {
        MacroMolLeavingGroup newLeavingGroup(leavingGroup);
        newTemplateBuilder->addLeavingGroup(std::move(newLeavingGroup));
      }
    }

    // now find the hb groups to be added according to the chemical 
    // structure of the template

    auto &mol = newTemplateBuilder->getMol();
    if (mol.needsUpdatePropertyCache()) {
      mol.updatePropertyCache(false);
    }
    
    for (const auto &querySmi : hbondQueries) {
      RDKit::SubstructMatchParameters params;
      auto match = SubstructMatch(mol, *querySmi.dp_mol, params);

      if (!match.empty()) {
        donorFlagsMap[templateSymbol] = std::vector<bool>();
        for (const auto &hBondConnection : querySmi.d_hBondConnections) {
          donorFlagsMap[templateSymbol].push_back(hBondConnection.d_isDonor);
          auto attachPt = ScsrTypeToAttachPointType["Hb" + std::to_string(donorFlagsMap[templateSymbol].size())];
          MacroMolLeavingGroup newLeavingGroup(std::vector<unsigned int>(),
                        match[0][hBondConnection.d_templateAtomIdx].second,
                        -1, attachPt);
          newTemplateBuilder->addLeavingGroup(std::move(newLeavingGroup));
        }          
        break; // only one match
      }   
    }
    newTemplateLibrary->addTemplate(newTemplateBuilder->build());
  }

 
  //Now find all bonds in the main mol that refer to atoms that
  // are based on BASE templates, and fix them


  for (auto &bond : macroMol->bonds()) {
    if (bond->getBondType() != Bond::BondType::HYDROGEN) { 
      continue;
    }

    if (isAtomABaseReference(*bond->getBeginAtom(), baseTemplates) && isAtomABaseReference(*bond->getEndAtom(), baseTemplates)) {
      // add a new macroBondInfo

      bond->setMacroBondInfo(new MacroBondInfo());
      auto macroBondInfo = bond->getMacroBondInfo();
      auto atom1 = bond->getBeginAtom();
      auto atom2 = bond->getEndAtom();
        
      makeHydrogenBonds(donorFlagsMap[atom1->getMacroAtomInfo()->getSymbol()],
                        donorFlagsMap[atom2->getMacroAtomInfo()->getSymbol()],
                        *macroBondInfo);

    }
  }
  
  macroMol->setLocalTemplateLibrary(std::move(newTemplateLibrary));
  
}
}  //namespace

void SCSRUtils::AdjustBaseHBonds(std::unique_ptr<MacroMol> &macroMol,  SCSRUtils::SCSRBaseHbondOptions options) {
  switch (options) {
    case SCSRUtils::SCSRBaseHbondOptions::RemoveHBonds:
      RemoveBaseHBonds(macroMol, true );  //remove HBond bonds from BASE 
      break;

    case SCSRUtils::SCSRBaseHbondOptions::RemovedHbondLeavingGroups:
      RemoveBaseHBonds(macroMol, false );  //keep the bonds, remove the Connections anbd SAPS
      break;

    case SCSRUtils::SCSRBaseHbondOptions::FullHBonds:
      AddBaseHBonds(macroMol);   // Add Hbonds and connections to MacroBondInfo and templates
      break;
  }
}


//------------------------------------------------
//
//  Read an SCSR molecule from a string
//
//------------------------------------------------
std::unique_ptr<RDKit::MacroMol> SCSRUtils::MacroMolFromSCSRBlock(
    const std::string &molBlock,
    const RDKit::v2::FileParsers::MolFileParserParams &params) {
  std::istringstream inStream(molBlock);
  unsigned int line = 0;
  return MacroMolFromSCSRDataStream(inStream, line, params);
}

//------------------------------------------------
//
//  Read an SCSR molecule from a file
//
//------------------------------------------------
std::unique_ptr<RDKit::MacroMol> SCSRUtils::MacroMolFromSCSRFile(
    const std::string &fName,
    const RDKit::v2::FileParsers::MolFileParserParams &params) {
  std::ifstream inStream(fName.c_str());
  if (!inStream || (inStream.bad())) {
    std::ostringstream errout;
    errout << "Bad input file " << fName;
    throw RDKit::BadFileException(errout.str());
  }
  if (!inStream.eof()) {
    unsigned int line = 0;
    return MacroMolFromSCSRDataStream(inStream, line, params);
  } else {
    throw RDKit::BadFileException("Empty file found");
  }
}

// std::unique_ptr<RDKit::RWMol> SCSRUtils::MolFromSCSRDataStream(
//     std::istream &inStream, unsigned int &line,
//     const RDKit::v2::FileParsers::MolFileParserParams &molFileParserParams,
//     const RDKit::MolFromMacroMolParams &molFromMacroMolParams) {
//   auto macroMol = SCSRUtils::MacroMolFromSCSRDataStream(
//       inStream, line, molFileParserParams);
//   return MolFromMacroMol(
//     *macroMol.get(), macroMol->getLocalTemplateLibrary(), molFromMacroMolParams);

// }

std::unique_ptr<RDKit::RWMol> SCSRUtils::MolFromSCSRBlock(
    const std::string &molBlock,
    const RDKit::v2::FileParsers::MolFileParserParams &molFileParserParams,
    const RDKit::MolFromMacroMolParams &molFromMacroMolParams,
    const RDKit::SCSRUtils::SCSRBaseHbondOptions scsrBaseHbondOptions) {
  auto macroMol = SCSRUtils::MacroMolFromSCSRBlock(molBlock, molFileParserParams);
  
  if (scsrBaseHbondOptions == RDKit::SCSRUtils::SCSRBaseHbondOptions::RemovedHbondLeavingGroups) {
    throw ValueErrorException("Full atomonistic output must have full Hbond treatment OR no hbonds for BASE tempates");
  }
  AdjustBaseHBonds(macroMol, scsrBaseHbondOptions);

  return MolFromMacroMol(
    *macroMol.get(), macroMol->getLocalTemplateLibrary(),  molFromMacroMolParams);

}


std::unique_ptr<RDKit::RWMol> SCSRUtils::MolFromSCSRFile(
    const std::string &fName,
    const RDKit::v2::FileParsers::MolFileParserParams &molFileParserParams,
    const RDKit::MolFromMacroMolParams &molFromMacroMolParams,
    const RDKit::SCSRUtils::SCSRBaseHbondOptions scsrBaseHbondOptions) {
  auto macroMol = SCSRUtils::MacroMolFromSCSRFile(fName, molFileParserParams);

  if (scsrBaseHbondOptions == RDKit::SCSRUtils::SCSRBaseHbondOptions::RemovedHbondLeavingGroups) {
    throw ValueErrorException("Full atomonistic output must have full Hbond treatment OR no hbonds for BASE tempates");
  }
  AdjustBaseHBonds(macroMol, scsrBaseHbondOptions);

  return MolFromMacroMol(
      *macroMol.get(), macroMol->getLocalTemplateLibrary(), molFromMacroMolParams);

 
}

std::string SCSRUtils::MacroMolToSCSRMolBlock(MacroMol &macroMol,
                                   const RDKit::MolWriterParams &params,
                                   int confId,
                                   bool useMacroNamesForAtomNames/* if set the macro names (eg 3 letter codes) are used for names */) {
  RDKit::Utils::LocaleSwitcher switcher;
  std::string res = "";
  auto localParams = params;
  localParams.forceV3000 = true;
  localParams.writeEndMolLine = false;

  // make a copy of the main macromol so we can modify it for output without
  // affecting the original

  auto localMol = std::unique_ptr<MacroMol>(new MacroMol(macroMol));

  AdjustBaseHBonds(localMol, SCSRUtils::SCSRBaseHbondOptions::RemovedHbondLeavingGroups);

  // move the MacroAtomInfo to be props of the atoms

  for (Atom *atom :localMol->atoms()) {
    auto macroAtomInfo = atom->getMacroAtomInfo();
    if (macroAtomInfo == nullptr) {
      continue;  // not a macro atom
    }

    auto TemplateForAtom = macroMol.getLocalTemplateLibrary().getBySymbol(
          macroAtomInfo->getMonomerClass(),macroAtomInfo->getSymbol());

    atom->setProp(common_properties::molAtomClass, getScsrClass(*TemplateForAtom));
    if (useMacroNamesForAtomNames) {
      atom->setProp(common_properties::dummyLabel, TemplateForAtom->getName());
    } else {
      atom->setProp(common_properties::dummyLabel, macroAtomInfo->getSymbol());
    }

    atom->setMacroAtomInfo(nullptr);
  }

  // move the connection point info from the bonds to the atoms for the main
  // mol

  for (auto &bond: localMol->bonds()) {
    auto macroBondInfo = bond->getMacroBondInfo();
    if (!macroBondInfo) {
      continue;
    }

    Atom *atom1 = bond->getBeginAtom();
    Atom *atom2 = bond->getEndAtom();

    std::vector<std::pair<unsigned int, std::string>> attchOrds;

    atom1->getPropIfPresent(RDKit::common_properties::molAttachOrderTemplate,
                              attchOrds);
                      
    for (const auto &bondInfo : macroBondInfo->getBonds()) {
      attchOrds.emplace_back(atom2->getIdx(), AttachPointTypeToSCSRType[(AttachPointType)bondInfo.beginAttachPt]);
    }

    atom1->setProp(RDKit::common_properties::molAttachOrderTemplate,
                            attchOrds);

    attchOrds.clear();
    atom2->getPropIfPresent(RDKit::common_properties::molAttachOrderTemplate,
                              attchOrds);
                      
    for (const auto &bondInfo : macroBondInfo->getBonds()) {
      attchOrds.emplace_back(atom1->getIdx(), AttachPointTypeToSCSRType[(AttachPointType)bondInfo.endAttachPt]);
    }

    atom2->setProp(RDKit::common_properties::molAttachOrderTemplate,
                            attchOrds);

    if ( macroBondInfo->getNumBonds() > 0) {
      bond->setBondType((Bond::BondType)macroBondInfo->getBond(0).bondType);
    }
    bond->setMacroBondInfo(nullptr);
  }

  res += RDKit::MolToMolBlock(*(localMol), localParams, confId);

  // for each template, create the main SUP and a SUP for each leaving group
  
  res += "M  V30 BEGIN TEMPLATE\n";
  unsigned int templateId = 0;
  for (const auto macroMolTemplate : macroMol.getLocalTemplateLibrary().getTemplates()) {
    // templates need to have the main atoms and the leaving groups converted to
    /// SGROUPS

    auto scsrTemplateClass = getScsrClass(*macroMolTemplate);
    auto templateSymbol = macroMolTemplate->getSymbol();
    auto templateName = macroMolTemplate->getName();

    
    res += "M  V30 TEMPLATE " + std::to_string(templateId + 1) + " " +
           scsrTemplateClass + "/" +  templateName + "/" +  templateSymbol + "/";
   

    for (auto &propName : macroMolTemplate->getMol().getPropList()) {
      if (propName != common_properties::molAtomClass &&
          propName != common_properties::templateNames && propName[0] != '_') {
        std::string propVal;
        macroMolTemplate->getMol().getProp<std::string>(propName, propVal);
        res += " " + propName + "=\"" + propVal + "\"";
      }
    }

    res += "\n";

    // copy template mol so we can modify it for output without affecting the
    // original, if needed

    // work on a copy: prepareMol mutates the mol, and the template is const
    std::unique_ptr<RWMol> tMol(new RWMol(macroMolTemplate->getMol()));
    RWMol *molToUse = tMol.get();

    // set the main sup group

    molToUse->clearSubstanceGroups();
    auto &sgroups = RDKit::getSubstanceGroups(*molToUse);
    sgroups.emplace_back(molToUse,"SUP");

    auto &sgroup = sgroups.back();

    for (auto atomId : macroMolTemplate->getMainAtomIdxs()) {
      sgroup.addAtomWithIdx(atomId);
    }

    sgroup.setProp("CLASS", scsrTemplateClass);
    sgroup.setProp("LABEL", templateSymbol);

    // for each leaving group, add the SAP to the main sgroup

    auto &saps = sgroup.getAttachPoints();
    for (auto &leavingGroup : macroMolTemplate->getLeavingGroups()) {
      RDKit::SubstanceGroup::AttachPoint attachPoint;
      attachPoint.aIdx = leavingGroup.attachAtomIdx;
      attachPoint.lvIdx = leavingGroup.leavingAtomIdx;
      attachPoint.id = AttachPointTypeToSCSRType[(AttachPointType)leavingGroup.attachPoint];\
      saps.push_back(attachPoint);

      // add a XBOND for each SAP, if the leaving group is not 0

      if (leavingGroup.leavingAtomIdx != -1) {
        auto bond = molToUse->getBondBetweenAtoms(attachPoint.aIdx, attachPoint.lvIdx);
        sgroup.addBondWithIdx(bond->getIdx());

         // also add a SUP for each SAP
        
        sgroups.emplace_back(molToUse,"SUP");
        auto &lvSgroup = sgroups.back();

        lvSgroup.setProp("CLASS", "LGRP");
        lvSgroup.setProp("LABEL", templateSymbol);
        for (auto atomId : leavingGroup.atomIdxs) {
          lvSgroup.addAtomWithIdx(atomId);
        }
        lvSgroup.addBondWithIdx(bond->getIdx());

      }
    }

    boost::dynamic_bitset<> aromaticBonds(molToUse->getNumBonds());
    RDKit::prepareMol(*molToUse, localParams, aromaticBonds);
    res += FileParserUtils::getV3000CTAB(*molToUse, aromaticBonds, confId,
                                         localParams.precision);
    ++templateId;
  }

  res += "M  V30 END TEMPLATE\n";
  res += "M  END\n";

  return res;
}

void SCSRUtils::MacroMolToSCSRMolFile(RDKit::MacroMol &macroMol, const std::string &fName,
                           const RDKit::MolWriterParams &params, int confId,
        bool useMacroNamesForAtomNames /* if set the macro names (eg 3 letter codes) are used for names */) {
  auto *outStream = new std::ofstream(fName.c_str());
  if (!(*outStream) || outStream->bad()) {
    delete outStream;
    std::ostringstream errout;
    errout << "Bad output file " << fName;
    throw BadFileException(errout.str());
  }
  std::string outString = SCSRUtils::MacroMolToSCSRMolBlock(macroMol, params, confId, useMacroNamesForAtomNames);
  *outStream << outString;
  delete outStream;
}

}  // namespace RDKit
