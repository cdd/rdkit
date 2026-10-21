//
//  Copyright (C) 2025 Tad Hurst and other RDKit contributors
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
#include <streambuf>

#include "RDGeneral/test.h"
#include <catch2/catch_all.hpp>
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>

#include "FileParsers.h"
#include "FileParserUtils.h"
#include "MolSGroupParsing.h"
#include <RDGeneral/StreamOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/SubstanceGroup.h>

#include <GraphMol/FileParsers/FileWriters.h>
#include <GraphMol/Atropisomers.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Atom.h>

#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/MacroMol.h>
#include <GraphMol/FileParsers/MolToMacroMol.h>
#include <GraphMol/FileParsers/SCSRUtils.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/lexical_cast.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <filesystem>
using namespace RDKit;
using namespace v2::FileParsers;

class ScsrMolTest {
 public:
 public:
  bool generateExpectedFiles = false;

  enum class ExpectedStatus {
    Success,
    FailsScsrParsing,
    FailsMolConversion
  };

  class ScsrTest {
   public:
    std::string fileName;
    unsigned int totalAtomCount;
    unsigned int totalBondCount;
    unsigned int sgroupCount;
    unsigned int totalQueryAtomCount;
    unsigned int totalQueryBondCount;
    unsigned int querySgroupCount;
    unsigned int scsrAtomCount;
    unsigned int scsrBondCount;
    ExpectedStatus expectedStatus;

    SCSRUtils::SCSRBaseHbondOptions scsrBaseHbondOptions;
    std::vector<std::pair<unsigned int, Atom::ChiralType>> chiralChecks;
    std::vector<std::pair<unsigned int, Atom::ChiralType>> chiralChecksQuery;

    ScsrTest(std::string fileNameInit, ExpectedStatus expectedStatusInit,
             SCSRUtils::SCSRBaseHbondOptions scsrBaseHbondOptions,
             unsigned int totalAtomCountInit, unsigned int totalBondCountInit,
             unsigned int sgroupCountInit, unsigned int totalQueryAtomCountInit,
             unsigned int totalQueryBondCountInit,
             unsigned int querySgroupCountInit,
             unsigned int scsrAtomCountInit = 0, unsigned int scsrBondCountInit = 0,  // 0 indicates that atom/bond count does not change on regeneration
             std::vector<std::pair<unsigned int, Atom::ChiralType>>
                 chiralChecksInit = {},
             std::vector<std::pair<unsigned int, Atom::ChiralType>>
                 chiralChecksQueryInit = {})
        : fileName(fileNameInit),
          totalAtomCount(totalAtomCountInit),
          totalBondCount(totalBondCountInit),
          sgroupCount(sgroupCountInit),
          totalQueryAtomCount(totalQueryAtomCountInit),
          totalQueryBondCount(totalQueryBondCountInit),
          querySgroupCount(querySgroupCountInit),
          scsrAtomCount(scsrAtomCountInit),
          scsrBondCount(scsrBondCountInit),
          expectedStatus(expectedStatusInit),
          scsrBaseHbondOptions(scsrBaseHbondOptions),
          chiralChecks(chiralChecksInit),
          chiralChecksQuery(chiralChecksQueryInit) {};
  };

  class ScsrMakeTest {
   public:
    std::string fileName;
    std::string templateFileName;
    unsigned int atomCount;
    unsigned int bondCount;
    unsigned int templateCount;
    bool expectedStatus;

    ScsrMakeTest(std::string fileNameInit, std::string templateFileNameInit,
                 unsigned int atomCountInit = 0, unsigned int bondCountInit = 0,
                 unsigned int templateCountInit = 0,
                 bool expectedStatusInit = true)
        : fileName(fileNameInit),
          templateFileName(templateFileNameInit),
          atomCount(atomCountInit),
          bondCount(bondCountInit),
          templateCount(templateCountInit),
          expectedStatus(expectedStatusInit) {};
  };

 public:
  std::list<ScsrMolTest::ScsrTest> ScsrTests;

  ScsrMolTest() {
    ScsrTests = {
        ScsrMolTest::ScsrTest("DnaBadPairs.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 84, 94, 10, 80,
                              90, 10),

        ScsrMolTest::ScsrTest("ValenceErrorScsr2.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 29, 29, 3, 26, 26, 3),

  
        ScsrMolTest::ScsrTest("DNASlurpErrorImport.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 81, 90, 11, 79, 88,
                              11, 0,0,
                              {{0, Atom::ChiralType::CHI_TETRAHEDRAL_CW},
                               {22, Atom::ChiralType::CHI_TETRAHEDRAL_CW},
                               {41, Atom::ChiralType::CHI_TETRAHEDRAL_CW},
                               {61, Atom::ChiralType::CHI_TETRAHEDRAL_CW}},
                              {{0, Atom::ChiralType::CHI_TETRAHEDRAL_CW},
                               {22, Atom::ChiralType::CHI_TETRAHEDRAL_CW},
                               {41, Atom::ChiralType::CHI_TETRAHEDRAL_CW},
                               {61, Atom::ChiralType::CHI_TETRAHEDRAL_CW}}),
        ScsrMolTest::ScsrTest("DNASlurpErrorSketch.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 84, 93, 12, 82, 91,
                              12, 0,  0,
                              {{51, Atom::ChiralType::CHI_TETRAHEDRAL_CCW},
                               {60, Atom::ChiralType::CHI_TETRAHEDRAL_CCW},
                               {68, Atom::ChiralType::CHI_TETRAHEDRAL_CCW},
                               {76, Atom::ChiralType::CHI_TETRAHEDRAL_CCW}},
                              {{51, Atom::ChiralType::CHI_TETRAHEDRAL_CCW},
                               {60, Atom::ChiralType::CHI_TETRAHEDRAL_CCW},
                               {68, Atom::ChiralType::CHI_TETRAHEDRAL_CCW},
                               {76, Atom::ChiralType::CHI_TETRAHEDRAL_CCW}}),

        ScsrMolTest::ScsrTest("ValenceErrorScsr.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 38, 39, 3, 35, 36, 3),
        ScsrMolTest::ScsrTest("ValenceErrorScsr2.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 29, 29, 3, 26, 26, 3),
        ScsrMolTest::ScsrTest("RiboseFullname.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 45, 49, 6, 43, 47, 6),
        ScsrMolTest::ScsrTest("testSCSR.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 64, 66, 8, 57,
                              59, 8),

        ScsrMolTest::ScsrTest("ModifiedPeptide2.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 438, 444, 50, 407,
                              413, 50),
        ScsrMolTest::ScsrTest("DnaBadPairs_NoCh.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 84, 94, 10, 80, 90,
                              10),
        ScsrMolTest::ScsrTest("DnaBadPairs.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 84, 94, 10, 80,
                              90, 10),
        ScsrMolTest::ScsrTest("DnaTest.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 254, 300, 34,
                              250, 296, 34),
        ScsrMolTest::ScsrTest("wobblePairs2.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 169, 196, 22,
                              165, 192, 22),
        ScsrMolTest::ScsrTest("wobblePairs.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 169, 196, 22,
                              165, 192, 22),
        ScsrMolTest::ScsrTest("KellanRNA.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 299, 353, 40,
                              295, 349, 40),
        ScsrMolTest::ScsrTest("DnaTest2.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 83, 97, 10, 79,
                              93, 10),
        ScsrMolTest::ScsrTest("DnaTest3.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 165, 194, 22,
                              161, 190, 22),
        ScsrMolTest::ScsrTest("KellanError.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 244, 263, 31,
                              236, 255, 31),

        ScsrMolTest::ScsrTest("TestRNA2_fixed.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 106, 117, 15,
                              104, 115, 15),


        ScsrMolTest::ScsrTest("DnaTest.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::RemoveHBonds, 254, 282, 34, 250,
                              278, 34),
    
        ScsrMolTest::ScsrTest("wobblePairs2.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::RemoveHBonds, 169, 188, 22,
                              165, 184, 22),
        ScsrMolTest::ScsrTest("wobblePairs.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::RemoveHBonds, 169, 188, 22, 165,
                              184, 22),
        ScsrMolTest::ScsrTest("DnaBadPairs.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::RemoveHBonds, 84, 92, 10, 80, 88,
                              10),
        ScsrMolTest::ScsrTest("DnaTest2.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::RemoveHBonds, 83, 91, 10, 79, 87,
                              10),
        ScsrMolTest::ScsrTest("DnaTest3.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::RemoveHBonds, 165, 183, 22, 161,
                              179, 22),
        ScsrMolTest::ScsrTest("KellanError.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::RemoveHBonds, 244, 263, 31, 236,
                              255, 31),
        ScsrMolTest::ScsrTest("TestRNA2_fixed.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::RemoveHBonds, 106, 117, 15, 104,
                              115, 15),

        ScsrMolTest::ScsrTest(
            "TrastuzumabMaxPlus3Register.mol", ExpectedStatus::Success,
            SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 7606, 7793, 925, 7080, 7267, 925),
        ScsrMolTest::ScsrTest(
            "TrastuzumabMaxRegister.mol", ExpectedStatus::Success,
            SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 7576, 7763, 922, 7053, 7240, 922),
        ScsrMolTest::ScsrTest("Mixed.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 51, 54, 8, 51,
                              54, 8),
        ScsrMolTest::ScsrTest("CrossLink.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 47, 48, 8, 45,
                              46, 8),
        ScsrMolTest::ScsrTest("cyclic.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 45, 47, 8, 45,
                              47, 8),

        ScsrMolTest::ScsrTest("Triplet.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 30, 30, 3, 27,
                              27, 3),
        ScsrMolTest::ScsrTest("FromBioviaDoc.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 27, 26, 5, 25,
                              24, 5),


       
        ScsrMolTest::ScsrTest(
            "badAtomName.mol", ExpectedStatus::FailsScsrParsing,
            SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 0, 0, 0, 0, 0, 0),
        ScsrMolTest::ScsrTest("badClass.mol", ExpectedStatus::FailsScsrParsing,
                              SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 0, 0, 0, 0, 0,
                              0),
        ScsrMolTest::ScsrTest(
            "badClassTemplate.mol", ExpectedStatus::FailsScsrParsing,
            SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 0, 0, 0, 0, 0, 0),
        ScsrMolTest::ScsrTest(
            "badMissingTemplate.mol", ExpectedStatus::FailsScsrParsing,
            SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 0, 0, 0, 0, 0, 0),
        ScsrMolTest::ScsrTest("obj3dTest.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 27, 26, 5, 25,
                              24, 5),
        ScsrMolTest::ScsrTest("obj3dTest.mol", ExpectedStatus::Success,
                              SCSRUtils::SCSRBaseHbondOptions::RemoveHBonds, 27, 26, 5, 25,
                              24, 5),
        ScsrMolTest::ScsrTest(
            "obj3dFoundTwice.mol", ExpectedStatus::FailsScsrParsing,
            SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 27, 26, 0, 27, 26, 0),
        ScsrMolTest::ScsrTest(
            "SgroupFoundTwice.mol", ExpectedStatus::FailsScsrParsing,
            SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 0, 0, 0, 0, 0, 0),
    };
  }

  void testScsrFiles(const ScsrTest *ScsrTest) {
    BOOST_LOG(rdInfoLog) << "testing scsr  files" << std::endl;

    INFO(ScsrTest->fileName);

    std::string rdbase = getenv("RDBASE");
    std::string fName = rdbase +
                        "/Code/GraphMol/FileParsers/test_data/macromols/" +
                        ScsrTest->fileName;
    std::string fnameOutBase = rdbase +
                           "/build/Testing/Temporary/" +
                           ScsrTest->fileName;
    if (ScsrTest->scsrBaseHbondOptions == SCSRUtils::SCSRBaseHbondOptions::RemoveHBonds) {
      fnameOutBase += ".NoHbonds";
    }
  
    std::string fOutName = fnameOutBase + ".out.mol";
    std::string fOutName2 = fnameOutBase + ".scsrout.mol";
    std::string fOutName3 = fnameOutBase + ".scsroutAnnotated.mol";
                            
    std::string fOutName4 = fnameOutBase + ".outNoLeavingGroups.mol";
    std::string fOutName5 = fnameOutBase + ".outNoSgroups.mol";

    RDKit::v2::FileParsers::MolFileParserParams pp;
    pp.sanitize = true;
    pp.removeHs = false;
    pp.strictParsing = true;

    MolFromMacroMolParams molFromMacroMolParams;
    molFromMacroMolParams.includeLeavingGroups = true;

    SCSRUtils::SCSRBaseHbondOptions scsrBaseHbondOptions = ScsrTest->scsrBaseHbondOptions;

    std::unique_ptr<RDKit::ROMol> mol;
    if (ScsrTest->expectedStatus == ExpectedStatus::Success) {
      REQUIRE_NOTHROW(mol = SCSRUtils::MolFromSCSRFile(fName, pp, molFromMacroMolParams,
                                            scsrBaseHbondOptions));
    } else {
      REQUIRE_THROWS(mol = SCSRUtils::MolFromSCSRFile(fName, pp, molFromMacroMolParams,
                                           scsrBaseHbondOptions));
      return;
    }

    RDKit::Chirality::removeNonExplicit3DChirality(*(mol.get()));

    CHECK(mol != nullptr);
    CHECK(mol->getNumAtoms() == ScsrTest->totalAtomCount);
    CHECK(mol->getNumBonds() == ScsrTest->totalBondCount);
    CHECK(getSubstanceGroups(*mol).size() == ScsrTest->sgroupCount);

    
    std::ifstream file(fName);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    std::stringstream buffer;
    buffer << file.rdbuf(); 
    auto molBlock = buffer.str();   
    mol = SCSRUtils::MolFromSCSRBlock(molBlock, pp, molFromMacroMolParams, scsrBaseHbondOptions);

    RDKit::Chirality::removeNonExplicit3DChirality(*(mol.get()));

    CHECK(mol != nullptr);
    CHECK(mol->getNumAtoms() == ScsrTest->totalAtomCount);
    CHECK(mol->getNumBonds() == ScsrTest->totalBondCount);
    CHECK(getSubstanceGroups(*mol).size() == ScsrTest->sgroupCount);


    MolWriterParams params;
    params.forceV3000= true;
    RDKit::MolToMolFile(*mol, fOutName, params, -1);

    auto macroMol = SCSRUtils::MacroMolFromSCSRFile(fName, pp);

    // count the hbonds which should be deleted if any

    unsigned int deletedHbondCount(0);
    if (scsrBaseHbondOptions == SCSRUtils::SCSRBaseHbondOptions::RemoveHBonds) {
      for (auto &bond : macroMol->bonds()) {
        auto macroBondInfo = bond->getMacroBondInfo();
        if (macroBondInfo) {
          if (macroBondInfo->getBond(0).bondType == Bond::BondType::HYDROGEN) {
            deletedHbondCount++;
          }
        } else if (bond->getBondType() == Bond::BondType::HYDROGEN) {
           deletedHbondCount++;
        }
      }
    }

  //   if (mol->needsUpdatePropertyCache()) {
  //     mol->updatePropertyCache(false);
  // }


    auto outMol = RDKit::MolToMacroMol(*(mol.get()), macroMol->getLocalTemplateLibrary());
    RDKit::MolToMolFile(*outMol, fOutName3, params, -1);

    auto outMacroMol = CollapseMacroMol(*outMol, macroMol->getLocalTemplateLibrary());

    CHECK(outMacroMol != nullptr);

    unsigned int macroAtomCount(ScsrTest->scsrAtomCount);
    unsigned int macroBondCount(ScsrTest->scsrBondCount);
    
    if (macroAtomCount == 0) {
      macroAtomCount = macroMol->getNumAtoms();
    }
    if (macroBondCount == 0) {
      macroBondCount = macroMol->getNumBonds();
    }
    CHECK(outMacroMol->getNumAtoms() == macroAtomCount);
    CHECK(outMacroMol->getNumBonds() + deletedHbondCount == macroBondCount);
    CHECK(outMacroMol->getLocalTemplateLibrary().getTemplates().size() == macroMol->getLocalTemplateLibrary().getTemplates().size());

    SCSRUtils::MacroMolToSCSRMolFile(*(outMacroMol.get()), fOutName2);

    std::unique_ptr<RDKit::ROMol> molReadBackIn = SCSRUtils::MolFromSCSRFile(
        fOutName2, pp, molFromMacroMolParams, scsrBaseHbondOptions);

    CHECK(molReadBackIn != nullptr);
    CHECK(molReadBackIn->getNumAtoms() == mol->getNumAtoms());
    CHECK(molReadBackIn->getNumBonds() == mol->getNumBonds());

    for (auto chiralCheck : ScsrTest->chiralChecks) {
      CHECK(mol->getAtomWithIdx(chiralCheck.first)->getChiralTag() ==
            chiralCheck.second);
    }

    // now make the expanded mol in "query" mode - not including any leaving
    // groups

    molFromMacroMolParams.includeLeavingGroups = false;
    std::unique_ptr<RDKit::RWMol> molNoLeavingGroups;
    if (ScsrTest->expectedStatus == ExpectedStatus::Success) {
      molNoLeavingGroups = SCSRUtils::MolFromSCSRFile(fName, pp, molFromMacroMolParams,
                                           scsrBaseHbondOptions);
    } else {
      molNoLeavingGroups = SCSRUtils::MolFromSCSRFile(fName, pp, molFromMacroMolParams,
                                           scsrBaseHbondOptions);
      return;
    }

    CHECK(molNoLeavingGroups != nullptr);
    CHECK(molNoLeavingGroups->getNumAtoms() == ScsrTest->totalQueryAtomCount);
    CHECK(molNoLeavingGroups->getNumBonds() == ScsrTest->totalQueryBondCount);
    CHECK(getSubstanceGroups(*molNoLeavingGroups).size() ==
          ScsrTest->querySgroupCount);

    for (auto chiralCheck : ScsrTest->chiralChecksQuery) {
      CHECK(mol->getAtomWithIdx(chiralCheck.first)->getChiralTag() ==
            chiralCheck.second);
    }

    RDKit::MolToMolFile(*molNoLeavingGroups, fOutName4, params, -1);

    /* test atomistic output WITHOUT sgroups */
    molFromMacroMolParams.outputSgroups = false;
    molFromMacroMolParams.includeLeavingGroups = true;

    if (ScsrTest->expectedStatus == ExpectedStatus::Success) {
      REQUIRE_NOTHROW(mol = SCSRUtils::MolFromSCSRFile(fName, pp, molFromMacroMolParams,
                                            scsrBaseHbondOptions));
    } else {
      REQUIRE_THROWS(mol = SCSRUtils::MolFromSCSRFile(fName, pp, molFromMacroMolParams,
                                           scsrBaseHbondOptions));
      return;
    }

    RDKit::Chirality::removeNonExplicit3DChirality(*(mol.get()));

    CHECK(mol != nullptr);
    CHECK(mol->getNumAtoms() == ScsrTest->totalAtomCount);
    CHECK(mol->getNumBonds() == ScsrTest->totalBondCount);
    CHECK(getSubstanceGroups(*mol).size() == 0);
  RDKit::MolToMolFile(*mol, fOutName5, params, -1);
  }

  void threeLetterCodeTest(const ScsrTest *ScsrTest) {
    BOOST_LOG(rdInfoLog) << "testing scsr  files with three letter codes"
                         << std::endl;

    std::string rdbase = getenv("RDBASE");
    std::string fName = rdbase +
                        "/Code/GraphMol/FileParsers/test_data/macromols/" +
                        ScsrTest->fileName;

    RDKit::v2::FileParsers::MolFileParserParams pp;
    pp.sanitize = false;
    pp.removeHs = false;
    pp.strictParsing = false;

    MolFromMacroMolParams molFromMacroMolParams;
    molFromMacroMolParams.includeLeavingGroups = true;
    SCSRUtils::SCSRBaseHbondOptions scsrBaseHbondOptions = ScsrTest->scsrBaseHbondOptions;

    std::unique_ptr<RDKit::RWMol> mol;
    if (ScsrTest->expectedStatus == ExpectedStatus::Success) {
      REQUIRE_NOTHROW(mol = SCSRUtils::MolFromSCSRFile(fName, pp, molFromMacroMolParams,
                                            scsrBaseHbondOptions));
    } else {
      REQUIRE_THROWS(mol = SCSRUtils::MolFromSCSRFile(fName, pp, molFromMacroMolParams,
                                           scsrBaseHbondOptions));
      return;
    }

    RDKit::Chirality::removeNonExplicit3DChirality(*(mol.get()));

    CHECK(mol);
    CHECK(mol->getNumAtoms() == ScsrTest->totalAtomCount);
    CHECK(mol->getNumBonds() == ScsrTest->totalBondCount);
    CHECK(getSubstanceGroups(*mol).size() == ScsrTest->sgroupCount);

    molFromMacroMolParams.includeLeavingGroups = true;
    std::unique_ptr<RWMol> mol2;
    if (ScsrTest->expectedStatus == ExpectedStatus::Success) {
      REQUIRE_NOTHROW(mol2 = SCSRUtils::MolFromSCSRFile(fName, pp, molFromMacroMolParams,
                                             scsrBaseHbondOptions));
    } else {
      REQUIRE_THROWS(mol2 = SCSRUtils::MolFromSCSRFile(fName, pp, molFromMacroMolParams,
                                            scsrBaseHbondOptions));
      return;
    }

    // const std::unique_ptr<RDKit::RWMol> mol;
    CHECK(mol2);
    CHECK(mol2->getNumAtoms() == ScsrTest->totalQueryAtomCount);
    CHECK(mol2->getNumBonds() == ScsrTest->totalQueryBondCount);
    CHECK(getSubstanceGroups(*mol2).size() == ScsrTest->querySgroupCount);

    auto parserParams = RDKit::v2::FileParsers::MolFileParserParams();
    auto writerParams = RDKit::MolWriterParams();
    auto macroMol = SCSRUtils::MacroMolFromSCSRFile(fName, parserParams);
    auto scsrBlock = SCSRUtils::MacroMolToSCSRMolBlock(*macroMol.get(), writerParams, -1, true);
    macroMol = SCSRUtils::MacroMolFromSCSRBlock(scsrBlock,parserParams);

    std::string templateName;
    macroMol->getAtomWithIdx(1)->getProp(common_properties::dummyLabel ,templateName);
    CHECK(templateName  == "Cys");

  };

};

TEST_CASE("scsrTests", "scsrTests") {
  SECTION("basics") {
    ScsrMolTest ScsrMolTest;
    for (auto ScsrTest : ScsrMolTest.ScsrTests) {
      BOOST_LOG(rdInfoLog) << "Test: " << ScsrTest.fileName << std::endl;

      ScsrMolTest.testScsrFiles(&ScsrTest);
    }
  }
}

TEST_CASE("nestedParens", "nestedParens") {
  SECTION("basics") {
    BOOST_LOG(rdInfoLog) << "testing names with parens" << std::endl;
    std::string filename = "Ugly.mol";
    INFO(filename);

    std::string rdbase = getenv("RDBASE");
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/macromols/" + filename;

    RDKit::v2::FileParsers::MolFileParserParams pp;
    pp.sanitize = true;
    pp.removeHs = false;
    pp.strictParsing = true;

    MolFromMacroMolParams molFromMacroMolParams;
    molFromMacroMolParams.includeLeavingGroups = true;
    SCSRUtils::SCSRBaseHbondOptions scsrBaseHbondOptions = SCSRUtils::SCSRBaseHbondOptions::FullHBonds;

    std::unique_ptr<RDKit::RWMol> mol;

    REQUIRE_NOTHROW(mol = SCSRUtils::MolFromSCSRFile(fName, pp, molFromMacroMolParams,
                                          scsrBaseHbondOptions));

    CHECK(mol != nullptr);
    CHECK(mol->getNumAtoms() == 140);
    CHECK(mol->getNumBonds() == 151);
    CHECK(getSubstanceGroups(*mol).size() == 11);

    // check that the macro atoms parsed have the correct names.  In the
    // output mol these names will not appear in the SGROUP that defined the
    // macroatom. we willcheck just a couple of them

    std::string sgroupName;
    getSubstanceGroups(*mol)[1].getProp("LABEL", sgroupName);
    std::string expected = "((cPr)O(2S-Me)Et)NGly";
    CHECK(sgroupName.substr(0, expected.length()) == expected);

    getSubstanceGroups(*mol)[5].getProp("LABEL", sgroupName);
    expected = "Phe(b-Me2)";
    CHECK(sgroupName.substr(0, expected.length()) == expected);

    getSubstanceGroups(*mol)[7].getProp("LABEL", sgroupName);
    expected = "(PhO(2S-Me)Et)NGly";
    CHECK(sgroupName.substr(0, expected.length()) == expected);
  }
}

TEST_CASE("threeLetterCodeTest", "threeLetterCodeTest") {
  SECTION("basics") {
    std::list<ScsrMolTest::ScsrTest> scsrTests{
        ScsrMolTest::ScsrTest(
            "PepTla.mol", ScsrMolTest::ExpectedStatus::Success,
            SCSRUtils::SCSRBaseHbondOptions::FullHBonds, 26, 25, 4, 26, 25, 4),

    };
    ScsrMolTest scsrMolTest;

    for (auto scsrTest : scsrTests) {
      BOOST_LOG(rdInfoLog) << "Test: " << scsrTest.fileName << std::endl;

      scsrMolTest.threeLetterCodeTest(&scsrTest);
    }
  }
}


TEST_CASE("MissingTemplate") {
  SECTION("basic") {

    std::string rdbase = getenv("RDBASE");
    std::string fName = rdbase +
                        "/Code/GraphMol/FileParsers/test_data/macromols/FromBioviaDoc.mol";
    std::string templateFName = rdbase +
                        "/Code/GraphMol/FileParsers/test_data/macromols/TwoHelmTemplates.mol";
    
    RDKit::v2::FileParsers::MolFileParserParams pp;
    pp.sanitize = true;
    pp.removeHs = false;
    pp.strictParsing = true;

    std::unique_ptr<RDKit::RWMol> mol;
    std::unique_ptr<RDKit::MacroMol> templates;
    RDKit::v2::FileParsers::MolFileParserParams params;
    RDKit::MolFromMacroMolParams molFromMacroMolParams;
    //mol = SCSRUtils::MolFromSCSRBlock(molOut, params, molFromMacroMolParams);
    molFromMacroMolParams.includeLeavingGroups = true;

    REQUIRE_NOTHROW(mol = SCSRUtils::MolFromSCSRFile(fName, params,molFromMacroMolParams,  SCSRUtils::SCSRBaseHbondOptions::FullHBonds));

    REQUIRE_NOTHROW(templates = SCSRUtils::MacroMolFromSCSRFile(templateFName, params));


    auto groupedMol = MolToMacroMol( *mol.get(),  templates->getLocalTemplateLibrary());

    auto macroMol = CollapseMacroMol(*groupedMol.get(), templates->getLocalTemplateLibrary());

    auto molOut = SCSRUtils::MacroMolToSCSRMolBlock(*macroMol);
    mol = SCSRUtils::MolFromSCSRBlock(molOut, params, molFromMacroMolParams);
     
    CHECK(mol->getNumAtoms() == 27);
    CHECK(mol->getNumBonds() == 26);
    
  }

}
