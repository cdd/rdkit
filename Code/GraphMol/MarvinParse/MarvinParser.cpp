//
//  Copyright (C) 2002-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

// this file parses MRV file for molecules and reactions


#include <string>
#include <exception>
#include <iostream>
#include <fstream>
#include <vector>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSGroupParsing.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include "MarvinParser.h"
#include "MarvinDefs.h"
#include <GraphMol/Conformer.h>

#include <GraphMol/RDKitQueries.h>
#include <GraphMol/StereoGroup.h>
#include <GraphMol/SubstanceGroup.h>

#include "MarvinParser.h"

#include <RDGeneral/StreamOps.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/LocaleSwitcher.h>




#ifdef RDKIT_USE_BOOST_REGEX
#include <boost/regex.hpp>
using boost::regex;
using boost::regex_match;
using boost::smatch;
//using boost::algorithm>;

#else
#include <regex>
using std::regex;
using std::regex_match;
using std::smatch;
#endif


using namespace RDKit::SGroupParsing;

namespace RDKit
{

  /*
	Imports the Marvin-specific dialect of CML (Chemical Markup Language) and converts it to datastructures
	that are compatible with Molfile, RXNfile, and Molfile complemented with canvas objects.
  */
  class MarvinCMLReader
  {
    
    public:
      RWMol *mol;
      ChemicalReaction *rxn;

   public:
    MarvinCMLReader()
      : mol(NULL)
      , rxn(NULL)
    {
    };

    ~MarvinCMLReader()
    { 
    };




    public:
  
    //this routine does the work of parsing.  It returns either an RWMol * or a ChemicalStructure *
    // either way it is cast to a void *

    void *parse(std::istream &is, bool &isReaction, bool sanitize=false, bool removeHs=false)
    {
      // Create empty property tree object
      using boost::property_tree::ptree;
      ptree tree;

      // Parse the XML into the property tree.

      read_xml(is, tree);

      // find a molecule or a reaction

      bool molFlag = false;
      boost::property_tree::ptree molOrRxn;
      try
      {
        molOrRxn = tree.get_child("cml.MDocument.MChemicalStruct.molecule");
        molFlag = true;
      }
      catch(const std::exception& e)
      {
      // do nothing 

      }
    
      if (molFlag)
      {
        mol = (RWMol *) parseMolecule(molOrRxn, sanitize, removeHs);
        isReaction = false;
        return (void *)mol;
      }
      
      try
      {
        molOrRxn = tree.get_child("cml.MDocument.MChemicalStruct.reaction");
      }
      catch(const std::exception& e)
      {
        // check for sparse marvin data - no reaction nor molecule parts
        // just return an empty RWMol *

        try
        {
          molOrRxn = tree.get_child("cml.MDocument.MChemicalStruct");
          return new RWMol();
        }
        catch(const std::exception& e)
        {
          try
          {
          molOrRxn = tree.get_child("cml.MDocument");
          return new RWMol();
          }
          catch(const std::exception& e)
          {
          try
          {
            molOrRxn = tree.get_child("cml");
            return new RWMol();
          }
          catch(const std::exception& e)
          {
            throw FileParseException("Expected \"molecule\" or \"reaction\" in MRV file");
          }
          }
        }
      }

      rxn = parseReaction(molOrRxn, tree.get_child("cml.MDocument"), sanitize, removeHs); 
      isReaction = true;
      return (void *)rxn;
    }

  
    bool getCleanDouble(std::string strToParse, double &outDouble)
    {
    
      if (boost::algorithm::trim_copy(strToParse) != strToParse)   // should be no white space
      return false;
      try
      {
        outDouble = boost::lexical_cast<double>(strToParse);
      }
      catch(const std::exception& e)
      {
        return false;
      }


      return true;

    }
    
    bool getCleanInt(std::string strToParse, int &outInt)
    {
    
      if (boost::algorithm::trim_copy(strToParse) != strToParse)  // should be no white space
      return false;
      try
      {
        outInt = boost::lexical_cast<int>(strToParse);
      }
      catch(const std::exception& e)
      {
        return false;
      }

      return true;
    }
    
    

    Atom *molAtomFromMarvinAtom(MarvinAtom *marvinAtom)
    {
      Atom *res = NULL;

      try
      {
        std::string symb = marvinAtom->elementType;
        boost::trim(symb);
        
        if (symb == "*") 
        {
          auto *query = new QueryAtom(0);   // wHAT IS AUTO *
          res = query;
          query->setQuery(makeAtomNullQuery());
          query->setProp(common_properties::dummyLabel, "*");
          // queries have no implicit Hs:
          res->setNoImplicit(true);
        } 
        else if (symb  == "R") 
        {
          // It must have an rgroupRef/>
          
          if ( marvinAtom->rgroupRef < 0)
            throw FileParseException("Expected an R group atom to have an rgroupRef designation");
          
          auto *query = new QueryAtom(0);   // wHAT IS AUTO *
          res = query;
      
          query->setProp(common_properties::_MolFileRLabel, marvinAtom->rgroupRef);
          std::string dLabel = "R" + std::to_string(marvinAtom->rgroupRef);
          query->setProp(common_properties::dummyLabel, dLabel);
          if (marvinAtom->mrvAlias != "")
            query->setProp(common_properties::molFileAlias,marvinAtom->mrvAlias);
          query->setIsotope(marvinAtom->rgroupRef);
          query->setQuery(makeAtomNullQuery());
        }

        else if (symb.size() <= 2) 
        {
          // must be a regular atom


          if (symb.size() == 2 &&  symb[1] >= 'A' && symb[1] <= 'Z')
            symb[1] = static_cast<char>(tolower(symb[1]));
          res = new Atom();
        
          try 
          {
            res->setAtomicNum(PeriodicTable::getTable()->getAtomicNumber(symb));
          } 
          catch (const Invar::Invariant &e) 
          {
            delete res;
            throw FileParseException(e.what());
          }
        }
        else
          throw FileParseException("Unrecognized atom type in MRV file");

      
        //res->setPos(marvinAtom.x2,marvinAtom->y2,0.0);
        
        if (marvinAtom->formalCharge != 0) 
          res->setFormalCharge(marvinAtom->formalCharge);

        
        if (marvinAtom->isotope != 0)
        {
          res->setIsotope(marvinAtom->isotope);
          res->setProp(common_properties::_hasMassQuery, true);
        }

        if (marvinAtom->radical != "")
          res->setNumRadicalElectrons(marvinRadicalToRadicalElectrons.at(marvinAtom->radical));

        //  we no not get parity (chirality) frommarvin file
        // res->setProp(common_properties::molParity, parity);
        // res->setProp(common_properties::molStereoCare, stereoCare);

        // total valence is not parsed from marvin file
        // res->setProp(common_properties::molTotValence, totValence);

        // rxnRole is not parsed from marvin file
        // res->setProp(common_properties::molRxnRole, rxnRole);
        // res->setProp(common_properties::molRxnComponent, rxnComponent);

        if (marvinAtom->mrvMap != 0)
          res->setProp(common_properties::molAtomMapNumber, marvinAtom->mrvMap);

      // not set
      //res->setProp(common_properties::molInversionFlag, inversionFlag);
      // res->setProp("molExactChangeFlag", exactChangeFlag);
        
        return res;
      }
      catch(const std::exception& e)
      {
        delete res;
        throw;
      }
    }

    void molBondFromMarvinBond(MarvinBond *marvinBond, MarvinMol *marvinMol, RWMol *mol, bool &chiralityPossible)
    {  
      unsigned int bType, stereo;
      Bond *bond = nullptr;

      try
      {
        int idx1 = marvinMol->getAtomIndex(marvinBond->atomRefs2[0]);
        int idx2 = marvinMol->getAtomIndex(marvinBond->atomRefs2[1]);

        Bond::BondType type;

        // if there is a query type, that takes precidence over the bond type
        std::string temp = boost::algorithm::to_upper_copy(marvinBond->queryType);
        if (temp != "")
        {
          if (temp == "SD")
          {
            bType = 5;
            type = Bond::UNSPECIFIED;
            bond = new QueryBond;
            bond->setQuery(makeSingleOrDoubleBondQuery());
          }
          else if (temp == "SA")
          {
            bType = 6;
            type = Bond::UNSPECIFIED;
            bond = new QueryBond;
            bond->setQuery(makeSingleOrAromaticBondQuery());
          }
          else if (temp == "DA")
          {
            bType = 7;
            type = Bond::UNSPECIFIED;
            bond = new QueryBond;
            bond->setQuery(makeDoubleOrAromaticBondQuery());
          }
          else
          {
            std::ostringstream err;
            err      << "unrecognized query bond type " << marvinBond->queryType << " in MRV File ";
            throw FileParseException(err.str());
          } 
        }
        else   // if no query type, the bond order takes over
        {
          temp = boost::algorithm::to_upper_copy(marvinBond->order);
          if (temp == "1")
          {
            type = Bond::SINGLE;
            bond = new Bond;
            bType = 1;
          }
          else if (temp == "2")
          {
            type = Bond::DOUBLE;
            bond = new Bond;
            bType = 2;
          }
          else if (temp == "3")
          {
            type = Bond::TRIPLE;
            bond = new Bond;
            bType = 3;
          }
          else if (temp == "A")
          {
            bType = 4;
            type = Bond::AROMATIC;
            bond = new Bond;
          }
          else
          {
            std::ostringstream err;
            err      << "unrecognized bond type " << marvinBond->order << " in MRV File ";
            throw FileParseException(err.str());
          } 
        }      


        bond->setBeginAtomIdx(idx1);
        bond->setEndAtomIdx(idx2);
        bond->setBondType(type);
        bond->setProp(common_properties::_MolFileBondType, bType);
        
        temp = boost::algorithm::to_lower_copy(marvinBond->bondStereo);
        if (temp == "")  
        {         
          bond->setBondDir(Bond::NONE);
          stereo = 0;
        }
        else if (temp == "w")  
        {         
          bond->setBondDir(Bond::BEGINWEDGE);
          chiralityPossible = true;
          stereo = 1;
        } 
        else if (temp == "h")  
        {         
          stereo = 6;
          chiralityPossible = true;
          bond->setBondDir(Bond::BEGINDASH);
        }
        else
        {
          std::ostringstream err;
          err      << "unrecognized bond stereo " << marvinBond->bondStereo << " in MRV File ";
          throw FileParseException(err.str());
        }
        
        bond->setProp(common_properties::_MolFileBondStereo, stereo);
        
        // if we got an aromatic bond set the flag on the bond and the connected
        // atoms
        if (bond->getBondType() == Bond::AROMATIC) 
          bond->setIsAromatic(true);

        // v2k has no way to set stereoCare on bonds, so set the property if both
        // the beginning and end atoms have it set:
        int care1 = 0;
        int care2 = 0;
        if (!bond->hasProp(common_properties::molStereoCare) &&
          mol->getAtomWithIdx(bond->getBeginAtomIdx())
            ->getPropIfPresent(common_properties::molStereoCare, care1) &&
          mol->getAtomWithIdx(bond->getEndAtomIdx())
            ->getPropIfPresent(common_properties::molStereoCare, care2)) 
        {
          if (care1 && care2) 
          {
            bond->setProp(common_properties::molStereoCare, 1);
          }
        }
        mol->addBond(bond, true);
      }
      catch(const std::exception& e)
      {
        delete bond;
        throw;
      }
    }   
  
    RWMol *parseMolecule(MarvinMol *marvinMol, bool sanitize=false, bool removeHs=false)
    {
      
      std::vector<MarvinStereoGroup *> stereoGroups;
      Conformer *conf = NULL;
      RWMol *mol = NULL;

      try
      {
        mol = new RWMol();

        //mol->setProp(common_properties::_Name, marvinMol->molID);
        mol->setProp("_MolFileComments", "Generated by RDKit");

        // set the atoms
        
        unsigned int i;
        std::vector<MarvinAtom *>::const_iterator atomIter;

        conf = new Conformer(marvinMol->atoms.size());
        conf->set3D(false);

        for (i = 0, atomIter = marvinMol->atoms.begin() ; atomIter != marvinMol->atoms.end();  ++i, ++atomIter)
        {		  
          Atom *atom = molAtomFromMarvinAtom(*atomIter);
          unsigned int aid = mol->addAtom(atom, false, true);

          RDGeom::Point3D pos;
          pos.x = (*atomIter)->x2;
          pos.y = (*atomIter)->y2;
          pos.z = 0.0;
          
          conf->setAtomPos(aid, pos);
          //mol->setAtomBookmark(atom, i+1);


          //also collect the stereo groups here

          if ((*atomIter)->mrvStereoGroup != "")
          {
            RDKit::StereoGroupType groupType;
            int groupNumber;

            // get the group parts
            std::string temp = boost::algorithm::to_lower_copy((*atomIter)->mrvStereoGroup);
            if (boost::starts_with(temp, "abs"))
            {
              groupType  = RDKit::StereoGroupType::STEREO_ABSOLUTE;
              if (! getCleanInt(temp.substr(3), groupNumber))
              throw FileParseException("Group Number must be an integer in a stereo group in a MRV file"); 
            }
            else if (boost::starts_with(temp, "and"))
            {
              groupType  = RDKit::StereoGroupType::STEREO_AND;
              if (! getCleanInt(temp.substr(3),groupNumber))
              throw FileParseException("Group Number must be an integer in a stereo group in a MRV file"); 
            }
            else if (boost::starts_with(temp, "or"))
            {
              groupType  = RDKit::StereoGroupType::STEREO_OR;
              if (! getCleanInt(temp.substr(2), groupNumber))
              throw FileParseException("Group Number must be an integer in a stereo group in a MRV file"); 
            }
            else
              throw FileParseException("Unrecognized group definition"); 


            // see if the group already exists

            MarvinStereoGroup *marvinStereoGroup;
            auto groupIter = find_if(stereoGroups.begin(), stereoGroups.end(), [groupType, groupNumber](const MarvinStereoGroup *arg) { 
                            return arg->groupNumber ==groupNumber && arg->groupType == groupType; });
            if (groupIter != stereoGroups.end())
              marvinStereoGroup = *groupIter;
            else 
            {
              marvinStereoGroup = new MarvinStereoGroup(groupType, groupNumber);
              stereoGroups.push_back(marvinStereoGroup);
            }

            // add this atom to the group

            marvinStereoGroup->atoms.push_back((unsigned int)marvinMol->getAtomIndex((*atomIter)->id));
          }
        }

        mol->addConformer(conf, true);
        conf = NULL;  // conf now owned by mol
        


        // set the bonds

        std::vector<MarvinBond *>::const_iterator bondIter;
        bool chiralityPossible = false;

        for (i = 0, bondIter = marvinMol->bonds.begin() ; bondIter != marvinMol->bonds.end();  ++i, ++bondIter)
          molBondFromMarvinBond(*bondIter, marvinMol, mol, chiralityPossible);
          
        //  add the stereo groups

        std::vector<StereoGroup> groups;

        for (std::vector<MarvinStereoGroup *>::const_iterator groupIter = stereoGroups.begin() ; groupIter != stereoGroups.end() ; ++groupIter)
        {
          std::vector<Atom *> atoms;
          for (std::vector<unsigned int>::const_iterator it = (*groupIter)->atoms.begin() ;  it != (*groupIter)->atoms.end();  ++it)
          atoms.push_back(mol->getAtomWithIdx(*it));

          groups.emplace_back((*groupIter)->groupType, std::move(atoms));
          if (!groups.empty()) 
          mol->setStereoGroups(std::move(groups));
        }
        // add the super atoms records

        int sequenceId;
        std::vector<MarvinSuperInfo *>::iterator superIter;
        for ( sequenceId= 1, superIter = marvinMol->superInfos.begin() ; superIter != marvinMol->superInfos.end() ; ++superIter, ++sequenceId)
        {
          std::string typ = "SUP";         
          SubstanceGroup sgroup = SubstanceGroup(mol, typ);
          sgroup.setProp<unsigned int>("index", sequenceId);

          void (SubstanceGroup::*sGroupAddIndexedElement)(const unsigned int) = nullptr;
          sGroupAddIndexedElement = &SubstanceGroup::addAtomWithIdx;


          std::vector<std::string>::const_iterator atomIter;
          for (atomIter = (*superIter)->atoms.begin(); atomIter != (*superIter)->atoms.end(); ++atomIter)
          (sgroup.*sGroupAddIndexedElement)(marvinMol->getAtomIndex(*atomIter));

          sgroup.setProp("LABEL", (*superIter)->title);    

          if (sgroup.getIsValid())
          addSubstanceGroup(*mol, sgroup);
        }
      
        // now the SruGroups

        // note: sequence continues counting from the loop above
        for ( std::vector<MarvinSruSgroup *>::iterator sruIter = marvinMol->sruSgroups.begin() ; sruIter != marvinMol->sruSgroups.end() ; ++sruIter, ++sequenceId)
        {
          std::string typ = "SRU";         
          SubstanceGroup sgroup = SubstanceGroup(mol, typ);
          sgroup.setProp<unsigned int>("index", sequenceId);

          sgroup.setProp("CONNECT", (*sruIter)->connect);

          void (SubstanceGroup::*sGroupAddIndexedElement)(const unsigned int) = nullptr;
          sGroupAddIndexedElement = &SubstanceGroup::addAtomWithIdx;

          void (SubstanceGroup::*sGroupAddIndexedElementBond)(const unsigned int) = nullptr;
          sGroupAddIndexedElementBond = &SubstanceGroup::addBondWithIdx;

          std::vector<MarvinAtom *>::const_iterator atomIter;
          for (atomIter = (*sruIter)->atoms.begin(); atomIter != (*sruIter)->atoms.end(); ++atomIter)
          {
            int atomIndex = marvinMol->getAtomIndex((*atomIter)->id);
            (sgroup.*sGroupAddIndexedElement)(atomIndex);
          }
          std::vector<MarvinBond *>::const_iterator bondIter;
          for (bondIter = (*sruIter)->bonds.begin(); bondIter != (*sruIter)->bonds.end(); ++bondIter)
          {
            //int bondIndex = marvinMol->getBondIndex((*bondIter)) + 1;
            int bondIndex = marvinMol->getBondIndex((*bondIter)->id);
            (sgroup.*sGroupAddIndexedElementBond)(bondIndex);
          }
          sgroup.setProp("LABEL", (*sruIter)->title);

          if (sgroup.getIsValid())
            addSubstanceGroup(*mol, sgroup);
        }


        mol->clearAllAtomBookmarks();
        mol->clearAllBondBookmarks();

        // calculate explicit valence on each atom:
        for (RWMol::AtomIterator atomIt = mol->beginAtoms(); atomIt != mol->endAtoms(); ++atomIt) 
        {
          (*atomIt)->calcExplicitValence(false);
          (*atomIt)->calcImplicitValence(false);
        }


        // update the chirality and stereo-chemistry
        //
        // NOTE: we detect the stereochemistry before sanitizing/removing
        // hydrogens because the removal of H atoms may actually remove
        // the wedged bond from the molecule.  This wipes out the only
        // sign that chirality ever existed and makes us sad... so first
        // perceive chirality, then remove the Hs and sanitize.
        //

        const Conformer &conf2 = mol->getConformer();
        if (chiralityPossible) 
          DetectAtomStereoChemistry(*mol, &conf2);

        if (sanitize) 
        {
          if (removeHs)
            MolOps::removeHs(*mol, false, false);
          else 
            MolOps::sanitizeMol(*mol);
      
          // now that atom stereochem has been perceived, the wedging
          // information is no longer needed, so we clear
          // single bond dir flags:
          ClearSingleBondDirFlags(*mol);

          // unlike DetectAtomStereoChemistry we call detectBondStereochemistry
          // here after sanitization because we need the ring information:
          MolOps::detectBondStereochemistry(*mol);
          
          MolOps::assignStereochemistry(*mol, true, true, true);
        } 
        else 
        {
          // we still need to do something about double bond stereochemistry
          // (was github issue 337)
          // now that atom stereochem has been perceived, the wedging
          // information is no longer needed, so we clear
          // single bond dir flags:
          ClearSingleBondDirFlags(*mol);
          MolOps::detectBondStereochemistry(*mol);
        }

        if (mol->hasProp(common_properties::_NeedsQueryScan)) {
          mol->clearProp(common_properties::_NeedsQueryScan);
          QueryOps::completeMolQueries(mol);
        }

        // clean up

        for (std::vector<MarvinStereoGroup *>::iterator it = stereoGroups.begin() ; it != stereoGroups.end(); ++it)
          delete *it;
        
        return mol;
      }

      catch(const std::exception& e)
      {
        delete mol;

        delete conf;

        for (std::vector<MarvinStereoGroup *>::iterator it = stereoGroups.begin() ; it != stereoGroups.end(); ++it)
            delete *it;
        
        throw;
      }
    }

    

    RWMol *parseMolecule(boost::property_tree::ptree molTree, bool sanitize=false, bool removeHs=false)
    {
      MarvinMol *marvinMol = NULL;

      try
      {
        marvinMol = (MarvinMol *)parseMarvinMolecule(molTree);
        marvinMol->convertFromSuperAtoms();

        RWMol* mol = parseMolecule(marvinMol , sanitize, removeHs);

        delete marvinMol;

        return mol;
      }
      catch(const std::exception& e)
      {
        delete mol;

        delete marvinMol;
        
        throw;
      }
    }

    MarvinMolBase *parseMarvinMolecule(boost::property_tree::ptree molTree, MarvinMol *parentMol=NULL)  // parent is for sub-mols
    {
      MarvinMolBase *res=NULL;

      try
      {
        std::string molID = molTree.get<std::string>("<xmlattr>.molID", "");

        if (molID == "")
          throw FileParseException("Expected a molID in MRV file");
        
        std::string role = "";

        if (parentMol == NULL)\
        {
          res = new MarvinMol();
        } 
        else // is is a sub-mol - used for super atoms and sruGroups
        {
          role = molTree.get<std::string>("<xmlattr>.role", "");
          if (role == "")
          throw FileParseException("Expected a role for a sub-molecule in MRV file");

          if (role == "SuperatomSgroup")
          {
            MarvinSuperatomSgroup *marvinSuperatomSgroup = new MarvinSuperatomSgroup();
            res = marvinSuperatomSgroup;  

            marvinSuperatomSgroup->id = molTree.get<std::string>("<xmlattr>.id", "");

            marvinSuperatomSgroup->title = molTree.get<std::string>("<xmlattr>.title", "");
            if (marvinSuperatomSgroup->title == "")
            throw FileParseException("Expected  title for a SuperatomSgroup definition in MRV file");
          
          } 
          else if (role == "SruSgroup")
          {
            //      <molecule molID="m2" id="sg1" role="SruSgroup" atomRefs="a1 a3 a4" title="n" connect="ht" correspondence="" bondList=""/>
            //      <molecule molID="m3" id="sg2" role="SruSgroup" atomRefs="a5 a6 a7 a8" title="n" connect="hh" correspondence="" bondList=""/>
            //      <molecule molID="m4" id="sg3" role="SruSgroup" atomRefs="a10 a11" title="n" connect="eu" correspondence="" bondList=""/></molecule>

            MarvinSruSgroup *marvinSruSgroup = new MarvinSruSgroup();
            res = marvinSruSgroup;
            
            marvinSruSgroup->id = molTree.get<std::string>("<xmlattr>.id", "");

            std::string atomRefsStr = molTree.get<std::string>("<xmlattr>.atomRefs", "");
            if (atomRefsStr == "")
              throw FileParseException("Expected  atomRefs for a SruSgroup definition in MRV file");


            std::vector<std::string> atomList;
            boost::algorithm::split(atomList, atomRefsStr, boost::algorithm::is_space());
            for (std::vector<std::string>::const_iterator it = atomList.begin() ;  it != atomList.end(); ++it)
            {
              auto atomIter = find_if(parentMol->atoms.begin(), parentMol->atoms.end(), [it](const MarvinAtom *arg) { 
                            return arg->id == *it; });
              if (atomIter == parentMol->atoms.end())
              throw FileParseException("AtomRef specification for an SRU group definition was not found in the atom array in MRV file");
              marvinSruSgroup->atoms.push_back(*atomIter);
            }


            // boost::algorithm::split(marvinSruSgroup->atomRefs, atomRefsStr, boost::algorithm::is_space());
            // for (std::vector<std::string>::const_iterator it = marvinSruSgroup->atomRefs.begin() ;  it != marvinSruSgroup->atomRefs.end(); ++it)
            //   if (!boost::algorithm::contains(parentMol->atoms, std::vector<std::string>{*it}, MarvinMol::atomRefInAtoms ))
            //     throw FileParseException("All of the AtomRefs in the sruSgroup must in the parent molecule");

            marvinSruSgroup->title = molTree.get<std::string>("<xmlattr>.title", "");
            if (marvinSruSgroup->title == "")
              throw FileParseException("Expected  title for a SruSgroup definition in MRV file");

            marvinSruSgroup->connect = molTree.get<std::string>("<xmlattr>.connect", "");
            if (!boost::algorithm::contains(sruSgroupConnectChoices,std::vector<std::string>{marvinSruSgroup->connect}))
            {
              std::ostringstream err;
              err << "Expected  a connect  string of \"" << boost::algorithm::join(sruSgroupConnectChoices, ", ") << "\" for a SruSgroup definition in MRV file";
              throw FileParseException(err.str());
            }

            marvinSruSgroup->correspondence = molTree.get<std::string>("<xmlattr>.correspondence", "");

            std::string bondListStr = molTree.get<std::string>("<xmlattr>.bondList", "");
            if (bondListStr != "")
            {
              std::vector<std::string> bondList;
              boost::algorithm::split(bondList, bondListStr, boost::algorithm::is_space());
              for (std::vector<std::string>::const_iterator it = bondList.begin() ;  it != bondList.end(); ++it)
              {
              auto bondIter = find_if(parentMol->bonds.begin(), parentMol->bonds.end(), [it](const MarvinBond *arg) { 
                            return arg->id == *it; });
              if (bondIter == parentMol->bonds.end())
                throw FileParseException("Bond specification for an SRU group definition was not found in the bond array in MRV file");
              marvinSruSgroup->bonds.push_back(*bondIter);
              }

              // boost::algorithm::split(marvinSruSgroup->bondList, bondListStr, boost::algorithm::is_space());
              // for (std::vector<std::string>::const_iterator it = marvinSruSgroup->bondList.begin() ;  it != marvinSruSgroup->bondList.end(); ++it)
              //   if (!boost::algorithm::contains(parentMol->bonds, std::vector<std::string>{*it}, MarvinMol::bondRefInBonds ))
              //     throw FileParseException("All of the BondList entries in the sruSgroup must in the parent molecule");
            }
          }
          else if (role == "MultipleSgroup")
          {
            throw FileParseException("MultipleSgroup in not yet implemented in MRV file");
          }
          else
            throw FileParseException("Unexpected role " + role + " in MRV file");
        }
        
        res->molID = molID;

        // get atoms if this is NOT a role == "SruSgroup"

        if (role != "SruSgroup")
        {
          boost::property_tree::ptree atomArray = molTree.get_child("atomArray");

          // there are two types of atom arrays:
          // <atomArray atomID="a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11" elementType="C C C C C C Cl C N O O" formalCharge="0 0 0 0 0 0 0 0 1 0 -1" lonePair="0 0 0 0 0 0 3 0 0 2 3" x2="-4.3334 -5.6670 -5.6670 -4.3334 -2.9997 -2.9997 -4.3335 -1.6660 -7.0007 -1.6660 -0.3323" y2="1.7693 0.9993 -0.5409 -1.3109 -0.5409 0.9993 3.3093 -1.3109 -1.3109 -2.8509 -0.5410"></atomArray>
          //  AND
          // <atomArray>
          //       <atom id="a1" elementType="C" x2="-9.4583" y2="1.9358" mrvStereoGroup="and1"/>
          //       <atom id="a2" elementType="C" x2="-10.7921" y2="1.1658"/>
          //       <atom id="a3" elementType="C" x2="-10.7921" y2="-0.3744"/>
          //       <atom id="a8" elementType="O" x2="-12.1257" y2="-1.1444" lonePair="2"/>
          //   </atomArray>
          
          // See which one we have

          
          std::string atomID = atomArray.get<std::string>("<xmlattr>.atomID", "");
          if (atomID == "")
          {
            // long form - each atom on a line

            BOOST_FOREACH( boost::property_tree::ptree::value_type &v, molTree.get_child("atomArray"))
            {
              MarvinAtom *mrvAtom = new MarvinAtom();   
              res->atoms.push_back(mrvAtom);

              mrvAtom->id = v.second.get<std::string>("<xmlattr>.id", "");
              mrvAtom->elementType = v.second.get<std::string>("<xmlattr>.elementType", "");
        
              if (mrvAtom->id == "" || mrvAtom->elementType == "")
                throw FileParseException("Expected id, elementType for an atom definition in MRV file");
              
              std::string x2 = v.second.get<std::string>("<xmlattr>.x2", "");
              std::string y2 = v.second.get<std::string>("<xmlattr>.y2", "");

              // x2 and y2 are doubles
              
              if (x2 != "" && y2 != ""  && (!getCleanDouble(x2, mrvAtom->x2) || !getCleanDouble(y2, mrvAtom->y2)))
                throw FileParseException("The values x2 and y2 must be large floats in MRV file");        

              std::string formalCharge = v.second.get<std::string>("<xmlattr>.formalCharge", "");
              if (formalCharge != "")
              {
                if (!getCleanInt(formalCharge, mrvAtom->formalCharge) )
                  throw FileParseException("The value for formalCharge must be an integer in MRV file"); 
              }
              else
                mrvAtom->formalCharge = 0;
              

              mrvAtom->radical = v.second.get<std::string>("<xmlattr>.radical", "");
              if (mrvAtom->radical != "")
              {
                if (!boost::algorithm::contains(marvinRadicalVals, std::vector<std::string>{mrvAtom->radical} ))
                {
                  std::ostringstream err;
                  err << "The value for radical must be one of " << boost::algorithm::join(marvinRadicalVals,", ") << " in MRV file";
                  throw FileParseException(err.str()); 
                }
              }
              else 
                mrvAtom->radical = "";

              std::string isotopeStr = v.second.get<std::string>("<xmlattr>.isotope", "");
              if (isotopeStr != "")
              {
                if (!getCleanInt(isotopeStr, mrvAtom->isotope) || mrvAtom->isotope <=0)
                {
                  std::ostringstream err;
                  err << "The value for isotope must be a positive number in MRV file";
                  throw FileParseException(err.str()); 
                }
              }
              else 
                mrvAtom->isotope = 0;

              mrvAtom->mrvAlias = v.second.get<std::string>("<xmlattr>.mrvAlias", "");

              mrvAtom->rgroupRef = v.second.get<int>("<xmlattr>.rgroupRef", -1);

              mrvAtom->mrvStereoGroup = v.second.get<std::string>("<xmlattr>.mrvStereoGroup", "");
              if (mrvAtom->mrvStereoGroup == "0")
                mrvAtom->mrvStereoGroup = "";

              std::string mrvMap = v.second.get<std::string>("<xmlattr>.mrvMap", "");
              if (mrvMap != "")
                {
                  if (!getCleanInt(mrvMap, mrvAtom->mrvMap) || mrvAtom->mrvMap <=0)
                    throw FileParseException("The value for mrvMap must be an non-=negative integer in MRV file"); 
                }
              else
                mrvAtom->mrvMap = 0;

              mrvAtom->sgroupRef = v.second.get<std::string>("<xmlattr>.sgroupRef", "");

              if (role== "SuperatomSgroup")
              {
                mrvAtom->sgroupAttachmentPoint = v.second.get<std::string>("<xmlattr>.sgroupAttachmentPoint", "");

                //atom->setProp(common_properties::molAttachPoint, ival);
              }           
            }
          }
          else  // single line form of atoms
          {
            // <atomArray atomID="a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11" elementType="C C C C C C Cl C N O O" formalCharge="0 0 0 0 0 0 0 0 1 0 -1" lonePair="0 0 0 0 0 0 3 0 0 2 3" x2="-4.3334 -5.6670 -5.6670 -4.3334 -2.9997 -2.9997 -4.3335 -1.6660 -7.0007 -1.6660 -0.3323" y2="1.7693 0.9993 -0.5409 -1.3109 -0.5409 0.9993 3.3093 -1.3109 -1.3109 -2.8509 -0.5410"></atomArray>
          
            std::vector<std::string> atomIds;
            boost::algorithm::split(atomIds,atomID,boost::algorithm::is_space());
            size_t atomCount = atomIds.size();

            std::vector<std::string> elementTypes;
            std::string elementType = atomArray.get<std::string>("<xmlattr>.elementType", "");
            boost::algorithm::split(elementTypes,elementType,boost::algorithm::is_space());
          
            std::vector<std::string> x2s;
            std::string x2 = atomArray.get<std::string>("<xmlattr>.x2", "");
            boost::algorithm::split(x2s,x2,boost::algorithm::is_space());

            std::vector<std::string> y2s;
            std::string y2 = atomArray.get<std::string>("<xmlattr>.y2", "");
            boost::algorithm::split(y2s,y2,boost::algorithm::is_space());

            std::vector<std::string> formalCharges;
            std::string formalCharge = atomArray.get<std::string>("<xmlattr>.formalCharge", "");
            boost::algorithm::split(formalCharges,formalCharge,boost::algorithm::is_space());

            std::vector<std::string> isotopes;
            std::string isotope = atomArray.get<std::string>("<xmlattr>.isotope", "");
            boost::algorithm::split(isotopes,isotope,boost::algorithm::is_space());

            std::vector<std::string> radicals;
            std::string radical = atomArray.get<std::string>("<xmlattr>.radical", "");
            boost::algorithm::split(radicals,radical,boost::algorithm::is_space());

            std::vector<std::string> mrvAliases;
            std::string mrvAlias = atomArray.get<std::string>("<xmlattr>.mrvAlias", "");
            boost::algorithm::split(mrvAliases,mrvAlias,boost::algorithm::is_space());


            std::vector<std::string> rgroupRefs;
            std::string rgroupRef = atomArray.get<std::string>("<xmlattr>.rgroupRef", "");
            boost::algorithm::split(rgroupRefs,rgroupRef,boost::algorithm::is_space());
        
            std::vector<std::string> mrvStereoGroups;
            std::string mrvStereoGroup = atomArray.get<std::string>("<xmlattr>.mrvStereoGroup", "");
            boost::algorithm::split(mrvStereoGroups,mrvStereoGroup,boost::algorithm::is_space());

            std::vector<std::string> mrvMaps;
            std::string mrvMap = atomArray.get<std::string>("<xmlattr>.mrvMap", "");
            boost::algorithm::split(mrvMaps,mrvMap,boost::algorithm::is_space());

            std::vector<std::string> sgroupRefs;
            std::string sgroupRef = atomArray.get<std::string>("<xmlattr>.sgroupRef", "");
            boost::algorithm::split(sgroupRefs,sgroupRef,boost::algorithm::is_space());

            std::vector<std::string> sgroupAttachmentPoints;
            std::string sgroupAttachmentPoint = atomArray.get<std::string>("<xmlattr>.sgroupAttachmentPoint", "");
            boost::algorithm::split(sgroupAttachmentPoints,sgroupAttachmentPoint,boost::algorithm::is_space());

            if (atomID == "" ||elementType == "" ||  elementTypes.size() < atomCount )
              throw FileParseException("Expected id, and elementType arrays for an atomArray definition in MRV file");
            if (elementTypes.size() < atomCount )
              throw FileParseException("There must be an element type for each atom id");
            
            for (size_t i = 0 ; i < atomCount ; ++i)
            {
              MarvinAtom *mrvAtom = new MarvinAtom();   
              res->atoms.push_back(mrvAtom);

              mrvAtom->id = atomIds[i];

              mrvAtom->elementType = elementTypes[i];
              
              if (x2 != ""  && y2 != "" && x2s.size() > i && y2s.size() > i)
                if (!getCleanDouble(x2s[i], mrvAtom->x2) || !getCleanDouble(y2s[i], mrvAtom->y2))
                  throw FileParseException("The values x2 and y2 must be large floats in MRV file");        

              if (formalCharge != "" && formalCharges.size() > i)
              {
                if (!getCleanInt(formalCharges[i], mrvAtom->formalCharge) )
                  throw FileParseException("The value for formalCharge must be an integer in MRV file"); 
              }
              else
                mrvAtom->formalCharge = 0;

              if (isotope != "" && isotopes.size() > i)
              {
                if (!getCleanInt(isotopes[i], mrvAtom->isotope) )
                  throw FileParseException("The value for formalCharge must be an integer in MRV file"); 
              }
              else
                mrvAtom->isotope = 0;
            
            
              if (radical != "" && radicals.size() > i)
              {
                mrvAtom->radical = radicals[i];
                if (!boost::algorithm::contains(marvinRadicalVals, std::vector<std::string>{mrvAtom->radical} ))
                {
                  std::ostringstream err;
                  err << "The value for radical must be one of " << boost::algorithm::join(marvinRadicalVals,", ") << " in MRV file";
                  throw FileParseException(err.str()); 
                }
              }
              else 
                mrvAtom->radical = "";

              if (mrvAlias != "" && mrvAliases.size() > i)
                mrvAtom->mrvAlias = mrvAliases[i];
              else
                mrvAtom->mrvAlias = "";


              if (rgroupRef != "" && rgroupRefs.size() > i)
              {
                if (!getCleanInt(rgroupRefs[i],mrvAtom->rgroupRef))
                  throw FileParseException("rgroupRef value must be an integer in MRV file");
              }
              else
                mrvAtom->rgroupRef = (-1);                  
              
              if (mrvStereoGroup != "" && mrvStereoGroups.size() > i && mrvStereoGroups[i] != "0")  // "0" is NOT a stereo group
              {
                mrvAtom->mrvStereoGroup = mrvStereoGroups[i];
              }
              else
                mrvAtom->mrvStereoGroup = "";       
              
              if (mrvMap != "" && mrvMaps.size() > i)
              {
                if (!getCleanInt(mrvMaps[i], mrvAtom->mrvMap) || mrvAtom->mrvMap <=0)
                  throw FileParseException("The value for mrvMap must be an non-=negative integer in MRV file"); 
              }
              else
              mrvAtom->mrvMap = 0;       

              if (sgroupRef != "" && sgroupRefs.size() > i)
                mrvAtom->sgroupRef = sgroupRefs[i];
              else
                mrvAtom->sgroupRef = "";
              
              if (role== "SuperatomSgroup" && sgroupAttachmentPoint != "" && sgroupAttachmentPoints.size() > i)
              {
                mrvAtom->sgroupAttachmentPoint = sgroupAttachmentPoints[i];
              }
              else
                mrvAtom->sgroupAttachmentPoint = "";
            }
          }
        }

        // get bonds if this is NOT a role == "SruSgroup"

        if (role != "SruSgroup")
        {
          BOOST_FOREACH(
            boost::property_tree::ptree::value_type &v, molTree.get_child("bondArray"))
          {
          MarvinBond *mrvBond = new MarvinBond();   
          res->bonds.push_back(mrvBond);

          mrvBond->id = v.second.get<std::string>("<xmlattr>.id", "");
          if (mrvBond->id == ""  )
            throw FileParseException("Expected id for an bond definition in MRV file");
          
          
          std::string atomRefs2 = v.second.get<std::string>("<xmlattr>.atomRefs2", "");

          std::vector<std::string> atomRefs2s;
          boost::algorithm::split(atomRefs2s,atomRefs2,boost::algorithm::is_space());
          mrvBond->atomRefs2[0] =  atomRefs2s[0];
          mrvBond->atomRefs2[1] =  atomRefs2s[1];
          if (atomRefs2s.size() != 2 
            || !boost::algorithm::contains(res->atoms, std::vector<std::string>{mrvBond->atomRefs2[0]}, MarvinMol::atomRefInAtoms )
            || !boost::algorithm::contains(res->atoms, std::vector<std::string>{mrvBond->atomRefs2[1]}, MarvinMol::atomRefInAtoms ) )
            throw FileParseException("atomRefs2 must contain two atom refs that must appear in the atoms array in MRV file");


          mrvBond->order = v.second.get<std::string>("<xmlattr>.order", "");
          if (!boost::algorithm::contains(marvinBondOrders, std::vector<std::string>{mrvBond->order} ))
          {
            std::ostringstream err;
            err << "Expected one of  " << boost::algorithm::join(marvinBondOrders,", ") << " for order for an bond definition in MRV file";
            throw FileParseException(err.str());
          }

          mrvBond->queryType = v.second.get<std::string>("<xmlattr>.queryType", "");
          if (mrvBond->queryType != "")
          {
            if (!boost::algorithm::contains(marvinQueryBondsTypes, std::vector<std::string>{mrvBond->queryType} ))
            {
              std::ostringstream err;
              err << "Expected one of  " << boost::algorithm::join(marvinQueryBondsTypes,", ") << " for queryType for an bond definition in MRV file";
              throw FileParseException(err.str());
            }
          }

          mrvBond->bondStereo = v.second.get<std::string>("bondStereo", "");
          if (mrvBond->bondStereo == "" || boost::algorithm::to_lower_copy(mrvBond->bondStereo) == "w" || boost::algorithm::to_lower_copy(mrvBond->bondStereo) == "h")
          {
            // do nothing  - this is OK
          }
          else if (boost::algorithm::to_lower_copy(mrvBond->bondStereo) == "c" || boost::algorithm::to_lower_copy(mrvBond->bondStereo) == "t")
            mrvBond->bondStereo = "";  // cis and trans are ignored
          else
            throw FileParseException("The value for bondStereo must be \"H\", \"W\", \"C\" or \"T\" in MRV file (\"C\" and \"T\" are ignored)");
          }
        }

        if (role== "SuperatomSgroup")
        {

          // see if there is an AttachmentPointArray
          bool found;
          
          try 
          {
          boost::property_tree::ptree AttachmentPointArrayTree = molTree.get_child("AttachmentPointArray");
          found = true;

          }
          catch(const std::exception& e)
          {
          found = false;
          }

          if (found)
          {
          BOOST_FOREACH(
            boost::property_tree::ptree::value_type &v, molTree.get_child("AttachmentPointArray"))
            {
              MarvinAttachmentPoint *marvinAttachmentPoint = new MarvinAttachmentPoint();   
              ((MarvinSuperatomSgroup *)res)->attachmentPoints.push_back(marvinAttachmentPoint);

              marvinAttachmentPoint->atom = v.second.get<std::string>("<xmlattr>.atom", "");


              marvinAttachmentPoint->order = v.second.get<std::string>("<xmlattr>.order", "");
              marvinAttachmentPoint->bond = v.second.get<std::string>("<xmlattr>.bond", "");
              if (marvinAttachmentPoint->atom == "" || marvinAttachmentPoint->order == "" || marvinAttachmentPoint->bond == ""  )
                throw FileParseException("Expected atom, order and bond,  for an AttachmentPoint definition in MRV file");

              // atom must be found in the atoms vector

              if (!boost::algorithm::contains(res->atoms, std::vector<std::string>{marvinAttachmentPoint->atom}, MarvinMol::atomRefInAtoms ))
                throw FileParseException("Atom specification for an AttachmentPoint definition must be in the parent's atom array in MRV file");

              // bond must be found in the  bonds vector

              if (!boost::algorithm::contains(parentMol->bonds, std::vector<std::string>{marvinAttachmentPoint->bond}, MarvinMol::bondRefInBonds ))
                throw FileParseException("Bond specification for an AttachmentPoint definition must be in the bond array in MRV file");

            }
          }
        }

        if (parentMol == NULL)  // is this is a parent type mol - it has no parent
        {
          BOOST_FOREACH(
            boost::property_tree::ptree::value_type &v, molTree)
          {
            if(v.first == "molecule")
            {
              MarvinMolBase *subMol = (MarvinSuperatomSgroup *)parseMarvinMolecule(v.second, (MarvinMol *)res);

              if (subMol->role() == "SuperatomSgroup")
                ((MarvinMol *)res)->superatomSgroups.push_back((MarvinSuperatomSgroup *)subMol);
              else if (subMol->role() == "SruSgroup")
                ((MarvinMol *)res)->sruSgroups.push_back((MarvinSruSgroup *)subMol);
              else
                throw FileParseException("Unexpected role while parsing sub-molecues in MRV file");
            }
          }
        }

        return res;
      }
      catch(const std::exception& e)
      {
        delete res;

        throw;
      }
    };

    ChemicalReaction *parseReaction(boost::property_tree::ptree rxnTree, boost::property_tree::ptree documentTree, bool sanitize=false, bool removeHs=false)
    {
      ChemicalReaction *rxn = NULL;
      MarvinReaction *marvinReaction = NULL;
      RWMol *mol = NULL;
      
      try
      {

        rxn = new ChemicalReaction();

        marvinReaction= parseMarvinReaction(rxnTree, documentTree);
        marvinReaction->convertFromSuperAtoms();

        // get each reactant

        std::vector<MarvinMol *>::iterator molIter;
        for (molIter = marvinReaction->reactants.begin(); molIter != marvinReaction->reactants.end(); ++molIter)
        {
          mol = parseMolecule((*molIter) , sanitize, removeHs);
          ROMol *roMol = new ROMol(*mol);
          delete mol;

          rxn->addReactantTemplate(ROMOL_SPTR(roMol));   //roMol  now owned by rxn;
        }

        // get each agent

        for (molIter = marvinReaction->agents.begin(); molIter != marvinReaction->agents.end(); ++molIter)
        {
          mol = parseMolecule((*molIter) , sanitize, removeHs);
          ROMol *roMol = new ROMol(*mol);
          delete mol;
          
          rxn->addAgentTemplate(ROMOL_SPTR(roMol)); //roMol  now owned by rxn;
        }

        // get each product

        for (molIter = marvinReaction->products.begin(); molIter != marvinReaction->products.end(); ++molIter)
        {
          mol = parseMolecule((*molIter) , sanitize, removeHs);
          ROMol *roMol = new ROMol(*mol);
          delete mol;

          rxn->addProductTemplate(ROMOL_SPTR(roMol)); //roMol  now owned by rxn;
        }

        // convert atoms to queries:
        for (MOL_SPTR_VECT::const_iterator iter = rxn->beginReactantTemplates(); iter != rxn->endReactantTemplates(); ++iter) 
        {
          // to write the mol block, we need ring information:
          for (ROMol::AtomIterator atomIt = (*iter)->beginAtoms();
            atomIt != (*iter)->endAtoms(); ++atomIt) 
          {
            QueryOps::replaceAtomWithQueryAtom((RWMol *)iter->get(), (*atomIt));
          }
        }
        for (MOL_SPTR_VECT::const_iterator iter = rxn->beginProductTemplates();
          iter != rxn->endProductTemplates(); ++iter) 
        {
          // to write the mol block, we need ring information:
          for (ROMol::AtomIterator atomIt = (*iter)->beginAtoms(); atomIt != (*iter)->endAtoms(); ++atomIt) 
          {
            QueryOps::replaceAtomWithQueryAtom((RWMol *)iter->get(), (*atomIt));
          }
        }
        //updateProductsStereochem(rxn);

        // RXN-based reactions do not have implicit properties
        //rxn->setImplicitPropertiesFlag(false);

        delete marvinReaction;
        return rxn;
      }
      catch(const std::exception& e)
      {
        delete marvinReaction;
        delete rxn;

        throw;
      } 
    }


    MarvinReaction  *parseMarvinReaction(boost::property_tree::ptree rxnTree, boost::property_tree::ptree documentTree)
    {
      MarvinReaction *res = new MarvinReaction();

      try
      {
      
        boost::property_tree::ptree childTree;
        bool foundChild=false;
        try
        {
          childTree = rxnTree.get_child("reactantList");
          foundChild = true;
        }
        catch(const std::exception& e)
        {
          foundChild = false;
        }
        
        if (foundChild)
        {
          BOOST_FOREACH(
            boost::property_tree::ptree::value_type &v, childTree)
          res->reactants.push_back((MarvinMol *)parseMarvinMolecule(v.second));
        }

        try
        {
          childTree = rxnTree.get_child("agentList");
          foundChild = true;
        }
        catch(const std::exception& e)
        {
          foundChild = false;
        }
        if (foundChild)
        {
          BOOST_FOREACH(
            boost::property_tree::ptree::value_type &v, childTree)
          res->agents.push_back((MarvinMol *)parseMarvinMolecule(v.second));
        }

        try
        {
          childTree = rxnTree.get_child("productList");
          foundChild = true;
        }
        catch(const std::exception& e)
        {
          foundChild = false;
        }
        if (foundChild)
        {
          BOOST_FOREACH(
            boost::property_tree::ptree::value_type &v, childTree)
          res->products.push_back((MarvinMol *)parseMarvinMolecule(v.second));
        }

        // <arrow type="DEFAULT" x1="-11.816189911577812" y1="-10.001443743444021" x2="-8.401759471454618" y2="-10.001443743444021"/>
        boost::property_tree::ptree arrow = rxnTree.get_child("arrow");
        res->arrow.type = arrow.get<std::string>("<xmlattr>.type", "");
        if (!getCleanDouble( arrow.get<std::string>("<xmlattr>.x1", ""), res->arrow.x1)
            ||!getCleanDouble( arrow.get<std::string>("<xmlattr>.y1", ""), res->arrow.y1)
            ||!getCleanDouble( arrow.get<std::string>("<xmlattr>.x2", ""), res->arrow.x2)
            ||!getCleanDouble( arrow.get<std::string>("<xmlattr>.y2", ""), res->arrow.y1))
          throw FileParseException("Arrow coordinates must all be large floating point numbers in MRV file"); 
          
        BOOST_FOREACH(
          boost::property_tree::ptree::value_type &v, documentTree)
        {

          if (v.first != "MReactionSign")
          continue;
          MarvinPlus *marvinPlus = new MarvinPlus();
          res->pluses.push_back(marvinPlus);  
          marvinPlus->id =  v.second.get<std::string>("<xmlattr>.id", "");
          int pointCount = 0;
          BOOST_FOREACH(
            boost::property_tree::ptree::value_type &v2, v.second)
          {
            if (v2.first == "MPoint")
            {

              double x;
              double y; 
              if (!getCleanDouble( v2.second.get<std::string>("<xmlattr>.x", ""),x)
                || !getCleanDouble( v2.second.get<std::string>("<xmlattr>.y", ""),y))
              throw FileParseException("Plus sign  coordinates must all be large floating point numbers in MRV file"); 
                
              switch( pointCount)
              {
              case 0:  //  first point - x1 and y1 are set
                marvinPlus->x1 = x;
                marvinPlus->y1 = y;
                break;
              case 1:   // x2 is set, y1 is checked
                marvinPlus->x2 = x;
                if (marvinPlus->y1 != y)
                  throw FileParseException("Plus sign coordinate Y in 2nd MPoint must be the same as that from the 1st MPoint in MRV file"); 
                break;
              case 2:   // y2 is set, x2 is checked
                marvinPlus->y2 = y;
                if (marvinPlus->x2 != x)
                  throw FileParseException("Plus sign coordinate X in 3rd MPoint must be the same as that from the 2nd MPoint in MRV file"); 
                break;
              case 3:   // x2 and y2 are checked
                if (marvinPlus->x1 != x)
                  throw FileParseException("Plus sign coordinate X in 4th MPoint must be the same as that from the 1st MPoint in MRV file"); 

                if (marvinPlus->y2 != y)
                  throw FileParseException("Plus sign coordinate Y in 4th MPoint must be the same as that from the 3rd MPoint in MRV file"); 
                break;

              default:
                throw FileParseException("Plus sign coordinate must have 4 MPoints in MRV file"); 
              }
              ++pointCount;
            }
          }
        }

        BOOST_FOREACH(
            boost::property_tree::ptree::value_type &v, documentTree)
        {
          if (v.first != "MTextBox")
          continue;

          MarvinCondition *marvinCondition = new MarvinCondition();
          res->conditions.push_back(marvinCondition);  
          marvinCondition->id =  v.second.get<std::string>("<xmlattr>.id", "");
          marvinCondition->halign =  v.second.get<std::string>("<xmlattr>.halign", "");
          marvinCondition->valign =  v.second.get<std::string>("<xmlattr>.valign", "");
          double fontScale;
          std::string fontScaleStr = v.second.get<std::string>("<xmlattr>.fontScale", "");
          if (fontScaleStr != "")
          {
            if (!getCleanDouble(fontScaleStr ,fontScale))
              throw FileParseException(
                  "Condition font scale must be a positive integer in MRV file"); 
          }
          else
            fontScale = 0.0;
          marvinCondition->fontScale = fontScale;

          marvinCondition->text = v.second.get<std::string>("Field", "");

          int pointCount = 0;
          BOOST_FOREACH(
            boost::property_tree::ptree::value_type &v2, v.second)
          {
            if (v2.first == "MPoint")
            {
              double x,y;
              std::string xStr = v2.second.get<std::string>("<xmlattr>.x", "");
              std::string yStr = v2.second.get<std::string>("<xmlattr>.y", "");
              if (! getCleanDouble(xStr,  x) ||  ! getCleanDouble(yStr,  y))
                throw FileParseException("Condition coordinate must valid integers in MRV file"); 
            
              switch( pointCount)
              {
              case 0:  //  first point - x1 and y1 are set
                marvinCondition->x = x;
                marvinCondition->y = y;
                break;
              case 1:   
              case 2:
              case 3:   // x and Y are checked - must be the same as point 1
                if ( marvinCondition->x != x || marvinCondition->y != y)
                  throw FileParseException("Condition coordinates must be the same as those from the 1st MPoint in MRV file"); 
                break;

              default:
                throw FileParseException("Condition defs must have 4 MPoints in MRV file"); 
              }
              ++pointCount;
            }
          }
        }
    
        return res;
      }
      catch(const std::exception& e)
      {
        delete res;

        throw;
      }
    }
  };    

  //------------------------------------------------
  //
  //  Read a molecule from a stream
  //
  //------------------------------------------------
  
  void *MrvDataStreamParser(std::istream *inStream, bool &isReaction, bool sanitize, bool removeHs)
  {
    PRECONDITION(inStream, "no stream");
    std::string tempStr;
    void *res = NULL;
    MarvinCMLReader marvinCML;

    res = marvinCML.parse(*inStream, isReaction, sanitize, removeHs);

    // if (marvinCML.isReaction())
    //   res = (void *)marvinCML.rxn;
    // else
    //   res = (void *)marvinCML.mol;

    return res;
  }

  void *MrvDataStreamParser(std::istream &inStream, bool &isReaction, bool sanitize, bool removeHs)
  {
    return MrvDataStreamParser(&inStream, isReaction, sanitize, removeHs);
  }
  //------------------------------------------------
  //
  //  Read a molecule from a string
  //
  //------------------------------------------------
  void *MrvBlockParser(const std::string &molmrvText, bool &isReaction, bool sanitize, bool removeHs)
  {
    std::istringstream inStream(molmrvText);
    // unsigned int line = 0;
    return MrvDataStreamParser(inStream, isReaction, sanitize, removeHs);
  }

  //------------------------------------------------
  //
  //  Read a RWMOL from a file
  //
  //------------------------------------------------
  void *MrvFileParser(const std::string &fName, bool &isReaction, bool sanitize, bool removeHs)
  {
    std::ifstream inStream(fName.c_str());
    if (!inStream || (inStream.bad()))
    {
      std::ostringstream errout;
      errout << "Bad input file " << fName;
      throw BadFileException(errout.str());
    }
    void *res = nullptr;
    if (!inStream.eof())
    {
      res = MrvDataStreamParser(inStream,  isReaction, sanitize, removeHs);
    }
    return res;
  }

  //------------------------------------------------
  //
  //  Read a RWMol from a stream 
  //
  //------------------------------------------------
  RWMol *MrvMolDataStreamParser(std::istream *inStream, bool sanitize, bool removeHs)
  {
  
    void *res = NULL;
    
    bool isReaction=false;
    res = MrvDataStreamParser(inStream, isReaction, sanitize, removeHs);
    if (isReaction)
    {
      delete (ChemicalReaction *)res;
      throw FileParseException("The file parsed as a reaction, not a molecule"); 
    }

    return (RWMol *)res;
  }
  //------------------------------------------------
  //
  //  Read a RWMol from a stream reference
  //
  //------------------------------------------------
  RWMol *MrvMolDataStreamParser(std::istream &inStream, bool sanitize, bool removeHs)
  {
    return MrvMolDataStreamParser(&inStream, sanitize, removeHs);
  }
  //------------------------------------------------
  //
  //  Read a RWMol from a string
  //
  //------------------------------------------------
  RWMol *MrvMolStringParser(const std::string &molmrvText, bool sanitize, bool removeHs)
  {
    std::istringstream inStream(molmrvText);
    // unsigned int line = 0;
    return MrvMolDataStreamParser(inStream, sanitize, removeHs);
  }

  //------------------------------------------------
  //
  //  Read an RWMol from a file
  //
  //------------------------------------------------
  RWMol *MrvMolFileParser(const std::string &fName, bool sanitize, bool removeHs)
  {
    std::ifstream inStream(fName.c_str());
    if (!inStream || (inStream.bad()))
    {
      std::ostringstream errout;
      errout << "Bad input file " << fName;
      throw BadFileException(errout.str());
    }
    RWMol *res = nullptr;
    if (!inStream.eof())
      res = MrvMolDataStreamParser(inStream, sanitize, removeHs);
    return res;
  }

  //------------------------------------------------
  //
  //  Read a ChemicalReaction from a stream 
  //
  //------------------------------------------------
  ChemicalReaction *MrvRxnDataStreamParser(std::istream *inStream, bool sanitize, bool removeHs)
  {
    void *res = NULL;
    
    bool isReaction=false;
    res = MrvDataStreamParser(inStream, isReaction, sanitize, removeHs);
    if (!isReaction)
    {
      delete (RWMol *)res;
      throw FileParseException("The file parsed as a molecule, not a reaction"); 
    }

    return (ChemicalReaction *)res;
  }

  //------------------------------------------------
  //
  //  Read a ChemicalReaction from a stream reference
  //
  //------------------------------------------------
  ChemicalReaction *MrvRxnDataStreamParser(std::istream &inStream, bool sanitize, bool removeHs)
  {
    return MrvRxnDataStreamParser(&inStream, sanitize, removeHs);
  }
  //------------------------------------------------
  //
  //  Read a ChemicalReaction from a string
  //
  //------------------------------------------------
  ChemicalReaction *MrvRxnStringParser(const std::string &molmrvText, bool sanitize, bool removeHs)
  {
    std::istringstream inStream(molmrvText);
    // unsigned int line = 0;
    return MrvRxnDataStreamParser(inStream, sanitize, removeHs);
  }

  //------------------------------------------------
  //
  //  Read a ChemicalReaction from a file
  //
  //------------------------------------------------
  ChemicalReaction *MrvRxnFileParser(const std::string &fName, bool sanitize, bool removeHs)
  {
    std::ifstream inStream(fName.c_str());
    if (!inStream || (inStream.bad()))
    {
      std::ostringstream errout;
      errout << "Bad input file " << fName;
      throw BadFileException(errout.str());
    }
    ChemicalReaction *res = MrvRxnDataStreamParser(inStream, sanitize, removeHs);
    return res;
  }
}// namespace RDKit
