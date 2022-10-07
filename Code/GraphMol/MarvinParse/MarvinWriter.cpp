//
//  Copyright (C) 2002-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

// this file WRITES MRV file for molecules and reactions

#include <string>
#include <exception>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cstdio>

#include <GraphMol/SubstanceGroup.h>
#include <RDGeneral/Ranking.h>
#include <RDGeneral/LocaleSwitcher.h>
#include <RDGeneral/Invariant.h>


#include <boost/format.hpp>
#include <boost/dynamic_bitset.hpp>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/GenericGroups/GenericGroups.h>

#include <boost/algorithm/string.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSGroupWriting.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include "MarvinParser.h"
#include "MarvinDefs.h"


#include <GraphMol/RDKitQueries.h>
#include <GraphMol/StereoGroup.h>
#

#include <RDGeneral/StreamOps.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/LocaleSwitcher.h>

using namespace RDKit::SGroupWriting;

namespace RDKit 
{

  bool hasComplexQuery(const Atom *atom) 
  {
    PRECONDITION(atom, "bad atom");
    bool res = false;
    if (atom->hasQuery()) 
    {
      res = true;
      // counter examples:
      //  1) atomic number
      //  2) the smarts parser inserts AtomAnd queries
      //     for "C" or "c":
      //
      std::string descr = atom->getQuery()->getDescription();
      if (descr == "AtomAtomicNum") 
        res = false;
      else if (descr == "AtomAnd") 
      {
        if ((*atom->getQuery()->beginChildren())->getDescription() == "AtomAtomicNum") 
          res = false;
      }
    }
    return res;
  }

  const std::string GetMarvinAtomSymbol(
      const Atom *atom) 
  {
    PRECONDITION(atom, "");

    std::string res;
    if (atom->hasProp(common_properties::_MolFileRLabel)) 
    {
      unsigned int lbl;
      if (atom->getPropIfPresent(common_properties::_MolFileRLabel, lbl) && lbl > 0)
        res = "R" + std::to_string(lbl);
      else 
        res = "R";
    } 
    else if (atom->getAtomicNum()) 
    {
      res = atom->getSymbol();
    } 
    else 
    {
      if (!atom->hasProp(common_properties::dummyLabel)) 
      {
        if (atom->hasQuery() &&
            (atom->getQuery()->getTypeLabel() == "A" ||
            (atom->getQuery()->getNegation() &&
              atom->getQuery()->getDescription() == "AtomAtomicNum" &&
              static_cast<ATOM_EQUALS_QUERY *>(atom->getQuery())->getVal() ==
                  1))) 
        {
          res = "*";
        } 
        else if (atom->hasQuery() &&
                  (atom->getQuery()->getTypeLabel() == "Q" ||
                    (atom->getQuery()->getNegation() &&
                    atom->getQuery()->getDescription() == "AtomOr" &&
                    atom->getQuery()->endChildren() -
                            atom->getQuery()->beginChildren() ==
                        2 &&
                    (*atom->getQuery()->beginChildren())->getDescription() ==
                        "AtomAtomicNum" &&
                    static_cast<ATOM_EQUALS_QUERY *>(
                        (*atom->getQuery()->beginChildren()).get())
                            ->getVal() == 6 &&
                    (*++(atom->getQuery()->beginChildren()))->getDescription() ==
                        "AtomAtomicNum" &&
                    static_cast<ATOM_EQUALS_QUERY *>(
                        (*++(atom->getQuery()->beginChildren())).get())
                            ->getVal() == 1))) 
        {
          throw MarvinWriterException("Query atoms are not supported for MarvinWriter");
        } 
        else if (atom->hasQuery() && atom->getQuery()->getTypeLabel() == "X") 
        {
          throw MarvinWriterException("Query atoms are not supported for MarvinWriter");
        } 
        else if (atom->hasQuery() && atom->getQuery()->getTypeLabel() == "M") 
        {
          throw MarvinWriterException("Query atoms are not supported for MarvinWriter");
        } 
        else if (atom->hasQuery() && atom->getQuery()->getTypeLabel() == "AH") 
        {
          throw MarvinWriterException("Query atoms are not supported for MarvinWriter");
        } 
        else if (atom->hasQuery() && atom->getQuery()->getTypeLabel() == "QH") 
        {
          throw MarvinWriterException("Query atoms are not supported for MarvinWriter");
        } 
        else if (atom->hasQuery() && atom->getQuery()->getTypeLabel() == "XH") 
        {
          throw MarvinWriterException("Query atoms are not supported for MarvinWriter");
        } 
        else if (atom->hasQuery() && atom->getQuery()->getTypeLabel() == "MH") 
        {
          throw MarvinWriterException("Query atoms are not supported for MarvinWriter");
        } 
        else if (hasComplexQuery(atom)) 
        {
          throw MarvinWriterException("Query atoms are not supported for MarvinWriter");
        } 
        else 
        {
          res = "R";
        }
      } 
      else 
      {
        std::string symb;
        atom->getProp(common_properties::dummyLabel, symb);
        if (symb == "*") 
          res = "R";
        else if (symb == "X") 
          res = "R";
        else if (symb == "Xa") 
          res = "R1";
        else if (symb == "Xb") 
          res = "R2";
        else if (symb == "Xc") 
          res = "R3";
        else if (symb == "Xd") 
          res = "R4";
        else if (symb == "Xf") 
          res = "R5";
        else if (symb == "Xg") 
          res = "R6";
        else if (symb == "Xh") 
          res = "R7";
        else if (symb == "Xi") 
          res = "R8";
        else if (symb == "Xj") 
          res = "R9";
        else 
          throw MarvinWriterException("Query atoms are not supported for MarvinWriter");
      }
    
    }
    
    return res;
  }

  bool isQueryBondInRing(const Bond *bond) 
  {
    PRECONDITION(bond, "no bond");
    PRECONDITION(bond->hasQuery(), "no query");
    Bond::QUERYBOND_QUERY *qry = bond->getQuery();
    // start by catching combined bond order + bond topology queries

    if (qry->getDescription() == "BondAnd" && !qry->getNegation() &&
        qry->endChildren() - qry->beginChildren() == 2) 
    {
      auto child1 = qry->beginChildren();
      auto child2 = child1 + 1;
      if (((*child1)->getDescription() == "BondInRing") != ((*child2)->getDescription() == "BondInRing")) 
      {
        if ((*child1)->getDescription() != "BondInRing") 
          qry = child2->get();
        else
          qry = child1->get();
      }
    }
    if (qry->getDescription() == "BondInRing") 
      return true;

    return false;
  }

  std::string getMarvinQueryBondSymbol(const Bond *bond) 
  {
    PRECONDITION(bond, "no bond");
    PRECONDITION(bond->hasQuery(), "no query");

    Bond::QUERYBOND_QUERY *qry = bond->getQuery();
    if (qry->getDescription() == "BondOrder" || isQueryBondInRing(bond)) 
      return "";
    else 
    {
      // start by catching combined bond order + bond topology queries
      if (qry->getDescription() == "BondAnd" && !qry->getNegation() &&
          qry->endChildren() - qry->beginChildren() == 2) 
      {
        auto child1 = qry->beginChildren();
        auto child2 = child1 + 1;
        if ((*child2)->getDescription() == "BondInRing") 
          qry = child1->get();
        else if ((*child1)->getDescription() == "BondInRing") 
          qry = child2->get();
      }
      if (qry->getDescription() == "BondOr" && !qry->getNegation()) 
      {
        if (qry->endChildren() - qry->beginChildren() == 2) 
        {
          auto child1 = qry->beginChildren();
          auto child2 = child1 + 1;
          if ((*child1)->getDescription() == "BondOrder" &&
              !(*child1)->getNegation() &&
              (*child2)->getDescription() == "BondOrder" &&
              !(*child2)->getNegation()) 
          {
            // ok, it's a bond query we have a chance of dealing with
            int t1 = static_cast<BOND_EQUALS_QUERY *>(child1->get())->getVal();
            int t2 = static_cast<BOND_EQUALS_QUERY *>(child2->get())->getVal();
            if (t1 > t2) 
              std::swap(t1, t2);
            if (t1 == Bond::SINGLE && t2 == Bond::DOUBLE) 
              return "SD";
            else if (t1 == Bond::SINGLE && t2 == Bond::AROMATIC) 
              return "SA";
            else if (t1 == Bond::DOUBLE && t2 == Bond::AROMATIC) 
              return "DA";
          }
        }
      } 
      else if (qry->getDescription() == "SingleOrAromaticBond" && !qry->getNegation()) 
       return "SA";
      else if (qry->getDescription() == "SingleOrDoubleBond" && !qry->getNegation()) 
        return "SD";
      else if (qry->getDescription() == "DoubleOrAromaticBond" && !qry->getNegation()) 
        return "DA";
  
    }
    throw MarvinWriterException("Only SA, DA, and SD query bond are supported for MarvinWriter");
  }

 

  void GetMarvinBondSymbol(const Bond *bond, std::string &order, std::string &queryType)
  {
    PRECONDITION(bond, "");
    
    order = "1";
    if (bond->hasQuery())
    {
      order = "1";
      queryType = getMarvinQueryBondSymbol(bond);
      if (queryType == "")
         throw MarvinWriterException("Only 1,2,3,Aromatic, and query bonds SA, DA, and SD are supported for MarvinWriter");
      return;
    }
    
    queryType = "";   // not s query

    switch (bond->getBondType()) 
    {
      case Bond::SINGLE:
        if (bond->getIsAromatic()) 
           queryType = "A";
        else
          queryType =  "1";
        break;

      case Bond::DOUBLE:
        if (bond->getIsAromatic()) 
           queryType = "A";
        else
          queryType =  "2";
        break;
      case Bond::TRIPLE:
        queryType = "3";
        break;

      case Bond::AROMATIC:
        queryType =  "A";
        break;

      default:
        throw MarvinWriterException("Only 1,2,3,Aromatic, and query bonds SA, DA, and SD are supported for MarvinWriter");
    }
  }

  
  void GetMarvinBondStereoInfo(const Bond *bond, const INT_MAP_INT &wedgeBonds,
                                const Conformer *conf,  Bond::BondDir &dir,
                                bool &reverse) 
  {
    PRECONDITION(bond, "");
    reverse = false;
    dir = Bond::NONE;
    if (bond->getBondType() == Bond::SINGLE) 
    {
      // single bond stereo chemistry
      dir = DetermineBondWedgeState(bond, wedgeBonds, conf);

      // if this bond needs to be wedged it is possible that this
      // wedging was determined by a chiral atom at the end of the
      // bond (instead of at the beginning). In this case we need to
      // reverse the begin and end atoms for the bond when we write
      // the mol file

      if ((dir == Bond::BEGINDASH) || (dir == Bond::BEGINWEDGE)) 
      {
        auto wbi = wedgeBonds.find(bond->getIdx());
        if (wbi != wedgeBonds.end() && static_cast<unsigned int>(wbi->second) != bond->getBeginAtomIdx()) 
          reverse = true;
      }
    } 
  }



  MarvinMol *MolToMarvinMol(ROMol *mol, int &molCount, int &atomCount, int &bondCount, int &sgCount)
  {
    //molCount is the starting and ending molCount - used when called from a rxn
    
    MarvinMol *marvinMol = NULL;
    const Conformer *conf = NULL;

    try
    {
      marvinMol = new MarvinMol();

      marvinMol->molID = 'm' + std::to_string(++molCount);

      // get a 2D conformer

      for (unsigned int confId = 0; confId < mol->getNumConformers(); ++confId) 
      {
        auto testConf = mol->getConformer(confId);
        if (!testConf.is3D())
        { 
          conf = &testConf;
          break;
        }
      }

      for (auto atom : mol->atoms())
      {
        auto marvinAtom = new MarvinAtom();
        marvinMol->atoms.push_back(marvinAtom);

        marvinAtom->id = 'a' + std::to_string(++atomCount);

        marvinAtom->elementType = GetMarvinAtomSymbol(atom);
        
        marvinAtom->formalCharge = atom->getFormalCharge();
        marvinAtom->isotope = atom->getIsotope();
        
        if (!atom->getPropIfPresent(common_properties::molAtomMapNumber, marvinAtom->mrvMap))
          marvinAtom->mrvMap=0;

        if (conf != NULL)
        {
          const RDGeom::Point3D pos = conf->getAtomPos(atom->getIdx());
          marvinAtom->x2 = pos.x;
          marvinAtom->y2 = pos.y;
        }
        else
        {
          marvinAtom->x2 = 0.0;
          marvinAtom->y2 = 0.0;
        }
   
        unsigned int nRadEs = atom->getNumRadicalElectrons();
        if (nRadEs != 0)
           throw MarvinWriterException("Radicals are not handled by MarvinWriter"); 
       
        // atom maps for rxns

        if (!atom->getPropIfPresent(common_properties::molAtomMapNumber, marvinAtom->mrvMap)) 
          marvinAtom->mrvMap = 0;
      }

      INT_MAP_INT wedgeBonds = pickBondsToWedge(*mol);
      
      for (auto bond : mol->bonds())
      {
        auto marvinBond = new MarvinBond();
        marvinMol->bonds.push_back(marvinBond);

        marvinBond->id = 'b' + std::to_string(++bondCount);

        GetMarvinBondSymbol(bond, marvinBond->order, marvinBond->queryType );
        

        Bond::BondDir bondDirection;
        bool reverse;
        GetMarvinBondStereoInfo(bond, wedgeBonds, conf, bondDirection, reverse);

       
        if (reverse) 
        {
          // switch the begin and end atoms on the bond line
          marvinBond->atomRefs2[0] = marvinMol->atoms[bond->getEndAtomIdx()]->id;
          marvinBond->atomRefs2[1] = marvinMol->atoms[bond->getBeginAtomIdx()]->id;
        } 
        else 
        {
          marvinBond->atomRefs2[0] = marvinMol->atoms[bond->getBeginAtomIdx()]->id;
          marvinBond->atomRefs2[1] = marvinMol->atoms[bond->getEndAtomIdx()]->id;
        }

        switch (bondDirection)
        {
          case Bond::NONE:
            marvinBond->bondStereo = "";
            break;
          case Bond::BEGINWEDGE:
              marvinBond->bondStereo = "w";
              break;
          case Bond::BEGINDASH:
                      marvinBond->bondStereo = "d";
              break;
          default:
            throw MarvinWriterException("Unsupported bond direction");
        }       
      }

       // get all stereoGroups - add them to the correct atoms
      int orCount=0;
      int andCount=0;
      int absCount=0;  // really should be only one, but maybe ...
      int *currentCount;

      for (const StereoGroup group : mol->getStereoGroups())
      {
        std::string stereoGroupType;

        switch (group.getGroupType()) {
          case RDKit::StereoGroupType::STEREO_ABSOLUTE:
            stereoGroupType = "ABS";
            currentCount = &absCount;
            break;
          case RDKit::StereoGroupType::STEREO_OR:
            stereoGroupType = "OR";
            currentCount = &orCount;
            break;
          case RDKit::StereoGroupType::STEREO_AND:
              stereoGroupType = "AND";
              currentCount = &andCount;
            break;
          default:
            throw MarvinWriterException("Unrecognized stereo group type"); 
        }
        for (auto &&atom : group.getAtoms()) 
          marvinMol->atoms[atom->getIdx()]->mrvStereoGroup = stereoGroupType + std::to_string(++(*currentCount));
      }

      int sruMolCount = 0, sruSgCount=0;
      for (const SubstanceGroup &sgroup:  getSubstanceGroups(*mol))
      {
        std::string type = sgroup.getProp<std::string>("TYPE");
        if (type == "SRU")
        {
                            
          auto marvinSruSgroup =new MarvinSruSgroup();
          marvinMol->sruSgroups.push_back(marvinSruSgroup);

          marvinSruSgroup->title = sgroup.getProp<std::string>(std::string("LABEL"));
            
          marvinSruSgroup->connect = sgroup.getProp<std::string>("CONNECT");
          marvinSruSgroup->id = "sg" + std::to_string(++sruSgCount);
          marvinSruSgroup->molID  = 'm' + std::to_string(++sruMolCount);


          for (auto atomIndex : sgroup.getAtoms())
          {
              marvinSruSgroup->atomRefs.push_back(marvinMol->atoms[atomIndex]->id);
              marvinMol->atoms[atomIndex]->sgroupRef =  marvinSruSgroup->id;
          }
          
          for (auto bondIndex : sgroup.getBonds())
          {
              marvinSruSgroup->bondList.push_back(marvinMol->bonds[bondIndex]->id);
          }
        }

        else if (type == "SUP")
        {
          auto marvinSuperInfo =new MarvinSuperInfo();
          marvinMol->superInfos.push_back(marvinSuperInfo);

          marvinSuperInfo->title = sgroup.getProp<std::string>(std::string("LABEL"));
            
          for (auto atomIndex : sgroup.getAtoms())
              marvinSuperInfo->atoms.push_back(marvinMol->atoms[atomIndex]->id);
        }
          
      }

      // convert the superInfos to supergroups

      marvinMol->convertFromSuperAtoms();
      marvinMol->cleanUpNumbering(molCount, atomCount, bondCount, sgCount);

      return marvinMol;
    }
    catch(const std::exception& e)
    {
      delete marvinMol;
      throw;
    }
  }

  MarvinMol *MolToMarvinMol(ROMol *mol) 
  {
      int molCount=0, atomCount = 0, bondCount = 0, sgCount = 0;
      return MolToMarvinMol(mol, molCount, atomCount, bondCount, sgCount);
  }


  bool compareRowsOfRectangles(std::vector<MarvinRectangle> &v1, std::vector<MarvinRectangle> &v2)
  {
    return MarvinRectangle::compareRectanglesByY(v1.front(), v2.front());  // just compare the first one in each row
  }


  void AddMarvinPluses(MarvinReaction &rxn, std::vector<MarvinMol *>molList, int &arrowPlusCount)
  {
    // dividing the mols into rows and sorted by y value

    std::vector<std::vector<MarvinRectangle>> rowsOfRectangles;
    for (auto mol : molList)
    {
      // see if there is horizontal overlap with any existing row

      MarvinRectangle molRect(mol->atoms);
      bool foundRow = false;
      for (std::vector<MarvinRectangle> row : rowsOfRectangles)
      {
          for ( MarvinRectangle rect : row)
          {
            if (molRect.overlapsVertically(rect))
            {
              foundRow = true;
              row.push_back(molRect);
            }
          }
      }

      if (!foundRow)  // no overlap with a current row, so make a new one
      {
        std::vector<MarvinRectangle> newRow;
        newRow.push_back(molRect);
        rowsOfRectangles.push_back(newRow);
      }
    }

    // sort the members of each  row by X

    for (std::vector<MarvinRectangle> row : rowsOfRectangles)
        std::sort(row.begin(), row.end(), MarvinRectangle::compareRectanglesByX);
        
    //sort the rows by Y

    std::sort(rowsOfRectangles.begin(), rowsOfRectangles.end(), compareRowsOfRectangles);

    // make a plus between each rect on each row

    for (auto rowPtr = rowsOfRectangles.begin(); rowPtr != rowsOfRectangles.end() ; ++ rowPtr)
    {
      for (auto rect1 = rowPtr->begin() ; rect1 != rowPtr->end() ; ++rect1)
      {
          auto rect2 = rect1 + 1;
          if (rect2 == rowPtr->end())
            break;
          
          double x = (rect2->lowerRight.x + rect1->upperLeft.x)/2.0;
          double y;

          // see if there is room between for the +
          if (rect2->lowerRight.x - rect1->upperLeft.x >= 0.5)
              y = (rect2->getCenter().y + rect1->getCenter().y)/2.0;
          else  // put it under the two rectangles
            y = std::min<double>(rect1->lowerRight.y, rect2->lowerRight.y) - 0.25;

          auto newMarvinPlus = new MarvinPlus();
          rxn.pluses.push_back(newMarvinPlus);

          newMarvinPlus->id = "o" + std::to_string(++arrowPlusCount);
          newMarvinPlus->x1 = x - 0.25;
          newMarvinPlus->y1 = y - 0.25;
          newMarvinPlus->x2 = x + 0.25;
          newMarvinPlus->y2 = y + 0.25;
      }

    // for each row after the first, add plus

      if (rowPtr != rowsOfRectangles.begin())   // 2nd and subsequent rows
      {
        auto newMarvinPlus = new MarvinPlus();
        rxn.pluses.push_back(newMarvinPlus);

        double x = rowPtr->front().upperLeft.x;
        double y = rowPtr->front().getCenter().y;
        newMarvinPlus->id = "o" + std::to_string(++arrowPlusCount);
        newMarvinPlus->x1 = x - 0.25;
        newMarvinPlus->y1 = y - 0.25;
        newMarvinPlus->x2 = x + 0.25;
        newMarvinPlus->y2 = y + 0.25;
      }
    }



  }

  MarvinReaction *MolToMarvinRxn(ChemicalReaction *rxn)
  {
    MarvinReaction *marvinReaction = NULL;
    try
    {
      auto marvinReaction = new MarvinReaction();
      int molCount=0, atomCount = 0, bondCount = 0, sgCount = 0;
      for (auto mol : rxn->getReactants())
        marvinReaction->reactants.push_back(MolToMarvinMol(&(*mol), molCount, atomCount, bondCount, sgCount));
      for (auto mol : rxn->getAgents())
        marvinReaction->agents.push_back(MolToMarvinMol(&(*mol), molCount, atomCount, bondCount, sgCount));
      for (auto mol : rxn->getProducts())
        marvinReaction->products.push_back(MolToMarvinMol(&(*mol), molCount, atomCount, bondCount, sgCount));

      // make up some arrows

      int arrowPlusCount = 0;
      AddMarvinPluses(*marvinReaction, marvinReaction->reactants, arrowPlusCount);
      AddMarvinPluses(*marvinReaction, marvinReaction->products, arrowPlusCount);

      // add a reaction arrow


      // get the overall rectangle for the reactants and the one for the products

      if (marvinReaction->reactants.size() > 0 && marvinReaction->products.size() > 0)
      {
        marvinReaction->arrow.type = "DEFAULT";

        MarvinRectangle reactantRect(marvinReaction->reactants.front()->atoms);
        for (auto reactantPtr = marvinReaction->reactants.begin()+1 ; reactantPtr != marvinReaction->reactants.end(); ++reactantPtr )
          reactantRect.extend(MarvinRectangle((*reactantPtr)->atoms));

        MarvinRectangle productRect(marvinReaction->products.front()->atoms);
        for (auto productPtr = marvinReaction->products.begin()+1 ; productPtr != marvinReaction->products.end(); ++productPtr)
          productRect.extend(MarvinRectangle((*productPtr)->atoms));

        // if there is room between the reactants and products, put the arrow there

        if (productRect.upperLeft.x - reactantRect.lowerRight.x  > 2.0)
        {
          marvinReaction->arrow.x1 = reactantRect.lowerRight.x + 0.5;
          marvinReaction->arrow.x2 = productRect.upperLeft.x - 0.5;
          marvinReaction->arrow.y1 =(reactantRect.getCenter().y + productRect.getCenter().y)/2.0;
          marvinReaction->arrow.y2 = marvinReaction->arrow.y1;
        }
        // if not enough room between the reactants and product horizontally, try vertically

        else if (productRect.upperLeft.y - reactantRect.lowerRight.y  > 2.0)
        {
          marvinReaction->arrow.x1 = (reactantRect.getCenter().x + productRect.getCenter().x)/2.0;
          marvinReaction->arrow.x2 = marvinReaction->arrow.x1;
          marvinReaction->arrow.y1 =reactantRect.lowerRight.y + 0.5;
          marvinReaction->arrow.y2 =productRect.upperLeft.y - 0.5;
        }

        // if not good horizontal nor vertical place, just put it between the centers (hack)

        else if ((reactantRect.getCenter() - productRect.getCenter()).length() > 2.0)
        {
          marvinReaction->arrow.x1 = reactantRect.getCenter().x;
          marvinReaction->arrow.x2 = productRect.getCenter().x;
          marvinReaction->arrow.y1 = reactantRect.getCenter().y;
          marvinReaction->arrow.y2 = productRect.getCenter().y;
        }

        // really no good place for the arrow

        else
        {
          marvinReaction->arrow.x1 = reactantRect.lowerRight.x + 0.5;
          marvinReaction->arrow.x2 = reactantRect.lowerRight.x + 1.5;
          marvinReaction->arrow.y1 =(reactantRect.getCenter().y + productRect.getCenter().y)/2.0;
          marvinReaction->arrow.y2 = marvinReaction->arrow.y1;
        }

      }

      return marvinReaction;
    }
    catch(const std::exception& e)
    {
      delete marvinReaction;
      throw;
    }

  }
}
