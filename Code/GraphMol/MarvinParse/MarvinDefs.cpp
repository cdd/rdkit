//
//  Copyright (C) 2002-2022 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "MarvinDefs.h"
// #include <GraphMol/RDKitBase.h>
// #include <GraphMol/ChemReactions/Reaction.h>
// #include <RDGeneral/FileParseException.h>
// #include <RDGeneral/BadFileException.h>
// #include <RDGeneral/LocaleSwitcher.h>

// #include <string>
// #include <iostream>

namespace RDKit 
{
  std::string MarvinArrow::toString() const
  {
    std::ostringstream out;
    out << "<arrow type=\"" << type << "\" x1=\"" << x1<< "\" y1=\"" << y1 << "\" x2=\"" << x2<< "\" y2=\"" << y2 << "\"/>";
  
    return out.str();
  }
   
  std::string MarvinPlus::toString() const
  {
    std::ostringstream out;

    out << "<MReactionSign id=\"" << id << "\" toptions=\"NOROT\" fontScale=\"14.0\" halign=\"CENTER\" valign=\"CENTER\" autoSize=\"true\">"
      "<Field name=\"text\"><![CDATA[ {D font=SansSerif,size=18,bold}+ ]]></Field>"
      "<MPoint x=\"" << x1 << "\" y=\"" << y1 << "\"/>"
      "<MPoint x=\"" << x2 << "\" y=\"" << y1 << "\"/>"
      "<MPoint x=\"" << x2 << "\" y=\"" << y2 << "\"/>"
      "<MPoint x=\"" << x1 << "\" y=\"" << y2 << "/>"
      "</MReactionSign>";

    return out.str();
  }

  std::string MarvinCondition::toString() const
  {
    std::ostringstream out;

    out << "<MTextBox id=\"" << id << "\" toption=\"NOROT\" halign=\"" << halign << "\" valign=\"" << valign << "\" autoSize=\"true\"";
    if (fontScale > 0)
      out << " fontScale=\"" << fontScale << "\"";
  
    out << ">"
      "<Field name=\"text\">" << text << "</Field>"
      "<MPoint x=\"" << x << "\" y=\"" << y << "\"/>"
      "<MPoint x=\"" << x << "\" y=\"" << y << "\"/>"
      "<MPoint x=\"" << x << "\" y=\"" << y << "\"/>"
      "<MPoint x=\"" << x << "\" y=\"" << y << "\"/>"
      "</MTextBox>";



    return out.str();
  }

  std::string MarvinAttachmentPoint::toString() const
  {
    std::ostringstream out;

    out << "<attachmentPoint atom=\"" << atom << "\" order=\"" << order << "\" bond=\"" << bond << "\"/>";

    return out.str();
  }

  MarvinAtom::MarvinAtom()
    : x2(DBL_MAX)
    , y2(DBL_MAX)
    , formalCharge(0)
    , mrvValence(-1)
    , hydrogenCount(-1)
    , mrvMap(0)
    , rgroupRef (-1)   // indicates that it was not specified

  {
  }

  bool MarvinAtom::operator==(const MarvinAtom& rhs) const
  {
    return this->id == rhs.id;
  }
  
  bool MarvinAtom::operator==(const MarvinAtom *rhs) const
  {
    return this->id == rhs->id;
  }

  bool MarvinAtom::isElement() const
  {
    return this->elementType != "R" && this->elementType != "*";
  }

  std::string MarvinAtom::toString() const
  {
    // <atom id="a7" elementType="C" x2="15.225" y2="-8.3972" sgroupAttachmentPoint="1"/>

    std::ostringstream out;
    out << "<atom id=\"" << id << "\" elementType=\"" << elementType << "\"";
    
    if (x2 != DBL_MAX && y2 != DBL_MAX)
       out << " x2=\"" << x2 << "\" y2=\"" << y2 << "\"";

    if (formalCharge !=0)
      out << " formalCharge=\"" << formalCharge << "\"";

    if (radical != "")
      out << " radical=\"" << radical << "\"";

    if (isElement() &&  isotope != 0)
      out << " isotope=\"" << isotope << "\"";

    if (mrvAlias != "")
      out << " mrvAlias=\"" << mrvAlias << "\"";

    if (mrvValence >= 0)
      out << " mrvValence=\"" << mrvValence << "\"";

    if (hydrogenCount > 0)
      out << " hydrogenCount=\"" << hydrogenCount << "\"";

    if (mrvStereoGroup != "")
      out << " mrvStereoGroup=\"" << mrvStereoGroup << "\"";

    if (mrvMap != 0)
      out << " mrvMap=\"" << mrvMap << "\"";

    if (sgroupRef != "")
      out << " sgroupRef=\"" << sgroupRef << "\"";

    if (rgroupRef >= 0) 
      out << " rgroupRef=\"" << rgroupRef << "\"";

    if (sgroupAttachmentPoint != "")
      out << " sgroupAttachmentPoint=\"" << sgroupAttachmentPoint << "\"";

    out << "/>";

    return out.str();
  }     

  const std::string MarvinBond::getBondType() const
  {
    std::string tempQueryType = boost::algorithm::to_upper_copy(queryType);
    std::string tempOrder = boost::algorithm::to_upper_copy(order);
    std::string tempConvention = boost::algorithm::to_upper_copy(convention);

    if (tempQueryType != "")
    {
      if (tempQueryType == "SD" || tempQueryType == "SA"|| tempQueryType == "DA" || tempQueryType == "ANY")
        return tempQueryType;
      else
      {
        std::ostringstream err;
        err      << "unrecognized query bond type " << queryType << " in MRV File ";
        throw FileParseException(err.str());
      } 
    }
    else if (tempConvention != "")// if no query type, check for convention 
    {
      if (tempConvention == "CXN:COORD")
        return "DATIVE";
      else
      {
          std::ostringstream err;
          err      << "unrecognized convention " << convention << " in MRV File ";
          throw FileParseException(err.str());
      } 
    }      
    else if (tempOrder != "") // if no query type not conventtion,  so check for order
    {
      if (tempOrder == "1" || tempOrder == "2" || tempOrder == "3" || tempOrder == "A")
      {
        return tempOrder;      
      }
      else
      {
          std::ostringstream err;
          err      << "unrecognized bond type " << order << " in MRV File ";
          throw FileParseException(err.str());
      } 
    }
    else
    {
        std::ostringstream err;
        err      << "bond must have one of:  order, queryType, or convention in MRV File ";
        throw FileParseException(err.str());
    }       
}

  bool MarvinBond::isEqual(const MarvinAtom& other) const
  {
    return this->id == other.id;
  }

  bool MarvinBond::operator==(const MarvinAtom& rhs) const
  {
    return this->isEqual(rhs);
  }

  std::string MarvinBond::toString() const
  {
  // <bond id="b8" atomRefs2="a1 a7" order="1">

    std::ostringstream out;
    
    out << "<bond id=\"" << id << "\" atomRefs2=\"" << atomRefs2[0] << " " << atomRefs2[1] << "\"";
  
    if (order != "")
      out << " order=\"" << order << "\"";

    if (queryType != "")
      out << " queryType=\"" << queryType << "\"";

    if (convention != "")
      out << " convention=\"" << convention << "\"";
    
    if (bondStereo != "")
      out << "><bondStereo>" << bondStereo << "</bondStereo></bond>";
    else
      out << "/>";          
    return out.str();
  }

  MarvinMolBase::~MarvinMolBase()
  {   
  }
  
  int MarvinMolBase::getAtomIndex(std::string id)
  {
    auto atomIter = find_if(atoms.begin(), atoms.end(), [id](const MarvinAtom *arg) { return arg->id == id; });
    if (atomIter != atoms.end())
      return atomIter - atoms.begin();
    else 
      return -1;
  }

  int MarvinMolBase::getBondIndex(std::string id)
  {
    auto bondIter = find_if(bonds.begin(), bonds.end(), [id](const MarvinBond *arg) { return arg->id == id; });
    if (bondIter != bonds.end())
      return bondIter - bonds.begin();
    else 
      return -1;
  }

  const std::vector<std::string> MarvinMolBase::getBondList() const
  {
    std::vector<std::string> bondList;
    for (auto bond : bonds)
      bondList.push_back(bond->id);

    return bondList;
  }
  
  const std::vector<std::string> MarvinMolBase::getAtomList() const
  {
    std::vector<std::string> atomList;
    for (auto atom : atoms)
      atomList.push_back(atom->id);

    return atomList;
  }

  bool MarvinMolBase::hasCoords() const
  {
    for (auto atom: atoms)
      if (atom->x2 == DBL_MAX || atom->y2 == DBL_MAX)
        return false;

    return true;
  }

  void MarvinMolBase::removeCoords()
  {
    for (auto atom: atoms)
    {
      atom->x2 = DBL_MAX;
      atom->y2 = DBL_MAX;
    }
  }

  int MarvinMolBase::getExplicitValence(const MarvinAtom &marvinAtom) const
  {
    unsigned int resTimes10=0;   // calculated as 10 * the actual value so we can use int match, and have 1.5 order bonds

    for (auto bondPtr : bonds)
    {
      if (bondPtr->atomRefs2[0] != marvinAtom.id && bondPtr->atomRefs2[1] != marvinAtom.id)
        continue;  // this bond is NOT to the atom

      std::string tempConvention = boost::algorithm::to_upper_copy(bondPtr->convention);
      std::string marvinBondType = bondPtr->getBondType();

      if (marvinBondType == "SD" || marvinBondType == "SA" || marvinBondType == "DA")
      {
        resTimes10 += 15;   // really 1.5 order bond
      }
      else if (marvinBondType == "ANY")
      {
        resTimes10 += 10;   // no good answer for Any bonds - treat as a single bond
      }
      else if (marvinBondType == "DATIVE")
      {
        // if (bondPtr->atomRefs2[1] == marvinAtom.id) //second atom of dative bond count as 1 (first atom is zero)
        //   resTimes10 += 10;  // really 1 order bond
      }
      else if (marvinBondType == "1")
      {
          resTimes10 += 10; // really 1 order bond
      }
      else if (marvinBondType == "2")
      {
         resTimes10 += 20;   // really 2 order bond 
      }
      else if (marvinBondType == "3")
      {
        resTimes10 += 30;  // really 3 order bond
      }
      else if (marvinBondType == "A")
      {
         resTimes10 += 15;   // really 1.5 order bond
      }
    }

    return resTimes10/10;
  }

 
  std::string MarvinSruSgroup::toString() const
  {
    std::ostringstream out;

    out << "<molecule molID=\"" << molID << "\" id=\"" << id << "\" role=\"SruSgroup\" atomRefs=\"" << boost::algorithm::join(getAtomList()," ") << "\" title=\"" << title 
    <<"\" connect=\"" << connect << "\" correspondence=\"" << correspondence << "\" bondList=\"" << boost::algorithm::join(getBondList()," ") << "\"/>";

    return out.str();
  }

  std::string MarvinSruSgroup::role() const
  {
    return std::string("SruSgroup");
  } 

  bool MarvinSruSgroup::hasAtomBondBlocks() const
  {
    return false;
  } 
  
  
  std::string MarvinDataSgroup::toString() const
  {
    std::ostringstream out;

    out << "<molecule molID=\"" << molID << "\" id=\"" << id << "\" role=\"DataSgroup\" atomRefs=\"" << boost::algorithm::join(getAtomList()," ") << "\" context=\"" << context 
    <<"\" fieldName=\"" << fieldName << "\" placement=\"" << placement << "\" unitsDisplayed=\"" << unitsDisplayed << "\" fieldData=\"" << fieldData;
  
    if (units != "")
      out << "\" units=\"" << units;
    
    if (queryType != "" && queryOp != "")
      out << "\" queryType=\"" << queryType << "\" queryOp=\"" << boost::property_tree::xml_parser::encode_char_entities(queryOp);

    out << "\" x=\"" << x << "\" y=\"" << y << "\"/>";

    return out.str();
  }

  std::string MarvinDataSgroup::role() const
  {
    return std::string("DataSgroup");
  } 

  bool MarvinDataSgroup::hasAtomBondBlocks() const
  {
    return false;
  } 
  

  std::string MarvinMultipleSgroup::role() const
  {
    return std::string("MultipleSgroup");
  } 

  bool MarvinMultipleSgroup::hasAtomBondBlocks() const
  {
    return false;
  } 

  std::string MarvinMultipleSgroup::toString() const
  {
    std::ostringstream out;

    out << "<molecule molID=\"" << molID << "\" id=\"" << id << "\" role=\"MultipleSgroup\" atomRefs=\"" << boost::algorithm::join(getAtomList()," ") << "\" title=\"" << title 
    << "\"/>";

    return out.str();
  }

  MarvinSuperatomSgroupExpanded::~MarvinSuperatomSgroupExpanded()
  {
    for ( std::vector<MarvinAttachmentPoint *>::iterator it = attachmentPoints.begin(); it != attachmentPoints.end(); ++it)
      delete(*it);
  
  }

  std::string MarvinSuperatomSgroupExpanded::toString() const
  {
    std::ostringstream out;

    out << "<molecule molID=\"" << molID << "\" id=\"" << id << "\" role=\"SuperatomSgroup\" atomRefs=\"" << boost::algorithm::join(getAtomList()," ") << "\" title=\"" << title 
    << "\"/>";

    return out.str();
  }

  std::string MarvinSuperatomSgroupExpanded::role() const
  {
    return std::string("SuperatomSgroupExpanded");
  } 

  bool MarvinSuperatomSgroupExpanded::hasAtomBondBlocks() const
  {
    return false;
  } 

  MarvinSuperatomSgroup::~MarvinSuperatomSgroup()
  {
    for ( std::vector<MarvinAttachmentPoint *>::iterator it = attachmentPoints.begin(); it != attachmentPoints.end(); ++it)
      delete(*it);
  
    for ( std::vector<MarvinAtom *>::iterator it = atoms.begin(); it != atoms.end(); ++it)
    {
      delete(*it);
      *it = NULL;
    }
    for ( std::vector<MarvinBond *>::iterator it = bonds.begin(); it != bonds.end(); ++it)
    {
      delete(*it);
      *it = NULL;
    }
  }

  std::string MarvinSuperatomSgroup::role() const
  {
    return std::string("SuperatomSgroup");
  } 

  bool MarvinSuperatomSgroup::hasAtomBondBlocks() const
  {
    return true;
  } 

  std::string MarvinSuperatomSgroup::toString() const
  {
    std::ostringstream out;

    out << "<molecule molID=\"" << molID << "\" id=\"" << id << "\" role=\"SuperatomSgroup\" title=\"" << title << "\">";

    out <<"<atomArray>";
    for (std::vector<MarvinAtom *>::const_iterator it = atoms.begin();  it != atoms.end(); ++it)
      out << (*it)->toString();
    out << "</atomArray>";

    out <<"<bondArray>";
    for (std::vector<MarvinBond *>::const_iterator it = bonds.begin();  it != bonds.end(); ++it)
      out << (*it)->toString();
    out << "</bondArray>";

    if (attachmentPoints.size() > 0)
    {
      out <<"<AttachmentPointArray>";
      for (std::vector<MarvinAttachmentPoint *>::const_iterator it = attachmentPoints.begin();  it != attachmentPoints.end(); ++it)
        out << (*it)->toString();
      out << "</AttachmentPointArray>";
    }
    out << "</molecule>";
    
    return out.str();
  }   

  MarvinMol::~MarvinMol()
  {
    for ( std::vector<MarvinAtom *>::iterator it = atoms.begin(); it != atoms.end(); ++it)
    {
      delete(*it);
      *it = NULL;
    }
    for ( std::vector<MarvinBond *>::iterator it = bonds.begin(); it != bonds.end(); ++it)
    {
      delete(*it);
      *it = NULL;
    }
    for (std::vector<MarvinSuperatomSgroup *>::iterator it = superatomSgroups.begin(); it != superatomSgroups.end(); ++it)
      delete(*it);
    for (std::vector<MarvinSruSgroup *>::iterator it = sruSgroups.begin(); it != sruSgroups.end(); ++it)
      delete(*it);
    for (std::vector<MarvinMultipleSgroup *>::iterator it = multipleSgroups.begin(); it != multipleSgroups.end(); ++it)
      delete(*it);
    for (std::vector<MarvinSuperatomSgroupExpanded *>::iterator it = superatomSgroupsExpanded.begin(); it != superatomSgroupsExpanded.end(); ++it)
      delete(*it);

    for (std::vector<MarvinDataSgroup *>::iterator it =dataSgroups.begin(); it != dataSgroups.end(); ++it)
      delete(*it);

    for (std::vector<MarvinSuperInfo *>::iterator it = superInfos.begin(); it != superInfos.end(); ++it)
      delete(*it);
  }

  std::string MarvinMol::role() const
  {
    return std::string("Parent");
  } 

  bool MarvinMol::hasAtomBondBlocks() const
  {
    return true;
  } 

  bool MarvinMol::atomRefInAtoms(MarvinAtom *a, std::string b )
  { 
    return a->id == b; 
  }

  bool MarvinMol::bondRefInBonds(MarvinBond *a, std::string b )
  { 
    return a->id == b; 
  }

  void MarvinMol::cleanUpNumbering(int &molCount    // this is the starting mol count, and receives the ending mol count - THis is used when MarvinMol->convertToSuperAtaoms is called multiple times from a RXN
  , int &atomCount  // starting and ending atom count
  , int &bondCount  // starting and ending bond count
  , int &sgCount)  // starting and ending sg  count)
  {
    // clean and renumber the sgroups

    std::map<std::string,std::string> sgMap;
    for (MarvinSuperatomSgroup *marvinSuperatomSgroup : this->superatomSgroups)
    {
      std::string newId = "sg" + std::to_string(++sgCount);
      sgMap[marvinSuperatomSgroup->id] = newId;
      marvinSuperatomSgroup->id = newId;
    }

    for (MarvinSuperatomSgroupExpanded *marvinSuperatomSgroupExpanded : this->superatomSgroupsExpanded)         
    {
      std::string newId = "sg" + std::to_string(++sgCount);
      sgMap[marvinSuperatomSgroupExpanded->id] = newId;
      marvinSuperatomSgroupExpanded->id = newId;
    }

    for (MarvinMultipleSgroup *marvinMultipleSgroup : this->multipleSgroups)         
    {
      std::string newId = "sg" + std::to_string(++sgCount);
      sgMap[marvinMultipleSgroup->id] = newId;
      marvinMultipleSgroup->id = newId;
    }

    for (MarvinSruSgroup *marvinSruSGgroup : this->sruSgroups)         
    {
      std::string newId = "sg" + std::to_string(++sgCount);
      sgMap[marvinSruSGgroup->id] = newId;
      marvinSruSGgroup->id = newId;
    }

    for (MarvinDataSgroup *marvinDataSGgroup : this->dataSgroups)         
    {
      std::string newId = "sg" + std::to_string(++sgCount);
      sgMap[marvinDataSGgroup->id] = newId;
      marvinDataSGgroup->id = newId;
    }


    // clean up the mol ids, the atomIds and bondIds.  make  map of the old to new atom ids and bond ids

    this->molID = "m" + std::to_string(++molCount);
    
    std::map<std::string,std::string> atomMap;
    std::map<std::string,std::string> bondMap;
    
    for (auto atomPtr : this->atoms)
    {
      std::string newId = "a" + std::to_string(++atomCount);
      atomMap[atomPtr->id] =newId;
      atomPtr->id = newId;

      // fix the sgroupRef 
      if (atomPtr->sgroupRef != "")
        atomPtr->sgroupRef = sgMap[atomPtr->sgroupRef];
    }

    for (auto bondPtr : this->bonds)
    {
      std::string newId = "b" + std::to_string(++bondCount);
      bondMap[bondPtr->id] =newId;
      bondPtr->id = newId;

      // fix the bond's references to atoms

      bondPtr->atomRefs2[0] = atomMap[bondPtr->atomRefs2[0]];
      bondPtr->atomRefs2[1] = atomMap[bondPtr->atomRefs2[1]];
    }

    //Now the atoms and bonds for the superatoms

    for (MarvinSuperatomSgroup *marvinSuperatomSgroup : this->superatomSgroups)
    {
      marvinSuperatomSgroup->molID = "m" + std::to_string(++molCount);
      for (auto atomPtr : marvinSuperatomSgroup->atoms)
      {
        std::string newId = "a" + std::to_string(++atomCount);
        atomMap[atomPtr->id] =newId;
        atomPtr->id = newId;
      }

      for (auto bondPtr : marvinSuperatomSgroup->bonds)
      {
        std::string newId = "b" + std::to_string(++bondCount);
        bondMap[bondPtr->id] =newId;
        bondPtr->id = newId;

        // fix the bond's references to atoms

        bondPtr->atomRefs2[0] = atomMap[bondPtr->atomRefs2[0]];
        bondPtr->atomRefs2[1] = atomMap[bondPtr->atomRefs2[1]];
      }

      for (auto attachmentPoint : marvinSuperatomSgroup->attachmentPoints)
      {
        attachmentPoint->atom = atomMap[attachmentPoint->atom];
        attachmentPoint->bond = bondMap[attachmentPoint->bond];  // bond is actually in the parent
      }
    }
  }

  void MarvinMol::convertFromSuperAtoms()
  {
  // the mol-style super atoms are significanyly different than the Marvin super atoms.  The mol style has all atoms and bonds in the main mol, and parameter lines that
  // indicate the name of the super atom and the atom indices affected.
  //
  // There are 3 general forms in MRV:
  //    1) SUP atoms   
  //    2) MultipleSgroup
  //    3) SruGroup
  //
  //    The SUP form can appear in conracted or expanded form in MRV blocks:
  //   
  // In the contracted form, the MRV block can have one or more dummy atoms with the super atom name in the main molecule, and a separate sub-mol with the atoms and bonds of the superatom.  
  // It also can have one or more separate records 
  // called attachmentPoints that specifiy which atom(s) in the super atom sub-mol that replaces the dummy atom(s), and also a bond pointer (in the parent mol).
  // 
  // In the expanded form, all of the atoms are in the parent molecule, and the sub-mol only refers to them.  The attachement points refer to the atoms in the parent mol.
  // The MultipleSgroup and SruGroup seem very much alike.   In both, all atoms are in the parent mol, and the sub-mol refers to a group of atoms in that parent.
  // The SruGroup specifies the name and also the connection informat (head-to-tail, head-head, and unspecified).  The Multiple S-group specifies the report count for the group
  //
  // This routine deals with only the contracted form of the Supergroups, and copies the atoms and bonds from the sub=mol to the parent mol, and deleted the dummy atom form the parent mol.  It also saves the infor needed to make a mol-file type
  // superatom in the MarvinSuperInfo array.


    for(std::vector<MarvinSuperatomSgroup *>::iterator subMolIter =  superatomSgroups.begin() ; subMolIter != superatomSgroups.end() ; ++subMolIter)
    {
      //save the name of the superatom
      
      auto marvinSuperInfo = new MarvinSuperInfo();
      this->superInfos.push_back(marvinSuperInfo);

      marvinSuperInfo->title = (*subMolIter)->title;

      //  remove and delete the dummy atoms from the parent.

      bool coordExist = (*subMolIter)->hasCoords() && hasCoords();
      auto dummyAtomIter = find_if(atoms.begin(), atoms.end(), [subMolIter](const MarvinAtom *arg) { 
                        return arg->sgroupRef == (*subMolIter)->id; });
      if (dummyAtomIter == atoms.end())
        throw FileParseException("No contracted atom found for a superatomSgroup");
    
      // get a list of the bonds that will need to be fixed after the copy of atoms and bonds

      std::vector<MarvinAttachmentPoint *> orphanedBonds;
      RDGeom::Point3D centerOfGroup;
      for (auto orphanedBond = bonds.begin(); orphanedBond != bonds.end() ; ++orphanedBond)
      {
        std::string attachedAtomId = "";
        if ((*orphanedBond)->atomRefs2[0] == (*dummyAtomIter)->id)
          attachedAtomId = (*orphanedBond)->atomRefs2[1];
        else if ((*orphanedBond)->atomRefs2[1] == (*dummyAtomIter)->id)
          attachedAtomId = (*orphanedBond)->atomRefs2[0];
       
        if (attachedAtomId != "")
        {
          auto attachmentPoint = find_if((*subMolIter)->attachmentPoints.begin(), (*subMolIter)->attachmentPoints.end(), [orphanedBond](const MarvinAttachmentPoint *arg) { 
                        return arg->bond == (*orphanedBond)->id; });
          if (attachmentPoint == (*subMolIter)->attachmentPoints.end())
             throw FileParseException("No attachment point found for bond to the conedensed atom in a superatomSgroup");
          orphanedBonds.push_back(*attachmentPoint);
          
          auto subMolAtomIndex = (*subMolIter)->getAtomIndex((*attachmentPoint)->atom);
          if (subMolAtomIndex < 0)
            throw FileParseException("Attachment atom not found in the superatomsGroup");
          MarvinAtom *attachedAtom = (*subMolIter)->atoms[subMolAtomIndex];
          if (coordExist)
          {
            centerOfGroup.x += attachedAtom->x2;
            centerOfGroup.y += attachedAtom->y2;
          }
        }
      }

      if (coordExist && orphanedBonds.size() > 0)
      {
        centerOfGroup.x /= orphanedBonds.size();
        centerOfGroup.y /= orphanedBonds.size();
        RDGeom::Point3D offset;
        offset.x = (*dummyAtomIter)->x2 - centerOfGroup.x;
        offset.y = (*dummyAtomIter)->y2 - centerOfGroup.y;

        for (auto atom : (*subMolIter)->atoms)
        {
          atom->x2 += offset.x;
          atom->y2 += offset.y;
        }
      }

      delete *dummyAtomIter;  // get rid of the MolAtom
      atoms.erase(dummyAtomIter);   // get rid of the atoms pointer to the old dummy atom


      // add the atoms and bonds from the super group to the parent

      for (std::vector<MarvinAtom *>::iterator subAtomIter = (*subMolIter)->atoms.begin() ; subAtomIter != (*subMolIter)->atoms.end() ; ++subAtomIter)
      {
        atoms.push_back( *subAtomIter);
        marvinSuperInfo->atoms.push_back((*subAtomIter)->id);

        // remove the sgroupRef from the atom if it has one
        (*subAtomIter)->sgroupAttachmentPoint = "";

        *subAtomIter = NULL;  // so that it is not in two places  - prevents double deleting
      }
      for (std::vector<MarvinBond *>::iterator subBondIter = (*subMolIter)->bonds.begin()  ; subBondIter != (*subMolIter)->bonds.end() ; ++subBondIter)
      {
        bonds.push_back( *subBondIter);
        *subBondIter = NULL;  // so that it is not in two places  - prevents double deleting
      }

      // process the attachment points - fix the bond that was made wrong by deleting the dummy atom(s)

      for (std::vector<MarvinAttachmentPoint *>::iterator attachIter = (*subMolIter) -> attachmentPoints.begin() ;  attachIter != (*subMolIter) -> attachmentPoints.end(); ++attachIter)
      {
        // find the bond in the parent
        
        auto bondIter = find_if(bonds.begin(), bonds.end(), [attachIter](const MarvinBond *arg) { 
                        return arg->id == (*attachIter)->bond; });
        if (bondIter == bonds.end())
          throw FileParseException("Bond specification for an AttachmentPoint definition was not found in the bond array in MRV file");
        
        marvinSuperInfo->bonds.push_back((*attachIter)->bond);

        // one of the two atoms in the bond is NOT in the mol - we deleted the dummy atom.

        int atomIndex;
        for (atomIndex = 0 ; atomIndex < 2 ; ++atomIndex)
        {
          if (!boost::algorithm::contains(atoms, std::vector<std::string>{(*bondIter)->atomRefs2[atomIndex]}, atomRefInAtoms ))
          {
            (*bondIter)->atomRefs2[atomIndex] = (*attachIter)->atom;   // the attach atom
            break;
          }
        }
        if (atomIndex == 2)  // not found?
        {
          std::ostringstream err;
          err << "Bond " << (*attachIter)->bond.c_str() << " from attachment point did not have a missing atom in MRV file";
          throw FileParseException(err.str());
        }

        delete *attachIter;         
      }
      
      (*subMolIter)->attachmentPoints.clear();
      (*subMolIter)->atoms.clear();
      (*subMolIter)->bonds.clear();
      delete *subMolIter;
    }
    this->superatomSgroups.clear();
  }

  void MarvinMol::convertToSuperAtoms()
  {
    // the mol-style super atoms are significatnly different than the Marvin super atoms.  The mol style has all atoms and bonds in the main mol, and parameter lines that
    // indicate the name of the super atom and the atom indices affected.
    //
    // The Marvin as a dummy atom with the super atom name in the main molecule, and a separate sub-mol with the atoms and bonds of the superatom.  It also has a separate record 
    // called attachmentPoint that atom in the super atom sub=mol that is replaces the dummy atom, and also a bond pointer and bond order in case the bond order changes.
    //
    //  Takes information from a MarvinSuperInfos and converts the Marvin structure  to have the super atoms defined

    int superGroupsAdded = 0;    
    for (auto marvinSuperInfo : this->superInfos)
    {
      // make a new sub mol

      auto marvinSuperatomSgroup = new MarvinSuperatomSgroup();
      this->superatomSgroups.push_back(marvinSuperatomSgroup);
      superGroupsAdded++;
      std::string newAtomName =  "NA" + std::to_string( superGroupsAdded);

      marvinSuperatomSgroup->molID = "NM" + std::to_string( superGroupsAdded);;
      marvinSuperatomSgroup->title = marvinSuperInfo->title;
      marvinSuperatomSgroup->id = "nsg" + std::to_string(superGroupsAdded);  // n in nsg ensures no colllision with other sgs already in the mol 

      bool coordsExist = hasCoords();  // we have to check before we add the dummy atom  - it will not have coords (yet)

      // add the dummy atom into the parent

      MarvinAtom *dummyParentAtom = new MarvinAtom();
      this->atoms.push_back(dummyParentAtom);
      dummyParentAtom->elementType = "R";
      dummyParentAtom->id = newAtomName;
      dummyParentAtom->sgroupRef = marvinSuperatomSgroup->id;

      for (auto atom : marvinSuperInfo->atoms)
      {
        // Find the atom in the main structure

        int index = this->getAtomIndex(atom);
        MarvinAtom *atomToMove = this->atoms[index];
        marvinSuperatomSgroup->atoms.push_back(atomToMove);

        this->atoms.erase(this->atoms.begin() + index);  
      }

      // move the bonds of the group
      
      int attachmentPointsAdded = 0;
      MarvinAtom *atomPtr = NULL;
      RDGeom::Point3D centerOfGroup;

      for (auto bond : this->bonds)
      {
        bool atom1InGroup = boost::algorithm::contains(marvinSuperatomSgroup->atoms, std::vector<std::string>{(bond)->atomRefs2[0]}, atomRefInAtoms );
        bool atom2InGroup = boost::algorithm::contains(marvinSuperatomSgroup->atoms, std::vector<std::string>{(bond)->atomRefs2[1]}, atomRefInAtoms );
        if (atom1InGroup && atom2InGroup )  // both are in, so move the bond
          marvinSuperatomSgroup->bonds.push_back(bond);

        else if(atom1InGroup || atom2InGroup ) // one is in so this is a connection point
        {
          // fix the bonds to the dummy atom and add an attachment point

          if (atom1InGroup)
          {
            atomPtr = marvinSuperatomSgroup->atoms[marvinSuperatomSgroup->getAtomIndex((bond)->atomRefs2[0])];       
            (bond)->atomRefs2[0] = newAtomName;

          }
          else 
          {
            atomPtr = marvinSuperatomSgroup->atoms[marvinSuperatomSgroup->getAtomIndex((bond)->atomRefs2[1])];       
            (bond)->atomRefs2[1] = newAtomName;
          }
              
          atomPtr->sgroupAttachmentPoint = std::to_string(++attachmentPointsAdded);

          // add an attachentPoint structure

          auto marvinAttachmentPoint = new MarvinAttachmentPoint();
          marvinSuperatomSgroup->attachmentPoints.push_back(marvinAttachmentPoint);
          marvinAttachmentPoint->atom = atomPtr->id;
          marvinAttachmentPoint->bond = bond->id;;
          marvinAttachmentPoint->order = std::to_string(attachmentPointsAdded); 

          if (coordsExist)
          {
            centerOfGroup.x += atomPtr->x2;
            centerOfGroup.y += atomPtr->y2;
          }       
        }
      }

      // now remove the bonds that were moved to the superGroup from the parent
      
      for (auto bond : marvinSuperatomSgroup->bonds)
      {
        int index = this->getBondIndex(bond->id);
        this->bonds.erase(this->bonds.begin() + index);
      }

    
      if(coordsExist)   
      {
        if (marvinSuperatomSgroup->attachmentPoints.size() > 0) // Any attachment points?  if so we use the center of the attached atoms in the supergroup
        {
          // put the new dummy atom at the center of the removed group

          dummyParentAtom->x2 = centerOfGroup.x / marvinSuperatomSgroup->attachmentPoints.size();
          dummyParentAtom->y2 = centerOfGroup.y / marvinSuperatomSgroup->attachmentPoints.size();
        }
        else
        {
          // No attachments to the supergroup - (probably all atoms in the mol are in the supergroup) - doesn't matter where we put the dummy atom
          dummyParentAtom->x2 = 0.0;
          dummyParentAtom->y2 =  0.0;
        }
      }
    }
    
    for (auto marvinSuperInfo : this->superInfos)
      delete marvinSuperInfo;

    this->superInfos.clear();
  }

  std::string MarvinMol::toString() const
  {
    std::ostringstream out;

    out << "<molecule molID=\"" << molID << "\">";

    out <<"<atomArray>";
    for (std::vector<MarvinAtom *>::const_iterator it = atoms.begin();  it != atoms.end(); ++it)
      out << (*it)->toString();
    out <<"</atomArray>";

    out <<"<bondArray>";
    for (std::vector<MarvinBond *>::const_iterator it = bonds.begin();  it != bonds.end(); ++it)
      out << (*it)->toString();
    out << "</bondArray>";

    for (std::vector<MarvinSuperatomSgroup *>::const_iterator it = superatomSgroups.begin();  it != superatomSgroups.end(); ++it)
      out << (*it)->toString();
    for (std::vector<MarvinSuperatomSgroupExpanded *>::const_iterator it = superatomSgroupsExpanded.begin();  it != superatomSgroupsExpanded.end(); ++it)
      out << (*it)->toString();
    for (std::vector<MarvinMultipleSgroup *>::const_iterator it = multipleSgroups.begin();  it != multipleSgroups.end(); ++it)
      out << (*it)->toString();
    for (std::vector<MarvinSruSgroup *>::const_iterator it = sruSgroups.begin();  it != sruSgroups.end(); ++it)
      out << (*it)->toString();
    for (std::vector<MarvinDataSgroup *>::const_iterator it = dataSgroups.begin();  it != dataSgroups.end(); ++it)
      out << (*it)->toString();
    
    out <<"</molecule>";
    
    return out.str();
  }   

  std::string MarvinMol::generateMolString()
  {
    std::ostringstream out;

    out << "<cml xmlns=\"http://www.chemaxon.com\"  xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
      "xsi:schemaLocation=\"http://www.chemaxon.com http://www.chemaxon.com/marvin/schema/mrvSchema_20_20_0.xsd\">"
      "<MDocument><MChemicalStruct>";

    out << toString();
    
    out << "</MChemicalStruct></MDocument></cml>";
    return out.str();
  }

  MarvinReaction::~MarvinReaction()
  {
    for (std::vector<MarvinMol *>::iterator it = reactants.begin(); it != reactants.end(); ++it)
      delete(*it);      
    for (std::vector<MarvinMol *>::iterator it = agents.begin(); it != agents.end(); ++it)
      delete(*it);      
    for (std::vector<MarvinMol *>::iterator it = products.begin(); it != products.end(); ++it)
      delete(*it);
    for (std::vector<MarvinPlus *>::iterator it = pluses.begin(); it != pluses.end(); ++it)
      delete(*it);
    for (std::vector<MarvinCondition *>::iterator it = conditions.begin(); it != conditions.end(); ++it)
      delete(*it);
  }

  void MarvinReaction::convertFromSuperAtoms()
  {
    //This routine converts all the mols in the rxn to be ready for conversion to RDKIT mols

    for (std::vector<MarvinMol *>::iterator molIter = reactants.begin() ; molIter != reactants.end() ; ++molIter)
      (*molIter)->convertFromSuperAtoms();
    for (std::vector<MarvinMol *>::iterator molIter = agents.begin() ; molIter != agents.end() ; ++molIter)
      (*molIter)->convertFromSuperAtoms();
    for (std::vector<MarvinMol *>::iterator molIter = products.begin() ; molIter != products.end() ; ++molIter)
      (*molIter)->convertFromSuperAtoms();
  }

  void MarvinReaction::convertToSuperAtoms()
  {
    //This routine converts all the mols in the rxn to be ready for conversion to RDKIT mols
    int molCount=0, atomCount=0, bondCount=0, sgCount=0;

    for (std::vector<MarvinMol *>::iterator molIter = reactants.begin() ; molIter != reactants.end() ; ++molIter)
    {
      (*molIter)->convertToSuperAtoms();
      (*molIter)->cleanUpNumbering(molCount,atomCount, bondCount, sgCount);
    }
    for (std::vector<MarvinMol *>::iterator molIter = agents.begin() ; molIter != agents.end() ; ++molIter)
    {
      (*molIter)->convertToSuperAtoms();
      (*molIter)->cleanUpNumbering(molCount,atomCount, bondCount, sgCount);
    }       
    for (std::vector<MarvinMol *>::iterator molIter = products.begin() ; molIter != products.end() ; ++molIter)
    {
      (*molIter)->convertToSuperAtoms();
      (*molIter)->cleanUpNumbering(molCount,atomCount, bondCount, sgCount);
    }       
  }

  std::string MarvinReaction::toString()
  {
    std::ostringstream out;

    out << "<cml xmlns=\"http://www.chemaxon.com\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\""
        " xsi:schemaLocation=\"http://www.chemaxon.com http://www.chemaxon.com/marvin/schema/mrvSchema_20_20_0.xsd\">"
        "<MDocument><MChemicalStruct><reaction>";

    out <<"<reactantList>";
    for (std::vector<MarvinMol *>::const_iterator it = reactants.begin();  it != reactants.end(); ++it)
      out << (*it)->toString();
    out <<"</reactantList>";
    out <<"<agentList>";
    for (std::vector<MarvinMol *>::const_iterator it = agents.begin();  it != agents.end(); ++it)
      out << (*it)->toString();
    out <<"</agentList>";
    out <<"<productList>";
    for (std::vector<MarvinMol *>::const_iterator it = products.begin();  it != products.end(); ++it)
      out << (*it)->toString();
    out <<"</productList>";

    out << arrow.toString();

    out << "</reaction></MChemicalStruct>";

    for (std::vector<MarvinPlus *>::const_iterator it = pluses.begin();  it != pluses.end(); ++it)
      out << (*it)->toString();
    for (std::vector<MarvinCondition *>::const_iterator it = conditions.begin();  it != conditions.end(); ++it)
      out << (*it)->toString();

    
    out << "</MDocument></cml>";
    return out.str();
  }

  MarvinRectangle::MarvinRectangle(double left, double right, double top, double bottom)
  {
    upperLeft.x = left;
    upperLeft.y = top;
    lowerRight.x = right;
    lowerRight.y = bottom;
    centerIsStale = true;
  }

  MarvinRectangle::MarvinRectangle(const RDGeom::Point3D &upperLeftInit, const RDGeom::Point3D &lowerRightInit)
  {
    upperLeft = upperLeftInit;
    lowerRight = lowerRightInit;
    centerIsStale = true;
  }

  MarvinRectangle::MarvinRectangle(const std::vector<MarvinAtom *> atoms)
  {
    centerIsStale = true;

    if (atoms.size() == 0)
      return;
    upperLeft.x = DBL_MAX;
    upperLeft.y = -DBL_MAX;
    lowerRight.x = -DBL_MAX;
    lowerRight.y = DBL_MAX;


    for (auto atom : atoms)
    {
      if (atom->x2 < upperLeft.x)
        upperLeft.x = atom->x2;
      if (atom->x2 > lowerRight.x)
        lowerRight.x = atom->x2;

      
      if (atom->y2 > upperLeft.y)
        upperLeft.y = atom->y2;
      if (atom->y2 < lowerRight.y)
        lowerRight.y = atom->y2;
    }
  }

  void MarvinRectangle::extend(const MarvinRectangle &otherRectangle)
  {  
    if (otherRectangle.upperLeft.x < upperLeft.x)
      upperLeft.x = otherRectangle.upperLeft.x;
    if (otherRectangle.lowerRight.x > lowerRight.x)
      lowerRight.x = otherRectangle.lowerRight.x;

    
    if (otherRectangle.upperLeft.y > upperLeft.y)
      upperLeft.y = otherRectangle.upperLeft.y;
    if (otherRectangle.lowerRight.y < lowerRight.y)
      lowerRight.y = otherRectangle.lowerRight.y;
    
    centerIsStale = true;
  }

  RDGeom::Point3D &MarvinRectangle::getCenter()
  {
    if (centerIsStale)
    {
      center.x = (lowerRight.x + upperLeft.x)/2.0;
      center.y = (lowerRight.y + upperLeft.y)/2.0;
      centerIsStale = false;
    }
    return center;
  }

  bool MarvinRectangle::overlapsVertically(const MarvinRectangle &otherRectangle) const
  {
    if (otherRectangle.upperLeft.y < lowerRight.y || otherRectangle.lowerRight.y > upperLeft.y)
      return false;
    return true;
  }

  bool MarvinRectangle::overlapsVHorizontally(const MarvinRectangle &otherRectangle) const
  {
    if (otherRectangle.upperLeft.x > lowerRight.x || otherRectangle.lowerRight.x < upperLeft.x)
      return false;
    return true;
  }

  bool MarvinRectangle::compareRectanglesByX(MarvinRectangle &r1, MarvinRectangle &r2)
  {
    return (r1.getCenter().x < r2.getCenter().x);
  }

  bool MarvinRectangle::compareRectanglesByY(MarvinRectangle &r1, MarvinRectangle &r2)
  {
    return (r1.getCenter().y < r2.getCenter().y);
  }

  MarvinStereoGroup::MarvinStereoGroup(StereoGroupType grouptypeInit, int groupNumberInit )
  {
    groupType = grouptypeInit;
    groupNumber = groupNumberInit;
  }
}
