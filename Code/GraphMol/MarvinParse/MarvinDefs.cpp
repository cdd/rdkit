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
      "<Field name=\"text\"><![CDATA[ {D font=SansSerif,size=18,bold}+]]></Field>"
      "<MPoint x=\"" << x1 << "\" y=\"" << y1 << "\"/>"
      "<MPoint x=\"" << x2 << "\" y=\"" << y1 << "\"/>"
      "<MPoint x=\"" << x2 << "\" y=\"" << y2 << "\"/>"
      "<MPoint x=\"" << x1 << "\" y=\"" << y2 << "\"/>"
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

   MarvinAtom::MarvinAtom(const MarvinAtom &atomToCopy, std::string newId)
  :    id(newId)
    , elementType(atomToCopy.elementType)
    , x2(atomToCopy.x2)
    , y2(atomToCopy.y2)
    , formalCharge(atomToCopy.formalCharge)
    , radical(atomToCopy.radical)
    , isotope(atomToCopy.isotope)
    , mrvValence(atomToCopy.mrvValence)
    , hydrogenCount(atomToCopy.hydrogenCount)
    , mrvAlias(atomToCopy.mrvAlias)
    , mrvStereoGroup(atomToCopy.mrvStereoGroup)
    , mrvMap(atomToCopy.mrvMap)
    , sgroupRef(atomToCopy.sgroupRef)
    , sgroupAttachmentPoint(atomToCopy.sgroupAttachmentPoint)
    , rgroupRef(atomToCopy.rgroupRef)   // indicates that it was not specified
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

  MarvinBond::MarvinBond(const MarvinBond &bondToCopy, std::string newId, std::string atomRef1, std::string atomRef2)
    : id(newId)
    , order(bondToCopy.order)
    , bondStereo(bondToCopy.bondStereo)
    , queryType(bondToCopy.queryType)
    , convention(bondToCopy.convention)
  {
    atomRefs2[0] = atomRef1;
    atomRefs2[1] = atomRef2;
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

  std::string MarvinMulticenterSgroup::role() const
  {
    return std::string("MulticenterSgroup");
  } 

  bool MarvinMulticenterSgroup::hasAtomBondBlocks() const
  {
    return false;
  } 

  std::string MarvinMulticenterSgroup::toString() const
  {
     // <molecule molID="m2" id="sg1" role="MulticenterSgroup" atomRefs="a2 a6 a5 a4 a3" center="a18"/>
    std::ostringstream out;

    out << "<molecule molID=\"" << molID << "\" id=\"" << id << "\" role=\"MulticenterSgroup\" atomRefs=\"" << boost::algorithm::join(getAtomList()," ") << "\" center=\"" 
    << center->id << "\"/>";

    return out.str();
  }
  
  std::string MarvinGenericSgroup::role() const
  {
    return std::string("GenericSgroup");
  } 

  bool MarvinGenericSgroup::hasAtomBondBlocks() const
  {
    return false;
  } 

  std::string MarvinGenericSgroup::toString() const
  {
     // <molecule molID="m2" id="sg1" role="MulticenterSgroup" atomRefs="a2 a6 a5 a4 a3" center="a18"/>
    std::ostringstream out;

    out << "<molecule molID=\"" << molID << "\" id=\"" << id << "\" role=\"GenericSgroup\" atomRefs=\"" << boost::algorithm::join(getAtomList()," ")  << "\"/>";

    return out.str();
  }

  std::string MarvinMonomerSgroup::role() const
  {
    return std::string("MonomerSgroup");
  } 

  bool MarvinMonomerSgroup::hasAtomBondBlocks() const
  {
    return false;
  } 

  std::string MarvinMonomerSgroup::toString() const
  {
// <molecule id="sg1" role="MonomerSgroup" title="mon" charge="onAtoms" molID="m2" atomRefs="a2 a1 a3 a4">
      //     <MBracket type="SQUARE" orientation="DOUBLE">
      //         <MPoint x="-0.8726666666666667" y="1.078"></MPoint>
      //         <MPoint x="1.2833333333333334" y="1.078"></MPoint>
      //         <MPoint x="1.2833333333333334" y="-1.078"></MPoint>
      //         <MPoint x="-0.8726666666666667" y="-1.078"></MPoint>
      //     </MBracket>
      // </molecule> 
      
    std::ostringstream out;

    out << "<molecule molID=\"" << molID << "\" id=\"" << id << "\" role=\"MonomerSgroup\" atomRefs=\"" << boost::algorithm::join(getAtomList()," ")  
    << "\" title=\"" << title << "\" charge=\"" << charge << "\">" << bracket.toString() << "</molecule>";

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
    for (std::vector<MarvinDataSgroup *>::iterator it = dataSgroups.begin(); it != dataSgroups.end(); ++it)
      delete(*it);
    for (std::vector<MarvinMulticenterSgroup *>::iterator it = multicenterSgroups.begin(); it != multicenterSgroups.end(); ++it)
      delete(*it);
    for (std::vector<MarvinGenericSgroup *>::iterator it = genericSgroups.begin(); it != genericSgroups.end(); ++it)
      delete(*it);
    for (std::vector<MarvinMonomerSgroup *>::iterator it = monomerSgroups.begin(); it != monomerSgroups.end(); ++it)
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

    for (MarvinMulticenterSgroup *marvinMulticenterSgroup : this->multicenterSgroups)         
    {
      std::string newId = "sg" + std::to_string(++sgCount);
      sgMap[marvinMulticenterSgroup->id] = newId;
      marvinMulticenterSgroup->id = newId;
    }

    for (MarvinGenericSgroup *marvinGenericSgroup : this->genericSgroups)         
    {
      std::string newId = "sg" + std::to_string(++sgCount);
      sgMap[marvinGenericSgroup->id] = newId;
      marvinGenericSgroup->id = newId;
    }

    for (MarvinMonomerSgroup *marvinMonomerSgroup : this->monomerSgroups)         
    {
      std::string newId = "sg" + std::to_string(++sgCount);
      sgMap[marvinMonomerSgroup->id] = newId;
      marvinMonomerSgroup->id = newId;
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
  // the mol-style super atoms are significanTly different than the Marvin super atoms.  The mol style has all atoms and bonds in the main mol, and parameter lines that
  // indicate the name of the super atom and the atom indices affected.
  //
  // There are several general forms of super atom groups (Sgroups): SuperatomSgroup (SUP), MultipleSgroup (MUL), SruSGroup (SRU), DataSgroup (DAT)
  //      , MultiCenterSgroup (no RDKIT equiv),  GenericSgroup (GEN), and MonomerSgrop (MON).
  //  THis routine deals only with SuperatomSgrups. 
  //
  //  The SUP form can appear in contracted or expanded form in MRV blocks:
  //   
  // In the contracted form, the MRV block can have one or more dummy atoms with the super atom name in the main molecule, and a separate sub-mol with the atoms and bonds of the superatom.  
  // It also can have one or more separate records 
  // called attachmentPoints that specifiy which atom(s) in the super atom sub-mol that replaces the dummy atom(s), and also a bond pointer (in the parent mol).
  // 
  // In the expanded form, all of the atoms are in the parent molecule, and the sub-mol only refers to them.  The attachement points refer to the atoms in the parent mol.
  // The MultipleSgroup and SruGroup seem very much alike.   In both, all atoms are in the parent mol, and the sub-mol refers to a group of atoms in that parent.
  // The SruGroup specifies the name and also the connection informat (head-to-tail, head-head, and unspecified).  The Multiple S-group specifies the report count for the group
  //
  // This routine deals with only the contracted form of the Supergroups, and copies the atoms and bonds from the sub=mol to the parent mol, and deletes the dummy atom form the parent mol.
  //  It also saves the info needed to make a mol-file type superatom in the MarvinSuperInfo array.


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
      RDGeom::Point3D centerOfAttachmentPoints;
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
             throw FileParseException("No attachment point found for bond to the condensed atom in a superatomSgroup");
          orphanedBonds.push_back(*attachmentPoint);
          
          auto subMolAtomIndex = (*subMolIter)->getAtomIndex((*attachmentPoint)->atom);
          if (subMolAtomIndex < 0)
            throw FileParseException("Attachment atom not found in the superatomsGroup");
          MarvinAtom *attachedAtom = (*subMolIter)->atoms[subMolAtomIndex];
          if (coordExist)
          {
            centerOfAttachmentPoints.x += attachedAtom->x2;
            centerOfAttachmentPoints.y += attachedAtom->y2;
          }
        }
      }

      if (coordExist)
      {
        RDGeom::Point3D offset;

        if (orphanedBonds.size() > 0)
        {
          centerOfAttachmentPoints.x /= orphanedBonds.size();
          centerOfAttachmentPoints.y /= orphanedBonds.size();
          offset.x = (*dummyAtomIter)->x2 - centerOfAttachmentPoints.x;
          offset.y = (*dummyAtomIter)->y2 - centerOfAttachmentPoints.y;    
        }
        else // use the center of all atoms in the group
        {
          RDGeom::Point3D centerOfGroup;

          for (auto atom : (*subMolIter)->atoms)
          {
            centerOfGroup.x += atom->x2;
            centerOfGroup.y += atom->y2;
          }
          centerOfGroup.x /=  (*subMolIter)->atoms.size();
          centerOfGroup.y /=  (*subMolIter)->atoms.size();
          offset.x = (*dummyAtomIter)->x2 - centerOfGroup.x;
          offset.y = (*dummyAtomIter)->y2 - centerOfGroup.y;
        }

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

      RDGeom::Point3D centerOfGroup;

      for (auto atom : marvinSuperInfo->atoms)
      {
        // Find the atom in the main structure

        int index = this->getAtomIndex(atom);
        MarvinAtom *atomToMove = this->atoms[index];
        marvinSuperatomSgroup->atoms.push_back(atomToMove);

        this->atoms.erase(this->atoms.begin() + index);  
        if (coordsExist)  // get the center of all atoms in the group - we might use this if there are no attachement points
        {
          centerOfGroup.x += atomToMove->x2;
          centerOfGroup.y += atomToMove->y2;
        }     
      }

      // move the bonds of the group
      
      int attachmentPointsAdded = 0;
      MarvinAtom *atomPtr = NULL;
      RDGeom::Point3D centerOfAttachmentPoints;

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
            centerOfAttachmentPoints.x += atomPtr->x2;
            centerOfAttachmentPoints.y += atomPtr->y2;
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

          dummyParentAtom->x2 = centerOfAttachmentPoints.x / marvinSuperatomSgroup->attachmentPoints.size();
          dummyParentAtom->y2 = centerOfAttachmentPoints.y / marvinSuperatomSgroup->attachmentPoints.size();
        }
        else
        {
          // No attachments to the supergroup - (probably all atoms in the mol are in the supergroup) -use the center of all atoms in the group
          dummyParentAtom->x2 = centerOfGroup.x / marvinSuperatomSgroup->atoms.size();
          dummyParentAtom->y2 = centerOfGroup.y / marvinSuperatomSgroup->atoms.size();
        }
      }
    }
    
    for (auto marvinSuperInfo : this->superInfos)
      delete marvinSuperInfo;

    this->superInfos.clear();
  }

  bool MarvinMolBase::AnyOverLappingAtoms(const MarvinMolBase *otherMol) const
  {
    std::vector<std::string> firstList = getAtomList();
    std::sort(firstList.begin(), firstList.end());
    std::vector<std::string> secondList = otherMol->getAtomList();
    std::sort(secondList.begin(), secondList.end());
    std::vector<std::string> intersection; 
    std::set_intersection(firstList.begin(), firstList.end(), secondList.begin(), secondList.end(), back_inserter(intersection));

    return intersection.size() > 0;
  }


  void MarvinMol::expandMultipleSgroups()
  {
    // MulitpleSgrups are handled differently in Marvin and RDKit/mol file format
    //
    // In Marvin, the atoms of the group are indicated, and the "title" contains the number of replicates
    // In mol files, the atoms of the group are actually replicated in the atomBlock, and the MUL contructs indicates the list of
    // ALL such atoms and the list of the original parent atoms from which the others were replicated.  It also indicate the two bonds from
    // the group atoms to atoms NOT in the group (although these could be derived
    //
    //  This routine takes a MarvinMol and  prepares it for generation of an RDKit mol.  Each MultipleSgroup is expanded by:
    //    1) expanding the atom block - the new replicate sets of atoms are places just after the last atom in the original (parent) set
    //    2) records the two bonds to atoms NOT in the expanded group - one is from a "parent" atom and one from a replicated atom.
    //    3) bonds the parent and replicated sets to each other in a head to tail fashion.

      // make sure no two groups share any atoms

    for(std::vector<MarvinMultipleSgroup *>::const_iterator subMolIter = multipleSgroups.begin() ; subMolIter != multipleSgroups.end() ; ++subMolIter)
      for(std::vector<MarvinMultipleSgroup *>::const_iterator subMolIter2 = subMolIter + 1 ; subMolIter2 != multipleSgroups.end() ; ++subMolIter2)
        if ((*subMolIter)->AnyOverLappingAtoms(*subMolIter))
          throw FileParseException("Error - overlapping MultipleSgroups were detected");



    for(MarvinMultipleSgroup *subMolPtr :  multipleSgroups)
    {
      // make sure it is not already expanded - this should not happen

      if (subMolPtr->isExpanded)
         continue;

      // find the two bonds (or zero bonds) to atoms NOT in the group


      std::vector<std::string> originalConnectorAtoms;   // the bonds inthe parent that were to outside atoms BEFORE replication
      std::string connectorAtoms[2];           // working atoms to connect the replicates
      std::string outsideConnectorAtom;
      std::vector<MarvinBond *> bondsInGroup;
      subMolPtr->bondsToAtomsNotInExpandedGroup.clear();

      for (MarvinBond *bondPtr : bonds)
      {
        bool atom1InSet = boost::algorithm::contains(subMolPtr->atoms, std::vector<std::string>{bondPtr->atomRefs2[0]}, atomRefInAtoms);  
        bool atom2InSet = boost::algorithm::contains(subMolPtr->atoms, std::vector<std::string>{bondPtr->atomRefs2[1]}, atomRefInAtoms); 
        if (atom1InSet && atom2InSet)
          bondsInGroup.push_back(bondPtr);
        else if (atom1InSet)  
        {
            originalConnectorAtoms.push_back(bondPtr->atomRefs2[0]);
            subMolPtr->bondsToAtomsNotInExpandedGroup.push_back(bondPtr);
            outsideConnectorAtom = bondPtr->atomRefs2[1];  // last one set wins - when done should be the outside atom for the 2nd bond to the outside
        }
        else if (atom2InSet)  
        {
            originalConnectorAtoms.push_back(bondPtr->atomRefs2[1]);
            subMolPtr->bondsToAtomsNotInExpandedGroup.push_back(bondPtr);
            outsideConnectorAtom = bondPtr->atomRefs2[0];   // last one set wins - when done should be the outside atom for the 2nd bond to the outside
        }
      }

      //should be two bonds to the outside
      if (subMolPtr->bondsToAtomsNotInExpandedGroup.size() != 2 && subMolPtr->bondsToAtomsNotInExpandedGroup.size() != 0)
        throw FileParseException("Error - there must be zero or two bonds from the group to atoms not in the group");

      // find the hightest index of an atom in the Group - we will insert replicated atoms after that one

      std::vector<MarvinAtom *>::const_iterator lastAtomInGroupPtr = atoms.begin();
      for (MarvinAtom *atomPtr : subMolPtr->atoms)
      {
        std::vector<MarvinAtom *>::const_iterator parentAtomPtr = std::find(atoms.begin(), atoms.end(), atomPtr);
        if (parentAtomPtr > lastAtomInGroupPtr)
          lastAtomInGroupPtr = parentAtomPtr;
      }

      // ready to make copies of the group

      int copyCount = std::stoi(subMolPtr->title);

      std::vector<MarvinAtom *> atomsToAdd;
      std::vector<MarvinBond *> bondsToAdd;
      if (subMolPtr->bondsToAtomsNotInExpandedGroup.size() == 2)
        connectorAtoms[1]  = originalConnectorAtoms[1];  // this one does NOT have an _# on it

      for (int copyIndex = 2 ; copyIndex <= copyCount ; ++ copyIndex)  // start at 1, we already have the first copy
      {
        for (auto parentAtomPtr : subMolPtr->parentAtoms)
        {
          MarvinAtom *copyAtom =  new MarvinAtom( *parentAtomPtr, parentAtomPtr->id + "_" + std::to_string(copyIndex));
          atomsToAdd.push_back(copyAtom);
          subMolPtr->atoms.push_back(copyAtom);
        }

        // add the bonds for this copy

        for (auto parentBondPtr : bondsInGroup)
        {
          MarvinBond *copyBond =  new MarvinBond( *parentBondPtr, parentBondPtr->id + "_" + std::to_string(copyIndex)
              , parentBondPtr->atomRefs2[0] + "_" + std::to_string(copyIndex)
              , parentBondPtr->atomRefs2[1] + "_" + std::to_string(copyIndex));
          bondsToAdd.push_back(copyBond);
        }

        // add the bond from the last group to the new one
        
        if (subMolPtr->bondsToAtomsNotInExpandedGroup.size() == 2)
        {
          connectorAtoms[0] = originalConnectorAtoms[0] + "_" + std::to_string(copyIndex);
          MarvinBond *connectorBond =  new MarvinBond( *subMolPtr->bondsToAtomsNotInExpandedGroup[0]
              , subMolPtr->bondsToAtomsNotInExpandedGroup[0]->id + "_" + std::to_string(copyIndex) 
              , connectorAtoms[0]
              , connectorAtoms[1]);
          bondsToAdd.push_back(connectorBond);

          // update connectorAtom[1] for the next pass (or the final bond to the outside atom)

          connectorAtoms[1] = originalConnectorAtoms[1] + "_" + std::to_string(copyIndex);        
        }
      }

      // fix the bond to the second outside atom to reflect that it comes from the last replicate

      if (subMolPtr->bondsToAtomsNotInExpandedGroup.size() == 2)
      {
        subMolPtr->bondsToAtomsNotInExpandedGroup[1]->atomRefs2[0] = outsideConnectorAtom;
        subMolPtr->bondsToAtomsNotInExpandedGroup[1]->atomRefs2[1] = connectorAtoms[1];
      }

      //now add the new atoms and bonds

      atoms.insert(lastAtomInGroupPtr + 1, atomsToAdd.begin(), atomsToAdd.end());
      bonds.insert(bonds.end(), bondsToAdd.begin(), bondsToAdd.end());
   
      subMolPtr->isExpanded = true;   
    }
  }

 void MarvinMol::contractMultipleSgroups()
  {
    //this routine takes the epxpanded MultipleSgroup (which comes from the RDKit version, or is ready to product the RDKit version), and contracts it
    // to the Marbin format version.
    // the replicates are deleted (atoms and bonds), and the two orphaned bonds are  repaired 

    for(MarvinMultipleSgroup *subMolPtr :  multipleSgroups)
    {
      // make sure it is not already contracted - this should not happen

      if (!subMolPtr->isExpanded)
         continue;

    // find the atoms to be deleted
  
      std::vector<MarvinAtom *> atomsToDelete;   
      std::vector<MarvinBond *> bondsToDelete;   
      std::vector<MarvinBond *> orphanedBonds;   
      std::vector<std::string> orphanedAtomIds;   
      
      for (MarvinAtom * atomPtr : subMolPtr->atoms)
      {
        // Atoms in the group but NOT in the parent part are to be deleted
        if (!boost::algorithm::contains(subMolPtr->parentAtoms, std::vector<std::string>{atomPtr->id}, atomRefInAtoms))
          atomsToDelete.push_back(atomPtr);
      }

      for (MarvinBond *bondPtr : bonds)
      {
        bool atom1InSet = boost::algorithm::contains(atomsToDelete, std::vector<std::string>{bondPtr->atomRefs2[0]}, atomRefInAtoms);  
        bool atom2InSet = boost::algorithm::contains(atomsToDelete, std::vector<std::string>{bondPtr->atomRefs2[1]}, atomRefInAtoms); 
        if (atom1InSet && atom2InSet)
          bondsToDelete.push_back(bondPtr);
        else if (atom1InSet) 
        {  
          orphanedBonds.push_back(bondPtr);
          orphanedAtomIds.push_back(bondPtr->atomRefs2[1]);  // second atom is orphaned
        }
        else if (atom2InSet)   
        {
          orphanedBonds.push_back(bondPtr);
          orphanedAtomIds.push_back(bondPtr->atomRefs2[0]);  // first atom is orphaned

        }
      }

      // remove the atoms

      for (MarvinAtom *atomPtr : atomsToDelete)
      {
          atoms.erase(std::find(atoms.begin(), atoms.end(), atomPtr));
          subMolPtr->atoms.erase(std::find(subMolPtr->atoms.begin(), subMolPtr->atoms.end(), atomPtr));  // also remove from the submol
          delete atomPtr;
      }
      atomsToDelete.clear();

      //remove the bonds

      for (MarvinBond *bondPtr : bondsToDelete)
      {
          bonds.erase(std::find(bonds.begin(), bonds.end(), bondPtr));
          delete bondPtr;
      }
      bondsToDelete.clear();

      // now fix the orphaned bond - The first gets the atoms from the second is was NOT removed.
      // the second orphaned bond is deleted

      orphanedBonds[0]->atomRefs2[0] = orphanedAtomIds[0];
      orphanedBonds[0]->atomRefs2[1] = orphanedAtomIds[1];
      auto bondToDelete = std::find(bonds.begin(), bonds.end(), orphanedBonds[1]);
      delete *bondToDelete;
      bonds.erase(bondToDelete);  // delete the second orphaned atom
   
      subMolPtr->isExpanded = false;   
    }
  }


  void MarvinMol::processMulticenterSgroups()
  {
    //This routine removes the MultiCenter groups - RDKit does not have a place to put them
    // also removes the X atom (the center atom) and its bonds

    for(std::vector<MarvinMulticenterSgroup *>::iterator subMolIter =  multicenterSgroups.begin() ; subMolIter != multicenterSgroups.end() ; ++subMolIter)
    {
  
      //delete the bonds to the dummy atom
      std::vector<MarvinBond *> orphanedBonds;  // list of bonds to delete
      for (auto orphanedBond = bonds.begin(); orphanedBond != bonds.end() ; ++orphanedBond)
      {
        std::string attachedAtomId = "";
        if ((*orphanedBond)->atomRefs2[0] == (*subMolIter)->center->id || (*orphanedBond)->atomRefs2[1] == (*subMolIter)->center->id)
        {
         orphanedBonds.push_back(*orphanedBond);
        }
      }
      for (auto orphanedBond = orphanedBonds.begin(); orphanedBond != orphanedBonds.end() ; ++orphanedBond)
      {
        delete *orphanedBond;
        auto orphanedBondIter = std::find(bonds.begin(), bonds.end(), *orphanedBond);
        bonds.erase(orphanedBondIter);
      }
      
      auto centerIter = std::find(atoms.begin(), atoms.end(), (*subMolIter)->center);
      delete (*subMolIter)->center;  // get rid of the MolAtom
      atoms.erase(centerIter);   // get rid of the atoms pointer to the old dummy atom

      delete *subMolIter;
    }
    this->multicenterSgroups.clear();
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
    for (std::vector<MarvinMulticenterSgroup *>::const_iterator it = multicenterSgroups.begin();  it != multicenterSgroups.end(); ++it)
      out << (*it)->toString();    
    for (std::vector<MarvinGenericSgroup *>::const_iterator it = genericSgroups.begin();  it != genericSgroups.end(); ++it)
      out << (*it)->toString();   
    for (std::vector<MarvinMonomerSgroup *>::const_iterator it = monomerSgroups.begin();  it != monomerSgroups.end(); ++it)
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

  void MarvinReaction::expandMultipleSgroups()
  {
    //This routine epxans MultipleSgrpoup to have full replicates of the repeating group in preparation for creating an RDKut version

    for (std::vector<MarvinMol *>::iterator molIter = reactants.begin() ; molIter != reactants.end() ; ++molIter)
    {
      (*molIter)->expandMultipleSgroups();
    }
    for (std::vector<MarvinMol *>::iterator molIter = agents.begin() ; molIter != agents.end() ; ++molIter)
    {
      (*molIter)->expandMultipleSgroups();
    }       
    for (std::vector<MarvinMol *>::iterator molIter = products.begin() ; molIter != products.end() ; ++molIter)
    {
      (*molIter)->expandMultipleSgroups();
    }       
  }

  
  void MarvinReaction::processMulticenterSgroups()
  {
    //This routine removes all multiCenter Sgroups

    for (std::vector<MarvinMol *>::iterator molIter = reactants.begin() ; molIter != reactants.end() ; ++molIter)
      (*molIter)->processMulticenterSgroups();
    for (std::vector<MarvinMol *>::iterator molIter = agents.begin() ; molIter != agents.end() ; ++molIter)
      (*molIter)->processMulticenterSgroups();
    for (std::vector<MarvinMol *>::iterator molIter = products.begin() ; molIter != products.end() ; ++molIter)
      (*molIter)->processMulticenterSgroups();
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

  MarvinRectangle::MarvinRectangle()
  {
    upperLeft.x = 0.0;
    upperLeft.y = 0.0;
    lowerRight.x = 0.0;
    lowerRight.y = 0.0;
    centerIsStale = true;
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

  MarvinRectangle::MarvinRectangle(const std::vector<MarvinRectangle> rects)
  {
    centerIsStale = true;

    if (rects.size() == 0)
      return;
    upperLeft.x = DBL_MAX;
    upperLeft.y = -DBL_MAX;
    lowerRight.x = -DBL_MAX;
    lowerRight.y = DBL_MAX;


    for (auto rect : rects)
    {
      this->extend(rect);
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

  bool MarvinRectangle::compareRectanglesByYReverse(MarvinRectangle &r1, MarvinRectangle &r2)
  {
    return (r1.getCenter().y > r2.getCenter().y);
  }

  MarvinStereoGroup::MarvinStereoGroup(StereoGroupType grouptypeInit, int groupNumberInit )
  {
    groupType = grouptypeInit;
    groupNumber = groupNumberInit;
  }
}
