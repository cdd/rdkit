//
//  Copyright (C) 2002-2022 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RD_MARVINDEFS_H
#define RD_MARVINDEFS_H

#include <GraphMol/RDKitBase.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/LocaleSwitcher.h>

#include <string>
#include <iostream>

namespace RDKit 
{

    class MarvinWriterException : public std::runtime_error 
    {
        public:
        explicit MarvinWriterException(std::string message)
            : std::runtime_error(message){};
    };

    class MarvinArrow
    {
        public:
        std::string type;
        double x1;
        double y1;
        double x2;
        double y2;

        std::string toString() const
        {

            std::ostringstream out;
            out << "<arrow type=\"" << type << "\" x1=\"" << x1<< "\" y1=\"" << y1 << "\" x2=\"" << x2<< "\" y2=\"" << y2 << "\"/>";
        
            return out.str();
        }
    };

    class MarvinPlus
    {
        public:
        std::string id;
        double x1;
        double y1;
        double x2;
        double y2;


        std::string toString() const
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
    };

    class MarvinCondition
    {
        public:
        std::string id;
        std::string text;
        double x;
        double y;
        double fontScale;

        std::string halign;
        std::string valign;

        std::string toString() const
        {
            std::ostringstream out;

            out << "<MTextBox id=\"" << id << "\" toption=\"NOROT\" fontScale=\"" << fontScale << "\" halign=\"" << halign << "\" valign=\"" << valign << "\" autoSize=\"true\">"
                "<Field name=\"text\">" << text << "</Field>"
                "<MPoint x=\"" << x << "\" y=\"" << y << "\"/>"
                "<MPoint x=\"" << x << "\" y=\"" << y << "\"/>"
                "<MPoint x=\"" << x << "\" y=\"" << y << "\"/>"
                "<MPoint x=\"" << x << "\" y=\"" << y << "\"/>"
                "</MTextBox>";

            return out.str();
        }
    };

    
    class MarvinAttachmentPoint
    {
        public:

        // <AttachmentPoint atom="a7" order="1" bond="b6"/>
        std::string atom;
        std::string bond;
        std::string order;

        std::string toString() const
        {
            std::ostringstream out;

            out << "<AttachmentPoint atom=\"" << atom << "\" order=\"" << order << "\" bond=\"" << bond << "\"/>";

            return out.str();
        }
    };

    class MarvinAtom
    {
        public:
        std::string id;
        std::string elementType;
        double x2;
        double y2;
        int formalCharge;
        std::string radical;
        int isotope;
        std::string mrvAlias;
        std::string mrvStereoGroup;
        int mrvMap;
        std::string sgroupRef;
        std::string sgroupAttachmentPoint;

        MarvinAtom()
            : x2(DBL_MIN)
            , y2(DBL_MIN)
            , formalCharge(0)
            , mrvMap(0)
        {
        }

    
        bool operator==(const MarvinAtom& rhs) const
        {
            return this->id == rhs.id;
        }
        
        bool operator==(const MarvinAtom *rhs) const
        {
            return this->id == rhs->id;
        }

        std::string toString() const
        {
            // <atom id="a7" elementType="C" x2="15.225" y2="-8.3972" sgroupAttachmentPoint="1"/>

            std::ostringstream out;
            out << "<atom id=\"" << id << "\" elementType=\"" << elementType << "\" x2=\"" << x2 << "\" y2=\"" << y2;

            if (formalCharge !=0)
            out << " formalCharge=\"" << formalCharge << "\"";

            if (radical != "")
            out << " radical=\"" << radical << "\"";

            if (isotope != 0)
            out << " isotope=\"" << isotope << "\"";

            if (mrvAlias != "")
            out << " mrvAlias=\"" << mrvAlias << "\"";

            if (mrvStereoGroup != "")
            out << " mrvStereoGroup=\"" << mrvStereoGroup << "\"";

            if (mrvMap != 0)
            out << " mrvMap=\"" << mrvMap << "\"";

            if (sgroupRef != "")
            out << " sgroupRef=\"" << sgroupRef << "\"";

            if (sgroupAttachmentPoint != "")
            out << " sgroupAttachmentPoint=\"" << sgroupAttachmentPoint << "\"";

            out << "/>";

            return out.str();
        }     
    };

    class MarvinBond
    {
        public:
        std::string id;
        std::string atomRefs2[2];
        std::string order;
        std::string bondStereo;
        std::string queryType;
    
        bool isEqual(const MarvinAtom& other) const
        {
            return this->id == other.id;
        }

        bool operator==(const MarvinAtom& rhs) const
        {
            return this->isEqual(rhs);
        }

        std::string toString() const
        {
        // <bond id="b8" atomRefs2="a1 a7" order="1">

            std::ostringstream out;
            
            out << "<bond id=\"" << id << "\" atomRefs2=\"" << atomRefs2[0] << " " << atomRefs2[1] << "\" order=\"" << order << "\"";
        

            if (queryType != "")
            out << " queryType=\"" << queryType << "\"";
            

            if (bondStereo != "")
            out << "><bondStereo>" << bondStereo << "</bondStereo></bond>";
            else
            out << "/>";          
            return out.str();
    }
    };

    class MarvinMolBase
    {
        public:
        std::string molID;
        std::vector<MarvinAtom *> atoms;
        std::vector<MarvinBond *> bonds;      

        
        virtual std::string role() = 0;

        MarvinMolBase()
        {

        }

        virtual ~MarvinMolBase()
        {
            for ( std::vector<MarvinAtom *>::iterator it = atoms.begin(); it != atoms.end(); ++it)
            delete(*it);
            for ( std::vector<MarvinBond *>::iterator it = bonds.begin(); it != bonds.end(); ++it)
            delete(*it);
        }
        
        int getAtomIndex(std::string id)
        {
            auto atomIter = find_if(atoms.begin(), atoms.end(), [id](const MarvinAtom *arg) { return arg->id == id; });
            if (atomIter != atoms.end())
                return atomIter - atoms.begin();
            else 
                return -1;
        }

        int getBondIndex(std::string id)
        {
            auto bondIter = find_if(bonds.begin(), bonds.end(), [id](const MarvinBond *arg) { return arg->id == id; });
            if (bondIter != bonds.end())
                return bondIter - bonds.begin();
            else 
                return -1;
        }


    };

    class MarvinSruSgroup : public MarvinMolBase
    {
        public:
        std::string id;
        std::vector<std::string> atomRefs;
        std::string title;
        std::string connect;
        std::string correspondence;
        std::vector<std::string> bondList;

        MarvinSruSgroup()
        {
            
        }
        ~MarvinSruSgroup()
        {
            
        }

        std::string toString() const
        {
            std::ostringstream out;

            out << "<molecule molID=\"" << molID << "\" id=\"" << id << "\" role=\"SruSgroup\" atomRefs=\"" << boost::algorithm::join(atomRefs,",") << "\" title=\"" << title 
            <<"\" connect=\"" << connect << "\" correspondence=\"" << correspondence << "\" bondList=\"" << boost::algorithm::join(bondList,",") << "\"/>";

            return out.str();
        }

        std::string role()
        {
            return std::string("SruSgroup");
        } 

    };

    
    class MarvinSuperatomSgroup : public MarvinMolBase
    {
        public:
        std::string id;
        std::string title;
        std::vector<MarvinAttachmentPoint *> attachmentPoints;

        MarvinSuperatomSgroup()
        {

        }

      ~MarvinSuperatomSgroup()
      {
        for ( std::vector<MarvinAttachmentPoint *>::iterator it = attachmentPoints.begin(); it != attachmentPoints.end(); ++it)
            delete(*it);
      }

      std::string role()
      {
        return std::string("SuperatomSgroup");
      } 

      std::string toString() const
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
    };

    class MarvinSuperInfo   // used in converting superatoms group to mol stule groups
    {
      public:
      std::string title;
      std::vector<std::string> atoms;
    };
    
    class MarvinMol : public MarvinMolBase
    {
      public:

      std::vector<MarvinSuperatomSgroup *> superatomSgroups;
      std::vector<MarvinSruSgroup *>  sruSgroups;
      std::vector<MarvinSuperInfo *> superInfos;  // used in convertng superatomSgroups to mol-type CTs

      MarvinMol()
      {

      }

      ~MarvinMol()
      {
        for (std::vector<MarvinSuperatomSgroup *>::iterator it = superatomSgroups.begin(); it != superatomSgroups.end(); ++it)
            delete(*it);
        for (std::vector<MarvinSruSgroup *>::iterator it = sruSgroups.begin(); it != sruSgroups.end(); ++it)
            delete(*it);
        for (std::vector<MarvinSuperInfo *>::iterator it = superInfos.begin(); it != superInfos.end(); ++it)
            delete(*it);
      }

      std::string role()
      {
            return std::string("Parent");
      } 

      static bool atomRefInAtoms(MarvinAtom *a, std::string b )
      { 
            return a->id == b; 
      }

      static bool bondRefInBonds(MarvinBond *a, std::string b )
      { 
            return a->id == b; 
      }

      void cleanUpNumbering(int &molCount    // this is the starting mol count, and receives the ending mol count - THis is used when MarvinMol->convertToSuperAtaoms is called multiple times from a RXN
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

        for (MarvinSruSgroup *marvinmarvinSruSGgroup : this->sruSgroups)         
        {
            std::string newId = "sg" + std::to_string(++sgCount);
            sgMap[marvinmarvinSruSGgroup->id] = newId;
            marvinmarvinSruSGgroup->id = newId;
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
            atomMap[bondPtr->id] =newId;
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
              atomMap[bondPtr->id] =newId;
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

          for (MarvinSruSgroup *marvinmarvinSruSGgroup : this->sruSgroups)
          {
              for (auto atomRefPtr : marvinmarvinSruSGgroup->atomRefs)
                  atomRefPtr = atomMap[atomRefPtr];

              for (auto bondListPtr : marvinmarvinSruSGgroup->bondList)
                  bondListPtr = bondMap[bondListPtr];
          }
        }
      }

      void convertFromSuperAtoms()
      {
        // the mol-style super atoms are significatnly different than the Marvin super atoms.  The mol style has all atoms and bonds in the main mol, and parameter lines that
        // indicate the name of the super atom and the atom indices affected.
        //
        // The Marvin as a dummy atom with the super atom name in the main molecule, and a separate sub-mol with the atoms and bonds of the superatom.  It also has a separate record 
        // called attachmentPoint that atom in the super atom sub=mol that is replaces the dummy atom, and also a bond pointer and bond order in case the bond order changes.
        //
        //This routine copies the aboms and bonds from the sub=mol to the parent mol, and deleted the dummy atom form the parent mol.  It also saves the infor needed to make a mol-file type
        // superatom in the MarvinSuperInfo array.

        for(std::vector<MarvinSuperatomSgroup *>::iterator subMolIter =  superatomSgroups.begin() ; subMolIter != superatomSgroups.end() ; ++subMolIter)
        {
            //save the name of the superatom
            
            auto marvinSuperInfo = new MarvinSuperInfo();
            this->superInfos.push_back(marvinSuperInfo);

            marvinSuperInfo->title = (*subMolIter)->title;

            //  remove and delete the dummy atom from the parent.

            auto dummyAtomIter = find_if(atoms.begin(), atoms.end(), [subMolIter](const MarvinAtom *arg) { 
                                                return arg->sgroupRef == (*subMolIter)->id; });
            if (dummyAtomIter != atoms.end())
            {
                delete *dummyAtomIter;  // get rid of the MolAtom
                atoms.erase(dummyAtomIter);   // get rid of the atoms pointer to the old dummy atom
            }

            // add the atoms and bonds from the super group to the parent

            for (std::vector<MarvinAtom *>::iterator subAtomIter = (*subMolIter)->atoms.begin() ; subAtomIter != (*subMolIter)->atoms.end() ; ++subAtomIter)
            {
                atoms.push_back( *subAtomIter);
                marvinSuperInfo->atoms.push_back((*subAtomIter)->id);

                // remove the sgroupRef from the atom (only one will have it)
                (*subAtomIter)->sgroupAttachmentPoint = "";

            }
            for (std::vector<MarvinBond *>::iterator subBondIter = (*subMolIter)->bonds.begin()  ; subBondIter != (*subMolIter)->bonds.end() ; ++subBondIter)
            {
                bonds.push_back( *subBondIter);
            }

            // process the attachment points - fix the bond that was made wrong by deleting the dummy atom

            for (std::vector<MarvinAttachmentPoint *>::iterator attachIter = (*subMolIter) -> attachmentPoints.begin() ;  attachIter != (*subMolIter) -> attachmentPoints.end(); ++attachIter)
            {
                // find the bond in the parent
                
                auto bondIter = find_if(bonds.begin(), bonds.end(), [attachIter](const MarvinBond *arg) { 
                                                return arg->id == (*attachIter)->bond; });
                if (bondIter == bonds.end())
                throw FileParseException("Bond specification for an AttachmentPoint definition was not found in the bond array in MRV file");

                // one of the two atoms in the bond is NOT in the mol - we deleted the dummy atom.

                int atomIndex;
                for (atomIndex = 0 ; atomIndex < 2 ; ++atomIndex)
                {
                    if (!boost::algorithm::contains(atoms, std::vector<std::string>{(*bondIter)->atomRefs2[atomIndex]}, atomRefInAtoms ))
                    {
                        (*bondIter)->atomRefs2[atomIndex] = (*attachIter)->atom;   // the attach atom
                        (*bondIter)->order = (*attachIter)->order;  // fix the bond order
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
            
            // clean up - delete the super atom mols and the attach points

            (*subMolIter)->attachmentPoints.clear();
            (*subMolIter)->atoms.clear();
            (*subMolIter)->bonds.clear();
            delete *subMolIter;
        }
        this->superatomSgroups.clear();
      }

      void convertToSuperAtoms()
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

          marvinSuperatomSgroup->molID = "AS";
          marvinSuperatomSgroup->title = marvinSuperInfo->title;
          marvinSuperatomSgroup->id = "nsg" + std::to_string(superGroupsAdded);  // n in nsg ensures no colllision with other sgs already in the mol 

          for (auto atom : marvinSuperInfo->atoms)
          {
              // Find the atom in the main structure

              int index = this->getAtomIndex(atom);
              MarvinAtom *atomToMove = this->atoms[index];
              marvinSuperatomSgroup->atoms.push_back(atomToMove);
              this->atoms.erase(this->atoms.begin() + index);
          }

          // move the bonds of the group
          
          MarvinBond *attachmentBondInParent = NULL;
          std::string atomInGroup = "";
          
          for (auto bond : this->bonds)
          {
            bool atom1IsInGroup = boost::algorithm::contains(marvinSuperatomSgroup->atoms, std::vector<std::string>{(bond)->atomRefs2[0]}, atomRefInAtoms );
            bool atom2IsInGroup = boost::algorithm::contains(marvinSuperatomSgroup->atoms, std::vector<std::string>{(bond)->atomRefs2[1]}, atomRefInAtoms );


            if (atom1IsInGroup && atom2IsInGroup )  // both are in, so move the bond
            {
                int index = this->getBondIndex(bond->id);
                MarvinBond *bondToMove = this->bonds[index];
                marvinSuperatomSgroup->bonds.push_back(bondToMove);
                this->bonds.erase(this->bonds.begin() + index);
            }

            // see if one atom of the bond is in the group to be created.  
            else if (atom1IsInGroup)
            {
                if (atomInGroup != "")  // already found - we aonly alow one attachment point
                    throw MarvinWriterException("Multiple attachment points for SuperGroup is not supported");
                
                
                atomInGroup = bond->atomRefs2[0];
                bond->atomRefs2[0] = newAtomName;  // fix the bond that was attached to the moved attachment atom
                attachmentBondInParent = bond;

            }
            else if  (atom2IsInGroup)
            {
                if (atomInGroup != "")  // already found - we aonly alow one attachment point
                    throw  MarvinWriterException("Multiple attachment points for SuperGroup is not supported");
                
                atomInGroup = bond->atomRefs2[1];
                bond->atomRefs2[1] = newAtomName;  // fix the bond that was attached to the moved attachment atom
                attachmentBondInParent = bond;
            }

            // now if we found one atom of a bond that was in the group, we have an attachment point.

            if (atomInGroup != "")
            {
                // add an attachment  atom to the parent 

                MarvinAtom *attachAtom = new MarvinAtom();
                this->atoms.push_back(attachAtom);
                attachAtom->elementType = "R";
                attachAtom->id = newAtomName;
                attachAtom->sgroupAttachmentPoint = marvinSuperInfo->title;
                MarvinAtom *atomPtr = marvinSuperatomSgroup->atoms[marvinSuperatomSgroup->getAtomIndex(atomInGroup)];
                attachAtom->x2 = atomPtr->x2;
                attachAtom->y2 = atomPtr->y2;
                atomPtr->sgroupAttachmentPoint = "1";

                // add an attachentPoint structure

                auto marvinAttachmentPoint = new MarvinAttachmentPoint();
                marvinSuperatomSgroup->attachmentPoints.push_back(marvinAttachmentPoint);
                marvinAttachmentPoint->atom = atomInGroup;
                marvinAttachmentPoint->bond = attachmentBondInParent->id;
                marvinAttachmentPoint->order = attachmentBondInParent->order;
            }
          } 
        }
          
        this->superatomSgroups.clear();
      }

     
  
      std::string toString() const
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
        for (std::vector<MarvinSruSgroup *>::const_iterator it = sruSgroups.begin();  it != sruSgroups.end(); ++it)
          out << (*it)->toString();
        
        out <<"</molecule>";
        
        return out.str();


      }   

      std::string generateMolString()
      {
        std::ostringstream out;

        out << "<cml xmlns=\"http://www.chemaxon.com\" version=\"ChemAxon file format v20.20.0, generated by RDKit\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
            "xsi:schemaLocation=\"http://www.chemaxon.com http://www.chemaxon.com/marvin/schema/mrvSchema_20_20_0.xsd\">"
            "<MDocument><MChemicalStruct>";

        out << toString();
      
        out << "</MChemicalStruct></MDocument></cml>";
        return out.str();
      }

    };

  

    class MarvinReaction
    {
      public:

      std::vector<MarvinMol *> reactants;
      std::vector<MarvinMol *> agents;
      std::vector<MarvinMol *> products;

      MarvinArrow arrow;
      std::vector<MarvinPlus *> pluses;
      std::vector<MarvinCondition *> conditions;

      ~MarvinReaction()
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

      void convertFromSuperAtoms()
      {
        //This routine converts all the mols in the rxn to be ready for conversion to RDKIT mols

        for (std::vector<MarvinMol *>::iterator molIter = reactants.begin() ; molIter != reactants.end() ; ++molIter)
          (*molIter)->convertFromSuperAtoms();
        for (std::vector<MarvinMol *>::iterator molIter = agents.begin() ; molIter != agents.end() ; ++molIter)
          (*molIter)->convertFromSuperAtoms();
        for (std::vector<MarvinMol *>::iterator molIter = products.begin() ; molIter != products.end() ; ++molIter)
          (*molIter)->convertFromSuperAtoms();
        
      }

      void convertToSuperAtoms()
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

      std::string toString()
      {
        std::ostringstream out;

        out << "<cml xmlns=\"http://www.chemaxon.com\" version=\"ChemAxon file format v20.20.0, generated by RDKit\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\""
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
    };

    class MarvinRectangle
    {
      private:
      RDGeom::Point3D *center;
      
      public:
      RDGeom::Point3D upperLeft;
      RDGeom::Point3D lowerRight;

      MarvinRectangle(double left, double right, double top, double bottom)
      {
        upperLeft.x = left;
        upperLeft.y = top;
        lowerRight.x = right;
        lowerRight.y = bottom;
        center = NULL;
      }

      MarvinRectangle(const RDGeom::Point3D &upperLeftInit, const RDGeom::Point3D &lowerRightInit)
      {
        upperLeft = upperLeftInit;
        lowerRight = lowerRightInit;
        center = NULL;
      }

      MarvinRectangle(const std::vector<MarvinAtom *> atoms)
      {
        center = NULL;

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

      void extend(const MarvinRectangle &otherRectangle)
      {  
        if (otherRectangle.upperLeft.x < upperLeft.x)
          upperLeft.x = otherRectangle.upperLeft.x;
        if (otherRectangle.lowerRight.x > lowerRight.x)
          lowerRight.x = otherRectangle.lowerRight.x;

        
        if (otherRectangle.upperLeft.y > upperLeft.y)
          upperLeft.y = otherRectangle.upperLeft.y;
        if (otherRectangle.lowerRight.y < lowerRight.y)
          lowerRight.y = otherRectangle.lowerRight.y;
        
        center = NULL;
      }

      RDGeom::Point3D getCenter()
      {
        if (center == NULL)
        {
          center = new RDGeom::Point3D();
          center->x = (lowerRight.x + upperLeft.x)/2.0;
          center->y = (lowerRight.y + upperLeft.y)/2.0;
        }
        return *center;
      }

      bool overlapsVertically(const MarvinRectangle &otherRectangle) const
      {
        if (otherRectangle.upperLeft.y < lowerRight.y || otherRectangle.lowerRight.y > upperLeft.y)
          return false;
        return true;
      }

      
      bool overlapsVHorizontally(const MarvinRectangle &otherRectangle) const
      {
        if (otherRectangle.upperLeft.x > lowerRight.x || otherRectangle.lowerRight.x < upperLeft.x)
          return false;
        return true;
      }

      static bool compareRectanglesByX(MarvinRectangle &r1, MarvinRectangle &r2)
      {
          return (r1.getCenter().x < r2.getCenter().x);
      }

      static bool compareRectanglesByY(MarvinRectangle &r1, MarvinRectangle &r2)
      {
          return (r1.getCenter().y < r2.getCenter().y);
      }
    };

    class MarvinStereoGroup
    {
      public:

      StereoGroupType groupType;  // one of ABS AND OR
      int groupNumber;
      std::vector<unsigned int> atoms;

      MarvinStereoGroup(StereoGroupType grouptypeInit, int groupNumberInit )
      {
        groupType = grouptypeInit;
        groupNumber = groupNumberInit;
      }
    };


}

#endif //RD_MARVINDEFS_H

