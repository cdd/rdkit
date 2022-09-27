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
  class MarvinCML
  {
    
    public:
      RWMol *mol;
      ChemicalReaction *rxn;

    const std::vector<std::string> sruSgroupConnectChoices  { "hh","ht","eu"};
    const std::vector<std::string> marvinBondOrders  { "1","2","3","A"};
    const std::vector<std::string> marvinQueryBondsTypes  {"SD","SA","DA"};
    const std::vector<std::string> marvinRadicalVals  {"monovalent",   "divalent1", "trivalent4", "4"};
    

    public:
    MarvinCML()
        : mol(NULL)
        , rxn(NULL)
    {
    };

    ~MarvinCML()
    { 
    };


    private:

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

      void convertSuperAtoms()
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

      void convertSuperAtoms()
      {
        //This routine converts all the mols in the rxn to be ready for conversion to RDKIT mols

        for (std::vector<MarvinMol *>::iterator molIter = reactants.begin() ; molIter != reactants.end() ; ++molIter)
          (*molIter)->convertSuperAtoms();
        for (std::vector<MarvinMol *>::iterator molIter = agents.begin() ; molIter != agents.end() ; ++molIter)
          (*molIter)->convertSuperAtoms();
        for (std::vector<MarvinMol *>::iterator molIter = products.begin() ; molIter != products.end() ; ++molIter)
          (*molIter)->convertSuperAtoms();
        
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


    public:
 
    //this routine does the work of parsing.  It returns either an RWMol * or a ChemicalStructure *
    // either way it is cast to a void *

    void *parse(std::istream &is, bool &isReaction)
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
        // do nothing try a reaction
      }
      if (molFlag)
      {
        mol = (RWMol *) parseMolecule(molOrRxn, true);
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

      rxn = parseReaction(molOrRxn, tree.get_child("cml.MDocument.MChemicalStruct.reaction")); 
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
        
        if (symb == "R") 
        {
            auto *query = new QueryAtom(0);
            res = query;
            query->setQuery(makeAtomNullQuery());
            // queries have no implicit Hs:
            res->setNoImplicit(true);
        } 
        else if (symb[0] == 'R' && symb >= "R0" && symb <= "R99") 
        {
          auto *query = new QueryAtom(0);
          res = query;
          query->setQuery(makeAtomNullQuery());

          if ( marvinAtom->isotope == 0)
          {
            std::string rlabel = "";
            rlabel = symb.substr(1, symb.length() - 1);
            int rnumber;
            try 
            {
              rnumber = boost::lexical_cast<int>(rlabel);
            } catch (boost::bad_lexical_cast &) {
              rnumber = -1;
            }
            if (rnumber >= 0) 
            {
              res->setIsotope(rnumber);
            }
          }
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

        //  we no not get parity (chirality) frommarvin file
        //res->setProp(common_properties::molParity, parity);
        //res->setProp(common_properties::molStereoCare, stereoCare);
        
        // total valence is not parsed from marvin file
        //res->setProp(common_properties::molTotValence, totValence);
      
        // rxnRole is not parsed from marvin file
        //res->setProp(common_properties::molRxnRole, rxnRole);
        //res->setProp(common_properties::molRxnComponent, rxnComponent);
        
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
              bType = 6;
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
        } else if (temp == "h")  
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
                ->getPropIfPresent(common_properties::molStereoCare, care2)) {
          if (care1 && care2) {
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
 
    RWMol *parseMolecule(MarvinMol *marvinMol, bool sanitize=true)
    {
      
      std::vector<MarvinStereoGroup *> stereoGroups;
      Conformer *conf = NULL;
      RWMol *mol = NULL;

      try
      {

        mol = new RWMol();

        mol->setProp(common_properties::_Name, marvinMol->molID);
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

          std::vector<std::string>::const_iterator atomIter;
          for (atomIter = (*sruIter)->atomRefs.begin(); atomIter != (*sruIter)->atomRefs.end(); ++atomIter)
          {
            int atomIndex = marvinMol->getAtomIndex((*atomIter));
            (sgroup.*sGroupAddIndexedElement)(atomIndex);
          }
          std::vector<std::string>::const_iterator bondIter;
          for (bondIter = (*sruIter)->bondList.begin(); bondIter != (*sruIter)->bondList.end(); ++bondIter)
          {
            //int bondIndex = marvinMol->getBondIndex((*bondIter)) + 1;
            int bondIndex = marvinMol->getBondIndex((*bondIter));
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
          (*atomIt)->calcExplicitValence(false);


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
          
          MolOps::removeHs(*mol, false, false);
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
        printf("Caught error in parseMolecule");
        delete mol;

        delete conf;

        for (std::vector<MarvinStereoGroup *>::iterator it = stereoGroups.begin() ; it != stereoGroups.end(); ++it)
          delete *it;
        
        throw;
      }
    }


    RWMol *parseMolecule(boost::property_tree::ptree molTree, bool sanitize=true)
    {
      MarvinMol *marvinMol = NULL;

      try
      {


        marvinMol = (MarvinMol *)parseMarvinMolecule(molTree);
        //printf("%s\n", marvinMol->generateMolString().c_str());
        marvinMol->convertSuperAtoms();
        //printf("%s\n", marvinMol->generateMolString().c_str());



        RWMol* mol = parseMolecule(marvinMol , sanitize);

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
            boost::algorithm::split(marvinSruSgroup->atomRefs, atomRefsStr, boost::algorithm::is_space());
            for (std::vector<std::string>::const_iterator it = marvinSruSgroup->atomRefs.begin() ;  it != marvinSruSgroup->atomRefs.end(); ++it)
              if (!boost::algorithm::contains(parentMol->atoms, std::vector<std::string>{*it}, MarvinMol::atomRefInAtoms ))
                throw FileParseException("All of the AtomRefs in the sruSgroup must in the parent molecule");

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
              boost::algorithm::split(marvinSruSgroup->bondList, bondListStr, boost::algorithm::is_space());
              for (std::vector<std::string>::const_iterator it = marvinSruSgroup->bondList.begin() ;  it != marvinSruSgroup->bondList.end(); ++it)
                if (!boost::algorithm::contains(parentMol->bonds, std::vector<std::string>{*it}, MarvinMol::bondRefInBonds ))
                  throw FileParseException("All of the BondList entries in the sruSgroup must in the parent molecule");
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

            BOOST_FOREACH(
                boost::property_tree::ptree::value_type &v, molTree.get_child("atomArray"))
            {
              MarvinAtom *mrvAtom = new MarvinAtom();   
              res->atoms.push_back(mrvAtom);

              mrvAtom->id = v.second.get<std::string>("<xmlattr>.id", "");
              mrvAtom->elementType = v.second.get<std::string>("<xmlattr>.elementType", "");
              std::string x2 = v.second.get<std::string>("<xmlattr>.x2", "");
              std::string y2 = v.second.get<std::string>("<xmlattr>.y2", "");
              if (mrvAtom->id == "" || mrvAtom->elementType == "" || x2 == "" || y2 == "" )
                  throw FileParseException("Expected id, elementType. x2 and y2 for an atom definition in MRV file");

              // x2 and y2 are doubles

              if (!getCleanDouble(x2, mrvAtom->x2) || !getCleanDouble(y2, mrvAtom->y2))
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

              mrvAtom->mrvStereoGroup = v.second.get<std::string>("<xmlattr>.mrvStereoGroup", "");

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

            if (atomID == "" ||elementType == "" || x2 == "" || y2 == ""
                    ||  elementTypes.size() != atomCount || x2s.size() != atomCount || y2s.size() != atomCount )
              throw FileParseException("Expected id, elementType. x2 and y2 arrays for an atomArray definition in MRV file, and the counts of each must the same");

            if (formalCharge != "" && formalCharges.size() != atomCount)
              throw FileParseException("formalCharges, if specified, must have the same count as the atomIDs for an atomArray definition in MRV file");

            if (radical != "" && radicals.size() != atomCount)
              throw FileParseException("radicals, if specified, must have the same count as the atomIDs for an atomArray definition in MRV file");
  
            if (isotope != "" && isotopes.size() != atomCount)
              throw FileParseException("isotopes, if specified, must have the same count as the atomIDs for an atomArray definition in MRV file");

            if (mrvAlias != "" && mrvAliases.size() != atomCount)
              throw FileParseException("mrvAliases, if specified, must have the same count as the atomIDs for an atomArray definition in MRV file");

            if (mrvStereoGroup != ""  && mrvStereoGroups.size() != atomCount)
              throw FileParseException("mrvStereoGroups, if specified, must have the same count as the atomIDs for an atomArray definition in MRV file");

            if (mrvMap != "" && mrvMaps.size() != atomCount)
              throw FileParseException("mrvMaps, if specified, must have the same count as the atomIDs for an atomArray definition in MRV file");

            if (sgroupRef != ""  && sgroupRefs.size() != atomCount)
              throw FileParseException("sgroupRefs, if specified, must have the same count as the atomIDs for an atomArray definition in MRV file");

            if (sgroupAttachmentPoint != ""  && sgroupAttachmentPoints.size() != atomCount)
              throw FileParseException("sgroupAttachmentPoint, if specified, must have the same count as the atomIDs for an atomArray definition in MRV file");

            for (size_t i = 0 ; i < atomCount ; ++i)
            {
              MarvinAtom *mrvAtom = new MarvinAtom();   
              res->atoms.push_back(mrvAtom);

              mrvAtom->id = atomIds[i];

              mrvAtom->elementType = elementTypes[i];
              
              if (!getCleanDouble(x2s[i], mrvAtom->x2) || !getCleanDouble(y2s[i], mrvAtom->y2))
                throw FileParseException("The values x2 and y2 must be large floats in MRV file");        

              if (formalCharge != "")
              {
                if (!getCleanInt(formalCharges[i], mrvAtom->formalCharge) )
                    throw FileParseException("The value for formalCharge must be an integer in MRV file"); 
              }
              else
                mrvAtom->formalCharge = 0;

              if (isotope != "")
              {
                if (!getCleanInt(isotopes[i], mrvAtom->isotope) )
                    throw FileParseException("The value for formalCharge must be an integer in MRV file"); 
              }
              else
                mrvAtom->isotope = 0;
            
            
              if (radical != "")
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

              if (mrvAlias != "")
                mrvAtom->mrvAlias = mrvAliases[i];
              else
                mrvAtom->mrvAlias = "";
              
              if (mrvStereoGroup != "")
                mrvAtom->mrvStereoGroup = mrvStereoGroups[i];
              else
                mrvAtom->mrvStereoGroup = "";       
              
              if (mrvMap != "")
              {
                if (!getCleanInt(mrvMap, mrvAtom->mrvMap) || mrvAtom->mrvMap <=0)
                    throw FileParseException("The value for mrvMap must be an non-=negative integer in MRV file"); 
              }
              else
                mrvAtom->mrvMap = 0;       

              if (sgroupRef != "")
                mrvAtom->sgroupRef = sgroupRefs[i];
              else
                mrvAtom->sgroupRef = "";
              
              if (role== "SuperatomSgroup" && sgroupAttachmentPoint != "")
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
            if (mrvBond->bondStereo != "" &&  boost::algorithm::to_lower_copy(mrvBond->bondStereo) != "w" && boost::algorithm::to_lower_copy(mrvBond->bondStereo) != "h")
                throw FileParseException("The value for bondStereo must be \"h\" or \"w\" in MRV file");

          
          }

        }

        if (role== "SuperatomSgroup")
        {

          // see if there is an AttachmentPointArray
          bool found;
          
          try 
          {
            BOOST_FOREACH(
                boost::property_tree::ptree::value_type &v, molTree.get_child("AttachmentPointArray"))
            {
              break;
            }
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

    ChemicalReaction *parseReaction(boost::property_tree::ptree rxnTree, boost::property_tree::ptree documentTree)
    {
      ChemicalReaction *rxn = NULL;
      MarvinReaction *marvinReaction = NULL;
      RWMol *mol = NULL;
      
      try
      {

        rxn = new ChemicalReaction();

        marvinReaction= parseMarvinReaction(rxnTree, documentTree);
        //printf("%s\n", marvinReaction->toString().c_str());
        marvinReaction->convertSuperAtoms();
        //printf("%s\n", marvinReaction->toString().c_str());

        // get each reactant

        std::vector<MarvinMol *>::iterator molIter;
        for (molIter = marvinReaction->reactants.begin(); molIter != marvinReaction->reactants.end(); ++molIter)
        {
          mol = parseMolecule((*molIter) , true);
          ROMol *roMol = new ROMol(*mol);
          delete mol;

          rxn->addReactantTemplate(ROMOL_SPTR(roMol));   //roMol  now owned by rxn;
        }

       // get each agent

        for (molIter = marvinReaction->agents.begin(); molIter != marvinReaction->agents.end(); ++molIter)
        {
          mol = parseMolecule((*molIter) , true);
          ROMol *roMol = new ROMol(*mol);
          delete mol;
          
          rxn->addAgentTemplate(ROMOL_SPTR(roMol)); //roMol  now owned by rxn;
        }

        // get each product

        for (molIter = marvinReaction->products.begin(); molIter != marvinReaction->products.end(); ++molIter)
        {
          mol = parseMolecule((*molIter) , true);
          ROMol *roMol = new ROMol(*mol);
          delete mol;

          rxn->addProductTemplate(ROMOL_SPTR(roMol)); //roMol  now owned by rxn;
        }

        // convert atoms to queries:
        for (MOL_SPTR_VECT::const_iterator iter = rxn->beginReactantTemplates();
            iter != rxn->endReactantTemplates(); ++iter) {
          // to write the mol block, we need ring information:
          for (ROMol::AtomIterator atomIt = (*iter)->beginAtoms();
              atomIt != (*iter)->endAtoms(); ++atomIt) {
            QueryOps::replaceAtomWithQueryAtom((RWMol *)iter->get(), (*atomIt));
          }
        }
        for (MOL_SPTR_VECT::const_iterator iter = rxn->beginProductTemplates();
            iter != rxn->endProductTemplates(); ++iter) {
          // to write the mol block, we need ring information:
          for (ROMol::AtomIterator atomIt = (*iter)->beginAtoms();
              atomIt != (*iter)->endAtoms(); ++atomIt) {
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
       
      

        BOOST_FOREACH(
              boost::property_tree::ptree::value_type &v, rxnTree.get_child("reactantList"))
            res->reactants.push_back((MarvinMol *)parseMarvinMolecule(v.second));

        BOOST_FOREACH(
              boost::property_tree::ptree::value_type &v, rxnTree.get_child("agentList"))
            res->agents.push_back((MarvinMol *)parseMarvinMolecule(v.second));

        BOOST_FOREACH(
              boost::property_tree::ptree::value_type &v, rxnTree.get_child("productList"))
            res->products.push_back((MarvinMol *)parseMarvinMolecule(v.second));

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
              if (!getCleanDouble( v.second.get<std::string>("<xmlattr>.x", ""),x)
                  || !getCleanDouble( v.second.get<std::string>("<xmlattr>.y", ""),y))
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
          int fontScale;
          if (! getCleanInt(v.second.get<std::string>("<xmlattr>.fontScale", ""),  fontScale))
            throw FileParseException("Condition font scale must be an integersin MRV file"); 
          marvinCondition->fontScale = fontScale;
          marvinCondition->text = v.second.get<std::string>("Field", "");

          int pointCount = 0;
          BOOST_FOREACH(
              boost::property_tree::ptree::value_type &v2, v.second)
          {
            if (v2.first == "MPoint")
            {
              int x,y;
              if (! getCleanInt(v.second.get<std::string>("<xmlattr>.x", ""),  x)
              ||  ! getCleanInt(v.second.get<std::string>("<xmlattr>.x", ""),  y))
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
  
  void *MrvDataStreamParser(std::istream *inStream, bool &isReaction)
  {
    PRECONDITION(inStream, "no stream");
    std::string tempStr;
    void *res = NULL;
    MarvinCML marvinCML;

    res = marvinCML.parse(*inStream, isReaction);

    // if (marvinCML.isReaction())
    //   res = (void *)marvinCML.rxn;
    // else
    //   res = (void *)marvinCML.mol;

    return res;
  }

  void *MrvDataStreamParser(std::istream &inStream, bool &isReaction)
  {
    return MrvDataStreamParser(&inStream, isReaction);
  }
  //------------------------------------------------
  //
  //  Read a molecule from a string
  //
  //------------------------------------------------
  void *MrvStringParser(const std::string &molmrvText, bool &isReaction)
  {
    std::istringstream inStream(molmrvText);
    // unsigned int line = 0;
    return MrvDataStreamParser(inStream, isReaction);
  }

  //------------------------------------------------
  //
  //  Read a RWMOL from a file
  //
  //------------------------------------------------
  void *MrvFileParser(const std::string &fName, bool &isReaction)
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
      res = MrvDataStreamParser(inStream,  isReaction);
    }
    return res;
  }

  //------------------------------------------------
  //
  //  Read a RWMol from a stream 
  //
  //------------------------------------------------
  RWMol *MrvMolDataStreamParser(std::istream *inStream)
  {
  
    void *res = NULL;
    
    bool isReaction;
    res = MrvDataStreamParser(inStream, isReaction);
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
  RWMol *MrvMolDataStreamParser(std::istream &inStream)
  {
    return MrvMolDataStreamParser(&inStream);
  }
  //------------------------------------------------
  //
  //  Read a RWMol from a string
  //
  //------------------------------------------------
  RWMol *MrvMolStringParser(const std::string &molmrvText)
  {
    std::istringstream inStream(molmrvText);
    // unsigned int line = 0;
    return MrvMolDataStreamParser(inStream);
  }

  //------------------------------------------------
  //
  //  Read an RWMol from a file
  //
  //------------------------------------------------
  RWMol *MrvMolFileParser(const std::string &fName)
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
      res = MrvMolDataStreamParser(inStream);
    return res;
  }

  //------------------------------------------------
  //
  //  Read a ChemicalReaction from a stream 
  //
  //------------------------------------------------
  ChemicalReaction *MrvRxnDataStreamParser(std::istream *inStream)
  {
  
    void *res = NULL;
    
    bool isReaction;
    res = MrvDataStreamParser(inStream, isReaction);
    if (!isReaction)
    {
        delete (RWMol *)res;
        throw FileParseException("The file parsed as a reaction, not a molecule"); 
    }

    return (ChemicalReaction *)res;
  }

  //------------------------------------------------
  //
  //  Read a ChemicalReaction from a stream reference
  //
  //------------------------------------------------
  ChemicalReaction *MrvRxnDataStreamParser(std::istream &inStream)
  {
    return MrvRxnDataStreamParser(&inStream);
  }
  //------------------------------------------------
  //
  //  Read a ChemicalReaction from a string
  //
  //------------------------------------------------
  ChemicalReaction *MrvRxnStringParser(const std::string &molmrvText)
  {
    std::istringstream inStream(molmrvText);
    // unsigned int line = 0;
    return MrvRxnDataStreamParser(inStream);
  }

  //------------------------------------------------
  //
  //  Read a ChemicalReaction from a file
  //
  //------------------------------------------------
  ChemicalReaction *MrvRxnFileParser(const std::string &fName)
  {
    std::ifstream inStream(fName.c_str());
    if (!inStream || (inStream.bad()))
    {
      std::ostringstream errout;
      errout << "Bad input file " << fName;
      throw BadFileException(errout.str());
    }
    ChemicalReaction *res = MrvRxnDataStreamParser(inStream);
    return res;
  }
}// namespace RDKit
