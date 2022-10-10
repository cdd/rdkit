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

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <RDGeneral/BoostEndInclude.h>


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

        std::string toString() const;        
    };

    class MarvinPlus
    {
        public:
        std::string id;
        double x1;
        double y1;
        double x2;
        double y2;


        std::string toString() const;
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

        std::string toString() const;       
    };

    
    class MarvinAttachmentPoint
    {
        public:

        // <AttachmentPoint atom="a7" order="1" bond="b6"/>
        std::string atom;
        std::string bond;
        std::string order;

        std::string toString() const;
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

        MarvinAtom();
    
        bool operator==(const MarvinAtom& rhs) const;
       
        bool operator==(const MarvinAtom *rhs) const;

        std::string toString() const;
    };

    class MarvinBond
    {
        public:
        std::string id;
        std::string atomRefs2[2];
        std::string order;
        std::string bondStereo;
        std::string queryType;
    
        bool isEqual(const MarvinAtom& other) const;
       

        bool operator==(const MarvinAtom& rhs) const;
        
        std::string toString() const;
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

        virtual ~MarvinMolBase();
              
        int getAtomIndex(std::string id);
      
        int getBondIndex(std::string id);       
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

        std::string toString() const;
       
        std::string role();
       
    };

    
    class MarvinSuperatomSgroup : public MarvinMolBase
    {
        public:
        std::string id;
        std::string title;
        std::vector<MarvinAttachmentPoint *> attachmentPoints;

     
      ~MarvinSuperatomSgroup();
      
      std::string role();
    
      std::string toString() const;

    };

    class MarvinSuperInfo   // used in converting superatoms group to mol style groups
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

      ~MarvinMol();
    
      std::string role();
    
      static bool atomRefInAtoms(MarvinAtom *a, std::string b );

      static bool bondRefInBonds(MarvinBond *a, std::string b );
    
      void cleanUpNumbering(int &molCount    // this is the starting mol count, and receives the ending mol count - THis is used when MarvinMol->convertToSuperAtaoms is called multiple times from a RXN
        , int &atomCount  // starting and ending atom count
        , int &bondCount  // starting and ending bond count
        , int &sgCount);  // starting and ending sg  count)
      
      void convertFromSuperAtoms();
      
      void convertToSuperAtoms();
       
      std::string toString() const;
      
      std::string generateMolString();   
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

      ~MarvinReaction();
      
      void convertFromSuperAtoms();
      

      void convertToSuperAtoms();
     
      std::string toString();
    };

    class MarvinRectangle
    {
      private:
      RDGeom::Point3D *center;
      
      public:
      RDGeom::Point3D upperLeft;
      RDGeom::Point3D lowerRight;

      MarvinRectangle(double left, double right, double top, double bottom);
      MarvinRectangle(const RDGeom::Point3D &upperLeftInit, const RDGeom::Point3D &lowerRightInit);
      MarvinRectangle(const std::vector<MarvinAtom *> atoms);
      
      void extend(const MarvinRectangle &otherRectangle);
     
      RDGeom::Point3D getCenter();
      
      bool overlapsVertically(const MarvinRectangle &otherRectangle) const;
      
      bool overlapsVHorizontally(const MarvinRectangle &otherRectangle) const;

      static bool compareRectanglesByX(MarvinRectangle &r1, MarvinRectangle &r2);
     
      static bool compareRectanglesByY(MarvinRectangle &r1, MarvinRectangle &r2);
    };

    class MarvinStereoGroup
    {
      public:

      StereoGroupType groupType;  // one of ABS AND OR
      int groupNumber;
      std::vector<unsigned int> atoms;

      MarvinStereoGroup(StereoGroupType grouptypeInit, int groupNumberInit );
    };


}

#endif //RD_MARVINDEFS_H

