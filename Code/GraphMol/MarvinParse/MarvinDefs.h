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

  const std::vector<std::string> sruSgroupConnectChoices{"hh", "ht", "eu"};
  const std::vector<std::string> marvinBondOrders{"1", "2", "3", "A"};
  const std::vector<std::string> marvinQueryBondsTypes{"SD", "SA", "DA", "Any"};
  const std::vector<std::string> marvinConventionTypes{"cxn:coord"};
  const std::vector<std::string> marvinRadicalVals{"monovalent", "divalent", "divalent1", "divalent3", "trivalent","trivalent2","trivalent4", "4"};
  const std::map<std::string, int> marvinRadicalToRadicalElectrons{
      {"monovalent", 1},
      {"divalent", 2},
      {"divalent1", 2},
      {"divalent3", 2},
      {"trivalent", 3},
      {"trivalent2", 3},
      {"trivalent4", 3}
      };

  const std::map<int, std::string> radicalElectronsToMarvinRadical{
      {1, "monovalent"}, {2, "divalent"}, {3, "trivalent4"}, {4, "4"}};

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
    double fontScale=0.0;

    std::string halign;
    std::string valign;

    std::string toString() const;       
  };


  class MarvinAttachmentPoint
  {
    public:

    // <attachmentPoint atom="a7" order="1" bond="b6"/>
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
    int mrvValence;
    int hydrogenCount;
    std::string mrvAlias;
    std::string mrvStereoGroup;
    int mrvMap;
    std::string sgroupRef;
    std::string sgroupAttachmentPoint;
    int rgroupRef;

    MarvinAtom();
    MarvinAtom(const MarvinAtom &atomToCopy, std::string newId);


    bool operator==(const MarvinAtom& rhs) const;
    
    bool operator==(const MarvinAtom *rhs) const;

    bool isElement() const;

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
    std::string convention;

    MarvinBond()
    {
    }

    MarvinBond(const MarvinBond &bondToCopy, std::string newId, std::string atomRef1, std::string atomRef2);

    bool isEqual(const MarvinAtom& other) const;
    
    bool operator==(const MarvinAtom& rhs) const;

    const std::string getBondType() const;

    std::string toString() const;
  };

  class MarvinRectangle
  {
    protected:
    RDGeom::Point3D center;
    bool centerIsStale;
    
    public:
    RDGeom::Point3D upperLeft;
    RDGeom::Point3D lowerRight;

    MarvinRectangle();
    MarvinRectangle(double left, double right, double top, double bottom);
    MarvinRectangle(const RDGeom::Point3D &upperLeftInit, const RDGeom::Point3D &lowerRightInit);
    MarvinRectangle(const std::vector<MarvinAtom *> atoms);
    MarvinRectangle(const std::vector<MarvinRectangle> rects);
    
    void extend(const MarvinRectangle &otherRectangle);
    
    RDGeom::Point3D &getCenter();
    
    bool overlapsVertically(const MarvinRectangle &otherRectangle) const;
    
    bool overlapsVHorizontally(const MarvinRectangle &otherRectangle) const;

    static bool compareRectanglesByX(MarvinRectangle &r1, MarvinRectangle &r2);
    
    static bool compareRectanglesByYReverse(MarvinRectangle &r1, MarvinRectangle &r2);
  };

  class MarvinBracket : public MarvinRectangle
  {
    public:

    MarvinBracket() {};
    MarvinBracket(double left, double right, double top, double bottom) : MarvinRectangle(left, right, top, bottom) {};
    MarvinBracket(const RDGeom::Point3D &upperLeftInit, const RDGeom::Point3D &lowerRightInit) : MarvinRectangle(upperLeftInit, lowerRightInit) {} ;
    MarvinBracket(const std::vector<MarvinAtom *> atoms) : MarvinRectangle(atoms) {};

    std::string toString() const
    {
      std::ostringstream out;

      out << "<MBracket type=\"SQUARE\" orientation=\"DOUBLE\"><MPoint x=\""
      << upperLeft.x << "\" y=\"" << upperLeft.y << "\"></MPoint><MPoint x=\""
      << lowerRight.x << "\" y=\"" << upperLeft.y << "\"></MPoint><MPoint x=\""
      << upperLeft.x << "\" y=\"" << lowerRight.y << "\"></MPoint><MPoint x=\""
      << lowerRight.x <<"\" y=\"" << lowerRight.y << "\"></MPoint></MBracket>";

      return out.str();
    }
  };

  class MarvinMolBase
  {
    public:
    std::string molID;
    std::vector<MarvinAtom *> atoms;
    std::vector<MarvinBond *> bonds;    
    MarvinBracket bracket;    //  only used for derived classes that do not have their own atoms 
  
    virtual std::string role() const = 0;
    virtual bool hasAtomBondBlocks() const = 0;


    int getExplicitValence(const MarvinAtom &marvinAtom) const;

    MarvinMolBase()
    {
    }

    virtual ~MarvinMolBase();
          
    int getAtomIndex(std::string id);
  
    int getBondIndex(std::string id);  

    const std::vector<std::string> getBondList() const;
    const std::vector<std::string> getAtomList() const;
    bool AnyOverLappingAtoms(const MarvinMolBase *otherMol) const;

    bool hasCoords() const;  
    void removeCoords();
  };

  class MarvinSruSgroup : public MarvinMolBase
  {
    public:
    std::string id;
    std::string title;
    std::string connect;
    std::string correspondence;

    std::string toString() const;
    
    std::string role() const;    
    bool hasAtomBondBlocks() const;
  };

  class MarvinDataSgroup : public MarvinMolBase
  {
    public:
    std::string id;
    std::string context;
    std::string fieldName;
    std::string placement;
    std::string unitsDisplayed;
    std::string queryType;
    std::string queryOp;
    std::string fieldData;
    std::string units;
    double x;
    double y;

    std::string toString() const;
    
    std::string role() const;    
    bool hasAtomBondBlocks() const;
  };



  class MarvinSuperatomSgroupExpanded : public MarvinMolBase
  {
    public:
    std::string id;
    std::string title;

    std::vector<MarvinAttachmentPoint *> attachmentPoints;
    
    ~MarvinSuperatomSgroupExpanded();
    
    std::string toString() const;
    
    std::string role() const;
    bool hasAtomBondBlocks() const;
  };

  class MarvinMultipleSgroup : public MarvinMolBase
  {
    public:

    std::string id;
    std::string title;
    bool isExpanded = false;
    std::vector<MarvinAtom *>parentAtoms; 
    std::vector<MarvinBond *>bondsToAtomsNotInExpandedGroup;      // only when expanded

    std::string toString() const;
    
    std::string role() const;
    bool hasAtomBondBlocks() const;
  };

  class MarvinMulticenterSgroup : public MarvinMolBase
  {
    // <molecule molID="m2" id="sg1" role="MulticenterSgroup" atomRefs="a2 a6 a5 a4 a3" center="a18"/>
    public:
    std::string id;
        
    std::string toString() const;
    MarvinAtom *center;
    std::string role() const;
    bool hasAtomBondBlocks() const;
  };

  class MarvinGenericSgroup : public MarvinMolBase
  {
    // <molecule molID="m2" id="sg1" role="GenericSgroup" atomRefs="a1 a2 a3 a4 a5 a6 a7 a8 a9 a13 a10 a11 a12" charge="onAtoms"/></molecule>
    public:
    std::string id;
    std::string charge;   // onAtoms or onBrackets
    std::string toString() const;
    std::string role() const;
    bool hasAtomBondBlocks() const;
  };

  class MarvinMonomerSgroup : public MarvinMolBase
  {
      // <molecule id="sg1" role="MonomerSgroup" title="mon" charge="onAtoms" molID="m2" atomRefs="a2 a1 a3 a4">
      //     <MBracket type="SQUARE" orientation="DOUBLE">
      //         <MPoint x="-0.8726666666666667" y="1.078"></MPoint>
      //         <MPoint x="1.2833333333333334" y="1.078"></MPoint>
      //         <MPoint x="1.2833333333333334" y="-1.078"></MPoint>
      //         <MPoint x="-0.8726666666666667" y="-1.078"></MPoint>
      //     </MBracket>
      // </molecule> 
    public:
    std::string id;
    std::string title;
    std::string charge;   // onAtoms or onBrackets
    std::string toString() const;
    std::string role() const;
    bool hasAtomBondBlocks() const;
  };

  class MarvinSuperatomSgroup : public MarvinMolBase
  {
    public:
    std::string id;
    std::string title;
    std::vector<MarvinAttachmentPoint *> attachmentPoints;
    
    ~MarvinSuperatomSgroup();
    
    std::string role() const;
    bool hasAtomBondBlocks() const;
  
    std::string toString() const;
  };

  class MarvinSuperInfo   // used in converting superatoms group to mol style groups
  {
    public:
    std::string title;
    std::vector<std::string> atoms;
    std::vector<std::string> bonds;
  };
  
  class MarvinMol : public MarvinMolBase
  {
    public:

    std::vector<MarvinSuperatomSgroup *> superatomSgroups;
    std::vector<MarvinSruSgroup *>  sruSgroups;
    std::vector<MarvinSuperatomSgroupExpanded *>  superatomSgroupsExpanded;
    std::vector<MarvinMultipleSgroup *>  multipleSgroups;
    std::vector<MarvinDataSgroup *> dataSgroups;
    std::vector<MarvinMulticenterSgroup *> multicenterSgroups;
    std::vector<MarvinGenericSgroup *>genericSgroups;
    std::vector<MarvinMonomerSgroup *>monomerSgroups;
    std::vector<MarvinSuperInfo *> superInfos;  // used in converting superatomSgroups to mol-type CTs


    ~MarvinMol();
  
    std::string role() const;
    bool hasAtomBondBlocks() const;

    static bool atomRefInAtoms(MarvinAtom *a, std::string b );

    static bool bondRefInBonds(MarvinBond *a, std::string b );
  
    void cleanUpNumbering(int &molCount    // this is the starting mol count, and receives the ending mol count - THis is used when MarvinMol->convertToSuperAtaoms is called multiple times from a RXN
      , int &atomCount  // starting and ending atom count
      , int &bondCount  // starting and ending bond count
      , int &sgCount);  // starting and ending sg  count)
    
    void convertFromSuperAtoms();
    
    void convertToSuperAtoms();

    void processMulticenterSgroups();
    void expandMultipleSgroups();
    void contractMultipleSgroups();
      
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

    void expandMultipleSgroups();
    

    void processMulticenterSgroups();
    
    std::string toString();
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
