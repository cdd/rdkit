//
//  Copyright (C) 2002-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/BoostStartInclude.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/format.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <string>
#include <set>
#include <exception>
#include <iostream>

namespace pt = boost::property_tree;

#include "FileParsers.h"
#include "FileParserUtils.h"
#include "MolSGroupParsing.h"
#include "MolFileStereochem.h"

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/StereoGroup.h>
#include <GraphMol/SubstanceGroup.h>
#include <RDGeneral/StreamOps.h>
#include <RDGeneral/RDLog.h>

#include <fstream>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/LocaleSwitcher.h>
#include <typeinfo>
#include <exception>
#include <charconv>

#ifdef RDKIT_USE_BOOST_REGEX
#include <boost/regex.hpp>
using boost::regex;
using boost::regex_match;
using boost::smatch;
#else
#include <regex>
using std::regex;
using std::regex_match;
using std::smatch;
#endif
#include <sstream>
#include <locale>
#include <cstdlib>
#include <cstdio>
#include <string_view>
#include <list>

using namespace RDKit::SGroupParsing;

namespace RDKit
{

  /*
    Imports the Marvin-specific dialect of CML (Chemical Markup Language) and converts it to datastructures
    that are compatible with Molfile, RXNfile, and Molfile complemented with canvas objects.
  */

  class MarvinArrow
  {
  public:
    double x1;
    double y1;
    double x2;
    double Y2;
  };

  class MarvinPlus
  {
  public:
    double x;
    double y;
  };

  class MarvinText
  {
    std::string text;
    double x;
    double y;
    double fontScale;

    std::string halign;
    std::string valign;
  };

  class Watermark
  {
    double molID;
    double atomID;
    double bondID;
    double sgroupID;
    double objectID;
  };

  // const safeAttr = (el: Element, attr: string): string => el.getAttribute(attr) ?? '';

  class MarvinCML
  {
    private:
      std::string xml;

    public:
      RWMol *mol;

    std::list<MarvinText> textBoxes;

    std::list<RWMol *> reactants;
    std::list<RWMol *> agents;
    std::list<RWMol *> products;

    MarvinArrow arrow;

    std::list<MarvinPlus> plusSigns;

    MarvinCML()
        : mol(NULL)
        // ,reactants(new std::list<RWMol>())
        // ,agents(new std::list<RWMol>())
        // ,products(new std::list<RWMol>())
    {

    };

    void parse(std::istream &is)
    {
      // Create empty property tree object
      using boost::property_tree::ptree;
      ptree tree;

      // Parse the XML into the property tree.
      //pt::read_xml(filename, tree);

      read_xml(is, tree);

      // loop over all mols in the record

      BOOST_FOREACH(
          boost::property_tree::ptree::value_type &mol, tree.get_child("cml.MDocument.MChemicalStruct"))
      {
          std::string id = mol.second.get<std::string>("<xmlattr>.molID", "0");
          printf("id: %s\n", id.c_str());
      
        BOOST_FOREACH(
            boost::property_tree::ptree::value_type &v, mol.second.get_child("atomArray"))
        {
          std::string atomId = v.second.get<std::string>("<xmlattr>.id", "0");
          printf("atomId: %s\n", atomId.c_str());
          std::string elementType = v.second.get<std::string>("<xmlattr>.elementType", "0");
          printf("elementType: %s\n", elementType.c_str());
          std::string x2 = v.second.get<std::string>("<xmlattr>.x2", "0");
          printf("x2: %s\n", x2.c_str());
          std::string y2 = v.second.get<std::string>("<xmlattr>.y2", "0");
          printf("y2: %s\n", y2.c_str());
          std::string formalCharge = v.second.get<std::string>("<xmlattr>.formalCharge", "0");
          printf("formalCharge: %s\n", formalCharge.c_str());
        }
      }

      // Use the throwing version of get to find the debug filename.
      // If the path cannot be resolved, an exception is thrown.
      // m_file = tree.get<std::string>("debug.filename");

      // Use the default-value version of get to find the debug level.
      // Note that the default value is used to deduce the target type.
      // m_level = tree.get("debug.level", 0);

      // Use get_child to find the node containing the modules, and iterate over
      // its children. If the path cannot be resolved, get_child throws.
      // A C++11 for-range loop would also work.
      // BOOST_FOREACH (pt::ptree::value_type &v, tree.get_child("debug.modules"))
      // {
      //   // The data function is used to access the data stored in a node.
      //   m_modules.insert(v.second.data());
      // }

      // qqq this.xml = WebMolKit.XML.parseXML(this.cml);
      // if (this.xml == null)
      //  throw new Error('Unable to read CML document, XML parsing error.');

      //     const heads : Element[] = [];
      //     this.recursiveHeadFind(heads, this.xml.documentElement);

      //     // look for text labels
      //     for (const head of heads)
      //     {
      //       if (head.tagName == 'MTextBox')
      //       {
      //         this.extractTextBox(head);
      //       }
      //     }

      //     // if there are any reaction units, parse the first one and call it complete
      //     for (const head of heads)
      //     {
      //       if (head.tagName == 'reaction')
      //       {
      //         this.parseReaction(head);

      //         this.plusSigns = [];
      //         for (const head of heads)
      //         {
      //           if (head.tagName == 'MReactionSign')
      //           {
      //             this.extractPlusSign(head);
      //           }
      //         }
      //         this.postProcessPluses();

      //         return;
      //       }
      //     }

      //     // read all of the molecule(s) into a single object
      //     for (const head of heads)
      //     {
      //       const frag = this.parseMolecule(head);
      //       if (this.mol)
      //       {
      //         this.mol.append(frag);
      //       }
      //       else
      //       {
      //         this.mol = frag;
      //       }
      //     }

      //     if (!this.mol)
      //       throw new Error('No molecules or reactions found.');
    }

      // public
      //   getMolfile() : string
      //   {
      //     return new WebMolKit.MDLMOLWriter(this.mol).write();
      //   }

      // public
      //   isReaction() : boolean
      //   {
      //     return this.reactants != null && this.agents != null && this.products != null;
      //   }

      // public
      //   getRXNfile() : string
      //   {
      //     if (!this.isReaction())
      //       throw new Error('Source content is not a reaction.');
      //     const pad = (v
      //                  : number) = >
      //     {
      //       const str = v.toString();
      //       return ' '.repeat(3 - str.length) + str;
      //     };
      //     let rxn = '$RXN\n\nRXNfile\n\n' + pad(this.reactants.length) + pad(this.products.length) + pad(this.agents.length) + '\n';
      //     for (const mol of this.reactants)
      //     {
      //       rxn += '$MOL\n' + new WebMolKit.MDLMOLWriter(mol).write() + '\n';
      //     }
      //     for (const mol of this.products)
      //     {
      //       rxn += '$MOL\n' + new WebMolKit.MDLMOLWriter(mol).write() + '\n';
      //     }
      //     for (const mol of this.agents)
      //     {
      //       rxn += '$MOL\n' + new WebMolKit.MDLMOLWriter(mol).write() + '\n';
      //     }
      //     return rxn;
      //   }

      // private
      //   recursiveHeadFind(heads
      //                     : Element[], parent
      //                     : Element) : void
      //   {
      //     for (const child of WebMolKit.XML.childElements(parent))
      //     {
      //       if ([ 'molecule', 'reaction', 'MTextBox', 'MReactionSign' ].includes(child.tagName))
      //       {
      //         heads.push(child);
      //       }
      //       else
      //       {
      //         this.recursiveHeadFind(heads, child);
      //       }
      //     }
      //   }

      //   // extracts a single molecule from the element; for CML sources containing just a molecule, this is the primary workhorse; it
      //   // also gets called from the datasheet-extracting methods
      // private
      //   parseMolecule(elMol
      //                 : Element, attpoint
      //                 : string[] = null) : WebMolKit.Molecule
      //   {
      //     const mol = new WebMolKit.Molecule();
      //     const atomID : string[] = [];
      //     const sgroupID : string[] = [];

      //     // scan for atoms first
      //     for (const elAtomArray of WebMolKit.XML.findChildElements(elMol, 'atomArray'))
      //     {
      //       this.extractAtomOrArray(elAtomArray, mol, atomID, sgroupID, attpoint);
      //       for (const elAtom of WebMolKit.XML.findChildElements(elAtomArray, 'atom'))
      //       {
      //         this.extractAtomOrArray(elAtom, mol, atomID, sgroupID, attpoint);
      //       }
      //     }

      //     // scan for bonds (best to do it separately, since they refer to atoms by ID)
      //     for (const elBondArray of WebMolKit.XML.findChildElements(elMol, 'bondArray'))
      //     {
      //       for (const elBond of WebMolKit.XML.findChildElements(elBondArray, 'bond'))
      //       {
      //         this.extractBond(elBond, mol, atomID);
      //       }
      //     }

      //     // scan for sub-molecules that can be converted into abbreviations
      //     for (const elSubMol of WebMolKit.XML.findChildElements(elMol, 'molecule'))
      //     {
      //       const role = elSubMol.getAttribute('role');
      //       if (role == 'SuperatomSgroup')
      //       {
      //         const sgID = elSubMol.getAttribute('id');
      //         if (!sgID)
      //           continue;
      //         const headAtom = sgroupID.indexOf(sgID) + 1;
      //         if (headAtom > 0)
      //           this.affixAbbreviation(elSubMol, mol, headAtom);
      //       }
      //       else if (role == 'SruSgroup')
      //       {
      //         const atomList = safeAttr(elSubMol, 'atomRefs').split(' ').map((aid) = > atomID.indexOf(aid) + 1);
      //         if (atomList.some((atom) = > atom == 0))
      //           continue;
      //         const unit = new WebMolKit.PolymerBlockUnit(atomList);
      //         const connect = elSubMol.getAttribute('connect');
      //         if (connect == 'ht')
      //         {
      //           unit.connect = WebMolKit.PolymerBlockConnectivity.HeadToTail;
      //         }
      //         else if (connect == 'hh')
      //         {
      //           unit.connect = WebMolKit.PolymerBlockConnectivity.HeadToHead;
      //         }
      //         else if (connect == 'eu')
      //         {
      //           unit.connect = WebMolKit.PolymerBlockConnectivity.Random;
      //         }
      //         new WebMolKit.PolymerBlock(mol).createUnit(unit);
      //       }
      //       else if (role == 'MultipleSgroup')
      //       {
      //         /*
      //           ignoring for now: there's a story for implementing this; it needs to be stored using an intermediate form of
      //           "polymer" representation with the # of units retained; exporting this to & from Ketcher is easy enough, but
      //           going to Molfile format requires actual expansion of the atoms, then annotation of which ones were the original,
      //           which is a bit more tricky

      //         const atomList = safeAttr(elSubMol, 'atomRefs').split(' ').map((aid) => atomID.indexOf(aid) + 1);
      //         const count = parseInt(safeAttr(elSubMol, 'title'));
      //         const WITHIN_REASON = 30; // duplicating atoms, don't want to do anything crazy
      //         for (let n = 1; n < count && n < WITHIN_REASON && mol.numAtoms + count < 1000; n++) {
      //           this.replicateMultiBlock(mol, atomList);
      //         }
      //         */
      //       }
      //     }

      //     return mol;
      //   }

      //   // pulls out the content of an <atomArray> or <atom> node
      // private
      //   extractAtomOrArray(elAtom
      //                      : Element, mol
      //                      : WebMolKit.Molecule, atomID
      //                      : string[], sgroupID
      //                      : string[], attpoint
      //                      : string[]) : void
      //   {
      //     const id = elAtom.getAttribute('atomID') || elAtom.getAttribute('id');
      //     if (!id)
      //       return;

      //     const bitsID = id.split(' ');
      //     const bitsElement = safeAttr(elAtom, 'elementType').split(' ');
      //     const bitsX = safeAttr(elAtom, 'x2').split(' ');
      //     const bitsY = safeAttr(elAtom, 'y2').split(' ');
      //     const bitsCharge = safeAttr(elAtom, 'formalCharge').split(' ');
      //     const bitsRadical = safeAttr(elAtom, 'radical').split(' ');
      //     const bitsIsotope = safeAttr(elAtom, 'isotope').split(' ');
      //     const bitsSGroup = safeAttr(elAtom, 'sgroupRef').split(' ');
      //     const bitsAttPoint = (safeAttr(elAtom, 'attachmentPoint') || elAtom.getAttribute('sgroupAttachmentPoint') || '').split(' ');
      //     const bitsMapNum = safeAttr(elAtom, 'mrvMap').split(' ');
      //     const bitsQuery = safeAttr(elAtom, 'mrvQueryProps').split(' ');
      //     const bitsAlias = safeAttr(elAtom, 'mrvAlias').split(' ');
      //     for (let n = 0; n < bitsID.length; n++)
      //     {
      //       atomID.push(bitsID[n]);
      //       sgroupID.push(n < bitsSGroup.length ? bitsSGroup[n] : '');
      //       const x = Math.round(WebMolKit.safeFloat(bitsX[n]) * 10000) * 0.0001;
      //       const y = Math.round(WebMolKit.safeFloat(bitsY[n]) * 10000) * 0.0001;
      //       const a = mol.addAtom(bitsElement[n], x, y);
      //       if (n < bitsCharge.length)
      //         mol.setAtomCharge(a, WebMolKit.safeInt(bitsCharge[n]));
      //       if (bitsRadical[n] == 'monovalent')
      //       {
      //         mol.setAtomUnpaired(a, 1);
      //       }
      //       else if (bitsRadical[n] == 'divalent' || bitsRadical[n] == 'divalent1')
      //       {
      //         mol.setAtomUnpaired(a, 2);
      //       }
      //       else if (bitsRadical[n] == 'trivalent4')
      //       {
      //         mol.setAtomUnpaired(a, 3);
      //       }
      //       else if (bitsRadical[n] == '4')
      //       {
      //         mol.setAtomUnpaired(a, 4);
      //       }
      //       if (n < bitsIsotope.length)
      //       {
      //         const iso = parseInt(bitsIsotope[n]);
      //         if (iso > 0)
      //           mol.setAtomIsotope(a, iso);
      //       }
      //       if (n < bitsMapNum.length)
      //       {
      //         const mapnum = parseInt(bitsMapNum[n]);
      //         if (mapnum > 0)
      //           mol.setAtomMapNum(a, mapnum);
      //       }
      //       if (n < bitsQuery.length && bitsQuery[n] != null && bitsQuery[n] != '' && bitsQuery[n] != '0')
      //       {
      //         mol.setAtomElement(a, 'A');
      //       }
      //       if (n < bitsAlias.length && bitsAlias[n] != null && bitsAlias[n] != '' && bitsAlias[n] != '0')
      //       {
      //         mol.setAtomElement(a, bitsAlias[n]);
      //       }
      //       if (attpoint != null)
      //         attpoint.push(bitsAttPoint[n]); // (push'd value can be null)
      //     }

      //     // optionally keep other stuff, since CML has no limits on what can be included
      //   }

      //   // pulls out the content of a <bondArray> or <bond> node
      // private
      //   extractBond(elBond
      //               : Element, mol
      //               : WebMolKit.Molecule, atomID
      //               : string[]) : void
      //   {
      //     if (!elBond.hasAttribute('atomRefs2'))
      //       return;

      //     const aref = elBond.getAttribute('atomRefs2').split(' ');
      //     const a1 = atomID.indexOf(aref[0]), a2 = atomID.indexOf(aref[1]);
      //     if (a1 < 0 || a2 < 0)
      //       return;
      //     const strOrder = safeAttr(elBond, 'order');
      //     const order = WebMolKit.safeInt(strOrder) || 0;

      //     const elStereo = WebMolKit.XML.findElement(elBond, 'bondStereo');
      //     let btype = WebMolKit.Molecule.BONDTYPE_NORMAL;
      //     if (elStereo != null)
      //     {
      //       let stereo = WebMolKit.XML.nodeText(elStereo);
      //       if (stereo == 'W')
      //         btype = WebMolKit.Molecule.BONDTYPE_INCLINED;
      //       else if (stereo == 'H')
      //         btype = WebMolKit.Molecule.BONDTYPE_DECLINED;

      //       if (elStereo.hasAttribute('dictRef'))
      //       {
      //         stereo = elStereo.getAttribute('dictRef');
      //         if (stereo == 'cml:W')
      //           btype = WebMolKit.Molecule.BONDTYPE_INCLINED;
      //         else if (stereo == 'cml:H')
      //           btype = WebMolKit.Molecule.BONDTYPE_DECLINED;
      //       }
      //       if (elStereo.getAttribute('convention') == 'MDL')
      //       {
      //         const mdl = WebMolKit.safeInt(elStereo.getAttribute('conventionValue'));
      //         if (mdl == 1)
      //           btype = WebMolKit.Molecule.BONDTYPE_INCLINED;
      //         else if (mdl == 6)
      //           btype = WebMolKit.Molecule.BONDTYPE_DECLINED;
      //         else if (mdl == 3 || mdl == 4)
      //           btype = WebMolKit.Molecule.BONDTYPE_UNKNOWN;
      //       }
      //     }

      //     const b = mol.addBond(a1 + 1, a2 + 1, order, btype);
      //     if (strOrder == 'A')
      //     {
      //       WebMolKit.QueryUtil.setQueryBondOrders(mol, b, [-1]);
      //     }
      //   }

      //   // links an atom to a child-element that contains an Sgroup molecule; super groups can be defined as MMI-flavour, which
      //   // is round-trip compatible with the abbreviation mechanism, or ChemAxon-flavour, which is more similar to the MDL approach
      // private
      //   affixAbbreviation(elAbbrev
      //                     : Element, mol
      //                     : WebMolKit.Molecule, atom
      //                     : number) : void
      //   {
      //     const attpoint : string[] = [];
      //     const smol : WebMolKit.Molecule = this.parseMolecule(elAbbrev, attpoint);
      //     if (WebMolKit.MolUtil.isBlank(smol))
      //       return;

      //     const title = elAbbrev.getAttribute('title');
      //     if (title)
      //       mol.setAtomElement(atom, title);

      //     // define the abbreviation, starting with the guide atom; use the special attribute to define its position, if available
      //     const abv = new WebMolKit.Molecule();
      //     abv.addAtom(WebMolKit.MolUtil.ABBREV_ATTACHMENT, mol.atomX(atom), mol.atomY(atom));
      //     abv.append(smol);

      //     for (let n = 0; n < attpoint.length; n++)
      //     {
      //       if (attpoint[n] == '1')
      //         abv.addBond(1, n + 2, 1);
      //     }
      //     this.projectAttachmentPoint(abv);

      //     // if no connection to guide atom, this is no good
      //     let hasConn = false;
      //     for (let n = 1; n <= abv.numBonds; n++)
      //     {
      //       if (abv.bondFrom(n) == 1 || abv.bondTo(n) == 1)
      //       {
      //         hasConn = true;
      //         break;
      //       }
      //     }
      //     if (!hasConn)
      //     {
      //       // just create the thing outright, and rotate the first new atom into the old atom's spot, so that it
      //       // doesn't disrupt the overall process too much
      //       abv.deleteAtomAndBonds(1);
      //       const idx = mol.numAtoms + 1;
      //       mol.append(abv);
      //       mol.swapAtoms(atom, idx);
      //       mol.deleteAtomAndBonds(idx);
      //       return;
      //     }

      //     if (mol.atomAdjCount(atom) == 1)
      //       WebMolKit.MolUtil.setAbbrev(mol, atom, abv);
      //   }

      //   // for an abbreviation, where the alignment atom (#1) has been fabricated on account of there being no data for it,
      //   // determine a suitable position which is based on projecting outward from the attached atoms
      // private
      //   projectAttachmentPoint(mol
      //                          : WebMolKit.Molecule) : void
      //   {
      //     let numProj = 0;
      //     let projX = 0, projY = 0;

      //     const frag = mol.clone(); // this is without the phony placeholder atom
      //     frag.deleteAtomAndBonds(1);

      //     for (let n = 2; n <= mol.numAtoms; n++)
      //     {
      //       if (mol.findBond(1, n) > 0 && mol.atomAdjCount(n) > 0)
      //       {
      //         let ang = WebMolKit.SketchUtil.calculateNewBondAngles(frag, n - 1, 1);
      //         if (ang == null)
      //           ang = WebMolKit.SketchUtil.exitVectors(frag, n - 1);
      //         if (WebMolKit.Vec.len(ang) == 0)
      //           continue;

      //         let bestScore = 0, bestX = 0, bestY = 0;
      //         for (let i = 0; i < ang.length; i++)
      //         {
      //           const x = mol.atomX(n) + WebMolKit.Molecule.IDEALBOND * Math.cos(ang[i]);
      //           const y = mol.atomY(n) + WebMolKit.Molecule.IDEALBOND * Math.sin(ang[i]);
      //           const score = ang.length == 1 ? 0 : WebMolKit.CoordUtil.congestionPoint(frag, x, y);
      //           if (i == 0 || score < bestScore)
      //           {
      //             bestScore = score;
      //             bestX = x;
      //             bestY = y;
      //           }
      //         }
      //         numProj++;
      //         projX += bestX;
      //         projY += bestY;
      //       }
      //     }

      //     if (numProj == 0)
      //       return;
      //     if (numProj > 0)
      //     {
      //       const inv = 1.0 / numProj;
      //       projX *= inv;
      //       projY *= inv;
      //     }
      //     mol.setAtomPos(1, projX, projY);
      //   }

      //   /* partial implementation, which is not quite correct; see comments in calling code
      //   private replicateMultiBlock(mol: WebMolKit.Molecule, atomList: number[]): void {
      //     const atomReps: number[] = [];
      //     for (const i of atomList) {
      //       const j = mol.addAtom(mol.atomElement(i), mol.atomX(i), mol.atomY(i), mol.atomCharge(i), mol.atomUnpaired(i));
      //       mol.setAtomHExplicit(j, mol.atomHExplicit(i));
      //       mol.setAtomIsotope(j, mol.atomIsotope(i));
      //       mol.setAtomMapNum(j, mol.atomMapNum(i));
      //       mol.setAtomExtra(j, mol.atomExtra(i));
      //       atomReps.push(j);
      //     }

      //     for (let n = 1, num = mol.numBonds; n <= num; n++) {
      //       let bfr = mol.bondFrom(n), bto = mol.bondTo(n);
      //       const i = atomList.indexOf(bfr), j = atomList.indexOf(bto);
      //       if (i < 0 && j < 0) continue;
      //       if (i >= 0) bfr = atomReps[i];
      //       if (j >= 0) bto = atomReps[j];
      //       mol.addBond(bfr, bto, mol.bondOrder(n), mol.bondType(n));
      //     }
      //   } */

      // private
      //   parseReaction(elRxn
      //                 : Element) : void
      //   {
      //     this.mol = new WebMolKit.Molecule();
      //     this.reactants = [];
      //     this.agents = [];
      //     this.products = [];

      //     for (const child of WebMolKit.XML.childElements(elRxn))
      //     {
      //       if (child.tagName == 'reactantList')
      //       {
      //         this.extractComponents(this.reactants, child);
      //       }
      //       else if (child.tagName == 'productList')
      //       {
      //         this.extractComponents(this.products, child);
      //       }
      //       else if (child.tagName == 'agentList')
      //       {
      //         this.extractComponents(this.agents, child);
      //       }
      //       else if (child.tagName == 'arrow')
      //       {
      //         const arrow : MarvinArrow = {
      //           x1 : parseFloat(safeAttr(child, 'x1')),
      //           y1 : parseFloat(safeAttr(child, 'y1')),
      //           x2 : parseFloat(safeAttr(child, 'x2')),
      //           y2 : parseFloat(safeAttr(child, 'y2')),
      //         };
      //         if (!Number.isNaN(arrow.x1) && !Number.isNaN(arrow.y1) && !Number.isNaN(arrow.x2) && !Number.isNaN(arrow.y2))
      //         {
      //           this.arrow = arrow;
      //         }
      //       }
      //     }
      //   }

      // private
      //   extractComponents(molList
      //                     : WebMolKit.Molecule[], elList
      //                     : Element) : void
      //   {
      //     for (const elMol of WebMolKit.XML.findChildElements(elList, 'molecule'))
      //     {
      //       const mol = this.parseMolecule(elMol);
      //       molList.push(mol);
      //       this.mol.append(mol);
      //     }
      //   }

      // private
      //   extractTextBox(elText
      //                  : Element) : void
      //   {
      //     const elField = WebMolKit.XML.findElement(elText, 'Field');
      //     if (!elField)
      //       return;
      //     const text = WebMolKit.XML.nodeText(elField);

      //     let cx = 0, cy = 0, num = 0;
      //     for (const child of WebMolKit.XML.findChildElements(elText, 'MPoint'))
      //     {
      //       const x = parseFloat(safeAttr(child, 'x'));
      //       const y = parseFloat(safeAttr(child, 'y'));
      //       if (Number.isNaN(x) || Number.isNaN(y))
      //         return;
      //       cx += x;
      //       cy += y;
      //       num++;
      //     }
      //     if (num == 0)
      //       return;

      //     this.textBoxes.push({
      //       text,
      //       x : cx / num,
      //       y : cy / num,
      //       fontScale : parseFloat(elText.getAttribute('fontScale')),
      //       halign : elText.getAttribute('halign'),
      //       valign : elText.getAttribute('valign'),
      //     });
      //   }

      // private
      //   extractPlusSign(elPlus
      //                   : Element) : void
      //   {
      //     let cx = 0, cy = 0, num = 0;
      //     for (const child of WebMolKit.XML.findChildElements(elPlus, 'MPoint'))
      //     {
      //       const x = parseFloat(safeAttr(child, 'x'));
      //       const y = parseFloat(safeAttr(child, 'y'));
      //       if (Number.isNaN(x) || Number.isNaN(y))
      //         return;
      //       cx += x;
      //       cy += y;
      //       num++;
      //     }
      //     if (num > 0)
      //     {
      //       this.plusSigns.push({x : cx / num, y : cy / num});
      //     }
      //   }

      //   // sometimes the provided pluses do not match the reagent/product lists, so it may be necessary to create or delete
      // private
      //   postProcessPluses() : void
      //   {
      //     const residual = this.plusSigns.slice(0);
      //     this.plusSigns = [];

      //     const isHorizontal = Math.abs(this.arrow.x2 - this.arrow.x1) > Math.abs(this.arrow.y2 - this.arrow.y1);

      //     const processSide = (molList
      //                          : WebMolKit.Molecule[]) : void = >
      //     {
      //       const cx = molList.map((mol) = > WebMolKit.Vec.sum(WebMolKit.MolUtil.arrayAtomX(mol)) / mol.numAtoms);
      //       const idx = WebMolKit.Vec.idxSort(cx);
      //       for (let n = 0; n < idx.length - 1; n++)
      //       {
      //         const i1 = idx[n], i2 = idx[n + 1];
      //         const b1 = molList[i1].boundary(), b2 = molList[i2].boundary();
      //         const x1 = b1.maxX(), x2 = b2.minX();
      //         if (x1 > x2)
      //           continue; // the molecule positions overlap, so no arrow
      //         let hit = false;
      //         for (let i = 0; i < residual.length; i++)
      //         {
      //           if (residual[i].x >= x1 && residual[i].x <= x2)
      //           {
      //             this.plusSigns.push(residual[i]);
      //             residual.splice(i, 1);
      //             hit = true;
      //             break;
      //           }
      //         }
      //         if (!hit)
      //         {
      //           const x = 0.5 * (x1 + x2);
      //           const y = isHorizontal ? 0.5 * (this.arrow.y1 + this.arrow.y2) : 0.25 * (b1.minY() + b1.maxY() + b2.minY() + b2.maxY());
      //           this.plusSigns.push({x, y});
      //         }
      //       }
      //     };

      //     processSide(this.reactants);
      //     processSide(this.products);
      //   }

      // public
      //   serialize() : string
      //   {
      //     const xml = WebMolKit.XML.parseXML(
      //         '<cml xmlns="http://www.chemaxon.com" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"' +
      //         ' version="ChemAxon file format v20.20.0, generated by vunknown"' +
      //         ' xsi:schemaLocation="http://www.chemaxon.com http://www.chemaxon.com/marvin/schema/mrvSchema_20_20_0.xsd"/>');
      //     const elDoc = WebMolKit.XML.appendElement(xml.documentElement, 'MDocument');
      //     const elChem = WebMolKit.XML.appendElement(elDoc, 'MChemicalStruct');

      //     const watermark : Watermark = {
      //       molID : 0,
      //       atomID : 0,
      //       bondID : 0,
      //       sgroupID : 0,
      //       objectID : 0,
      //     };

      //     if (!this.isReaction())
      //     {
      //       const elMol = WebMolKit.XML.appendElement(elChem, 'molecule');
      //       elMol.setAttribute(
      //           'molID', `m$ { ++watermark.molID }`);
      //       this.encodeMolecule(elMol, this.mol, watermark);
      //     }
      //     else
      //     {
      //       const elReaction = WebMolKit.XML.appendElement(elChem, 'reaction');

      //       const sectionMols : [ string, WebMolKit.Molecule[] ][] = [ [ 'reactantList', this.reactants ], [ 'agentList', this.agents ], [ 'productList', this.products ] ];
      //       for (const[section, molList] of sectionMols)
      //       {
      //         const elList = WebMolKit.XML.appendElement(elReaction, section);
      //         for (const mol of molList)
      //         {
      //           const elMol = WebMolKit.XML.appendElement(elList, 'molecule');
      //           elMol.setAttribute(
      //               'molID', `m$ { ++watermark.molID }`);
      //           this.encodeMolecule(elMol, mol, watermark);
      //         }
      //       }

      //       const elArrow = WebMolKit.XML.appendElement(elReaction, 'arrow');
      //       elArrow.setAttribute('type', 'DEFAULT');
      //       elArrow.setAttribute('x1', this.arrow.x1.toString());
      //       elArrow.setAttribute('y1', this.arrow.y1.toString());
      //       elArrow.setAttribute('x2', this.arrow.x2.toString());
      //       elArrow.setAttribute('y2', this.arrow.y2.toString());

      //       for (const plus of this.plusSigns)
      //       {
      //         const elSign = WebMolKit.XML.appendElement(elDoc, 'MReactionSign');
      //         elSign.setAttribute(
      //             'id', `o$ { ++watermark.objectID }`);
      //         for (const[k, v] of[[ 'toptions', 'NOROT' ], [ 'fontScale', '14.0' ], [ 'halign', 'CENTER' ], [ 'valign', 'CENTER' ], [ 'autoSize', 'true' ]])
      //         {
      //           elSign.setAttribute(k, v);
      //         }

      //         const elField = WebMolKit.XML.appendElement(elSign, 'Field');
      //         elField.setAttribute('name', 'text');
      //         WebMolKit.XML.appendText(elField, ' {D font=SansSerif,size=18,bold}+ ', true);
      //         for (const[dx, dy] of[[ -0.25, -0.25 ], [ 0.25, -0.25 ], [ 0.25, 0.25 ], [ -0.25, 0.25 ]])
      //         {
      //           const elPoint = WebMolKit.XML.appendElement(elSign, 'MPoint');
      //           elPoint.setAttribute('x', (plus.x + dx).toString());
      //           elPoint.setAttribute('y', (plus.y + dy).toString());
      //         }
      //       }
      //     }

      //     for (const box of this.textBoxes)
      //     {
      //       const elBox = WebMolKit.XML.appendElement(elDoc, 'MTextBox');
      //       elBox.setAttribute(
      //           'id', `o$ { ++watermark.objectID }`);
      //       elBox.setAttribute('toption', 'NOROT');
      //       elBox.setAttribute('fontScale', box.fontScale.toString());
      //       elBox.setAttribute('halign', 'LEFT');
      //       elBox.setAttribute('valign', 'TOP');
      //       elBox.setAttribute('autoSize', 'true');

      //       const elField = WebMolKit.XML.appendElement(elBox, 'Field');
      //       elField.setAttribute('name', 'text');
      //       WebMolKit.XML.setText(elField, box.text);
      //       for (let n = 0; n < 4; n++)
      //       {
      //         const elPoint = WebMolKit.XML.appendElement(elBox, 'MPoint');
      //         elPoint.setAttribute('x', box.x.toString());
      //         elPoint.setAttribute('y', box.y.toString());
      //       }
      //     }

      //     return WebMolKit.XML.toString(xml);
      //   }

      // private
      //   encodeMolecule(elMol
      //                  : Element, mol
      //                  : WebMolKit.Molecule, watermark
      //                  : Watermark) : string[]
      //   {
      //     const elAtomArray = WebMolKit.XML.appendElement(elMol, 'atomArray');

      //     const atomIDList : string[] = [], bondIDList : string[] = [];

      //     const SGIDPFX = 'xSGroupAttachmentPoint';
      //     const attachments : {sgid : string, atom : number}[] = [];

      //     const polyUnits = new WebMolKit.PolymerBlock(mol).getUnits();
      //     const polyRef : string[] = [];
      //     for (let n = 0; n < polyUnits.length; n++)
      //       polyRef.push(`sg$ { ++watermark.sgroupID }`);

      //     for (let n = 1; n <= mol.numAtoms; n++)
      //     {
      //       const id = `a$ { ++watermark.atomID }
      //       `;
      //       atomIDList.push(id);

      //       let label = mol.atomElement(n), sgid : string = null, alias : string = null;
      //       if (WebMolKit.MolUtil.hasAbbrev(mol, n) && mol.atomAdjCount(n) == 1)
      //       {
      //         label = 'R';
      //         sgid = `sg$ { ++watermark.sgroupID }
      //         `;
      //         attachments.push({sgid, atom : n});
      //       }
      //       else if (mol.atomicNumber(n) == 0 && label != '*' && label != 'R')
      //       {
      //         alias = label;
      //         label = 'C';
      //       }

      //       const elAtom = WebMolKit.XML.appendElement(elAtomArray, 'atom');
      //       elAtom.setAttribute('id', id);
      //       elAtom.setAttribute('elementType', label);
      //       elAtom.setAttribute('x2', mol.atomX(n).toString());
      //       elAtom.setAttribute('y2', mol.atomY(n).toString());
      //       if (mol.atomCharge(n) != 0)
      //         elAtom.setAttribute('formalCharge', mol.atomCharge(n).toString());
      //       if (mol.atomUnpaired(n) == 1)
      //       {
      //         elAtom.setAttribute('radical', 'monovalent');
      //       }
      //       else if (mol.atomUnpaired(n) == 2)
      //       {
      //         elAtom.setAttribute('radical', 'divalent1');
      //       }
      //       else if (mol.atomUnpaired(n) == 3)
      //       {
      //         elAtom.setAttribute('radical', 'trivalent4');
      //       }
      //       if (mol.atomIsotope(n) > 0)
      //       {
      //         elAtom.setAttribute('isotope', mol.atomIsotope(n).toString());
      //       }
      //       if (mol.atomMapNum(n) > 0)
      //       {
      //         elAtom.setAttribute('mrvMap', mol.atomMapNum(n).toString());
      //       }
      //       if (alias != null)
      //       {
      //         elAtom.setAttribute('mrvAlias', alias);
      //       }
      //       const refids : string[] = [];
      //       if (sgid)
      //         refids.push(sgid);
      //       for (let i = 0; i < polyUnits.length; i++)
      //       {
      //         if (polyUnits[i].atoms.includes(n))
      //           refids.push(polyRef[i]);
      //       }
      //       if (refids.length > 0)
      //       {
      //         elAtom.setAttribute('sgroupRef', refids.join(' '));
      //       }
      //       if (mol.atomExtra(n).includes(SGIDPFX))
      //       {
      //         elAtom.setAttribute('sgroupAttachmentPoint', '1');
      //       }
      //     }

      //     const elBondArray = WebMolKit.XML.appendElement(elMol, 'bondArray');

      //     for (let n = 1; n <= mol.numBonds; n++)
      //     {
      //       const id = `b$ { ++watermark.bondID }
      //       `;
      //       bondIDList.push(id);
      //       const bfr = mol.bondFrom(n), bto = mol.bondTo(n), order = mol.bondOrder(n), type = mol.bondType(n);

      //       const elBond = WebMolKit.XML.appendElement(elBondArray, 'bond');
      //       elBond.setAttribute('id', id);
      //       elBond.setAttribute(
      //           'atomRefs2', `${atomIDList[bfr - 1]} $ { atomIDList[bto - 1] }`);
      //       elBond.setAttribute('order', order.toString());

      //       if (type > 0)
      //       {
      //         const elStereo = WebMolKit.XML.appendElement(elBond, 'bondStereo');
      //         elStereo.setAttribute('convention', 'MDL');
      //         if (type == WebMolKit.Molecule.BONDTYPE_INCLINED)
      //         {
      //           elStereo.setAttribute('conventionValue', '1');
      //         }
      //         else if (type == WebMolKit.Molecule.BONDTYPE_DECLINED)
      //         {
      //           elStereo.setAttribute('conventionValue', '6');
      //         }
      //         else if (type == WebMolKit.Molecule.BONDTYPE_UNKNOWN)
      //         {
      //           if (order == 1)
      //           {
      //             elStereo.setAttribute('conventionValue', '4');
      //           }
      //           else if (order == 2)
      //           {
      //             elStereo.setAttribute('conventionValue', '3');
      //           }
      //         }
      //       }
      //     }

      //     for (const attachment of attachments)
      //     {
      //       const elSgroup = WebMolKit.XML.appendElement(elMol, 'molecule');
      //       elSgroup.setAttribute(
      //           'molID', `m$ { ++watermark.molID }`);
      //       elSgroup.setAttribute('id', attachment.sgid);
      //       elSgroup.setAttribute('role', 'SuperatomSgroup');
      //       elSgroup.setAttribute('title', mol.atomElement(attachment.atom));

      //       const abvmol = WebMolKit.MolUtil.getAbbrev(mol, attachment.atom);
      //       const sgmol = abvmol.clone();
      //       for (const a of sgmol.atomAdjList(1))
      //       {
      //         sgmol.setAtomExtra(a, [... sgmol.atomExtra(a), SGIDPFX ]);
      //       }
      //       sgmol.deleteAtomAndBonds(1);
      //       const sgAtomList = this.encodeMolecule(elSgroup, sgmol, watermark);

      //       const elAttArray = WebMolKit.XML.appendElement(elSgroup, 'AttachmentPointArray');
      //       for (const a of abvmol.atomAdjList(1))
      //       {
      //         const elAttPoint = WebMolKit.XML.appendElement(elAttArray, 'AttachmentPoint');
      //         elAttPoint.setAttribute('atom', sgAtomList[a - 2]);
      //         elAttPoint.setAttribute('order', abvmol.bondOrder(abvmol.findBond(1, a)).toString());
      //         elAttPoint.setAttribute('bond', bondIDList[mol.atomAdjBonds(attachment.atom)[0] - 1]);
      //       }
      //     }

      //     for (let i = 0; i < polyUnits.length; i++)
      //     {
      //       const atomRefs = polyUnits[i].atoms.map((a) = > atomIDList[a - 1]);

      //       const elSgroup = WebMolKit.XML.appendElement(elMol, 'molecule');
      //       elSgroup.setAttribute(
      //           'molID', `m$ { ++watermark.molID }`);
      //       elSgroup.setAttribute('id', polyRef[i]);
      //       elSgroup.setAttribute('role', 'SruSgroup');
      //       elSgroup.setAttribute('atomRefs', atomRefs.join(' '));
      //       elSgroup.setAttribute('title', 'n');
      //       if (polyUnits[i].connect == WebMolKit.PolymerBlockConnectivity.HeadToHead)
      //       {
      //         elSgroup.setAttribute('connect', 'hh');
      //       }
      //       else if (polyUnits[i].connect == WebMolKit.PolymerBlockConnectivity.Random)
      //       {
      //         elSgroup.setAttribute('connect', 'eu');
      //       }
      //       else
      //       {
      //         elSgroup.setAttribute('connect', 'ht');
      //       }
      //       elSgroup.setAttribute('correspondence', '');
      //       elSgroup.setAttribute('bondList', '');
      //     }
      //     return atomIDList;
      //   }
      //}
  };

  //------------------------------------------------
  //
  //  Read a molecule from a stream
  //
  //------------------------------------------------
  
  RWMol *MrvDataStreamToMol(std::istream *inStream)
  {
    PRECONDITION(inStream, "no stream");
    std::string tempStr;
    RWMol *res = NULL;
    MarvinCML marvinCML;

    marvinCML.parse(*inStream);

    printf("Reactant Count: %lu\n",marvinCML.reactants.size());
    printf("Agent Count: %lu\n",marvinCML.agents.size());
    printf("Products Count: %lu\n",marvinCML.products.size());


    res = marvinCML.mol;

    return res;
  }

  RWMol *MrvDataStreamToMol(std::istream &inStream)
  {
    return MrvDataStreamToMol(&inStream);
  }
  //------------------------------------------------
  //
  //  Read a molecule from a string
  //
  //------------------------------------------------
  RWMol *MrvBlockToMol(const std::string &molBlock)
  {
    std::istringstream inStream(molBlock);
    // unsigned int line = 0;
    return MrvDataStreamToMol(inStream);
  }

  //------------------------------------------------
  //
  //  Read a molecule from a file
  //
  //------------------------------------------------
  RWMol *MrvFileToMol(const std::string &fName)
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
    {
      res = MrvDataStreamToMol(inStream);
    }
    return res;
  }
}// namespace RDKit
