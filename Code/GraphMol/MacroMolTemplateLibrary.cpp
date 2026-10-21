//
//  Copyright (C) 2026 Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//
#include "MacroMolTemplateLibrary.h"

#include <array>
#include <cstddef>
#include <memory>
#include <mutex>
#include <GraphMol/Chirality.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <cstdlib>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace RDKit {
namespace {

// Generated from the validated built-in template graphs; no parser is used at runtime.
struct DirectAtomDef { unsigned atomicNum; int formalCharge; unsigned isotope; Atom::ChiralType chiralTag; unsigned explicitHs; bool noImplicit; const char *pdbName; };
struct DirectBondDef { unsigned beginAtomIdx; unsigned endAtomIdx; Bond::BondType type; bool aromatic; Bond::BondStereo stereo; std::vector<unsigned int> stereoAtoms; };
struct BuiltinTemplateDef { const char *symbol; const char *name; MonomerClass monomerClass; const char *originalData; std::vector<DirectAtomDef> atoms; std::vector<DirectBondDef> bonds; std::vector<unsigned int> mainAtomIdxs; std::vector<MacroMolLeavingGroup> leavingGroups; };


const std::array<BuiltinTemplateDef, 22> builtinTemplates{{
    {"A", "ALA", MonomerClass::AminoAcid, "C[C@H](N[H])C(=O)[OH] |atomProp:0.pdbName. CB :1.pdbName. CA :2.pdbName. N :3.pdbName. H :4.pdbName. C :5.pdbName. O :6.pdbName. OXT|",
     {
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 4, 5}, {
      {{3}, 2, 3, 1},
      {{6}, 4, 6, 2},
     }},
    {"R", "ARG", MonomerClass::AminoAcid, "N=C(N)NCCC[C@H](N[H])C(=O)[OH] |atomProp:0.pdbName. NH2:1.pdbName. CZ :2.pdbName. NH1:3.pdbName. NE :4.pdbName. CD :5.pdbName. CG :6.pdbName. CB :7.pdbName. CA :8.pdbName. N :9.pdbName. H :10.pdbName. C :11.pdbName. O :12.pdbName. OXT|",
     {
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "NH2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CZ "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "NH1"},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "NE "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {6, 7, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {7, 8, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {8, 9, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {7, 10, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {10, 11, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {10, 12, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11}, {
      {{9}, 8, 9, 1},
      {{12}, 10, 12, 2},
     }},
    {"N", "ASN", MonomerClass::AminoAcid, "NC(=O)C[C@H](N[H])C(=O)[OH] |atomProp:0.pdbName. ND2:1.pdbName. CG :2.pdbName. OD1:3.pdbName. CB :4.pdbName. CA :5.pdbName. N :6.pdbName. H :7.pdbName. C :8.pdbName. O :9.pdbName. OXT|",
     {
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "ND2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "OD1"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 7, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {7, 8, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {7, 9, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 7, 8}, {
      {{6}, 5, 6, 1},
      {{9}, 7, 9, 2},
     }},
    {"D", "ASP", MonomerClass::AminoAcid, "O=C([C@H](CC(=O)[OH])N[H])[OH] |atomProp:0.pdbName. O :1.pdbName. C :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. OD1:6.pdbName. OD2:7.pdbName. N :8.pdbName. H :9.pdbName. OXT|",
     {
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "OD1"},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OD2"},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 7, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {7, 8, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 9, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 7}, {
      {{8}, 7, 8, 1},
      {{9}, 1, 9, 2},
      {{6}, 4, 6, 3},
     }},
    {"C", "CYS", MonomerClass::AminoAcid, "O=C([C@H](CS[H])N[H])[OH] |atomProp:0.pdbName. O :1.pdbName. C :2.pdbName. CA :3.pdbName. CB :4.pdbName. SG :5.pdbName. HG :6.pdbName. N :7.pdbName. H :8.pdbName. OXT|",
     {
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {16, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "SG "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "HG "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {6, 7, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 8, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 6}, {
      {{7}, 6, 7, 1},
      {{8}, 1, 8, 2},
      {{5}, 4, 5, 3},
     }},
    {"Q", "GLN", MonomerClass::AminoAcid, "NC(=O)CC[C@H](N[H])C(=O)[OH] |atomProp:0.pdbName. NE2:1.pdbName. CD :2.pdbName. OE1:3.pdbName. CG :4.pdbName. CB :5.pdbName. CA :6.pdbName. N :7.pdbName. H :8.pdbName. C :9.pdbName. O :10.pdbName. OXT|",
     {
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "NE2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "OE1"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {6, 7, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 8, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {8, 9, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {8, 10, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 6, 8, 9}, {
      {{7}, 6, 7, 1},
      {{10}, 8, 10, 2},
     }},
    {"E", "GLU", MonomerClass::AminoAcid, "O=C([C@H](CCC(=O)[OH])N[H])[OH] |atomProp:0.pdbName. O :1.pdbName. C :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. CD :6.pdbName. OE1:7.pdbName. OE2:8.pdbName. N :9.pdbName. H :10.pdbName. OXT|",
     {
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "OE1"},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OE2"},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 7, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 8, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {8, 9, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 10, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 6, 8}, {
      {{9}, 8, 9, 1},
      {{10}, 1, 10, 2},
      {{7}, 5, 7, 3},
     }},
    {"G", "GLY", MonomerClass::AminoAcid, "O=C(CN[H])[OH] |atomProp:0.pdbName. O :1.pdbName. C :2.pdbName. CA :3.pdbName. N :4.pdbName. H :5.pdbName. OXT|",
     {
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CA "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3}, {
      {{4}, 3, 4, 1},
      {{5}, 1, 5, 2},
     }},
    {"H", "HIS", MonomerClass::AminoAcid, "O=C([C@H](Cc1cnc[nH]1)N[H])[OH] |atomProp:0.pdbName. O :1.pdbName. C :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. CD2:6.pdbName. NE2:7.pdbName. CE1:8.pdbName. ND1:9.pdbName. N :10.pdbName. H :11.pdbName. OXT|",
     {
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD2"},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "NE2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CE1"},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 1, false, "ND1"},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {6, 7, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {7, 8, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {2, 9, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {9, 10, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 11, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {8, 4, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, {
      {{10}, 9, 10, 1},
      {{11}, 1, 11, 2},
     }},
    {"I", "ILE", MonomerClass::AminoAcid, "CC[C@H](C)[C@H](N[H])C(=O)[OH] |atomProp:0.pdbName. CD1:1.pdbName. CG1:2.pdbName. CB :3.pdbName. CG2:4.pdbName. CA :5.pdbName. N :6.pdbName. H :7.pdbName. C :8.pdbName. O :9.pdbName. OXT|",
     {
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD1"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG1"},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 7, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {7, 8, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {7, 9, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 7, 8}, {
      {{6}, 5, 6, 1},
      {{9}, 7, 9, 2},
     }},
    {"L", "LEU", MonomerClass::AminoAcid, "CC(C)C[C@H](N[H])C(=O)[OH] |atomProp:0.pdbName. CD1:1.pdbName. CG :2.pdbName. CD2:3.pdbName. CB :4.pdbName. CA :5.pdbName. N :6.pdbName. H :7.pdbName. C :8.pdbName. O :9.pdbName. OXT|",
     {
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD1"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 7, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {7, 8, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {7, 9, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 7, 8}, {
      {{6}, 5, 6, 1},
      {{9}, 7, 9, 2},
     }},
    {"K", "LYS", MonomerClass::AminoAcid, "O=C([C@H](CCCCN[H])N[H])[OH] |atomProp:0.pdbName. O :1.pdbName. C :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. CD :6.pdbName. CE :7.pdbName. NZ :8.pdbName. HZ1:9.pdbName. N :10.pdbName. H :11.pdbName. OXT|",
     {
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CE "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "NZ "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "HZ1"},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {6, 7, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {7, 8, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 9, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {9, 10, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 11, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 6, 7, 9}, {
      {{10}, 9, 10, 1},
      {{11}, 1, 11, 2},
      {{8}, 7, 8, 3},
     }},
    {"M", "MET", MonomerClass::AminoAcid, "CSCC[C@H](N[H])C(=O)[OH] |atomProp:0.pdbName. CE :1.pdbName. SD :2.pdbName. CG :3.pdbName. CB :4.pdbName. CA :5.pdbName. N :6.pdbName. H :7.pdbName. C :8.pdbName. O :9.pdbName. OXT|",
     {
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CE "},
      {16, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "SD "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 7, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {7, 8, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {7, 9, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 7, 8}, {
      {{6}, 5, 6, 1},
      {{9}, 7, 9, 2},
     }},
    {"F", "PHE", MonomerClass::AminoAcid, "O=C([C@H](Cc1ccccc1)N[H])[OH] |atomProp:0.pdbName. O :1.pdbName. C :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. CD1:6.pdbName. CE1:7.pdbName. CZ :8.pdbName. CE2:9.pdbName. CD2:10.pdbName. N :11.pdbName. H :12.pdbName. OXT|",
     {
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD1"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CE1"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CZ "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CE2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD2"},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {6, 7, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {7, 8, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {8, 9, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {2, 10, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {10, 11, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 12, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {9, 4, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, {
      {{11}, 10, 11, 1},
      {{12}, 1, 12, 2},
     }},
    {"P", "PRO", MonomerClass::AminoAcid, "O=C([C@@H]1CCCN1[H])[OH] |atomProp:0.pdbName. O :1.pdbName. C :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. CD :6.pdbName. N :7.pdbName. H :8.pdbName. OXT|",
     {
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {6, 7, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 8, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {6, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 6}, {
      {{7}, 6, 7, 1},
      {{8}, 1, 8, 2},
     }},
    {"S", "SER", MonomerClass::AminoAcid, "O=C([C@H](CO)N[H])[OH] |atomProp:0.pdbName. O :1.pdbName. C :2.pdbName. CA :3.pdbName. CB :4.pdbName. OG :5.pdbName. N :6.pdbName. H :7.pdbName. OXT|",
     {
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "OG "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 7, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5}, {
      {{6}, 5, 6, 1},
      {{7}, 1, 7, 2},
     }},
    {"T", "THR", MonomerClass::AminoAcid, "C[C@@H](O)[C@H](N[H])C(=O)[OH] |atomProp:0.pdbName. CG2:1.pdbName. CB :2.pdbName. OG1:3.pdbName. CA :4.pdbName. N :5.pdbName. H :6.pdbName. C :7.pdbName. O :8.pdbName. OXT|",
     {
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(1), 1, true, "CB "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "OG1"},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {6, 7, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {6, 8, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 6, 7}, {
      {{5}, 4, 5, 1},
      {{8}, 6, 8, 2},
     }},
    {"W", "TRP", MonomerClass::AminoAcid, "O=C([C@H](Cc1c[nH]c2ccccc12)N[H])[OH] |atomProp:0.pdbName. O :1.pdbName. C :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. CD1:6.pdbName. NE1:7.pdbName. CE2:8.pdbName. CZ2:9.pdbName. CH2:10.pdbName. CZ3:11.pdbName. CE3:12.pdbName. CD2:13.pdbName. N :14.pdbName. H :15.pdbName. OXT|",
     {
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD1"},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 1, false, "NE1"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CE2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CZ2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CH2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CZ3"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CE3"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD2"},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {6, 7, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {7, 8, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {8, 9, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {9, 10, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {10, 11, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {11, 12, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {2, 13, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {13, 14, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 15, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {12, 4, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {12, 7, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13}, {
      {{14}, 13, 14, 1},
      {{15}, 1, 15, 2},
     }},
    {"Y", "TYR", MonomerClass::AminoAcid, "O=C([C@H](Cc1ccc(O)cc1)N[H])[OH] |atomProp:0.pdbName. O :1.pdbName. C :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. CD1:6.pdbName. CE1:7.pdbName. CZ :8.pdbName. OH :9.pdbName. CE2:10.pdbName. CD2:11.pdbName. N :12.pdbName. H :13.pdbName. OXT|",
     {
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD1"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CE1"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CZ "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "OH "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CE2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD2"},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {6, 7, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {7, 8, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {7, 9, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {9, 10, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {2, 11, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {11, 12, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 13, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {10, 4, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}, {
      {{12}, 11, 12, 1},
      {{13}, 1, 13, 2},
     }},
    {"V", "VAL", MonomerClass::AminoAcid, "CC(C)[C@H](N[H])C(=O)[OH] |atomProp:0.pdbName. CG1:1.pdbName. CB :2.pdbName. CG2:3.pdbName. CA :4.pdbName. N :5.pdbName. H :6.pdbName. C :7.pdbName. O :8.pdbName. OXT|",
     {
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG1"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {6, 7, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {6, 8, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 6, 7}, {
      {{5}, 4, 5, 1},
      {{8}, 6, 8, 2},
     }},
    {"U", "SEC", MonomerClass::AminoAcid, "O=C([C@H](C[SeH])N[H])[OH] |atomProp:0.pdbName. O :1.pdbName. C :2.pdbName. CA :3.pdbName. CB :4.pdbName. SE :5.pdbName. N :6.pdbName. H :7.pdbName. OXT|",
     {
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {34, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "SE "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 7, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5}, {
      {{6}, 5, 6, 1},
      {{7}, 1, 7, 2},
     }},
    {"O", "PYL", MonomerClass::AminoAcid, "C[C@@H]1CC=N[C@H]1C(=O)NCCCC[C@H](N[H])C(=O)[OH] |atomProp:0.pdbName. CB2:1.pdbName. CG2:2.pdbName. CD2:3.pdbName. CE2:4.pdbName. N2 :5.pdbName. CA2:6.pdbName. C2 :7.pdbName. O2 :8.pdbName. NZ :9.pdbName. CE :10.pdbName. CD :11.pdbName. CG :12.pdbName. CB :13.pdbName. CA :14.pdbName. N :15.pdbName. H :16.pdbName. C :17.pdbName. O :18.pdbName. OXT|",
     {
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CG2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CE2"},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N2 "},
      {6, 0, 0, static_cast<Atom::ChiralType>(1), 1, true, "CA2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C2 "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O2 "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "NZ "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CE "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {6, 7, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {6, 8, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {8, 9, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {9, 10, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {10, 11, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {11, 12, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {12, 13, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {13, 14, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {14, 15, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {13, 16, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {16, 17, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {16, 18, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 1, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17}, {
      {{15}, 14, 15, 1},
      {{18}, 16, 18, 2},
     }},
}};

std::unique_ptr<MacroMolTemplate> makeBuiltinTemplate(
    const BuiltinTemplateDef &def) {
  RWMol mol;
  for (const auto &atomDef : def.atoms) {
    auto *atom = new Atom(atomDef.atomicNum);
    atom->setFormalCharge(atomDef.formalCharge);
    atom->setIsotope(atomDef.isotope);
    atom->setChiralTag(atomDef.chiralTag);
    atom->setNumExplicitHs(atomDef.explicitHs);
    atom->setNoImplicit(atomDef.noImplicit);
    if (atomDef.pdbName) atom->setProp("pdbName", std::string(atomDef.pdbName));
    mol.addAtom(atom, true);
  }
  for (const auto &bondDef : def.bonds) {
    auto *bond = new Bond(bondDef.type);
    bond->setIsAromatic(bondDef.aromatic);
    if (bondDef.stereoAtoms.size() == 2) bond->setStereoAtoms(bondDef.stereoAtoms[0], bondDef.stereoAtoms[1]);
    bond->setStereo(bondDef.stereo);
    bond->setBeginAtomIdx(bondDef.beginAtomIdx);
    bond->setEndAtomIdx(bondDef.endAtomIdx);
    mol.addBond(bond, true);
  }
  MacroMolTemplateBuilder builder(mol, def.monomerClass, def.name, def.symbol, def.originalData);
  auto mainAtomIdxs = def.mainAtomIdxs;
  builder.setMainGroup(std::move(mainAtomIdxs));
  for (const auto &leavingGroup : def.leavingGroups) {
    auto leavingGroupCopy = leavingGroup;
    builder.addLeavingGroup(std::move(leavingGroupCopy));
  }
  return builder.build();
}
}
namespace detail {
// Defined in DefaultMonomers.cpp, generated from Data/DefaultMonomers.json
std::string_view getEmbeddedDefaultMonomerDefinitions();
}  // namespace detail

namespace {

constexpr const char *monomerLibraryEnvVar = "RDKIT_MONOMER_LIBRARY";

struct MonomerDefinitions {
  std::string json;
  std::string source;  // describes where the definitions came from
};

// Uses the JSON file named by $RDKIT_MONOMER_LIBRARY if it is set, otherwise
// the default definitions embedded at build time.
MonomerDefinitions readMonomerDefinitions() {
  const auto *fileName = std::getenv(monomerLibraryEnvVar);
  if (!fileName || !*fileName) {
    return {std::string(detail::getEmbeddedDefaultMonomerDefinitions()),
            "embedded default monomer definitions"};
  }
  std::ifstream input(fileName);
  if (!input) {
    throw std::runtime_error(std::string("Could not open monomer library '") +
                             fileName + "' specified by " +
                             monomerLibraryEnvVar);
  }
  return {std::string(std::istreambuf_iterator<char>(input),
                      std::istreambuf_iterator<char>()),
          std::string("monomer library '") + fileName + "' specified by " +
              monomerLibraryEnvVar};
}

std::string getString(const boost::property_tree::ptree &object,
                      const std::string &key) {
  // ptree::get<std::string>() would return an empty string for an object or
  // array, so check that the value has no children. An empty object can't be
  // distinguished from an empty string, so empty values are rejected too.
  const auto value = object.get_child_optional(key);
  if (!value || !value->empty() || value->data().empty()) {
    throw std::runtime_error("Missing or invalid '" + key +
                             "' in monomer definition");
  }
  return value->data();
}

MonomerClass monomerClassFromPolymerType(const std::string &polymerType) {
  if (polymerType == "PEPTIDE") {
    return MonomerClass::AminoAcid;
  }
  if (polymerType == "RNA") {
    return MonomerClass::NucleicAcid;
  }
  throw std::runtime_error("Unsupported polymer type in monomer definition: " +
                           polymerType);
}

std::unique_ptr<RWMol> parseCXSmiles(const std::string &smiles,
                                     const std::string &what,
                                     const std::string &name) {
  v2::SmilesParse::SmilesParserParams params;
  params.removeHs = false;
  params.parseName = false;
  auto mol = v2::SmilesParse::MolFromSmiles(smiles, params);
  if (!mol) {
    throw std::runtime_error("Could not parse " + what + " for monomer " +
                             name);
  }
  // E/Z labels depend on CIP ranks, which change when dummy atoms are
  // replaced, so express double bond stereo relative to the stereo atoms
  for (const auto bond : mol->bonds()) {
    if ((bond->getStereo() == Bond::STEREOE ||
         bond->getStereo() == Bond::STEREOZ) &&
        bond->getStereoAtoms().size() == 2) {
      bond->setStereo(Chirality::translateEZLabelToCisTrans(bond->getStereo()));
    }
  }
  return mol;
}

// Returns n if the atom is a dummy atom labelled _Rn, otherwise 0.
int getRGroupNumber(const Atom *atom) {
  std::string label;
  if (atom->getAtomicNum() != 0 ||
      !atom->getPropIfPresent(common_properties::atomLabel, label) ||
      label.size() < 3 || label.compare(0, 2, "_R") != 0 ||
      label.find_first_not_of("0123456789", 2) != std::string::npos ||
      label[2] == '0') {
    return 0;
  }
  return std::stoi(label.substr(2));
}

// Replaces the core dummy atom at coreDummyIdx with the R-group, and returns
// the resulting leaving group. The R-group atom bonded to the R-group's dummy
// atom takes the place of the core dummy atom, keeping its index and bonds so
// that the stereochemistry of the attachment atom is preserved. Any further
// R-group atoms are appended to the molecule.
MacroMolLeavingGroup attachRGroup(RWMol &mol, unsigned int coreDummyIdx,
                                  int rGroupNumber, const ROMol &rGroup,
                                  const std::string &name) {
  const auto rGroupLabel = "R" + std::to_string(rGroupNumber);
  const Atom *rGroupDummy = nullptr;
  for (const auto atom : rGroup.atoms()) {
    if (atom->getAtomicNum() == 0) {
      if (rGroupDummy || getRGroupNumber(atom) != rGroupNumber) {
        throw std::runtime_error(
            "R-group " + rGroupLabel + " for monomer " + name +
            " must contain exactly one dummy atom, labelled _" + rGroupLabel);
      }
      rGroupDummy = atom;
    }
  }
  if (!rGroupDummy || rGroupDummy->getDegree() != 1) {
    throw std::runtime_error("The dummy atom in R-group " + rGroupLabel +
                             " for monomer " + name +
                             " must have exactly one neighbor");
  }
  const auto rGroupLeavingAtom = *rGroup.atomNeighbors(rGroupDummy).begin();

  const auto coreDummy = mol.getAtomWithIdx(coreDummyIdx);
  if (coreDummy->getDegree() != 1) {
    throw std::runtime_error("Attachment point _" + rGroupLabel +
                             " in monomer " + name +
                             " must have exactly one neighbor");
  }
  const auto attachAtomIdx = (*mol.atomNeighbors(coreDummy).begin())->getIdx();

  std::vector<int> newIdxs(rGroup.getNumAtoms(), -1);
  newIdxs[rGroupLeavingAtom->getIdx()] = coreDummyIdx;
  mol.replaceAtom(coreDummyIdx, rGroupLeavingAtom);
  std::vector<unsigned int> atomIdxs{coreDummyIdx};
  for (const auto atom : rGroup.atoms()) {
    if (atom != rGroupDummy && atom != rGroupLeavingAtom) {
      newIdxs[atom->getIdx()] = mol.addAtom(atom->copy(), false, true);
      atomIdxs.push_back(newIdxs[atom->getIdx()]);
    }
  }
  // the R-group's dummy atom corresponds to the core attachment atom
  const auto mapIdx = [&](unsigned int rGroupIdx) {
    return rGroupIdx == rGroupDummy->getIdx() ? attachAtomIdx
                                              : newIdxs[rGroupIdx];
  };
  for (const auto bond : rGroup.bonds()) {
    if (bond->getBeginAtom() == rGroupDummy ||
        bond->getEndAtom() == rGroupDummy) {
      continue;  // replaced by the core bond to the attachment atom
    }
    auto newBond = bond->copy();
    newBond->setBeginAtomIdx(mapIdx(bond->getBeginAtomIdx()));
    newBond->setEndAtomIdx(mapIdx(bond->getEndAtomIdx()));
    for (auto &stereoAtomIdx : newBond->getStereoAtoms()) {
      stereoAtomIdx = mapIdx(stereoAtomIdx);
    }
    mol.addBond(newBond, true);
  }

  // The leaving atom's chiral tag refers to the order of its bonds in the
  // R-group, which differs from the order in the combined molecule (where the
  // bond to the attachment atom comes first), so adjust it for the permutation.
  const auto leavingAtom = mol.getAtomWithIdx(coreDummyIdx);
  const auto chiralTag = leavingAtom->getChiralTag();
  if (chiralTag != Atom::CHI_UNSPECIFIED) {
    if (chiralTag != Atom::CHI_TETRAHEDRAL_CW &&
        chiralTag != Atom::CHI_TETRAHEDRAL_CCW) {
      throw std::runtime_error(
          "Unsupported chirality on the atom bonded to "
          "the attachment point of R-group " +
          rGroupLabel + " for monomer " + name);
    }
    INT_LIST bondOrder;
    for (const auto bond : rGroup.atomBonds(rGroupLeavingAtom)) {
      const auto otherIdx =
          mapIdx(bond->getOtherAtomIdx(rGroupLeavingAtom->getIdx()));
      bondOrder.push_back(
          mol.getBondBetweenAtoms(coreDummyIdx, otherIdx)->getIdx());
    }
    if (leavingAtom->getPerturbationOrder(bondOrder) % 2) {
      leavingAtom->invertChirality();
    }
  }
  return {std::move(atomIdxs), attachAtomIdx, coreDummyIdx, rGroupNumber};
}

// Builds a template from CORE_SMILES, a CXSMILES in which the attachment points
// are dummy atoms labelled _R1, _R2, etc., and RGROUPS, an object mapping the
// labels R1, R2, etc. to CXSMILES for the leaving groups. Each leaving group
// contains a dummy atom with the same label marking its attachment point.
std::unique_ptr<MacroMolTemplate> makeTemplate(
    const boost::property_tree::ptree &definition) {
  const auto symbol = getString(definition, "SYMBOL");
  const auto name = getString(definition, "PDBCODE");
  const auto coreSmiles = getString(definition, "CORE_SMILES");
  const auto monomerClass =
      monomerClassFromPolymerType(getString(definition, "POLYMER_TYPE"));
  const auto subclass = getString(definition, "MONOMER_TYPE");
  const auto rGroupDefinitions = definition.get_child_optional("RGROUPS");
  if (!rGroupDefinitions) {
    throw std::runtime_error("Missing 'RGROUPS' in definition of monomer " +
                             name);
  }

  auto mol = parseCXSmiles(coreSmiles, "CORE_SMILES", name);
  std::vector<std::pair<unsigned int, int>> coreDummies;
  for (const auto atom : mol->atoms()) {
    if (const auto rGroupNumber = getRGroupNumber(atom)) {
      coreDummies.emplace_back(atom->getIdx(), rGroupNumber);
    }
  }
  if (coreDummies.size() != rGroupDefinitions->size()) {
    throw std::runtime_error(
        "The attachment points in CORE_SMILES of monomer " + name +
        " do not match its RGROUPS");
  }

  std::vector<MacroMolLeavingGroup> leavingGroups;
  for (const auto &[coreDummyIdx, rGroupNumber] : coreDummies) {
    const auto rGroupLabel = "R" + std::to_string(rGroupNumber);
    const auto rGroupSmiles =
        rGroupDefinitions->get_child_optional(rGroupLabel);
    if (!rGroupSmiles || !rGroupSmiles->empty()) {
      throw std::runtime_error("Missing or invalid R-group " + rGroupLabel +
                               " in RGROUPS of monomer " + name);
    }
    const auto rGroup =
        parseCXSmiles(rGroupSmiles->data(), "R-group " + rGroupLabel, name);
    leavingGroups.push_back(
        attachRGroup(*mol, coreDummyIdx, rGroupNumber, *rGroup, name));
  }
  // Double bond stereo is fully described by the bond stereo and stereo atoms.
  // The bond directions from the SMILES are incomplete where dummy atoms were
  // replaced, so regenerate them before perceiving the stereochemistry (which
  // uses them with the legacy stereo perception), as parsing SMILES would.
  for (const auto bond : mol->bonds()) {
    if (bond->getBondDir() == Bond::ENDUPRIGHT ||
        bond->getBondDir() == Bond::ENDDOWNRIGHT) {
      bond->setBondDir(Bond::NONE);
    }
  }
  MolOps::sanitizeMol(*mol);
  MolOps::setDoubleBondNeighborDirections(*mol);
  const bool cleanIt = true;
  const bool force = true;
  MolOps::assignStereochemistry(*mol, cleanIt, force);

  std::vector<bool> isLeavingAtom(mol->getNumAtoms(), false);
  for (const auto &leavingGroup : leavingGroups) {
    for (const auto atomIdx : leavingGroup.atomIdxs) {
      isLeavingAtom[atomIdx] = true;
    }
  }
  std::vector<unsigned int> mainAtomIdxs;
  for (unsigned int atomIdx = 0; atomIdx < mol->getNumAtoms(); ++atomIdx) {
    if (!isLeavingAtom[atomIdx]) {
      mainAtomIdxs.push_back(atomIdx);
    }
  }

  // The original data is the monomer's JSON definition. Note that the boost
  // JSON writer escapes '/' as \/, which is valid JSON.
  std::ostringstream originalData;
  boost::property_tree::write_json(originalData, definition, false);
  auto originalDataStr = originalData.str();
  if (!originalDataStr.empty() && originalDataStr.back() == '\n') {
    originalDataStr.pop_back();
  }

  MacroMolTemplateBuilder builder(*mol, monomerClass, name, symbol,
                                  std::move(originalDataStr));
  builder.setMainGroup(std::move(mainAtomIdxs));
  builder.setSubclass(subclass);
  for (auto &leavingGroup : leavingGroups) {
    builder.addLeavingGroup(std::move(leavingGroup));
  }
  return builder.build();
}

const std::array<std::unique_ptr<const MacroMolTemplate>, builtinTemplates.size()> &getBuiltinTemplates() {
  static const auto templates = [] {
    std::array<std::unique_ptr<const MacroMolTemplate>, builtinTemplates.size()> result{};
    for (std::size_t i = 0; i < builtinTemplates.size(); ++i) result[i] = makeBuiltinTemplate(builtinTemplates[i]);
    return result;
  }();
  return templates;
}
// Builds a library from the monomer definitions. Any error, including
// duplicate monomers, is reported along with the source of the definitions.
MacroMolTemplateLibrary loadMacroMolTemplates() {
  const auto definitions = readMonomerDefinitions();
  try {
    std::istringstream input(definitions.json);
    boost::property_tree::ptree tree;
    boost::property_tree::read_json(input, tree);
    MacroMolTemplateLibrary result;
    for (const auto &definition : tree) {
      result.addTemplate(makeTemplate(definition.second));
    }
    return result;
  } catch (const std::exception &e) {
    throw std::runtime_error("Error loading " + definitions.source + ": " +
                             e.what());
  }
}

}  // namespace

void addBuiltinMacroMolTemplates(MacroMolTemplateLibrary &templates) {
  for (const auto *builtinTemplate :
       getGlobalMacroMolTemplateLibrary().getTemplates()) {
    templates.addTemplate(std::make_unique<MacroMolTemplate>(*builtinTemplate));
  }
}

const MacroMolTemplateLibrary &getGlobalMacroMolTemplateLibrary() {
  // If loading throws, the static is left uninitialized and the next call
  // retries from scratch, so a partially populated library is never exposed.
  static const MacroMolTemplateLibrary templates = loadMacroMolTemplates();
  return templates;
}

}  // namespace RDKit
