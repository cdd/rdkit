//
//  Copyright (C) 2021 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <string>
#include <cstring>
#include <iostream>
#include <regex>

#include <RDGeneral/versions.h>
#include <atomic>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesJSONParsers.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/MolInterchange/MolInterchange.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Descriptors/Property.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Fingerprints/AtomPairs.h>
#include <GraphMol/Fingerprints/MACCS.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/CIPLabeler/CIPLabeler.h>
#include <GraphMol/Abbreviations/Abbreviations.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionPickler.h>
#include <GraphMol/Chirality.h>
#include <DataStructs/BitOps.h>

#include "common.h"

#include <sstream>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <RDGeneral/BoostEndInclude.h>

#ifdef RDK_BUILD_INCHI_SUPPORT
#include <INCHI-API/inchi.h>
#endif

#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>
#include "cffiwrapper.h"

namespace rj = rapidjson;

using namespace RDKit;

#if (defined(__GNUC__) || defined(__GNUG__))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion-null"
#endif

namespace {
char *str_to_c(const std::string &str, size_t *len = nullptr) {
  if (len) {
    *len = 0;
  }
  char *res;
  res = (char *)malloc(str.size() + 1);
  if (res) {
    if (len) {
      *len = str.size();
    }
    memcpy(res, str.c_str(), str.size());
    res[str.size()] = '\0';
  }
  return res;
}
char *str_to_c(const char *str) {
  char *res;
  res = (char *)malloc(strlen(str) + 1);
  if (res) {
    strcpy(res, str);
  }
  return res;
}
}  // namespace

void mol_to_pkl(
    const ROMol &mol, char **mol_pkl, size_t *mol_pkl_sz,
    unsigned int propFlags = PicklerOps::PropertyPickleOptions::AllProps ^
                             PicklerOps::PropertyPickleOptions::ComputedProps) {
  std::string pkl;
  MolPickler::pickleMol(mol, pkl, propFlags);
  free(*mol_pkl);
  *mol_pkl = str_to_c(pkl, mol_pkl_sz);
}

RWMol mol_from_pkl(const char *pkl, size_t pkl_sz) {
  if (!pkl || !pkl_sz) {
    return RWMol();
  }
  std::string mol_pkl(pkl, pkl_sz);
  RWMol res(mol_pkl);
  res.setProp(common_properties::_StereochemDone, 1, true);
  return res;
}

ChemicalReaction rxn_from_pkl(const char *pkl, size_t pkl_sz) {
  if (!pkl || !pkl_sz) {
    return ChemicalReaction();
  }
  std::string rxn_pkl(pkl, pkl_sz);
  ChemicalReaction res(rxn_pkl);
  return res;
}

#ifdef PT_OPT_GET
#undef PT_OPT_GET
#endif
#define PT_OPT_GET(opt) opt = pt.get(#opt, opt);

namespace {
SmilesWriteParams smiles_helper(const char *details_json) {
  SmilesWriteParams params;
  updateSmilesWriteParamsFromJSON(params, details_json);
  return params;
}
std::string cxsmiles_helper(const char *pkl, size_t pkl_sz,
                            const char *details_json) {
  if (!pkl || !pkl_sz) {
    return "";
  }
  auto params = smiles_helper(details_json);
  auto mol = mol_from_pkl(pkl, pkl_sz);
  std::uint32_t cxSmilesFields = SmilesWrite::CXSmilesFields::CX_ALL;
  unsigned int restoreBondDirs = RestoreBondDirOptionClear;
  updateCXSmilesFieldsFromJSON(cxSmilesFields, restoreBondDirs, details_json);
  return MolToCXSmiles(mol, params, cxSmilesFields,
                       static_cast<RestoreBondDirOption>(restoreBondDirs));
}
}  // namespace
extern "C" char *get_smiles(const char *pkl, size_t pkl_sz,
                            const char *details_json) {
  if (!pkl || !pkl_sz) {
    return nullptr;
  }
  auto params = smiles_helper(details_json);
  auto mol = mol_from_pkl(pkl, pkl_sz);
  auto data = MolToSmiles(mol, params);
  return str_to_c(data);
}
extern "C" char *get_smarts(const char *pkl, size_t pkl_sz,
                            const char *details_json) {
  if (!pkl || !pkl_sz) {
    return nullptr;
  }
  auto params = smiles_helper(details_json);
  auto mol = mol_from_pkl(pkl, pkl_sz);
  auto data = MolToSmarts(mol, params.doIsomericSmiles, params.rootedAtAtom);
  return str_to_c(data);
}
extern "C" char *get_cxsmiles(const char *pkl, size_t pkl_sz,
                              const char *details_json) {
  if (!pkl || !pkl_sz) {
    return nullptr;
  }
  auto data = cxsmiles_helper(pkl, pkl_sz, details_json);
  return str_to_c(data);
}
extern "C" char *get_cxsmarts(const char *pkl, size_t pkl_sz,
                              const char *details_json) {
  if (!pkl || !pkl_sz) {
    return nullptr;
  }
  auto params = smiles_helper(details_json);
  auto mol = mol_from_pkl(pkl, pkl_sz);
  auto data = MolToCXSmarts(mol, params.doIsomericSmiles);
  return str_to_c(data);
}
extern "C" char *get_molblock(const char *pkl, size_t pkl_sz,
                              const char *details_json) {
  if (!pkl || !pkl_sz) {
    return nullptr;
  }
  auto mol = mol_from_pkl(pkl, pkl_sz);
  auto data = MinimalLib::molblock_helper(mol, details_json, false);
  return str_to_c(data);
}
extern "C" char *get_v3kmolblock(const char *pkl, size_t pkl_sz,
                                 const char *details_json) {
  if (!pkl || !pkl_sz) {
    return nullptr;
  }
  auto mol = mol_from_pkl(pkl, pkl_sz);
  auto data = MinimalLib::molblock_helper(mol, details_json, true);
  return str_to_c(data);
}
extern "C" char *get_json(const char *pkl, size_t pkl_sz, const char *) {
  if (!pkl || !pkl_sz) {
    return nullptr;
  }
  auto mol = mol_from_pkl(pkl, pkl_sz);
  auto data = MolInterchange::MolToJSONData(mol);
  return str_to_c(data);
}
extern "C" void free_ptr(char *ptr) {
  if (ptr) {
    free(ptr);
  }
}

extern "C" char *get_svg(const char *pkl, size_t pkl_sz,
                         const char *details_json) {
  if (!pkl || !pkl_sz) {
    return nullptr;
  }
  auto mol = mol_from_pkl(pkl, pkl_sz);
  unsigned int width = MinimalLib::d_defaultWidth;
  unsigned int height = MinimalLib::d_defaultHeight;
  return str_to_c(MinimalLib::mol_to_svg(mol, width, height, details_json));
}

extern "C" char *get_rxn_svg(const char *pkl, size_t pkl_sz,
                             const char *details_json) {
  if (!pkl || !pkl_sz) {
    return nullptr;
  }
  auto rxn = rxn_from_pkl(pkl, pkl_sz);
  unsigned int width = MinimalLib::d_defaultWidth;
  unsigned int height = MinimalLib::d_defaultHeight;
  return str_to_c(MinimalLib::rxn_to_svg(rxn, width, height, details_json));
}

#ifdef RDK_BUILD_INCHI_SUPPORT
extern "C" char *get_inchi(const char *pkl, size_t pkl_sz,
                           const char *details_json) {
  if (!pkl || !pkl_sz) {
    return nullptr;
  }
  auto mol = mol_from_pkl(pkl, pkl_sz);
  ExtraInchiReturnValues rv;
  auto options = MinimalLib::parse_inchi_options(details_json);
  return str_to_c(
      MolToInchi(mol, rv, !options.empty() ? options.c_str() : nullptr));
}

extern "C" char *get_inchi_for_molblock(const char *ctab,
                                        const char *details_json) {
  if (!ctab) {
    return str_to_c("");
  }
  ExtraInchiReturnValues rv;
  auto options = MinimalLib::parse_inchi_options(details_json);
  return str_to_c(
      MolBlockToInchi(ctab, rv, !options.empty() ? options.c_str() : nullptr));
}

extern "C" char *get_inchikey_for_inchi(const char *inchi) {
  if (!inchi) {
    return str_to_c("");
  }
  return str_to_c(InchiToInchiKey(inchi));
}
#endif

extern "C" char *get_mol(const char *input, size_t *pkl_sz,
                         const char *details_json) {
  std::unique_ptr<RWMol> mol{MinimalLib::mol_from_input(input, details_json)};
  if (!mol) {
    *pkl_sz = 0;
    return nullptr;
  }
  unsigned int propFlags = PicklerOps::PropertyPickleOptions::AllProps ^
                           PicklerOps::PropertyPickleOptions::ComputedProps;
  MinimalLib::updatePropertyPickleOptionsFromJSON(propFlags, details_json);
  std::string pkl;
  MolPickler::pickleMol(*mol, pkl, propFlags);
  return str_to_c(pkl, pkl_sz);
}

extern "C" char *get_qmol(const char *input, size_t *pkl_sz,
                          const char *details_json) {
  std::unique_ptr<RWMol> mol{MinimalLib::qmol_from_input(input, details_json)};
  if (!mol) {
    *pkl_sz = 0;
    return nullptr;
  }
  static const unsigned int propFlags =
      PicklerOps::PropertyPickleOptions::AllProps ^
      PicklerOps::PropertyPickleOptions::ComputedProps;
  std::string pkl;
  MolPickler::pickleMol(*mol, pkl, propFlags);
  return str_to_c(pkl, pkl_sz);
}

extern "C" char *get_rxn(const char *input, size_t *pkl_sz,
                         const char *details_json) {
  std::unique_ptr<ChemicalReaction> rxn{
      MinimalLib::rxn_from_input(input, details_json)};
  if (!rxn) {
    *pkl_sz = 0;
    return nullptr;
  }
  unsigned int propFlags = PicklerOps::PropertyPickleOptions::AllProps ^
                           PicklerOps::PropertyPickleOptions::ComputedProps;
  std::string pkl;
  ReactionPickler::pickleReaction(*rxn, pkl, propFlags);
  return str_to_c(pkl, pkl_sz);
}

extern "C" char **get_mol_frags(const char *pkl, size_t pkl_sz,
                                size_t **frags_pkl_sz_array, size_t *num_frags,
                                const char *details_json,
                                char **mappings_json) {
  if (!pkl || !pkl_sz || !frags_pkl_sz_array || !num_frags) {
    return nullptr;
  }
  *frags_pkl_sz_array = nullptr;
  *num_frags = 0;
  auto mol = mol_from_pkl(pkl, pkl_sz);
  std::vector<int> frags;
  std::vector<std::vector<int>> fragsMolAtomMapping;
  bool sanitizeFrags = true;
  bool copyConformers = true;
  if (details_json) {
    std::string json = details_json;
    MinimalLib::get_mol_frags_details(json, sanitizeFrags, copyConformers);
  }
  std::vector<ROMOL_SPTR> molFrags;
  try {
    molFrags = MolOps::getMolFrags(mol, sanitizeFrags, &frags,
                                   &fragsMolAtomMapping, copyConformers);
  } catch (...) {
  }
  if (molFrags.empty()) {
    return nullptr;
  }
  char **molPklArray = (char **)malloc(sizeof(char *) * molFrags.size());
  if (!molPklArray) {
    return nullptr;
  }
  *frags_pkl_sz_array = (size_t *)malloc(sizeof(size_t) * molFrags.size());
  if (!*frags_pkl_sz_array) {
    free(molPklArray);
    return nullptr;
  }
  memset(molPklArray, 0, sizeof(char *) * molFrags.size());
  *num_frags = molFrags.size();
  for (size_t i = 0; i < molFrags.size(); ++i) {
    mol_to_pkl(*molFrags[i], &molPklArray[i], &(*frags_pkl_sz_array)[i]);
  }
  if (mappings_json) {
    auto res = MinimalLib::get_mol_frags_mappings(frags, fragsMolAtomMapping);
    *mappings_json = str_to_c(res);
  }
  return molPklArray;
}

extern "C" char *version() { return str_to_c(rdkitVersion); }

extern "C" char *get_substruct_match(const char *mol_pkl, size_t mol_pkl_sz,
                                     const char *query_pkl, size_t query_pkl_sz,
                                     const char *options_json) {
  if (!mol_pkl || !mol_pkl_sz || !query_pkl || !query_pkl_sz) {
    return nullptr;
  }
  auto mol = mol_from_pkl(mol_pkl, mol_pkl_sz);
  auto query = mol_from_pkl(query_pkl, query_pkl_sz);

  SubstructMatchParameters params;
  if (options_json) {
    std::string json(options_json);
    updateSubstructMatchParamsFromJSON(params, json);
  }
  params.maxMatches = 1;

  std::string res = "{}";
  auto matches = SubstructMatch(mol, query, params);
  if (!matches.empty()) {
    auto match = matches[0];
    rj::Document doc;
    doc.SetObject();
    MinimalLib::get_sss_json(mol, query, match, doc, doc);
    rj::StringBuffer buffer;
    rj::Writer<rj::StringBuffer> writer(buffer);
    doc.Accept(writer);
    res = buffer.GetString();
  }

  return str_to_c(res);
}
extern "C" char *get_substruct_matches(const char *mol_pkl, size_t mol_pkl_sz,
                                       const char *query_pkl,
                                       size_t query_pkl_sz,
                                       const char *options_json) {
  if (!mol_pkl || !mol_pkl_sz || !query_pkl || !query_pkl_sz) {
    return nullptr;
  }
  auto mol = mol_from_pkl(mol_pkl, mol_pkl_sz);
  auto query = mol_from_pkl(query_pkl, query_pkl_sz);

  SubstructMatchParameters params;
  if (options_json) {
    std::string json(options_json);
    updateSubstructMatchParamsFromJSON(params, json);
  }

  std::string res = "{}";
  auto matches = SubstructMatch(mol, query, params);
  if (!matches.empty()) {
    rj::Document doc;
    doc.SetArray();

    for (const auto &match : matches) {
      rj::Value rjMatch(rj::kObjectType);
      MinimalLib::get_sss_json(mol, query, match, rjMatch, doc);
      doc.PushBack(rjMatch, doc.GetAllocator());
    }

    rj::StringBuffer buffer;
    rj::Writer<rj::StringBuffer> writer(buffer);
    doc.Accept(writer);
    res = buffer.GetString();
  }

  return str_to_c(res);
}

extern "C" char *get_descriptors(const char *mol_pkl, size_t mol_pkl_sz) {
  if (!mol_pkl || !mol_pkl_sz) {
    return nullptr;
  }
  auto mol = mol_from_pkl(mol_pkl, mol_pkl_sz);
  return str_to_c(MinimalLib::get_descriptors(mol));
}

extern "C" char *get_morgan_fp(const char *mol_pkl, size_t mol_pkl_sz,
                               const char *details_json) {
  if (!mol_pkl || !mol_pkl_sz) {
    return nullptr;
  }
  try {
    auto fp = MinimalLib::morgan_fp_as_bitvect(
        mol_from_pkl(mol_pkl, mol_pkl_sz), details_json);
    auto res = BitVectToText(*fp);
    return str_to_c(res);
  } catch (...) {
    return nullptr;
  }
}

extern "C" char *get_morgan_fp_as_bytes(const char *mol_pkl, size_t mol_pkl_sz,
                                        size_t *nbytes,
                                        const char *details_json) {
  if (!mol_pkl || !mol_pkl_sz) {
    return nullptr;
  }
  try {
    auto fp = MinimalLib::morgan_fp_as_bitvect(
        mol_from_pkl(mol_pkl, mol_pkl_sz), details_json);
    auto res = BitVectToBinaryText(*fp);
    return str_to_c(res, nbytes);
  } catch (...) {
    return nullptr;
  }
}

extern "C" char *get_rdkit_fp(const char *mol_pkl, size_t mol_pkl_sz,
                              const char *details_json) {
  if (!mol_pkl || !mol_pkl_sz) {
    return nullptr;
  }
  try {
    auto fp = MinimalLib::rdkit_fp_as_bitvect(mol_from_pkl(mol_pkl, mol_pkl_sz),
                                              details_json);
    auto res = BitVectToText(*fp);
    return str_to_c(res);
  } catch (...) {
    return nullptr;
  }
}

extern "C" char *get_rdkit_fp_as_bytes(const char *mol_pkl, size_t mol_pkl_sz,
                                       size_t *nbytes,
                                       const char *details_json) {
  if (!mol_pkl || !mol_pkl_sz) {
    return nullptr;
  }
  try {
    auto fp = MinimalLib::rdkit_fp_as_bitvect(mol_from_pkl(mol_pkl, mol_pkl_sz),
                                              details_json);
    auto res = BitVectToBinaryText(*fp);
    return str_to_c(res, nbytes);
  } catch (...) {
    return nullptr;
  }
}

extern "C" char *get_pattern_fp(const char *mol_pkl, size_t mol_pkl_sz,
                                const char *details_json) {
  if (!mol_pkl || !mol_pkl_sz) {
    return nullptr;
  }
  try {
    auto fp = MinimalLib::pattern_fp_as_bitvect(
        mol_from_pkl(mol_pkl, mol_pkl_sz), details_json);
    auto res = BitVectToText(*fp);
    return str_to_c(res);
  } catch (...) {
    return nullptr;
  }
}

extern "C" char *get_pattern_fp_as_bytes(const char *mol_pkl, size_t mol_pkl_sz,
                                         size_t *nbytes,
                                         const char *details_json) {
  if (!mol_pkl || !mol_pkl_sz) {
    return nullptr;
  }
  try {
    auto fp = MinimalLib::pattern_fp_as_bitvect(
        mol_from_pkl(mol_pkl, mol_pkl_sz), details_json);
    auto res = BitVectToBinaryText(*fp);
    return str_to_c(res, nbytes);
  } catch (...) {
    return nullptr;
  }
}

extern "C" char *get_topological_torsion_fp(const char *mol_pkl,
                                            size_t mol_pkl_sz,
                                            const char *details_json) {
  if (!mol_pkl || !mol_pkl_sz) {
    return nullptr;
  }
  try {
    auto fp = MinimalLib::topological_torsion_fp_as_bitvect(
        mol_from_pkl(mol_pkl, mol_pkl_sz), details_json);
    auto res = BitVectToText(*fp);
    return str_to_c(res);
  } catch (...) {
    return nullptr;
  }
}

extern "C" char *get_topological_torsion_fp_as_bytes(const char *mol_pkl,
                                                     size_t mol_pkl_sz,
                                                     size_t *nbytes,
                                                     const char *details_json) {
  if (!mol_pkl || !mol_pkl_sz) {
    return nullptr;
  }
  try {
    auto fp = MinimalLib::topological_torsion_fp_as_bitvect(
        mol_from_pkl(mol_pkl, mol_pkl_sz), details_json);
    auto res = BitVectToBinaryText(*fp);
    return str_to_c(res, nbytes);
  } catch (...) {
    return nullptr;
  }
}

extern "C" char *get_atom_pair_fp(const char *mol_pkl, size_t mol_pkl_sz,
                                  const char *details_json) {
  if (!mol_pkl || !mol_pkl_sz) {
    return nullptr;
  }
  try {
    auto fp = MinimalLib::atom_pair_fp_as_bitvect(
        mol_from_pkl(mol_pkl, mol_pkl_sz), details_json);
    auto res = BitVectToText(*fp);
    return str_to_c(res);
  } catch (...) {
    return nullptr;
  }
}

extern "C" char *get_atom_pair_fp_as_bytes(const char *mol_pkl,
                                           size_t mol_pkl_sz, size_t *nbytes,
                                           const char *details_json) {
  if (!mol_pkl || !mol_pkl_sz) {
    return nullptr;
  }
  try {
    auto fp = MinimalLib::atom_pair_fp_as_bitvect(
        mol_from_pkl(mol_pkl, mol_pkl_sz), details_json);
    auto res = BitVectToBinaryText(*fp);
    return str_to_c(res, nbytes);
  } catch (...) {
    return nullptr;
  }
}

extern "C" char *get_maccs_fp(const char *mol_pkl, size_t mol_pkl_sz) {
  if (!mol_pkl || !mol_pkl_sz) {
    return nullptr;
  }
  try {
    auto fp =
        MinimalLib::maccs_fp_as_bitvect(mol_from_pkl(mol_pkl, mol_pkl_sz));
    auto res = BitVectToText(*fp);
    return str_to_c(res);
  } catch (...) {
    return nullptr;
  }
}

extern "C" char *get_maccs_fp_as_bytes(const char *mol_pkl, size_t mol_pkl_sz,
                                       size_t *nbytes) {
  if (!mol_pkl || !mol_pkl_sz) {
    return nullptr;
  }
  try {
    auto fp =
        MinimalLib::maccs_fp_as_bitvect(mol_from_pkl(mol_pkl, mol_pkl_sz));
    auto res = BitVectToBinaryText(*fp);
    return str_to_c(res, nbytes);
  } catch (...) {
    return nullptr;
  }
}

#ifdef RDK_BUILD_AVALON_SUPPORT
extern "C" char *get_avalon_fp(const char *mol_pkl, size_t mol_pkl_sz,
                               const char *details_json) {
  if (!mol_pkl || !mol_pkl_sz) {
    return nullptr;
  }
  try {
    auto fp = MinimalLib::avalon_fp_as_bitvect(
        mol_from_pkl(mol_pkl, mol_pkl_sz), details_json);
    auto res = BitVectToText(*fp);
    return str_to_c(res);
  } catch (...) {
    return nullptr;
  }
}

extern "C" char *get_avalon_fp_as_bytes(const char *mol_pkl, size_t mol_pkl_sz,
                                        size_t *nbytes,
                                        const char *details_json) {
  if (!mol_pkl || !mol_pkl_sz) {
    return nullptr;
  }
  try {
    auto fp = MinimalLib::avalon_fp_as_bitvect(
        mol_from_pkl(mol_pkl, mol_pkl_sz), details_json);
    auto res = BitVectToBinaryText(*fp);
    return str_to_c(res, nbytes);
  } catch (...) {
    return nullptr;
  }
}
#endif

extern "C" void prefer_coordgen(short val) {
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  RDDepict::preferCoordGen = val;
#endif
};

extern "C" short has_coords(const char *mol_pkl, size_t mol_pkl_sz) {
  short res = 0;
  if (mol_pkl && mol_pkl_sz) {
    auto mol = mol_from_pkl(mol_pkl, mol_pkl_sz);
    res = (mol.getNumConformers() > 0);
    if (res) {
      res = mol.getConformer().is3D() ? 3 : 2;
    }
  }
  return res;
}

extern "C" short set_2d_coords(char **mol_pkl, size_t *mol_pkl_sz) {
  if (!mol_pkl || !mol_pkl_sz || !*mol_pkl || !*mol_pkl_sz) {
    return 0;
  }
  auto mol = mol_from_pkl(*mol_pkl, *mol_pkl_sz);
  RDDepict::compute2DCoords(mol);

  mol_to_pkl(mol, mol_pkl, mol_pkl_sz);
  return 1;
}

extern "C" short set_2d_coords_aligned(char **mol_pkl, size_t *mol_pkl_sz,
                                       const char *template_pkl,
                                       size_t template_sz,
                                       const char *details_json,
                                       char **match_json) {
  if (match_json) {
    *match_json = nullptr;
  }
  if (!mol_pkl || !mol_pkl_sz || !*mol_pkl || !*mol_pkl_sz || !template_pkl ||
      !template_sz || !template_pkl || !template_sz) {
    return 0;
  }
  auto templ = mol_from_pkl(template_pkl, template_sz);
  if (!templ.getNumConformers()) {
    return 0;
  }
  auto mol = mol_from_pkl(*mol_pkl, *mol_pkl_sz);
  auto match = MinimalLib::generate_aligned_coords(mol, templ, details_json);
  if (match.empty()) {
    return 0;
  } else {
    mol_to_pkl(mol, mol_pkl, mol_pkl_sz);
    if (match_json) {
      *match_json = str_to_c(match);
    }
    return 1;
  }
};

extern "C" short set_3d_coords(char **mol_pkl, size_t *mol_pkl_sz,
                               const char *params_json) {
  if (!mol_pkl || !mol_pkl_sz || !*mol_pkl || !*mol_pkl_sz) {
    return 0;
  }
  auto mol = mol_from_pkl(*mol_pkl, *mol_pkl_sz);
  std::string json;
  if (params_json) {
    json = params_json;
  }
  DGeomHelpers::EmbedParameters ps = DGeomHelpers::srETKDGv3;
  if (!json.empty()) {
    DGeomHelpers::updateEmbedParametersFromJSON(ps, json);
  }
  int res = DGeomHelpers::EmbedMolecule(mol, ps);
  if (res >= 0) {
    ++res;
  }
  // if we have a coordMap then be sure to clear up the memory that
  // updateEmbedParametersFromJSON() allocated for it
  if (ps.coordMap) {
    delete ps.coordMap;
  }
  mol_to_pkl(mol, mol_pkl, mol_pkl_sz);
  return (short)res;
}

extern "C" short add_hs(char **mol_pkl, size_t *mol_pkl_sz) {
  if (!mol_pkl || !mol_pkl_sz || !*mol_pkl || !*mol_pkl_sz) {
    return 0;
  }
  auto mol = mol_from_pkl(*mol_pkl, *mol_pkl_sz);
  MolOps::addHs(mol);
  // we don't need the properties that sets:
  for (auto atom : mol.atoms()) {
    if (atom->getAtomicNum() == 1) {
      atom->clearProp(common_properties::isImplicit);
    }
  }

  mol_to_pkl(mol, mol_pkl, mol_pkl_sz);
  return 1;
}

extern "C" short remove_all_hs(char **mol_pkl, size_t *mol_pkl_sz) {
  if (!mol_pkl || !mol_pkl_sz || !*mol_pkl || !*mol_pkl_sz) {
    return 0;
  }
  auto mol = mol_from_pkl(*mol_pkl, *mol_pkl_sz);
  MolOps::removeAllHs(mol);

  mol_to_pkl(mol, mol_pkl, mol_pkl_sz);
  return 1;
}

extern "C" short remove_hs(char **mol_pkl, size_t *mol_pkl_sz,
                           const char *details_json) {
  if (!mol_pkl || !mol_pkl_sz || !*mol_pkl || !*mol_pkl_sz) {
    return 0;
  }
  auto mol = mol_from_pkl(*mol_pkl, *mol_pkl_sz);
  MolOps::RemoveHsParameters ps;
  bool sanitize = true;
  MinimalLib::updateRemoveHsParametersFromJSON(ps, sanitize, details_json);
  MolOps::removeHs(mol, ps, sanitize);

  mol_to_pkl(mol, mol_pkl, mol_pkl_sz);
  return 1;
}

// standardization
namespace {
template <typename T>
short standardize_func(char **mol_pkl, size_t *mol_pkl_sz,
                       const char *details_json, T func) {
  if (!mol_pkl || !mol_pkl_sz || !*mol_pkl || !*mol_pkl_sz) {
    return 0;
  }
  auto mol = mol_from_pkl(*mol_pkl, *mol_pkl_sz);
  std::string json;
  if (details_json) {
    json = details_json;
  }
  std::unique_ptr<RWMol> res(func(mol, json));

  mol_to_pkl(*res, mol_pkl, mol_pkl_sz);
  return 1;
}
}  // namespace
extern "C" short cleanup(char **mol_pkl, size_t *mol_pkl_sz,
                         const char *details_json) {
  return standardize_func(mol_pkl, mol_pkl_sz, details_json,
                          MinimalLib::do_cleanup);
};
extern "C" short normalize(char **mol_pkl, size_t *mol_pkl_sz,
                           const char *details_json) {
  return standardize_func(mol_pkl, mol_pkl_sz, details_json,
                          MinimalLib::do_normalize);
};
extern "C" short canonical_tautomer(char **mol_pkl, size_t *mol_pkl_sz,
                                    const char *details_json) {
  return standardize_func(mol_pkl, mol_pkl_sz, details_json,
                          MinimalLib::do_canonical_tautomer);
};
extern "C" short charge_parent(char **mol_pkl, size_t *mol_pkl_sz,
                               const char *details_json) {
  return standardize_func(mol_pkl, mol_pkl_sz, details_json,
                          MinimalLib::do_charge_parent);
};
extern "C" short reionize(char **mol_pkl, size_t *mol_pkl_sz,
                          const char *details_json) {
  return standardize_func(mol_pkl, mol_pkl_sz, details_json,
                          MinimalLib::do_reionize);
};
extern "C" short neutralize(char **mol_pkl, size_t *mol_pkl_sz,
                            const char *details_json) {
  return standardize_func(mol_pkl, mol_pkl_sz, details_json,
                          MinimalLib::do_neutralize);
};
extern "C" short fragment_parent(char **mol_pkl, size_t *mol_pkl_sz,
                                 const char *details_json) {
  return standardize_func(mol_pkl, mol_pkl_sz, details_json,
                          MinimalLib::do_fragment_parent);
};

// chirality
extern "C" short use_legacy_stereo_perception(short value) {
  short was = Chirality::getUseLegacyStereoPerception();
  Chirality::setUseLegacyStereoPerception(value);
  return was;
}

extern "C" short allow_non_tetrahedral_chirality(short value) {
  short was = Chirality::getAllowNontetrahedralChirality();
  Chirality::setAllowNontetrahedralChirality(value);
  return was;
}

std::unique_ptr<MinimalLib::LoggerStateSingletons>
    MinimalLib::LoggerStateSingletons::d_instance;

extern "C" short enable_logging() {
  return MinimalLib::LogHandle::enableLogging();
}
extern "C" short enable_logger(const char *log_name) {
  return MinimalLib::LogHandle::enableLogging(log_name);
}

extern "C" short disable_logging() {
  return MinimalLib::LogHandle::disableLogging();
}
extern "C" short disable_logger(const char *log_name) {
  return MinimalLib::LogHandle::disableLogging(log_name);
}

extern "C" void *set_log_tee(const char *log_name) {
  return MinimalLib::LogHandle::setLogTee(log_name);
}

extern "C" void *set_log_capture(const char *log_name) {
  return MinimalLib::LogHandle::setLogCapture(log_name);
}

extern "C" short destroy_log_handle(void **log_handle) {
  if (!log_handle || !*log_handle) {
    return 0;
  }
  auto lh = reinterpret_cast<MinimalLib::LogHandle *>(*log_handle);
  delete lh;
  *log_handle = nullptr;
  return 1;
}

extern "C" char *get_log_buffer(void *log_handle) {
  return log_handle
             ? str_to_c(reinterpret_cast<MinimalLib::LogHandle *>(log_handle)
                            ->getBuffer())
             : nullptr;
}

extern "C" short clear_log_buffer(void *log_handle) {
  if (log_handle) {
    reinterpret_cast<MinimalLib::LogHandle *>(log_handle)->clearBuffer();
    return 1;
  }
  return 0;
}

extern "C" short has_prop(const char *mol_pkl, size_t mol_pkl_sz,
                          const char *key) {
  auto mol = mol_from_pkl(mol_pkl, mol_pkl_sz);
  return mol.hasProp(key);
}

extern "C" char **get_prop_list(const char *mol_pkl, size_t mol_pkl_sz,
                                short includePrivate, short includeComputed) {
  auto mol = mol_from_pkl(mol_pkl, mol_pkl_sz);
  auto propList = mol.getPropList(includePrivate, includeComputed);
  std::string propNames;
  for (const auto &prop : propList) {
    propNames += prop + ",";
  }
  auto resLen = sizeof(char *) * (propList.size() + 1);
  char **res = (char **)malloc(resLen);
  if (!res) {
    return nullptr;
  }
  memset(res, 0, resLen);
  for (size_t i = 0; i < propList.size(); ++i) {
    res[i] = strdup(propList.at(i).c_str());
    if (!res[i]) {
      while (i--) {
        free(res[i]);
      }
      return nullptr;
    }
  }
  return res;
}

extern "C" void set_prop(char **mol_pkl, size_t *mol_pkl_sz, const char *key,
                         const char *val, short computed) {
  auto mol = mol_from_pkl(*mol_pkl, *mol_pkl_sz);
  std::string valAsString(val);
  mol.setProp(key, valAsString, computed);
  mol_to_pkl(mol, mol_pkl, mol_pkl_sz);
}

extern "C" char *get_prop(const char *mol_pkl, size_t mol_pkl_sz,
                          const char *key) {
  auto mol = mol_from_pkl(mol_pkl, mol_pkl_sz);
  if (!mol.hasProp(key)) {
    return nullptr;
  }
  std::string val;
  mol.getProp(key, val);
  return strdup(val.c_str());
}

extern "C" short clear_prop(char **mol_pkl, size_t *mol_pkl_sz,
                            const char *key) {
  auto mol = mol_from_pkl(*mol_pkl, *mol_pkl_sz);
  short res = mol.hasProp(key);
  if (res) {
    mol.clearProp(key);
    mol_to_pkl(mol, mol_pkl, mol_pkl_sz);
  }
  return res;
}

extern "C" void keep_props(char **mol_pkl, size_t *mol_pkl_sz,
                           const char *details_json) {
  auto mol = mol_from_pkl(*mol_pkl, *mol_pkl_sz);
  unsigned int propFlags = PicklerOps::PropertyPickleOptions::AllProps ^
                           PicklerOps::PropertyPickleOptions::ComputedProps;
  MinimalLib::updatePropertyPickleOptionsFromJSON(propFlags, details_json);
  mol_to_pkl(mol, mol_pkl, mol_pkl_sz, propFlags);
}

extern "C" short add_mol_to_png_blob(char **png_blob, size_t *png_blob_sz,
                                     const char *pkl, size_t pkl_sz,
                                     const char *details_json) {
  if (!png_blob || !*png_blob || !png_blob_sz || !*png_blob_sz || !pkl ||
      !pkl_sz) {
    return 0;
  }
  PNGMetadataParams params;
  MinimalLib::updatePNGMetadataParamsFromJSON(params, details_json);
  std::string pngString(*png_blob, *png_blob_sz);
  try {
    auto mol = mol_from_pkl(pkl, pkl_sz);
    auto updatedPngString = addMolToPNGString(mol, pngString, params);
    auto updated_png_blob =
        static_cast<char *>(malloc(updatedPngString.size()));
    if (!updated_png_blob) {
      return 0;
    }
    memcpy(updated_png_blob, updatedPngString.data(), updatedPngString.size());
    free(*png_blob);
    *png_blob = updated_png_blob;
    *png_blob_sz = updatedPngString.size();
  } catch (...) {
    return 0;
  }
  return 1;
}

extern "C" short get_mol_from_png_blob(const char *png_blob, size_t png_blob_sz,
                                       char **pkl, size_t *pkl_sz,
                                       const char *details_json) {
  if (!png_blob || !png_blob_sz || !pkl || !pkl_sz) {
    return 0;
  }
  std::string pngString(png_blob, png_blob_sz);
  auto mols = MinimalLib::get_mols_from_png_blob_internal(pngString, true,
                                                          details_json);
  if (mols.empty()) {
    return 0;
  }
  char *pkl_local = nullptr;
  size_t pkl_sz_local = 0;
  mol_to_pkl(*mols.front(), &pkl_local, &pkl_sz_local);
  if (pkl_local && pkl_sz_local) {
    *pkl = pkl_local;
    *pkl_sz = pkl_sz_local;
    return 1;
  }
  return 0;
}

extern "C" short get_mols_from_png_blob(const char *png_blob,
                                        size_t png_blob_sz, char ***pkl_array,
                                        size_t **pkl_sz_array,
                                        const char *details_json) {
  if (!png_blob || !png_blob_sz || !pkl_array || !pkl_sz_array) {
    return 0;
  }
  std::string pngString(png_blob, png_blob_sz);
  auto mols = MinimalLib::get_mols_from_png_blob_internal(pngString, false,
                                                          details_json);
  if (mols.empty()) {
    return 0;
  }
  char **pkl_array_local = nullptr;
  size_t *pkl_sz_array_local = nullptr;
  size_t mol_array_len = mols.size() + 1;
  pkl_array_local = (char **)malloc(mol_array_len * sizeof(char *));
  if (pkl_array_local) {
    memset(pkl_array_local, 0, mol_array_len * sizeof(char *));
    pkl_sz_array_local = (size_t *)malloc(mol_array_len * sizeof(size_t));
  }
  if (pkl_sz_array_local) {
    memset(pkl_sz_array_local, 0, mol_array_len * sizeof(size_t));
    short i = 0;
    for (const auto &mol : mols) {
      mol_to_pkl(*mol, &pkl_array_local[i], &pkl_sz_array_local[i]);
      if (pkl_array_local[i] && pkl_sz_array_local[i]) {
        ++i;
      } else {
        break;
      }
    }
    if (i == static_cast<short>(mols.size())) {
      *pkl_array = pkl_array_local;
      *pkl_sz_array = pkl_sz_array_local;
      return i;
    }
  }
  free_mol_array(&pkl_array_local, &pkl_sz_array_local);
  return 0;
}

extern "C" void free_mol_array(char ***pkl_array, size_t **pkl_sz_array) {
  if (pkl_array && *pkl_array) {
    for (size_t i = 0; (*pkl_array)[i]; ++i) {
      free((*pkl_array)[i]);
      (*pkl_array)[i] = NULL;
    }
    free(*pkl_array);
    *pkl_array = NULL;
  }
  if (pkl_sz_array && *pkl_sz_array) {
    free(*pkl_sz_array);
    *pkl_sz_array = NULL;
  }
}

#if (defined(__GNUC__) || defined(__GNUG__))
#pragma GCC diagnostic pop
#endif
