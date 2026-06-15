//
//  Copyright (C) 2025 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#pragma once

#include <functional>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <GraphMol/MacroMol.h>

#include <RDGeneral/export.h>

namespace RDKit {

class ROMol;

//! Represents a monomer definition in the library
struct RDKIT_MONOMERMOL_EXPORT MonomerEntry {
    std::string symbol;           // e.g., "A" for Alanine
    std::string original_data;    // Original definition (SMILES, SDF, etc.)
    std::shared_ptr<ROMol> mol;   // Parsed molecule
    std::string monomer_class;    // "PEPTIDE", "RNA", etc.
    std::string pdb_code;         // e.g., "ALA" (optional)
};


//! Return type for getMonomerInfo: {symbol, original_data, monomer_class}
using monomer_info_t = std::optional<std::tuple<std::string, std::string, std::string>>;


// --- Query operations ---

//! Get original data (SMILES/SDF/etc) for a monomer by its symbol and class
[[nodiscard]] std::optional<std::string> RDKIT_MONOMERMOL_EXPORT getMonomerData(
    const MacroMolTemplateLib &macroMolTemplateLib,
    const std::string& monomer_id,
    const std::string& monomer_class);

//! Get full monomer info from a three-letter PDB code
/*!
    \param pdb_code The three-letter PDB residue code (e.g., "ALA")
    \return tuple of {symbol, original_data, monomer_class} or nullopt if not found
*/
[[nodiscard]] monomer_info_t RDKIT_MONOMERMOL_EXPORT
getMonomerInfo(const MacroMolTemplateLib &macroMolTemplateLib, const std::string& pdb_code, const std::string monomer_class);

//! Get PDB code for a monomer symbol
[[nodiscard]] std::optional<std::string> RDKIT_MONOMERMOL_EXPORT
getPdbCode(const MacroMolTemplateLib &macroMolTemplateLib, const std::string& symbol, const std::string& monomer_class);

[[nodiscard]] std::optional<std::vector<std::string>> RDKIT_MONOMERMOL_EXPORT
getPdbCodes(const MacroMolTemplateLib &macroMolTemplateLib, const std::string& symbol, const std::string& monomer_class);

// MacroMolTemplateLib &getMacroMolTemplateLib() {
//     return d_macroMolTemplateLib;
// }
// const MacroMolTemplateLib &getMacroMolTemplateLib() const {
//     return d_macroMolTemplateLib;
// }

//! Get parsed molecule for a monomer
//[[nodiscard]] std::shared_ptr<ROMol> getMonomer(
// [[nodiscard]] MacroMolTemplate *getMonomer(
//     const std::string& monomer_id,
//     const std::string& monomer_class) const;

// --- Mutation operations (for instance libraries) ---

//! Add a monomer from a SMILES string
void RDKIT_MONOMERMOL_EXPORT addMonomerFromSmiles(MacroMolTemplateLib &macroMolTemplateLib, 
                            const std::string& smiles,
                            const std::string& symbol,
                            const std::string& monomer_class,
                            const std::string& pdb_code = "");

//! Add a monomer from an SDF/molblock string
void RDKIT_MONOMERMOL_EXPORT addMonomerFromSDF(MacroMolTemplateLib &macroMolTemplateLib,
                        const std::string& sdf_data,
                        const std::string& symbol,
                        const std::string& monomer_class,
                        const std::string& pdb_code = "");

//! Add a monomer with a pre-parsed molecule
void RDKIT_MONOMERMOL_EXPORT addMonomer(MacroMolTemplateLib &macroMolTemplateLib, 
                std::unique_ptr<RWMol> &mol,
                const std::string& data,
                const std::string& symbol,
                const std::string& monomer_class,
                const std::string& pdb_code = "");

// //! Check if a monomer exists
// [[nodiscard]] bool RDKIT_MONOMERMOL_EXPORT hasMonomer(const MacroMolTemplateLib &macroMolTemplateLib, const std::string& symbol,
//                                 const std::string& monomer_class);

// --- Global library configuration ---

//! Get the global singleton instance

//static MonomerLibrary *getGlobalLibrary();

   
class RDKIT_MONOMERMOL_EXPORT GlobalMonomerLibrary{
public:
    static MacroMolTemplateLib *getGlobalLibrary();

 private:
    //! Load built-in monomer definitions
    //! TODO: this may load in a preset json or sqlite DB
    static void loadBuiltinDefinitions();

    //! Static state for global mode
    //static bool s_useGlobalLibrary
    static std::unique_ptr<MacroMolTemplateLib> s_globalLibrary;
    static std::once_flag s_globalLibraryOnce;



};

} // namespace RDKit
