/* -------------------------------------------------------------------------
 * Declares tools for Monomeric MacroMols
 *
 * A "MonomerMol" is an MacroMols with adornments for HELM like monomeric molecules. Chains
 * are represented via the AtomMonomerInfo structs on atoms, and linkages are
 * stord as a LINKAGE property on bonds in the form of RX-RY, where X is the attachment
 * point used on the begin monomer and Y is the attachment point used on the end monomer.
 *
 * Linkages between an atom and a monomer are represented similarly, but with just RX where
 * X is the attachment point on the monomer.
 *
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#pragma once

#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include <RDGeneral/BetterEnums.h>
#include <RDGeneral/export.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/MACROMol.h>
#include <GraphMol/MonomerMol/MonomerLibrary.h>

namespace RDKit
{

// Forward declarations
class Atom;
class Bond;
class SubstanceGroup;

const std::string LINKAGE{"attachmentPoints"};
const std::string EXTRA_LINKAGE{"extraAttachmentPoints"};
const std::string ATOM_LABEL{"atomLabel"};

// Some default linkage options
const std::string BRANCH_LINKAGE{"R3-R1"};
const std::string BACKBONE_LINKAGE{"R2-R1"};
const std::string CROSS_LINKAGE{"R3-R3"};
const std::string BRANCH2_LINKAGE{"R1-R3"};
const std::string BACKBONE2_LINKAGE{"R1-R2"};
const std::string MONOMERATOM_LINKAGE{"R1-ATOM"};
const std::string MONOMERATOM2_LINKAGE{"R2-ATOM"};
const std::string MONOMERATOM3_LINKAGE{"R3-ATOM"};
const std::string MONOMERATOM4_LINKAGE{"RTOM-R1"};
const std::string MONOMERATOM5_LINKAGE{"RTOM-R2"};
const std::string MONOMERATOM6_LINKAGE{"ATOM}-R3"};

const std::string HYDROGEN_LINKAGE{"pair-pair"};

// Monomer properties stored on atoms
const std::string BRANCH_MONOMER{"isBranchMonomer"};
const std::string SMILES_MONOMER{"isSmilesMonomer"};

// Substance group property to indicate an annotation on a chain
const std::string ANNOTATION{"annotation"};

// acts as a view to a polymer chain
struct Chain {
    std::vector<unsigned int> atoms;
    std::vector<unsigned int> bonds;
    std::string annotation;
};

// Types of polymer chains
enum class ChainType { PEPTIDE, RNA, DNA, CHEM, OTHER };
enum class MonomerType { REGULAR, SMILES };

// Free utility functions for working with monomer atoms

// Returns true if the atom represents a monomer in a MonomerMol, false otherwise
RDKIT_MONOMERMOL_EXPORT bool isMonomer(const Atom* atom);

// Get the polymer/chain ID for an atom
RDKIT_MONOMERMOL_EXPORT std::string getPolymerId(const Atom* atom);

// Get the residue number for an atom
RDKIT_MONOMERMOL_EXPORT unsigned int getResidueNumber(const Atom* atom);


void RDKIT_MONOMERMOL_EXPORT addGlobalLibrary(MACROMol &macroMol);

 
/*
  * Add a monomer to the molecule
  *
  * @param name The name of the monomer
  * @param residue_number The residue number of the monomer
  * @param chain_id The chain ID of the monomer, defaults to "A"
  * @param monomer_type The type of monomer to add
  *
  * @return The index of the added monomer
  */
size_t RDKIT_MONOMERMOL_EXPORT addMonomer(MACROMol &macroMol, std::string_view name, int residue_number,
                  std::string_view monomer_class,
                  std::string_view chain_id = "A",
                  MonomerType monomer_type = MonomerType::REGULAR);

/*
  * Add a monomer to the molecule. Overload that uses the last monomer
  * added to the molecule to determine the chain ID, residue number, and monomer class.
  *
  * @param name The name of the monomer
  * @param monomer_type The type of monomer to add
  *
  * @return The index of the added monomer
  */
size_t RDKIT_MONOMERMOL_EXPORT addMonomer(MACROMol &macroMol, std::string_view name,
                  MonomerType monomer_type = MonomerType::REGULAR);

  // ---- Connection Operations ----

/*
  * Add a connection between two monomers in the molecule. The connection has
  * directionality that starts at monomer1 and ends at monomer2.
  *
  * @param monomer1 The index of the first monomer
  * @param monomer2 The index of the second monomer
  * @param connection_type The type of connection to add
  */
void RDKIT_MONOMERMOL_EXPORT addConnection(MACROMol &macroMol, size_t monomer1, size_t monomer2,
                    const std::string &linkage, Bond::BondType bond_type =
                                Bond::BondType::SINGLE);


  // ---- Query Operations ----

  [[nodiscard]] Chain RDKIT_MONOMERMOL_EXPORT getPolymer(const MACROMol &macroMol, std::string_view polymer_id);

  [[nodiscard]] std::vector<std::string> RDKIT_MONOMERMOL_EXPORT getPolymerIds(MACROMol &macroMol);

  // ---- Chain Assignment ----

  // Discards existing chains and reassigns monomers to sequential chains where monomers
  // are reordered based on connectivity.
  void RDKIT_MONOMERMOL_EXPORT assignChains(MACROMol &macroMol);

  

} // namespace RDKit
