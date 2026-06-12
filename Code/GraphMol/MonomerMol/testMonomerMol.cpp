//
//  Copyright (C) 2025 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <sstream>
#include "RDGeneral/test.h"
#include <catch2/catch_all.hpp>
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MACROMolUtils.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/FileParseException.h>
#include <boost/algorithm/string.hpp>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/Chirality.h>

#include <GraphMol/FileParsers/SequenceParsers.h>

#include <GraphMol/MonomerMol/Conversions.h>
#include <GraphMol/MonomerMol/MonomerMol.h>
#include <GraphMol/MonomerMol/MonomerLibrary.h>
#include <GraphMol/MonomerInfo.h>

#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <string>

using namespace RDKit;

/*
 * Temporary, simple FASTA parser to show how to use the MonomerMol
*/
static std::string to_fasta(const MACROMol& mol)
{
    // make a set of all the chains in the molecule
    std::set<std::string> chains;
    for (auto atom : mol.atoms()) {
        chains.insert(getPolymerId(atom));
    }
    std::stringstream fasta;
    bool first = true;
    for (auto chain_id : chains) {
        if (!first) {
            fasta << "\n";
        }
        first = false;
        fasta << ">Chain " << chain_id << "\n"; // add title
        auto chain = getPolymer(mol,chain_id);
        for (auto atom : chain.atoms) {
            // Really there should be a lookup to ensure that these are all 1 letter
            // And a sort by residue number & connectivity
            const auto* res_info = static_cast<const RDKit::AtomPDBResidueInfo*>(mol.getAtomWithIdx(atom)->getMonomerInfo());
            fasta << res_info->getResidueName();
        }
    }
    return fasta.str();
}

static void neutralize_atoms(RDKit::ROMol& mol)
{
    // Algorithm for neutralizing molecules from
    // https://www.rdkit.org/docs/Cookbook.html#neutralizing-molecules by Noel
    // O'Boyle Will neutralize the molecule by adding or removing hydrogens as
    // needed. This will ensure SMILES can be used to match atomistic structures
    // to the correct monomer.
    static const std::unique_ptr<RDKit::RWMol> neutralize_query(
        RDKit::SmartsToMol(
            "[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]"));
    for (const auto& match : RDKit::SubstructMatch(mol, *neutralize_query)) {
        auto atom = mol.getAtomWithIdx(match[0].second);
        auto chg = atom->getFormalCharge();
        auto hcount = atom->getTotalNumHs();
        atom->setFormalCharge(0);
        atom->setNumExplicitHs(hcount - chg);
        atom->updatePropertyCache();
    }
}

static std::unique_ptr<RDKit::ROMol> resolve_his(const RDKit::ROMol& mol)
{
    // Some structures may contain different protonation states for histidine,
    // but we currently map all of them to the same single letter code 'H' in
    // the monomeric representation. Since we want to test the roundtrip
    // conversion against the original, we need to resolve the histidine
    // protonation state to what is in the monomer database.
    std::string smiles = RDKit::MolToSmiles(mol);

    std::vector<std::string> targets = {"Cc1c[nH]cn1"};
    std::string replace_with = "Cc1cnc[nH]1";
    for (const auto& target : targets) {
        size_t pos = 0;
        while ((pos = smiles.find(target, pos)) != std::string::npos) {
            smiles.replace(pos, target.length(), replace_with);
            pos += replace_with.length();
        }
    }
    return std::unique_ptr<RDKit::ROMol>(RDKit::SmilesToMol(smiles));
}

static void remove_solvents(RDKit::RWMol& rwmol)
{
    std::vector<unsigned int> atoms_to_remove;
    for (const auto& atom : rwmol.atoms()) {
        const auto res_info = static_cast<const RDKit::AtomPDBResidueInfo*>(
            atom->getMonomerInfo());

        if (res_info != nullptr) {
            std::string residue_name = res_info->getResidueName();
            if (residue_name == "HOH" || residue_name == "S04") {
                atoms_to_remove.push_back(atom->getIdx());
            }
        }
    }
    rwmol.beginBatchEdit();
    for (unsigned int atom_idx : atoms_to_remove) {
        rwmol.removeAtom(atom_idx);
    }
    rwmol.commitBatchEdit();
}

static bool same_roundtrip_mol(const RDKit::ROMol& original,
                               const RDKit::ROMol& roundtrip,
                               unsigned int num_chains) {
    // Verify the roundtrip produces a molecule with expected number of atoms,
    // roundtripped structure has additional oxygen atoms at the end of the
    // chains so there may be some extra atoms in the roundtrip (up to 2 *
    // number of chains)
    RDKit::SubstructMatchParameters params;
    params.maxMatches = 1;
    auto match = RDKit::SubstructMatch(roundtrip,
                                       original, params);
    if (match.empty()) {
      return false;
    }
    if (match[0].size() != original.getNumAtoms()) {
      return false;
    }
    return (roundtrip.getNumAtoms() - original.getNumAtoms()) <= (2 * num_chains);
}

TEST_CASE("FASTAConversions") {
  SECTION("SIMPLE") {
    // Build MonomerMol with a single chain
    MACROMol macroMol;
    addGlobalLibrary(macroMol);

    addMonomer(macroMol, "R", 1, "PEPTIDE", "PEPTIDE1");
    addMonomer(macroMol, "D");
    addMonomer(macroMol, "K");
    addMonomer(macroMol, "I");
    addMonomer(macroMol, "T");

    addConnection(macroMol,0, 1, BACKBONE_LINKAGE);
    addConnection(macroMol,1, 2, BACKBONE_LINKAGE);
    addConnection(macroMol,2, 3, BACKBONE_LINKAGE);
    addConnection(macroMol,3, 4, BACKBONE_LINKAGE);

    CHECK(macroMol.getNumAtoms() == 5);
    CHECK(macroMol.getNumBonds() == 4);
    CHECK(std::string(">Chain PEPTIDE1\nRDKIT") == to_fasta(macroMol));
  }

  SECTION("MultipleChains") {
    // Build MonomerMol with two chains
    MACROMol macroMol;
    addGlobalLibrary(macroMol);

    auto midx1 = addMonomer(macroMol, "R", 1, "PEPTIDE", "A");
    auto midx2 = addMonomer(macroMol, "D");

    auto midx3 = addMonomer(macroMol, "K", 1, "PEPTIDE", "B");
    auto midx4 = addMonomer(macroMol, "I");
    auto midx5 = addMonomer(macroMol, "T");

    addConnection(macroMol, midx1, midx2, BACKBONE_LINKAGE);
    addConnection(macroMol, midx3, midx4, BACKBONE_LINKAGE);
    addConnection(macroMol, midx4, midx5, BACKBONE_LINKAGE);

    CHECK(macroMol.getNumAtoms() == 5);
    CHECK(macroMol.getNumBonds() == 3);

    CHECK(std::string(">Chain A\nRD\n>Chain B\nKIT") == to_fasta(macroMol));
  }

  SECTION("UsingSequenceReader") {
    std::string seq = "CGCGAATTACCGCG";
    // the sequence parser creates an atomistic molecule
    auto atomistic_mol = SequenceToMol(seq);
    CHECK(atomistic_mol);

    // PDB info is used to convert the atomistic sturcture into a MonomerMol
    auto macroMol = toMonomeric(*atomistic_mol);

    CHECK(std::string(">Chain PEPTIDE1\nCGCGAATTACCGCG") == to_fasta(*macroMol));
    delete atomistic_mol;
  }
}

TEST_CASE("Conversions") {
  SECTION("toAtomistic") {
    std::string seq = "CGCGA";
    auto atomistic_mol = SequenceToMol(seq);

    MACROMol macroMol;
    addGlobalLibrary(macroMol);
    auto midx1 = addMonomer(macroMol, "C", 1, "PEPTIDE", "PEPTIDE1");
    auto midx2 = addMonomer(macroMol, "G");
    auto midx3 = addMonomer(macroMol, "C");
    auto midx4 = addMonomer(macroMol, "G");
    auto midx5 = addMonomer(macroMol, "A");

    addConnection(macroMol, midx1, midx2, BACKBONE_LINKAGE);
    addConnection(macroMol, midx2, midx3, BACKBONE_LINKAGE);
    addConnection(macroMol, midx3, midx4, BACKBONE_LINKAGE);
    addConnection(macroMol, midx4, midx5, BACKBONE_LINKAGE);
    auto atomistic_mol2 = toAtomistic(macroMol);
    
    // atomistic structure is same as using sequence parser
    std::string smi1 = MolToSmiles(*atomistic_mol);
    std::string smi2 = MolToSmiles(*atomistic_mol2);
    CHECK(smi1 == smi2);

    delete atomistic_mol;
  }

  SECTION("toAtomisticWithBranch") {
    // This is equivalent to HELM string "PEPTIDE1{A.D(C)P}$$$$"
    MACROMol macroMol;
    addGlobalLibrary(macroMol);
    auto midx1 = addMonomer(macroMol, "A", 1, "PEPTIDE", "PEPTIDE1");
    auto midx2 = addMonomer(macroMol, "D");
    auto midx3 = addMonomer(macroMol, "C");
    auto midx4 = addMonomer(macroMol, "P");

    addConnection(macroMol, midx1, midx2, BACKBONE_LINKAGE);
    addConnection(macroMol, midx2, midx3, BRANCH_LINKAGE);
    addConnection(macroMol, midx2, midx4, BACKBONE_LINKAGE);

    std::string smi = MolToSmiles(*toAtomistic(macroMol));
    CHECK(smi == "C[C@H](N)C(=O)N[C@@H](CC(=O)N[C@@H](CS)C(=O)O)C(=O)N1CCC[C@H]1C(=O)O");
  }


  SECTION("toAtomisticWithDisulfide") {
    std::string helm = "PEPTIDE1{C.A.A.A.C}$PEPTIDE1,PEPTIDE1,1:R3-5:R3$$$V2.0";
    auto atomistic_mol = HELMToMol(helm);

    MACROMol macroMol;
    addGlobalLibrary(macroMol);
    
    auto midx1 = addMonomer(macroMol, "C", 1, "PEPTIDE", "PEPTIDE1");
    auto midx2 = addMonomer(macroMol, "A");
    auto midx3 = addMonomer(macroMol, "A");
    auto midx4 = addMonomer(macroMol, "A");
    auto midx5 = addMonomer(macroMol, "C");

    addConnection(macroMol, midx1, midx2, BACKBONE_LINKAGE);
    addConnection(macroMol, midx2, midx3, BACKBONE_LINKAGE);
    addConnection(macroMol, midx3, midx4, BACKBONE_LINKAGE);
    addConnection(macroMol, midx4, midx5, BACKBONE_LINKAGE);
    addConnection(macroMol, midx1, midx5, CROSS_LINKAGE);

    std::string smi1 = MolToSmiles(*atomistic_mol);
    std::string smi2 = MolToSmiles(*toAtomistic(macroMol));
    CHECK(smi1 == smi2);
    delete atomistic_mol;
  }

  SECTION("toAtomisticCyclicPeptide") {
    std::string helm = "PEPTIDE1{F.Y.K.A.R.L}$PEPTIDE1,PEPTIDE1,6:R2-1:R1$$$V2.0";
    auto atomistic_mol = HELMToMol(helm);

    MACROMol macroMol;
    addGlobalLibrary(macroMol);
    auto midx1 = addMonomer(macroMol, "F", 1, "PEPTIDE", "PEPTIDE1");
    auto midx2 = addMonomer(macroMol, "Y");
    auto midx3 = addMonomer(macroMol, "K");
    auto midx4 = addMonomer(macroMol, "A");
    auto midx5 = addMonomer(macroMol, "R");
    auto midx6 = addMonomer(macroMol, "L");
    addConnection(macroMol, midx1, midx2, BACKBONE_LINKAGE);
    addConnection(macroMol, midx2, midx3, BACKBONE_LINKAGE);
    addConnection(macroMol, midx3, midx4, BACKBONE_LINKAGE);
    addConnection(macroMol, midx4, midx5, BACKBONE_LINKAGE);
    addConnection(macroMol, midx5, midx6, BACKBONE_LINKAGE);
    addConnection(macroMol, midx6, midx1, BACKBONE_LINKAGE);

    CHECK(macroMol.getNumAtoms() == 6);
    CHECK(macroMol.getNumBonds() == 6);


    
    std::string smi1 = MolToSmiles(*atomistic_mol);
    auto mol2 = toAtomistic(macroMol);
    std::string smi2 = MolToSmiles(*(mol2.get()));
    CHECK(smi1 == smi2);
    delete atomistic_mol;
  }

  SECTION("toAtomisticCyclicPeptideOutOfOrder") {
    std::string helm = "PEPTIDE1{F.Y.K.A.R.L}$PEPTIDE1,PEPTIDE1,6:R2-1:R1$$$V2.0";
    auto atomistic_mol = HELMToMol(helm);


    // toAtomistic should still work even when monomers are added out of order
    MACROMol macroMol;
    addGlobalLibrary(macroMol);
    auto midx1 = addMonomer(macroMol, "A", 1, "PEPTIDE", "PEPTIDE1");
    auto midx2 = addMonomer(macroMol, "K");
    auto midx3 = addMonomer(macroMol, "Y");
    auto midx4 = addMonomer(macroMol, "F");
    auto midx5 = addMonomer(macroMol, "L");
    auto midx6 = addMonomer(macroMol, "R");
    addConnection(macroMol, midx4, midx3, BACKBONE_LINKAGE);
    addConnection(macroMol, midx3, midx2, BACKBONE_LINKAGE);
    addConnection(macroMol, midx2, midx1, BACKBONE_LINKAGE);
    addConnection(macroMol, midx1, midx6, BACKBONE_LINKAGE);
    addConnection(macroMol, midx6, midx5, BACKBONE_LINKAGE);
    addConnection(macroMol, midx5, midx4, BACKBONE_LINKAGE);
    assignChains(macroMol);

    CHECK(macroMol.getNumAtoms() == 6);
    CHECK(macroMol.getNumBonds() == 6);

    std::string smi1 = MolToSmiles(*atomistic_mol);
    std::string smi2 = MolToSmiles(*toAtomistic(macroMol));
    CHECK(smi1 == smi2);
    delete atomistic_mol;
  }

  SECTION("toAtomisticSmilesMonomer") {
    // This is equivelent to HELM string "PEPTIDE1{A.[O=C([C@H]1CCCN1[*:1])[*:2]].P}$$$$"
    MACROMol macroMol;
    addGlobalLibrary(macroMol);
    auto midx1 = addMonomer(macroMol, "A", 1, "PEPTIDE", "PEPTIDE1");
    auto midx2 = addMonomer(macroMol, "O=C([C@H]1CCCN1[*:1])[*:2]", MonomerType::SMILES);
    auto midx3 = addMonomer(macroMol, "P");
    addConnection(macroMol, midx1, midx2, BACKBONE_LINKAGE);
    addConnection(macroMol, midx2, midx3, BACKBONE_LINKAGE);
    assignChains(macroMol);
    CHECK(macroMol.getNumAtoms() == 3);
    CHECK(macroMol.getNumBonds() == 2);
    auto atomistic_mol = toAtomistic(macroMol);
    std::string smi1 = MolToSmiles(*atomistic_mol);
    CHECK(smi1 == "C[C@H](N)C(=O)N1CCC[C@@H]1C(=O)N1CCC[C@H]1C(=O)O");
  }

  SECTION("toMonomeric1DNG") {
    std::string pdbfile = getenv("RDBASE");
    pdbfile += "/Code/GraphMol/MonomerMol/test_data/1dng.pdb";
    auto mol = std::unique_ptr<RWMol>(PDBFileToMol(pdbfile));
    auto monomer_mol = toMonomeric(*mol);

    CHECK(std::string(">Chain PEPTIDE1\nQAPAYEEAAEELAKS") == to_fasta(*monomer_mol));
    CHECK(monomer_mol->getNumAtoms() == 15);
    CHECK(monomer_mol->getNumBonds() == 14);

    auto roundtrip = toAtomistic(*monomer_mol);
    neutralize_atoms(*mol);
    neutralize_atoms(*roundtrip);
    CHECK(same_roundtrip_mol(*mol, *roundtrip, getPolymerIds(*monomer_mol.get()).size()));
  }

  SECTION("toMonomeric2N65") {
    // Example with multiple chains and a disulfide bond between chains
    std::string pdbfile = getenv("RDBASE");
    pdbfile += "/Code/GraphMol/MonomerMol/test_data/2n65.pdb";
    auto mol = std::unique_ptr<RWMol>(PDBFileToMol(pdbfile));
    auto monomer_mol = toMonomeric(*mol);

    CHECK(std::string(">Chain PEPTIDE1\nVARGWKRKCPLFGKGG\n>Chain PEPTIDE2\nVARGWKRKCPLFGKGG") == to_fasta(*monomer_mol));
    CHECK(monomer_mol->getNumAtoms() == 32);
    CHECK(monomer_mol->getNumBonds() == 31);

    auto roundtrip = toAtomistic(*monomer_mol);
    neutralize_atoms(*mol);
    auto resolved_his = resolve_his(*mol);
    neutralize_atoms(*roundtrip);
    CHECK(same_roundtrip_mol(*resolved_his, *roundtrip, getPolymerIds(*monomer_mol).size()));
  }

  SECTION("toMonomeric4QAF") {
    // Example with disulfide bond between chains and within chains
    std::string pdbfile = getenv("RDBASE");
    pdbfile += "/Code/GraphMol/MonomerMol/test_data/4qaf.pdb";
    auto mol = std::unique_ptr<RWMol>(PDBFileToMol(pdbfile));
    remove_solvents(*mol); // we don't care about water or S04

    auto monomer_mol = toMonomeric(*(mol.get()));

    // There should be 4 chains
    auto polymer_ids = getPolymerIds(*monomer_mol);
    CHECK(polymer_ids.size() == 4);

    CHECK(monomer_mol->getNumAtoms() == 395);
    CHECK(monomer_mol->getNumBonds() == 392);

    auto roundtrip = toAtomistic(*monomer_mol);
    neutralize_atoms(*mol);
    neutralize_atoms(*roundtrip);
    CHECK(same_roundtrip_mol(*mol, *roundtrip, polymer_ids.size()));
  }

  SECTION("toMonomeric5VAV") {
    // Example with disulfide bond and R2-R1 cycle closure
    std::string pdbfile = getenv("RDBASE");
    pdbfile += "/Code/GraphMol/MonomerMol/test_data/5vav.pdb";
    auto mol = std::unique_ptr<RWMol>(PDBFileToMol(pdbfile));
    auto monomer_mol = toMonomeric(*mol);
    
    remove_solvents(*mol); // we don't care about water or S04

    // Checking sequence doesn't make sense due to SMILES monomer
    // Two extra bonds than atoms because this structure has two cycle closures, one
    // by a disulfide bond and another from an R2-R1 closure
    CHECK(monomer_mol->getNumAtoms() == 14);
    CHECK(monomer_mol->getNumBonds() == 15);

    auto roundtrip = toAtomistic(*monomer_mol);
    neutralize_atoms(*mol);
    neutralize_atoms(*roundtrip);
    CHECK(same_roundtrip_mol(*mol, *roundtrip, getPolymerIds(*monomer_mol).size()));
  }
}
