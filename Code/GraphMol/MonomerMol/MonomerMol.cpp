//
//  Copyright (C) 2025 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "MonomerMol.h"
#include "MonomerLibrary.h"

#include <GraphMol/QueryAtom.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SubstanceGroup.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/MolPickler.h>

#include <boost/functional/hash.hpp>

namespace RDKit
{

// Helper functions
namespace {

std::pair<unsigned int, unsigned int> getAttchpts(const std::string& linkage)
{
    // in form RX-RY, returns {X, Y}
    auto dash = linkage.find('-');
    if (dash == std::string::npos) {
        throw std::runtime_error("Invalid linkage format: " + linkage);
    }
    return {std::stoi(linkage.substr(1, dash - 1)),
            std::stoi(linkage.substr(dash + 2))};
}

} // anonymous namespace


size_t addMonomer(MacroMol &macroMol, std::string_view name, int residue_number, std::string_view monomer_class,
                              std::string_view chain_id, MonomerType monomer_type)

{
    // see if this is an "immediate mode" macro ref (smiles type)
    //if so, create the template on the fly and assign the once-name to the new reference


    std::string n{name};
    std::string monomerClass(monomer_class);
    //if (is_smiles) {
    if (monomer_type == MonomerType::SMILES) {

      // new name
      std::string smiles(name);

      n = "SM" + std::to_string(macroMol.getTemplateLibrary()->size() + 1);

      addMonomerFromSmiles(*macroMol.getTemplateLibrary(), smiles, n, monomerClass); // owned by the local lib
    }


    //auto a = std::make_unique<::RDKit::Atom>();
    auto newAtomId = macroMol.addMacroAtom(monomerClass, n);
    auto a = macroMol.getAtomWithIdx(newAtomId);
    
    // Allows monomer names show up in image renderings
    a->setProp(RDKit::common_properties::atomLabel, n);
    // Allows monomer names to be written to SMILES
    a->setProp(RDKit::common_properties::smilesSymbol, n);
    // Always start with BRANCH_MONOMER as false, will be set to
    // true if branch linkage is made to this monomer. Will be used
    // for rendering and monomer ordering purposes.
    //a->setProp(BRANCH_MONOMER, false);
    // Provided monomer can be an unnamed monomer represented by a SMILES string,
    // indicated by is_smiles=True
    //a->setProp(SMILES_MONOMER, is_smiles);

    // Give canonicalization to monomers based on name
    // TODO: Canonicalize on name and class?
    static boost::hash<std::string> hasher;
    a->setIsotope(hasher(n));
    a->setNoImplicit(true);

    // This is how the information about this monomer is represented
    auto* monomer_info = new ::RDKit::AtomMonomerInfo();
    monomer_info->setResidueNumber(residue_number);
    monomer_info->setResidueName(n);
    monomer_info->setChainId(std::string{chain_id});
    monomer_info->setMonomerClass(monomerClass);
    a->setMonomerInfo(monomer_info);
    return newAtomId;
}

void addConnection(MacroMol &macroMol, size_t monomer1, size_t monomer2,
                               const std::string& linkage, Bond::BondType bond_type)
{   

    std::vector<std::string> connectionLabels;
 
    auto atom1 = macroMol.getAtomWithIdx(monomer1);
    auto atom2 = macroMol.getAtomWithIdx(monomer2);

    // if either atom is a monomer, there must be a linkage parameter

    if (isMonomer(atom1 )|| isMonomer(atom2)) {

        // bond must be single or hydrogen if it is between two monomers or is a monomer-atom bond
        if (bond_type != Bond::BondType::SINGLE && bond_type == Bond::BondType::HYDROGEN) {
            throw std::runtime_error(
                    "Bond type to any Monomer is not SINGLE or HYDROGEN for bond between atom=" +
                    std::to_string(monomer1) + " and atom=" + std::to_string(monomer2));
        };

        // linkage must be of the form. Ln1-ln2 .  For now one of:  R2-R1, R3-R1, R3-R3.
        // For a hydrogen bond, the bond type is HYDROGEN.

        // get the connect points from the linkage

        boost::algorithm::split(connectionLabels, linkage, boost::algorithm::is_any_of("-"));

        if (connectionLabels.size() != 2) {
                throw std::runtime_error(
                    "Linkage " + linkage + "is not of the form \"L1-L2\" for bond between atom=" +
                    std::to_string(monomer1) + " and atom=" + std::to_string(monomer2));
        }

        // make sure there are not two bonds with the same connection labels
        
        for (auto bondI : boost::make_iterator_range(macroMol.getAtomBonds(atom1))) {
            auto bond = macroMol[bondI];

            if (bond->getOtherAtom(atom1) != atom2) {
                continue;  // only looking at bonds between the same two atoms
            }

            std::string tempProp1="", tempProp2= "";
            if ((bond->getPropIfPresent(common_properties::_MolFileBondAttachPt1,tempProp1) &&
                (tempProp1 == connectionLabels[0] || tempProp1 == connectionLabels[1])) ||
                (bond->getPropIfPresent(common_properties::_MolFileBondAttachPt2,tempProp2) &&
                (tempProp1 == connectionLabels[0] || tempProp1 == connectionLabels[1]))) {
                throw std::runtime_error(
                    "Bond with linkage " + connectionLabels[0] + " or "  + connectionLabels[1] + " already found for bond between atom=" +
                    std::to_string(monomer1) + " and atom=" + std::to_string(monomer2));
            }
        }
    }

    const auto new_total = macroMol.addBond(monomer1, monomer2, bond_type);
    auto bond = macroMol.getBondWithIdx(new_total - 1);

    //if (isMonomer(macroMol.getAtomWithIdx(monomer1)) && connectionLabels[0] != "ATOM") {
        bond->setProp(common_properties::_MolFileBondAttachPt1, connectionLabels[0]);
    //}
    //if (isMonomer(macroMol.getAtomWithIdx(monomer2)) && connectionLabels[1] != "ATOM") {
        bond->setProp(common_properties::_MolFileBondAttachPt2, connectionLabels[1]);
    //}
}

std::vector<std::string> getPolymerIds(MacroMol &macroMol)
{
    std::vector<std::string> polymer_ids;
    for (auto atom : macroMol.atoms()) {
        auto id = getPolymerId(atom);
        // in vector to preseve order of polymers
        if (std::find(polymer_ids.begin(), polymer_ids.end(), id) ==
            polymer_ids.end()) {
            polymer_ids.push_back(id);
        }
    }
    return polymer_ids;
}

size_t addMonomer(MacroMol &macroMol,std::string_view name, MonomerType monomer_type)
{
    if (macroMol.getNumAtoms() == 0) {
        throw std::invalid_argument(
            "No atoms in molecule to determine chain ID");
    }
    const auto* last_monomer = macroMol.getAtomWithIdx(macroMol.getNumAtoms() - 1);
    const auto chain_id = getPolymerId(last_monomer);
    const auto residue_number = getResidueNumber(last_monomer) + 1;
    const auto monomer_class = last_monomer->getMonomerInfo()->getMonomerClass();
    return addMonomer(macroMol, name, residue_number, monomer_class, chain_id, monomer_type);
}

Chain getPolymer(const MacroMol &macroMol, std::string_view polymer_id)
{
    std::vector<unsigned int> chain_atoms;
    for (auto atom : macroMol.atoms()) {
        if (getPolymerId(atom) == polymer_id) {
            chain_atoms.push_back(atom->getIdx());
        }
    }
    std::sort(chain_atoms.begin(), chain_atoms.end(),
              [macroMol](unsigned int a, unsigned int b) {
                  return getResidueNumber(macroMol.getAtomWithIdx(a)) <
                         getResidueNumber(macroMol.getAtomWithIdx(b));
              });
    std::vector<unsigned int> chain_bonds;
    for (auto bond : macroMol.bonds()) {
        if (getPolymerId(bond->getBeginAtom()) == polymer_id &&
            getPolymerId(bond->getEndAtom()) == polymer_id) {
            chain_bonds.push_back(bond->getIdx());
        }
    }

    // Annotations stored as a COP substance group ("heavy chain", "light chain")
    std::string annotation{};
    for (const auto& sg : ::RDKit::getSubstanceGroups(macroMol)) {
        if ((sg.getProp<std::string>("TYPE") != "COP") ||
            !sg.hasProp(ANNOTATION) || !sg.hasProp("ID")) {
            continue;
        }
        if (sg.getProp<std::string>("ID") == polymer_id) {
            annotation = sg.getProp<std::string>(ANNOTATION);
            break;
        }
    }
    return {chain_atoms, chain_bonds, annotation};
}

void addGlobalLibrary(MacroMol &macroMol) {
    
    auto lobalLibrary = GlobalMonomerLibrary::getGlobalLibrary();
    macroMol.getTemplateLibrary()->copyTemplateLib(*lobalLibrary);

  }

} // namespace RDKit
