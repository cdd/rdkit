//
//  Copyright (C) 2025 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD licenseFget
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "RDGeneral/test.h"
#include <catch2/catch_all.hpp>

#include <GraphMol/RWMol.h>
#include <GraphMol/MonomerMol/MonomerLibrary.h>
#include <GraphMol/MonomerMol/MonomerMol.h>
#include <GraphMol/MonomerMol/Conversions.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace RDKit;

TEST_CASE("MonomerLibrary") {
  SECTION("GlobalLibraryBasics") {
    // The global library is available by default
    //CHECK(MonomerLibrary::isUsingGlobalLibrary() == true);

    // Get a library with the globals loaded
    MACROMolTemplateLib *globalLib = GlobalMonomerLibrary::getGlobalLibrary();

    // Check that built-in monomers are available
    CHECK(globalLib->contains("PEPTIDE", "A"));
    CHECK(globalLib->contains("PEPTIDE", "G"));
    CHECK(globalLib->contains("PEPTIDE", "C"));

    // Check that non-existent monomers return false
    CHECK(globalLib->contains("PEPTIDE", "XYZ") == false);
    CHECK(globalLib->contains("RNA", "A") == false);

    // Get monomer data
    auto alanineData = getMonomerData(*globalLib, "A", "PEPTIDE");
    REQUIRE(alanineData.has_value());
    CHECK(alanineData->find("C[C@H]") != std::string::npos);  // Alanine SMILES

    // Get monomer info from PDB code
    auto alaInfo = getMonomerInfo(*globalLib, "ALA", "PEPTIDE");
    REQUIRE(alaInfo.has_value());
    CHECK(std::get<0>(*alaInfo) == "A");  // symbol
    CHECK(std::get<2>(*alaInfo) == "PEPTIDE");  // class

    // Get PDB code from symbol
    auto pdbCode = getPdbCode(*globalLib, "A", "PEPTIDE");
    REQUIRE(pdbCode.has_value());
    CHECK(*pdbCode == "ALA");
  }

  SECTION("GlobalLibraryWithMonomerMol") {
    // When no custom library is set, MonomerMol uses the global library
    MACROMol mol; 
    addGlobalLibrary(mol);

    // getGlobalLibrary() returns the global library
    //CHECK(lib.contains("A", "PEPTIDE"));

    // Build a simple peptide using the global library
    addMonomer(mol, "A", 1, "PEPTIDE", "PEPTIDE1");
    addMonomer(mol, "G");
    addMonomer(mol, "C");
    addConnection(mol, 0, 1, BACKBONE_LINKAGE);
    addConnection(mol, 1, 2, BACKBONE_LINKAGE);

    CHECK(mol.getNumAtoms() == 3);
    CHECK(mol.getNumBonds() == 2);

    // Convert to atomistic - uses global library
    auto atomistic = toAtomistic(mol);
    CHECK(atomistic->getNumAtoms() > 0);
  }

  SECTION("CustomLibraryViaConstructor") {
    // Create a custom library with built-ins plus a non-standard monomer
    auto customLib = std::make_shared<MACROMolTemplateLib>(); 

    // Add a custom monomer (fictitious amino acid "X")
    addMonomerFromSmiles(*customLib,
        "CC(C)(C)[C@H](N[H:1])C(=O)[OH:2]",  // tert-butyl glycine
        "X",
        "PEPTIDE",
        "TBG"
    );
    CHECK(customLib->contains( "PEPTIDE", "X"));

    // Create MonomerMol with the custom library
    MACROMol mol;
    mol.getTemplateLibrary()->copyTemplateLib(*customLib.get());
    //mol.addGlobalLibrary();

    auto &internalLib = *mol.getTemplateLibrary();
    addMonomerFromSmiles(internalLib,
        "CC(C)[C@H](N[H:1])C(=O)[OH:2]",  // i-propyl glycine
        "XX",
        "PEPTIDE",
        "IPG"
    );

    // CHECK(mol.hasCustomLibrary() == true);
    // CHECK(mol.usingGlobalLibrary() == true);
    // CHECK(mol.hasLocalTemplates() == true);

    // The custom monomer is accessible
    auto lib = mol.getTemplateLibrary();
    CHECK(lib->contains("PEPTIDE", "X"));
    CHECK(internalLib.contains("PEPTIDE", "XX"));

    // Build a peptide with the custom monomer
    addMonomer(mol, "XX", 1, "PEPTIDE", "PEPTIDE1");
    addMonomer(mol, "X");  // custom monomer
    addMonomer(mol, "XX");  // customer monomer from internal library
    addMonomer(mol, "X");
    addConnection(mol, 0, 1, BACKBONE_LINKAGE);
    addConnection(mol, 1, 2, BACKBONE_LINKAGE);
    addConnection(mol, 2, 3, BACKBONE_LINKAGE);

    CHECK(mol.getNumAtoms() == 4);

    // Convert to atomistic - uses custom library
    auto atomistic = toAtomistic(mol);
    CHECK(atomistic->getNumAtoms() > 0);

    // Verify the custom monomer was expanded correctly
    std::string smiles = MolToSmiles(*atomistic);
    CHECK(smiles.find("C(C)(C)") != std::string::npos);  // tert-butyl group
  }

  SECTION("CustomLibraryViaSetMonomerLibrary") {
    // Create a MonomerMol first (uses global by default)
    MACROMol mol;
    //CHECK(mol.hasLocalTemplates() == false);

    // Create and set a custom library
    auto customLib = std::make_shared<MACROMolTemplateLib>();  // load built-ins
    addMonomerFromSmiles(*customLib.get(),
        "NCCC[C@H](N[H:1])C(=O)[OH:2]",  // ornithine-like
        "Z",
        "PEPTIDE",
        "ORN"
    );


    mol.getTemplateLibrary()->copyTemplateLib(*customLib.get());
    //mol.addGlobalLibrary();


    // CHECK(mol.hasLocalTemplates() == false);
    // CHECK(mol.hasCustomLibrary() == true);

    // Now the custom monomer is available
    CHECK(mol.getTemplateLibrary()->contains("PEPTIDE", "Z"));

    // Build peptide with custom monomer
    addMonomer(mol, "Z", 1, "PEPTIDE", "PEPTIDE1");
    //mol.addMonomer("A");
    addMonomer(mol, "Z");
    addConnection(mol, 0, 1, BACKBONE_LINKAGE);

    auto atomistic = toAtomistic(mol);
    CHECK(atomistic->getNumAtoms() > 0);
  }

  // SECTION("ClearCustomLibrary") {
  //   auto customLib = std::make_shared<MonomerLibrary>();
  //   customLib->addMonomerFromSmiles("CC", "Y", "PEPTIDE");

  //   MonomerMol mol(customLib.get());
  //   CHECK(mol.hasCustomLibrary() == true);
  //   CHECK(mol.getCustomLibrary()->getMACROMolTemplateLib().getNumTemplates() > 0);

  //   // Clear 
  //   mol.getCustomLibrary()->getMACROMolTemplateLib().clearTemplateLib();
  //   CHECK(mol.hasCustomLibrary() == true);
  //   CHECK(mol.getCustomLibrary()->getMACROMolTemplateLib().getNumTemplates() == 0);

  //   // Now uses global library again
  //   mol.addGlobalLibrary();
  //   CHECK(mol.usingGlobalLibrary());
  // }

  // SECTION("SharedLibraryBetweenMolecules") {
  //   // Multiple MonomerMols can share the same custom library
  //   auto externalLib = std::make_shared<MonomerLibrary>();
  //   externalLib->addMonomerFromSmiles("CC", "Q1", "PEPTIDE", "QQ1");

  //   MonomerMol mol1(externalLib.get());
  //   MonomerMol mol2(externalLib.get());

  //   CHECK(mol1.hasLocalTemplates() == false);
  //   CHECK(mol2.hasLocalTemplates() == false);

  //   // Adding to the shared library affects both
  //   externalLib->addMonomerFromSmiles("CCC", "Q2", "PEPTIDE", "QQ2");
  //   CHECK(mol1.getCustomLibrary()->contains("Q2", "PEPTIDE") ==  true);
  //   CHECK(mol2.getCustomLibrary()->contains("Q2", "PEPTIDE") ==  true);


  //   MonomerMol mol3(externalLib.get());
  //   CHECK(mol3.hasLocalTemplates() == false);
  //   CHECK(mol3.getCustomLibrary()->contains("Q2", "PEPTIDE") ==  true);

  // }

  SECTION("CopyPreservesLibrary") {
    auto customLib = std::make_unique<MACROMolTemplateLib>();
    addMonomerFromSmiles(*customLib, "CC", "W1", "PEPTIDE");

    MACROMol mol1;
    mol1.getTemplateLibrary()->copyTemplateLib(*customLib.get());
    addMonomer(mol1, "W1", 1, "PEPTIDE", "PEPTIDE1");
    
    // Copy constructor preserves the library
    MACROMol mol2(mol1);

    addMonomerFromSmiles(*customLib, "CCCC", "W2", "PEPTIDE");

    CHECK(mol2.getTemplateLibrary()->size() == mol1.getTemplateLibrary()->size());
    CHECK(mol2.getTemplateLibrary()->contains("PEPTIDE", "W1")== true);
    CHECK(mol2.getTemplateLibrary()->contains("PEPTIDE", "W2")== false);  // added after copied to mol
    
    mol2.getTemplateLibrary()->copyTemplateLib(*(customLib.get()));
    CHECK(mol2.getTemplateLibrary()->contains("PEPTIDE", "W2")== true);
  }


  SECTION("GetMonomer") {
    auto lib = GlobalMonomerLibrary::getGlobalLibrary();

    // Get a parsed molecule for alanine
    auto alaMol = lib->find("PEPTIDE", "A");
    auto test = alaMol != nullptr;
    CHECK(test);
  }

  SECTION("AddMonomerWithMol") {
    auto customLib = std::make_unique<MACROMolTemplateLib>();

    // Pre-parse a molecule
    auto mol = std::unique_ptr<RWMol>(SmilesToMol("CC(N)C(=O)O"));
    auto origAtomCount = mol->getNumAtoms();

    // Add with pre-parsed mol (no original data needed)
    addMonomer(*customLib.get(), mol, "CC(N)C(=O)O",  "TST" ,"PEPTIDE");


    // getMol returns the same pre-parsed molecule
    auto retrieved = customLib->find("PEPTIDE", "TST");
    CHECK(retrieved->getNumAtoms() == origAtomCount);
  }

  SECTION("AddMonomerFromSDF") {
    auto customLib = std::make_unique<MACROMolTemplateLib>();

    // Simple alanine-like structure in SDF format
    std::string sdfData = R"(
     RDKit          3D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  2  0  0  0  0  0  0
    2.2500    1.2990    0.0000 N   0  0  0  0  0  1  0  0  0  0  0  0
    2.2500   -1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  2  4  2  0
M  END
)";

    addMonomerFromSDF(*customLib.get(), sdfData, "SDF1", "PEPTIDE", "SD1");

    CHECK(customLib->contains("PEPTIDE", "SDF1"));

    // getMol should parse the SDF and return a molecule
    auto mol = customLib->find("PEPTIDE", "SDF1");
    REQUIRE(mol != nullptr);
    CHECK(mol->getNumAtoms() == 4);
    CHECK(mol->getNumBonds() == 3);
  }

  // SECTION("GlobalLibraryConfiguration") {
  //   // Save original state
  //   bool originalState = MonomerLibrary::isUsingGlobalLibrary();

  //   // Can toggle global library mode
  //   MonomerLibrary::useGlobalLibrary(false);
  //   CHECK(MonomerLibrary::isUsingGlobalLibrary() == false);

  //   MonomerLibrary::useGlobalLibrary(true);
  //   CHECK(MonomerLibrary::isUsingGlobalLibrary() == true);

  //   // Restore original state
  //   MonomerLibrary::useGlobalLibrary(originalState);
  // }

  SECTION("EmptyLibrary") {
    // Create an empty library (default behavior)
    auto emptyLib = std::make_unique<MACROMolTemplateLib>();

    // Should not have any built-in monomers
    CHECK(emptyLib->contains("A", "PEPTIDE") == false);
    CHECK(emptyLib->contains("G", "PEPTIDE") == false);
    CHECK(emptyLib->contains("C", "PEPTIDE") == false);

    // Can still add custom monomers
    addMonomerFromSmiles(*emptyLib.get(),
        "CC[C@H](N[H:1])C(=O)[OH:2]",
        "CUSTOM",
        "PEPTIDE",
        "CUS"
    );
    CHECK(emptyLib->contains("PEPTIDE", "CUSTOM"));

    // Use with MACROMol
    MACROMol mol;
    mol.getTemplateLibrary()->copyTemplateLib(*emptyLib.get());
    addMonomer(mol, "CUSTOM", 1, "PEPTIDE", "PEPTIDE1");
    CHECK(mol.getNumAtoms() == 1);

    // Convert to atomistic works with custom monomer
    auto atomistic = toAtomistic(mol);
    CHECK(atomistic->getNumAtoms() > 0);
  }

  SECTION("LibraryWithBuiltins") {
    // Explicitly create library with built-ins
    auto libWithBuiltins(GlobalMonomerLibrary::getGlobalLibrary());

    // Should have built-in monomers
    CHECK(libWithBuiltins->contains("PEPTIDE", "A"));
    CHECK(libWithBuiltins->contains("PEPTIDE","G"));
  }
}
