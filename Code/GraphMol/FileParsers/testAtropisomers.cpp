//
//  Copyright (C) 2002-2021 Collaboartive Drug Discovery and other RDKit
//  contributors
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/SequenceParsers.h>
#include <GraphMol/FileParsers/SequenceWriters.h>
#include <GraphMol/MolFileStereochem.h>
#include <GraphMol/Depictor/RDDepictor.h>

#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Atropisomers.h>

#include <string>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <filesystem>
using namespace RDKit;

enum LoadAs { LoadAsMolOrRxn, LoadAsMol, LoadAsRxn };

class MolOrRxnTest {
 public:
  std::string fileName;
  bool expectedResult;
  LoadAs loadAs;

  MolOrRxnTest(std::string fileNameInit, bool expectedResultInit,
               LoadAs loadAsInit)
      : fileName(fileNameInit),
        expectedResult(expectedResultInit),
        loadAs(loadAsInit){};

  virtual bool isRxnTest() const = 0;
};

class MolTest : public MolOrRxnTest {
 public:
  unsigned int atomCount;
  unsigned int bondCount;

  MolTest(std::string fileNameInit, bool expectedResultInit, LoadAs loadAsInit,
          int atomCountInit, int bondCountInit)
      : MolOrRxnTest(fileNameInit, expectedResultInit, loadAsInit),
        atomCount(atomCountInit),
        bondCount(bondCountInit){};

  bool isRxnTest() const override { return false; }
};

void testLookForAtropisomersInSDdfFiles(std::string fileName,
                                        unsigned int expectedHits,
                                        unsigned int expectedMisses) {
  BOOST_LOG(rdInfoLog) << "Looking for atropisomers in " << fileName
                       << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName =
      rdbase + "/Code/GraphMol/FileParsers/test_data/atropisomers/" + fileName;

  std::ifstream in;
  in.open(fName);
  std::string line;
  unsigned int foundCount = 0;
  unsigned int notFoundCount = 0;
  while (!in.eof()) {
    std::string molBlock = "";
    while (std::getline(in, line)) {
      if (line.find("$$$$") != std::string::npos) {
        break;
      }

      molBlock += line + "\n";
    }

    if (molBlock.length() < 10) {
      continue;  // try for another;
    }

    std::unique_ptr<RWMol> mol(
        MolBlockToMol(molBlock, false, false, false, true));
    TEST_ASSERT(mol != nullptr);

    auto hasAtropisomers = doesMolHaveAtropisomers(*mol);

    if (hasAtropisomers) {
      BOOST_LOG(rdInfoLog) << "Found atropisomers in " << fileName << std::endl;
      foundCount++;
      printf("Atropisomers- %d hits   %d misses\r", foundCount, notFoundCount);
      std::flush(std::cout);
      std::ofstream out;
      out.open(fName + "_" + std::to_string(foundCount) + ".sdf");
      out << molBlock << std::endl;
    } else {
      notFoundCount++;
      if (notFoundCount % 100 == 0) {
        printf("Atropisomers- %d hits   %d misses\r", foundCount,
               notFoundCount);
        std::flush(std::cout);
      }
    }
  }
  printf("\nFinal results:\nFound atropisomers in %s - %d hits   %d misses\n",
         fileName.c_str(), foundCount, notFoundCount);

  TEST_ASSERT(foundCount == expectedHits);
  TEST_ASSERT(notFoundCount == expectedMisses);
}

void testMolFiles(const MolTest *molFileTest) {
  BOOST_LOG(rdInfoLog) << "testing mol files with atropisomers" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName = rdbase +
                      "/Code/GraphMol/FileParsers/test_data/atropisomers/" +
                      molFileTest->fileName;

  try {
    std::unique_ptr<RWMol> mol(MolFileToMol(fName, false, false, false, true));

    TEST_ASSERT(mol != nullptr);
    TEST_ASSERT(mol->getNumAtoms() == molFileTest->atomCount)
    TEST_ASSERT(mol->getNumBonds() == molFileTest->bondCount)

    {
      std::string expectedFileName = fName + ".expected.sdf";
      std::string outMolStr = "";
      MolOps::Kekulize(*mol);
      reapplyMolBlockWedging(*mol);
      outMolStr = MolToMolBlock(*mol, true, 0, true, true);

      // code to create the expected files for new or changed tests

      // {
      //   std::ofstream out;
      //   out.open(fName + ".NEW.sdf");
      //   out << outMolStr;
      // }

      std::stringstream expectedMolStr;
      std::ifstream in;
      in.open(expectedFileName);
      expectedMolStr << in.rdbuf();
      std::string expectedStr = expectedMolStr.str();

      TEST_ASSERT(expectedStr == outMolStr);
    }

    // 2nd pass without reapplying the mol block wedging -
    //     the itropisomers will be marked automatically

    mol = std::unique_ptr<RWMol>(MolFileToMol(fName, true, false, false, true));

    {
      std::string expectedFileName = fName + ".expected2.sdf";
      std::string outMolStr = "";
      outMolStr = MolToMolBlock(*mol, true, 0, true, true);

      // code to create the expected files for new or changed tests

      // {
      //   std::ofstream out;
      //   out.open(fName + ".NEW2.sdf");
      //   out << outMolStr;
      // }

      std::stringstream expectedMolStr;
      std::ifstream in;
      in.open(expectedFileName);
      expectedMolStr << in.rdbuf();
      std::string expectedStr = expectedMolStr.str();

      TEST_ASSERT(expectedStr == outMolStr);
    }

    // CXSMILES

    mol =
        std::unique_ptr<RWMol>(MolFileToMol(fName, false, false, false, true));

    {
      std::string expectedFileName = fName + ".expected.cxsmi";
      SmilesWriteParams ps;
      ps.canonical = false;
      ps.doKekule = true;

      unsigned int flags = SmilesWrite::CXSmilesFields::CX_COORDS |
                           SmilesWrite::CXSmilesFields::CX_MOLFILE_VALUES |
                           SmilesWrite::CXSmilesFields::CX_ATOM_PROPS |
                           SmilesWrite::CXSmilesFields::CX_BOND_CFG |
                           SmilesWrite::CXSmilesFields::CX_ENHANCEDSTEREO
          //| SmilesWrite::CXSmilesFields::CX_ALL
          ;

      std::string smilesOut =
          MolToCXSmiles(*mol, ps, flags, RestoreBondDirOptionTrue);

      // code to create the expected files for new or changed tests

      // {
      //   std::ofstream out;
      //   out.open(fName + ".NEW.cxsmi");
      //   out << smilesOut;
      // }

      std::stringstream expectedMolStr;
      std::ifstream in;
      in.open(expectedFileName);
      expectedMolStr << in.rdbuf();
      std::string expectedStr = expectedMolStr.str();

      TEST_ASSERT(expectedStr == smilesOut);
    }

    BOOST_LOG(rdInfoLog) << "done" << std::endl;
  } catch (const std::exception &e) {
    if (molFileTest->expectedResult != false) {
      throw;
    }
    return;
  }

  TEST_ASSERT(molFileTest->expectedResult == true);

  return;
}

void RunTests() {
  // the molecule tests

  std::list<MolTest> sdfTests{
      MolTest("AtropTest.sdf", true, LoadAsMolOrRxn, 38, 41),
      MolTest("AtropManyChiralsEnhanced.sdf", true, LoadAsMolOrRxn, 20, 20),
      MolTest("AtropManyChiralsEnhanced2.sdf", true, LoadAsMolOrRxn, 20, 20),
      MolTest("AtropManyChirals.sdf", true, LoadAsMolOrRxn, 20, 20),
      MolTest("BMS-986142.sdf", true, LoadAsMolOrRxn, 42, 47),
      MolTest("BMS-986142_3d_chiral.sdf", true, LoadAsMolOrRxn, 72, 77),
      MolTest("BMS-986142_3d.sdf", true, LoadAsMolOrRxn, 72, 77),
      MolTest("BMS-986142_atrop1.sdf", true, LoadAsMolOrRxn, 42, 47),
      MolTest("BMS-986142_atrop2.sdf", true, LoadAsMolOrRxn, 42, 47),
      MolTest("BMS-986142_atrop3.sdf", true, LoadAsMolOrRxn, 42, 47),
      MolTest("BMS-986142_atrop4.sdf", true, LoadAsMolOrRxn, 42, 47),
      MolTest("BMS-986142_atrop5.sdf", true, LoadAsMolOrRxn, 42, 47),
      MolTest("BMS-986142_atrop6.sdf", true, LoadAsMolOrRxn, 42, 47),
      MolTest("BMS-986142_atrop7.sdf", true, LoadAsMolOrRxn, 42, 47),
      MolTest("BMS-986142_atrop8.sdf", true, LoadAsMolOrRxn, 42, 47),
      MolTest("BMS-986142_atropBad2.sdf", true, LoadAsMolOrRxn, 42, 47),
      MolTest("JDQ443.sdf", true, LoadAsMolOrRxn, 38, 44),
      MolTest("JDQ443_3d.sdf", true, LoadAsMolOrRxn, 66, 72),
      MolTest("JDQ443_atrop1.sdf", true, LoadAsMolOrRxn, 38, 44),
      MolTest("JDQ443_atrop2.sdf", true, LoadAsMolOrRxn, 38, 44),
      MolTest("JDQ443_atrop3.sdf", true, LoadAsMolOrRxn, 38, 44),
      MolTest("JDQ443_atropBad1.sdf", true, LoadAsMolOrRxn, 38, 44),
      MolTest("RP-6306.sdf", true, LoadAsMolOrRxn, 24, 26),
      MolTest("RP-6306_atrop1.sdf", true, LoadAsMolOrRxn, 24, 26),
      MolTest("RP-6306_atrop2.sdf", true, LoadAsMolOrRxn, 24, 26),
      MolTest("RP-6306_atrop3.sdf", true, LoadAsMolOrRxn, 24, 26),
      MolTest("RP-6306_atrop4.sdf", true, LoadAsMolOrRxn, 24, 26),
      MolTest("RP-6306_atrop5.sdf", true, LoadAsMolOrRxn, 24, 26),
      MolTest("RP-6306_atropBad1.sdf", true, LoadAsMolOrRxn, 24, 26),
      MolTest("RP-6306_atropBad2.sdf", true, LoadAsMolOrRxn, 24, 26),
      // note the rp-6306_3d.sdf is backwards from the 2D versions
      // the 2D version were based on images from drug hunter
      // the 3D version came from PUBCHEM
      MolTest("RP-6306_3d.sdf", true, LoadAsMolOrRxn, 44, 46),
      MolTest("Sotorasib.sdf", true, LoadAsMolOrRxn, 41, 45),
      MolTest("Sotorasib_atrop1.sdf", true, LoadAsMolOrRxn, 41, 45),
      MolTest("Sotorasib_atrop2.sdf", true, LoadAsMolOrRxn, 41, 45),
      MolTest("Sotorasib_atrop3.sdf", true, LoadAsMolOrRxn, 41, 45),
      MolTest("Sotorasib_atrop4.sdf", true, LoadAsMolOrRxn, 41, 45),
      MolTest("Sotorasib_atrop5.sdf", true, LoadAsMolOrRxn, 41, 45),
      MolTest("Sotorasib_atropBad1.sdf", true, LoadAsMolOrRxn, 41, 45),
      MolTest("Sotorasib_atropBad2.sdf", true, LoadAsMolOrRxn, 41, 45),
      // note the sotorasib_3d.sdf is backwards from the 2D versions
      // the 2D version were based on images from drug hunter
      // the 3D version came from PUBCHEM
      MolTest("Sotorasib_3d.sdf", true, LoadAsMolOrRxn, 71, 75),
      MolTest("ZM374979.sdf", true, LoadAsMolOrRxn, 45, 49),
      MolTest("ZM374979_atrop1.sdf", true, LoadAsMolOrRxn, 45, 49),
      MolTest("ZM374979_atrop2.sdf", true, LoadAsMolOrRxn, 45, 49),
      MolTest("ZM374979_atrop3.sdf", true, LoadAsMolOrRxn, 45, 49),
      MolTest("ZM374979_atropBad1.sdf", true, LoadAsMolOrRxn, 45, 49),
      MolTest("Mrtx1719.sdf", true, LoadAsMolOrRxn, 33, 37),
      // note the Mrtx1719_3d.sdf is backwards from the 2D versions
      // the 2D version were based on images from drug hunter
      // the 3D version came from PUBCHEM
      MolTest("Mrtx1719_3d.sdf", true, LoadAsMolOrRxn, 51, 55),
      MolTest("Mrtx1719_atrop1.sdf", true, LoadAsMolOrRxn, 33, 37),
      MolTest("Mrtx1719_atrop2.sdf", true, LoadAsMolOrRxn, 33, 37),
      MolTest("Mrtx1719_atrop3.sdf", true, LoadAsMolOrRxn, 33, 37),
      MolTest("Mrtx1719_atropBad1.sdf", true, LoadAsMolOrRxn, 33, 37),
  };

  for (auto sdfTest : sdfTests) {
    BOOST_LOG(rdInfoLog) << "Test: " << sdfTest.fileName << std::endl;

    printf("Test\n\n %s\n\n", sdfTest.fileName.c_str());
    testMolFiles(&sdfTest);
  }
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;
  //  std::locale::global(std::locale("de_DE.UTF-8"));

  RDLog::InitLogs();
  BOOST_LOG(rdInfoLog) << " ---- Running with POSIX locale ----- " << std::endl;

  RunTests();  // run with C locale

  testLookForAtropisomersInSDdfFiles("TestMultInSDF.sdf", 1, 4);
  // testLookForAtropisomersInSDdfFiles("AllVendors.sdf", 76, 1523042);

  return 0;
}
