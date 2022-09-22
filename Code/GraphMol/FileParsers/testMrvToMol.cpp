//
//  Copyright (C) 2002-2021 Greg Landrum and other RDKit contributors
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/test.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Canon.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/MolPickler.h>
#include "FileParsers.h"
#include "SequenceParsers.h"
#include "SequenceWriters.h"
#include "MolFileStereochem.h"
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/FileParsers/ProximityBonds.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/LocaleSwitcher.h>
#include <clocale>
#include <cstdlib>


#include <GraphMol/QueryOps.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>

#include <string>
#include <fstream>
#include <boost/lexical_cast.hpp>

using namespace RDKit;

class MolTest
{
  public:
  std::string fileName;
  int atomCount;
  int bondCount;

  MolTest(std::string fileNameInit, int atomCountInit, int bondCountInit)
    : fileName(fileNameInit), atomCount(atomCountInit), bondCount(bondCountInit)
  {};

};

class RxnTest
{
  public:
  std::string fileName;
  int recactantCount;
  int agentCount;
  int productCount;

  RxnTest(std::string fileNameInit, int recactantCountInit, int agentCountInit, int productCountInit)
    : fileName(fileNameInit), recactantCount(recactantCountInit), agentCount(agentCountInit), productCount(productCountInit)
  {};

};

void testMrvToMol(const MolTest molTest)
{
  BOOST_LOG(rdInfoLog) << "testing simple mol parsing" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName =
      rdbase + "/Code/GraphMol/test_data/" + molTest.fileName;
  RWMol *m = (RWMol *)MrvFileParser(fName);
  // MolOps::sanitizeMol(*m);
  TEST_ASSERT(m != NULL);

  //TEST_ASSERT(m.atoms.count() == molTest.numberOfAtoms)
  //TEST_ASSERT(m.bonds.count() == molTest.numberOfbonds)
  
  delete m;
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testMrvToRxn(const RxnTest rxnTest)
{
  BOOST_LOG(rdInfoLog) << "testing simple mol parsing" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName =
      rdbase + "/Code/GraphMol/test_data/" + rxnTest.fileName;
  ChemicalReaction *r = (ChemicalReaction *)MrvFileParser(fName);
  // MolOps::sanitizeMol(*m);
  TEST_ASSERT(r != NULL);

  //TEST_ASSERT(r.reagents.count() == molTest.numberOfAtoms)
  //TEST_ASSERT(r.agents.count() == molTest.numberOfAtoms)
  //TEST_ASSERT(m.products.count() == molTest.numberOfbonds)
  
  delete r;
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void RunTests()
{
  
  // first the molecule tests

  std::list<MolTest> molFileNames
  {
     MolTest("ketback01.mrv",11,11)
    ,MolTest("ketback02.mrv",9,9)
    ,MolTest("ketback07.mrv",12,11)
    ,MolTest("ketback10.mrv",10,10)
    ,MolTest("marvin06.mrv",11,11)
    ,MolTest("ketback12.mrv",10,10)

  };

  for (std::list<MolTest>::const_iterator it = molFileNames.begin() ; it != molFileNames.end(); ++it)
  {
    printf("Test\n\n %s\n\n", it->fileName.c_str());
    testMrvToMol(*it);
  }


// now the reactions

 std::list<RxnTest> rxnFileNames
 {
     RxnTest("ketback03.mrv",1,1,1)
    ,RxnTest("ketback04.mrv",2,1,2)
    ,RxnTest("ketback08.mrv",2,4,2)
    ,RxnTest("ketback09.mrv",2,3,2)
    ,RxnTest("ketback11.mrv",2,0,1)
    ,RxnTest("marvin05.mrv",2,1,1)
  };

  for (std::list<RxnTest>::const_iterator it = rxnFileNames.begin() ; it != rxnFileNames.end(); ++it)
  {
    printf("Test\n\n %s\n\n", it->fileName.c_str());
    testMrvToRxn(*it);
  }

};

int main(int argc, char *argv[])
{
  (void)argc;
  (void)argv;
  //  std::locale::global(std::locale("de_DE.UTF-8"));

  RDLog::InitLogs();
  BOOST_LOG(rdInfoLog) << " ---- Running with POSIX locale ----- " << std::endl;
  RunTests(); // run with C locale

  return 0;
}
