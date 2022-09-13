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

#include <string>
#include <fstream>
#include <boost/lexical_cast.hpp>

using namespace RDKit;

void testMrvToMol()
{
  BOOST_LOG(rdInfoLog) << "testing atom query parsing" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName =
      rdbase + "/Code/GraphMol/test_data/ketback01.mrv";
  RWMol *m = MrvFileToMol(fName);
  // MolOps::sanitizeMol(*m);
  TEST_ASSERT(m);

  // TEST_ASSERT(m->getNumAtoms() == 6);
  // std::string smi = MolToSmiles(*m);
  // TEST_ASSERT(smi == "C1=CC=CC=C1");

  // m->updatePropertyCache();
  // smi = MolToSmarts(*m);
  // TEST_ASSERT(smi == "[#6]1=[#6]-[#6]=[#6]-[#6]=[#6,#7,#15]-1");

  // smi = "C1=CC=CC=C1";
  // RWMol *m2 = SmilesToMol(smi, false, false);
  // MatchVectType mv;
  // TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  // TEST_ASSERT(mv.size() == 6);
  // // sanitize it, which will aromatize the bonds... we will not match:
  // MolOps::sanitizeMol(*m2);
  // TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  // TEST_ASSERT(mv.size() == 0);

  // delete m2;
  // smi = "N1=CC=CC=C1";
  // m2 = SmilesToMol(smi, false, false);
  // TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  // TEST_ASSERT(mv.size() == 6);
  // delete m2;
  // smi = "S1=CC=CC=C1";
  // m2 = SmilesToMol(smi, false, false);
  // TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  // TEST_ASSERT(mv.size() == 0);
  // delete m2;
  // smi = "P1=CC=CC=C1";
  // m2 = SmilesToMol(smi, false, false);
  // TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  // TEST_ASSERT(mv.size() == 6);
  // delete m2;

  // delete m;
  // fName = rdbase + "/Code/GraphMol/FileParsers/test_data/not-list-query.mol";
  // m = MolFileToMol(fName);
  // TEST_ASSERT(m);

  // smi = "CC(=C)C";
  // m2 = SmilesToMol(smi, false, false);
  // TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  // TEST_ASSERT(mv.size() == 4);
  // delete m2;

  // smi = "CC(=O)C";
  // m2 = SmilesToMol(smi, false, false);
  // TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  // delete m2;

  // smi = "CC(=N)C";
  // m2 = SmilesToMol(smi, false, false);
  // TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  // delete m2;

  // smi = "CC(=O)C(=C)C";
  // m2 = SmilesToMol(smi, false, false);
  // TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  // TEST_ASSERT(mv.size() == 4);
  // delete m2;

  // smi = "C(=C)C";
  // m2 = SmilesToMol(smi, false, false);
  // TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  // delete m2;

  // // make sure new-style atom lists override old-style atom lists:
  // delete m;
  // fName = rdbase +
  //         "/Code/GraphMol/FileParsers/test_data/conflicting-list-query.mol";
  // m = MolFileToMol(fName);
  // TEST_ASSERT(m);

  // smi = "CC(=C)C";
  // m2 = SmilesToMol(smi, false, false);
  // TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  // TEST_ASSERT(mv.size() == 4);
  // delete m2;

  // smi = "CC(=O)C";
  // m2 = SmilesToMol(smi, false, false);
  // TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  // delete m2;

  // // longer list queries, this was issue 2413431:
  // delete m;
  // fName = rdbase + "/Code/GraphMol/FileParsers/test_data/list-query-long.mol";
  // m = MolFileToMol(fName);
  // TEST_ASSERT(m);
  // TEST_ASSERT(m->getAtomWithIdx(14)->hasQuery());

  // smi = "C1COC2=CC3=CC4=C(C=CC=C4)C=C3C=C2C1";
  // m2 = SmilesToMol(smi);
  // TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  // delete m2;
  // smi = "C1C[Se]C2=CC3=CC4=C(C=CC=C4)C=C3C=C2C1";
  // m2 = SmilesToMol(smi);
  // TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  // delete m2;
  // smi = "C1C[Te]C2=CC3=CC4=C(C=CC=C4)C=C3C=C2C1";
  // m2 = SmilesToMol(smi);
  // TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  // delete m2;
  // smi = "C1C[As]C2=CC3=CC4=C(C=CC=C4)C=C3C=C2C1";
  // m2 = SmilesToMol(smi);
  // TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  // delete m2;

  delete m;
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void RunTests()
{
#if 1
  testMrvToMol();
#endif
}

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
