//
//  Copyright (C) 2002-2021 Greg Landrum and other RDKit contributors
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
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include "MarvinParser.h"


#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>

#include <string>
#include <fstream>
#include <boost/lexical_cast.hpp>

using namespace RDKit;

class MolOrRxnTest
{
  public:
  std::string fileName;
  bool expectedResult;

  MolOrRxnTest(std::string fileNameInit, bool expectedResultInit)
    : fileName(fileNameInit)
    , expectedResult(expectedResultInit)
  {};

  virtual bool isRxnTest() const =0;
};

class MolTest : public MolOrRxnTest
{
  public:
  unsigned int atomCount;
  unsigned int bondCount;

  MolTest(std::string fileNameInit, bool expectedResultInit, int atomCountInit, int bondCountInit)
    : MolOrRxnTest(fileNameInit, expectedResultInit)
    ,atomCount(atomCountInit)
    , bondCount(bondCountInit)
  {};

  bool isRxnTest() const
  {
    return false;
  }

};

class RxnTest : public MolOrRxnTest
{
  public:
  unsigned int reactantCount;
  unsigned int agentCount;
  unsigned int productCount;
  unsigned int warnings;
  unsigned int errors;

  RxnTest(std::string fileNameInit, bool expectedResultInit, int reactantCountInit, int agentCountInit, int productCountInit, int warnInit, int errorInit)
    : MolOrRxnTest(fileNameInit, expectedResultInit)
    , reactantCount(reactantCountInit)
    , agentCount(agentCountInit)
    , productCount(productCountInit)
    ,warnings(warnInit), errors(errorInit)
  {};

  bool isRxnTest() const
  {
    return true;
  }
};



void testMrvToMol(RWMol *mol,const MolTest *molTest)
{
  // MolOps::sanitizeMol(*m);
  TEST_ASSERT(mol != NULL);

  if (mol->getNumAtoms() != molTest->atomCount)
    printf("mol->getNumAtoms(): %d    molTest->atomCount: %d\n" ,mol->getNumAtoms(), molTest->atomCount);
  if (mol->getNumBonds() != molTest->bondCount)
    printf("mol->getNumBonds(): %d    molTest->bondCount: %d\n" ,mol->getNumBonds(), molTest->bondCount);
  TEST_ASSERT(mol->getNumAtoms() == molTest->atomCount)
  TEST_ASSERT(mol->getNumBonds() == molTest->bondCount)
  
  delete mol;
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testMrvToRxn(ChemicalReaction *rxn, const RxnTest *rxnTest)
{
  unsigned int nWarn, nError;
           
  TEST_ASSERT(rxn!= NULL);

  printf("rxnTest->reactantCount): %d rxnTest->productCount: %d rxnTest->agentCount: %d\n", rxnTest->reactantCount, rxnTest->productCount, rxnTest->agentCount);
  printf("rxn->getNumReactantTemplates(): %d >getNumProductTemplates(): %d rxn->getNumAgentTemplates(): %d\n",
       rxn->getNumReactantTemplates(), rxn->getNumProductTemplates(), rxn->getNumAgentTemplates());

  TEST_ASSERT(rxn->getNumReactantTemplates() == rxnTest->reactantCount);
  TEST_ASSERT(rxn->getNumProductTemplates() == rxnTest->productCount);
  TEST_ASSERT(rxn->getNumAgentTemplates() == rxnTest->agentCount);
  rxn->initReactantMatchers();
  TEST_ASSERT(rxn->validate(nWarn, nError, false));

  printf("nWarn: %d nError: %d\n", nWarn,nError);
  TEST_ASSERT(nWarn == rxnTest->warnings);
  TEST_ASSERT(nError == rxnTest->errors);

  
  delete rxn;
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testMarvin(const MolOrRxnTest *molOrRxnTest)
{
  BOOST_LOG(rdInfoLog) << "testing simple mol parsing" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName =
      rdbase + "/Code/GraphMol/MarvinParse/test_data/" + molOrRxnTest->fileName;
  std::string ofName = fName + ".sdf";

  bool isReaction;
  void *molOrRxn = (RWMol *)MrvFileParser(fName, isReaction);

  if (isReaction != molOrRxnTest->isRxnTest())
  {
    printf("Wrong type of MRV file\n");
    TEST_ASSERT(molOrRxnTest->expectedResult == false);
    return; 
  }

  try
  {
    if (isReaction)
    {
      // reaction test

      auto rxn = (ChemicalReaction *)molOrRxn;
      std::string ofName = fName + ".rxn";
      std::string outMolStr = ChemicalReactionToRxnBlock(*rxn, false, true);
      {
        std::ofstream  out;
        out.open(ofName);
        out << outMolStr << "\n";
      }      

      testMrvToRxn(rxn,(RxnTest *) molOrRxnTest);       
    }
    else
    {
      // mol test

      auto mol = (RWMol *)molOrRxn;
      std::string ofName = fName + ".sdf";
      std::string outMolStr =  MolToMolBlock(*mol, true, 0, true, true);
        {                    
          std::ofstream  out;
          out.open(ofName);
          out << outMolStr << "\n";
        }      

      testMrvToMol(mol,(MolTest *) molOrRxnTest);

    }
  }
  catch(const std::exception& e)
  {
    TEST_ASSERT(molOrRxnTest->expectedResult == false);
    return; 
  }
  
}

void RunTests()
{
  
  // first the molecule tests

  std::list<MolTest> molFileNames
  {
     MolTest("ketback01.mrv", true,11,11)
    ,MolTest("ketback02.mrv",true,9,9)
    ,MolTest("ketback07.mrv",true,12,11)
    ,MolTest("ketback10.mrv",true,10,10)
    ,MolTest("marvin06.mrv",true,11,11)
    ,MolTest("ketback12.mrv",true,31,33)
    ,MolTest("ketback03.mrv",false,31,33)  // should fail - thisis a reaction

  };

  for (std::list<MolTest>::const_iterator it = molFileNames.begin() ; it != molFileNames.end(); ++it)
  {
    printf("Test\n\n %s\n\n", it->fileName.c_str());
    testMarvin(&*it);
  }


// now the reactions

 std::list<RxnTest> rxnFileNames
 {
     RxnTest("ketback03.mrv",true,1,1,1,2,0)
    ,RxnTest("ketback04.mrv",true,2,1,2,4,0)
    ,RxnTest("ketback08.mrv",true,2,3,2,4,0)
    ,RxnTest("ketback09.mrv",true,2,3,2,4,0)
    ,RxnTest("ketback11.mrv",true,2,0,1,0,0)
    ,RxnTest("marvin05.mrv",true,2,1,1,3,0)
    ,RxnTest("ketback01.mrv",false,2,1,1,3,0)  // should fail - this is a reaction file
  };

  for (std::list<RxnTest>::const_iterator it = rxnFileNames.begin() ; it != rxnFileNames.end(); ++it)
  {
    printf("Test\n\n %s\n\n", it->fileName.c_str());
    testMarvin(&*it);
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
