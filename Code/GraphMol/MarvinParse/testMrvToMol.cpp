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


enum LoadAs
{
  LoadAsMolOrRxn
  ,LoadAsMol
  ,LoadAsRxn
};


class MolOrRxnTest
{
  public:



  std::string fileName;
  bool expectedResult;
  LoadAs loadAs;

  MolOrRxnTest(std::string fileNameInit, bool expectedResultInit, LoadAs loadAsInit)
    : fileName(fileNameInit)
    , expectedResult(expectedResultInit)
    , loadAs(loadAsInit)
  {};

  virtual bool isRxnTest() const =0;
};

class MolTest : public MolOrRxnTest
{
  public:
  unsigned int atomCount;
  unsigned int bondCount;

  MolTest(std::string fileNameInit, bool expectedResultInit, LoadAs loadAsInit, int atomCountInit, int bondCountInit)
    : MolOrRxnTest(fileNameInit, expectedResultInit, loadAsInit)
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

  RxnTest(std::string fileNameInit, bool expectedResultInit, LoadAs loadAsInit, int reactantCountInit, int agentCountInit, int productCountInit, int warnInit, int errorInit)
    : MolOrRxnTest(fileNameInit, expectedResultInit, loadAsInit)
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





void *GetMolOrReaction(const MolOrRxnTest *molOrRxnTest, bool &isReaction)
{
  std::string rdbase = getenv("RDBASE");
  std::string fName =
      rdbase + "/Code/GraphMol/MarvinParse/test_data/" + molOrRxnTest->fileName;
  std::string ofName = fName + ".sdf";


  switch(molOrRxnTest->loadAs)
  {
    case LoadAsMolOrRxn:
        return MrvFileParser(fName, isReaction);

    case LoadAsMol:
        isReaction = false;
        return (void *)MrvMolFileParser(fName);  

    case LoadAsRxn:
        isReaction = true;
        return (void *)MrvRxnFileParser(fName);

  }
  
  throw   BadFileException("Should never happen");
  
}

void testMarvin(const MolOrRxnTest *molOrRxnTest)
{
  BOOST_LOG(rdInfoLog) << "testing marvin parsing" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName =
      rdbase + "/Code/GraphMol/MarvinParse/test_data/" + molOrRxnTest->fileName;
  std::string ofName = fName + ".sdf";

  class LocalVars  // protext agains mem leak on error
  {
    public:
    void *molOrRxn = NULL;
    bool isReaction=false;

    LocalVars(){};

    ~LocalVars()
    {  

      //printf ("deleting molOrRxn\n");
      if (isReaction)
        delete (ChemicalReaction *)molOrRxn;
      else
        delete (RWMol *) molOrRxn;
    }
  } localVars;

  try
  {  

    localVars.molOrRxn = GetMolOrReaction(molOrRxnTest, localVars.isReaction);
    if (localVars.isReaction != molOrRxnTest->isRxnTest())
    {
      //printf("Wrong type of MRV file\n");
      TEST_ASSERT(molOrRxnTest->expectedResult == false);
      //printf("Expected failure!\n");
      return;
    }

 
    if (localVars.isReaction)
    {
      // reaction test

      auto rxn = (ChemicalReaction *)localVars.molOrRxn;
      std::string ofName = fName + ".rxn";
      std::string outMolStr = ChemicalReactionToRxnBlock(*rxn, false, true);
      {
        std::ofstream  out;
        out.open(ofName);
        out << outMolStr << "\n";
      }  
      
      ofName =  fName + ".OUT.mrv";
      outMolStr = ChemicalReactionToMrvBlock(*rxn);
      {
        std::ofstream  out;
        out.open(ofName);
        out << outMolStr << "\n";
      }      

      unsigned int nWarn=0, nError=0;
        
      TEST_ASSERT(rxn!= NULL);

      //printf("rxnTest->reactantCount): %d rxnTest->productCount: %d rxnTest->agentCount: %d\n", rxnTest->reactantCount, rxnTest->productCount, rxnTest->agentCount);
      //printf("rxn->getNumReactantTemplates(): %d >getNumProductTemplates(): %d rxn->getNumAgentTemplates(): %d\n",
      //    rxn->getNumReactantTemplates(), rxn->getNumProductTemplates(), rxn->getNumAgentTemplates());

      auto rxnTest = (RxnTest *)molOrRxnTest;
      TEST_ASSERT(rxn->getNumReactantTemplates() == rxnTest->reactantCount);
      TEST_ASSERT(rxn->getNumProductTemplates() == rxnTest->productCount);
      TEST_ASSERT(rxn->getNumAgentTemplates() == rxnTest->agentCount);
      rxn->initReactantMatchers(true);


      if (rxn->getNumReactantTemplates() > 0 && rxn->getNumProductTemplates() > 0)
      {
        TEST_ASSERT(rxn->validate(nWarn, nError, true));
      }
      else
      {
        nWarn=0;
        nError=0;
      }

      //printf("nWarn: %d nError: %d\n", nWarn,nError);
      TEST_ASSERT(nWarn == rxnTest->warnings);
      TEST_ASSERT(nError == rxnTest->errors);

      BOOST_LOG(rdInfoLog) << "done" << std::endl;     
    }
    else
    {
      // mol test

      auto mol = (RWMol *)localVars.molOrRxn;

      std::string ofName = fName + ".sdf";
      std::string outMolStr =  MolToMolBlock(*mol, true, 0, true, true);
      {                    
        std::ofstream  out;
        out.open(ofName);
        out << outMolStr << "\n";
      }      

      ofName =  fName + ".OUT.mrv";
      outMolStr = MolToMrvBlock(*mol,true, -1, true);
      {
        std::ofstream  out;
        out.open(ofName);
        out << outMolStr << "\n";
      }      

      // MolOps::sanitizeMol(*m);
      TEST_ASSERT(mol != NULL);

      auto molTest = (MolTest *)molOrRxnTest;

      // if (mol->getNumAtoms() != molTest->atomCount)
      //   printf("mol->getNumAtoms(): %d    molTest->atomCount: %d\n" ,mol->getNumAtoms(), molTest->atomCount);
      // if (mol->getNumBonds() != molTest->bondCount)
      //   printf("mol->getNumBonds(): %d    molTest->bondCount: %d\n" ,mol->getNumBonds(), molTest->bondCount);
      TEST_ASSERT(mol->getNumAtoms() == molTest->atomCount)
      TEST_ASSERT(mol->getNumBonds() == molTest->bondCount)
      
      BOOST_LOG(rdInfoLog) << "done" << std::endl;    
    }
  }
  catch(const std::exception& e)
  {


    if(molOrRxnTest->expectedResult != false)
        throw;
    return;
    
  }

  TEST_ASSERT(molOrRxnTest->expectedResult == true);
   
  return; 

}

void RunTests()
{
  
  // first the molecule tests

  std::list<MolTest> molFileNames
  {
     MolTest("ketback01.mrv", true, LoadAsMolOrRxn, 11,11)
    ,MolTest("ketback01.mrv", true, LoadAsMol, 11,11)
    ,MolTest("ketback01.mrv", false, LoadAsRxn, 11,11)   // should fail
    ,MolTest("ketback02.mrv",true,LoadAsMolOrRxn, 9,9)
    ,MolTest("ketback07.mrv",true,LoadAsMolOrRxn, 12,11)
    ,MolTest("ketback10.mrv",true,LoadAsMolOrRxn, 10,10)
    ,MolTest("marvin06.mrv",true,LoadAsMolOrRxn, 11,11)
    ,MolTest("ketback12.mrv",true,LoadAsMolOrRxn, 31,33)
    ,MolTest("EmptyMol.mrv",true,LoadAsMolOrRxn, 0,0)
    ,MolTest("Sparse.mrv",true,LoadAsMolOrRxn, 0,0)
    ,MolTest("Sparse2.mrv",true,LoadAsMolOrRxn, 0,0)
    ,MolTest("Sparse3.mrv",true,LoadAsMolOrRxn, 0,0)
    ,MolTest("MarvinNoCoords.mrv",true,LoadAsMolOrRxn, 6,6)
    ,MolTest("aspirin.mrv",true, LoadAsMolOrRxn, 13,13)  
    ,MolTest("MarvinStereoGroupsZeros.mrv",true, LoadAsMolOrRxn, 8,8)  
    ,MolTest("MarvinMissingX2.mrv",true,LoadAsMolOrRxn, 12,11)  // should fail - missing atom X2 coord
    ,MolTest("MarvinMissingY2.mrv",true,LoadAsMolOrRxn, 12,11)  // should fail - missing atom Y2 coord
    ,MolTest("ketback03.mrv",false,LoadAsMolOrRxn, 31,33)  // should fail - this is a reaction
    ,MolTest("MarvinBadMissingMolID.mrv",false,LoadAsMolOrRxn, 12,11)  // should fail - no molId
    ,MolTest("MarvinBadMissingAtomID.mrv",false,LoadAsMolOrRxn, 12,11)  // should fail - missing atom Id
    ,MolTest("MarvinBadX2.mrv",false,LoadAsMolOrRxn, 12,11)  // should fail - 
    ,MolTest("MarvinBadY2.mrv",false,LoadAsMolOrRxn, 12,11)  // should fail - 
    ,MolTest("MarvinBadElementType.mrv",false, LoadAsMolOrRxn, 12,11)  // should fail - 
    ,MolTest("MarvinBadMissingBondID.mrv",false, LoadAsMolOrRxn, 12,11)  // should fail - 
    ,MolTest("MarvinBadMissingBondAtomRefs",false, LoadAsMolOrRxn, 12,11)  // should fail - 
    ,MolTest("MarvinBadMissingBondOrder.mrv",false, LoadAsMolOrRxn, 12,11)  // should fail - 
    ,MolTest("MarvinBadMissingSruMolID.mrv",false, LoadAsMolOrRxn, 12,11)  // should fail - 
    ,MolTest("MarvinBadMissingSruID.mrv",false, LoadAsMolOrRxn, 12,11)  // should fail - 
    ,MolTest("MarvinBadMissingSruRole.mrv",false, LoadAsMolOrRxn, 12,11)  // should fail - 
    ,MolTest("MarvinBadMissingSruAtomRef.mrv",false, LoadAsMolOrRxn, 12,11)  // should fail - 
    ,MolTest("MarvinBadMissingSruTitle.mrv",false, LoadAsMolOrRxn, 12,11)  // should fail - 
    ,MolTest("MarvinBadSruAtomRef.mrv",false, LoadAsMolOrRxn, 12,11)  // should fail - 
    ,MolTest("MarvinBadSruID.mrv",false, LoadAsMolOrRxn, 12,11)  // should fail - 
    ,MolTest("MarvinBadSruRole.mrv",false, LoadAsMolOrRxn, 12,11)  // should fail - 
    ,MolTest("MarvinBadSruAtomRef.mrv",false, LoadAsMolOrRxn, 12,11)  // should fail - 
    ,MolTest("MarvinBadSruAtomRef.mrv",false, LoadAsMolOrRxn, 12,11)  // should fail - 
    ,MolTest("MarvinBadSruConnect.mrv",false, LoadAsMolOrRxn, 12,11)  // should fail - 
    ,MolTest("MarvinBadSupAttachAtom.mrv",false, LoadAsMolOrRxn, 9,9)  // should fail -
    ,MolTest("MarvinBadSupAttachBond.mrv",false, LoadAsMolOrRxn, 9,9)  // should fail -
    ,MolTest("MarvinBadSupAttachOrder.mrv",false, LoadAsMolOrRxn, 9,9)  // should fail -
    ,MolTest("MarvinBadSupAttachAtom.mrv",false, LoadAsMolOrRxn, 9,9)  // should fail -
    ,MolTest("MarvinBadSupAttachAtom.mrv",false, LoadAsMolOrRxn, 9,9)  // should fail -
    ,MolTest("MarvinBadSupAttachAtom.mrv",false, LoadAsMolOrRxn, 9,9)  // should fail -
    ,MolTest("MarvinBadSupMissingAttachBond.mrv",false, LoadAsMolOrRxn, 9,9)  // should fail -
    ,MolTest("MarvinBadSupMissingAttachOrder.mrv",false, LoadAsMolOrRxn, 9,9)  // should fail -
  };

  for (std::list<MolTest>::const_iterator it = molFileNames.begin() ; it != molFileNames.end(); ++it)
  {
    BOOST_LOG(rdInfoLog) << "Test: " << it->fileName << std::endl;

    printf("Test\n\n %s\n\n", it->fileName.c_str());
    testMarvin(&*it);
  }


// now the reactions

 std::list<RxnTest> rxnFileNames
 {
     RxnTest("ketback03.mrv",true, LoadAsMolOrRxn,1,1,1,2,0)
    ,RxnTest("ketback03.mrv",true, LoadAsRxn,1,1,1,2,0)
    ,RxnTest("ketback03.mrv",false, LoadAsMol,1,1,1,2,0)   // should fail
    ,RxnTest("ketback04.mrv",true, LoadAsMolOrRxn,2,1,2,4,0)
    ,RxnTest("ketback08.mrv",true, LoadAsMolOrRxn,2,3,2,4,0)
    ,RxnTest("ketback09.mrv",true, LoadAsMolOrRxn,2,3,2,4,0)
    ,RxnTest("ketback11.mrv",true, LoadAsMolOrRxn,2,0,1,0,0)
    ,RxnTest("marvin05.mrv",true, LoadAsMolOrRxn,2,1,1,3,0)
    ,RxnTest("EmptyRxn.mrv",true, LoadAsMolOrRxn,0,0,0,0,0)
    ,RxnTest("ketback01.mrv",false,LoadAsMolOrRxn,2,1,1,3,0)  // should fail - this is a mol file
    ,RxnTest("aspirineSynthesisWithAttributes.mrv",true,LoadAsMolOrRxn,2,0,1,3,0)  // should fail - this is a mol file
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
