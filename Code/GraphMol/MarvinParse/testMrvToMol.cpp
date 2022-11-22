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
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <string>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <filesystem>

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

class SmilesTest
{
  public:
  std::string name;
  std::string smiles;
  bool expectedResult;
  unsigned int atomCount;
  unsigned int bondCount;

  SmilesTest(std::string nameInit, std::string smilesInit, bool expectedResultInit, int atomCountInit, int bondCountInit)
    : name(nameInit)
    , smiles(smilesInit)
    , expectedResult(expectedResultInit)
    ,atomCount(atomCountInit)
    , bondCount(bondCountInit)
  {};

  bool isRxnTest() const
  {
    return false;
  }

};

void *GetMolOrReaction(const MolOrRxnTest *molOrRxnTest, bool &isReaction)
{
  std::string rdbase = getenv("RDBASE");
  std::string fName =
      rdbase + "/Code/GraphMol/MarvinParse/test_data/" + molOrRxnTest->fileName;

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

void testSmilesToMarvin(const SmilesTest *smilesTest)
{
  BOOST_LOG(rdInfoLog) << "testing smiles to marin " << std::endl;
  std::string rdbase = getenv("RDBASE");
  std::string fName = rdbase + "/Code/GraphMol/MarvinParse/test_data/" + smilesTest->name;

  class LocalVars  // protext against mem leak on error
  {
   public:
    RWMol *smilesMol;

    LocalVars(){};

    ~LocalVars() 
    {     
      delete smilesMol;
    }
  } localVars;

  try 
  {
    SmilesParserParams smilesParserParams;
    smilesParserParams.sanitize = true;

    localVars.smilesMol = SmilesToMol(smilesTest->smiles, smilesParserParams);
    std::string mrvBlock = MolToMrvBlock(*localVars.smilesMol, true, -1, true);

    delete localVars.smilesMol;
    localVars.smilesMol = NULL;

    localVars.smilesMol = MrvMolStringParser(mrvBlock);

    TEST_ASSERT(localVars.smilesMol->getNumAtoms() == smilesTest->atomCount);
    TEST_ASSERT(localVars.smilesMol->getNumBonds() == smilesTest->bondCount);

    {
      std::string expectedSdfName =fName + ".expected.sdf";
      std::string outMolStr="";
      try
      {
        outMolStr =  MolToMolBlock(*localVars.smilesMol, true, 0, true, true);
      }
      catch (const RDKit::KekulizeException &e)
      {
        outMolStr = "";
      }
      catch(...)
      {
        throw;  // re-trhow the error if not a kekule error
      }
      if (outMolStr == "")
        outMolStr =  MolToMolBlock(*localVars.smilesMol, true, 0,false, true);  // try without kekule'ing

      std::stringstream  expectedMolStr;
      std::ifstream  in;
      in.open(expectedSdfName);
      expectedMolStr << in.rdbuf();
      std::string expectedStr = expectedMolStr.str();

      TEST_ASSERT(expectedStr == outMolStr); 
    }

    {
      std::string expectedMrvName =fName + ".expected.mrv";
      std::string outMolStr="";
      try
      {
        outMolStr = MolToMrvBlock(*localVars.smilesMol,true, -1, true);
      }
      catch (const RDKit::KekulizeException &e)
      {
        outMolStr = "";
      }
      catch(...)
      {
        throw;  // re-trhow the error if not a kekule error
      }
      if (outMolStr == "")
        outMolStr = MolToMrvBlock(*localVars.smilesMol,true, -1,false);  // try without kekule'ing

      std::stringstream  expectedMolStr;
      std::ifstream  in;
      in.open(expectedMrvName);
      expectedMolStr << in.rdbuf();
      std::string expectedStr = expectedMolStr.str();

      TEST_ASSERT(expectedStr == outMolStr); 
    }
    BOOST_LOG(rdInfoLog) << "done" << std::endl;   
  } 
  catch (const std::exception &e) 
  {
    if(smilesTest->expectedResult != false)
        throw;
    return;

  }

  TEST_ASSERT(smilesTest->expectedResult == true);
}

void testMarvin(const MolOrRxnTest *molOrRxnTest)
{
  BOOST_LOG(rdInfoLog) << "testing marvin parsing" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName = rdbase + "/Code/GraphMol/MarvinParse/test_data/" + molOrRxnTest->fileName;

  class LocalVars  // protext against mem leak on error
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
      auto rxnTest = (RxnTest *)molOrRxnTest;

      // check for errors

      unsigned int nWarn=0, nError=0;

      TEST_ASSERT(rxn!= NULL);

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

      TEST_ASSERT(nWarn == rxnTest->warnings);
      TEST_ASSERT(nError == rxnTest->errors);

      {
        std::string outMolStr = ChemicalReactionToRxnBlock(*rxn, false, true);
        std::string expectedRxnName = fName + ".expected.rxn";

        std::stringstream  expectedMolStr;
        std::ifstream  in;
        in.open(expectedRxnName);
        expectedMolStr << in.rdbuf();
        std::string expectedStr = expectedMolStr.str();

        TEST_ASSERT(expectedStr == outMolStr); 
      }

      {
        std::string outMolStr = ChemicalReactionToMrvBlock(*rxn);
        std::string expectedRxnName = fName + ".expected.mrv";
        std::stringstream  expectedMolStr;
        std::ifstream  in;
        in.open(expectedRxnName);
        expectedMolStr << in.rdbuf();
        std::string expectedStr = expectedMolStr.str();

        TEST_ASSERT(expectedStr == outMolStr); 
      }

      BOOST_LOG(rdInfoLog) << "done" << std::endl;
    }
    else
    {   
      // mol  test

      auto mol = (RWMol *)localVars.molOrRxn;
      auto molTest = (MolTest *)molOrRxnTest;
      TEST_ASSERT(mol != NULL);

      TEST_ASSERT(mol->getNumAtoms() == molTest->atomCount)
      TEST_ASSERT(mol->getNumBonds() == molTest->bondCount)

     
      {
        std::string expectedMrvName =fName + ".expected.sdf";
        std::string outMolStr="";
        try
        {
          outMolStr =  MolToMolBlock(*mol, true, 0, true, true);
        }
        catch (const RDKit::KekulizeException &e)
        {
          outMolStr = "";
        }
        catch(...)
        {
          throw;  // re-trhow the error if not a kekule error
        }
        if (outMolStr == "")
          outMolStr =  MolToMolBlock(*mol, true, 0,false, true);  // try without kekule'ing

        std::stringstream  expectedMolStr;
        std::ifstream  in;
        in.open(expectedMrvName);
        expectedMolStr << in.rdbuf();
        std::string expectedStr = expectedMolStr.str();

        TEST_ASSERT(expectedStr == outMolStr); 
      }

      {
        std::string expectedMrvName =fName + ".expected.mrv";

        std::string outMolStr="";
        try
        {
          outMolStr = MolToMrvBlock(*mol,true, -1, true);
        }
        catch (const RDKit::KekulizeException &e)
        {
          outMolStr = "";
        }
        catch(...)
        {
          throw;  // re-trhow the error if not a kekule error
        }
        if (outMolStr == "")
          outMolStr = MolToMrvBlock(*mol,true, -1,false);  // try without kekule'ing
  
        std::stringstream  expectedMolStr;
        std::ifstream  in;
        in.open(expectedMrvName);
        expectedMolStr << in.rdbuf();
        std::string expectedStr = expectedMolStr.str();

        TEST_ASSERT(expectedStr == outMolStr); 
      }

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
    MolTest("radical_value.mrv", true, LoadAsMolOrRxn, 3, 2)
    ,  MolTest("emptyOneLineAtomList.mrv", true, LoadAsMolOrRxn, 0, 0)
    ,  MolTest("mrvValence_value.mrv", true, LoadAsMolOrRxn, 3, 2)
    ,  MolTest("ChiralTest2.mrv", true, LoadAsMolOrRxn, 46, 47)
    ,  MolTest("ChiralTest.mrv", true, LoadAsMolOrRxn, 8, 7)
    ,  MolTest("SnCl2.mrv", true, LoadAsMolOrRxn, 3, 2)
    ,  MolTest("SnH2Cl2.mrv", true, LoadAsMolOrRxn, 3, 2)
    ,  MolTest("marvin01.mrv", true, LoadAsMolOrRxn, 11, 11)
    ,  MolTest("marvin01.mrv", true, LoadAsMol, 11, 11)
    ,  MolTest("marvin01.mrv", false, LoadAsRxn, 11, 11)  // should fail
    ,  MolTest("marvin02.mrv", true, LoadAsMolOrRxn, 9, 9)
    ,  MolTest("marvin07.mrv", true, LoadAsMolOrRxn, 12, 11)
    ,  MolTest("marvin10.mrv", true, LoadAsMolOrRxn, 10, 10)
    ,  MolTest("marvin06.mrv", true, LoadAsMolOrRxn, 11, 11)
    ,  MolTest("marvin12.mrv", true, LoadAsMolOrRxn, 31, 33)
    ,  MolTest("EmptyMol.mrv", true, LoadAsMolOrRxn, 0, 0)
    ,  MolTest("Sparse.mrv", true, LoadAsMolOrRxn, 0, 0)
    ,  MolTest("Sparse2.mrv", true, LoadAsMolOrRxn, 0, 0)
    ,  MolTest("Sparse3.mrv", true, LoadAsMolOrRxn, 0, 0)
    ,  MolTest("MarvinNoCoords.mrv", true, LoadAsMolOrRxn, 6, 6)
    ,  MolTest("aspirin.mrv", true, LoadAsMolOrRxn, 13, 13)
    ,  MolTest("MarvinStereoGroupsZeros.mrv", true, LoadAsMolOrRxn, 8, 8)
    ,  MolTest("MarvinStereoGroupsAbs.mrv", true, LoadAsMolOrRxn, 8, 8)
    ,  MolTest("triphenylphosphine.mrv", true, LoadAsMolOrRxn, 19, 21)
    ,  MolTest("MarvinOldSuperGroupTest.mrv", true, LoadAsMolOrRxn, 89, 93)
    ,  MolTest("RadicalTests.mrv", true, LoadAsMolOrRxn, 9, 9)
    ,  MolTest("AnyBond.mrv", true, LoadAsMolOrRxn, 4, 3)
    ,  MolTest("cisBenzene.mrv", true, LoadAsMolOrRxn, 6, 6)
    ,  MolTest("DativeBond.mrv", true, LoadAsMolOrRxn, 6, 5)
    ,  MolTest("MultipleSgroup.mrv", true, LoadAsMolOrRxn, 75, 74)
    ,  MolTest("SgroupExpanded.mrv", true, LoadAsMolOrRxn, 5, 4)
    ,  MolTest("SgroupMultAttach.mrv", true, LoadAsMolOrRxn, 44, 45)  
    ,  MolTest("MarvinMissingX2.mrv", true, LoadAsMolOrRxn, 12, 11)  
    ,  MolTest("MarvinMissingY2.mrv", true, LoadAsMolOrRxn, 12, 11)  
    ,  MolTest("DataSgroup.mrv", true, LoadAsMolOrRxn, 7, 6)
    ,  MolTest("MulticenterSgroup.mrv", true, LoadAsMolOrRxn, 17, 16)
    ,  MolTest("GenericSgroup.mrv", true, LoadAsMolOrRxn, 13, 13)
    ,  MolTest("MonomerSgroup.mrv", true, LoadAsMolOrRxn, 4, 3)
    ,  MolTest("marvin03.mrv", false, LoadAsMolOrRxn, 31,33)  // should fail - this is a reaction
    ,  MolTest("MarvinBadMissingMolID.mrv", false, LoadAsMolOrRxn, 12, 11)  // should fail - no molId
    ,  MolTest("MarvinBadMissingAtomID.mrv", false, LoadAsMolOrRxn, 12, 11)  // should fail - missing atom Id
    ,  MolTest("MarvinBadX2.mrv", false, LoadAsMolOrRxn, 12, 11)  // should fail -
    ,  MolTest("MarvinBadY2.mrv", false, LoadAsMolOrRxn, 12, 11)  // should fail -
    ,  MolTest("MarvinBadStereoGroupsAbs.mrv", false, LoadAsMolOrRxn, 8, 8) // should fail -
    ,  MolTest("MarvinBadElementType.mrv", false, LoadAsMolOrRxn, 12, 11)  // should fail -
    ,  MolTest("MarvinBadMissingBondID.mrv", false, LoadAsMolOrRxn, 12, 11)  // should fail -
    ,  MolTest("MarvinBadMissingBondAtomRefs", false, LoadAsMolOrRxn, 12, 11)  // should fail -
    ,  MolTest("MarvinBadMissingBondOrder.mrv", false, LoadAsMolOrRxn, 12, 11)  // should fail -
    ,  MolTest("MarvinBadMissingSruMolID.mrv", false, LoadAsMolOrRxn, 12, 11)  // should fail -
    ,  MolTest("MarvinBadMissingSruID.mrv", false, LoadAsMolOrRxn, 12, 11)  // should fail -
    ,  MolTest("MarvinBadMissingSruRole.mrv", false, LoadAsMolOrRxn, 12, 11)  // should fail -
    ,  MolTest("MarvinBadMissingSruAtomRef.mrv", false, LoadAsMolOrRxn, 12, 11)  // should fail -
    ,  MolTest("MarvinBadMissingSruTitle.mrv", false, LoadAsMolOrRxn, 12, 11)  // should fail -
    ,  MolTest("MarvinBadSruAtomRef.mrv", false, LoadAsMolOrRxn, 12, 11)  // should fail -
    ,  MolTest("MarvinBadSruID.mrv", false, LoadAsMolOrRxn, 12, 11)  // should fail -
    ,  MolTest("MarvinBadSruRole.mrv", false, LoadAsMolOrRxn, 12, 11)  // should fail -
    ,  MolTest("MarvinBadSruAtomRef.mrv", false, LoadAsMolOrRxn, 12, 11)  // should fail -
    ,  MolTest("MarvinBadSruAtomRef.mrv", false, LoadAsMolOrRxn, 12, 11)  // should fail -
    ,  MolTest("MarvinBadSruConnect.mrv", false, LoadAsMolOrRxn, 12, 11)  // should fail -
    ,  MolTest("MarvinBadSupAttachAtom.mrv", false, LoadAsMolOrRxn, 9,  9)  // should fail -
    ,  MolTest("MarvinBadSupAttachBond.mrv", false, LoadAsMolOrRxn, 9,  9)  // should fail -
    ,  MolTest("MarvinBadSupAttachOrder.mrv", false, LoadAsMolOrRxn, 9,  9)  // should fail -
    ,  MolTest("MarvinBadSupAttachAtom.mrv", false, LoadAsMolOrRxn, 9,  9)  // should fail -
    ,  MolTest("MarvinBadSupAttachAtom.mrv", false, LoadAsMolOrRxn, 9,  9)  // should fail -
    ,  MolTest("MarvinBadSupAttachAtom.mrv", false, LoadAsMolOrRxn, 9,  9)  // should fail -
    ,  MolTest("MarvinBadSupMissingAttachBond.mrv", false, LoadAsMolOrRxn, 9,  9)  // should fail -
    ,  MolTest("MarvinBadSupMissingAttachOrder.mrv", false, LoadAsMolOrRxn,9,  9)  // should fail -
  };

  for (std::list<MolTest>::const_iterator it = molFileNames.begin();
        it != molFileNames.end(); ++it) 
  {
    BOOST_LOG(rdInfoLog) << "Test: " << it->fileName << std::endl;

    printf("Test\n\n %s\n\n", it->fileName.c_str());
    testMarvin(&*it);
  }

  // now the reactions

  std::list<RxnTest> rxnFileNames
  {
      RxnTest("bondArray_node.mrv", true, LoadAsMolOrRxn, 2, 4, 1, 3, 0),
      RxnTest("marvin03.mrv", true, LoadAsMolOrRxn, 1, 1, 1, 2, 0),
      RxnTest("marvin03.mrv", true, LoadAsRxn, 1, 1, 1, 2, 0),
      RxnTest("marvin03.mrv", false, LoadAsMol, 1, 1, 1, 2, 0),  // should fail
      RxnTest("marvin04.mrv", true, LoadAsMolOrRxn, 2, 1, 2, 4, 0),
      RxnTest("marvin08.mrv", true, LoadAsMolOrRxn, 2, 3, 2, 4, 0),
      RxnTest("marvin09.mrv", true, LoadAsMolOrRxn, 2, 3, 2, 4, 0),
      RxnTest("marvin11.mrv", true, LoadAsMolOrRxn, 2, 0, 1, 0, 0),
      RxnTest("marvin05.mrv", true, LoadAsMolOrRxn, 2, 1, 1, 3, 0),
      RxnTest("EmptyRxn.mrv", true, LoadAsMolOrRxn, 0, 0, 0, 0, 0),
      RxnTest("condition_coordinates_mpoint.mrv", true, LoadAsMolOrRxn, 1, 0, 1, 0, 0), 
      RxnTest("marvin01.mrv", false, LoadAsMolOrRxn, 2, 1, 1, 3, 0)  // should fail - this is a mol file

      ,
      RxnTest("aspirineSynthesisWithAttributes.mrv", true, LoadAsMolOrRxn,
              2, 0, 1, 3, 0)  // should fail - this is a mol file
  };

  for (std::list<RxnTest>::const_iterator it = rxnFileNames.begin();
        it != rxnFileNames.end(); ++it) 
  {
    printf("Test\n\n %s\n\n", it->fileName.c_str());
    testMarvin(&*it);
  }

  // now smiles tests

   std::list<SmilesTest> smiFileNames
   {
      SmilesTest("Smiles1","N[C@@H]([O-])c1cc[13c]cc1",true,9,9)
   };

  for (std::list<SmilesTest>::const_iterator it = smiFileNames.begin(); it != smiFileNames.end(); ++it) 
  {
    printf("Test\n\n %s\n\n", it->name.c_str());
    testSmilesToMarvin(&*it);
  }
}


int main(int argc, char *argv[]) 
{
  (void)argc;
  (void)argv;
  //  std::locale::global(std::locale("de_DE.UTF-8"));


  RDLog::InitLogs();
  BOOST_LOG(rdInfoLog) << " ---- Running with POSIX locale ----- "
                        << std::endl;

  RunTests();  // run with C locale

  return 0;
}
