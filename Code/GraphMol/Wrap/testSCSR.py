#
#  Copyright (C) 2024 Tad Hurst
#         All Rights Reserved
#
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

#
import os
import unittest

from rdkit import Chem
from rdkit.Chem import RDConfig


class TestCase(unittest.TestCase):

  def setUp(self):
    pass

  def testSCSR(self):
    """Test the SCSR system"""

    ofile = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'macromols', 'Triplet.mol')
    with open(ofile) as inf:
      scsrBlock = inf.read()

    molFromMacroMolParams = Chem.MolFromMacroMolParams()
    molFromMacroMolParams.includeLeavingGroups = True
    molFromMacroMolParams.outputSgroups = True

  
    for mol in (Chem.MolFromSCSRBlock(scsrBlock, False, False, molFromMacroMolParams),
                Chem.MolFromSCSRFile(ofile, False, False, molFromMacroMolParams)):

      self.assertTrue(mol.GetNumAtoms() == 30)
      sgs = Chem.GetMolSubstanceGroups(mol)
      self.assertEqual(len(sgs), 3)

    # check defaults:
    for mol in (Chem.MolFromSCSRBlock(scsrBlock), Chem.MolFromSCSRFile(ofile)):
      self.assertTrue(mol.GetNumAtoms() == 30)
      sgs = Chem.GetMolSubstanceGroups(mol)
      self.assertEqual(len(sgs), 3)

  def testSCSRRna(self):
    """Test the SCSR system with and RNA double strand"""

    ofile = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'macromols', 'DnaTest.mol')
    with open(ofile) as inf:
      scsrBlock = inf.read()

    molFromMacroMolParams = Chem.MolFromMacroMolParams()
    molFromMacroMolParams.includeLeavingGroups = True
    molFromMacroMolParams.outputSgroups = True

    mol = Chem.MolFromSCSRBlock(scsrBlock, False, False, molFromMacroMolParams,True)

    self.assertTrue(mol.GetNumAtoms() == 254)
    self.assertTrue(mol.GetNumBonds() == 300)
    sgs = Chem.GetMolSubstanceGroups(mol)
    print ("len(sgs): " + str(len(sgs)))
    self.assertTrue(len(sgs) == 34)

    mol = Chem.MolFromSCSRBlock(scsrBlock, False, False, molFromMacroMolParams, False)

    self.assertTrue(mol.GetNumAtoms() == 254)
    self.assertTrue(mol.GetNumBonds() == 282)
    sgs = Chem.GetMolSubstanceGroups(mol)
    self.assertTrue(len(sgs) == 34)

  def testThreeLetterCodes(self):
    """Test the SCSR system with three letter codes"""

    ofile = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'macromols', 'PepTla.mol')
    with open(ofile) as inf:
      scsrBlock = inf.read()

    molFromMacroMolParams = Chem.MolFromMacroMolParams()
    molFromMacroMolParams.includeLeavingGroups = True
    molFromMacroMolParams.outputSgroups = True

    mol = Chem.MolFromSCSRBlock(scsrBlock, False, False, molFromMacroMolParams)

    self.assertEqual(mol.GetNumAtoms(), 26)
    self.assertEqual(mol.GetNumBonds(), 25)
    sgs = Chem.GetMolSubstanceGroups(mol)

    self.assertEqual(len(sgs), 4)
    sgs[0].GetProp('LABEL')
    self.assertEqual(sgs[0].GetProp('LABEL'), 'Ala')
          
if __name__ == '__main__':
  unittest.main()
