//
//
//  Copyright (C) 2020 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <string>

#include <RDBoost/Wrap.h>
#include <RDBoost/python.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/CIPLabeler/CIPLabeler.h>
#include <GraphMol/FileParsers/FileParsers.h>


namespace python = boost::python;
using RDKit::CIPLabeler::assignCIPLabels;

void assignCIPLabelsWrapHelper(RDKit::ROMol &mol,
                               const python::object &atomsToLabel,
                               const python::object &bondsToLabel,
                               unsigned int maxRecursiveIterations) {
  auto atoms = pythonObjectToDynBitset(atomsToLabel, mol.getNumAtoms());
  auto bonds = pythonObjectToDynBitset(bondsToLabel, mol.getNumBonds());

  // If both atoms and bonds are None, assign all the mol.
  if (!atomsToLabel && !bondsToLabel) {
    atoms.set();
    bonds.set();
  }

  assignCIPLabels(mol, atoms, bonds,maxRecursiveIterations);
}

BOOST_PYTHON_MODULE(rdCIPLabeler) {
  python::scope().attr("__doc__") =
      "Module containing a function to assign stereochemical labels based "
      "on an accurate CIP rules implementation. This algoritm is a port "
      "of https://github.com/SiMolecule/centres, which was originally "
      "written by John Mayfield. The original algorithm is described in:\n\n"
      "Hanson, R. M., Musacchio, S., Mayfield, J. W., Vainio, M. J., Yerin, "
      "A., Redkin, D.\nAlgorithmic Analysis of Cahn--Ingold--Prelog Rules of "
      "Stereochemistry:\nProposals for Revised Rules and a Guide for Machine "
      "Implementation.\nJ. Chem. Inf. Model. 2018, 58, 1755-1765.\n";

  std::string docString =
      "New implementation of Stereo assignment using a true CIP ranking.\n\
       On return:  The molecule to contains CIP flags\n\
       Errors:  when maxRecursiveIterations is exceeded, throws a MaxIterationsExceeded error\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - atomsToLabel: (optional) list of atoms to label\n\
    - bondsToLabel: (optional) list of bonds to label\n\
    - maxRecursiveIterations: (optional) protects against pseudo-infinite recursiion for highly symmetical structures.\n\
       A value of 1,250,000 take about 1 second.  Most strucutres requires less than 10,000 iterations (0 = default - no limit)\n";

  python::def(
      "AssignCIPLabels", assignCIPLabelsWrapHelper,
      (python::arg("mol"), 
       python::arg("atomsToLabel") = python::object(),
       python::arg("bondsToLabel") = python::object(),
       python::arg("maxRecursiveIterations") = 0),
      docString.c_str());
}
