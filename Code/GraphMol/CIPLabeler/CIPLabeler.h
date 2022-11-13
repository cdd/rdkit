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
#pragma once

#include <boost/dynamic_bitset.hpp>

#include <RDGeneral/export.h>

#include <chrono>

namespace RDKit {

class ROMol;

namespace CIPLabeler {

RDKIT_CIPLABELER_EXPORT std::chrono::system_clock::time_point &getCipTimeOut();

/*
  Some very symmetrical mols can cause pseudo infinite processing (e.g. bukcy-ball)
  To avoid this a timeout can be set by the caller with setTimeOutInSeconds()
  If that timeout is exceeded, the folloowing error is thrown
*/

class CipLabelerTimeoutException : public std::runtime_error {
 public:
  explicit CipLabelerTimeoutException()
      : std::runtime_error("CipLabelerTimeoutException"){};
};



/**
 * Calculate Stereochemical labels based on an accurate implementation
 * of the CIP rules.
 *
 * This is a C++ port of https://github.com/SiMolecule/centres, which was
 * originally written by John Mayfield in Java. The original algorithm was
 * described in:
 *
 * Hanson, R. M., Musacchio, S., Mayfield, J. W., Vainio, M. J., Yerin, A.,
 * Redkin, D. Algorithmic Analysis of Cahn--Ingold--Prelog Rules of
 * Stereochemistry: Proposals for Revised Rules and a Guide for Machine
 * Implementation. J. Chem. Inf. Model. 2018, 58, 1755-1765.
 *
 *   \param mol - the molecule to be labelled.
 *
 *   \note only atoms with chiral tags and double bonds with proper
 *          bond directions will be labelled.
 *   \note Labels will be stored under the common_properties::_CIPCode
 *          property of the relevant atoms/bonds.
 */
RDKIT_CIPLABELER_EXPORT void assignCIPLabels(ROMol &mol, unsigned int timeOutInSeconds = 0);

/**
 * Overload that allows selecting which atoms and/or bonds will be labeled.
 *
 *   \param mol - the molecule to be labelled.
 *
 *   \param atoms - bitset with the atom indexes to be labeled.
 *
 *   \param bonds - bitset with the bond indexes to be labeled.
 *
 */
RDKIT_CIPLABELER_EXPORT void assignCIPLabels(
    ROMol &mol, const boost::dynamic_bitset<> &atoms,
    const boost::dynamic_bitset<> &bonds, 
    unsigned int timeOutInSeconds = 0);
}  // namespace CIPLabeler
}  // namespace RDKit
