//
//  Copyright (C) 2004-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_MOL_FILE_STEREOCHEM_H
#define RD_MOL_FILE_STEREOCHEM_H

#include <GraphMol/RDKitBase.h>

namespace RDKit {

RDKIT_FILEPARSERS_EXPORT void DetectAtomStereoChemistry(RWMol &mol,
                                                        const Conformer *conf);
//! deprecated, please use MolOps::detectBondStereoChemistry instead
RDKIT_FILEPARSERS_EXPORT void DetectBondStereoChemistry(ROMol &mol,
                                                        const Conformer *conf);

RDKIT_FILEPARSERS_EXPORT void WedgeMolBonds(ROMol &mol, const Conformer *conf);
RDKIT_FILEPARSERS_EXPORT void WedgeBond(Bond *bond, unsigned int fromAtomIdx,
                                        const Conformer *conf);

struct RDKIT_FILEPARSERS_EXPORT StereoBondThresholds {
  const static unsigned DBL_BOND_NO_STEREO =
      1000;  //!< neighboring double bond without stereo info
  const static unsigned DBL_BOND_SPECIFIED_STEREO =
      10000;  //!< neighboring double bond with stereo specified
  const static unsigned CHIRAL_ATOM =
      100000;  //!< atom with specified chirality
  const static unsigned DIRECTION_SET =
      1000000;  //!< single bond with the direction already set
};
//! set wavy bonds around double bonds with STEREOANY stereo
/*!
 \param mol molecule to be modified
 \param clearDoubleBondFlags when this is true flags for unknown double bond
   stereo will also be removed.
 \param addWhenImpossible if nonzero a neighboring single bond will be made
 wavy even if it connects to a chiral center or double bond with specified
 stereo. one example of this would be the middle double bond in
 C/C=C/C=C/C=C/C (if that's set to STEREOANY after constructing the molecule)
 Otherwise, no wavy bond will be set
*/
RDKIT_FILEPARSERS_EXPORT void addWavyBondsForStereoAny(
    ROMol &mol, bool clearDoubleBondFlags = true,
    unsigned addWhenImpossible = StereoBondThresholds::DBL_BOND_NO_STEREO);

//! picks the bonds which should be wedged
/// \returns a map from bond idx -> controlling atom idx
RDKIT_FILEPARSERS_EXPORT INT_MAP_INT pickBondsToWedge(const ROMol &mol);
//! deprecated, please use MolOps::clearSingleBondDirFlags instead
RDKIT_FILEPARSERS_EXPORT void ClearSingleBondDirFlags(ROMol &mol);
RDKIT_FILEPARSERS_EXPORT Bond::BondDir DetermineBondWedgeState(
    const Bond *bond, unsigned int fromAtomIdx, const Conformer *conf);
RDKIT_FILEPARSERS_EXPORT Bond::BondDir DetermineBondWedgeState(
    const Bond *bond, const INT_MAP_INT &wedgeBonds, const Conformer *conf);
//! Remove MolBlock bond wedging information from molecule.
/*!
 \param mol: molecule to modify
 */
RDKIT_FILEPARSERS_EXPORT void clearMolBlockWedgingInfo(ROMol &mol);
//! Invert bond wedging information read from a mol block (if present).
/*!
 \param mol: molecule to modify
 */
RDKIT_FILEPARSERS_EXPORT void invertMolBlockWedgingInfo(ROMol &mol);

//! Set double bonds with unspecified stereo to STEREOANY and add wavy bonds
//! to
///  potential stereocenters with unspecified chirality
RDKIT_FILEPARSERS_EXPORT void markUnspecifiedStereoAsUnknown(ROMol &mol,
                                                             int confId = -1);
//! gets stereo info for a bond
/*!
 \param bond: bond to check
 \param wedgeBonds - the list of bonds to have wedges
 \param conf -  Conformer to use
 \param dirCode - receives the dircode for the bond
 \param reverse - receives the reverse flag
 \param explicitUnknownDoubleBondOnly - is set, the code for unknown stereo is
 only returned if it was exlicility set witha wiggle bond
 */

void GetMolFileBondStereoInfo(const Bond *bond, const INT_MAP_INT &wedgeBonds,
                              const Conformer *conf, int &dirCode,
                              bool &reverse,
                              bool explicitUnknownDoubleBondOnly = false);

//! Forces use of atom wedging from MolBlock if present.
/*
 \param mol: molecule to have its wedges altered
 */

RDKIT_FILEPARSERS_EXPORT void reapplyMolBlockWedging(
    RWMol &mol, bool throwAromaticError = false);

//! used when search for an wedged bond was found to be aromatic
//! the cxsmiles uses this
class RDKIT_RDGENERAL_EXPORT AromaticException : public std::runtime_error {
 public:
  //! construct with an error message
  explicit AromaticException(const char *msg)
      : std::runtime_error("AromaticException"), _msg(msg) {}
  //! construct with an error message
  explicit AromaticException(const std::string msg)
      : std::runtime_error("AromaticException"), _msg(msg) {}
  //! get the error message
  const char *what() const noexcept override { return _msg.c_str(); }
  ~AromaticException() noexcept override = default;

 private:
  std::string _msg;
};

}  // namespace RDKit
#endif