//
//  Copyright (C) 2004-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//
#include <list>
#include <RDGeneral/RDLog.h>
#include <GraphMol/Chirality.h>
#include "MolFileStereochem.h"
#include <Geometry/point.h>
#include <boost/dynamic_bitset.hpp>
#include <algorithm>
#include <RDGeneral/Ranking.h>
#include <RDGeneral/FileParseException.h>

constexpr double REALLY_SMALL_BOND_LEN = 0.0000001;

namespace RDKit {
void GetMolFileBondStereoInfo(const Bond *bond, const INT_MAP_INT &wedgeBonds,
                              const Conformer *conf, int &dirCode,
                              bool &reverse);

typedef std::list<double> DOUBLE_LIST;

void WedgeBond(Bond *bond, unsigned int fromAtomIdx, const Conformer *conf) {
  PRECONDITION(bond, "no bond");
  PRECONDITION(conf, "no conformer");
  PRECONDITION(&conf->getOwningMol() == &bond->getOwningMol(),
               "bond and conformer do not belong to same molecule");
  if (bond->getBondType() != Bond::SINGLE) {
    return;
  }
  Bond::BondDir dir = DetermineBondWedgeState(bond, fromAtomIdx, conf);
  if (dir == Bond::BEGINWEDGE || dir == Bond::BEGINDASH) {
    bond->setBondDir(dir);
  }
}

void WedgeMolBonds(ROMol &mol, const Conformer *conf) {
  PRECONDITION(conf, "no conformer");
  auto wedgeBonds = pickBondsToWedge(mol);
  for (auto bond : mol.bonds()) {
    if (bond->getBondType() == Bond::SINGLE) {
      Bond::BondDir dir = DetermineBondWedgeState(bond, wedgeBonds, conf);
      if (dir == Bond::BEGINWEDGE || dir == Bond::BEGINDASH) {
        bond->setBondDir(dir);

        // it is possible that this
        // wedging was determined by a chiral atom at the end of the
        // bond (instead of at the beginning). In this case we need to
        // reverse the begin and end atoms for the bond
        auto wbi = wedgeBonds.find(bond->getIdx());
        if (wbi != wedgeBonds.end() &&
            static_cast<unsigned int>(wbi->second) != bond->getBeginAtomIdx()) {
          auto tmp = bond->getBeginAtomIdx();
          bond->setBeginAtomIdx(bond->getEndAtomIdx());
          bond->setEndAtomIdx(tmp);
        }
      }
    }
  }
}

std::tuple<unsigned int, unsigned int, unsigned int, unsigned int>
getDoubleBondPresence(const ROMol &mol, const Atom &atom) {
  unsigned int hasDouble = 0;
  unsigned int hasKnownDouble = 0;
  unsigned int hasAnyDouble = 0;
  unsigned int hasAtropIsomer = 0;
  for (const auto bond : mol.atomBonds(&atom)) {
    if (bond->getBondType() == Bond::BondType::DOUBLE) {
      ++hasDouble;
      if (bond->getStereo() == Bond::BondStereo::STEREOANY) {
        ++hasAnyDouble;
      } else if (bond->getStereo() > Bond::BondStereo::STEREOANY) {
        ++hasKnownDouble;
      }
    } else if (bond->getBondType() == Bond::BondType::SINGLE &&
               (bond->getStereo() == Bond::BondStereo::STEREOATROPCW ||
                bond->getStereo() == Bond::BondStereo::STEREOATROPCCW)) {
      ++hasAtropIsomer;
    }
  }
  return std::make_tuple(hasDouble, hasKnownDouble, hasAnyDouble,
                         hasAtropIsomer);
}

std::pair<bool, INT_VECT> countChiralNbrs(const ROMol &mol, int noNbrs) {
  // we need ring information; make sure findSSSR has been called before
  // if not call now
  if (!mol.getRingInfo()->isInitialized()) {
    MolOps::findSSSR(mol);
  }

  INT_VECT nChiralNbrs(mol.getNumAtoms(), noNbrs);

  // start by looking for bonds that are already wedged
  for (const auto bond : mol.bonds()) {
    if (bond->getBondDir() == Bond::BEGINWEDGE ||
        bond->getBondDir() == Bond::BEGINDASH ||
        bond->getBondDir() == Bond::UNKNOWN) {
      if (bond->getBeginAtom()->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
          bond->getBeginAtom()->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) {
        nChiralNbrs[bond->getBeginAtomIdx()] = noNbrs + 1;
      } else if (bond->getEndAtom()->getChiralTag() ==
                     Atom::CHI_TETRAHEDRAL_CW ||
                 bond->getEndAtom()->getChiralTag() ==
                     Atom::CHI_TETRAHEDRAL_CCW) {
        nChiralNbrs[bond->getEndAtomIdx()] = noNbrs + 1;
      }
    }
  }

  // now rank atoms by the number of chiral neighbors or Hs they have:
  bool chiNbrs = false;
  for (const auto at : mol.atoms()) {
    if (nChiralNbrs[at->getIdx()] > noNbrs) {
      // std::cerr << " SKIPPING1: " << at->getIdx() << std::endl;
      continue;
    }
    Atom::ChiralType type = at->getChiralTag();
    if (type != Atom::CHI_TETRAHEDRAL_CW && type != Atom::CHI_TETRAHEDRAL_CCW) {
      continue;
    }
    nChiralNbrs[at->getIdx()] = 0;
    chiNbrs = true;
    for (const auto nat : mol.atomNeighbors(at)) {
      if (nat->getAtomicNum() == 1) {
        // special case: it's an H... we weight these especially high:
        nChiralNbrs[at->getIdx()] -= 10;
        continue;
      }
      type = nat->getChiralTag();
      if (type != Atom::CHI_TETRAHEDRAL_CW &&
          type != Atom::CHI_TETRAHEDRAL_CCW) {
        continue;
      }
      nChiralNbrs[at->getIdx()] -= 1;
    }
  }
  return std::pair<bool, INT_VECT>(chiNbrs, nChiralNbrs);
}

// picks a bond for atom that we will wedge when we write the mol file
// returns idx of that bond.
int pickBondToWedge(const Atom *atom, const ROMol &mol,
                    const INT_VECT &nChiralNbrs, const INT_MAP_INT &resSoFar,
                    int noNbrs) {
  PRECONDITION(atom, "no atom");

  // here is what we are going to do
  // - at each chiral center look for a bond that is begins at the atom and
  //   is not yet picked to be wedged for a different chiral center,
  //   preferring bonds to Hs
  // - if we do not find a bond that begins at the chiral center - we will
  // take
  //   the first bond that is not yet picked by any other chiral centers
  // we use the orders calculated above to determine which order to do the
  // wedging

  // If we call wedgeMolBonds() on a fragment, it can happen that we end up
  // with atoms that don't have enough neighbors. Those are going to cause
  // problems, so just bail here.
  std::vector<std::pair<int, int>> nbrScores;
  for (const auto bond : mol.atomBonds(atom)) {
    // can only wedge single bonds:
    if (bond->getBondType() != Bond::SINGLE) {
      continue;
    }

    int bid = bond->getIdx();
    if (resSoFar.find(bid) == resSoFar.end()) {
      // very strong preference for Hs:
      if (bond->getOtherAtom(atom)->getAtomicNum() == 1) {
        nbrScores.emplace_back(-1000000,
                               bid);  // lower than anything else can be
        continue;
      }
      // prefer lower atomic numbers with lower degrees and no specified
      // chirality:
      const Atom *oatom = bond->getOtherAtom(atom);
      int nbrScore = oatom->getAtomicNum() + 100 * oatom->getDegree() +
                     1000 * ((oatom->getChiralTag() != Atom::CHI_UNSPECIFIED));
      // prefer neighbors that are nonchiral or have as few chiral neighbors
      // as possible:
      int oIdx = oatom->getIdx();
      if (nChiralNbrs[oIdx] < noNbrs) {
        // the counts are negative, so we have to subtract them off
        nbrScore -= 100000 * nChiralNbrs[oIdx];
      }
      // prefer bonds to non-ring atoms:
      nbrScore += 10000 * mol.getRingInfo()->numAtomRings(oIdx);
      // prefer non-ring bonds;
      nbrScore += 20000 * mol.getRingInfo()->numBondRings(bid);
      // prefer bonds to atoms which don't have a double bond from them
      unsigned int hasDoubleBond;       // is a double bond there?
      unsigned int hasKnownDoubleBond;  // is specified stereo there?
      unsigned int hasAnyDoubleBond;    // is STEREOANY there?
      unsigned int hasAtropisomer;      // on end of an atropisomer
      std::tie(hasDoubleBond, hasKnownDoubleBond, hasAnyDoubleBond,
               hasAtropisomer) = getDoubleBondPresence(mol, *oatom);
      nbrScore += 11000 * hasDoubleBond;
      nbrScore += 12000 * hasAtropisomer;
      nbrScore += 12000 * hasKnownDoubleBond;
      nbrScore += 23000 * hasAnyDoubleBond;

      // std::cerr << "    nrbScore: " << idx << " - " << oIdx << " : "
      //           << nbrScore << " nChiralNbrs: " << nChiralNbrs[oIdx]
      //           << std::endl;
      nbrScores.emplace_back(nbrScore, bid);
    }
  }
  // There's still one situation where this whole thing can fail: an unlucky
  // situation where all neighbors of all neighbors of an atom are chiral and
  // that atom ends up being the last one picked for stereochem assignment.
  // This also happens in cases where the chiral atom doesn't have all of its
  // neighbors (like when working with partially sanitized fragments)
  //
  // We'll bail here by returning -1
  if (nbrScores.empty()) {
    return -1;
  }
  std::sort(nbrScores.begin(), nbrScores.end(), Rankers::pairLess);
  return nbrScores[0].second;
}

INT_MAP_INT pickBondsToWedge(const ROMol &mol) {
  // returns map of bondIdx -> bond begin atom for those bonds that
  // need wedging.
  std::vector<unsigned int> indices(mol.getNumAtoms());
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    indices[i] = i;
  }
  static int noNbrs = 100;
  std::pair<bool, INT_VECT> retVal = countChiralNbrs(mol, noNbrs);
  bool chiNbrs = retVal.first;
  INT_VECT nChiralNbrs = retVal.second;
  if (chiNbrs) {
    std::sort(indices.begin(), indices.end(), [&](auto i1, auto i2) {
      return nChiralNbrs[i1] < nChiralNbrs[i2];
    });
  }
#if 0
  std::cerr << "  nbrs: ";
  std::copy(nChiralNbrs.begin(), nChiralNbrs.end(),
            std::ostream_iterator<int>(std::cerr, " "));
  std::cerr << std::endl;
  std::cerr << "  order: ";
  std::copy(indices.begin(), indices.end(),
            std::ostream_iterator<int>(std::cerr, " "));
  std::cerr << std::endl;
#endif
  INT_MAP_INT res;
  for (auto idx : indices) {
    if (nChiralNbrs[idx] > noNbrs) {
      // std::cerr << " SKIPPING2: " << idx << std::endl;
      continue;  // already have a wedged bond here
    }
    const Atom *atom = mol.getAtomWithIdx(idx);
    Atom::ChiralType type = atom->getChiralTag();
    // the indices are ordered such that all chiral atoms come first. If
    // this has no chiral flag, we can stop the whole loop:
    if (type != Atom::CHI_TETRAHEDRAL_CW && type != Atom::CHI_TETRAHEDRAL_CCW) {
      break;
    }
    int bnd = pickBondToWedge(atom, mol, nChiralNbrs, res, noNbrs);
    if (bnd >= 0) {
      res[bnd] = idx;
    }
  }
  return res;
}

std::vector<Bond *> getBondNeighbors(ROMol &mol, const Bond &bond) {
  std::vector<Bond *> res;
  for (auto nbrBond : mol.atomBonds(bond.getBeginAtom())) {
    if (nbrBond == &bond) {
      continue;
    }
    res.push_back(nbrBond);
  }
  for (auto nbrBond : mol.atomBonds(bond.getEndAtom())) {
    if (nbrBond == &bond) {
      continue;
    }
    res.push_back(nbrBond);
  }
  return res;
}

const Atom *getNonsharedAtom(const Bond &bond1, const Bond &bond2) {
  if (bond1.getBeginAtomIdx() == bond2.getBeginAtomIdx() ||
      bond1.getBeginAtomIdx() == bond2.getEndAtomIdx()) {
    return bond1.getEndAtom();
  } else if (bond1.getEndAtomIdx() == bond2.getBeginAtomIdx() ||
             bond1.getEndAtomIdx() == bond2.getEndAtomIdx()) {
    return bond1.getBeginAtom();
  }
  POSTCONDITION(0, "bonds don't share an atom");
}

const unsigned StereoBondThresholds::DBL_BOND_NO_STEREO;
const unsigned StereoBondThresholds::DBL_BOND_SPECIFIED_STEREO;
const unsigned StereoBondThresholds::CHIRAL_ATOM;
const unsigned StereoBondThresholds::DIRECTION_SET;

// a note on the way the StereoBondThresholds are used:
//  the penalties are all 1/10th of the corresponding threshold, so
//  the penalty for being connected to a chiral atom is
//  StereoBondThresholds::CHIRAL_ATOM/10
//  This allows us to just add up the penalties for a particular
//  single bond and still use one set of thresholds - an individual
//  single bond will never have any particular penalty term applied
//  more than a couple of times
//
void addWavyBondsForStereoAny(ROMol &mol, bool clearDoubleBondFlags,
                              unsigned addWhenImpossible) {
  std::vector<int> singleBondScores(mol.getNumBonds(), 0);
  // used to store the double bond neighbors, if any, of each single bond
  std::map<unsigned, std::vector<unsigned>> singleBondNeighbors;
  boost::dynamic_bitset<> doubleBondsToSet(mol.getNumBonds());
  // mark single bonds adjacent to double bonds
  for (const auto dblBond : mol.bonds()) {
    if (dblBond->getBondType() != Bond::BondType::DOUBLE) {
      continue;
    }
    if (dblBond->getStereo() == Bond::BondStereo::STEREOANY) {
      doubleBondsToSet.set(dblBond->getIdx());
    }
    for (auto singleBond : getBondNeighbors(mol, *dblBond)) {
      if (singleBond->getBondType() != Bond::BondType::SINGLE) {
        continue;
      }
      // NOTE: we could make this canonical by initializing scores to the
      // canonical atom ranks
      int score = singleBondScores[singleBond->getIdx()];
      ++score;

      // penalty for having a direction already set
      if (singleBond->getBondDir() != Bond::BondDir::NONE) {
        score += StereoBondThresholds::DIRECTION_SET / 10;
      }

      // penalties from the double bond itself:

      // penalize being adjacent to a double bond with empty stereo:
      if (dblBond->getStereo() == Bond::BondStereo::STEREONONE) {
        score += StereoBondThresholds::DBL_BOND_NO_STEREO / 10;
      } else if (dblBond->getStereo() > Bond::BondStereo::STEREOANY) {
        // penalize being adjacent to a double bond with specified stereo:
        score += StereoBondThresholds::DBL_BOND_SPECIFIED_STEREO / 10;
      }

      // atom-related penalties
      auto otherAtom = getNonsharedAtom(*singleBond, *dblBond);
      // favor atoms with smaller numbers of neighbors:
      score += 10 * otherAtom->getDegree();
      // penalty for being adjacent to an atom with specified stereo
      if (otherAtom->getChiralTag() != Atom::ChiralType::CHI_UNSPECIFIED &&
          otherAtom->getChiralTag() != Atom::ChiralType::CHI_OTHER) {
        score += StereoBondThresholds::CHIRAL_ATOM / 10;
      }
      singleBondScores[singleBond->getIdx()] = score;
      if (dblBond->getStereo() == Bond::BondStereo::STEREOANY) {
        singleBondNeighbors[singleBond->getIdx()].push_back(dblBond->getIdx());
      }
    }
  }
  std::vector<std::tuple<int, unsigned int, size_t>> sortedScores;
  for (size_t i = 0; i < mol.getNumBonds(); ++i) {
    auto score = singleBondScores[i];
    if (!score) {
      continue;
    }
    sortedScores.push_back(
        std::make_tuple(-1 * singleBondNeighbors[i].size(), score, i));
  }
  std::sort(sortedScores.begin(), sortedScores.end());
  for (const auto &tpl : sortedScores) {
    // FIX: check if dir is already set
    for (auto dblBondIdx : singleBondNeighbors[std::get<2>(tpl)]) {
      if (doubleBondsToSet[dblBondIdx]) {
        if (addWhenImpossible) {
          if (std::get<1>(tpl) > addWhenImpossible) {
            continue;
          }
        } else if (std::get<1>(tpl) >
                   StereoBondThresholds::DBL_BOND_NO_STEREO) {
          BOOST_LOG(rdWarningLog)
              << "Setting wavy bond flag on bond " << std::get<2>(tpl)
              << " which may make other stereo info ambiguous" << std::endl;
        }
        mol.getBondWithIdx(std::get<2>(tpl))
            ->setBondDir(Bond::BondDir::UNKNOWN);
        if (clearDoubleBondFlags) {
          auto dblBond = mol.getBondWithIdx(dblBondIdx);
          if (dblBond->getBondDir() == Bond::BondDir::EITHERDOUBLE) {
            dblBond->setBondDir(Bond::BondDir::NONE);
          }
          dblBond->setStereo(Bond::BondStereo::STEREONONE);
        }
        doubleBondsToSet.reset(dblBondIdx);
      }
    }
  }
  if (addWhenImpossible) {
    if (doubleBondsToSet.count()) {
      std::stringstream sstr;
      sstr << " unable to set wavy bonds for double bonds:";
      for (size_t i = 0; i < mol.getNumBonds(); ++i) {
        if (doubleBondsToSet[i]) {
          sstr << " " << i;
        }
      }
      BOOST_LOG(rdWarningLog) << sstr.str() << std::endl;
    }
  }
}
//
// Determine bond wedge state
///
Bond::BondDir DetermineBondWedgeState(const Bond *bond,
                                      unsigned int fromAtomIdx,
                                      const Conformer *conf) {
  PRECONDITION(bond, "no bond");
  PRECONDITION(bond->getBondType() == Bond::SINGLE,
               "bad bond order for wedging");
  const ROMol *mol = &(bond->getOwningMol());
  PRECONDITION(mol, "no mol");

  Bond::BondDir res = bond->getBondDir();
  if (!conf) {
    return res;
  }

  Atom *atom, *bondAtom;  // = bond->getBeginAtom();
  if (bond->getBeginAtom()->getIdx() == fromAtomIdx) {
    atom = bond->getBeginAtom();
    bondAtom = bond->getEndAtom();
  } else {
    atom = bond->getEndAtom();
    bondAtom = bond->getBeginAtom();
  }

  Atom::ChiralType chiralType = atom->getChiralTag();
  CHECK_INVARIANT(chiralType == Atom::CHI_TETRAHEDRAL_CW ||
                      chiralType == Atom::CHI_TETRAHEDRAL_CCW,
                  "");

  // if we got this far, we really need to think about it:
  INT_LIST neighborBondIndices;
  DOUBLE_LIST neighborBondAngles;
  RDGeom::Point3D centerLoc, tmpPt;
  centerLoc = conf->getAtomPos(atom->getIdx());
  tmpPt = conf->getAtomPos(bondAtom->getIdx());
  centerLoc.z = 0.0;
  tmpPt.z = 0.0;
  RDGeom::Point3D refVect = centerLoc.directionVector(tmpPt);

  neighborBondIndices.push_back(bond->getIdx());
  neighborBondAngles.push_back(0.0);

  ROMol::OEDGE_ITER beg, end;
  boost::tie(beg, end) = mol->getAtomBonds(atom);
  while (beg != end) {
    const Bond *nbrBond = (*mol)[*beg];
    Atom *otherAtom = nbrBond->getOtherAtom(atom);
    if (nbrBond != bond) {
      tmpPt = conf->getAtomPos(otherAtom->getIdx());
      tmpPt.z = 0.0;
      RDGeom::Point3D tmpVect = centerLoc.directionVector(tmpPt);
      double angle = refVect.signedAngleTo(tmpVect);
      if (angle < 0.0) {
        angle += 2. * M_PI;
      }
      auto nbrIt = neighborBondIndices.begin();
      auto angleIt = neighborBondAngles.begin();
      // find the location of this neighbor in our angle-sorted list
      // of neighbors:
      while (angleIt != neighborBondAngles.end() && angle > (*angleIt)) {
        ++angleIt;
        ++nbrIt;
      }
      neighborBondAngles.insert(angleIt, angle);
      neighborBondIndices.insert(nbrIt, nbrBond->getIdx());
    }
    ++beg;
  }

  // at this point, neighborBondIndices contains a list of bond
  // indices from the central atom.  They are arranged starting
  // at the reference bond in CCW order (based on the current
  // depiction).
  int nSwaps = atom->getPerturbationOrder(neighborBondIndices);

  // in the case of three-coordinated atoms we may have to worry about
  // the location of the implicit hydrogen - Issue 209
  // Check if we have one of these situation
  //
  //      0        1 0 2
  //      *         \*/
  //  1 - C - 2      C
  //
  // here the hydrogen will be between 1 and 2 and we need to add an
  // additional swap
  if (neighborBondAngles.size() == 3) {
    // three coordinated
    auto angleIt = neighborBondAngles.begin();
    ++angleIt;  // the first is the 0 (or reference bond - we will ignoire
                // that
    double angle1 = (*angleIt);
    ++angleIt;
    double angle2 = (*angleIt);
    if (angle2 - angle1 >= (M_PI - 1e-4)) {
      // we have the above situation
      nSwaps++;
    }
  }

#ifdef VERBOSE_STEREOCHEM
  BOOST_LOG(rdDebugLog) << "--------- " << nSwaps << std::endl;
  std::copy(neighborBondIndices.begin(), neighborBondIndices.end(),
            std::ostream_iterator<int>(BOOST_LOG(rdDebugLog), " "));
  BOOST_LOG(rdDebugLog) << std::endl;
  std::copy(neighborBondAngles.begin(), neighborBondAngles.end(),
            std::ostream_iterator<double>(BOOST_LOG(rdDebugLog), " "));
  BOOST_LOG(rdDebugLog) << std::endl;
#endif
  if (chiralType == Atom::CHI_TETRAHEDRAL_CCW) {
    if (nSwaps % 2 == 1) {  // ^ reverse) {
      res = Bond::BEGINDASH;
    } else {
      res = Bond::BEGINWEDGE;
    }
  } else {
    if (nSwaps % 2 == 1) {  // ^ reverse) {
      res = Bond::BEGINWEDGE;
    } else {
      res = Bond::BEGINDASH;
    }
  }

  return res;
}
Bond::BondDir DetermineBondWedgeState(const Bond *bond,
                                      const INT_MAP_INT &wedgeBonds,
                                      const Conformer *conf) {
  PRECONDITION(bond, "no bond");
  int bid = bond->getIdx();
  auto wbi = wedgeBonds.find(bid);
  if (wbi == wedgeBonds.end()) {
    return bond->getBondDir();
  }

  unsigned int waid = wbi->second;
  return DetermineBondWedgeState(bond, waid, conf);
}

// handles stereochem markers set by the Mol file parser and
// converts them to the RD standard:

void DetectAtomStereoChemistry(RWMol &mol, const Conformer *conf) {
  PRECONDITION(conf, "no conformer");
  PRECONDITION(&(conf->getOwningMol()) == &mol,
               "conformer does not belong to molecule");
  MolOps::assignChiralTypesFromBondDirs(mol, conf->getId(), true);
}

void GetAtropisomerAtomsAndBonds(const Bond *bond, Atom *atoms[2],
                                 std::vector<Bond *> bonds[2],
                                 const ROMol &mol) {
  PRECONDITION(bond, "no bond");
  atoms[0] = bond->getBeginAtom();
  atoms[1] = bond->getEndAtom();

  // get the one or two bonds on each end

  for (int bondAtomIndex = 0; bondAtomIndex < 2; ++bondAtomIndex) {
    for (const auto nbrBond : mol.atomBonds(atoms[bondAtomIndex])) {
      if (nbrBond == bond) {
        continue;  // a bond is NOT its own neighbor
      }
      bonds[bondAtomIndex].push_back(nbrBond);
    }

    // make sure the bond with this lowest atom is is first

    if (bonds[bondAtomIndex].size() == 2 &&
        bonds[bondAtomIndex][1]->getOtherAtom(atoms[bondAtomIndex])->getIdx() <
            bonds[bondAtomIndex][0]
                ->getOtherAtom(atoms[bondAtomIndex])
                ->getIdx()) {
      std::swap(bonds[bondAtomIndex][0], bonds[bondAtomIndex][1]);
    }
  }
}

bool getBondFrameOfReference(const Atom *const atoms[2], const Conformer *conf,
                             RDGeom::Point3D &xAxis, RDGeom::Point3D &yAxis,
                             RDGeom::Point3D &zAxis) {
  // create a frame of reference that has its X-axis along the atrop bond
  // for 2D confs, the yAxis is in the 2D plane and the zAxis is perpendicular
  // to that plane) for 3D confs  the yAxis and the zAxis are arbitrary.

  PRECONDITION(atoms[0], "bad atom");
  PRECONDITION(atoms[1], "bad atom");

  xAxis = conf->getAtomPos(atoms[1]->getIdx()) -
          conf->getAtomPos(atoms[0]->getIdx());
  if (xAxis.length() < REALLY_SMALL_BOND_LEN) {
    return false;  // bond len is xero
  }
  xAxis.normalize();
  if (!conf->is3D()) {
    yAxis = RDGeom::Point3D(-xAxis.y, xAxis.x, 0);
    yAxis.normalize();
    zAxis = RDGeom::Point3D(0.0, 0.0, 1.0);
    return true;
  }

  // here for 3D conf

  if (xAxis.x > REALLY_SMALL_BOND_LEN || xAxis.y > REALLY_SMALL_BOND_LEN) {
    zAxis = RDGeom::Point3D(-xAxis.y, xAxis.x,
                            0);  // temp z axis - used to find yAxis
  } else {
    zAxis = RDGeom::Point3D(xAxis.z, xAxis.z,
                            0);  // temp z axis - used to find yAxis
  }

  yAxis = zAxis.crossProduct(xAxis);
  zAxis = xAxis.crossProduct(yAxis);
  yAxis.normalize();
  zAxis.normalize();

  return true;
}

Bond::BondDir getBondDirForAtropisomer2d(RDGeom::Point3D bondVecs[2],
                                         Bond::BondStereo bondStereo,
                                         unsigned int whichEnd,
                                         unsigned int whichBond) {
  PRECONDITION(whichEnd <= 1, "whichEnd must be 0 or 1");
  PRECONDITION(whichBond <= 1, "whichBond must be 0 or 1");
  PRECONDITION(bondStereo == Bond::BondStereo::STEREOATROPCW ||
                   bondStereo == Bond::BondStereo::STEREOATROPCCW,
               "bondStereo must be BondAtropisomerCW or BondAtropisomerCCW");

  int flips = 0;
  if (bondStereo == Bond::BondStereo::STEREOATROPCCW) {
    ++flips;
  }
  if (whichBond == 1) {
    ++flips;
  }
  if (whichEnd == 1) {
    ++flips;
  }
  if (bondVecs[1 - whichEnd].y < 0) {
    ++flips;
  }
  // if the OTHER end if negative for the low index bond vec, it
  // is a flip

  return flips % 2 ? Bond::BEGINDASH : Bond::BEGINWEDGE;
}
Bond::BondDir getBondDirForAtropisomer3d(Bond *whichBond,
                                         const Conformer *conf) {
  // for 3D we mark it as wedge or hash depending on the z-value of the bond
  // vector
  //  IT really doesn't matter since we ignore these except as MARKERS for
  //  which bonds are atropisomer bonds
  if ((conf->getAtomPos(whichBond->getEndAtom()->getIdx()).z -
       conf->getAtomPos(whichBond->getBeginAtom()->getIdx()).z) >
      REALLY_SMALL_BOND_LEN) {
    return Bond::BondDir::BEGINWEDGE;
  } else {
    return Bond::BondDir::BEGINDASH;
  }
}

bool getAtropIsomerEndVect(const Atom *mainBondAtom,
                           const std::vector<Bond *> endBonds,
                           const RDGeom::Point3D &yAxis,
                           const RDGeom::Point3D &zAxis, const Conformer *conf,
                           RDGeom::Point3D &bondVec) {
  PRECONDITION(mainBondAtom, "bad bond");
  PRECONDITION(endBonds.size() > 0 && endBonds.size() < 3, "bad bond size");
  PRECONDITION(endBonds[0], "bad first bond");
  PRECONDITION(endBonds.size() == 1 || endBonds[1], "bad second bond");

  bondVec =
      conf->getAtomPos(endBonds[0]->getOtherAtom(mainBondAtom)->getIdx()) -
      conf->getAtomPos(mainBondAtom->getIdx());  // in old frame of reference

  bondVec = RDGeom::Point3D(0.0, bondVec.dotProduct(yAxis),
                            bondVec.dotProduct(zAxis));  // in new frame

  // make sure the other atom is on the other side

  if (endBonds.size() == 2) {
    RDGeom::Point3D otherVec =
        conf->getAtomPos(endBonds[1]->getOtherAtom(mainBondAtom)->getIdx()) -
        conf->getAtomPos(mainBondAtom->getIdx());  // in old frame of reference
    otherVec = RDGeom::Point3D(0.0, otherVec.dotProduct(yAxis),
                               otherVec.dotProduct(zAxis));  // in new frame

    if (bondVec.length() < REALLY_SMALL_BOND_LEN) {
      bondVec = -otherVec;  // put it on the other side of otherVec
    } else if (bondVec.dotProduct(otherVec) > REALLY_SMALL_BOND_LEN) {
      // the product of dotproducts (y-values) should be
      // negative (or at least zero)
      BOOST_LOG(rdWarningLog)
          << "Both bonds on one end of an atropisomer are on the same side - atoms is : "
          << mainBondAtom->getIdx() << std::endl;
      return false;
    }
  }
  if (bondVec.length() < REALLY_SMALL_BOND_LEN) {
    BOOST_LOG(rdWarningLog)
        << "Could not find a bond on one end of an atropisomer that is not co-linear - atoms are : "
        << mainBondAtom->getIdx() << std::endl;
    return false;
  }

  bondVec.normalize();
  return true;
}
void DetectAtropisomerChiralityOneBond(Bond *bond, ROMol &mol,
                                       const Conformer *conf) {
  // the approach is this:
  // we will view the system along the line from the potential atropisomer
  // bond, from atom1 to atom 2 and we do a coordinate transformation to
  // the plane of reference where that vector, from a1 to a2, is the x-AXIS.
  // For 2D, the y axis is in the 2D plane, and the zaxis is perpendicaul to
  // the 2D plane For 3D, the Y and Z axes are taken arbitrarily to form
  // a right-handed system with the X-axis.
  //  atoms 1 and 2 each have one or two bonds out from the main potential
  //  atrop bond. for each end of the main bond, we find a vector to reprent
  //  the neighbor atom with the smallest index as its projection onto the
  //  x=0 plane.
  // (In 2d, this projection is on the y-AXIS for the end that does NOT have a
  // wedge/hash bond, and  on the z axis - out of the plane - for the end that
  // does have a wedge/hash). The chirality is recorded as the direction we
  // rotate from, atom 1's projection to atom2's proejection - either clockwise
  // or counter clockwise

  PRECONDITION(bond, "bad bond");

  Atom *atoms[2];
  // std::vector<std::vector<Bond *>> bonds;
  std::vector<Bond *> bonds[2];  // one vector for each end - each one
                                 // should end up with 1 ro 2 entries

  GetAtropisomerAtomsAndBonds(bond, atoms, bonds, mol);

  // make sure we do not have wiggle bonds

  for (auto endBonds : bonds) {
    for (auto endBond : endBonds) {
      if (endBond->getBondDir() == Bond::UNKNOWN) {
        return;  // not an atropisomer
      }
    }
  }

  // create a frame of reference that has its X-axis along the atrop bond

  RDGeom::Point3D xAxis, yAxis, zAxis;

  if (!getBondFrameOfReference(atoms, conf, xAxis, yAxis, zAxis)) {
    // connot percieve atroisomer
    BOOST_LOG(rdWarningLog)
        << "Failed to get a frame of reference along an atropisomer bond - atoms are: "
        << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx() << std::endl;
    return;
  }
  RDGeom::Point3D bondVecs[2];  // one bond vector from each end of the
                                // potential atropisomer bond

  bool foundWedgeOrHash = false;
  for (int bondAtomIndex = 0; bondAtomIndex < 2; ++bondAtomIndex) {
    // if the conf is 2D, we use the wedge bonds to set the coords for the
    // projected vector onto the xAxis perpendicular plane (looking down
    // the atrop bond )

    if (!conf->is3D()) {
      // get the wedge dir for this end of the bond
      // if the first bond 1 has a bondDir, use it
      // if the second bond has a bond dir use the opposite of if
      // if both bonds have a dir, make sure they are different

      auto bond1Dir = bonds[bondAtomIndex][0]->getBondDir();
      if (bond1Dir != Bond::BEGINWEDGE && bond1Dir != Bond::BEGINDASH) {
        bond1Dir = Bond::NONE;  //  we dont care if it any thing else
      }
      auto bond2Dir = bonds[bondAtomIndex].size() == 2
                          ? bonds[bondAtomIndex][1]->getBondDir()
                          : Bond::NONE;
      if (bond2Dir != Bond::BEGINWEDGE && bond2Dir != Bond::BEGINDASH) {
        bond2Dir = Bond::NONE;
      }

      // if both are set to a direction, they must NOT be the same - one
      // must be a dash and the other a hash

      if (bond1Dir != Bond::NONE && bond2Dir != Bond::NONE &&
          bond1Dir == bond2Dir) {
        BOOST_LOG(rdWarningLog)
            << "The bonds on one end of an atropisomer are both UP or both DOWN - atoms are: "
            << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
            << std::endl;
        return;
      }

      if (bond1Dir != Bond::NONE || bond2Dir != Bond::NONE) {
        if (foundWedgeOrHash) {
          // already found a wedge on the other end of the bond - this is
          // indeterminate
          BOOST_LOG(rdWarningLog)
              << "An atropisomer has a wedge or hash bond on both ends - atoms are: "
              << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
              << std::endl;
          return;
        }
        foundWedgeOrHash = true;

        if (bond1Dir == Bond::BEGINWEDGE || bond2Dir == Bond::BEGINDASH) {
          bondVecs[bondAtomIndex] = RDGeom::Point3D(0.0, 0.0, 1.0);
        } else if (bond1Dir == Bond::BEGINDASH ||
                   bond2Dir == Bond::BEGINWEDGE) {
          bondVecs[bondAtomIndex] = RDGeom::Point3D(0.0, 0.0, -1.0);
        }
      } else {
        // the coords are taken as the part of the 2d coord that is
        // in the the x=0 plane.  (for 2d, that means the y-value). The
        // result is normalized
        if (!getAtropIsomerEndVect(atoms[bondAtomIndex], bonds[bondAtomIndex],
                                   yAxis, zAxis, conf,
                                   bondVecs[bondAtomIndex])) {
          return;
        }
      }
    } else {  // the conf is 3D
      // to be considered, one or more neighbor bonds must have a wedge or
      // hash

      // find the projection of the bond(s) on this end in the frame of
      // reference's  x=0  plane
      RDGeom::Point3D tempBondVec =
          conf->getAtomPos(bonds[bondAtomIndex][0]
                               ->getOtherAtom(atoms[bondAtomIndex])
                               ->getIdx()) -
          conf->getAtomPos(atoms[bondAtomIndex]->getIdx());
      bondVecs[bondAtomIndex] = RDGeom::Point3D(
          0.0, tempBondVec.dotProduct(yAxis), tempBondVec.dotProduct(zAxis));

      if (bonds[bondAtomIndex].size() == 2) {
        tempBondVec = conf->getAtomPos(bonds[bondAtomIndex][1]
                                           ->getOtherAtom(atoms[bondAtomIndex])
                                           ->getIdx()) -
                      conf->getAtomPos(atoms[bondAtomIndex]->getIdx());

        // get the projection of the 2nd bond on the x=0 plane

        RDGeom::Point3D otherBondVec = RDGeom::Point3D(
            0.0, tempBondVec.dotProduct(yAxis), tempBondVec.dotProduct(zAxis));

        // if the first atom is co-linear with the main atrop bond, use
        // the opposite of the 2nd atom

        if (bondVecs[bondAtomIndex].length() < REALLY_SMALL_BOND_LEN) {
          bondVecs[bondAtomIndex] =
              -otherBondVec;  // note - it might still be co-linear- this
                              // is checked below
        } else if (bondVecs[bondAtomIndex].dotProduct(otherBondVec) >
                   REALLY_SMALL_BOND_LEN) {
          BOOST_LOG(rdWarningLog)
              << "Both bonds on one end of an atropisomer are on the same side - atoms are: "
              << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
              << std::endl;
          return;
        }
      }

      if (bondVecs[bondAtomIndex].length() < REALLY_SMALL_BOND_LEN) {
        BOOST_LOG(rdWarningLog)
            << "Failed to find a bond on one end of an atropisomer that is NOT co-linear - atoms are: "
            << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
            << std::endl;
        return;
      }
    }
  }

  auto crossProduct = bondVecs[0].crossProduct(bondVecs[1]);

  if (crossProduct.x > REALLY_SMALL_BOND_LEN) {
    bond->setStereo(Bond::BondStereo::STEREOATROPCCW);
  } else if (crossProduct.x < -REALLY_SMALL_BOND_LEN) {
    bond->setStereo(Bond::BondStereo::STEREOATROPCW);
  } else {
    BOOST_LOG(rdWarningLog)
        << "The 2 defining bonds for an atropisomer are co-planar - atoms are: "
        << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx() << std::endl;
    return;
  }
}

void DetectAtropisomerChirality(ROMol &mol, const Conformer *conf) {
  PRECONDITION(conf, "no conformer");
  PRECONDITION(&(conf->getOwningMol()) == &mol,
               "conformer does not belong to molecule");

  if (!mol.getRingInfo()->isInitialized()) {
    MolOps::symmetrizeSSSR(mol);
  }

  std::set<Bond *> bondsToTry;

  for (auto bond : mol.bonds()) {
    if (bond->getBondType() == Bond::SINGLE &&
        (bond->getBondDir() == Bond::BondDir::BEGINDASH ||
         bond->getBondDir() == Bond::BondDir::BEGINWEDGE)) {
      for (const auto &nbrBond : mol.atomBonds(bond->getBeginAtom())) {
        if (nbrBond == bond) {
          continue;  // a bond is NOT its own neighbor
        }
        bondsToTry.insert(nbrBond);
      }
    }
  }

  // const RingInfo *ri = bond->getOwningMol().getRingInfo();
  const RingInfo *ri = mol.getRingInfo();
  for (auto bondToTry : bondsToTry) {
    if (ri->numBondRings(bondToTry->getIdx()) > 0) {
      continue;
    }

    if (bondToTry->getBeginAtom()->getImplicitValence() == -1) {
      bondToTry->getBeginAtom()->calcExplicitValence(false);
      bondToTry->getBeginAtom()->calcImplicitValence(false);
    }
    if (bondToTry->getEndAtom()->getImplicitValence() == -1) {
      bondToTry->getEndAtom()->calcExplicitValence(false);
      bondToTry->getEndAtom()->calcImplicitValence(false);
    }
    if (bondToTry->getBondType() != Bond::SINGLE ||
        bondToTry->getStereo() == Bond::BondStereo::STEREOANY ||
        bondToTry->getBeginAtom()->getTotalDegree() < 2 ||
        bondToTry->getEndAtom()->getTotalDegree() < 2 ||
        bondToTry->getBeginAtom()->getTotalDegree() > 3 ||
        bondToTry->getEndAtom()->getTotalDegree() > 3) {
      continue;
    }

    DetectAtropisomerChiralityOneBond(bondToTry, mol, conf);
  }
}

void WedgeBondFromAtropisomerOneBond2d(Bond *bond, const ROMol &mol,
                                       const Conformer *conf,
                                       const INT_MAP_INT &wedgeBonds) {
  PRECONDITION(bond, "no bond");
  Atom *atoms[2];
  // std::vector<std::vector<Bond *>> bonds;
  std::vector<Bond *> bonds[2];  // one vector for each end - each one
                                 // should end up with 1 ro 2 entries

  GetAtropisomerAtomsAndBonds(bond, atoms, bonds, mol);

  //  make sure we do not have wiggle bonds

  for (auto endBonds : bonds) {
    for (auto endBond : endBonds) {
      if (endBond->getBondDir() == Bond::UNKNOWN) {
        return;  // not an atropisomer)
      }
    }
  }

  // create a frame of reference that has its X-axis along the atrop bond

  RDGeom::Point3D xAxis, yAxis, zAxis;

  if (!getBondFrameOfReference(atoms, conf, xAxis, yAxis, zAxis)) {
    // connot percieve atroisomer bond

    BOOST_LOG(rdWarningLog)
        << "Cound not get a frame of reference for an atropisomer bond - atoms are: "
        << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx() << std::endl;
    return;
  }

  RDGeom::Point3D bondVecs[2];  // one bond vector from each end of the
                                // potential atropisome bond

  for (int bondAtomIndex = 0; bondAtomIndex < 2; ++bondAtomIndex) {
    // find a vector to represent the lowest numbered atom on each end
    // this vector is NOT the bond vector, but is y-value in the bond
    // frame or reference

    if (!getAtropIsomerEndVect(atoms[bondAtomIndex], bonds[bondAtomIndex],
                               yAxis, zAxis, conf, bondVecs[bondAtomIndex])) {
      return;
    }

    if (bondVecs[bondAtomIndex].length() < REALLY_SMALL_BOND_LEN) {
      // did not find a non-colinear bond

      BOOST_LOG(rdWarningLog)
          << "Failed to get a representative vector for the defining bond of an atropisomer - atoms are: "
          << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
          << std::endl;
      return;
    }
  }

  // first see if any candidate bond is already set to a wedge or hash
  // if so, we will use that bond as a wedge or hash

  std::vector<int> useBondsAtEnd[2];
  bool foundBondDir = false;

  for (unsigned int whichEnd = 0; whichEnd < 2; ++whichEnd) {
    for (unsigned int whichBond = 0; whichBond < bonds[whichEnd].size();
         ++whichBond) {
      auto bondDir = bonds[whichEnd][whichBond]->getBondDir();

      // see if it is a wedge or hash and its origin is the atom in the
      // main bond

      if ((bondDir == Bond::BEGINWEDGE || bondDir == Bond::BEGINDASH) &&
          bonds[whichEnd][whichBond]->getBeginAtom() == atoms[whichEnd] &&
          bonds[whichEnd][whichBond]->getBondType() == Bond::SINGLE) {
        useBondsAtEnd[whichEnd].push_back(whichBond);
        foundBondDir = true;
      }
    }
  }

  if (foundBondDir) {
    for (unsigned int whichEnd = 0; whichEnd < 2; ++whichEnd) {
      for (unsigned int whichBondIndex = 0;
           whichBondIndex < useBondsAtEnd[whichEnd].size(); ++whichBondIndex) {
        bonds[whichEnd][useBondsAtEnd[whichEnd][whichBondIndex]]->setBondDir(
            getBondDirForAtropisomer2d(
                bondVecs, bond->getStereo(), whichEnd,
                useBondsAtEnd[whichEnd][whichBondIndex]));
      }
    }

    return;
  }

  // did not find a good bond dir - pick one to use
  // we would like to have one that is not in a ring, and will be a dash

  const RingInfo *ri = bond->getOwningMol().getRingInfo();

  int bestBondEnd = -1, bestBondNumber = -1;
  unsigned int bestRingCount = INT_MAX;
  Bond::BondDir bestBondDir = Bond::BondDir::NONE;
  for (unsigned int whichEnd = 0; whichEnd < 2; ++whichEnd) {
    for (unsigned int whichBond = 0; whichBond < bonds[whichEnd].size();
         ++whichBond) {
      auto bondToTry = bonds[whichEnd][whichBond];

      if (bondToTry->getBondType() != Bond::BondType::SINGLE ||
          wedgeBonds.find(bondToTry->getIdx()) != wedgeBonds.end()) {
        continue;  // must be a single bond and not already spoken for by a
                   // chiral center
      }

      if (bondToTry->getBondDir() != Bond::BondDir::NONE) {
        if (bondToTry->getBeginAtom()->getIdx() == atoms[whichEnd]->getIdx()) {
          BOOST_LOG(rdWarningLog)
              << "Wedge or hash bond found on atropisomer where not expected - atoms are: "
              << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
              << std::endl;
          return;
        } else {
          continue;  // wedge or hash bond affecting the OTHER atom =
                     // perhaps a chiral center
        }
      }
      auto ringCount = ri->numBondRings(bondToTry->getIdx());
      if (ringCount > bestRingCount) {
        continue;
      }

      else if (ringCount < bestRingCount) {
        bestBondEnd = whichEnd;
        bestBondNumber = whichBond;
        bestRingCount = ringCount;
        bestBondDir = getBondDirForAtropisomer2d(bondVecs, bond->getStereo(),
                                                 whichEnd, whichBond);
      } else {
        auto bondDir = getBondDirForAtropisomer2d(bondVecs, bond->getStereo(),
                                                  whichEnd, whichBond);
        if (bestBondDir == Bond::BondDir::NONE ||
            (bestBondDir == Bond::BondDir::BEGINDASH &&
             bondDir == Bond::BondDir::BEGINWEDGE)) {
          bestBondEnd = whichEnd;
          bestBondNumber = whichBond;
          bestRingCount = ringCount;
          bestBondDir = bondDir;
        }
      }
    }
  }

  if (bestBondEnd >= 0)  // we found a good one
  {
    // make sure the atoms on the bond are in the right order for the
    // wedge/hash the atom on the end of the main bond must be listed
    // first for the wedge/has bond

    auto bestBond = bonds[bestBondEnd][bestBondNumber];
    if (bestBond->getBeginAtom() != atoms[bestBondEnd]) {
      bestBond->setEndAtom(bestBond->getBeginAtom());
      bestBond->setBeginAtom(atoms[bestBondEnd]);
    }
    bonds[bestBondEnd][bestBondNumber]->setBondDir(bestBondDir);
  } else {
    BOOST_LOG(rdWarningLog)
        << "Failed to find a good bond to set as UP or DOWN for an atropisomer - atoms are: "
        << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx() << std::endl;
    return;
  }

  return;
}

void WedgeBondFromAtropisomerOneBond3d(Bond *bond, const ROMol &mol,
                                       const Conformer *conf,
                                       const INT_MAP_INT &wedgeBonds) {
  PRECONDITION(bond, "bad bond");
  Atom *atoms[2];
  // std::vector<std::vector<Bond *>> bonds;
  std::vector<Bond *> bonds[2];  // one vector for each end - each one
                                 // should end up with 1 ro 2 entries

  GetAtropisomerAtomsAndBonds(bond, atoms, bonds, mol);

  //  make sure we do not have wiggle bonds

  for (auto endBonds : bonds) {
    for (auto endBond : endBonds) {
      if (endBond->getBondDir() == Bond::UNKNOWN) {
        return;  // not an atropisomer)
      }
    }
  }

  // first see if any candidate bond is already set to a wedge or hash
  // if so, we will use that bond as a wedge or hash

  std::vector<Bond *> useBonds;

  for (unsigned int whichEnd = 0; whichEnd < 2; ++whichEnd) {
    for (unsigned int whichBond = 0; whichBond < bonds[whichEnd].size();
         ++whichBond) {
      auto bond = bonds[whichEnd][whichBond];
      auto bondDir = bond->getBondDir();

      // see if it is a wedge or hash and its origin is the atom in the
      // main bond

      if ((bondDir == Bond::BEGINWEDGE || bondDir == Bond::BEGINDASH) &&
          bond->getBeginAtom() == atoms[whichEnd] &&
          bond->getBondType() == Bond::SINGLE) {
        useBonds.push_back(bond);
      }
    }
  }

  // the following may seem redundat, since we just found the useBonds based on
  // their bond dir PRESENCE, but this endures that the values are correct.

  if (useBonds.size() > 0) {
    for (auto useBond : useBonds) {
      useBond->setBondDir(getBondDirForAtropisomer3d(useBond, conf));
    }

    return;
  }

  // did not find a used bond dir - pick one to use
  // we would like to have one that is not in a ring, and will be a dash

  const RingInfo *ri = bond->getOwningMol().getRingInfo();

  Bond *bestBond = nullptr;
  int bestBondEnd = -1;
  unsigned int bestRingCount = INT_MAX;
  Bond::BondDir bestBondDir = Bond::BondDir::NONE;
  for (unsigned int whichEnd = 0; whichEnd < 2; ++whichEnd) {
    for (unsigned int whichBond = 0; whichBond < bonds[whichEnd].size();
         ++whichBond) {
      auto bondToTry = bonds[whichEnd][whichBond];

      // cannot use a bond that is not single, nor if it is already slated
      // to be used for a chiral center

      if (bondToTry->getBondType() != Bond::BondType::SINGLE ||
          wedgeBonds.find(bond->getIdx()) != wedgeBonds.end()) {
        continue;  // must be a single bond and not already spoken for by a
                   // chiral center
      }

      // make sure the atoms on the bond are in the right order for the
      // wedge/hash the atom on the end of the main bond must be listed
      // first

      if (bondToTry->getBeginAtom() != atoms[whichEnd]) {
        bondToTry->setEndAtom(bondToTry->getBeginAtom());
        bondToTry->setBeginAtom(atoms[whichEnd]);
      }

      if (bondToTry->getBondDir() != Bond::BondDir::NONE) {
        if (bondToTry->getBeginAtom()->getIdx() == atoms[whichEnd]->getIdx()) {
          BOOST_LOG(rdWarningLog)
              << "Wedge or hash bond found on atropisomer where not expected - atoms are: "
              << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
              << std::endl;
          return;
        } else {
          continue;  // wedge or hash bond affecting the OTHER atom =
                     // perhaps a chiral center
        }
      }
      auto ringCount = ri->numBondRings(bondToTry->getIdx());
      if (ringCount > bestRingCount) {
        continue;
      }

      else if (ringCount < bestRingCount) {
        bestBond = bondToTry;
        bestBondEnd = whichEnd;
        bestRingCount = ringCount;
        bestBondDir = getBondDirForAtropisomer3d(bondToTry, conf);
      } else {
        auto bondDir = getBondDirForAtropisomer3d(bondToTry, conf);
        if (bestBondDir == Bond::BondDir::NONE ||
            (bestBondDir == Bond::BondDir::BEGINDASH &&
             bondDir == Bond::BondDir::BEGINWEDGE)) {
          bestBond = bondToTry;
          bestBondEnd = whichEnd;
          bestRingCount = ringCount;
          bestBondDir = bondDir;
        }
      }
    }
  }

  if (bestBond != nullptr) {  // we found a good one

    // make sure the atoms on the bond are in the right order for the
    // wedge/hash the atom on the end of the main bond must be listed
    // first for the wedge/has bond

    if (bestBond->getBeginAtom() != atoms[bestBondEnd]) {
      bestBond->setEndAtom(bestBond->getBeginAtom());
      bestBond->setBeginAtom(atoms[bestBondEnd]);
    }
    bestBond->setBondDir(bestBondDir);
  } else {
    BOOST_LOG(rdWarningLog)
        << "Failed to find a good bond to set as UP or DOWN for an atropisomer - atoms are: "
        << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx() << std::endl;
    return;
  }

  return;
}

void WedgeBondsFromAtropisomers(const ROMol &mol, const Conformer *conf,
                                const INT_MAP_INT &wedgeBonds) {
  PRECONDITION(conf, "no conformer");
  PRECONDITION(&(conf->getOwningMol()) == &mol,
               "conformer does not belong to molecule");

  if (!mol.getRingInfo()->isInitialized()) {
    MolOps::fastFindRings(mol);
  }

  auto *ri = mol.getRingInfo();
  for (auto bond : mol.bonds()) {
    if (ri->numBondRings(bond->getIdx()) > 0) {
      continue;
    }

    auto bondStereo = bond->getStereo();

    if (bond->getBondType() != Bond::SINGLE ||
        (bondStereo != Bond::BondStereo::STEREOATROPCW &&
         bondStereo != Bond::BondStereo::STEREOATROPCCW) ||
        bond->getBeginAtom()->getTotalDegree() < 2 ||
        bond->getEndAtom()->getTotalDegree() < 2 ||
        bond->getBeginAtom()->getTotalDegree() > 3 ||
        bond->getEndAtom()->getTotalDegree() > 3) {
      continue;
    }

    if (conf->is3D()) {
      WedgeBondFromAtropisomerOneBond3d(bond, mol, conf, wedgeBonds);
    } else {
      WedgeBondFromAtropisomerOneBond2d(bond, mol, conf, wedgeBonds);
    }
  }
}

void ClearSingleBondDirFlags(ROMol &mol) {
  MolOps::clearSingleBondDirFlags(mol);
}

void DetectBondStereoChemistry(ROMol &mol, const Conformer *conf) {
  PRECONDITION(conf, "no conformer");
  PRECONDITION(&(conf->getOwningMol()) == &mol,
               "conformer does not belong to molecule");
  MolOps::detectBondStereochemistry(mol, conf->getId());
}

void reapplyMolBlockWedging(RWMol &mol) {
  if (mol.needsUpdatePropertyCache()) {
    mol.updatePropertyCache(false);
  }

  MolOps::clearAllBondDirFlags(mol);
  for (auto b : mol.bonds()) {
    int explicit_unknown_stereo = -1, molFileBondStereo = (-1);
    if (b->getBondType() == Bond::SINGLE &&
        b->getPropIfPresent<int>(common_properties::_UnknownStereo,
                                 explicit_unknown_stereo) &&
        explicit_unknown_stereo &&
        (b->getBondType() == Bond::SINGLE ||
         b->getBondType() == Bond::AROMATIC)) {
      b->setBondDir(Bond::UNKNOWN);
      continue;
    }
    if (b->getBondType() == Bond::DOUBLE &&
        b->getPropIfPresent<int>(common_properties::_MolFileBondStereo,
                                 molFileBondStereo) &&
        molFileBondStereo == 3) {
      b->setBondDir(Bond::BondDir::EITHERDOUBLE);
      // clearDoubleBondStereoNeighborsOfWiggle(b, mol);
      continue;
    }
    int bond_dir = -1;
    if (b->getPropIfPresent<int>(common_properties::_MolFileBondStereo,
                                 bond_dir) &&
        (b->getBondType() == Bond::SINGLE ||
         b->getBondType() == Bond::AROMATIC)) {
      if (bond_dir == 1) {
        b->setBondDir(Bond::BEGINWEDGE);
      } else if (bond_dir == 6) {
        b->setBondDir(Bond::BEGINDASH);
      }
      continue;
    }
    int cfg = -1;
    if (b->getPropIfPresent<int>(common_properties::_MolFileBondCfg, cfg)) {
      switch (cfg) {
        case 1:
          if ((b->getBondType() == Bond::SINGLE ||
               b->getBondType() == Bond::AROMATIC)) {
            b->setBondDir(Bond::BEGINWEDGE);
          }
          break;
        case 2:
          if ((b->getBondType() == Bond::SINGLE ||
               b->getBondType() == Bond::AROMATIC)) {
            b->setBondDir(Bond::UNKNOWN);
          } else if (b->getBondType() == Bond::DOUBLE) {
            b->setBondDir(Bond::EITHERDOUBLE);
            b->setStereo(Bond::STEREOANY);
          }
          break;
        case 3:
          if ((b->getBondType() == Bond::SINGLE ||
               b->getBondType() == Bond::AROMATIC)) {
            b->setBondDir(Bond::BEGINDASH);
          }
          break;
      }
    }
  }
}

void clearMolBlockWedgingInfo(ROMol &mol) {
  for (auto b : mol.bonds()) {
    b->clearProp(common_properties::_MolFileBondStereo);
    b->clearProp(common_properties::_MolFileBondCfg);
  }
}

void invertMolBlockWedgingInfo(ROMol &mol) {
  for (auto b : mol.bonds()) {
    int bond_dir = -1;
    if (b->getPropIfPresent<int>(common_properties::_MolFileBondStereo,
                                 bond_dir)) {
      if (bond_dir == 1) {
        b->setProp<int>(common_properties::_MolFileBondStereo, 6);
      } else if (bond_dir == 6) {
        b->setProp<int>(common_properties::_MolFileBondStereo, 1);
      }
    }
    int cfg = -1;
    if (b->getPropIfPresent<int>(common_properties::_MolFileBondCfg, cfg)) {
      if (cfg == 1) {
        b->setProp<int>(common_properties::_MolFileBondCfg, 3);
      } else if (cfg == 3) {
        b->setProp<int>(common_properties::_MolFileBondCfg, 1);
      }
    }
  }
}

void markUnspecifiedStereoAsUnknown(ROMol &mol, int confId) {
  INT_MAP_INT wedgeBonds = pickBondsToWedge(mol);
  const auto conf = mol.getConformer(confId);
  for (auto b : mol.bonds()) {
    if (b->getBondType() == Bond::DOUBLE) {
      int dirCode;
      bool reverse;
      GetMolFileBondStereoInfo(b, wedgeBonds, &conf, dirCode, reverse);
      if (dirCode == 3) {
        b->setStereo(Bond::STEREOANY);
      }
    }
  }
  static int noNbrs = 100;
  auto si = Chirality::findPotentialStereo(mol);
  if (si.size()) {
    std::pair<bool, INT_VECT> retVal = countChiralNbrs(mol, noNbrs);
    INT_VECT nChiralNbrs = retVal.second;
    for (auto i : si) {
      if (i.type == Chirality::StereoType::Atom_Tetrahedral &&
          i.specified == Chirality::StereoSpecified::Unspecified) {
        i.specified = Chirality::StereoSpecified::Unknown;
        auto atom = mol.getAtomWithIdx(i.centeredOn);
        INT_MAP_INT resSoFar;
        int bndIdx = pickBondToWedge(atom, mol, nChiralNbrs, resSoFar, noNbrs);
        auto bond = mol.getBondWithIdx(bndIdx);
        bond->setBondDir(Bond::UNKNOWN);
      }
    }
  }
}

}  // namespace RDKit
