#include <algorithm>
#include <utility>
#include "StereoGroup.h"
#include "Atom.h"
#include "Bond.h"

namespace RDKit {

StereoGroup::StereoGroup(StereoGroupType grouptype, std::vector<Atom *> &&atoms,
                         std::vector<Bond *> &&bonds)
    : d_grouptype(grouptype), d_atoms(atoms), d_bonds(bonds) {}
StereoGroup::StereoGroup(StereoGroupType grouptype,
                         const std::vector<Atom *> &atoms,
                         std::vector<Bond *> &bonds)
    : d_grouptype(grouptype),
      d_atoms(std::move(atoms)),
      d_bonds(std::move(bonds)) {}

StereoGroupType StereoGroup::getGroupType() const { return d_grouptype; }

const std::vector<Atom *> &StereoGroup::getAtoms() const { return d_atoms; }
const std::vector<Bond *> &StereoGroup::getBonds() const { return d_bonds; }

void removeGroupsWithAtom(const Atom *atom, std::vector<StereoGroup> &groups) {
  auto containsAtom = [atom](const StereoGroup &group) {
    return std::find(group.getAtoms().cbegin(), group.getAtoms().cend(),
                     atom) != group.getAtoms().cend();
  };
  groups.erase(std::remove_if(groups.begin(), groups.end(), containsAtom),
               groups.end());
}

void removeGroupsWithBond(const Bond *bond, std::vector<StereoGroup> &groups) {
  auto containsBond = [bond](const StereoGroup &group) {
    return std::find(group.getBonds().cbegin(), group.getBonds().cend(),
                     bond) != group.getBonds().cend();
  };
  groups.erase(std::remove_if(groups.begin(), groups.end(), containsBond),
               groups.end());
}

void removeGroupsWithAtoms(const std::vector<Atom *> &atoms,
                           std::vector<StereoGroup> &groups) {
  auto containsAnyAtom = [&atoms](const StereoGroup &group) {
    for (auto atom : atoms) {
      if (std::find(group.getAtoms().cbegin(), group.getAtoms().cend(), atom) !=
          group.getAtoms().cend()) {
        return true;
      }
    }
    return false;
  };
  groups.erase(std::remove_if(groups.begin(), groups.end(), containsAnyAtom),
               groups.end());
}

void removeGroupsWithBonds(const std::vector<Bond *> &bonds,
                           std::vector<StereoGroup> &groups) {
  auto containsAnyBond = [&bonds](const StereoGroup &group) {
    for (auto bond : bonds) {
      if (std::find(group.getBonds().cbegin(), group.getBonds().cend(), bond) !=
          group.getBonds().cend()) {
        return true;
      }
    }
    return false;
  };
  groups.erase(std::remove_if(groups.begin(), groups.end(), containsAnyBond),
               groups.end());
}

}  // namespace RDKit

std::ostream &operator<<(std::ostream &target, const RDKit::StereoGroup &stg) {
  switch (stg.getGroupType()) {
    case RDKit::StereoGroupType::STEREO_ABSOLUTE:
      target << "ABS";
      break;
    case RDKit::StereoGroupType::STEREO_OR:
      target << "OR";
      break;
    case RDKit::StereoGroupType::STEREO_AND:
      target << "AND";
      break;
  }
  target << " Atoms: { ";
  for (auto atom : stg.getAtoms()) {
    target << atom->getIdx() << ' ';
  }
  if (stg.getBonds().size() > 0) {
    target << " Bonds: { ";
    for (auto bond : stg.getBonds()) {
      target << bond->getIdx() << ' ';
    }
    target << '}';
  }
  target << '}';

  return target;
}
