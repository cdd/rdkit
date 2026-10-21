//
//  Copyright (C) 2026 Tad Hurst, Schrödinger and other RDKit contributors
//
#include <GraphMol/FileParsers/MolToMacroMol.h>

#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/MacroMol.h>
#include <GraphMol/MacroMolTemplate.h>
#include <GraphMol/FileParsers/MolToMacroMol.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SubstanceGroup.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/Invariant.h>

#include <boost/functional/hash.hpp>

#include <algorithm>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace RDKit {
namespace {

constexpr const char *MACROMOL_GROUP = "MACROMOL";

struct AttachmentPointSpec {
  unsigned int sourceAtomIdx;
  unsigned int leavingAtomIdx;
  int attachmentIdx;
};

struct GroupSpec {
  const MacroMolTemplate *templ = nullptr;
  std::vector<unsigned int> sourceAtomIndices;
  std::vector<AttachmentPointSpec> attachmentPoints;

struct MainGroupQuery {
  std::unique_ptr<RWMol> query;
  std::vector<unsigned int> queryToTemplateAtom;
};

MainGroupQuery makeMainGroupQuery(const MacroMolTemplate &templ) {
  std::vector<unsigned int> mainAtoms(templ.getMainAtomIdxs().begin(),
                                      templ.getMainAtomIdxs().end());
  std::sort(mainAtoms.begin(), mainAtoms.end());
  const std::unordered_set<unsigned int> mainAtomSet(mainAtoms.begin(),
                                                     mainAtoms.end());

  auto query = std::make_unique<RWMol>(templ.getMol());
  query->beginBatchEdit();
  for (auto atom : query->atoms()) {
    if (!mainAtomSet.count(atom->getIdx())) {
      query->removeAtom(atom);
    }
  }
  query->commitBatchEdit();
  query->updatePropertyCache(false);
  return {std::move(query), std::move(mainAtoms)};
}

std::vector<MatchVectType> findTemplateMatches(const ROMol &mol,
                                               const RWMol &query) {
  SubstructMatchParameters params;
  params.recursionPossible = false;
  params.useChirality = true;
  params.useQueryQueryMatches = false;
  params.uniquify = false;
  params.maxMatches = 0;
  auto matches = SubstructMatch(mol, query, params);

  auto matchKey = [](const MatchVectType &match) {
    std::vector<int> sourceAtoms;
    std::vector<std::pair<int, int>> queryToSource;
    for (const auto &pair : match) {
      sourceAtoms.push_back(pair.second);
      queryToSource.push_back(pair);
    }
    std::sort(sourceAtoms.begin(), sourceAtoms.end());
    std::sort(queryToSource.begin(), queryToSource.end());
    return std::make_pair(sourceAtoms, queryToSource);
  };
  std::sort(matches.begin(), matches.end(),
            [&matchKey](const auto &lhs, const auto &rhs) {
              return matchKey(lhs) < matchKey(rhs);
            });
  return matches;
}


std::optional<GroupSpec> makeGroupForMatch(
    const MatchVectType &match, const ROMol &mol,
    const MacroMolTemplate &templ, const MainGroupQuery &mainQuery,
    const std::vector<bool> &claimed) {
  GroupSpec group;
  group.templ = &templ;
  std::unordered_map<unsigned int, unsigned int> templateToSource;
  std::map<unsigned int, unsigned int> molToQueryMap;
  for (const auto &[queryIdx, sourceIdxInt] : match) {
    const auto sourceIdx = static_cast<unsigned int>(sourceIdxInt);
    if (claimed[sourceIdx]) {
      return std::nullopt;
    }
    molToQueryMap[sourceIdxInt] = queryIdx;
    group.sourceAtomIndices.push_back(sourceIdx);
    templateToSource[mainQuery.queryToTemplateAtom[queryIdx]] = sourceIdx;
  }

  // for each atom in the match, check its bonds.  If it has a bond to
  // another atom in the hit,
  //  make sure that bond is the query.  If not, skip this match.

  // If the atom is an attachment point for the template, see if there is
  // a bond to an atom NOT in the match. If so, record the attachment
  // order for later.

  // if the atom is not an attachment point and the bond is not a hydrogen bond,  it has an external bond skip this
  // match.

  auto &leavingGroups = templ.getLeavingGroups();
  std::vector<bool> leavingGroupUsed(leavingGroups.size(), false);
  for (const auto &[atomIdx, queryAtomIdx] : molToQueryMap) {
    auto atom = mol.getAtomWithIdx(atomIdx);
    auto templateAtomIdx =  mainQuery.queryToTemplateAtom[queryAtomIdx];

    for (const auto bond : mol.atomBonds(atom)) {
      unsigned int nbrAtomIdx = bond->getOtherAtomIdx(atom->getIdx());
      if (molToQueryMap.contains(nbrAtomIdx)) {
        // make sure the bond is also in the query mol - if not skip this
        // match

        auto queryBond = mainQuery.query->getBondBetweenAtoms(
            molToQueryMap[bond->getBeginAtomIdx()],
            molToQueryMap[bond->getEndAtomIdx()]);

        if (!queryBond) {
          return std::nullopt;;
        }

        // bond to another atom is in the query - ok
        continue;
      } else {
        // bond to an atom NOT in the main SUP - could be a leaving group

        bool leavingGroupFound = false;
        for (unsigned int leavingGroupIndex = 0;
             leavingGroupIndex < leavingGroups.size(); ++leavingGroupIndex) {
          auto leavingGroup = leavingGroups[leavingGroupIndex];
    
          if (leavingGroup.attachAtomIdx != templateAtomIdx) {
            continue;
          }

          // don't reuse a leaving group.  If there are two or more on the same attachAtomIdx, 
          // if will be random how there are assigned.  And this could lead to undesireable pairings.
          // (like R1-R1 or Al-Al)
          //the example that causes this is the Phospate group that has two leaving groups on the same P
          // atom.
          //
          // The macroMol reoutines are not supposed to know about the meanings of the attachment types,
          // so it is up the the routines that use it (SCSR and HELM) to detect and fix these
          if (leavingGroupUsed[leavingGroupIndex]) {
            continue;
          }

          // see if the attachment point is an H or OH, and the nbr atom is
          // also an H or OH if so, this leaving group can be considered part
          // of the template

          leavingGroupFound = true;
          leavingGroupUsed[leavingGroupIndex] = true;

          auto nbrAtom = mol.getAtomWithIdx(nbrAtomIdx);
          const Atom *templateNbrAtom = nullptr;
          if (leavingGroup.leavingAtomIdx >= 0) {
            templateNbrAtom = templ.getMol().getAtomWithIdx(leavingGroup.leavingAtomIdx);
          }
          // if the template nbr atom is a single atom and so is the matched
          // atom, and they have the same atomic num and numHs, treat this as
          // part of the template

          if (leavingGroup.leavingAtomIdx >= 0 && nbrAtom->getDegree() == 1 &&
              templateNbrAtom->getDegree() == 1 &&
              templateNbrAtom->getAtomicNum() == nbrAtom->getAtomicNum() &&
              templateNbrAtom->getTotalNumHs() == nbrAtom->getTotalNumHs()) {
            // ok - treat this as part of the template

            group.sourceAtomIndices.push_back(nbrAtomIdx);
          } else {
            // save the attach point relative to the full atom molecule being
            // converted

            group.attachmentPoints.push_back({atomIdx, nbrAtomIdx, leavingGroup.attachPoint});
          }
          break;
        }

        // if a sap for the external bond is not found, and it is not a Hydrogen bond, this is not a valid hit to the template query

        if (!leavingGroupFound && bond->getBondType() != Bond::BondType::HYDROGEN) {
          // external bond that is not allowed - skip this match
          return std::nullopt;;
        }
      }
    }
  }

  std::sort(group.sourceAtomIndices.begin(), group.sourceAtomIndices.end());
  return group;
}

std::vector<GroupSpec> findGroups(const ROMol &mol,
                                  const MacroMolTemplateLibrary &templates) {
  std::vector<GroupSpec> groups;
  std::vector<bool> claimed(mol.getNumAtoms(), false);
  for (const auto *templ : templates.getTemplates()) {
    const auto mainQuery = makeMainGroupQuery(*templ);

    for (const auto &match : findTemplateMatches(mol, *mainQuery.query)) {
      auto group = makeGroupForMatch(match, mol, *templ, mainQuery, claimed);
      if (!group) {
        continue;
      }
      for (const auto atomIdx : group->sourceAtomIndices) {
        claimed[atomIdx] = true;
      }
      groups.push_back(std::move(*group));
    }
  }
  return groups;
}

void addGroupMetadata(RWMol &result, const GroupSpec &group) {
  SubstanceGroup sgroup(&result, "SUP");
  sgroup.setProp<bool>(MACROMOL_GROUP, true);
  sgroup.setProp<std::string>("CLASS",
                              monomerClassToString(group.templ->getMonomerClass()));
  sgroup.setProp<std::string>("LABEL", group.templ->getSymbol());
  sgroup.setAtoms(group.sourceAtomIndices);
  for (const auto &attach : group.attachmentPoints) {
    sgroup.addAttachPoint(attach.sourceAtomIdx, attach.leavingAtomIdx,
                          std::to_string(attach.attachmentIdx));
  }
  addSubstanceGroup(result, std::move(sgroup));
}

std::unique_ptr<RWMol> annotate(const ROMol &mol,
                                const std::vector<GroupSpec> &groups) {
  auto result = std::make_unique<RWMol>(mol);
  result->clearSubstanceGroups();

  for (const auto &group : groups) {
    addGroupMetadata(*result, group);
  }
  return result;
}

class AttachAtoms {
 public:
  unsigned int aIdx;
  unsigned int lvIdx;

  AttachAtoms(unsigned int aIdxInit, unsigned int lvIdxInit)
      : aIdx(aIdxInit), lvIdx(lvIdxInit) {}

  bool operator==(const AttachAtoms &other) const {
    return aIdx == other.aIdx && lvIdx == other.lvIdx;
  }

  struct Hash {
    size_t operator()(const AttachAtoms &aa) const noexcept {
      size_t h = std::hash<unsigned int>{}(aa.aIdx);
      boost::hash_combine(h, aa.lvIdx);
      return h;
    }
  };
};

struct CollapseGroup {
  std::vector<unsigned int> atoms;
  std::unordered_map<AttachAtoms, int, AttachAtoms::Hash> attachmentByAtoms;
  const MacroMolTemplate *templ = nullptr;
  unsigned int outputAtom = std::numeric_limits<unsigned int>::max();
};

const MacroMolTemplate *findTemplate(const MacroMolTemplateLibrary &templates,
                                     const SubstanceGroup &sgroup) {
  const auto monomerClass = monomerClassFromString(
      sgroup.getProp<std::string>("CLASS"));
  const auto label = sgroup.getProp<std::string>("LABEL");
  auto *templ = templates.getBySymbol(monomerClass, label);
  if (!templ) {
    templ = templates.getByName(monomerClass, label);
  }
  return templ;
}

void copyRegularBond(MacroMol &result, const Bond &source,
                     unsigned int beginIdx, unsigned int endIdx,
                     const std::vector<unsigned int> &outputForAtom) {
  std::unique_ptr<Bond> bond(source.copy());
  bond->setBeginAtomIdx(beginIdx);
  bond->setEndAtomIdx(endIdx);
  for (auto &stereoAtom : bond->getStereoAtoms()) {
    if (stereoAtom >= 0 &&
        static_cast<unsigned int>(stereoAtom) < outputForAtom.size()) {
      stereoAtom = static_cast<int>(
          outputForAtom[static_cast<unsigned int>(stereoAtom)]);
    }
  }
  result.addBond(bond.get(), true);
  bond.release();
}

}  // namespace

std::unique_ptr<RWMol> MolToMacroMol(
  const ROMol &mol, const MacroMolTemplateLibrary &templates) {
  return annotate(mol, findGroups(mol, templates));
}

std::unique_ptr<MacroMol> CollapseMacroMol(
      const ROMol &groupedMol, const MacroMolTemplateLibrary &templates) {
  auto result = std::make_unique<MacroMol>();
  std::vector<int> groupForAtom(groupedMol.getNumAtoms(), -1);
  std::vector<CollapseGroup> groups;
  std::unordered_set<const MacroMolTemplate *> usedTemplates;

  for (const auto &sgroup : getSubstanceGroups(groupedMol)) {
    bool isMacroGroup = false;
    if (!sgroup.getPropIfPresent<bool>(MACROMOL_GROUP, isMacroGroup) ||
        !isMacroGroup) {
      continue;
    }
    CollapseGroup group;
    group.atoms = sgroup.getAtoms();
    group.templ = findTemplate(templates, sgroup);

    if (!group.templ || group.atoms.empty()) {
      throw ValueErrorException("invalid MacroMol SubstanceGroup");
    }
    usedTemplates.insert(group.templ);

    for (const auto atomIdx : group.atoms) {
      if (atomIdx >= groupedMol.getNumAtoms() || groupForAtom[atomIdx] != -1) {
        throw ValueErrorException("overlapping or invalid MacroMol groups");
      }
      groupForAtom[atomIdx] = static_cast<int>(groups.size());
    }
    for (const auto &attach : sgroup.getAttachPoints()) {
      if (attach.aIdx >= groupedMol.getNumAtoms() ||
          groupForAtom[attach.aIdx] != static_cast<int>(groups.size()) ||
          groupForAtom[attach.lvIdx] == static_cast<int>(groups.size())) {
        throw ValueErrorException("invalid MacroMol attachment point");
      }
      int attachmentIdx = 0;
      try {
        attachmentIdx = std::stoi(attach.id);
      } catch (const std::exception &) {
        throw ValueErrorException("invalid MacroMol attachment point ID");
      }
      group.attachmentByAtoms.emplace(AttachAtoms(attach.aIdx, attach.lvIdx), attachmentIdx);
    }
    groups.push_back(std::move(group));
  }

  std::vector<unsigned int> outputForAtom(groupedMol.getNumAtoms());
  for (const auto *atom : groupedMol.atoms()) {
    const auto atomIdx = atom->getIdx();
    const auto groupIdx = groupForAtom[atomIdx];
    if (groupIdx >= 0) {
      auto &group = groups[static_cast<unsigned int>(groupIdx)];
      if (group.outputAtom == std::numeric_limits<unsigned int>::max()) {
        group.outputAtom = result->addMacroAtom(group.templ->getSymbol(),
                                                group.templ->getMonomerClass());
      }
      outputForAtom[atomIdx] = group.outputAtom;
    } else {
      auto copy = std::unique_ptr<Atom>(atom->copy());
      outputForAtom[atomIdx] = result->addAtom(copy.get(), true, true);
      copy.release();
    }
  }

  for (const auto *bond : groupedMol.bonds()) {
    const auto begin = bond->getBeginAtomIdx();
    const auto end = bond->getEndAtomIdx();
    const auto beginGroup = groupForAtom[begin];
    const auto endGroup = groupForAtom[end];
    if (beginGroup >= 0 && beginGroup == endGroup) {
      continue;  // bond is inside the same group - ignore it
    }

    int beginAttach = -1;
    int endAttach = -1;

    if (beginGroup >= 0 ) {
      const auto &attachmentByAtoms = groups[beginGroup].attachmentByAtoms;

      auto attachment = attachmentByAtoms.find(AttachAtoms(begin,end));
      if (attachment != attachmentByAtoms.end()) {
        beginAttach = attachment->second;
      }
    }

    if (endGroup >= 0 ) {
      const auto &attachmentByAtoms = groups[endGroup].attachmentByAtoms;

      auto attachment = attachmentByAtoms.find(AttachAtoms(end, begin));
      if (attachment != attachmentByAtoms.end()) {
        endAttach = attachment->second;
      }
    }

    if (bond->getBondType() != Bond::BondType::HYDROGEN && ((beginGroup >= 0 && beginAttach == -1 ) ||
          (endGroup >= 0 && endAttach == -1))) {
      throw ValueErrorException(
          "MacroMol group has an unannotated crossing bond");
    }
                          

    const auto bondType = bond->getBondType();
    if (beginAttach >= 0 && endAttach >= 0) {
      result->addMacroBond(outputForAtom[begin], outputForAtom[end],
                            beginAttach, endAttach,
                            bondType);
    } else if (beginAttach >= 0) {
      result->addMacroAtomToAtomBond(outputForAtom[begin], outputForAtom[end],
                                      beginAttach, bondType);
    } else if (endAttach >= 0) {
      result->addAtomToMacroAtomBond(outputForAtom[begin], outputForAtom[end],
                                      endAttach, bondType);
    } else {
      if (result->getBondBetweenAtoms( outputForAtom[begin], outputForAtom[end]) == nullptr) {
        copyRegularBond(*result, *bond, outputForAtom[begin],
                      outputForAtom[end], outputForAtom);
      } else if (bondType != Bond::BondType::HYDROGEN) {
        throw ValueErrorException(
          "Attempt to add the same bond twice");
      }
    }
  }

  // add the used templates to the local library

  for (auto macroTemplate : usedTemplates) {
    auto newTemplate = std::unique_ptr<MacroMolTemplate>(new MacroMolTemplate(*macroTemplate));
    result->addLocalTemplate(std::move(newTemplate));
  }
  return result;
}

}  // namespace RDKit
