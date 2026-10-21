//
//  Copyright (C) 2026 Tad Hurst, Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/FileParsers/MolFromMacroMol.h>

#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SubstanceGroup.h>
#include <RDGeneral/FileParseException.h>
#include <GraphMol/Conformer.h>

#include <map>
#include <memory>
#include <set>
#include <string>

namespace RDKit {
namespace {

struct MacroAtomDetails {
  MonomerClass monomerClass = MonomerClass::Other;
  std::string label;
};

struct ExpandedMacroAtom {
  std::map<unsigned int, unsigned int> templateToResultAtom;
  std::map<int, unsigned int> attachPointToTemplateAtom;
  std::set<int> usedAttachPoints;
};

struct BuildState {
  RWMol &result;
  std::map<unsigned int, unsigned int> atomToResultAtom;
  std::map<unsigned int, ExpandedMacroAtom> macroAtoms;
  std::map<const MacroMolTemplate *, RDGeom::Point3D> templateCentroids;

  // map from final atom to an Sgroup Id
  std::map<const unsigned int,unsigned int> endAtomToSgroup; 
  double maxSize;
  const Conformer *conf;
  std::unique_ptr<Conformer> newConf = nullptr;

};

bool getMacroAtomDetails(const Atom &atom, MacroAtomDetails &info) {
  const auto *macroInfo = atom.getMacroAtomInfo();
  if (!macroInfo) {
    return false;
  }
  info.monomerClass = macroInfo->getMonomerClass();
  info.label = macroInfo->getSymbol();
  return !info.label.empty();
}

const MacroMolTemplate *findTemplate(
    const MacroMolTemplateLibrary &templates, const MacroAtomDetails &info) {
  const auto *templ = templates.getBySymbol(info.monomerClass, info.label);
  if (!templ) {
    templ = templates.getByName(info.monomerClass, info.label);
  }
  return templ;
}

unsigned int addBondCopy(RWMol &result, const Bond &oldBond,
                 unsigned int beginAtomIdx, unsigned int endAtomIdx,
                 Bond::BondType bondType) {
  auto *newBond = new Bond(bondType);
  newBond->setBeginAtomIdx(beginAtomIdx);
  newBond->setEndAtomIdx(endAtomIdx);
  newBond->setBondDir(oldBond.getBondDir());
  newBond->updateProps(oldBond, false);
  return result.addBond(newBond, true);
}

unsigned int addBondCopy(RWMol &result, const Bond &oldBond,
                 unsigned int beginAtomIdx, unsigned int endAtomIdx) {
  return addBondCopy(result, oldBond, beginAtomIdx, endAtomIdx,
              oldBond.getBondType());
}

void copyRegularAtom(BuildState &state, const Atom &atom) {
  auto *newAtom = new Atom(atom);
  state.atomToResultAtom[atom.getIdx()] =
      state.result.addAtom(newAtom, true, true);
  if (state.conf != nullptr) {
    state.newConf->setAtomPos(newAtom->getIdx(),
            state.conf->getAtomPos(atom.getIdx()) * state.maxSize);
  }

}

std::set<int> findUsedAttachPoints(const MacroMol &macroMol,
                                   const Atom &macroAtom) {
  std::set<int> usedAttachPoints;
  for (const auto *macroBond : macroMol.bonds()) {
    const auto *macroBondInfo = macroBond->getMacroBondInfo();
    if (!macroBondInfo) {
      continue;
    }

    for (const auto &bondInfo : macroBondInfo->getBonds()) {
      int attachPoint = -1;
      if (macroBond->getBeginAtomIdx() == macroAtom.getIdx()) {
        attachPoint = bondInfo.beginAttachPt;
      } else if (macroBond->getEndAtomIdx() == macroAtom.getIdx()) {
        attachPoint = bondInfo.endAttachPt;
      }
      if (attachPoint >= 0) {
        usedAttachPoints.insert(attachPoint);
      }
    }
  }
  return usedAttachPoints;
}

std::set<unsigned int> findLeavingAtomsToSkip(
    const MacroMolTemplate &templ, const std::set<int> &usedAttachPoints,
    const MolFromMacroMolParams &molFromMacroMolParams) {
  std::set<unsigned int> atomsToSkip;
  for (const auto &leavingGroup : templ.getLeavingGroups()) {
    if (!molFromMacroMolParams.includeLeavingGroups ||
      usedAttachPoints.find(leavingGroup.attachPoint) !=
        usedAttachPoints.end()) {
      atomsToSkip.insert(leavingGroup.atomIdxs.begin(),
                         leavingGroup.atomIdxs.end());
    }
  }
  return atomsToSkip;
}

void rememberTemplateAttachPoints(ExpandedMacroAtom &expanded,
                                  const MacroMolTemplate &templ) {
  for (const auto &leavingGroup : templ.getLeavingGroups()) {
    expanded.attachPointToTemplateAtom[leavingGroup.attachPoint] =
        leavingGroup.attachAtomIdx;
  }
}

void copyTemplateAtoms(BuildState &state, ExpandedMacroAtom &expanded,
                       const MacroMolTemplate &templ,
                       const std::set<unsigned int> &atomsToSkip,
                      RDGeom::Point3D  coordOffset) {

  if (state.newConf) {
    state.newConf->resize(state.newConf->getNumAtoms() + templ.getMol().getNumAtoms() - atomsToSkip.size());
  }

  for (const auto *oldAtom : templ.getMol().atoms()) {
    if (atomsToSkip.find(oldAtom->getIdx()) != atomsToSkip.end()) {
      continue;
    }

    auto *newAtom = new Atom(*oldAtom);
    expanded.templateToResultAtom[oldAtom->getIdx()] =
        state.result.addAtom(newAtom, true, true);
    if (state.newConf) {
      state.newConf->setAtomPos(
          newAtom->getIdx(),
          coordOffset +
              templ.getMol().getConformer().getAtomPos(oldAtom->getIdx()));
    }
  }
}

void copyTemplateBonds(BuildState &state, const ExpandedMacroAtom &expanded,
                       const MacroMolTemplate &templ) {
  for (const auto *oldBond : templ.getMol().bonds()) {
    const auto beginAtom =
        expanded.templateToResultAtom.find(oldBond->getBeginAtomIdx());
    const auto endAtom =
        expanded.templateToResultAtom.find(oldBond->getEndAtomIdx());
    if (beginAtom == expanded.templateToResultAtom.end() ||
        endAtom == expanded.templateToResultAtom.end()) {
      continue;
    }
    addBondCopy(state.result, *oldBond, beginAtom->second, endAtom->second);
  }
}

void copyMainSgroup(BuildState &state, const ExpandedMacroAtom &expanded,
                    const MacroMolTemplate &templ) {
  SubstanceGroup mainSgroup(&state.result, "SUP");
  mainSgroup.setProp(
      "CLASS", std::string(monomerClassToString(templ.getMonomerClass())));
  mainSgroup.setProp("LABEL", templ.getName());

  std::vector<unsigned int> resultMainAtoms;
  resultMainAtoms.reserve(templ.getMainAtomIdxs().size());
  for (const auto templateAtomIdx : templ.getMainAtomIdxs()) {
    auto resultAtomIdx = expanded.templateToResultAtom.at(templateAtomIdx);
    resultMainAtoms.push_back(resultAtomIdx);
    state.endAtomToSgroup[resultAtomIdx] = getSubstanceGroups(state.result).size();
  }

  for (const auto &leavingGroup : templ.getLeavingGroups()) {
    const auto leavingAtom =
        expanded.templateToResultAtom.find(leavingGroup.leavingAtomIdx);
    if (leavingAtom != expanded.templateToResultAtom.end()) {
      resultMainAtoms.push_back(leavingAtom->second);
      state.endAtomToSgroup[leavingAtom->second] = getSubstanceGroups(state.result).size();
    }
                              
  }


  mainSgroup.setAtoms(std::move(resultMainAtoms));

  addSubstanceGroup(state.result, std::move(mainSgroup));
}

void copyMacroAtom(BuildState &state, const MacroMol &macroMol,
                   const MacroMolTemplateLibrary &templates, const Atom &atom,
                   const MacroAtomDetails &info,
                   const MolFromMacroMolParams &molFromMacroMolParams) {
  const auto *templ = findTemplate(templates, info);
  if (!templ) {
    throw FileParseException("No template found for macro atom " +
                             std::to_string(atom.getIdx()));
  }

  ExpandedMacroAtom expanded;
  expanded.usedAttachPoints = findUsedAttachPoints(macroMol, atom);
  const auto atomsToSkip =
      findLeavingAtomsToSkip(*templ, expanded.usedAttachPoints, molFromMacroMolParams);

  RDGeom::Point3D coordOffset;
  if (state.conf) {
    coordOffset = (state.conf->getAtomPos(atom.getIdx()) * state.maxSize) -
                      state.templateCentroids[templ];
  }

  rememberTemplateAttachPoints(expanded, *templ);
  copyTemplateAtoms(state, expanded, *templ, atomsToSkip, coordOffset);
  copyTemplateBonds(state, expanded, *templ);

  if (molFromMacroMolParams.outputSgroups) {
    copyMainSgroup(state, expanded, *templ);
  }

  state.macroAtoms[atom.getIdx()] = expanded;
}

void copyAtoms(BuildState &state, const MacroMol &macroMol,
               const MacroMolTemplateLibrary &templates,
               const MolFromMacroMolParams &molFromMacroMolParams) {
  for (const auto *atom : macroMol.atoms()) {
    MacroAtomDetails info;
    if (getMacroAtomDetails(*atom, info)) {
      copyMacroAtom(state, macroMol, templates, *atom, info,
               molFromMacroMolParams);
    } else {
      copyRegularAtom(state, *atom);
    }
  }
}

unsigned int getResultAtomForBond(const BuildState &state, const Atom &atom,
                                  int attachPointId) {
  MacroAtomDetails info;
  if (!getMacroAtomDetails(atom, info)) {
    return state.atomToResultAtom.at(atom.getIdx());
  }

  const auto &expanded = state.macroAtoms.at(atom.getIdx());
  const auto attachPoint =
      expanded.attachPointToTemplateAtom.find(attachPointId);
  if (attachPoint == expanded.attachPointToTemplateAtom.end()) {
    throw FileParseException("Macro atom bond is missing an attachment point");
  }

  return expanded.templateToResultAtom.at(attachPoint->second);
}

void copyMacroMolBonds(BuildState &state, const MacroMol &macroMol, bool outputSgroups) {
  for (const auto *oldBond : macroMol.bonds()) {
    const auto *macroBondInfo = oldBond->getMacroBondInfo();
    if (!macroBondInfo) {
      const auto beginAtomIdx =
          getResultAtomForBond(state, *oldBond->getBeginAtom(), -1);
      const auto endAtomIdx =
          getResultAtomForBond(state, *oldBond->getEndAtom(), -1);
      addBondCopy(state.result, *oldBond, beginAtomIdx, endAtomIdx);
      continue;
    }

    for (const auto &bondInfo : macroBondInfo->getBonds()) {
      const auto beginAtomIdx = getResultAtomForBond(
          state, *oldBond->getBeginAtom(), bondInfo.beginAttachPt);
      const auto endAtomIdx = getResultAtomForBond(
          state, *oldBond->getEndAtom(), bondInfo.endAttachPt);
          
      auto newBondIndex = addBondCopy(state.result, *oldBond, beginAtomIdx, endAtomIdx,
                  static_cast<Bond::BondType>(bondInfo.bondType)) - 1;

      if (outputSgroups) {
        if (state.endAtomToSgroup.contains(beginAtomIdx)){
          auto sgroupId =  state.endAtomToSgroup[beginAtomIdx];
          getSubstanceGroups(state.result)[sgroupId].addBondWithIdx(newBondIndex);
        }

        if (state.endAtomToSgroup.contains(endAtomIdx)) {

          auto sgroupId =  state.endAtomToSgroup[endAtomIdx];
          getSubstanceGroups(state.result)[sgroupId].addBondWithIdx(newBondIndex);
        }
      }



    }
  }
}

}  // namespace

std::unique_ptr<RWMol> MolFromMacroMol(
    const MacroMol &macroMol, const MacroMolTemplateLibrary &templates, 
    const RDKit::MolFromMacroMolParams &molFromMacroMolParams) {
  auto result = std::make_unique<RWMol>();
  BuildState state{*result, {}, {}, {}, {}, 0.0, nullptr, nullptr};

   // first get some information from the templates to be used when
  // creating the coords for the new atoms. this is a dirty approach
  // that simply expands the orginal macro atom coords to be big
  // enough to hold any expanded macro atom. No attempt is made to
  // make this look nice, or to avoid overlaps.

  if (macroMol.getNumConformers() != 0) {
    state.conf = &macroMol.getConformer(0);
    state.newConf.reset(new Conformer(macroMol.getNumAtoms()));
    state.newConf->set3D(state.conf->is3D());

    // loop over all main mol atoms and make a list of the templates used
    // external libs could have thousands of tempaltes, most not used by the 
    // current macro mol

    std::set<const MacroMolTemplate *> templatesInUse;
    for (unsigned int atomIdx = 0 ; atomIdx != macroMol.getNumAtoms(); ++atomIdx) {
      auto atom = macroMol.getAtomWithIdx(atomIdx);
      auto macroAtomInfo = atom->getMacroAtomInfo();
      if (macroAtomInfo == nullptr) {
        continue;
      }

      auto monomerClass = macroAtomInfo->getMonomerClass();
      auto templateSymbol = macroAtomInfo->getSymbol();

      auto templatePtr = templates.getBySymbol(monomerClass, templateSymbol);

      if (templatePtr != nullptr && !templatesInUse.contains(templatePtr)) {
        templatesInUse.insert(templatePtr);
      }
    }

    for (auto templateMol : templatesInUse) {
      RDGeom::Point3D sumOfCoords;
      const RDKit::Conformer *templateConf = nullptr;
      auto confCount = templateMol->getMol().getNumConformers();
      if (confCount == 0) {
        state.conf = nullptr;   // all templates must have coords or we punt
        break;
      }
      templateConf = &templateMol->getMol().getConformer(0);
      RDGeom::Point3D maxCoord = templateConf->getAtomPos(0);
      RDGeom::Point3D minCoord = maxCoord;
      for (unsigned int atomIdx = 0; atomIdx < templateMol->getMol().getNumAtoms();
            ++atomIdx) {
        auto atomCoord = templateConf->getAtomPos(atomIdx);
        sumOfCoords += atomCoord;

        if (atomCoord.x > maxCoord.x) {
          maxCoord.x = atomCoord.x;
        }
        if (atomCoord.y > maxCoord.y) {
          maxCoord.y = atomCoord.y;
        }
        if (atomCoord.z > maxCoord.z) {
          maxCoord.z = atomCoord.z;
        }
        if (atomCoord.x < minCoord.x) {
          minCoord.x = atomCoord.x;
        }
        if (atomCoord.y < minCoord.y) {
          minCoord.y = atomCoord.y;
        }
        if (atomCoord.z < minCoord.z) {
          minCoord.z = atomCoord.z;
        }
      }
      state.templateCentroids[templateMol] = (sumOfCoords / templateMol->getMol().getNumAtoms());
      if (maxCoord.x - minCoord.x > state.maxSize) {
        state.maxSize = maxCoord.x - minCoord.x;
      }
      if (maxCoord.y - minCoord.y > state.maxSize) {
        state.maxSize = maxCoord.y - minCoord.y;
      }
      if (maxCoord.z - minCoord.z > state.maxSize) {
        state.maxSize = maxCoord.z - minCoord.z;
      }
    }
  }
  
  copyAtoms(state, macroMol, templates, molFromMacroMolParams);
  copyMacroMolBonds(state, macroMol, molFromMacroMolParams.outputSgroups);

  if (state.conf != nullptr) {
    state.newConf->resize(state.result.getNumAtoms());
    state.result.addConformer(state.newConf.release());
  }

  result->updatePropertyCache(false);
  return result;
}

std::unique_ptr<RWMol> MolFromMacroMol(const MacroMol &macroMol, 
  const RDKit::MolFromMacroMolParams &molFromMacroMolParams) {
  return MolFromMacroMol(macroMol, getGlobalMacroMolTemplateLibrary(), molFromMacroMolParams);
}

}  // namespace RDKit
