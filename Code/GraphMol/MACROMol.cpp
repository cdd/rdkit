#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include "FileParsers/FileParsers.h"
#include "FileParsers/FileParserUtils.h"

#include "MACROMol.h"
#include "Atom.h"

unsigned int RDKit::MACROMol::atomIdxToTemplateIdx(unsigned int atomIdx) {
  if (p_atomIdxToTemplateIdxIsStale) {
    // rebuild the map
    p_atomIdxToTemplateIdx.clear();

    for (auto atom : p_mol->atoms()) {
      std::string atomClass;
      std::string dummyLabel = "";

      if (!atom->getPropIfPresent(common_properties::dummyLabel, dummyLabel) ||
          dummyLabel != "" ||
          !atom->getPropIfPresent<std::string>(common_properties::molAtomClass,
                                               atomClass) ||
          atomClass == "") {
        p_atomIdxToTemplateIdx[atom->getIdx()] = UINT_MAX;
        continue;
      }

      for (unsigned int tIdx = 0; tIdx < getTemplateCount(); ++tIdx) {
        const ROMol *templateMol = getTemplate(tIdx);
        std::string templateAtomClass;
        std::vector<std::string> templateNames;
        templateMol->getPropIfPresent<std::string>(
            common_properties::molAtomClass, templateAtomClass);
        templateMol->getPropIfPresent<std::vector<std::string>>(
            common_properties::templateNames, templateNames);
        if (templateAtomClass == atomClass &&
            std::find(templateNames.begin(), templateNames.end(), dummyLabel) !=
                templateNames.end()) {
          p_atomIdxToTemplateIdx[atom->getIdx()] = tIdx;
        }
      }

      if (p_atomIdxToTemplateIdx.contains(atom->getIdx()) == false) {
        throw FileParseException("Template for macro atom not found");
      }
    }

    p_atomIdxToTemplateIdxIsStale = false;
  }
  return p_atomIdxToTemplateIdx[atomIdx];
}
