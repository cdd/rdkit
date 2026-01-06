//
//  Copyright (C) 2024 Tad Hurst and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RD_MACROMOL_H
#define RD_MACROMOL_H

#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include "FileParsers/FileParsers.h"
#include "FileParsers/FileParserUtils.h"

namespace RDKit {

enum class MACROTemplateNames {
  AsEntered,      //<! use the name of the temlate as entered in the MACRO Mol
  UseFirstName,   //<!Use the first name in the template
                  // def (For AA, the 3 letter code
  UseSecondName,  //<!use the second name in the tempate def (
                  // For AA, the 1 letter code)
  All             //<! use all names in the template def
};

struct RDKIT_FILEPARSERS_EXPORT MolFromMACROMolParams {
  bool includeLeavingGroups =
      true; /**< when true, leaving groups on atoms that are not exo-bonded
               are retained.  When false, no leaving groups are retained */
  MACROTemplateNames macroTemplateNames = MACROTemplateNames::All;
};

class RDKIT_GRAPHMOL_EXPORT MACROMol {
 private:
  std::unique_ptr<RWMol> p_mol;
  std::vector<std::unique_ptr<ROMol>> p_templates;

 public:
  MACROMol() {};
  MACROMol(const MACROMol &other) = delete;
  MACROMol(MACROMol &&other) noexcept = delete;
  MACROMol &operator=(MACROMol &&other) noexcept = delete;

  MACROMol &operator=(const MACROMol &) = delete;  // disable assignment
  ~MACROMol() {}

  void addTemplate(std::unique_ptr<ROMol> templateMol) {
    PRECONDITION(templateMol, "bad template molecule");
    p_templates.push_back(std::move(templateMol));
  }

  unsigned int getTemplateCount() const { return p_templates.size(); }

  ROMol *getTemplate(unsigned int index) { return p_templates[index].get(); };
  const ROMol *getTemplate(unsigned int index) const {
    return p_templates[index].get();
  };

  const ROMol *getMol() const { return p_mol.get(); }

  RWMol *getMol() { return p_mol.get(); }

  void setMol(std::unique_ptr<RWMol> mol) {
    PRECONDITION(mol, "bad molecule");
    p_mol = std::move(mol);
  }

  static std::unique_ptr<RDKit::RWMol> MolFromMacroRMol(
      MACROMol *scsrMol,
      const RDKit::v2::FileParsers::MolFileParserParams &molFileParserParams,
      const MolFromMACROMolParams &molFromMACROMolParams);
};

typedef boost::shared_ptr<MACROMol> MACROMol_SPTR;

}  // namespace RDKit

#endif
