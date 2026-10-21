//
//  Copyright (C) 2026 Tad Hurst, Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RD_MOLFROMMACROMOL_H
#define RD_MOLFROMMACROMOL_H

#include <GraphMol/MacroMol.h>
#include <GraphMol/MacroMolTemplate.h>
#include <GraphMol/MacroMolTemplateLibrary.h>
#include <RDGeneral/export.h>

#include <memory>

namespace RDKit {

struct RDKIT_FILEPARSERS_EXPORT MolFromMacroMolParams {
    bool includeLeavingGroups =
        true; /**< when true, leaving groups on atoms that are not exo-bonded
                are retained.  When false, no leaving groups are retained */
    bool outputSgroups = true;
};


RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RWMol> MolFromMacroMol(
    const MacroMol &macroMol, const RDKit::MolFromMacroMolParams &molFromMacroMolParams=RDKit::MolFromMacroMolParams());

RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RWMol> MolFromMacroMol(
    const MacroMol &macroMol, const MacroMolTemplateLibrary &templates, const RDKit::MolFromMacroMolParams &molFromMacroMolParams=RDKit::MolFromMacroMolParams());

}  // namespace RDKit

#endif
