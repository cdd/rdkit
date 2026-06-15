//
//  Copyright (C) 2025 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//


#include "Conversions.h"

#include <queue>
#include <unordered_map>
#include <memory>

#include <boost/filesystem/path.hpp>

#include <GraphMol/RWMol.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/MacroMolUtils.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include "MonomerLibrary.h"
#include "MonomerMol.h"

namespace RDKit
{


std::unique_ptr<RDKit::RWMol> toAtomistic(const MacroMol& monomer_mol) {
    RDKit::MolFromMacroMolParams molFromMacroMolParams;
    molFromMacroMolParams.includeLeavingGroups = true;
    molFromMacroMolParams.macroTemplateNames = MACROTemplateNames::UseFirstName;
    molFromMacroMolParams.outputSgroups = false;
    
    RDKit::v2::FileParsers::MolFileParserParams molFileParserParams;
    return RDKit::MolFromMacroMol(&monomer_mol, molFileParserParams, molFromMacroMolParams);
}


} // namespace RDKit