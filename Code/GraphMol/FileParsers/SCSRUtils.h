//
//  Copyright (C) 2010-2025 Tad Hurst and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_SCSRUTILS_H
#define RD_SCSRUTILS_H

#include <string>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include "FileParsers.h"
#include <string_view>

namespace RDKit {
class RWMol;
class Conformer;

enum class SCSRTemplateNames {
  AsEntered,      //<! use the name of the temlate as entered in the SCSR Mol
  UseFirstName,   //<!Use the first name in the template
                  // def (For AA, the 3 letter code
  UseSecondName,  //<!use the second name in the tempate def (
                  // For AA, the 1 letter code)
  All             //<! use all names in the template def
};

enum class SCSRBaseHbondOptions {
  Ignore,     //<! Do not include base Hbonds in expanded output
  UseSapAll,  //<!use all hbonds defined in SAPs
              // can be more than one per base
  UseSapOne,  //<!use only one SAP hbond per base
              // If multiple SAPs are defined, use the first
              // even if it is not the best
              //(this just maintains the relationship between
              // the to base pairs)
  Auto        //<!For bases that are C,G,A,T,U,In (and
              // derivatives) use the standard Watson-Crick
              // Hbonding.  No SAPs need to be defined, and if
              // defined, they are ignored.
};

struct RDKIT_FILEPARSERS_EXPORT MolFromSCSRParams {
  bool includeLeavingGroups =
      true; /**< when true, leaving groups on atoms that are not exo-bonded are
                retained.  When false, no leaving groups are retained */
  SCSRTemplateNames scsrTemplateNames = SCSRTemplateNames::All;

  SCSRBaseHbondOptions scsrBaseHbondOptions = SCSRBaseHbondOptions::UseSapAll;
};
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RWMol> MolFromMolDataStream(
    std::istream &inStream, unsigned int &line,
    const RDKit::v2::FileParsers::MolFileParserParams &params =
        RDKit::v2::FileParsers::MolFileParserParams());
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RWMol> MolFromMolBlock(
    const std::string &molBlock,
    const RDKit::v2::FileParsers::MolFileParserParams &params =
        RDKit::v2::FileParsers::MolFileParserParams());
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RWMol> MolFromMolFile(
    const std::string &fName,
    const RDKit::v2::FileParsers::MolFileParserParams &params =
        RDKit::v2::FileParsers::MolFileParserParams());

RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RDKit::RWMol> MolFromSCSRDataStream(
    std::istream &inStream, unsigned int &line,
    const RDKit::v2::FileParsers::MolFileParserParams &molFileParserParams =
        RDKit::v2::FileParsers::MolFileParserParams(),
    const MolFromSCSRParams &molFromSCSRParams = MolFromSCSRParams());
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RDKit::RWMol> MolFromSCSRBlock(
    const std::string &molBlock,
    const RDKit::v2::FileParsers::MolFileParserParams &molFileParserParams =
        RDKit::v2::FileParsers::MolFileParserParams(),
    const MolFromSCSRParams &molFromSCSRParams = MolFromSCSRParams());
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RDKit::RWMol> MolFromSCSRFile(
    const std::string &fName,
    const RDKit::v2::FileParsers::MolFileParserParams &molFileParserParams =
        RDKit::v2::FileParsers::MolFileParserParams(),
    const MolFromSCSRParams &molFromSCSRParams = MolFromSCSRParams());

RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RDKit::SCSRMol>
SCSRMolFromSCSRDataStream(
    std::istream &inStream, unsigned int &line,
    const RDKit::v2::FileParsers::MolFileParserParams &params =
        RDKit::v2::FileParsers::MolFileParserParams());

RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RDKit::SCSRMol> SCSRMolFromSCSRBlock(
    const std::string &molBlock,
    const RDKit::v2::FileParsers::MolFileParserParams &params =
        RDKit::v2::FileParsers::MolFileParserParams());

std::unique_ptr<RDKit::SCSRMol> SCSRMolFromSCSRFile(
    const std::string &fName,
    const RDKit::v2::FileParsers::MolFileParserParams &params =
        RDKit::v2::FileParsers::MolFileParserParams());

}  // namespace RDKit

#endif
