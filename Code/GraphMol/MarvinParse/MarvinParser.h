//
//  Copyright (C) 2002-2022 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_MARVINPARSER_H
#define RD_MARVINPARSER_H

#include <RDGeneral/types.h>
#include <GraphMol/RDKitBase.h>

#include <string>
#include <string_view>
#include <iostream>
#include <vector>
#include <exception>

#include <boost/shared_ptr.hpp>

namespace RDKit {


RDKIT_FILEPARSERS_EXPORT void *MrvFileParser(const std::string &fname, bool &isReaction);
RDKIT_FILEPARSERS_EXPORT void *MrvDataStreamParser(std::istream *inStream);
RDKIT_FILEPARSERS_EXPORT void *MrvDataStreamParser(std::istream &inStream);
RDKIT_FILEPARSERS_EXPORT void *MrvBlockParser(const std::string &molBlock);



}  // namespace RDKit

#endif
