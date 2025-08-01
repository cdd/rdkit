

if(NOT RDK_USE_BOOST_IOSTREAMS)
    find_package(ZLIB)
    set(STANDALONE_ZLIB_LIBRARY ${ZLIB_LIBRARIES})
    include_directories(${ZLIB_INCLUDE_DIRS})
    add_definitions("-DRDK_USE_STANDALONE_ZLIB")
endif(NOT RDK_USE_BOOST_IOSTREAMS)

if(RDK_BUILD_MAEPARSER_SUPPORT)
    set(MAEPARSER_LIB maeparser)
    set(maesupplier MaeMolSupplier.cpp MaeWriter.cpp)
    if ("${RDK_MAEPARSER_LIBS}" STREQUAL "")
        # if we're not building maeparser as part of the RDKit build), then find
        # the imported package
        find_package(ZLIB QUIET)
        find_package(${MAEPARSER_LIB})
    else()
        # if we built maeparser as part of the RDKit build, then we need to use
        # the headers from the source directory
        include_directories(${maeparser_INCLUDE_DIRS})
    endif()
endif()

rdkit_library(FileParsers
    CDXMLParser.cpp
    Mol2FileParser.cpp MolFileParser.cpp 
    SCSRMolFileParser.cpp
    MolSGroupParsing.cpp MolSGroupWriting.cpp
    MolFileStereochem.cpp MolFileWriter.cpp
    ForwardSDMolSupplier.cpp SDMolSupplier.cpp SDWriter.cpp
    SmilesMolSupplier.cpp 
    SmilesWriter.cpp
    TDTMolSupplier.cpp 
    TDTWriter.cpp
    TplFileParser.cpp TplFileWriter.cpp
    PDBParser.cpp PDBWriter.cpp 
    CMLWriter.cpp XYZFileWriter.cpp XYZFileParser.cpp
    ${maesupplier}
    ProximityBonds.cpp
    SequenceParsers.cpp SequenceWriters.cpp
    SVGParser.cpp PNGParser.cpp
    MultithreadedMolSupplier.cpp
    MultithreadedSmilesMolSupplier.cpp
    MultithreadedSDMolSupplier.cpp
    LINK_LIBRARIES GenericGroups Depictor SmilesParse ChemTransforms GraphMol SubstructMatch ${MAEPARSER_LIB} ${RDK_CHEMDRAW_LIBS} ${STANDALONE_ZLIB_LIBRARY}
)
if(STANDALONE_ZLIB_LIBRARY)
    target_include_directories(FileParsers PRIVATE ${ZLIB_INCLUDE_DIRS})
endif(STANDALONE_ZLIB_LIBRARY)
target_compile_definitions(FileParsers PRIVATE RDKIT_FILEPARSERS_BUILD)

rdkit_headers(CDXMLParser.h
    FileParsers.h
    FileParserUtils.h
    MolFileStereochem.h
    FileWriters.h
    MolSupplier.h MolSupplier.v1API.h
    MolWriters.h
    SequenceParsers.h SequenceWriters.h
    GeneralFileReader.h
    MultithreadedMolSupplier.h
    MultithreadedSmilesMolSupplier.h
    MultithreadedSDMolSupplier.h
    PNGParser.h
    DEST GraphMol/FileParsers)

rdkit_test(fileParsersTest1 test1.cpp
    LINK_LIBRARIES CIPLabeler FileParsers SubstructMatch)

rdkit_test(testAtropisomers testAtropisomers.cpp
    LINK_LIBRARIES FileParsers SubstructMatch)

rdkit_catch_test(testSquiggleDoubleBonds testSquiggleDoubleBonds.cpp
    LINK_LIBRARIES FileParsers SubstructMatch)

rdkit_test(testMolSupplier testMolSupplier.cpp
    LINK_LIBRARIES FileParsers RDStreams)

rdkit_test(testGeneralFileReader testGeneralFileReader.cpp
    LINK_LIBRARIES FileParsers RDStreams)

if(RDK_TEST_MULTITHREADED)
    rdkit_test(testMultithreadedMolSupplier testMultithreadedMolSupplier.cpp
        LINK_LIBRARIES FileParsers Fingerprints RDStreams)
endif(RDK_TEST_MULTITHREADED)

rdkit_test(testMolWriter testMolWriter.cpp LINK_LIBRARIES CIPLabeler FileParsers)

rdkit_test(testTplParser testTpls.cpp LINK_LIBRARIES FileParsers GraphMol)

rdkit_test(testMol2ToMol testMol2ToMol.cpp LINK_LIBRARIES FileParsers)

rdkit_test(testSequence testSequence.cpp LINK_LIBRARIES FileParsers)

rdkit_test(testExtendedStereoParsing testExtendedStereoParsing.cpp
    LINK_LIBRARIES FileParsers)

rdkit_catch_test(fileParsersCatchTest file_parsers_catch.cpp
    LINK_LIBRARIES FileParsers)

rdkit_catch_test(macromolsCatchTest macromols_catch.cpp
    LINK_LIBRARIES FileParsers)

rdkit_catch_test(testPropertyLists testPropertyLists.cpp
    LINK_LIBRARIES FileParsers)

rdkit_catch_test(XYZFileParserCatchTest XYZFileParserCatchTest.cpp
    LINK_LIBRARIES FileParsers)

rdkit_catch_test(cdxmlParserCatchTest cdxml_parser_catch.cpp
    LINK_LIBRARIES FileParsers)

rdkit_catch_test(molfileStereoCatchTest molfile_stereo_catch.cpp
    LINK_LIBRARIES FileParsers CIPLabeler Subgraphs)

rdkit_catch_test(connectTheDotsTest connectTheDots_catch.cpp
    LINK_LIBRARIES FileParsers)

rdkit_catch_test(v2MolSuppliers v2_suppliers_catch.cpp
    LINK_LIBRARIES FileParsers)

rdkit_catch_test(v2FileParsersCatchTest v2_file_parsers_catch.cpp
    LINK_LIBRARIES FileParsers)

rdkit_catch_test(atropisomersCatch atropisomers_catch.cpp
    LINK_LIBRARIES FileParsers)



if(RDK_TEST_MULTITHREADED AND RDK_BUILD_THREADSAFE_SSS)
rdkit_catch_test(multithreadedSupplierCatchTest multithreaded_supplier_catch.cpp
    LINK_LIBRARIES FileParsers)
endif()

