CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/versions.cpp.cmake
        ${CMAKE_CURRENT_SOURCE_DIR}/versions.cpp)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/versions.h.cmake
        ${CMAKE_CURRENT_SOURCE_DIR}/versions.h)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/RDConfig.h.cmake
        ${CMAKE_CURRENT_SOURCE_DIR}/RDConfig.h)

FetchContent_MakeAvailable(better_enums)
if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/enum.h)
    file(COPY ${better_enums_SOURCE_DIR}/enum.h DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
endif (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/enum.h)

rdkit_library(RDGeneral
        Invariant.cpp types.cpp utils.cpp RDGeneralExceptions.cpp RDLog.cpp
        LocaleSwitcher.cpp versions.cpp SHARED)
target_compile_definitions(RDGeneral PRIVATE RDKIT_RDGENERAL_BUILD)

if (RDK_USE_BOOST_STACKTRACE AND UNIX AND NOT APPLE)
    target_compile_definitions(RDGeneral PRIVATE RDK_USE_BOOST_STACKTRACE)
    set(EXTRA_STACKTRACE_LIBS dl)
endif ()

target_link_libraries(RDGeneral PUBLIC ${RDKit_THREAD_LIBS} ${EXTRA_STACKTRACE_LIBS})


rdkit_headers(Exceptions.h
        BadFileException.h
        BoostStartInclude.h
        BoostEndInclude.h
        ControlCHandler.h
        Dict.h
        FileParseException.h
        Invariant.h
        RDAny.h
        RDValue.h
        RDValue-doublemagic.h
        RDValue-taggedunion.h
        RDLog.h
        RDProps.h
        RDThreads.h
        StreamOps.h
        types.h
        utils.h
        versions.h
        RDConfig.h
        LocaleSwitcher.h
        Ranking.h
        hanoiSort.h
        RDExportMacros.h
        export.h
        test.h
        ConcurrentQueue.h
        BetterEnums.h
        enum.h
        DEST RDGeneral)

if (NOT RDK_INSTALL_INTREE)
    install(DIRECTORY hash DESTINATION ${RDKit_HdrDir}/RDGeneral)
endif (NOT RDK_INSTALL_INTREE)

rdkit_catch_test(testDict testDict.cpp LINK_LIBRARIES RDGeneral)
rdkit_catch_test(testRDValue testRDValue.cpp LINK_LIBRARIES RDGeneral)

if (RDK_BUILD_THREADSAFE_SSS)
    rdkit_catch_test(testConcurrentQueue testConcurrentQueue.cpp LINK_LIBRARIES RDGeneral)
endif (RDK_BUILD_THREADSAFE_SSS)

if (RDK_BUILD_CPP_TESTS)
    add_library(rdkitCatch catch_main.cpp)
    target_link_libraries(rdkitCatch PUBLIC rdkit_base Catch2::Catch2)
endif (RDK_BUILD_CPP_TESTS)


rdkit_catch_test(dictTestsCatch catch_dict.cpp
        LINK_LIBRARIES RDGeneral)

rdkit_catch_test(logTestsCatch catch_logs.cpp
        LINK_LIBRARIES RDGeneral)
