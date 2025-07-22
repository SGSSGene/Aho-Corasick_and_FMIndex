# SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
# SPDX-License-Identifier: CC0-1.0
cmake_minimum_required (VERSION 3.25)

project(aho_fm LANGUAGES CXX C
               DESCRIPTION "aho_fm -- benchmarking aho-corasick and fm index")


include(cmake/CPM.cmake)
CPMAddPackage("gh:SGSSGene/cpm.dependencies@1.0.1")
CPMLoadDependenciesFile("${CMAKE_CURRENT_SOURCE_DIR}/cpm.dependencies")

add_executable(aho_fm src/aho_fm/main.cpp)
target_compile_features(aho_fm PUBLIC cxx_std_23)
target_link_libraries(aho_fm PRIVATE
    cereal::cereal
    clice::clice
    fmindex-collection::fmindex-collection
    fmt::fmt
    ivio::ivio
    ivsigma::ivsigma
    channel::channel
)
