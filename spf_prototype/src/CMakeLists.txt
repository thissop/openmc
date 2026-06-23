cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(spf_source CXX)

# Build the spin-polarized fusion CompiledSource as a shared library.
#   cmake -B build -DCMAKE_PREFIX_PATH=$HOME/spf_venv .
#   cmake --build build
# (CMAKE_PREFIX_PATH must point at the OpenMC install prefix; see NOTES_orientation.md.)
add_library(polarized_fusion_source SHARED polarized_fusion_source.cpp)

find_package(OpenMC REQUIRED)
target_link_libraries(polarized_fusion_source OpenMC::libopenmc)

# spf_sampler.hpp lives next to this source
target_include_directories(polarized_fusion_source PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(polarized_fusion_source PROPERTIES CXX_STANDARD 17)
