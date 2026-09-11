cmake_minimum_required(VERSION 3.15)
project(cifi LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(ZLIB REQUIRED)

# Find htslib
find_library(HTSLIB hts HINTS /usr/lib /usr/local/lib /opt/homebrew/lib)
find_path(HTSLIB_INCLUDE htslib/sam.h HINTS /usr/include /usr/local/include /opt/homebrew/include)

if(NOT HTSLIB OR NOT HTSLIB_INCLUDE)
    message(FATAL_ERROR "htslib not found. Install with: apt install libhts-dev (Linux) or brew install htslib (macOS)")
endif()

# parallel_hashmap (header-only)
include_directories(${CMAKE_SOURCE_DIR}/extern/parallel_hashmap)

# Find nanobind
execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

# Build the extension module
nanobind_add_module(
    _core
    STABLE_ABI
    src/cifi_cpp/cifi.cpp
    src/cifi_cpp/core/enzyme.cpp
    src/cifi_cpp/core/digestion.cpp
    src/cifi_cpp/core/segment_name.cpp
    src/cifi_cpp/stats/statistics.cpp
    src/cifi_cpp/io/writer.cpp
    src/cifi_cpp/filter/bam_filter.cpp
    src/cifi_cpp/contacts/contacts.cpp
)

target_include_directories(_core PRIVATE 
    ${HTSLIB_INCLUDE}
    ${CMAKE_SOURCE_DIR}/extern
)

target_link_libraries(_core PRIVATE 
    ${HTSLIB}
    ZLIB::ZLIB
)

# Install
install(TARGETS _core LIBRARY DESTINATION cifi)
