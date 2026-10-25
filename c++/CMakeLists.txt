cmake_minimum_required(VERSION 3.16)
project(spectrum_intel VERSION 0.2.0 LANGUAGES CXX)

# ---------------------------------------------------------------- options
option(SI_BUILD_TESTS "Build unit tests" ON)
option(SI_BUILD_APPS "Build the si_analyze command-line tool" ON)
option(SI_NATIVE "Tune for the build machine (-mcpu/-march=native); off for cross builds" OFF)
option(SI_SHARED "Build libspectrum_intel as a shared library (C API for other languages)" OFF)
set(SI_RECEIVERS "" CACHE STRING
    "Receivers to build, ';'-separated (empty = every folder in src/receivers, none = no receiver). Removing one removes it from the binary.")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Release: -O3 without -ffast-math. The library must reproduce the Python reference's
# decisions (thresholds on IEEE results); fast-math reassociation would move them.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  add_compile_options(-Wall -Wextra -Wno-unused-parameter -Wno-sign-compare)
  # no fused multiply-add contraction: the same IEEE results on x86 (GCC contracts by default) and ARM
  add_compile_options(-ffp-contract=off)
  # errno is never read: lets GCC vectorise sqrt (clang on Apple already does); no effect on results
  add_compile_options(-fno-math-errno)
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
  if(SI_NATIVE)
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag(-mcpu=native SI_HAS_MCPU_NATIVE)
    if(SI_HAS_MCPU_NATIVE)
      add_compile_options(-mcpu=native)
    else()
      add_compile_options(-march=native)
    endif()
  endif()
endif()

find_package(Threads REQUIRED)

set(SI_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
set(SI_REPO ${CMAKE_CURRENT_SOURCE_DIR}/..)

# ---------------------------------------------------------------- shared contract (schema/)
# The JSON files in <repo>/schema are compiled into the library, so a deployed binary
# carries the same response schema / taxonomy / requirements as the Python package.
include(cmake/EmbedFiles.cmake)
si_embed_files(${CMAKE_CURRENT_BINARY_DIR}/generated/si_embedded_schema.cpp
               ${SI_REPO}/schema/technologies.json
               ${SI_REPO}/schema/requirements.json
               ${SI_REPO}/schema/response.schema.json
               ${SI_REPO}/schema/capabilities.json
               ${SI_REPO}/schema/descriptions.json)

# ---------------------------------------------------------------- third party (header-only)
add_library(si_third_party INTERFACE)
target_include_directories(si_third_party SYSTEM INTERFACE
  ${SI_ROOT}/third_party/pocketfft ${SI_ROOT}/third_party/nlohmann ${SI_ROOT}/third_party)
target_compile_definitions(si_third_party INTERFACE POCKETFFT_CACHE_SIZE=64)

# ---------------------------------------------------------------- core: numpy/scipy-equivalent DSP + records
file(GLOB SI_CORE_SRC CONFIGURE_DEPENDS ${SI_ROOT}/src/core/*.cpp ${SI_ROOT}/src/io/*.cpp)
add_library(si_core STATIC ${SI_CORE_SRC} ${CMAKE_CURRENT_BINARY_DIR}/generated/si_embedded_schema.cpp)
# fft.cpp wraps pocketfft, which numpy and scipy compile with the compiler's default FP
# contraction (clang: on, GCC: fast). Matching that makes the FFTs bit-identical to
# numpy/scipy on the same platform; everything else stays -ffp-contract=off.
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set_source_files_properties(${SI_ROOT}/src/core/fft.cpp PROPERTIES COMPILE_OPTIONS "-ffp-contract=on")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set_source_files_properties(${SI_ROOT}/src/core/fft.cpp PROPERTIES COMPILE_OPTIONS "-ffp-contract=fast")
endif()
target_include_directories(si_core PUBLIC ${SI_ROOT}/src ${SI_ROOT}/include)
target_link_libraries(si_core PUBLIC si_third_party Threads::Threads)

# ---------------------------------------------------------------- plugin interface
file(GLOB SI_DECODE_SRC CONFIGURE_DEPENDS ${SI_ROOT}/src/decode/*.cpp)
list(FILTER SI_DECODE_SRC EXCLUDE REGEX ".*/registry\\.cpp$")   # needs the generated receiver list: built with spectrum_intel
add_library(si_decode STATIC ${SI_DECODE_SRC})
target_link_libraries(si_decode PUBLIC si_core)

# ---------------------------------------------------------------- receivers (one target each)
if(SI_RECEIVERS STREQUAL "")
  file(GLOB _rx_dirs LIST_DIRECTORIES true ${SI_ROOT}/src/receivers/*)
  set(SI_RECEIVER_LIST "")
  foreach(d ${_rx_dirs})
    if(IS_DIRECTORY ${d} AND EXISTS ${d}/plugin.cpp)
      get_filename_component(n ${d} NAME)
      list(APPEND SI_RECEIVER_LIST ${n})
    endif()
  endforeach()
elseif(SI_RECEIVERS STREQUAL "none")
  set(SI_RECEIVER_LIST "")
else()
  set(SI_RECEIVER_LIST ${SI_RECEIVERS})
endif()
list(SORT SI_RECEIVER_LIST)
set(SI_RECEIVER_TARGETS "")
foreach(n ${SI_RECEIVER_LIST})
  file(GLOB _src CONFIGURE_DEPENDS ${SI_ROOT}/src/receivers/${n}/*.cpp)
  add_library(si_rx_${n} STATIC ${_src})
  # a receiver's DSP sees only core; plugin.cpp also sees decode/ (checked by tools/check_modularity.py)
  target_link_libraries(si_rx_${n} PUBLIC si_decode)
  list(APPEND SI_RECEIVER_TARGETS si_rx_${n})
endforeach()
message(STATUS "spectrum_intel receivers: ${SI_RECEIVER_LIST}")

# the registry calls si::rx::<name>::register_plugins() for every built receiver
set(_reg ${CMAKE_CURRENT_BINARY_DIR}/generated/si_receivers.cpp)
set(_decl "")
set(_call "")
foreach(n ${SI_RECEIVER_LIST})
  string(APPEND _decl "namespace si::rx::${n} { void register_plugins(si::PluginList&); }\n")
  string(APPEND _call "    {\"${n}\", &si::rx::${n}::register_plugins},\n")
endforeach()
file(WRITE ${_reg}.in
"// generated by CMake from SI_RECEIVERS: do not edit
#include \"decode/registry.hpp\"
${_decl}
namespace si {
const std::vector<ReceiverEntry>& built_receivers() {
  static const std::vector<ReceiverEntry> list = {
${_call}  };
  return list;
}
}  // namespace si
")
configure_file(${_reg}.in ${_reg} COPYONLY)

# ---------------------------------------------------------------- pipeline + API
file(GLOB SI_PIPE_SRC CONFIGURE_DEPENDS
  ${SI_ROOT}/src/detect/*.cpp ${SI_ROOT}/src/features/*.cpp ${SI_ROOT}/src/classify/*.cpp
  ${SI_ROOT}/src/report/*.cpp ${SI_ROOT}/src/api/*.cpp ${SI_ROOT}/src/*.cpp ${SI_ROOT}/src/decode/registry.cpp)
if(SI_SHARED)
  add_library(spectrum_intel SHARED ${SI_PIPE_SRC} ${_reg})
else()
  add_library(spectrum_intel STATIC ${SI_PIPE_SRC} ${_reg})
endif()
target_link_libraries(spectrum_intel PUBLIC si_decode ${SI_RECEIVER_TARGETS})
target_compile_definitions(spectrum_intel PRIVATE SI_VERSION="${PROJECT_VERSION}")

# ---------------------------------------------------------------- apps / tests
if(SI_BUILD_APPS AND EXISTS ${SI_ROOT}/apps/si_analyze.cpp)
  add_executable(si_analyze ${SI_ROOT}/apps/si_analyze.cpp)
  target_link_libraries(si_analyze PRIVATE spectrum_intel)
endif()

if(SI_BUILD_TESTS)
  enable_testing()
  file(GLOB SI_TEST_SRC CONFIGURE_DEPENDS ${SI_ROOT}/tests/unit/*.cpp)
  # test_rx_<name>*.cpp only for the receivers in this build
  foreach(t ${SI_TEST_SRC})
    get_filename_component(tn ${t} NAME_WE)
    if(tn MATCHES "^test_rx_([a-z0-9]+(_[a-z0-9]+)*)")
      set(_keep FALSE)
      foreach(n ${SI_RECEIVER_LIST})
        if(tn STREQUAL "test_rx_${n}" OR tn MATCHES "^test_rx_${n}_")
          set(_keep TRUE)
        endif()
      endforeach()
      if(NOT _keep)
        list(REMOVE_ITEM SI_TEST_SRC ${t})
      endif()
    endif()
  endforeach()
  add_executable(si_tests ${SI_TEST_SRC})
  target_link_libraries(si_tests PRIVATE spectrum_intel)
  target_compile_definitions(si_tests PRIVATE SI_TESTDATA="${SI_REPO}/testdata" SI_TESTVEC="${SI_ROOT}/tests/vectors")
  add_test(NAME unit COMMAND si_tests)
endif()
