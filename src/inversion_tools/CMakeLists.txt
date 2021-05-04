cmake_minimum_required(VERSION 3.18)
project(inversion_tools)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

if (MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")
else()
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
endif()

set(THIRD_PARTY_DIR ../../third_party)

if (NOT TARGET alglib)
    add_subdirectory(${THIRD_PARTY_DIR}/alglib ${CMAKE_BINARY_DIR}/alglib)
endif()

include_directories(
    ${CMAKE_BINARY_DIR}/alglib/alglib/cpp/src
    ${THIRD_PARTY_DIR}/spectra/include
    ${THIRD_PARTY_DIR}/eigen
)

add_executable(CAMBALA_inversion_tools main.cpp sequential.cpp)

target_link_libraries(CAMBALA_inversion_tools alglib)
