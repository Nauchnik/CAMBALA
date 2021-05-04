cmake_minimum_required(VERSION 3.18)
project(normal_modes)

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

add_library(normal_modes normal_modes.cpp $<TARGET_OBJECTS:alglib.o>)

add_executable(normal_modes_sample main.cpp)

target_link_libraries(normal_modes_sample normal_modes)
