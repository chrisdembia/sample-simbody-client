# Generic CMakeLists.txt for making a Simbody-using executable.
# This shows how to use the provided SimbodyConfig.cmake to locate a Simbody
# installation on your machine so you can use it from your own code.
# You will most likely want to copy some of these lines into your own
# CMakeLists.txt rather than use this one verbatim.
cmake_minimum_required(VERSION 2.8)
project(sample-simbody-client)

# This depends on SimbodyConfig.cmake being located somewhere predictable
# on your machine. If you have installed it somewhere that CMake won't be
# able to guess, you'll need to tell find_package where to look.
find_package(Simbody REQUIRED)

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU" OR
        "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

include_directories(${Simbody_INCLUDE_DIR})
link_directories(${Simbody_LIB_DIR})

add_executable(myexe main.cpp)
target_link_libraries(myexe ${Simbody_LIBRARIES})
