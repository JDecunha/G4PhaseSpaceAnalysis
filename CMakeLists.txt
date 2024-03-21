#CMakeLists file for G4PhaseSpaceAnalysis
cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project (G4PhaseSpaceAnalysis VERSION 1.0.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD "17")


### Find and include ROOT libraries
find_package(ROOT CONFIG REQUIRED)
include("${ROOT_USE_FILE}") #include the ROOT CMakelists which defines $ROOT_LIBRARIES

### Configure executable
add_executable(G4PhaseSpaceAnalysis G4PhaseSpaceAnalysis.cc)

#Add our source files
add_subdirectory(src)

#Link my project
target_link_libraries(G4PhaseSpaceAnalysis PRIVATE src ${ROOT_LIBRARIES})


