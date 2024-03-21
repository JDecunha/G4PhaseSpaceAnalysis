#locate our project source files and make a library
file(GLOB_RECURSE sources ${CMAKE_SOURCE_DIR}/src/*.cc)
add_library(src ${sources})

#set includes
target_include_directories(src PUBLIC ${CMAKE_SOURCE_DIR}/include) #set the include files


