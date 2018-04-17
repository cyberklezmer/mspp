project(mspptest)
cmake_minimum_required(VERSION 2.8)

set( CPLEX_VER "1263" )
set( CPLEX_HEADER_DIR "/opt/ibm/ILOG/CPLEX_Studio1263/cplex/include" )
set( CONCERT_HEADER_DIR "/opt/ibm/ILOG/CPLEX_Studio1263/concert/include" )
set( CPLEX_LIB_DIR "/opt/ibm/ILOG/CPLEX_Studio1263/cplex/lib/x86-64_linux/static_pic" )
set( CPLEX_BIN_DIR "/opt/ibm/ILOG/CPLEX_Studio1263/cplex/bin/x86-64_linux" )
set( CONCERT_LIB_DIR "/opt/ibm/ILOG/CPLEX_Studio1263/concert/lib/x86-64_linux/static_pic" )

file(GLOB SRC_LIST "*.cpp")

#CPLEX - force it to use STD headers like <iostream> instead of <iostream.h>
add_definitions( -DIL_STD )
#CPLEX - headers, put your directory
include_directories(${CPLEX_HEADER_DIR} )
include_directories(${CONCERT_HEADER_DIR} )
include_directories(..)

#CPLEX - libraries to link
find_library(ILOCPLEX libilocplex.a ${CPLEX_LIB_DIR} )
find_library(CPLEX libcplex.a ${CPLEX_LIB_DIR} )
find_library(BINCPLEX libcplex${CPLEX_VER}.so ${CPLEX_BIN_DIR} )
find_library(CONCERT libconcert.a ${CONCERT_LIB_DIR} )

#Threads for boost and CPLEX
find_package (Threads REQUIRED)


add_executable(${PROJECT_NAME} ${SRC_LIST})

#Link with specified libraries (required by CPLEX)
target_link_libraries(${PROJECT_NAME} ${ILOCPLEX} )
target_link_libraries(${PROJECT_NAME} ${CPLEX} )
target_link_libraries(${PROJECT_NAME} ${BINCPLEX} )
target_link_libraries(${PROJECT_NAME} ${CONCERT} )
target_link_libraries(${PROJECT_NAME} ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(${PROJECT_NAME} ${CMAKE_DL_LIBS} )

#aux_source_directory(. SRC_LIST)

