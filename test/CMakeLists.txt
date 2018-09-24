project(mspptest)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


###required for sddp

include_directories(../sddp/source)
set( COIN_HEADER_DIR "/usr/local/lib/coin/build/include/coin" )

#COIN-OR headers, put your directory
set( COIN_LIB_DIR "/usr/local/lib/coin/build/lib" )
include_directories(source ${COIN_HEADER_DIR} )
file(GLOB SOURCES "source/*.cpp" "source/alglib/*.cpp" "source/rngstream/*.cpp" )
find_library(COINCBC libCbc.so ${COIN_LIB_DIR} )
find_library(COINCBCSOLVER libCbcSolver.so ${COIN_LIB_DIR} )
find_library(COINCLP libClp.so ${COIN_LIB_DIR} )
find_library(COINCLPSOLVER libClpSolver.so ${COIN_LIB_DIR} )
find_library(COINCGL libCgl.so ${COIN_LIB_DIR} )
find_library(COINUTILS libCoinUtils.so ${COIN_LIB_DIR} )
find_library(COINOSI libOsi.so ${COIN_LIB_DIR} )
find_library(COINOSICLP libOsiClp.so ${COIN_LIB_DIR} )
find_library(COINOSICBC libOsiCbc.so ${COIN_LIB_DIR} )

#Boost libraries
find_package(Boost 1.61 COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

#Armadillo library
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})


###required for cplex

set( CPLEX_VER "1263" )
set( CPLEX_HEADER_DIR "/opt/ibm/ILOG/CPLEX_Studio1263/cplex/include" )
set( CONCERT_HEADER_DIR "/opt/ibm/ILOG/CPLEX_Studio1263/concert/include" )
set( CPLEX_LIB_DIR "/opt/ibm/ILOG/CPLEX_Studio1263/cplex/lib/x86-64_linux/static_pic" )
set( CPLEX_BIN_DIR "/opt/ibm/ILOG/CPLEX_Studio1263/cplex/bin/x86-64_linux" )
set( CONCERT_LIB_DIR "/opt/ibm/ILOG/CPLEX_Studio1263/concert/lib/x86-64_linux/static_pic" )

#CPLEX - force it to use STD headers like <iostream> instead of <iostream.h>
add_definitions( -DIL_STD )
#CPLEX - headers, put your directory
include_directories(${CPLEX_HEADER_DIR} )
include_directories(${CONCERT_HEADER_DIR} )
include_directories(..)
include_directories(.)

#CPLEX - libraries to link
find_library(ILOCPLEX libilocplex.a ${CPLEX_LIB_DIR} )
find_library(CPLEX libcplex.a ${CPLEX_LIB_DIR} )
find_library(BINCPLEX libcplex${CPLEX_VER}.so ${CPLEX_BIN_DIR} )
find_library(CONCERT libconcert.a ${CONCERT_LIB_DIR} )

#Threads for boost and CPLEX
find_package (Threads REQUIRED)

# set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -lmcheck -pedantic")


file(GLOB SRC_LIST "src/*.cpp" "../sddp/source/*.cpp" "../sddp/source/alglib/*.cpp"  "../sddp/source/rngstream/*.cpp" )

add_executable(${PROJECT_NAME} ${SRC_LIST})

#Link with specified libraries (required by CPLEX)
target_link_libraries(${PROJECT_NAME} ${ILOCPLEX} )
target_link_libraries(${PROJECT_NAME} ${CPLEX} )
target_link_libraries(${PROJECT_NAME} ${BINCPLEX} )
target_link_libraries(${PROJECT_NAME} ${CONCERT} )
target_link_libraries(${PROJECT_NAME} ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(${PROJECT_NAME} ${CMAKE_DL_LIBS} )

#by sddp

target_link_libraries(${PROJECT_NAME} ${COINCBC} )
target_link_libraries(${PROJECT_NAME} ${COINCBCSOLVER} )
target_link_libraries(${PROJECT_NAME} ${COINCLP} )
target_link_libraries(${PROJECT_NAME} ${COINCLPSOLVER} )
target_link_libraries(${PROJECT_NAME} ${COINCGL} )
target_link_libraries(${PROJECT_NAME} ${COINUTILS} )
target_link_libraries(${PROJECT_NAME} ${COINOSI} )
target_link_libraries(${PROJECT_NAME} ${COINOSICLP} )
target_link_libraries(${PROJECT_NAME} ${COINOSICBC} )
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${ARMADILLO_LIBRARIES})



