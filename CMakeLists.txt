cmake_minimum_required(VERSION 2.8)
project( MPC_CACC )

set(Boost_USE_STATIC_LIBS ON)
set(CMAKE_CXX_STANDARD 11)

# ビルドは通るが，未定義動作
set(CMAKE_CXX_FLAGS "-fpermissive")

find_package(Boost REQUIRED)

# find_package(PkgConfig)
# pkg_check_modules(SIMPLEAMQPCLIENT libSimpleAmqpClient)
if(${Boost_FOUND})
message(STATUS "Boost include... ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost library path... ${Boost_LIBRARY_DIRS}")
message(STATUS "Boost library... ${Boost_LIBRARIES}")
else()
message(STATUS "Boost is NOT found")
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/cgmres)
include_directories(${Boost_INCLUDE_DIRS})

add_executable(mpc_cacc main.cpp)
target_link_libraries(mpc_cacc ${Boost_LIBRARIES})