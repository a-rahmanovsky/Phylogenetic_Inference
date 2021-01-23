cmake_minimum_required(VERSION 3.17)
project(KR_1_mag)

set(CMAKE_CXX_STANDARD 17)

add_executable(KR_1_mag main.cpp Matrix.h Matrix.cpp GenAlgo.cpp GenAlgo.h)

target_link_libraries(KR_1_mag stdc++fs)
target_link_libraries(KR_1_mag pthread)