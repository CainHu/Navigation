cmake_minimum_required(VERSION 3.16)
project(Navigation)

set(CMAKE_CXX_STANDARD 14)

include_directories("E:\\eigen3")
include_directories("${PROJECT_SOURCE_DIR}\\include\\eskf")

aux_source_directory(${PROJECT_SOURCE_DIR}/src/eskf ESKF)

add_executable(Navigation main.cpp ${ESKF})