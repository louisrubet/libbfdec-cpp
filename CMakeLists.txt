cmake_minimum_required(VERSION 3.27)
project(libbfdec-cpp)

add_library(libbf STATIC libbfdec++.cpp libbf/libbf.c libbf/cutils.c)
set_source_files_properties(libbf/libbf.c COMPILE_FLAGS "-DCONFIG_BIGNUM")
target_include_directories(libbf PRIVATE ${PROJECT_SOURCE_DIR} libbf quickjs-fake)

add_executable(libbfdec-cpp-test test/test.cpp)
target_include_directories(libbfdec-cpp-test PRIVATE ${PROJECT_SOURCE_DIR} libbf)
target_link_libraries(libbfdec-cpp-test libbf)
