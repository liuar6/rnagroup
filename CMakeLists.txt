cmake_minimum_required(VERSION 3.10)
project(RNAparser C)

set(CMAKE_C_STANDARD 99)

add_subdirectory(bioidx)
add_executable(rgfind rgfind.c file.c common.c)
add_executable(rgcount rgcount.c file.c common.c)
add_executable(rgmerge rgmerge.c file.c common.c)
target_link_libraries(rgfind hts m bioidx)
target_link_libraries(rgcount hts m bioidx)
target_link_libraries(rgmerge hts m bioidx)