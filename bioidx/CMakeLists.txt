cmake_minimum_required(VERSION 3.10)
project(bioidx C)

set(CMAKE_C_STANDARD 99)

add_library(bioidx SHARED bioidx.c binidx.c)
set_target_properties(bioidx PROPERTIES LIBRARY_OUTPUT_DIRECTORY lib)
install(TARGETS bioidx LIBRARY DESTINATION lib)
install(FILES bioidx.h DESTINATION include)

add_executable(bioidx_test bioidx_test.c)
target_link_libraries(bioidx_test bioidx)