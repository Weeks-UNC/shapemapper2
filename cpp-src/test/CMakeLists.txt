#[[----------------------------------------------------------------------#
#  This file is a part of ShapeMapper, and is licensed under the terms   #
#  of the MIT license. Copyright 2017 Steven Busan.                      #
#----------------------------------------------------------------------]]#

set(CMAKE_EXE_LINKER_FLAGS "-static-libstdc++")

enable_testing()
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR} ../src)

include_directories(${BAMTOOLS_INCLUDE_DIRS})

include_directories(${Boost_INCLUDE_DIRS})

add_library(libbamtools STATIC IMPORTED)
set_target_properties(libbamtools PROPERTIES IMPORTED_LOCATION ${BAMTOOLS_LIBRARIES})
add_dependencies(libbamtools bamtools)

add_executable(test_read_trimmer testReadTrimmer.cpp)
target_link_libraries(
        test_read_trimmer
        gtest gtest_main
        ${Boost_LIBRARIES}
        ${ZLIB_LIBRARIES}
)

add_executable(test_read_length_filter testReadLengthFilter.cpp)
target_link_libraries(
        test_read_length_filter
        gtest gtest_main
        ${Boost_LIBRARIES}
        ${ZLIB_LIBRARIES}
)


add_executable(test_mutation_parser testMutationParser.cpp)
target_link_libraries(
        test_mutation_parser
        libbamtools
        gtest gtest_main
        ${Boost_LIBRARIES}
        ${ZLIB_LIBRARIES}
)

add_executable(test_mutation_counter testMutationCounter.cpp)
target_link_libraries(
        test_mutation_counter
        libbamtools
        gtest gtest_main
        ${Boost_LIBRARIES}
        ${ZLIB_LIBRARIES}
)

add_test(run_all_unit_tests test_read_trimmer)
add_test(run_all_unit_tests test_mutation_parser)
add_test(run_all_unit_tests test_mutation_counter)
add_test(run_all_unit_tests test_read_length_filter)
