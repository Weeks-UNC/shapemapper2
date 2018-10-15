#[[----------------------------------------------------------------------#
#  This file is a part of ShapeMapper, and is licensed under the terms   #
#  of the MIT license. Copyright 2018 Steven Busan.                      #
#----------------------------------------------------------------------]]#

set(CMAKE_EXE_LINKER_FLAGS "-static-libstdc++")

enable_testing()
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR} ../src)

include_directories(${Boost_INCLUDE_DIRS})

add_executable(test_read_trimmer testReadTrimmer.cpp)
target_link_libraries(
        test_read_trimmer
        gtest gtest_main
        ${Boost_LIBRARIES}
        ${ZLIB_LIBRARIES}
)

add_executable(test_mutation_parser testMutationParser.cpp)
target_link_libraries(
        test_mutation_parser
        gtest gtest_main
        ${Boost_LIBRARIES}
        ${ZLIB_LIBRARIES}
)

add_executable(test_mutation_counter testMutationCounter.cpp)
target_link_libraries(
        test_mutation_counter
        gtest gtest_main
        ${Boost_LIBRARIES}
        ${ZLIB_LIBRARIES}
)

add_executable(test_histogram testHistogram.cpp)
target_link_libraries(
        test_histogram
        gtest gtest_main
)

add_test(run_all_unit_tests test_read_trimmer)
add_test(run_all_unit_tests test_mutation_parser)
add_test(run_all_unit_tests test_mutation_counter)
add_test(run_all_unit_tests test_histogram)

