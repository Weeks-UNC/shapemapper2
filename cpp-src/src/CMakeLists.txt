#[[----------------------------------------------------------------------#
#  This file is a part of ShapeMapper, and is licensed under the terms   #
#  of the MIT license. Copyright 2017 Steven Busan.                      #
#----------------------------------------------------------------------]]#

set(CMAKE_EXE_LINKER_FLAGS "-static-libstdc++")

include_directories(${BAMTOOLS_INCLUDE_DIRS})

add_library(libbamtools STATIC IMPORTED)
set_target_properties(libbamtools PROPERTIES IMPORTED_LOCATION ${BAMTOOLS_LIBRARIES})
add_dependencies(libbamtools bamtools)

include_directories(${Boost_INCLUDE_DIRS})

add_executable(shapemapper_read_trimmer ReadTrimmerExe.cpp)
target_link_libraries(
        shapemapper_read_trimmer
        ${Boost_LIBRARIES}
        ${ZLIB_LIBRARIES} # order is important, since this lib will be dynamically linked
)

add_executable(shapemapper_read_length_filter ReadLengthFilterExe.cpp)
target_link_libraries(
        shapemapper_read_length_filter
        ${Boost_LIBRARIES}
        ${ZLIB_LIBRARIES} # order is important, since this lib will be dynamically linked
)

add_executable(shapemapper_mutation_parser MutationParserExe.cpp)
target_link_libraries(
        shapemapper_mutation_parser
        libbamtools
        ${Boost_LIBRARIES}
        ${ZLIB_LIBRARIES} # order is important, since this lib will be dynamically linked
)

add_executable(shapemapper_mutation_counter MutationCounterExe.cpp)
target_link_libraries(
        shapemapper_mutation_counter
        ${Boost_LIBRARIES}
        ${ZLIB_LIBRARIES} # order is important, since this lib will be dynamically linked
)

