#[[----------------------------------------------------------------------#
#  This file is a part of ShapeMapper, and is licensed under the terms   #
#  of the MIT license. Copyright 2018 Steven Busan.                      #
#----------------------------------------------------------------------]]#

set(CMAKE_EXE_LINKER_FLAGS "-static-libstdc++")

include_directories(${Boost_INCLUDE_DIRS})

add_executable(shapemapper_read_trimmer ReadTrimmerExe.cpp)
target_link_libraries(
        shapemapper_read_trimmer
        ${Boost_LIBRARIES}
        ${ZLIB_LIBRARIES} # order is important, since this lib will be dynamically linked
)

add_executable(shapemapper_mutation_parser MutationParserExe.cpp)
target_link_libraries(
        shapemapper_mutation_parser
        ${Boost_LIBRARIES}
        ${ZLIB_LIBRARIES} # order is important, since this lib will be dynamically linked
)

add_executable(shapemapper_mutation_counter MutationCounterExe.cpp)
target_link_libraries(
        shapemapper_mutation_counter
        ${Boost_LIBRARIES}
        ${ZLIB_LIBRARIES} # order is important, since this lib will be dynamically linked
)
