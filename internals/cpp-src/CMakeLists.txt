#[[----------------------------------------------------------------------#
#  This file is a part of ShapeMapper, and is licensed under the terms   #
#  of the MIT license. Copyright 2018 Steven Busan.                      #
#----------------------------------------------------------------------]]#

cmake_minimum_required(VERSION 3.3) # 3.3 wanted by Clion
project(shapemapper)
include(ExternalProject)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../bin)

set(BUILD_TESTS ON) # Build unit tests for C++ modules.

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_C_FLAGS_CACHED ${CMAKE_C_FLAGS})
set(CMAKE_C_FLAGS "-Wall -g -O2")

# DEBUG print all cmake variables
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

find_package(Boost COMPONENTS filesystem program_options iostreams system REQUIRED)
if ( Boost_FOUND )
    message("BOOST include dirs at: " ${Boost_INCLUDE_DIRS})
    message("BOOST libraries at: " ${Boost_LIBRARIES})
else ()
    message("Could not find BOOST libs.")
endif ( Boost_FOUND )

find_package( ZLIB REQUIRED )

# DEBUG print all cmake variables
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

find_package(Doxygen)
if (DOXYGEN_FOUND)

    set(DOXY_IN ${SOURCE_FILES})
    set(DOXY_DIR ${CMAKE_BINARY_DIR}/doxygen_temp)
    set(DOXY_OUT ${DOXY_DIR}/doxy_out)
    set(BOOTSTRAP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/doxygen-bootstrapped)
    set(DOXY_CFG ${DOXY_DIR}/doxygen.cfg)
    add_custom_command(
        OUTPUT ${DOXY_OUT}

        COMMAND mkdir -p ${DOXY_DIR}
        COMMAND rm -fr ${DOXY_DIR}/*
    # generate doxygen config file
        COMMAND doxygen -g "${DOXY_CFG}"
    # replace parameters in config file (don't know why it won't just accept params on command line, gross
        COMMAND sed -i 's|^INPUT\\s*=|INPUT = src README.md|g' ${DOXY_CFG}
        COMMAND sed -i 's|^OUTPUT_DIRECTORY\\s*=|OUTPUT_DIRECTORY = ${DOXY_DIR}|g' ${DOXY_CFG}
        COMMAND sed -i 's|^PROJECT_NAME\\s*=.*|PROJECT_NAME = shapemapper|g' ${DOXY_CFG}
        # allow macro expansion so default parameters show up as actual numbers
        COMMAND sed -i 's|^MACRO_EXPANSION\\s*=.*|MACRO_EXPANSION = YES|g' ${DOXY_CFG}
        # disable sorting so functions are shown in declaration order
        COMMAND sed -i 's|^SORT_MEMBER_DOCS\\s*=.*|SORT_MEMBER_DOCS = NO|g' ${DOXY_CFG}
        # show namespaces even if not explicitly documented
        COMMAND sed -i 's|^EXTRACT_ALL\\s*=.*|EXTRACT_ALL = YES|g' ${DOXY_CFG}
        # don't show source code in docs
        COMMAND sed -i 's|^SOURCE_BROWSER\\s*=.*|SOURCE_BROWSER = NO|g' ${DOXY_CFG}
        # don't show include headers in docs
        #COMMAND sed -i 's|^VERBATIM_HEADERS\\s*=.*|VERBATIM_HEADERS = NO|g' ${DOXY_CFG}
        COMMAND sed -i 's|^SHOW_INCLUDE_FILES\\s*=.*|SHOW_INCLUDE_FILES = NO|g' ${DOXY_CFG}
        # misc
        #COMMAND sed -i 's|^MULTILINE_CPP_IS_BRIEF\\s*=.*|MULTILINE_CPP_IS_BRIEF = YES|g' ${DOXY_CFG}
        COMMAND sed -i 's|^USE_MDFILE_AS_MAINPAGE\\s*=.*|USE_MDFILE_AS_MAINPAGE = README.md|g' ${DOXY_CFG}

    # use prettier CSS formatting with twitter bootstrap
        COMMAND sed -i 's|^HTML_HEADER\\s*=|HTML_HEADER = ${BOOTSTRAP_DIR}/header.html|g' ${DOXY_CFG}
        COMMAND sed -i 's|^HTML_FOOTER\\s*=|HTML_FOOTER = ${BOOTSTRAP_DIR}/footer.html|g' ${DOXY_CFG}
        COMMAND sed -i 's|^HTML_EXTRA_STYLESHEET\\s*=|HTML_EXTRA_STYLESHEET = ${BOOTSTRAP_DIR}/customdoxygen.css|g' ${DOXY_CFG}
    # run doxygen
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND doxygen "${DOXY_CFG}"
    # copy bootstrap-required javascript file into doxygen html directory
        COMMAND cp ${BOOTSTRAP_DIR}/doxy-boot.js ${DOXY_DIR}/html/
    # move generated docs into doc folder for easier commit
        COMMAND rm -fr doc/html/*
        COMMAND mkdir -p doc/html
        COMMAND mv ${DOXY_DIR}/html/* doc/html/
            #COMMAND cp -f ${DOXY_CFG} doc/doxy.cfg
        COMMAND rm -fr ${DOXY_DIR}
        DEPENDS ${DOXY_IN}
    )
    add_custom_target(api_doc ALL DEPENDS ${DOXY_OUT})


else()
    message("Doxygen not found, skipping doc build.")

endif (DOXYGEN_FOUND)

message("BUILD_TESTS: " ${BUILD_TESTS})
if(BUILD_TESTS)
    message("Will build unit tests.")
    add_subdirectory(lib/gtest-1.7.0)
    add_subdirectory(test)
else()
    message("Will not build unit tests.")
endif(BUILD_TESTS)

add_subdirectory(src)
