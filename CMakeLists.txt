# Setup the project
project(restG4)

message("===============  ${PROJECT_NAME}  ==============")
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake
                      ${CMAKE_MODULE_PATH})

if (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    # this will only be accessed when build restG4 as a standalone package, not
    # when building the framework
    cmake_minimum_required(VERSION 3.16)
    find_package(REST REQUIRED)
    string(FIND ${REST_LIBRARIES} "RestGeant4" REST_GEANT4_FOUND)
    if (${REST_GEANT4_FOUND} EQUAL -1)
        message(
            SEND_ERROR
                "REST found but 'RestGeant4' library is not present. Please install REST with Geant4Lib"
        )
    endif ()
endif ()

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS RIO Geom)
execute_process(COMMAND root-config --cflags OUTPUT_VARIABLE ROOT_CFLAGS)
string(STRIP ${ROOT_CFLAGS} ROOT_CFLAGS)
message(
    STATUS
        "-- Found ROOT version: ${ROOT_VERSION} with compilation flags: ${ROOT_CFLAGS} and libraries: ${ROOT_LIBRARIES}"
)

# Find Geant4
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})
message("-- Found Geant4 version: ${Geant4_VERSION}")

# Check Geant4 C++ standard is correct
execute_process(COMMAND geant4-config --cxxstd OUTPUT_VARIABLE GEANT4_CXX_STD)
if (NOT ${GEANT4_CXX_STD} MATCHES "17")
    message(
        FATAL_ERROR
            "Geant4 installation was compiled with C++${GEANT4_CXX_STD} standard, but C++17 is required for REST"
    )
endif ()

# Fix for older Geant4 versions
if (${Geant4_VERSION} VERSION_LESS 11.0.0)
    add_compile_definitions(GEANT4_VERSION_LESS_11_0_0)
endif ()

# Fix to allow compatibility with older Geant4 versions that don't have
# G4RunManagerFactory
if (NOT EXISTS "${Geant4_INCLUDE_DIRS}/G4RunManagerFactory.hh")
    add_compile_definitions(GEANT4_WITHOUT_G4RunManagerFactory)
endif ()

if (NOT DEFINED CMAKE_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX ${REST_PATH})
endif ()

message(
    STATUS
        "Package ${PROJECT_NAME} will be installed in ${CMAKE_INSTALL_PREFIX}")

if (NOT CMAKE_CXX_FLAGS)
    set(CMAKE_CXX_FLAGS " -std=c++1y")
endif ()

# Set include and lib
set(INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS}
                 ${Geant4_INCLUDE_DIR} ${REST_INCLUDE_DIRS})

if (NOT ${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    set(INCLUDE_DIRS ${INCLUDE_DIRS} ${rest_include_dirs})
endif ()

set(LINK_LIBRARIES ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} RestFramework
                   RestGeant4)
string(STRIP "${LINK_LIBRARIES}" LINK_LIBRARIES)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cxx)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)

# Create a library, necessary to enable testing
set(LIBRARY RestRestG4) # this way the `ADD_LIBRARY_TEST` macro will pick the
                        # library name correctly
add_library(${LIBRARY} SHARED)

target_sources(${LIBRARY} PUBLIC ${sources})
target_link_libraries(${LIBRARY} ${LINK_LIBRARIES})
target_include_directories(${LIBRARY} SYSTEM PUBLIC ${INCLUDE_DIRS})

if (NOT ${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    add_library_test()
endif ()

# Add the executable
add_executable(${PROJECT_NAME} main.cxx)
target_link_libraries(${PROJECT_NAME} PUBLIC ${LIBRARY})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/examples
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

install(
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/examples/
    DESTINATION ./examples/restG4/
    COMPONENT install)

install(
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/mac
    DESTINATION macros
    COMPONENT install)

install(TARGETS ${PROJECT_NAME} DESTINATION bin)
if( CMAKE_SYSTEM_NAME MATCHES "Darwin")
    install(TARGETS ${LIBRARY} DESTINATION bin)
else()
    install(TARGETS ${LIBRARY} DESTINATION lib)
endif()

if (NOT ${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    set(package_added "${PROJECT_NAME} ")
    set(package_added
        ${package_added}
        PARENT_SCOPE)
endif ()
