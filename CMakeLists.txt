# Setup the project
project(restG4)

message("===============  ${PROJECT_NAME}  ==============")
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake ${CMAKE_MODULE_PATH})

# Find ROOT
set(ROOT_REQUIRED_LIBRARIES ${ROOT_REQUIRED_LIBRARIES} Gui Geom Gdml)
find_package(ROOT REQUIRED COMPONENTS ${ROOT_REQUIRED_LIBRARIES})

execute_process(COMMAND root-config --cflags OUTPUT_VARIABLE ROOT_CFLAGS)
string(STRIP ${ROOT_CFLAGS} ROOT_CFLAGS)
message(STATUS "-- Found ROOT version: ${ROOT_VERSION} with compilation flags: ${ROOT_CFLAGS} and libraries: ${ROOT_LIBRARIES}")

# Find Geant4
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})
message("-- Found Geant4 version: ${Geant4_VERSION}")

if (${Geant4_VERSION} VERSION_GREATER "10.2.9")
    add_definitions(-DG4104)
endif ()

# Fix to allow compatibility with older Geant4 versions that don't have G4RunManagerFactory
if (NOT EXISTS "${Geant4_INCLUDE_DIRS}/G4RunManagerFactory.hh")
    add_definitions(-DWITHOUT_G4RunManagerFactory)
endif ()

# Find REST
if (NOT DEFINED REST_PATH)
    if (DEFINED ENV{REST_PATH})
        set(REST_PATH $ENV{REST_PATH})
    else ()
        message(FATAL_ERROR "ERROR ::: REST_PATH must be defined as an environment variable and point to REST install directory")
        return()
    endif ()
endif ()

if (NOT DEFINED rest_include_dirs)
    set(rest_include_dirs ${REST_PATH}/include)
endif ()

if (NOT DEFINED CMAKE_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX ${REST_PATH})
endif ()
message(STATUS "Using REST_PATH: ${REST_PATH}")
message(STATUS "Using rest_include_dirs: ${rest_include_dirs}")
message(STATUS "Package ${PROJECT_NAME} will be installed in ${CMAKE_INSTALL_PREFIX}")

if (NOT CMAKE_CXX_FLAGS)
    set(CMAKE_CXX_FLAGS " -std=c++1y")
endif ()

# Set include and lib
set(INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS} ${rest_include_dirs} ${Geant4_INCLUDE_DIR})

set(LINK_LIBRARIES ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} RestGeant4)
string(STRIP "${LINK_LIBRARIES}" LINK_LIBRARIES)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cxx)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)

# Add the executable, and link it to the Geant4 libraries
add_executable(${PROJECT_NAME} restG4.cxx ${sources} ${headers})
target_link_libraries(${PROJECT_NAME} ${LINK_LIBRARIES})
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${INCLUDE_DIRS})

# Copy scripts to the build directory
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/examples/
        DESTINATION ./examples/restG4/
        COMPONENT install
        )

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/mac
        DESTINATION macros
        COMPONENT install
        )

# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
install(TARGETS ${PROJECT_NAME} DESTINATION bin)

set(package_added "${PROJECT_NAME} ")
set(package_added ${package_added} PARENT_SCOPE)
