project(restG4)

CMAKE_MINIMUM_REQUIRED(VERSION 3.5)
message("===============  restG4  ==============")
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake ${CMAKE_MODULE_PATH})

# Find REST
IF (NOT DEFINED REST_PATH)
    IF (DEFINED ENV{REST_PATH})
        set(REST_PATH $ENV{REST_PATH})
    ELSE ()
        message(FATAL_ERROR "ERROR ::: REST_PATH must be defined as an environment variable and point to REST install directory")
        return()
    ENDIF ()
ENDIF ()
set(CMAKE_INSTALL_PREFIX ${REST_PATH})
message(STATUS "Using REST PATH: " ${REST_PATH})

if (NOT CMAKE_CXX_FLAGS)
    SET(CMAKE_CXX_FLAGS "-std=c++1y")
endif ()

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake / cmake - gui
# to build a batch mode only executable
#----------------------------------------------------------------------------
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if (WITH_GEANT4_UIVIS)
    find_package(Geant4 REQUIRED ui_all vis_all)
else ()
    find_package(Geant4 REQUIRED)
endif ()
include(${Geant4_USE_FILE})
message("-- Found Geant4 version : ${Geant4_VERSION}")

# Find ROOT
find_package(ROOT REQUIRED)
message("-- Found ROOT version : ${ROOT_VERSION}")

# Set include and libs
set(external_libs ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} -L$ENV{ROOTSYS}/lib -lGui -lGeom -lGdml RestGeant4 RestFramework)
set(external_include_dirs ${PROJECT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIR} ${ROOT_INCLUDE_DIRS})

if (DEFINED rest_include_dirs)
    set(external_include_dirs ${external_include_dirs} ${rest_include_dirs})
else ()
    set(external_include_dirs ${external_include_dirs} ${REST_PATH}/include)
endif ()

include_directories(${incdir})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#----------------------------------------------------------------------------
#message( STATUS "external_libs G4: ${external_libs}" )
#message( STATUS "external_include_dirs G4: ${external_include_dirs}" )
#message( STATUS "sources G4: ${sources}" )
#message( STATUS "headers G4: ${headers}" )
add_executable(${PROJECT_NAME} restG4.cc ${sources} ${headers})
target_link_libraries(${PROJECT_NAME} ${external_libs})
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${external_include_dirs})

target_compile_options(${PROJECT_NAME} PRIVATE)


#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build sphereSimLab. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#----------------------------------------------------------------------------
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/examples
        DESTINATION .
        COMPONENT install
        )
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/mac
        DESTINATION macros
        COMPONENT install
        )

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#----------------------------------------------------------------------------
install(TARGETS ${PROJECT_NAME} DESTINATION bin)

set(package_added "${PROJECT_NAME}")
set(package_added ${package_added} PARENT_SCOPE)
