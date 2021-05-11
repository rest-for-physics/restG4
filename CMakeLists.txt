project(restG4)

CMAKE_MINIMUM_REQUIRED(VERSION 3.5)
message("==================================================  restG4  =================================================")
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake ${CMAKE_MODULE_PATH})

# Find REST
find_package(REST REQUIRED)

set(CMAKE_INSTALL_PREFIX ${REST_PATH})
message("-- Found REST version : ${REST_VERSION}")

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
set(external_libs ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} ${REST_LIBRARIES})
set(external_include_dirs ${PROJECT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIR} ${ROOT_INCLUDE_DIRS} ${REST_INCLUDE_DIRS})

include_directories(${incdir})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(${PROJECT_NAME} restG4.cc ${sources} ${headers})
target_link_libraries(${PROJECT_NAME} ${external_libs})
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${external_include_dirs})

target_compile_options(${PROJECT_NAME} PRIVATE)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build sphereSimLab. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#----------------------------------------------------------------------------
file(GLOB examples
        "examples/*"
        )
file(COPY ${examples} DESTINATION examples)

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
