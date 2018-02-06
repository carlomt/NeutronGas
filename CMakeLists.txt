#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(TestEm7)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Setup ROOT include directories and compile definitions
#
find_package(ROOT)
if(NOT ROOT_FOUND OR WITHOUT_ROOT)
  message(WARNING "ROOT not found or disabled by hand")
  add_definitions(-D__WITHOUT_ROOT__)
  list(REMOVE_ITEM sources ${PROJECT_SOURCE_DIR}/src/TreeManager.cc)
  list(REMOVE_ITEM sources ${PROJECT_SOURCE_DIR}/src/TreeManagerMessenger.cc)
  list(REMOVE_ITEM headers ${PROJECT_SOURCE_DIR}/include/TreeManager.hh)
  list(REMOVE_ITEM headers ${PROJECT_SOURCE_DIR}/include/TreeManagerMessenger.hh)
else()
  include_directories(${ROOT_INCLUDE_DIRS})
endif()

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(TestEm7 TestEm7.cc ${sources} ${headers})
if(ROOT_FOUND)
  target_link_libraries(TestEm7 ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
else()
  target_link_libraries(TestEm7 ${Geant4_LIBRARIES} )  
endif()
#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build TestEm7. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
file(GLOB TestEm7_SCRIPTS ${PROJECT_SOURCE_DIR}/*.mac)

set(TestEm7_PLOT_SCRIPTS
  TestEm7.in plotHisto.C plot.gp pres.txt 
  )

foreach(_script ${TestEm7_SCRIPTS})
  get_filename_component(barename ${_script} NAME)
  configure_file(
    ${_script}
    ${PROJECT_BINARY_DIR}/${barename}
    COPYONLY
    )
endforeach()

foreach(_script ${TestEm7_PLOT_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS TestEm7 DESTINATION bin)

