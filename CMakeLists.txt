#=========================================================
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
CMAKE_POLICY(SET CMP0012 NEW)
CMAKE_POLICY(SET CMP0025 NEW)

#=========================================================
PROJECT(Gate)

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake-modules ${CMAKE_MODULE_PATH})

#=========================================================
# Configure Dart testing support.  This should be done before any
# MESSAGE(FATAL_ERROR ...) commands are invoked.
OPTION(BUILD_TESTING "Build the testing tree." OFF)
INCLUDE(CTest)

IF(BUILD_TESTING)
  ENABLE_TESTING()
  SET(BUILDNAME "${BUILDNAME}" CACHE STRING "Name of build on the dashboard")
  MARK_AS_ADVANCED(BUILDNAME)
ENDIF(BUILD_TESTING)
#=========================================================

#=========================================================
# Help function to debug CMake
MACRO (DD in)
  MESSAGE(${in}=[${${in}}])
ENDMACRO(DD)

#=========================================================
# Default build type
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

#=========================================================
# Option to use visualisation feature of G4
OPTION(GATE_USE_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
IF(GATE_USE_GEANT4_UIVIS)
  FIND_PACKAGE(Geant4 REQUIRED ui_all vis_all)
ELSE()
  FIND_PACKAGE(Geant4 REQUIRED)
ENDIF()

#=========================================================
# Refuse to build with MT geant4, allow dev override
IF(Geant4_multithreaded_FOUND)
    MESSAGE(FATAL_ERROR "GATE should be compiled with a non-multithreaded installation of Geant4")
ENDIF()
#OPTION(GATE_ALLOW_MT_GEANT4 "Allow multithreaded geant4 (for MT Gate development)" OFF)
#IF(Geant4_multithreaded_FOUND)
    #IF(GATE_ALLOW_MT_GEANT4)
        #MESSAGE(WARNING "GATE does not yet work correctly with multithreaded Geant4, proceed with caution...")
    #ELSE()
        #MESSAGE(FATAL_ERROR "GATE should be compiled with a non-multithreaded installation of Geant4")
    #ENDIF()
#ENDIF()

# Check if OpenGL headers are still available
IF(Geant4_qt_FOUND OR Geant4_vis_opengl_x11_FOUND)
  FIND_PACKAGE(OpenGL QUIET)
  INCLUDE_DIRECTORIES(${OPENGL_INCLUDE_DIR})
  SET(GATE_USE_OPENGL ${OPENGL_FOUND})
ENDIF()

#=========================================================
# Build variables for G4 Version

STRING(REGEX REPLACE "(.*)[.](.*)[.](.*)" "\\1" G4VERSION_MAJOR ${Geant4_VERSION})
STRING(REGEX REPLACE "(.*)[.](.*)[.](.*)" "\\2" G4VERSION_MINOR ${Geant4_VERSION})
STRING(REGEX REPLACE "(.*)[.](.*)[.](.*)" "\\3" G4VERSION_PATCH ${Geant4_VERSION})
#MESSAGE(${Geant4_VERSION})
#MESSAGE(${G4VERSION_MAJOR})
#MESSAGE(${G4VERSION_MINOR})
#MESSAGE(${G4VERSION_PATCH})
IF(NOT ${G4VERSION_MAJOR}.${G4VERSION_MINOR} EQUAL 10.4 AND NOT ${G4VERSION_MAJOR}.${G4VERSION_MINOR} EQUAL 10.5)
  MESSAGE("Warning! GATE version 8.2 is not validated for Geant4 ${G4VERSION_MAJOR}.${G4VERSION_MINOR} distribution. Please use Geant4 10.4 or 10.5 distribution instead.")
ENDIF()

#=========================================================
# Option for Root use
SET(GATE_USE_ROOT ON)
SET(GATE_ANALYSIS_USE_GENERAL ON)
SET(GATE_ANALYSIS_USE_FILE ON)

OPTION(GATE_USE_OPTICAL "Gate use optical photon" OFF)
OPTION(GATE_USE_XRAYLIB "Option for XRAYLIB use" OFF)
OPTION(GATE_USE_DAVIS "Gate use Davis LUT model" OFF)
OPTION(GATE_USE_LMF "Gate use LMF" OFF)
OPTION(GATE_USE_ECAT7 "Gate use ECAT7" OFF)
OPTION(GATE_USE_SYSTEM_CLHEP "If 'ON', Gate does not use the standard CLHEP of GEANT4. Use OFF if you compile G4 with embedded CLHEP" OFF)
OPTION(GATE_USE_GPU "Gate use GPU (voxelized source and tracking), need CUDA " OFF)

#=========================================================
# RTK / ITK
option(GATE_USE_RTK "Use the Reconstruction Toolkit (RTK, requires also ITK)" OFF)
option(GATE_USE_ITK "Use the Insight Toolkit (ITK, required by RTK, DICOM reader and thermal therapy capability)" OFF)

#=========================================================
# ROOT
IF(GATE_USE_ROOT)
  FIND_PACKAGE(ROOT REQUIRED)
  INCLUDE_DIRECTORIES(${ROOT_INCLUDE_DIR})
  SET(LINK_DIRECTORIES ${ROOT_LIBRARY_DIR})
  #MESSAGE("Found Root in ${ROOT_INCLUDE_DIR} ${ROOT_LIBRARY_DIR} ${ROOT_LIBRARIES}")
  #IF(ROOT_INCLUDE_DIR)
    #MESSAGE("Your current Root version is ${ROOT_VERSION}")
  #ELSE()
  IF(NOT ROOT_INCLUDE_DIR)
    MESSAGE("Root not found! Please install Root and make it available in PATH")
  ENDIF()
ENDIF()

#=========================================================
# Option for Optical Photon use
IF(GATE_USE_OPTICAL)
  FIND_PACKAGE(LIBXML2 REQUIRED)
  INCLUDE_DIRECTORIES(${LIBXML2_INCLUDE_DIR})
  SET(LINK_DIRECTORIES ${LIBXML2_LIBRARIES})
ENDIF(GATE_USE_OPTICAL)

#=========================================================
# Option for XRAYLIB use
IF(GATE_USE_XRAYLIB)
  FIND_PACKAGE(LIBXRL REQUIRED)
  INCLUDE_DIRECTORIES(${LIBXRL_INCLUDE_DIR})
  SET(LINK_DIRECTORIES ${LIBXRL_LIBRARIES})
  IF(NOT LIBXRL_INCLUDE_DIR)
    MESSAGE("xraylib not found! Please install xraylib and make it available in PATH or manually set the LIBXRL_LIBRARY and LIBXRL_INCLUDE variables.")
  ENDIF()
ENDIF(GATE_USE_XRAYLIB)

#=========================================================
# Option for Davis model
IF(GATE_USE_DAVIS)
  MESSAGE("Warning! You have to install a manually modified version of Geant4 to use the Davis LUT model. Please refer to http://wiki.opengatecollaboration.org/index.php/Enabling_LUT_Davis_Model for detailed step-by-step compiling instructions")
ENDIF(GATE_USE_DAVIS)

#=========================================================
# Option for Clhep use
IF(GATE_USE_SYSTEM_CLHEP)
  FIND_PACKAGE(CLHEP REQUIRED)
  INCLUDE_DIRECTORIES(${CLHEP_INCLUDE_DIR})
  SET(LINK_DIRECTORIES ${CLHEP_LIBRARIES})
  MESSAGE("Your current CLHEP version is ${CLHEP_VERSION}")
  MESSAGE("Please refer to http://geant4.web.cern.ch/geant4/support/ReleaseNotes4.${G4VERSION_MAJOR}.${G4VERSION_MINOR}.html to verify that your CLHEP version is supported by Geant4.")

  #MESSAGE("Found Clhep in ${CLHEP_INCLUDE_DIR} ${CLHEP_LIBRARY}")
ENDIF()

#=========================================================
# Option for GPU use
IF(GATE_USE_GPU)
  FIND_PACKAGE(CUDA)
  INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/source/gpu/include/)
  FILE(GLOB sourcesGPU
    ${PROJECT_SOURCE_DIR}/source/gpu/src/GateGPUManager.cu
    )
ENDIF(GATE_USE_GPU)

#=========================================================
# Option for LMF use
IF(GATE_USE_LMF)
  FIND_PATH(LMF_HOME lmf /usr /usr/local)
  SET(LMF_INCLUDE_DIR "${LMF_HOME}/includes")
  SET(LMF_LIBRARIES "${LMF_HOME}/lib")
  SET(LMF_LIBRARY "-L${LMF_LIBRARIES} -lLMF")
  INCLUDE_DIRECTORIES(${LMF_INCLUDE_DIR})
  SET(LINK_DIRECTORIES ${LMF_LIBRARIES})
  #MESSAGE("Found Lmf in ${LMF_INCLUDE_DIR} ${LMF_LIBRARIES} ${LMF_LIBRARY}")
ENDIF(GATE_USE_LMF)

#=========================================================
# Option for ECAT7 use
IF(GATE_USE_ECAT7)
  FIND_PATH(ECAT7_HOME ecat /usr /usr/local)
  SET(ECAT7_INCLUDE_DIR "${ECAT7_HOME}/include")
  SET(ECAT7_LIBRARIES "${ECAT7_HOME}/lib")
  SET(ECAT7_LIBRARY "-L${ECAT7_LIBRARIES} -lecat")
  INCLUDE_DIRECTORIES(${ECAT7_INCLUDE_DIR})
  SET(LINK_DIRECTORIES ${ECAT7_LIBRARIES})
  #MESSAGE("Found Ecat in ${ECAT7_INCLUDE_DIR} ${ECAT7_LIBRARIES} ${ECAT7_LIBRARY}")
ENDIF(GATE_USE_ECAT7)

#=========================================================
# Create configuration file
CONFIGURE_FILE(GateConfiguration.h.in ${PROJECT_BINARY_DIR}/GateConfiguration.h)
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR})

#=========================================================
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
INCLUDE(${Geant4_USE_FILE})
INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/source/arf/include
  ${PROJECT_SOURCE_DIR}/source/geometry/include
  ${PROJECT_SOURCE_DIR}/source/physics/include
  ${PROJECT_SOURCE_DIR}/source/digits_hits/include
  ${PROJECT_SOURCE_DIR}/source/general/include
  ${PROJECT_SOURCE_DIR}/source/externals/clhep/include)

#=========================================================
# Locate sources and headers for this project
FILE(GLOB sources
  ${PROJECT_SOURCE_DIR}/source/arf/src/*.cc
  ${PROJECT_SOURCE_DIR}/source/geometry/src/*.cc
  ${PROJECT_SOURCE_DIR}/source/digits_hits/src/*.cc
  ${PROJECT_SOURCE_DIR}/source/physics/src/*.cc
  ${PROJECT_SOURCE_DIR}/source/general/src/*.cc
  ${PROJECT_SOURCE_DIR}/source/externals/clhep/src/CLHEP/Matrix/*.cc
  ${PROJECT_SOURCE_DIR}/source/externals/clhep/src/CLHEP/RandomObjects/*.cc)

FILE(GLOB headers
  ${PROJECT_SOURCE_DIR}/source/arf/include/*.hh
  ${PROJECT_SOURCE_DIR}/source/geometry/include/*.hh
  ${PROJECT_SOURCE_DIR}/source/physics/include/*.hh
  ${PROJECT_SOURCE_DIR}/source/digits_hits/include/*.hh
  ${PROJECT_SOURCE_DIR}/source/general/include/*.hh
  ${PROJECT_SOURCE_DIR}/source/externals/clhep/include/*.hh
  ${PROJECT_SOURCE_DIR}/source/gpu/include/*.hh)

#=========================================================
# RTK
IF(GATE_USE_RTK OR GATE_USE_ITK)
  SET(ITK_MINIMAL_VERSION 4.9)
  SET(ITK_COMPONENTS ITKIOImageBase;ITKIOGDCM;ITKThresholding;ITKFFT;ITKImageCompose;ITKSmoothing)
  IF(GATE_USE_RTK)
    LIST(APPEND ITK_COMPONENTS RTK)
    SET(ITK_MINIMAL_VERSION 4.12)
  ENDIF()
  FIND_PACKAGE(ITK ${ITK_MINIMAL_VERSION} REQUIRED COMPONENTS ${ITK_COMPONENTS})
  IF(GATE_USE_RTK AND "${RTK_VERSION_MAJOR}.${RTK_VERSION_MINOR}" VERSION_LESS 2.0)
    MESSAGE(FATAL_ERROR "RTK version 2.0 is required")
  ENDIF()
  INCLUDE(${ITK_USE_FILE})
  STRING(REPLACE "-pedantic" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
ELSE()
  ADD_SUBDIRECTORY(source/externals/itk-mhd ${PROJECT_BINARY_DIR}/itk-mhd)
  SET(ITK_LIBRARIES MetaIO itkzlib)
  INCLUDE_DIRECTORIES(
    ${PROJECT_SOURCE_DIR}/source/externals/itk-mhd/
    ${PROJECT_SOURCE_DIR}/source/externals/itk-mhd/itkzlib
    ${PROJECT_SOURCE_DIR}/source/externals/itk-mhd/src
    ${PROJECT_BINARY_DIR}/itk-mhd/)
  SET(ITKMHD_LIBS MetaIO itkzlib)
ENDIF()



#=========================================================
# Add the executable, and link it to the Geant4/ROOT/CLHEP/ITK libraries
IF(GATE_USE_GPU AND CUDA_FOUND)
  SET(CUDA_NVCC_FLAGS "-G;-gencode arch=compute_30,code=sm_30;-gencode arch=compute_32,code=sm_32;-gencode arch=compute_35,code=sm_35;-gencode arch=compute_37,code=sm_37;-gencode arch=compute_50,code=sm_50;-gencode arch=compute_52,code=sm_52;-gencode arch=compute_53,code=sm_53;-gencode arch=compute_60,code=sm_60;-gencode arch=compute_61,code=sm_61;--use_fast_math;-w;--ptxas-options=-v")
  CUDA_ADD_EXECUTABLE(Gate Gate.cc ${sources} ${sourcesGPU} ${headers})
ELSE(GATE_USE_GPU AND CUDA_FOUND)
  ADD_EXECUTABLE(Gate Gate.cc ${sources} ${headers})
ENDIF(GATE_USE_GPU AND CUDA_FOUND)

TARGET_LINK_LIBRARIES(Gate ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} ${CLHEP_LIBRARIES} ${LIBXML2_LIBRARIES} ${LIBXRL_LIBRARIES} ${LMF_LIBRARY} ${ECAT7_LIBRARY} ${ITK_LIBRARIES} pthread)

#=========================================================
INSTALL(TARGETS Gate DESTINATION bin)

find_package(GSL REQUIRED)    # See below (2)
target_link_libraries(Gate GSL::gsl GSL::gslcblas)

#=========================================================
option(GATE_DOWNLOAD_BENCHMARKS_DATA "Download the missing binary data to run gate benchmarks" OFF)
IF(GATE_DOWNLOAD_BENCHMARKS_DATA OR GATE_DOWNLOAD_EXAMPLES_DATA)
  INCLUDE(cmake-modules/GateExternalData.cmake)
ENDIF(GATE_DOWNLOAD_BENCHMARKS_DATA OR GATE_DOWNLOAD_EXAMPLES_DATA)
IF(GATE_DOWNLOAD_BENCHMARKS_DATA)
  ADD_SUBDIRECTORY(benchmarks)
  ExternalData_add_target(GateBenchmarkData)
ENDIF(GATE_DOWNLOAD_BENCHMARKS_DATA)
#=========================================================


#=========================================================
# We remove the warning option "shadow", because there are tons of
# such warning related to clhep/g4 system of units. Additionally,
# we force the c++11 std. This is mandatory for Geant4 >= 10.2.
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-shadow")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#=========================================================

#=========================================================
# Display Geant4 version
MESSAGE("Your current Geant4 version is ${Geant4_VERSION}")

#=========================================================
# Add gprof to RELWITHDEBINFO target
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -pg")
SET(CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO "${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO} -pg")
SET(CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO "${CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO} -pg")
#=========================================================
