cmake_minimum_required(VERSION 2.8.11)

project(SMIV)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/build_utils/CMakeModules/")

if( ${CMAKE_HOST_SYSTEM_VERSION} MATCHES "^3.1" )
  message( "Hopefully, the Ubuntu 14.04 machine, with a modern compiler." )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )
endif()
if( ${CMAKE_HOST_SYSTEM_VERSION} MATCHES "^2.6." )
  message( "It looks like a trusty old Centos 6 machine with an ancient compiler." )
endif()

find_package(OEToolkits COMPONENTS oedepict oeiupac oechem oesystem oeplatform)
find_package(Qt5Widgets)
message( "Qt5Widgets_INCLUDE_DIRS : ${Qt5Widgets_INCLUDE_DIRS}" )
message( "Qt5Widgets_DEFINITIONS : ${Qt5Widgets_DEFINITIONS}" )
message( "Qt5Widgets_COMPILE_DEFINITIONS : ${Qt5Widgets_COMPILE_DEFINITIONS}" )
message( "Qt5Widgets_EXECUTABLE_COMPILE_FLAGS : ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS}" )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS} -Wall" )

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

find_package(Boost COMPONENTS program_options regex iostreams filesystem system date_time REQUIRED)

find_package(Cairo REQUIRED)
message( "CAIRO_INCLUDE_DIR : ${CAIRO_INCLUDE_DIR}" )

include_directories(${SMIV_SOURCE_DIR})
include_directories( SYSTEM
${OEToolkits_INCLUDE_DIRS}
${Qt5Widgets_INCLUDE_DIRS}
${Boost_INCLUDE_DIRS})

set(SMIV_SRCS
smiv_main.cc
SmiV.cc
SmivDataTable.cc
SmiVFindMoleculeDialog.cc
SmiVPanel.cc
SmiVRecord.cc
SmiVSettings.cc
apply_daylight_arom_model_to_oemol.cc
build_time.cc)

set(SMIV_INCS
SmiV.H
SmivDataTable.H
SmiVFindMoleculeDialog.H
SmiVSettings.H
SmiVPanel.H
SmiVRecord.H)

set(SMIV_DACLIB_SRCS
QTMolDisplay2D.cc
QTSmartsEditDialog.cc
QTSmartsIntPickDialog.cc
apply_daylight_arom_model_to_oemol.cc
check_oechem_licence.cc
create_oesubsearch.cc
draw_oemol_to_qimage.cc
extract_smarts_from_smirks.cc
read_smarts_file.cc
split_smiles_into_atom_bits.cc
QT4SelectItems.cc
QTSmilesEditDialog.cc)

set(SMIV_DACLIB_INCS
DACOEMolAtomIndex.H
DACOEMolBondIndex.H
FileExceptions.H
QTMolDisplay2D.H
QT4SelectItems.H
QTSmartsEditDialog.H
QTSmartsIntPickDialog.H
QTSmilesEditDialog.H
SMARTSExceptions.H
stddefs.H)

set(LIBS ${LIBS}
  ${OEToolkits_LIBRARIES}
  ${Boost_LIBRARIES}
  ${CAIRO_LIBRARY})

set(EXECUTABLE_OUTPUT_PATH ${SMIV_SOURCE_DIR}/exe_${CMAKE_BUILD_TYPE})

add_executable(smiv ${SMIV_SRCS} ${SMIV_DACLIB_SRCS}
  ${SMIV_INCS} ${SMIV_DACLIB_INCS})
target_link_libraries(smiv z ${LIBS} Qt5::Widgets z pthread)
