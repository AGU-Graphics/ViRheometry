include( CMakeSourceFiles.txt )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(CMAKE_CXX_STANDARD 17)
add_definitions( -std=c++17 -Wall )
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

string( TIMESTAMP BUILD_DATE_TIME "%Y-%m-%d %H:%M:%S" )
configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/CompileDefinitions.h.in" "${CMAKE_CURRENT_SOURCE_DIR}/CompileDefinitions.h" @ONLY )

add_executable( GLRender3d ${Headers} ${Sources} )

if( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" )
  target_compile_definitions( GLRender3d PUBLIC CMAKE_DETECTED_GNU_COMPILER )
elseif( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel" )
  target_compile_definitions( GLRender3d PUBLIC CMAKE_DETECTED_INTEL_COMPILER )
elseif( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" )
  target_compile_definitions( GLRender3d PUBLIC CMAKE_DETECTED_CLANG_COMPILER )
elseif( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang" )
  target_compile_definitions( GLRender3d PUBLIC CMAKE_DETECTED_CLANG_COMPILER )
else()
  message( FATAL_ERROR "Unsupported C++ compiler ${CMAKE_CXX_COMPILER_ID}. Please email the maintainer." )
endif()

find_package( Eigen3 3.2.9 REQUIRED )
target_include_directories( GLRender3d SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIR} )
target_compile_definitions( GLRender3d PUBLIC EIGEN_DONT_PARALLELIZE )
if( CMAKE_BUILD_TYPE STREQUAL "Debug" )
  target_compile_definitions( GLRender3d PUBLIC EIGEN_DONT_ALIGN EIGEN_INITIALIZE_MATRICES_BY_NAN )
endif()

find_package( OpenCV REQUIRED )
target_include_directories( GLRender3d SYSTEM PRIVATE ${OpenCV_INCLUDE_DIRS} )
target_link_libraries( GLRender3d ${OpenCV_LIBS} )

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} )
target_link_libraries(GLRender3d ${OPENGL_LIBRARIES} ${GLUT_glut_LIBRARY} )

