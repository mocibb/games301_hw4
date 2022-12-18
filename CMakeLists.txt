cmake_minimum_required(VERSION 3.7)
project(SurfaceFrameworkCmake)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE "Release")

#设置不生成MinSizeRel;RelWithDebInfo
set(CMAKE_CONFIGURATION_TYPES "Debug;Release")

add_definitions(/D_USE_MATH_DEFINES)

#complier
# if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
#     message(STATUS "编译器是gcc")
# elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
#     message(STATUS "编译器是msvc")
# elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
#     message(STATUS "编译器是AppleClang")
#     set(Qt5Core_DIR /opt/Qt/5.4/gcc_64/lib/cmake/Qt5Widgets)#设置QT的目录
#     set(Qt5Gui_DIR /opt/Qt/5.4/gcc_64/lib/cmake/Qt5Widgets)
#     set(Qt5OpenGL_DIR /opt/Qt/5.4/gcc_64/lib/cmake/Qt5Widgets)
#     set(Qt5Widgets_DIR /opt/Qt/5.4/gcc_64/lib/cmake/Qt5Widgets)
# endif()

#mac
# cmake .. -DCMAKE_PREFIX_PATH=/opt/homebrew/opt/qt5/

#QT
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)
find_package(Qt5Widgets REQUIRED)  
find_package(Qt5Core REQUIRED)
find_package(Qt5Gui REQUIRED)
find_package(Qt5OpenGL REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories(${OPENGL_INCLUDE_DIR})
include_directories(${GLUT_INCLUDE_DIR})
include_directories(/usr/include/eigen3/)
include_directories(/opt/homebrew/include/eigen3)
include_directories(/opt/homebrew/include/)
include_directories(src)


set(QRC_SOURCE_FILES src/surfacemeshprocessing.qrc)
qt5_add_resources(QRC_FILES ${QRC_SOURCE_FILES})

file (GLOB GUI_FILE
    src/main.cpp
    src/MeshViewer/MainViewerWidget.*
    src/MeshViewer/MeshParamWidget.*
    src/surfacemeshprocessing.h
    src/surfacemeshprocessing.cpp
)

file (GLOB MeshViewer_FILE
    src/MeshViewer/InteractiveViewerWidget.*
    src/MeshViewer/MeshViewerWidget.*
    src/MeshViewer/QGLViewerWidget.*
)

file (GLOB PolyMesh_FILE
    src/PolyMesh/IOManger.cpp
    src/PolyMesh/PolyMesh.cpp
    src/PolyMesh/PolyMesh_Base.cpp
    src/PolyMesh/PolyMeshIterators.cpp
    src/include/Math/Algorithms.h
    src/include/Math/Massage.h
    src/include/Math/MPoint3.h
    src/include/Math/Mrgb.h
    src/include/Math/MVector3.h
    src/include/PolyMesh/IOManger.h
    src/include/PolyMesh/MemoryPool.h
    src/include/PolyMesh/PolyMesh.h
    src/include/PolyMesh/PolyMesh_Base.h
    src/include/PolyMesh/PolyMeshIterators.h
)

file (GLOB TUTTE_FILE
    src/Tutte/Tutte.cpp
)

file (GLOB AES_FILE
    src/AES/AES.cpp
)

file (GLOB BFF_FILE
    src/BFF/BFF.cpp
)

add_executable(${PROJECT_NAME}
    ${GUI_FILE}
    ${MeshViewer_FILE}
    ${PolyMesh_FILE}
    ${TUTTE_FILE}
    ${AES_FILE}
    ${BFF_FILE}
    ${QRC_FILES}
)

source_group(GUI FILES ${GUI_FILE})
source_group(MeshViewer FILES ${MeshViewer_FILE})
source_group(PolyMesh FILES ${PolyMesh_FILE})
source_group(TUTTE FILES ${TUTTE_FILE})
source_group(AES FILES ${AES_FILE})
source_group(BFF FILES ${BFF_FILE})

target_link_libraries(${PROJECT_NAME} Qt5::Widgets Qt5::Core Qt5::Gui Qt5::OpenGL)

target_link_libraries(${PROJECT_NAME} 
    ${GLUT_LIBRARY} 
    ${OPENGL_LIBRARY}
)

