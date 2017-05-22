cmake_minimum_required(VERSION 3.7)
project(unfold)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        common.h
        craft.cpp
        dual.cpp
        feature.cpp
        gene.cpp
        gene.h
        main.cpp
        matching.cpp
        meshui.cpp
        mst.cpp
        project.cpp
        random.cpp
        random.h
        render.cpp
        scan.h
        search.cpp
        sheetui.cpp
        stitch.cpp
        timer.cpp
        timer.h
        ui.h)

find_package( OpenCV )
find_package(CGAL REQUIRED)
FIND_PACKAGE(OpenGL REQUIRED)
FIND_PACKAGE(GLUT REQUIRED)
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
include( ${CGAL_USE_FILE} )
INCLUDE_DIRECTORIES(${OPENGL_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${GLUT_INCLUDE_DIR})

add_executable(unfold ${SOURCE_FILES})

target_link_libraries( unfold ${OpenCV_LIBS} )
target_link_libraries( unfold ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
TARGET_LINK_LIBRARIES(unfold ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
target_link_libraries( unfold ${GSL_LIBRARIES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
