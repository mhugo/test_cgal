cmake_minimum_required(VERSION 2.6)

include( cmake/FindCGAL.cmake )

find_package(Qt4 REQUIRED)
find_package(CGAL REQUIRED)

include(${QT_USE_FILE})

include_directories( ${CGAL_INCLUDE_DIRS} ${QT_INCLUDE_DIRS} )
link_directories( ${CGAL_LIBRARY_DIRS} ${QT_LIBRARY_DIRS} )
add_definitions( -frounding-math )

add_executable( tracer src/ray.cc )
target_link_libraries( tracer ${CGAL_LIBRARIES} ${QT_LIBRARIES} )

add_executable( poly src/poly.cc )
target_link_libraries( poly ${CGAL_LIBRARIES} )
