cmake_minimum_required(VERSION 3.10)

project(
    FDoG_app
    VERSION 0.1
    DESCRIPTION "Practive Project"
    LANGUAGES CXX
)

find_package( OpenCV REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )

add_subdirectory(src)

add_executable (${PROJECT_NAME} main.cc)
target_link_libraries(${PROJECT_NAME} 
  FDoG
  ${OpenCV_LIBS}
)
add_executable (getCL main2.cc)
target_link_libraries(getCL 
  FDoG
  ${OpenCV_LIBS}
)
add_executable (getFBL main3.cc)
target_link_libraries(getFBL 
  FDoG
  ${OpenCV_LIBS}
)

set(CMAKE_CXX_FLAGS "-Wall -g")

target_compile_options(${PROJECT_NAME} PRIVATE)
target_compile_options(getCL PRIVATE) 
target_compile_options(getFBL PRIVATE) 