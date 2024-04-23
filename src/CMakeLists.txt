# src/CMakeLists.txt

file(GLOB_RECURSE SRC_FILES CONFIGURE_DEPENDS
  ${CMAKE_CURRENT_SOURCE_DIR}/*.cc
)

add_library(FDoG SHARED ${SRC_FILES})

target_include_directories(FDoG PUBLIC ${CMAKE_SOURCE_DIR}/include)

target_compile_options(FDoG PRIVATE)