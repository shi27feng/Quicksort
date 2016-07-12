project(Quicksort)

cmake_minimum_required(VERSION 2.8)

enable_testing()

set(bitness "x86")
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  set(bitness "x64")
endif()

add_definitions(-D_CRT_SECURE_NO_WARNINGS)
include_directories($ENV{INTELOCLSDKROOT}/include)
link_directories($ENV{INTELOCLSDKROOT}/lib/${bitness})
add_executable(Quicksort QuicksortMain.cpp ImageUtils.h OpenCLUtils.h Quicksort.h QuicksortKernels.cl)
target_link_libraries(Quicksort OpenCL)

add_test(NAME Quicksort_smoke
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
         COMMAND Quicksort 100 gpu intel 2048  2048  how_cl)
