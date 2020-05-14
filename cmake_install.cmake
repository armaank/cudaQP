# Install script for directory: /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/out/libosqp.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/osqp" TYPE FILE FILES
    "/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/include/algebra_vector.h"
    "/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/include/algebra_matrix.h"
    "/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/include/auxil.h"
    "/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/include/csc_type.h"
    "/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/include/error.h"
    "/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/include/glob_opts.h"
    "/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/include/lin_alg.h"
    "/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/include/osqp.h"
    "/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/include/osqp_api_constants.h"
    "/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/include/osqp_api_functions.h"
    "/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/include/osqp_api_types.h"
    "/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/include/osqp_configure.h"
    "/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/include/proj.h"
    "/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/include/scaling.h"
    "/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/include/types.h"
    "/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/include/util.h"
    "/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/include/polish.h"
    "/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/include/lin_sys.h"
    "/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/include/ctrlc.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libosqp.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libosqp.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libosqp.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/out/libosqp.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libosqp.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libosqp.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libosqp.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/osqp/osqp-targets.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/osqp/osqp-targets.cmake"
         "/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/CMakeFiles/Export/lib/cmake/osqp/osqp-targets.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/osqp/osqp-targets-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/osqp/osqp-targets.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/osqp" TYPE FILE FILES "/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/CMakeFiles/Export/lib/cmake/osqp/osqp-targets.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^()$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/osqp" TYPE FILE FILES "/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/CMakeFiles/Export/lib/cmake/osqp/osqp-targets-noconfig.cmake")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/osqp" TYPE FILE FILES "/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/osqp-config.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/src/cmake_install.cmake")
  include("/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/include/cmake_install.cmake")
  include("/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/algebra/cmake_install.cmake")
  include("/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/cmake_install.cmake")
  include("/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/svm/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
