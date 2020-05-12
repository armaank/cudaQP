# Minimum version required
cmake_minimum_required (VERSION 3.2)

# Project name
project (osqp)

# Export compile commands
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set the output folder where your program will be created
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/out)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/out)

# Some non-standard CMake modules
LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/configure/cmake)
INCLUDE(Utils)

# Set options
# ----------------------------------------------
option (COVERAGE "Perform code coverage" OFF)
# message(STATUS "Code coverage is ${COVERAGE}")
# Are unittests generated?
# option (UNITTESTS "Enable unittests generation" OFF)

# Is the code generated for embedded platforms?
#   1 :   Yes. Matrix update not allowed.
# #   2 :   Yes. Matrix update allowed.
# if (NOT DEFINED EMBEDDED)
#     message(STATUS "Embedded is OFF")
# else()
#     message(STATUS "Embedded is ${EMBEDDED}")
# endif()

# Is the CUDA support enabled?
option (CUDA_SUPPORT "Enable CUDA support" OFF)
message(STATUS "CUDA support is ${CUDA_SUPPORT}")

# Is printing enabled?
option (PRINTING "Enable solver printing" ON)
# if (DEFINED EMBEDDED)
#     message(STATUS "Disabling printing for embedded")
#     set(PRINTING OFF)
# endif()
message(STATUS "Printing is ${PRINTING}")


# Is profiling enabled?
option (PROFILING "Enable solver profiling (timing)" ON)
# if (DEFINED EMBEDDED)
#     message(STATUS "Disabling profiling for embedded")
#     set(PROFILING OFF)
# endif()
message(STATUS "Profiling is ${PROFILING}")

# Is user interrupt enabled?
option (CTRLC "Enable user interrupt (Ctrl-C)" ON)
# if (DEFINED EMBEDDED)
#     message(STATUS "Disabling user interrupt for embedded")
#     set(CTRLC OFF)
# endif()
message(STATUS "User interrupt is ${CTRLC}")

# Use floats instead of doubles
if (CUDA_SUPPORT)
    # Use floats with CUDA by default
    option (DFLOAT "Use float numbers instead of doubles" ON)
else ()
    option (DFLOAT "Use float numbers instead of doubles" OFF)
endif()
message(STATUS "Floats are ${DFLOAT}")

# Use long integers for indexing
option (DLONG "Use long integers (64bit) for indexing" ON)
if (NOT (CMAKE_SIZEOF_VOID_P EQUAL 8))
    message(STATUS "Disabling long integers (64bit) on 32bit machine")
    set(DLONG OFF)
elseif (CUDA_SUPPORT)
	message(STATUS "Disabling long integers (64bit) for CUDA")
	set(DLONG OFF)
endif()
message(STATUS "Long integers (64bit) are ${DLONG}")


option (DEBUG "Debug mode" OFF)
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
	set (DEBUG ON)
	message(STATUS "Debug mode is ${DEBUG}")
endif()

# Add code coverage
# option (COVERAGE "Perform code coverage" OFF)
# message(STATUS "Code coverage is ${COVERAGE}")


# CUDA support
# ----------------------------------------------
if (CUDA_SUPPORT)

    if (NOT CMAKE_CUDA_COMPILER)
        set(CMAKE_CUDA_COMPILER "/usr/local/cuda/bin/nvcc")
    endif()

    enable_language(CUDA)
    find_package(CUDA)

    # cublas and cusparse headers
    include_directories(${CUDA_INCLUDE_DIRS})

    # sm_75 -> Turing support
    if (DFLOAT)
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --gpu-architecture=compute_30 --gpu-code=compute_30,sm_75")
    else()
        # To use doubles we need compute capability 6.0 for atomic operations
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --gpu-architecture=compute_60 --gpu-code=compute_60,sm_75")
    endif()
endif()


# Memory allocators
# ----------------------------------------------

#Report on custom user header options.  This is intended to allow
#users to provide definitions of their own memory functions
# The header should define the functions as follows
#
# define c_malloc mymalloc
# define c_calloc mycalloc
# define c_realloc myrealloc
# define c_free myfree

if(OSQP_CUSTOM_MEMORY)
	message(STATUS "User custom memory management header: ${OSQP_CUSTOM_MEMORY}")
endif()

# Generate header file with the global options
# ---------------------------------------------
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/configure/osqp_configure.h.in
               ${CMAKE_CURRENT_SOURCE_DIR}/include/osqp_configure.h
               NEWLINE_STYLE LF)

# Set Compiler flags
# ----------------------------------------------
set(CMAKE_POSITION_INDEPENDENT_CODE ON)  # -fPIC
# if (COVERAGE)
#     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
# if(FORTRAN)
#     set(CMAKE_FORTRAN_FLAGS "${CMAKE_FORTRAN_FLAGS} --coverage")
# endif(FORTRAN)
# endif()

# if (DEBUG)
#     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g")
# else()
#     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
# endif()

set(CMAKE_C_STANDARD_LIBRARIES "${CMAKE_C_STANDARD_LIBRARIES} -lm")      # Include math
# Include real time library in linux
# if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
set(CMAKE_C_STANDARD_LIBRARIES "${CMAKE_C_STANDARD_LIBRARIES} -lrt -ldl")
# endif()

# Set OSQP core sources and includes
# ----------------------------------------------
add_subdirectory (src)
add_subdirectory (include)

#location of the osqp api includes
set(osqp_api_includes ${CMAKE_CURRENT_SOURCE_DIR}/include/)

# Algebra implementation and CSC support
# ----------------------------------------------
add_subdirectory(algebra)
set(osqp_objects $<TARGET_OBJECTS:osqp_algebra>)

# Create Static Library
# ----------------------------------------------

# Add linear system solvers cumulative library
add_subdirectory(lin_sys)

# Static library
add_library (osqpstatic STATIC ${osqp_src} ${osqp_objects} ${osqp_headers} ${linsys_solvers})

# Give same name to static library output
set_target_properties(osqpstatic PROPERTIES OUTPUT_NAME osqp)

# Include directories for linear system solvers
target_include_directories(osqpstatic PRIVATE
                           ${linsys_solvers_includes}
                           ${algebra_impl_includes}
                           ${csc_tools_includes})


# Declare include directories for the cmake exported target
target_include_directories(osqpstatic
                           PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
                                  "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}/osqp>")

# Install Static Library
# ----------------------------------------------

include(GNUInstallDirs)

install(TARGETS osqpstatic
        EXPORT  ${PROJECT_NAME}
        ARCHIVE       DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        LIBRARY       DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        RUNTIME       DESTINATION "${CMAKE_INSTALL_BINDIR}")


# Install Headers
# ----------------------------------------------

install(FILES ${osqp_headers} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/osqp")


if (CUDA_SUPPORT)
    target_link_libraries (osqpstatic cublas cusparse)
endif()

add_library (osqp SHARED ${osqp_src} ${osqp_headers} ${osqp_objects} ${linsys_solvers})

# Include directories
target_include_directories(osqp PRIVATE
                            ${linsys_solvers_includes}
                            ${csc_tools_includes}
                            ${algebra_impl_includes}
)

# Declare include directories for the cmake exported target
target_include_directories(osqp
                            PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
                                    "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}/osqp>")

if (CUDA_SUPPORT AND MSVC)
    # Visual Studio needs the cuda libaries to be linked as well
    target_link_libraries (osqp cublas cusparse)
endif (CUDA_SUPPORT AND MSVC)

# Install osqp shared library
install(TARGETS osqp
        EXPORT  ${PROJECT_NAME}
        LIBRARY       DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        ARCHIVE       DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        RUNTIME       DESTINATION "${CMAKE_INSTALL_BINDIR}")

# Create demo executable (linked to static library)
add_executable (osqp_demo ${PROJECT_SOURCE_DIR}/examples/osqp_demo.c)
target_link_libraries (osqp_demo osqpstatic)

# Create CMake packages for the build directory
# ----------------------------------------------
include(CMakePackageConfigHelpers)

export(EXPORT ${PROJECT_NAME}
  FILE "${CMAKE_CURRENT_BINARY_DIR}/osqp-targets.cmake"
  NAMESPACE osqp::)

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/osqp-config.cmake)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/osqp-config.cmake "include(\"\${CMAKE_CURRENT_LIST_DIR}/osqp-targets.cmake\")\n")
endif()


# Create CMake packages for the install directory
# ----------------------------------------------

set(ConfigPackageLocation lib/cmake/osqp)

install(EXPORT ${PROJECT_NAME}
        FILE osqp-targets.cmake
        NAMESPACE osqp::
        DESTINATION ${ConfigPackageLocation})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/osqp-config.cmake
        DESTINATION ${ConfigPackageLocation})



# Add uninstall command
# ----------------------------------------------
if(NOT TARGET uninstall)
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/configure/cmake/cmake_uninstall.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
        IMMEDIATE @ONLY)

    add_custom_target(uninstall
        COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()

