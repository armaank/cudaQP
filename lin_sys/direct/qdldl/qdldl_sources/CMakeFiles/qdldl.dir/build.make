# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP

# Include any dependencies generated for this target.
include lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/depend.make

# Include the progress variables for this target.
include lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/progress.make

# Include the compile flags for this target's objects.
include lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/flags.make

lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/src/qdldl.c.o: lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/flags.make
lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/src/qdldl.c.o: lin_sys/direct/qdldl/qdldl_sources/src/qdldl.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/src/qdldl.c.o"
	cd /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/direct/qdldl/qdldl_sources && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/qdldl.dir/src/qdldl.c.o   -c /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/direct/qdldl/qdldl_sources/src/qdldl.c

lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/src/qdldl.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/qdldl.dir/src/qdldl.c.i"
	cd /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/direct/qdldl/qdldl_sources && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/direct/qdldl/qdldl_sources/src/qdldl.c > CMakeFiles/qdldl.dir/src/qdldl.c.i

lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/src/qdldl.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/qdldl.dir/src/qdldl.c.s"
	cd /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/direct/qdldl/qdldl_sources && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/direct/qdldl/qdldl_sources/src/qdldl.c -o CMakeFiles/qdldl.dir/src/qdldl.c.s

lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/src/qdldl.c.o.requires:

.PHONY : lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/src/qdldl.c.o.requires

lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/src/qdldl.c.o.provides: lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/src/qdldl.c.o.requires
	$(MAKE) -f lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/build.make lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/src/qdldl.c.o.provides.build
.PHONY : lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/src/qdldl.c.o.provides

lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/src/qdldl.c.o.provides.build: lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/src/qdldl.c.o


# Object files for target qdldl
qdldl_OBJECTS = \
"CMakeFiles/qdldl.dir/src/qdldl.c.o"

# External object files for target qdldl
qdldl_EXTERNAL_OBJECTS =

lin_sys/direct/qdldl/qdldl_sources/out/libqdldl.so: lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/src/qdldl.c.o
lin_sys/direct/qdldl/qdldl_sources/out/libqdldl.so: lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/build.make
lin_sys/direct/qdldl/qdldl_sources/out/libqdldl.so: lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C shared library out/libqdldl.so"
	cd /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/direct/qdldl/qdldl_sources && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/qdldl.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/build: lin_sys/direct/qdldl/qdldl_sources/out/libqdldl.so

.PHONY : lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/build

lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/requires: lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/src/qdldl.c.o.requires

.PHONY : lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/requires

lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/clean:
	cd /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/direct/qdldl/qdldl_sources && $(CMAKE_COMMAND) -P CMakeFiles/qdldl.dir/cmake_clean.cmake
.PHONY : lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/clean

lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/depend:
	cd /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/direct/qdldl/qdldl_sources /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/direct/qdldl/qdldl_sources /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldl.dir/depend

