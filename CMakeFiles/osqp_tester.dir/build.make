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
include CMakeFiles/osqp_tester.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/osqp_tester.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/osqp_tester.dir/flags.make

svm/basic_qp/data.h: svm/generate_tests_data.py
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating unittests data files using Python"
	cd /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/svm && /usr/bin/python generate_tests_data.py

svm/basic_svm/data.h: svm/basic_qp/data.h
	@$(CMAKE_COMMAND) -E touch_nocreate svm/basic_svm/data.h

svm/lin_alg/data.h: svm/basic_qp/data.h
	@$(CMAKE_COMMAND) -E touch_nocreate svm/lin_alg/data.h

svm/non_cvx/data.h: svm/basic_qp/data.h
	@$(CMAKE_COMMAND) -E touch_nocreate svm/non_cvx/data.h

svm/primal_dual_infeasibility/data.h: svm/basic_qp/data.h
	@$(CMAKE_COMMAND) -E touch_nocreate svm/primal_dual_infeasibility/data.h

svm/primal_infeasibility/data.h: svm/basic_qp/data.h
	@$(CMAKE_COMMAND) -E touch_nocreate svm/primal_infeasibility/data.h

svm/solve_linsys/data.h: svm/basic_qp/data.h
	@$(CMAKE_COMMAND) -E touch_nocreate svm/solve_linsys/data.h

svm/unconstrained/data.h: svm/basic_qp/data.h
	@$(CMAKE_COMMAND) -E touch_nocreate svm/unconstrained/data.h

svm/update_matrices/data.h: svm/basic_qp/data.h
	@$(CMAKE_COMMAND) -E touch_nocreate svm/update_matrices/data.h

CMakeFiles/osqp_tester.dir/svm/osqp_tester.c.o: CMakeFiles/osqp_tester.dir/flags.make
CMakeFiles/osqp_tester.dir/svm/osqp_tester.c.o: svm/osqp_tester.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/osqp_tester.dir/svm/osqp_tester.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/osqp_tester.dir/svm/osqp_tester.c.o   -c /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/svm/osqp_tester.c

CMakeFiles/osqp_tester.dir/svm/osqp_tester.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/osqp_tester.dir/svm/osqp_tester.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/svm/osqp_tester.c > CMakeFiles/osqp_tester.dir/svm/osqp_tester.c.i

CMakeFiles/osqp_tester.dir/svm/osqp_tester.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/osqp_tester.dir/svm/osqp_tester.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/svm/osqp_tester.c -o CMakeFiles/osqp_tester.dir/svm/osqp_tester.c.s

CMakeFiles/osqp_tester.dir/svm/osqp_tester.c.o.requires:

.PHONY : CMakeFiles/osqp_tester.dir/svm/osqp_tester.c.o.requires

CMakeFiles/osqp_tester.dir/svm/osqp_tester.c.o.provides: CMakeFiles/osqp_tester.dir/svm/osqp_tester.c.o.requires
	$(MAKE) -f CMakeFiles/osqp_tester.dir/build.make CMakeFiles/osqp_tester.dir/svm/osqp_tester.c.o.provides.build
.PHONY : CMakeFiles/osqp_tester.dir/svm/osqp_tester.c.o.provides

CMakeFiles/osqp_tester.dir/svm/osqp_tester.c.o.provides.build: CMakeFiles/osqp_tester.dir/svm/osqp_tester.c.o


# Object files for target osqp_tester
osqp_tester_OBJECTS = \
"CMakeFiles/osqp_tester.dir/svm/osqp_tester.c.o"

# External object files for target osqp_tester
osqp_tester_EXTERNAL_OBJECTS =

out/osqp_tester: CMakeFiles/osqp_tester.dir/svm/osqp_tester.c.o
out/osqp_tester: CMakeFiles/osqp_tester.dir/build.make
out/osqp_tester: out/libosqp.a
out/osqp_tester: CMakeFiles/osqp_tester.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C executable out/osqp_tester"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/osqp_tester.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/osqp_tester.dir/build: out/osqp_tester

.PHONY : CMakeFiles/osqp_tester.dir/build

CMakeFiles/osqp_tester.dir/requires: CMakeFiles/osqp_tester.dir/svm/osqp_tester.c.o.requires

.PHONY : CMakeFiles/osqp_tester.dir/requires

CMakeFiles/osqp_tester.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/osqp_tester.dir/cmake_clean.cmake
.PHONY : CMakeFiles/osqp_tester.dir/clean

CMakeFiles/osqp_tester.dir/depend: svm/basic_qp/data.h
CMakeFiles/osqp_tester.dir/depend: svm/basic_svm/data.h
CMakeFiles/osqp_tester.dir/depend: svm/lin_alg/data.h
CMakeFiles/osqp_tester.dir/depend: svm/non_cvx/data.h
CMakeFiles/osqp_tester.dir/depend: svm/primal_dual_infeasibility/data.h
CMakeFiles/osqp_tester.dir/depend: svm/primal_infeasibility/data.h
CMakeFiles/osqp_tester.dir/depend: svm/solve_linsys/data.h
CMakeFiles/osqp_tester.dir/depend: svm/unconstrained/data.h
CMakeFiles/osqp_tester.dir/depend: svm/update_matrices/data.h
	cd /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/CMakeFiles/osqp_tester.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/osqp_tester.dir/depend
