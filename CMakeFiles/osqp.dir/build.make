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
include CMakeFiles/osqp.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/osqp.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/osqp.dir/flags.make

CMakeFiles/osqp.dir/src/auxil.c.o: CMakeFiles/osqp.dir/flags.make
CMakeFiles/osqp.dir/src/auxil.c.o: src/auxil.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/osqp.dir/src/auxil.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/osqp.dir/src/auxil.c.o   -c /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/src/auxil.c

CMakeFiles/osqp.dir/src/auxil.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/osqp.dir/src/auxil.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/src/auxil.c > CMakeFiles/osqp.dir/src/auxil.c.i

CMakeFiles/osqp.dir/src/auxil.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/osqp.dir/src/auxil.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/src/auxil.c -o CMakeFiles/osqp.dir/src/auxil.c.s

CMakeFiles/osqp.dir/src/auxil.c.o.requires:

.PHONY : CMakeFiles/osqp.dir/src/auxil.c.o.requires

CMakeFiles/osqp.dir/src/auxil.c.o.provides: CMakeFiles/osqp.dir/src/auxil.c.o.requires
	$(MAKE) -f CMakeFiles/osqp.dir/build.make CMakeFiles/osqp.dir/src/auxil.c.o.provides.build
.PHONY : CMakeFiles/osqp.dir/src/auxil.c.o.provides

CMakeFiles/osqp.dir/src/auxil.c.o.provides.build: CMakeFiles/osqp.dir/src/auxil.c.o


CMakeFiles/osqp.dir/src/error.c.o: CMakeFiles/osqp.dir/flags.make
CMakeFiles/osqp.dir/src/error.c.o: src/error.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/osqp.dir/src/error.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/osqp.dir/src/error.c.o   -c /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/src/error.c

CMakeFiles/osqp.dir/src/error.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/osqp.dir/src/error.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/src/error.c > CMakeFiles/osqp.dir/src/error.c.i

CMakeFiles/osqp.dir/src/error.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/osqp.dir/src/error.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/src/error.c -o CMakeFiles/osqp.dir/src/error.c.s

CMakeFiles/osqp.dir/src/error.c.o.requires:

.PHONY : CMakeFiles/osqp.dir/src/error.c.o.requires

CMakeFiles/osqp.dir/src/error.c.o.provides: CMakeFiles/osqp.dir/src/error.c.o.requires
	$(MAKE) -f CMakeFiles/osqp.dir/build.make CMakeFiles/osqp.dir/src/error.c.o.provides.build
.PHONY : CMakeFiles/osqp.dir/src/error.c.o.provides

CMakeFiles/osqp.dir/src/error.c.o.provides.build: CMakeFiles/osqp.dir/src/error.c.o


CMakeFiles/osqp.dir/src/osqp_api.c.o: CMakeFiles/osqp.dir/flags.make
CMakeFiles/osqp.dir/src/osqp_api.c.o: src/osqp_api.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object CMakeFiles/osqp.dir/src/osqp_api.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/osqp.dir/src/osqp_api.c.o   -c /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/src/osqp_api.c

CMakeFiles/osqp.dir/src/osqp_api.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/osqp.dir/src/osqp_api.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/src/osqp_api.c > CMakeFiles/osqp.dir/src/osqp_api.c.i

CMakeFiles/osqp.dir/src/osqp_api.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/osqp.dir/src/osqp_api.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/src/osqp_api.c -o CMakeFiles/osqp.dir/src/osqp_api.c.s

CMakeFiles/osqp.dir/src/osqp_api.c.o.requires:

.PHONY : CMakeFiles/osqp.dir/src/osqp_api.c.o.requires

CMakeFiles/osqp.dir/src/osqp_api.c.o.provides: CMakeFiles/osqp.dir/src/osqp_api.c.o.requires
	$(MAKE) -f CMakeFiles/osqp.dir/build.make CMakeFiles/osqp.dir/src/osqp_api.c.o.provides.build
.PHONY : CMakeFiles/osqp.dir/src/osqp_api.c.o.provides

CMakeFiles/osqp.dir/src/osqp_api.c.o.provides.build: CMakeFiles/osqp.dir/src/osqp_api.c.o


CMakeFiles/osqp.dir/src/proj.c.o: CMakeFiles/osqp.dir/flags.make
CMakeFiles/osqp.dir/src/proj.c.o: src/proj.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object CMakeFiles/osqp.dir/src/proj.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/osqp.dir/src/proj.c.o   -c /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/src/proj.c

CMakeFiles/osqp.dir/src/proj.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/osqp.dir/src/proj.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/src/proj.c > CMakeFiles/osqp.dir/src/proj.c.i

CMakeFiles/osqp.dir/src/proj.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/osqp.dir/src/proj.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/src/proj.c -o CMakeFiles/osqp.dir/src/proj.c.s

CMakeFiles/osqp.dir/src/proj.c.o.requires:

.PHONY : CMakeFiles/osqp.dir/src/proj.c.o.requires

CMakeFiles/osqp.dir/src/proj.c.o.provides: CMakeFiles/osqp.dir/src/proj.c.o.requires
	$(MAKE) -f CMakeFiles/osqp.dir/build.make CMakeFiles/osqp.dir/src/proj.c.o.provides.build
.PHONY : CMakeFiles/osqp.dir/src/proj.c.o.provides

CMakeFiles/osqp.dir/src/proj.c.o.provides.build: CMakeFiles/osqp.dir/src/proj.c.o


CMakeFiles/osqp.dir/src/scaling.c.o: CMakeFiles/osqp.dir/flags.make
CMakeFiles/osqp.dir/src/scaling.c.o: src/scaling.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object CMakeFiles/osqp.dir/src/scaling.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/osqp.dir/src/scaling.c.o   -c /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/src/scaling.c

CMakeFiles/osqp.dir/src/scaling.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/osqp.dir/src/scaling.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/src/scaling.c > CMakeFiles/osqp.dir/src/scaling.c.i

CMakeFiles/osqp.dir/src/scaling.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/osqp.dir/src/scaling.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/src/scaling.c -o CMakeFiles/osqp.dir/src/scaling.c.s

CMakeFiles/osqp.dir/src/scaling.c.o.requires:

.PHONY : CMakeFiles/osqp.dir/src/scaling.c.o.requires

CMakeFiles/osqp.dir/src/scaling.c.o.provides: CMakeFiles/osqp.dir/src/scaling.c.o.requires
	$(MAKE) -f CMakeFiles/osqp.dir/build.make CMakeFiles/osqp.dir/src/scaling.c.o.provides.build
.PHONY : CMakeFiles/osqp.dir/src/scaling.c.o.provides

CMakeFiles/osqp.dir/src/scaling.c.o.provides.build: CMakeFiles/osqp.dir/src/scaling.c.o


CMakeFiles/osqp.dir/src/util.c.o: CMakeFiles/osqp.dir/flags.make
CMakeFiles/osqp.dir/src/util.c.o: src/util.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building C object CMakeFiles/osqp.dir/src/util.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/osqp.dir/src/util.c.o   -c /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/src/util.c

CMakeFiles/osqp.dir/src/util.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/osqp.dir/src/util.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/src/util.c > CMakeFiles/osqp.dir/src/util.c.i

CMakeFiles/osqp.dir/src/util.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/osqp.dir/src/util.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/src/util.c -o CMakeFiles/osqp.dir/src/util.c.s

CMakeFiles/osqp.dir/src/util.c.o.requires:

.PHONY : CMakeFiles/osqp.dir/src/util.c.o.requires

CMakeFiles/osqp.dir/src/util.c.o.provides: CMakeFiles/osqp.dir/src/util.c.o.requires
	$(MAKE) -f CMakeFiles/osqp.dir/build.make CMakeFiles/osqp.dir/src/util.c.o.provides.build
.PHONY : CMakeFiles/osqp.dir/src/util.c.o.provides

CMakeFiles/osqp.dir/src/util.c.o.provides.build: CMakeFiles/osqp.dir/src/util.c.o


CMakeFiles/osqp.dir/src/polish.c.o: CMakeFiles/osqp.dir/flags.make
CMakeFiles/osqp.dir/src/polish.c.o: src/polish.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building C object CMakeFiles/osqp.dir/src/polish.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/osqp.dir/src/polish.c.o   -c /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/src/polish.c

CMakeFiles/osqp.dir/src/polish.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/osqp.dir/src/polish.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/src/polish.c > CMakeFiles/osqp.dir/src/polish.c.i

CMakeFiles/osqp.dir/src/polish.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/osqp.dir/src/polish.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/src/polish.c -o CMakeFiles/osqp.dir/src/polish.c.s

CMakeFiles/osqp.dir/src/polish.c.o.requires:

.PHONY : CMakeFiles/osqp.dir/src/polish.c.o.requires

CMakeFiles/osqp.dir/src/polish.c.o.provides: CMakeFiles/osqp.dir/src/polish.c.o.requires
	$(MAKE) -f CMakeFiles/osqp.dir/build.make CMakeFiles/osqp.dir/src/polish.c.o.provides.build
.PHONY : CMakeFiles/osqp.dir/src/polish.c.o.provides

CMakeFiles/osqp.dir/src/polish.c.o.provides.build: CMakeFiles/osqp.dir/src/polish.c.o


CMakeFiles/osqp.dir/src/lin_sys.c.o: CMakeFiles/osqp.dir/flags.make
CMakeFiles/osqp.dir/src/lin_sys.c.o: src/lin_sys.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building C object CMakeFiles/osqp.dir/src/lin_sys.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/osqp.dir/src/lin_sys.c.o   -c /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/src/lin_sys.c

CMakeFiles/osqp.dir/src/lin_sys.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/osqp.dir/src/lin_sys.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/src/lin_sys.c > CMakeFiles/osqp.dir/src/lin_sys.c.i

CMakeFiles/osqp.dir/src/lin_sys.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/osqp.dir/src/lin_sys.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/src/lin_sys.c -o CMakeFiles/osqp.dir/src/lin_sys.c.s

CMakeFiles/osqp.dir/src/lin_sys.c.o.requires:

.PHONY : CMakeFiles/osqp.dir/src/lin_sys.c.o.requires

CMakeFiles/osqp.dir/src/lin_sys.c.o.provides: CMakeFiles/osqp.dir/src/lin_sys.c.o.requires
	$(MAKE) -f CMakeFiles/osqp.dir/build.make CMakeFiles/osqp.dir/src/lin_sys.c.o.provides.build
.PHONY : CMakeFiles/osqp.dir/src/lin_sys.c.o.provides

CMakeFiles/osqp.dir/src/lin_sys.c.o.provides.build: CMakeFiles/osqp.dir/src/lin_sys.c.o


CMakeFiles/osqp.dir/src/ctrlc.c.o: CMakeFiles/osqp.dir/flags.make
CMakeFiles/osqp.dir/src/ctrlc.c.o: src/ctrlc.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building C object CMakeFiles/osqp.dir/src/ctrlc.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/osqp.dir/src/ctrlc.c.o   -c /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/src/ctrlc.c

CMakeFiles/osqp.dir/src/ctrlc.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/osqp.dir/src/ctrlc.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/src/ctrlc.c > CMakeFiles/osqp.dir/src/ctrlc.c.i

CMakeFiles/osqp.dir/src/ctrlc.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/osqp.dir/src/ctrlc.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/src/ctrlc.c -o CMakeFiles/osqp.dir/src/ctrlc.c.s

CMakeFiles/osqp.dir/src/ctrlc.c.o.requires:

.PHONY : CMakeFiles/osqp.dir/src/ctrlc.c.o.requires

CMakeFiles/osqp.dir/src/ctrlc.c.o.provides: CMakeFiles/osqp.dir/src/ctrlc.c.o.requires
	$(MAKE) -f CMakeFiles/osqp.dir/build.make CMakeFiles/osqp.dir/src/ctrlc.c.o.provides.build
.PHONY : CMakeFiles/osqp.dir/src/ctrlc.c.o.provides

CMakeFiles/osqp.dir/src/ctrlc.c.o.provides.build: CMakeFiles/osqp.dir/src/ctrlc.c.o


CMakeFiles/osqp.dir/lin_sys/lib_handler.c.o: CMakeFiles/osqp.dir/flags.make
CMakeFiles/osqp.dir/lin_sys/lib_handler.c.o: lin_sys/lib_handler.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building C object CMakeFiles/osqp.dir/lin_sys/lib_handler.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/osqp.dir/lin_sys/lib_handler.c.o   -c /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/lib_handler.c

CMakeFiles/osqp.dir/lin_sys/lib_handler.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/osqp.dir/lin_sys/lib_handler.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/lib_handler.c > CMakeFiles/osqp.dir/lin_sys/lib_handler.c.i

CMakeFiles/osqp.dir/lin_sys/lib_handler.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/osqp.dir/lin_sys/lib_handler.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/lib_handler.c -o CMakeFiles/osqp.dir/lin_sys/lib_handler.c.s

CMakeFiles/osqp.dir/lin_sys/lib_handler.c.o.requires:

.PHONY : CMakeFiles/osqp.dir/lin_sys/lib_handler.c.o.requires

CMakeFiles/osqp.dir/lin_sys/lib_handler.c.o.provides: CMakeFiles/osqp.dir/lin_sys/lib_handler.c.o.requires
	$(MAKE) -f CMakeFiles/osqp.dir/build.make CMakeFiles/osqp.dir/lin_sys/lib_handler.c.o.provides.build
.PHONY : CMakeFiles/osqp.dir/lin_sys/lib_handler.c.o.provides

CMakeFiles/osqp.dir/lin_sys/lib_handler.c.o.provides.build: CMakeFiles/osqp.dir/lin_sys/lib_handler.c.o


# Object files for target osqp
osqp_OBJECTS = \
"CMakeFiles/osqp.dir/src/auxil.c.o" \
"CMakeFiles/osqp.dir/src/error.c.o" \
"CMakeFiles/osqp.dir/src/osqp_api.c.o" \
"CMakeFiles/osqp.dir/src/proj.c.o" \
"CMakeFiles/osqp.dir/src/scaling.c.o" \
"CMakeFiles/osqp.dir/src/util.c.o" \
"CMakeFiles/osqp.dir/src/polish.c.o" \
"CMakeFiles/osqp.dir/src/lin_sys.c.o" \
"CMakeFiles/osqp.dir/src/ctrlc.c.o" \
"CMakeFiles/osqp.dir/lin_sys/lib_handler.c.o"

# External object files for target osqp
osqp_EXTERNAL_OBJECTS = \
"/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/algebra/CMakeFiles/osqp_algebra.dir/default/algebra_libs.c.o" \
"/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/algebra/CMakeFiles/osqp_algebra.dir/default/matrix.c.o" \
"/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/algebra/CMakeFiles/osqp_algebra.dir/default/vector.c.o" \
"/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/algebra/CMakeFiles/osqp_algebra.dir/csc_tools/csc_math.c.o" \
"/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/algebra/CMakeFiles/osqp_algebra.dir/csc_tools/csc_utils.c.o" \
"/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/direct/kkt_common/CMakeFiles/kkt_common.dir/kkt.c.o" \
"/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/direct/qdldl/CMakeFiles/linsys_qdldl.dir/amd/src/amd_1.c.o" \
"/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/direct/qdldl/CMakeFiles/linsys_qdldl.dir/amd/src/amd_2.c.o" \
"/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/direct/qdldl/CMakeFiles/linsys_qdldl.dir/amd/src/amd_aat.c.o" \
"/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/direct/qdldl/CMakeFiles/linsys_qdldl.dir/amd/src/amd_control.c.o" \
"/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/direct/qdldl/CMakeFiles/linsys_qdldl.dir/amd/src/amd_defaults.c.o" \
"/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/direct/qdldl/CMakeFiles/linsys_qdldl.dir/amd/src/amd_info.c.o" \
"/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/direct/qdldl/CMakeFiles/linsys_qdldl.dir/amd/src/amd_order.c.o" \
"/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/direct/qdldl/CMakeFiles/linsys_qdldl.dir/amd/src/amd_post_tree.c.o" \
"/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/direct/qdldl/CMakeFiles/linsys_qdldl.dir/amd/src/amd_postorder.c.o" \
"/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/direct/qdldl/CMakeFiles/linsys_qdldl.dir/amd/src/amd_preprocess.c.o" \
"/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/direct/qdldl/CMakeFiles/linsys_qdldl.dir/amd/src/amd_valid.c.o" \
"/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/direct/qdldl/CMakeFiles/linsys_qdldl.dir/amd/src/SuiteSparse_config.c.o" \
"/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/direct/qdldl/CMakeFiles/linsys_qdldl.dir/qdldl_interface.c.o" \
"/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldlobject.dir/src/qdldl.c.o"

out/libosqp.so: CMakeFiles/osqp.dir/src/auxil.c.o
out/libosqp.so: CMakeFiles/osqp.dir/src/error.c.o
out/libosqp.so: CMakeFiles/osqp.dir/src/osqp_api.c.o
out/libosqp.so: CMakeFiles/osqp.dir/src/proj.c.o
out/libosqp.so: CMakeFiles/osqp.dir/src/scaling.c.o
out/libosqp.so: CMakeFiles/osqp.dir/src/util.c.o
out/libosqp.so: CMakeFiles/osqp.dir/src/polish.c.o
out/libosqp.so: CMakeFiles/osqp.dir/src/lin_sys.c.o
out/libosqp.so: CMakeFiles/osqp.dir/src/ctrlc.c.o
out/libosqp.so: CMakeFiles/osqp.dir/lin_sys/lib_handler.c.o
out/libosqp.so: algebra/CMakeFiles/osqp_algebra.dir/default/algebra_libs.c.o
out/libosqp.so: algebra/CMakeFiles/osqp_algebra.dir/default/matrix.c.o
out/libosqp.so: algebra/CMakeFiles/osqp_algebra.dir/default/vector.c.o
out/libosqp.so: algebra/CMakeFiles/osqp_algebra.dir/csc_tools/csc_math.c.o
out/libosqp.so: algebra/CMakeFiles/osqp_algebra.dir/csc_tools/csc_utils.c.o
out/libosqp.so: lin_sys/direct/kkt_common/CMakeFiles/kkt_common.dir/kkt.c.o
out/libosqp.so: lin_sys/direct/qdldl/CMakeFiles/linsys_qdldl.dir/amd/src/amd_1.c.o
out/libosqp.so: lin_sys/direct/qdldl/CMakeFiles/linsys_qdldl.dir/amd/src/amd_2.c.o
out/libosqp.so: lin_sys/direct/qdldl/CMakeFiles/linsys_qdldl.dir/amd/src/amd_aat.c.o
out/libosqp.so: lin_sys/direct/qdldl/CMakeFiles/linsys_qdldl.dir/amd/src/amd_control.c.o
out/libosqp.so: lin_sys/direct/qdldl/CMakeFiles/linsys_qdldl.dir/amd/src/amd_defaults.c.o
out/libosqp.so: lin_sys/direct/qdldl/CMakeFiles/linsys_qdldl.dir/amd/src/amd_info.c.o
out/libosqp.so: lin_sys/direct/qdldl/CMakeFiles/linsys_qdldl.dir/amd/src/amd_order.c.o
out/libosqp.so: lin_sys/direct/qdldl/CMakeFiles/linsys_qdldl.dir/amd/src/amd_post_tree.c.o
out/libosqp.so: lin_sys/direct/qdldl/CMakeFiles/linsys_qdldl.dir/amd/src/amd_postorder.c.o
out/libosqp.so: lin_sys/direct/qdldl/CMakeFiles/linsys_qdldl.dir/amd/src/amd_preprocess.c.o
out/libosqp.so: lin_sys/direct/qdldl/CMakeFiles/linsys_qdldl.dir/amd/src/amd_valid.c.o
out/libosqp.so: lin_sys/direct/qdldl/CMakeFiles/linsys_qdldl.dir/amd/src/SuiteSparse_config.c.o
out/libosqp.so: lin_sys/direct/qdldl/CMakeFiles/linsys_qdldl.dir/qdldl_interface.c.o
out/libosqp.so: lin_sys/direct/qdldl/qdldl_sources/CMakeFiles/qdldlobject.dir/src/qdldl.c.o
out/libosqp.so: CMakeFiles/osqp.dir/build.make
out/libosqp.so: CMakeFiles/osqp.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Linking C shared library out/libosqp.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/osqp.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/osqp.dir/build: out/libosqp.so

.PHONY : CMakeFiles/osqp.dir/build

CMakeFiles/osqp.dir/requires: CMakeFiles/osqp.dir/src/auxil.c.o.requires
CMakeFiles/osqp.dir/requires: CMakeFiles/osqp.dir/src/error.c.o.requires
CMakeFiles/osqp.dir/requires: CMakeFiles/osqp.dir/src/osqp_api.c.o.requires
CMakeFiles/osqp.dir/requires: CMakeFiles/osqp.dir/src/proj.c.o.requires
CMakeFiles/osqp.dir/requires: CMakeFiles/osqp.dir/src/scaling.c.o.requires
CMakeFiles/osqp.dir/requires: CMakeFiles/osqp.dir/src/util.c.o.requires
CMakeFiles/osqp.dir/requires: CMakeFiles/osqp.dir/src/polish.c.o.requires
CMakeFiles/osqp.dir/requires: CMakeFiles/osqp.dir/src/lin_sys.c.o.requires
CMakeFiles/osqp.dir/requires: CMakeFiles/osqp.dir/src/ctrlc.c.o.requires
CMakeFiles/osqp.dir/requires: CMakeFiles/osqp.dir/lin_sys/lib_handler.c.o.requires

.PHONY : CMakeFiles/osqp.dir/requires

CMakeFiles/osqp.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/osqp.dir/cmake_clean.cmake
.PHONY : CMakeFiles/osqp.dir/clean

CMakeFiles/osqp.dir/depend:
	cd /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP /home/armaan/Documents/cooperunion/spring_2020/ece453/final/cudaQP/CMakeFiles/osqp.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/osqp.dir/depend

