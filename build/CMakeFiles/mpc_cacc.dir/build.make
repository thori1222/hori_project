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
CMAKE_SOURCE_DIR = /home/thori/src/hori_project

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/thori/src/hori_project/build

# Include any dependencies generated for this target.
include CMakeFiles/mpc_cacc.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/mpc_cacc.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mpc_cacc.dir/flags.make

CMakeFiles/mpc_cacc.dir/main.cpp.o: CMakeFiles/mpc_cacc.dir/flags.make
CMakeFiles/mpc_cacc.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/thori/src/hori_project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/mpc_cacc.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mpc_cacc.dir/main.cpp.o -c /home/thori/src/hori_project/main.cpp

CMakeFiles/mpc_cacc.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mpc_cacc.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/thori/src/hori_project/main.cpp > CMakeFiles/mpc_cacc.dir/main.cpp.i

CMakeFiles/mpc_cacc.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mpc_cacc.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/thori/src/hori_project/main.cpp -o CMakeFiles/mpc_cacc.dir/main.cpp.s

CMakeFiles/mpc_cacc.dir/main.cpp.o.requires:

.PHONY : CMakeFiles/mpc_cacc.dir/main.cpp.o.requires

CMakeFiles/mpc_cacc.dir/main.cpp.o.provides: CMakeFiles/mpc_cacc.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/mpc_cacc.dir/build.make CMakeFiles/mpc_cacc.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/mpc_cacc.dir/main.cpp.o.provides

CMakeFiles/mpc_cacc.dir/main.cpp.o.provides.build: CMakeFiles/mpc_cacc.dir/main.cpp.o


# Object files for target mpc_cacc
mpc_cacc_OBJECTS = \
"CMakeFiles/mpc_cacc.dir/main.cpp.o"

# External object files for target mpc_cacc
mpc_cacc_EXTERNAL_OBJECTS =

mpc_cacc: CMakeFiles/mpc_cacc.dir/main.cpp.o
mpc_cacc: CMakeFiles/mpc_cacc.dir/build.make
mpc_cacc: CMakeFiles/mpc_cacc.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/thori/src/hori_project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable mpc_cacc"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mpc_cacc.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mpc_cacc.dir/build: mpc_cacc

.PHONY : CMakeFiles/mpc_cacc.dir/build

CMakeFiles/mpc_cacc.dir/requires: CMakeFiles/mpc_cacc.dir/main.cpp.o.requires

.PHONY : CMakeFiles/mpc_cacc.dir/requires

CMakeFiles/mpc_cacc.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/mpc_cacc.dir/cmake_clean.cmake
.PHONY : CMakeFiles/mpc_cacc.dir/clean

CMakeFiles/mpc_cacc.dir/depend:
	cd /home/thori/src/hori_project/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/thori/src/hori_project /home/thori/src/hori_project /home/thori/src/hori_project/build /home/thori/src/hori_project/build /home/thori/src/hori_project/build/CMakeFiles/mpc_cacc.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/mpc_cacc.dir/depend

