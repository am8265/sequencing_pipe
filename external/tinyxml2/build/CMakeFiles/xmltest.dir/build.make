# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

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
CMAKE_COMMAND = /nfs/goldstein/software/cmake-3.5.2/usr/local/bin/cmake

# The command to remove a file.
RM = /nfs/goldstein/software/cmake-3.5.2/usr/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/dh2880/pipeline/sequencing_pipe/testing/tinyxml2/tinyxml2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/dh2880/pipeline/sequencing_pipe/testing/tinyxml2/build

# Include any dependencies generated for this target.
include CMakeFiles/xmltest.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/xmltest.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/xmltest.dir/flags.make

CMakeFiles/xmltest.dir/xmltest.cpp.o: CMakeFiles/xmltest.dir/flags.make
CMakeFiles/xmltest.dir/xmltest.cpp.o: /home/dh2880/pipeline/sequencing_pipe/testing/tinyxml2/tinyxml2/xmltest.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/dh2880/pipeline/sequencing_pipe/testing/tinyxml2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/xmltest.dir/xmltest.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/xmltest.dir/xmltest.cpp.o -c /home/dh2880/pipeline/sequencing_pipe/testing/tinyxml2/tinyxml2/xmltest.cpp

CMakeFiles/xmltest.dir/xmltest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/xmltest.dir/xmltest.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dh2880/pipeline/sequencing_pipe/testing/tinyxml2/tinyxml2/xmltest.cpp > CMakeFiles/xmltest.dir/xmltest.cpp.i

CMakeFiles/xmltest.dir/xmltest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/xmltest.dir/xmltest.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dh2880/pipeline/sequencing_pipe/testing/tinyxml2/tinyxml2/xmltest.cpp -o CMakeFiles/xmltest.dir/xmltest.cpp.s

CMakeFiles/xmltest.dir/xmltest.cpp.o.requires:

.PHONY : CMakeFiles/xmltest.dir/xmltest.cpp.o.requires

CMakeFiles/xmltest.dir/xmltest.cpp.o.provides: CMakeFiles/xmltest.dir/xmltest.cpp.o.requires
	$(MAKE) -f CMakeFiles/xmltest.dir/build.make CMakeFiles/xmltest.dir/xmltest.cpp.o.provides.build
.PHONY : CMakeFiles/xmltest.dir/xmltest.cpp.o.provides

CMakeFiles/xmltest.dir/xmltest.cpp.o.provides.build: CMakeFiles/xmltest.dir/xmltest.cpp.o


# Object files for target xmltest
xmltest_OBJECTS = \
"CMakeFiles/xmltest.dir/xmltest.cpp.o"

# External object files for target xmltest
xmltest_EXTERNAL_OBJECTS =

xmltest: CMakeFiles/xmltest.dir/xmltest.cpp.o
xmltest: CMakeFiles/xmltest.dir/build.make
xmltest: libtinyxml2.so.6.2.0
xmltest: CMakeFiles/xmltest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/dh2880/pipeline/sequencing_pipe/testing/tinyxml2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable xmltest"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/xmltest.dir/link.txt --verbose=$(VERBOSE)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Configuring xmltest resources directory: /home/dh2880/pipeline/sequencing_pipe/testing/tinyxml2/build/resources"
	/nfs/goldstein/software/cmake-3.5.2/usr/local/bin/cmake -E copy_directory /home/dh2880/pipeline/sequencing_pipe/testing/tinyxml2/tinyxml2/resources /home/dh2880/pipeline/sequencing_pipe/testing/tinyxml2/build/resources
	/nfs/goldstein/software/cmake-3.5.2/usr/local/bin/cmake -E make_directory /home/dh2880/pipeline/sequencing_pipe/testing/tinyxml2/build/resources/out

# Rule to build all files generated by this target.
CMakeFiles/xmltest.dir/build: xmltest

.PHONY : CMakeFiles/xmltest.dir/build

CMakeFiles/xmltest.dir/requires: CMakeFiles/xmltest.dir/xmltest.cpp.o.requires

.PHONY : CMakeFiles/xmltest.dir/requires

CMakeFiles/xmltest.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/xmltest.dir/cmake_clean.cmake
.PHONY : CMakeFiles/xmltest.dir/clean

CMakeFiles/xmltest.dir/depend:
	cd /home/dh2880/pipeline/sequencing_pipe/testing/tinyxml2/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/dh2880/pipeline/sequencing_pipe/testing/tinyxml2/tinyxml2 /home/dh2880/pipeline/sequencing_pipe/testing/tinyxml2/tinyxml2 /home/dh2880/pipeline/sequencing_pipe/testing/tinyxml2/build /home/dh2880/pipeline/sequencing_pipe/testing/tinyxml2/build /home/dh2880/pipeline/sequencing_pipe/testing/tinyxml2/build/CMakeFiles/xmltest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/xmltest.dir/depend

