# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /mnt/c/DataCopied/Research/tree/source

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/c/DataCopied/Research/tree/build

# Include any dependencies generated for this target.
include CMakeFiles/merge_main.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/merge_main.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/merge_main.dir/flags.make

CMakeFiles/merge_main.dir/merge_main.cpp.o: CMakeFiles/merge_main.dir/flags.make
CMakeFiles/merge_main.dir/merge_main.cpp.o: /mnt/c/DataCopied/Research/tree/source/merge_main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/DataCopied/Research/tree/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/merge_main.dir/merge_main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/merge_main.dir/merge_main.cpp.o -c /mnt/c/DataCopied/Research/tree/source/merge_main.cpp

CMakeFiles/merge_main.dir/merge_main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/merge_main.dir/merge_main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/DataCopied/Research/tree/source/merge_main.cpp > CMakeFiles/merge_main.dir/merge_main.cpp.i

CMakeFiles/merge_main.dir/merge_main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/merge_main.dir/merge_main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/DataCopied/Research/tree/source/merge_main.cpp -o CMakeFiles/merge_main.dir/merge_main.cpp.s

# Object files for target merge_main
merge_main_OBJECTS = \
"CMakeFiles/merge_main.dir/merge_main.cpp.o"

# External object files for target merge_main
merge_main_EXTERNAL_OBJECTS =

libmerge_main.so: CMakeFiles/merge_main.dir/merge_main.cpp.o
libmerge_main.so: CMakeFiles/merge_main.dir/build.make
libmerge_main.so: CMakeFiles/merge_main.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/c/DataCopied/Research/tree/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library libmerge_main.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/merge_main.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/merge_main.dir/build: libmerge_main.so

.PHONY : CMakeFiles/merge_main.dir/build

CMakeFiles/merge_main.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/merge_main.dir/cmake_clean.cmake
.PHONY : CMakeFiles/merge_main.dir/clean

CMakeFiles/merge_main.dir/depend:
	cd /mnt/c/DataCopied/Research/tree/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/c/DataCopied/Research/tree/source /mnt/c/DataCopied/Research/tree/source /mnt/c/DataCopied/Research/tree/build /mnt/c/DataCopied/Research/tree/build /mnt/c/DataCopied/Research/tree/build/CMakeFiles/merge_main.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/merge_main.dir/depend

