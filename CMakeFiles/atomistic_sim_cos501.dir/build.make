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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/srijit/Dropbox/spaul/COS508_test/atomistic_sim_cos501

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/srijit/Dropbox/spaul/COS508_test/atomistic_sim_cos501

# Include any dependencies generated for this target.
include CMakeFiles/atomistic_sim_cos501.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/atomistic_sim_cos501.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/atomistic_sim_cos501.dir/flags.make

CMakeFiles/atomistic_sim_cos501.dir/main.cpp.o: CMakeFiles/atomistic_sim_cos501.dir/flags.make
CMakeFiles/atomistic_sim_cos501.dir/main.cpp.o: main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/srijit/Dropbox/spaul/COS508_test/atomistic_sim_cos501/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/atomistic_sim_cos501.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/atomistic_sim_cos501.dir/main.cpp.o -c /home/srijit/Dropbox/spaul/COS508_test/atomistic_sim_cos501/main.cpp

CMakeFiles/atomistic_sim_cos501.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/atomistic_sim_cos501.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/srijit/Dropbox/spaul/COS508_test/atomistic_sim_cos501/main.cpp > CMakeFiles/atomistic_sim_cos501.dir/main.cpp.i

CMakeFiles/atomistic_sim_cos501.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/atomistic_sim_cos501.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/srijit/Dropbox/spaul/COS508_test/atomistic_sim_cos501/main.cpp -o CMakeFiles/atomistic_sim_cos501.dir/main.cpp.s

CMakeFiles/atomistic_sim_cos501.dir/main.cpp.o.requires:

.PHONY : CMakeFiles/atomistic_sim_cos501.dir/main.cpp.o.requires

CMakeFiles/atomistic_sim_cos501.dir/main.cpp.o.provides: CMakeFiles/atomistic_sim_cos501.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/atomistic_sim_cos501.dir/build.make CMakeFiles/atomistic_sim_cos501.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/atomistic_sim_cos501.dir/main.cpp.o.provides

CMakeFiles/atomistic_sim_cos501.dir/main.cpp.o.provides.build: CMakeFiles/atomistic_sim_cos501.dir/main.cpp.o


# Object files for target atomistic_sim_cos501
atomistic_sim_cos501_OBJECTS = \
"CMakeFiles/atomistic_sim_cos501.dir/main.cpp.o"

# External object files for target atomistic_sim_cos501
atomistic_sim_cos501_EXTERNAL_OBJECTS =

atomistic_sim_cos501: CMakeFiles/atomistic_sim_cos501.dir/main.cpp.o
atomistic_sim_cos501: CMakeFiles/atomistic_sim_cos501.dir/build.make
atomistic_sim_cos501: CMakeFiles/atomistic_sim_cos501.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/srijit/Dropbox/spaul/COS508_test/atomistic_sim_cos501/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable atomistic_sim_cos501"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/atomistic_sim_cos501.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/atomistic_sim_cos501.dir/build: atomistic_sim_cos501

.PHONY : CMakeFiles/atomistic_sim_cos501.dir/build

CMakeFiles/atomistic_sim_cos501.dir/requires: CMakeFiles/atomistic_sim_cos501.dir/main.cpp.o.requires

.PHONY : CMakeFiles/atomistic_sim_cos501.dir/requires

CMakeFiles/atomistic_sim_cos501.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/atomistic_sim_cos501.dir/cmake_clean.cmake
.PHONY : CMakeFiles/atomistic_sim_cos501.dir/clean

CMakeFiles/atomistic_sim_cos501.dir/depend:
	cd /home/srijit/Dropbox/spaul/COS508_test/atomistic_sim_cos501 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/srijit/Dropbox/spaul/COS508_test/atomistic_sim_cos501 /home/srijit/Dropbox/spaul/COS508_test/atomistic_sim_cos501 /home/srijit/Dropbox/spaul/COS508_test/atomistic_sim_cos501 /home/srijit/Dropbox/spaul/COS508_test/atomistic_sim_cos501 /home/srijit/Dropbox/spaul/COS508_test/atomistic_sim_cos501/CMakeFiles/atomistic_sim_cos501.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/atomistic_sim_cos501.dir/depend

