# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

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
CMAKE_COMMAND = /home/aespielberg/Downloads/CLion-2018.2.2/clion-2018.2.2/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/aespielberg/Downloads/CLion-2018.2.2/clion-2018.2.2/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/aespielberg/ResearchCode/OpenFab

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/aespielberg/ResearchCode/OpenFab/cmake-build-debug

# Include any dependencies generated for this target.
include NetworkCommunication/CMakeFiles/NetworkCommunication.dir/depend.make

# Include the progress variables for this target.
include NetworkCommunication/CMakeFiles/NetworkCommunication.dir/progress.make

# Include the compile flags for this target's objects.
include NetworkCommunication/CMakeFiles/NetworkCommunication.dir/flags.make

NetworkCommunication/CMakeFiles/NetworkCommunication.dir/src/fossfab.pb.cc.o: NetworkCommunication/CMakeFiles/NetworkCommunication.dir/flags.make
NetworkCommunication/CMakeFiles/NetworkCommunication.dir/src/fossfab.pb.cc.o: ../NetworkCommunication/src/fossfab.pb.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/aespielberg/ResearchCode/OpenFab/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object NetworkCommunication/CMakeFiles/NetworkCommunication.dir/src/fossfab.pb.cc.o"
	cd /home/aespielberg/ResearchCode/OpenFab/cmake-build-debug/NetworkCommunication && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/NetworkCommunication.dir/src/fossfab.pb.cc.o -c /home/aespielberg/ResearchCode/OpenFab/NetworkCommunication/src/fossfab.pb.cc

NetworkCommunication/CMakeFiles/NetworkCommunication.dir/src/fossfab.pb.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NetworkCommunication.dir/src/fossfab.pb.cc.i"
	cd /home/aespielberg/ResearchCode/OpenFab/cmake-build-debug/NetworkCommunication && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/aespielberg/ResearchCode/OpenFab/NetworkCommunication/src/fossfab.pb.cc > CMakeFiles/NetworkCommunication.dir/src/fossfab.pb.cc.i

NetworkCommunication/CMakeFiles/NetworkCommunication.dir/src/fossfab.pb.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NetworkCommunication.dir/src/fossfab.pb.cc.s"
	cd /home/aespielberg/ResearchCode/OpenFab/cmake-build-debug/NetworkCommunication && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/aespielberg/ResearchCode/OpenFab/NetworkCommunication/src/fossfab.pb.cc -o CMakeFiles/NetworkCommunication.dir/src/fossfab.pb.cc.s

# Object files for target NetworkCommunication
NetworkCommunication_OBJECTS = \
"CMakeFiles/NetworkCommunication.dir/src/fossfab.pb.cc.o"

# External object files for target NetworkCommunication
NetworkCommunication_EXTERNAL_OBJECTS =

NetworkCommunication/libNetworkCommunicationd_.a: NetworkCommunication/CMakeFiles/NetworkCommunication.dir/src/fossfab.pb.cc.o
NetworkCommunication/libNetworkCommunicationd_.a: NetworkCommunication/CMakeFiles/NetworkCommunication.dir/build.make
NetworkCommunication/libNetworkCommunicationd_.a: NetworkCommunication/CMakeFiles/NetworkCommunication.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/aespielberg/ResearchCode/OpenFab/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libNetworkCommunicationd_.a"
	cd /home/aespielberg/ResearchCode/OpenFab/cmake-build-debug/NetworkCommunication && $(CMAKE_COMMAND) -P CMakeFiles/NetworkCommunication.dir/cmake_clean_target.cmake
	cd /home/aespielberg/ResearchCode/OpenFab/cmake-build-debug/NetworkCommunication && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/NetworkCommunication.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
NetworkCommunication/CMakeFiles/NetworkCommunication.dir/build: NetworkCommunication/libNetworkCommunicationd_.a

.PHONY : NetworkCommunication/CMakeFiles/NetworkCommunication.dir/build

NetworkCommunication/CMakeFiles/NetworkCommunication.dir/clean:
	cd /home/aespielberg/ResearchCode/OpenFab/cmake-build-debug/NetworkCommunication && $(CMAKE_COMMAND) -P CMakeFiles/NetworkCommunication.dir/cmake_clean.cmake
.PHONY : NetworkCommunication/CMakeFiles/NetworkCommunication.dir/clean

NetworkCommunication/CMakeFiles/NetworkCommunication.dir/depend:
	cd /home/aespielberg/ResearchCode/OpenFab/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/aespielberg/ResearchCode/OpenFab /home/aespielberg/ResearchCode/OpenFab/NetworkCommunication /home/aespielberg/ResearchCode/OpenFab/cmake-build-debug /home/aespielberg/ResearchCode/OpenFab/cmake-build-debug/NetworkCommunication /home/aespielberg/ResearchCode/OpenFab/cmake-build-debug/NetworkCommunication/CMakeFiles/NetworkCommunication.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : NetworkCommunication/CMakeFiles/NetworkCommunication.dir/depend

