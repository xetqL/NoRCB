# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

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
CMAKE_COMMAND = /var/lib/snapd/snap/clion/137/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /var/lib/snapd/snap/clion/137/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/xetql/projects/cpp/NoRCB

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/xetql/projects/cpp/NoRCB/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/norcb.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/norcb.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/norcb.dir/flags.make

CMakeFiles/norcb.dir/src/main.cpp.o: CMakeFiles/norcb.dir/flags.make
CMakeFiles/norcb.dir/src/main.cpp.o: ../src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xetql/projects/cpp/NoRCB/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/norcb.dir/src/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/norcb.dir/src/main.cpp.o -c /home/xetql/projects/cpp/NoRCB/src/main.cpp

CMakeFiles/norcb.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/norcb.dir/src/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/xetql/projects/cpp/NoRCB/src/main.cpp > CMakeFiles/norcb.dir/src/main.cpp.i

CMakeFiles/norcb.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/norcb.dir/src/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/xetql/projects/cpp/NoRCB/src/main.cpp -o CMakeFiles/norcb.dir/src/main.cpp.s

CMakeFiles/norcb.dir/src/norcb.cpp.o: CMakeFiles/norcb.dir/flags.make
CMakeFiles/norcb.dir/src/norcb.cpp.o: ../src/norcb.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xetql/projects/cpp/NoRCB/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/norcb.dir/src/norcb.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/norcb.dir/src/norcb.cpp.o -c /home/xetql/projects/cpp/NoRCB/src/norcb.cpp

CMakeFiles/norcb.dir/src/norcb.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/norcb.dir/src/norcb.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/xetql/projects/cpp/NoRCB/src/norcb.cpp > CMakeFiles/norcb.dir/src/norcb.cpp.i

CMakeFiles/norcb.dir/src/norcb.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/norcb.dir/src/norcb.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/xetql/projects/cpp/NoRCB/src/norcb.cpp -o CMakeFiles/norcb.dir/src/norcb.cpp.s

# Object files for target norcb
norcb_OBJECTS = \
"CMakeFiles/norcb.dir/src/main.cpp.o" \
"CMakeFiles/norcb.dir/src/norcb.cpp.o"

# External object files for target norcb
norcb_EXTERNAL_OBJECTS =

norcb: CMakeFiles/norcb.dir/src/main.cpp.o
norcb: CMakeFiles/norcb.dir/src/norcb.cpp.o
norcb: CMakeFiles/norcb.dir/build.make
norcb: /usr/local/lib/libCGAL.so.13.0.3
norcb: /usr/lib/libgmpxx.so
norcb: /usr/lib/libmpfr.so
norcb: /usr/lib/libgmp.so
norcb: CMakeFiles/norcb.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/xetql/projects/cpp/NoRCB/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable norcb"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/norcb.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/norcb.dir/build: norcb

.PHONY : CMakeFiles/norcb.dir/build

CMakeFiles/norcb.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/norcb.dir/cmake_clean.cmake
.PHONY : CMakeFiles/norcb.dir/clean

CMakeFiles/norcb.dir/depend:
	cd /home/xetql/projects/cpp/NoRCB/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/xetql/projects/cpp/NoRCB /home/xetql/projects/cpp/NoRCB /home/xetql/projects/cpp/NoRCB/cmake-build-debug /home/xetql/projects/cpp/NoRCB/cmake-build-debug /home/xetql/projects/cpp/NoRCB/cmake-build-debug/CMakeFiles/norcb.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/norcb.dir/depend

