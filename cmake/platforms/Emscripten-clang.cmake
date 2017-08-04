#-------------------------------------------------------------------
# Flags for Emscripten (javascript target)
#-------------------------------------------------------------------

# Shell script extension
set(SHELL_SUFFIX "sh")

find_path(EMSCRIPTEN_DIR
      emcc
      HINTS
        ENV EMSCRIPTEN
      PATHS
        "C:/Program Files/emscripten"
         /usr/lib/emscripten
)

set(CMAKE_C_COMPILER "")

set(CMAKE_CXX_COMPILER "")

include(${EMSCRIPTEN_DIR}/cmake/Modules/Platform/Emscripten.cmake)

set(GEOGRAM_WITH_EMSCRIPTEN TRUE)

# Warning flags
set(NORMAL_WARNINGS -Wall -Wextra)

set(FULL_WARNINGS
  -Weverything
    -Wno-disabled-macro-expansion # else we got a warning each time cout is used
    -Wno-padded # Disable generating a message each time padding is used
    -Wno-float-equal # Sometimes we compare floats (against 0.0 or 1.0 mainly)
    -Wno-global-constructors
    -Wno-exit-time-destructors
    -Wno-old-style-cast # Yes, old-style cast is sometime more legible...
    -Wno-format-nonliteral # Todo: use Laurent Alonso's trick
)

# Compile with full warnings by default
add_definitions(${FULL_WARNINGS})

# Run the static analyzer
if(VORPALINE_WITH_CLANGSA)
    add_definitions(--analyze)
endif()

#opengl flags
#-s FORCE_ALIGNED_MEMORY=1 
#https://kripken.github.io/emscripten-site/docs/optimizing/Optimizing-Code.html
#-s ALLOW_MEMORY_GROWTH=1
#-s TOTAL_MEMORY=256000000
set(EM_FLAGS -O2 -s USE_GLFW=3 -s TOTAL_MEMORY=256000000 )

# Add static and dynamic bounds checks (optimization required)
add_flags_no_remove_duplicates(CMAKE_CXX_FLAGS_RELEASE ${EM_FLAGS})
add_flags_no_remove_duplicates(CMAKE_C_FLAGS_RELEASE ${EM_FLAGS})

# Enable glibc parallel mode
#add_flags(CMAKE_CXX_FLAGS -D_GLIBCXX_PARALLEL)

# Generate debug information even in release mode
#add_flags(CMAKE_CXX_FLAGS_RELEASE -g)
#add_flags(CMAKE_C_FLAGS_RELEASE -g)

# Additional debug flags
# deactivated for now: I added bound checking in VOR::vector<>.
#add_flags(CMAKE_CXX_FLAGS_DEBUG -D_GLIBCXX_DEBUG)

# Compile and link with OpenMP ** NOT YET SUPPORTED in clang 3 **
#add_flags(CMAKE_CXX_FLAGS -fopenmp)
#add_flags(CMAKE_C_FLAGS -fopenmp)

# Profiler compilation flags
if(VORPALINE_WITH_GPROF)
    message(FATAL_ERROR "Profiling is not (yet) available with clang")
    message(STATUS "Building for code profiling")
    #add_flags(CMAKE_CXX_FLAGS -pg -DPROFILER)
    #add_flags(CMAKE_C_FLAGS -pg -DPROFILER)
endif()

# Code coverage compilation flags
if(VORPALINE_WITH_GCOV)
    message(STATUS "Building for coverage analysis")
    add_flags(CMAKE_CXX_FLAGS --coverage)
    add_flags(CMAKE_C_FLAGS --coverage)
endif()

# Compilation flags for Google's AddressSanitizer
# These flags can only be specified for dynamic builds
if(VORPALINE_WITH_ASAN)
    if(VORPALINE_BUILD_DYNAMIC)
        message(STATUS "Building with AddressSanitizer (debug only)")
        add_flags(CMAKE_CXX_FLAGS_DEBUG -fsanitize=address -fno-omit-frame-pointer)
        add_flags(CMAKE_C_FLAGS_DEBUG -fsanitize=address -fno-omit-frame-pointer)
    else()
        message(WARNING "AddressSanitizer can be used with dynamic builds only")
        set(VORPALINE_WITH_ASAN false)
    endif()
endif()
if(NOT VORPALINE_WITH_ASAN)
    # Use native GCC stack smash Protection and buffer overflow detection (debug only)
    add_flags(CMAKE_CXX_FLAGS_DEBUG -fstack-protector-all)
    add_flags(CMAKE_C_FLAGS_DEBUG -fstack-protector-all)
endif()


# Compilation flags for Google's ThreadSanitizer
# Does not work for the moment: cannot figure out how to link with library libtsan
if(VORPALINE_WITH_TSAN)
    message(STATUS "Building with ThreadSanitizer (debug only)")
    message(FATAL_ERROR "ThreadSanitizer is not available: cannot figure out how to link with library libtsan")
    add_flags(CMAKE_CXX_FLAGS_DEBUG -fsanitize=thread)
    add_flags(CMAKE_C_FLAGS_DEBUG -fsanitize=thread)
    if(NOT VORPALINE_BUILD_DYNAMIC)
        add_flags(CMAKE_EXE_LINKER_FLAGS -static-libtsan)
    endif()
endif()


# Reset the warning level for third parties
function(vor_reset_warning_level)
    remove_definitions(${FULL_WARNINGS})
    add_definitions(${NORMAL_WARNINGS})
endfunction()

macro(vor_add_executable)
    if(NOT VORPALINE_BUILD_DYNAMIC)
        # Create a statically linked executable
        # Link with static libraries
        add_flags(CMAKE_CXX_FLAGS -static)
        add_flags(CMAKE_C_FLAGS -static)
    endif()

    add_executable(${ARGN})

    if(NOT VORPALINE_BUILD_DYNAMIC AND DEFINED VORPALINE_WITH_DDT)
        # Static builds running with Allinea's DDT must be linked with a
        # special malloc library which replaces the malloc primitives of
        # the Glibc (We must allow multiple definitions)
        add_flags(CMAKE_EXE_LINKER_FLAGS -Wl,--allow-multiple-definition)

        if(VORPALINE_ARCH_64)
            link_directories(${VORPALINE_WITH_DDT}/lib/64)
        else()
            link_directories(${VORPALINE_WITH_DDT}/lib/32)
        endif()
        target_link_libraries(${ARGV0} dmallocthcxx)
    endif()

    if(UNIX)
        target_link_libraries(${ARGV0} m)
    endif()

endmacro()

