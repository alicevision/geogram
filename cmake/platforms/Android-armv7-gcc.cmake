#-------------------------------------------------------------------
# Flags for compiling with Android NDK
#-------------------------------------------------------------------

include(${GEOGRAM_SOURCE_DIR}/cmake/platforms/Linux.cmake)

# Set the Android compilers
set(CMAKE_CXX_COMPILER arm-linux-androideabi-g++)
set(CMAKE_C_COMPILER arm-linux-androideabi-gcc)

# No graphics (yet) for Android
set(GEOGRAM_WITH_GRAPHICS FALSE)

# Warning flags
set(NORMAL_WARNINGS -Wall -Wextra)
set(FULL_WARNINGS
    ${NORMAL_WARNINGS}
    -pedantic
    -Wno-long-long
    # Detect conversion problems (lot of warnings)
    -Wconversion
    -Wsign-conversion
    -Wdouble-promotion
)

# Compile with full warnings by default
add_definitions(${FULL_WARNINGS})

# Warn about missing virtual destructor (C++ only)
add_flags(CMAKE_CXX_FLAGS -Wnon-virtual-dtor)

# Add static and dynamic bounds checks (optimization required)
add_flags(CMAKE_CXX_FLAGS_RELEASE -D_FORTIFY_SOURCE=2)
add_flags(CMAKE_C_FLAGS_RELEASE -D_FORTIFY_SOURCE=2)

# Parameters for FPU instruction set and ABI:
# * -mfpu=neon-vfpv4: activate support for streaming vector floating point
#   (neon) and scalar floating point v4 (vfpv4). Note that vfpv4 supports
#   fma (fused multiply add), which is good, but which should not be 
#   generated automatically (-ffp-contract=off), since it will break
#   multi-precision routines.
# * -mfloat-abi=hard (or -mhard-float)
#  hard-float ABI (floating point values passed in vfp registers) needs 
#  a libstdc++ compiled with the same flags (NDK r9d has it, cool !)
# * define _NDK_MATH_NO_SOFTFP=1 and link with libm_hard
#  _NDK_MATH_NO_SOFTFP=1 ensures that math.h prototypes are not declared
#  as softfp ABI, and libm_hard is the version of libm compiled with hard ABI.

# This one works, but uses integer registers to return float/double results,
# a bit slower than what it could be... (however, it uses FP registers for
# the math functions, by declaring -D_NDK_MATH_NO_SOFTFP=1 and linking with m_hard).
set(ARCH_FLAGS -march=armv7-a -mfloat-abi=softfp -mfpu=neon-vfpv4 -ffp-contract=off -D_NDK_MATH_NO_SOFTFP=1)

# This one does not work everywhere (generates an executable that has an incoherent behavior,
#  maybe STL's C++ function do not have attribute((pcs("aapcs"))), this would explain....)
# It is not a problem anyway, performance seems to be the same.... (at least on the parts that work)
# set(ARCH_FLAGS -march=armv7-a -mhard-float -mfpu=neon-vfpv4 -ffp-contract=off -D_NDK_MATH_NO_SOFTFP=1)

add_flags(CMAKE_CXX_FLAGS ${ARCH_FLAGS})
add_flags(CMAKE_C_FLAGS ${ARCH_FLAGS})
add_flags(CMAKE_EXE_LINKER_FLAGS ${ARCH_FLAGS} -Wl,--fix-cortex-a8 -Wl,--no-warn-mismatch -lm_hard)

# Generate debug information even in release mode
#add_flags(CMAKE_CXX_FLAGS_RELEASE -g)
#add_flags(CMAKE_C_FLAGS_RELEASE -g)

# Additional debug flags
# deactivated for now: I added bound checking in VOR::vector<>.
#add_flags(CMAKE_CXX_FLAGS_DEBUG -D_GLIBCXX_DEBUG)


# Compile and link with OpenMP
add_flags(CMAKE_CXX_FLAGS -fopenmp -pthread)
add_flags(CMAKE_C_FLAGS -fopenmp -pthread)


# Profiler compilation flags
if(VORPALINE_WITH_GPROF)
    message(STATUS "Building for code profiling")
    add_flags(CMAKE_CXX_FLAGS -pg -DPROFILER)
    add_flags(CMAKE_C_FLAGS -pg -DPROFILER)
endif()


# Code coverage compilation flags
if(VORPALINE_WITH_GCOV)
    message(STATUS "Building for coverage analysis")
    add_flags(CMAKE_CXX_FLAGS --coverage)
    add_flags(CMAKE_C_FLAGS --coverage)
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
        # ... does not work with NDK 10.d
        #   (causes errors / multiply linked symbols)
#      add_flags(CMAKE_CXX_FLAGS -static-libstdc++ -static-libgcc -static)
#      add_flags(CMAKE_C_FLAGS -static-libgcc -static)
    endif()

    add_executable(${ARGN})
endmacro()

