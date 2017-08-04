#-------------------------------------------------------------------
# Flags common to all Linux based platforms with Intel compiler
#-------------------------------------------------------------------

include(${GEOGRAM_SOURCE_DIR}/cmake/platforms/Linux.cmake)

# Set the Intel compilers
set(CMAKE_CXX_COMPILER icpc)
set(CMAKE_C_COMPILER icc)
set(CMAKE_Fortran_COMPILER ifort)

# Warning flags
add_definitions(
    -Wall
    -Wno-deprecated
    -Wunused
    -diag-disable 1881
    # Disable annoying remarks for Intel C++ compiler
    # remark #304: access control not specified ("public" by default)
    -wd,383,981,304
)

# Set the C standard
add_flags(CMAKE_C_FLAGS -std=c99)

# Add Intel system includes
add_definitions(
    -isystem $ENV{INTEL}/include
    -isystem $ENV{INTEL}/include/intel64
    -isystem $ENV{INTEL}/mkl/include
    -isystem $ENV{INTEL}/ipp/include
)

# Compile and link with OpenMP
add_flags(CMAKE_CXX_FLAGS -openmp -restrict)
add_flags(CMAKE_C_FLAGS -openmp -restrict)

# Link flags to force link of shared libs to resolve all the symbols
add_flags(CMAKE_EXE_LINKER_FLAGS -z defs)

# Flags for the Release build type
#-axSSE3,SSSE3,SSE4.2,AVX -vec-report6
#add_flags(CMAKE_CXX_FLAGS_RELEASE -xHost -ip -axSSE3)
#add_flags(CMAKE_C_FLAGS_RELEASE -xHost -ip -axSSE3)


# icc options related with FPU:
# https://software.intel.com/en-us/node/522979
# default fp model is "fast" (not good !!)
# 'strict' implies both 'precise' and 'except'

add_flags(CMAKE_CXX_FLAGS -fp-model strict)
add_flags(CMAKE_C_FLAGS -fp-model strict)

# Reset the warning level for third parties
function(vor_reset_warning_level)
endfunction()

macro(vor_add_executable)
    if(NOT VORPALINE_BUILD_DYNAMIC)
        # Create a statically linked executable
        # Link with static libraries: temporarily deactivated (causes some errors in link phase)
        #add_flags(CMAKE_CXX_FLAGS -static -static-intel -static-libgcc)
        #add_flags(CMAKE_C_FLAGS -static -static-intel -static-libgcc)
    endif()

    add_executable(${ARGN})
endmacro()

