/*
 *  Copyright (c) 2004-2014, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     levy@loria.fr
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 */

#ifndef __GEOGRAM_BASIC__
#define __GEOGRAM_BASIC__

#include <math.h>
#include <limits>
#include <float.h>
#include <limits.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#ifndef M_PI
/**
 * \brief Value of the constant PI if not defined by the system
 */
#define M_PI 3.14159265358979323846
#endif

/**
 * \def GEO_DEBUG
 * \brief This macro is set when compiling in debug mode
 *
 * \def GEO_PARANOID
 * \brief This macro is set when compiling in debug mode
 *
 * \def GEO_OS_LINUX
 * \brief This macro is set on Linux systems (Android included).
 *
 * \def GEO_OS_UNIX
 * \brief This macro is set on Unix systems (Android included).
 *
 * \def GEO_OS_WINDOWS
 * \brief This macro is set on Windows systems.
 *
 * \def GEO_OS_APPLE
 * \brief This macro is set on Apple systems.
 *
 * \def GEO_OS_ANDROID
 * \brief This macro is set on Android systems (in addition to GEO_OS_LINUX
 * and GEO_OS_UNIX).
 *
 * \def GEO_OS_X11
 * \brief This macro is set on X11 is supported on the current system.
 *
 * \def GEO_ARCH_32
 * \brief This macro is set if the current system is a 32 bits architecture.
 *
 * \def GEO_ARCH_64
 * \brief This macro is set if the current system is a 64 bits architecture.
 *
 * \def GEO_OPENMP
 * \brief This macro is set if OpenMP is supported on the current system.
 *
 * \def GEO_COMPILER_GCC
 * \brief This macro is set if the source code is compiled with GNU's gcc.
 *
 * \def GEO_COMPILER_INTEL
 * \brief This macro is set if the source code is compiled with Intel's icc.
 *
 * \def GEO_COMPILER_MSVC
 * \brief This macro is set if the source code is compiled with Microsoft's
 * Visual C++.
 */

#ifdef NDEBUG
#undef GEO_DEBUG
#undef GEO_PARANOID
#else
#define GEO_DEBUG
#define GEO_PARANOID
#endif

// =============================== LINUX defines ===========================

#if defined(__ANDROID__)
#define GEO_OS_ANDROID
#endif

#if defined(__linux__)

#define GEO_OS_LINUX
#define GEO_OS_UNIX

#ifndef GEO_OS_ANDROID
#define GEO_OS_X11
#endif

#if defined(_OPENMP)
#  define GEO_OPENMP
#endif

#if defined(__INTEL_COMPILER)
#  define GEO_COMPILER_INTEL
#elif defined(__clang__)
#  define GEO_COMPILER_CLANG
#elif defined(__GNUC__)
#  define GEO_COMPILER_GCC
#else
#  error "Unsupported compiler"
#endif

// The following works on GCC and ICC
#if defined(__x86_64)
#  define GEO_ARCH_64
#else
#  define GEO_ARCH_32
#endif

// =============================== WINDOWS defines =========================

#elif defined(WIN32)

#define GEO_OS_WINDOWS

#if defined(_OPENMP)
#  define GEO_OPENMP
#endif

#if defined(_MSC_VER)
#  define GEO_COMPILER_MSVC
#else
#  error "Unsupported compiler"
#endif

#if defined(_WIN64)
#  define GEO_ARCH_64
#else
#  define GEO_ARCH_32
#endif

// =============================== APPLE defines ===========================

#elif defined(__APPLE__)

#define GEO_OS_APPLE
#define GEO_OS_UNIX

#if defined(_OPENMP)
#  define GEO_OPENMP
#endif

#if defined(__x86_64) || defined(__ppc64__)
#  define GEO_ARCH_64
#else
#  define GEO_ARCH_32
#endif

// =============================== Unsupported =============================
#else

#error "Unsupported operating system"

#endif

#ifdef DOXYGEN_ONLY
// Keep doxygen happy
#define GEO_OS_WINDOWS
#define GEO_OS_APPLE
#define GEO_OS_ANDROID
#define GEO_ARCH_32
#define GEO_COMPILER_INTEL
#define GEO_COMPILER_MSVC
#endif


namespace GEO {

    /**
     * \brief The type for storing and manipulating indices.
     * \internal
     * Vorpaline uses 32 bit indices (can be changed to 64 bits
     * if need be, but this will double memory consumption of
     * all combinatorial data structures).
     */
    typedef unsigned int index_t;

    /**
     * \brief The type for storing and manipulating indices differences.
     * \details Can be negative (for instance to indicate special values like
     * borders).
     */
    typedef int signed_index_t;

    /**
     * \brief The type for storing coordinate indices, and iterating on
     *  the coordinates of a point.
     */
    typedef unsigned char coord_index_t;


    /**
     * \brief Integer constants that represent the sign of a value
     */
    enum Sign {
        /** Value is negative */
        NEGATIVE = -1,
        /** Value is zero */
        ZERO = 0,
        /** Value is positive */
        POSITIVE = 1
    };

    /**
     * \brief Gets the sign of a value
     * \details Returns -1, 0, or 1 whether value \p x is resp. negative, zero
     * or positive. The function uses operator<() and operator>() to compare
     * the value agains to 0 (zero). The integer constant zero must make
     * senses for the type of the value, or T must be constructible from
     * integer constant zero.
     * \param[in] x the value to test
     * \tparam T the type of the value
     * \return the sign of the value
     * \see Sign
     */
    template <class T>
    inline Sign geo_sgn(const T& x) {
        return (x > 0) ? POSITIVE : (
            (x < 0) ? NEGATIVE : ZERO
        );
    }

    /**
     * \brief Gets the absolute value of a value
     * \details The function uses operator< to compare the value againts 0
     * (zero). The integer constant zero must make senses for the type of the
     * value, or T must be constructible from integer constant zero.
     * \param[in] x a value of type \p T
     * \tparam T the type of the value
     * \return -x if the value is negative, x otherwise
     */
    template <class T>
    inline T geo_abs(T x) {
        return (x < 0) ? -x : x;
    }
 
    /**
     * \brief Gets the square value of a value
     * \param[in] x a value of type \p T
     * \tparam T the type of the value
     * \return the sqre value of \p x
     */
    template <class T>
    inline T geo_sqr(T x) {
        return x * x;
    }


    /**
     * \brief Computes a two-by-two determinant.
     */
    inline double det2x2(
        double a11, double a12,                    
        double a21, double a22
    ) {                                 
        return a11*a22-a12*a21 ;
    }

    /**
     * \brief Computes a three-by-three determinant.
     */
    inline double det3x3(
        double a11, double a12, double a13,                
        double a21, double a22, double a23,                
        double a31, double a32, double a33
    ) {
    return
         a11*det2x2(a22,a23,a32,a33)   
        -a21*det2x2(a12,a13,a32,a33)   
        +a31*det2x2(a12,a13,a22,a23);
    }   


    /**
     * \brief Computes a four-by-four determinant.
     */
    inline double det4x4(
        double a11, double a12, double a13, double a14,
        double a21, double a22, double a23, double a24,               
        double a31, double a32, double a33, double a34,  
        double a41, double a42, double a43, double a44  
    ) {
        double m12 = a21*a12 - a11*a22;
        double m13 = a31*a12 - a11*a32;
        double m14 = a41*a12 - a11*a42;
        double m23 = a31*a22 - a21*a32;
        double m24 = a41*a22 - a21*a42;
        double m34 = a41*a32 - a31*a42;

        double m123 = m23*a13 - m13*a23 + m12*a33;
        double m124 = m24*a13 - m14*a23 + m12*a43;
        double m134 = m34*a13 - m14*a33 + m13*a43;
        double m234 = m34*a23 - m24*a33 + m23*a43;
        
        return (m234*a14 - m134*a24 + m124*a34 - m123*a44);
    }   

    namespace Numeric {

        /** Generic pointer type */
        typedef void* pointer;

        /** Integer type with a width of 8 bits */
        typedef int8_t int8;

        /** Integer type with a width of 16 bits */
        typedef int16_t int16;

        /** Integer type with a width of 32 bits */
        typedef int32_t int32;

        /** Integer type with a width of 64 bits */
        typedef int64_t int64;

        /** Unsigned integer type with a width of 8 bits */
        typedef uint8_t uint8;

        /** Unsigned integer type with a width of 16 bits */
        typedef uint16_t uint16;

        /** Unsigned integer type with a width of 32 bits */
        typedef uint32_t uint32;

        /** Unsigned integer type with a width of 64 bits */
        typedef uint64_t uint64;

        /** Floating point type with a width of 32 bits */
        typedef float float32;

        /** Floating point type with a width of 64 bits */
        typedef double float64;

        /** Floating point type with a width of 64 bits */
        typedef double float64;

        /**
         * \brief Gets 64 bits float maximum positive value
         */
        inline float64 max_float64() {
            return std::numeric_limits<float64>::max();
        }

        /**
         * \brief Gets 64 bits float minimum negative value
         */
        inline float64 min_float64() {
            // Note: numeric_limits<>::min() is not
            // what we want (it returns the smallest
            // positive non-denormal).
            return -max_float64();
        }

        /**
         * \brief Returns a 32 bits integer between 0 and RAND_MAX
         */
        inline int32 random_int32() {
#ifdef GEO_OS_WINDOWS
            return rand();
#else
            return int32(random() % std::numeric_limits<int32>::max());
#endif
        }

    }
    

#ifdef GEO_DEBUG
#define geo_assert(x) assert(x)
#define geo_debug_assert(x) assert(x)
#define geo_assert_not_reached assert(0)
#else
#define geo_assert(x) assert(x)
#define geo_debug_assert(x) 
#define geo_assert_not_reached assert(0)
#endif    

    /**
     * \brief Utilities for memory management.
     */
    namespace Memory {
        /** \brief Unsigned byte type */
        typedef unsigned char byte;

        /** \brief Pointer to unsigned byte(s) */        
        typedef byte* pointer;

        /**
         * \brief Clears a memory block
         * \details Clears (set to zero) the first \p size bytes of array \p
         * addr.
         * \param[in] addr an array of bytes
         * \param[in] size the number of bytes to clear
         */
        inline void clear(void* addr, size_t size) {
            ::memset(addr, 0, size);
        }

        /**
         * \brief Copies a memory block
         * \details Copies the first \p size bytes of array \p from to array
         * \p to. Note that this function has unpredictable results if the
         * memory areas pointed to by \p to and \p from overlap.
         * \param[in] to the destination array of bytes
         * \param[in] from the array of bytes to copy
         * \param[in] size the number of bytes to copy
         */
        inline void copy(void* to, const void* from, size_t size) {
            ::memcpy(to, from, size);
        }

        /**
         * \brief Value of the null pointer
         * \internal Should be replaces by nullptr in C++11
         */
#define nil 0
    }


    /**
     * \brief Suppresses compiler warnings about unused parameters
     * \details This function is meant to get rid of warnings
     * concerning non used parameters (for instance,
     * parameters to callbacks). The corresponding code
     * is supposed to be wiped out by the optimizer.
     */
    template <class T> 
    inline void geo_argused(const T&) {
    }

/**
 * \brief Placeholder for linkage (dllimport/dllexport). All
 *  classes in Geogram are declared with GEOGRAM_API.
 */
#define GEOGRAM_API

    using std::vector;

    /**
     * \brief Some geometric functions.
     */
    namespace Geom {

        /**
         * \brief Computes the squared distance between two nd points.
         * \param[in] p1 a pointer to the coordinates of the first point
         * \param[in] p2 a pointer to the coordinates of the second point
         * \param[in] dim dimension (number of coordinates of the points)
         * \return the squared distance between \p p1 and \p p2
         * \tparam COORD_T the numeric type of the point coordinates
         */
        template <class COORD_T>
        inline double distance2(
            const COORD_T* p1, const COORD_T* p2, coord_index_t dim
        ) {
            double result = 0.0;
            for(coord_index_t i = 0; i < dim; i++) {
                result += geo_sqr(double(p2[i]) - double(p1[i]));
            }
            return result;
        }
    }

}

#endif

