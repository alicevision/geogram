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
#include <vector>

#ifndef M_PI
/**
 * \brief Value of the constant PI if not defined by the system
 */
#define M_PI 3.14159265358979323846
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
    }

        /**
         * \brief Value of the null pointer
         * \internal Should be replaces by nullptr in C++11
         */
#define nil 0

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

}

#endif

