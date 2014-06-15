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

#ifndef __GEOGRAM_NUMERICS_PREDICATES__
#define __GEOGRAM_NUMERICS_PREDICATES__

#include <geogram/basic.h>

/**
 * \file geogram/numerics/predicates.h
 * \brief Filtered exact predicates.
 */

namespace GEO {

    /**
     * \brief PCK (Predicate Construction Kit) implements a set of
     *  geometric predicates. PCK uses arithmetic filters (Meyer and Pion),
     *  expansion arithmetics (Shewchuk) and simulation of simplicity 
     *  (Edelsbrunner).
     */
    namespace PCK {

        /**
         * \brief Computes the side of a point (given directly)
         *  relative to a bisector.
         * \details Computes the side of \f$ q0 \f$ relative to
         * \f$ \Pi(p0,p1) \f$.
         * Symbolic perturbation is applied whenever equality holds.
         * \param[in] p0 first extremity of the bisector
         * \param[in] p1 second extremity of the bisector
         * \param[in] q0 point to be tested
         * \param[in] DIM number of coordinates of the point
         * \retval POSITIVE if d(p0,q0) < d(p1,q0)
         * \retval NEGATIVE if d(p0,q0) > d(p1,q1)
         * \retval perturb() if f(p0,q0) = d(p1,q1),
         *  where \c perturb() denotes a globally
         *  consistent perturbation, that returns either POSITIVE or NEGATIVE
         * \note Only some specific dimensions are implemented (3,4,6 and 7)
         */
        Sign GEOGRAM_API side1_SOS(
            const double* p0, const double* p1,
            const double* q0,
            coord_index_t DIM
        );

        /**
         * \brief Computes the side of a point (given as the intersection
         *  between a segment and a bisector) relative to another bisector.
         * \details Computes the side of \f$ q = \Pi(p0,p1) \cap [q0,q1] \f$
         * relative to \f$ \Pi(p0,p2) \f$.
         * Symbolic perturbation is applied whenever equality holds.
         * \param[in] p0 first extremity of the bisectors
         * \param[in] p1 second extremity of the first bisector
         *  (that defines the intersection q)
         * \param[in] p2 second extremity of the second bisector
         *  (against which orientation is tested)
         * \param[in] q0 first extremity of the segment
         *  (that defines the intersection q)
         * \param[in] q1 second extremity of the segment
         *  (that defines the intersection q)
         * \retval POSITIVE if d(p0,q) < d(p2,q)
         * \retval NEGATIVE if d(p0,q) > d(p2,q)
         * \retval perturb() if d(p0,q) = d(p2,q),
         *  where \c perturb() denotes a globally
         *  consistent perturbation, that returns either POSITIVE or NEGATIVE
         * \note Only some specific dimensions are implemented (3,4,6 and 7)
         */
        Sign GEOGRAM_API side2_SOS(
            const double* p0, const double* p1, const double* p2,
            const double* q0, const double* q1,
            coord_index_t DIM
        );

        /**
         * \brief Computes the side of a point (given as the intersection
         *  between a facet and two bisectors) relative to another bisector.
         * \details Computes the side of
         *  \f$ q = \Pi(p0,p1) \cap Pi(p0,p2) \cap \Delta[q0,q1,q2] \f$
         * relative to \f$ \Pi(p0,p3) \f$.
         * Symbolic perturbation is applied whenever equality holds.
         * \param[in] p0 first extremity of the bisectors
         * \param[in] p1 second extremity of the first bisector
         *  (that defines the intersection q)
         * \param[in] p2 second extremity of the second bisector
         *  (that defines the intersection q)
         * \param[in] p3 second extremity of the third bisector
         *  (against which orientation is tested)
         * \param[in] q0 first vertex of the triangle
         *  (that defines the intersection q)
         * \param[in] q1 second vertex of the triangle
         *  (that defines the intersection q)
         * \param[in] q2 third vertex of the triangle
         *  (that defines the intersection q)
         * \retval POSITIVE if d(p0,q) < d(p3,q)
         * \retval NEGATIVE if d(p0,q) > d(p3,q)
         * \retval perturb() if d(p0,q) = d(p3,q),
         *  where \c perturb() denotes a globally
         *  consistent perturbation, that returns either POSITIVE or NEGATIVE
         * \note Only some specific dimensions are implemented (3,4,6 and 7)
         */
        Sign GEOGRAM_API side3_SOS(
            const double* p0, const double* p1, const double* p2, const double* p3,
            const double* q0, const double* q1, const double* q2,
            coord_index_t DIM
        );

        /**
         * \brief Computes the side of a point (given as the intersection
         *   between a tetrahedron and three bisectors) relative to
         *  another bisector.
         * \details Computes the side of
         *  \f$ q = \Pi(p0,p1) \cap Pi(p0,p2) \cap Pi(p0,p3)
         * \cap \Delta[q0,q1,q2,q3] \f$ relative to \f$ \Pi(p0,p4) \f$.
         * Symbolic perturbation is applied whenever equality holds.
         * \param[in] p0 first extremity of the bisectors
         * \param[in] p1 second extremity of the first bisector
         *  (that defines the intersection q)
         * \param[in] p2 second extremity of the second bisector
         *  (that defines the intersection q)
         * \param[in] p3 second extremity of the third bisector
         *  (that defines the intersection q)
         * \param[in] p4 second extremity of the fourth bisector
         *  (against which orientation is tested)
         * \param[in] q0 first vertex of the tetrahedron
         *  (that defines the intersection q)
         * \param[in] q1 second vertex of the tetrahedron
         *  (that defines the intersection q)
         * \param[in] q2 third vertex of the tetrahedron
         *  (that defines the intersection q)
         * \param[in] q3 third vertex of the tetrahedron
         *  (that defines the intersection q)
         * \retval POSITIVE if d(p0,q) < d(p4,q)
         * \retval NEGATIVE if d(p0,q) > d(p4,q)
         * \retval perturb() if d(p0,q) = d(p4,q),
         *  where \c perturb() denotes a globally
         *  consistent perturbation, that returns either POSITIVE or NEGATIVE
         * \note Only some specific dimensions are implemented (3,4,6 and 7)
         */
        Sign GEOGRAM_API side4_SOS(
            const double* p0,
            const double* p1, const double* p2, const double* p3, const double* p4,
            const double* q0, const double* q1, const double* q2, const double* q3,
            coord_index_t DIM
        );


        /**
         * \brief Computes the side of a point (given as the intersection
         *   between three bisectors) relative to another bisector.
         * \details Computes the side of
         *  \f$ q = \Pi(p0,p1) \cap \Pi(p0,p2) \cap \Pi(p0,p3) \f$
         * relative to \f$ Pi(p0,p4) \f$.
         * This version does not apply symbolic perturbation when equality
         * holds.
         * side4_3d() is a special case of side4(), where the ambient and
         * intrinsic dimensions coincide (therefore no embedding tetrahedron
         * is needed).
         * \param[in] p0 first extremity of the bisectors
         * \param[in] p1 second extremity of the first bisector
         *  (that defines the intersection q)
         * \param[in] p2 second extremity of the second bisector
         *  (that defines the intersection q)
         * \param[in] p3 second extremity of the third bisector
         *  (that defines the intersection q)
         * \param[in] p4 second extremity of the fourth bisector
         *  (against which orientation is tested)
         * \retval POSITIVE if d(p0,q) < d(p4,q)
         * \retval NEGATIVE if d(p0,q) > d(p4,q)
         * \retval ZERO if d(p0,q) = d(p4,q),
         */
        Sign GEOGRAM_API side4_3d(
            const double* p0,
            const double* p1, const double* p2, const double* p3, const double* p4
        );

        /**
         * \brief Computes the side of a point (given as the intersection
         *   between three bisectors) relative to another bisector.
         * \details Computes the side of
         *  \f$ q = \Pi(p0,p1) \cap \Pi(p0,p2) \cap \Pi(p0,p3) \f$
         * relative to \f$ Pi(p0,p4) \f$.
         * Symbolic perturbation is applied whenever equality holds.
         * side4_3d() is a special case of side4(), where the ambient and
         * intrinsic dimensions coincide (therefore no embedding tetrahedron
         * is needed).
         * \param[in] p0 first extremity of the bisectors
         * \param[in] p1 second extremity of the first bisector
         *  (that defines the intersection q)
         * \param[in] p2 second extremity of the second bisector
         *  (that defines the intersection q)
         * \param[in] p3 second extremity of the third bisector
         *  (that defines the intersection q)
         * \param[in] p4 second extremity of the fourth bisector
         *  (against which orientation is tested)
         * \retval POSITIVE if d(p0,q) < d(p4,q)
         * \retval NEGATIVE if d(p0,q) > d(p4,q)
         * \retval perturb() if d(p0,q) = d(p4,q),
         *  where \c perturb() denotes a globally
         *  consistent perturbation, that returns either POSITIVE or NEGATIVE
         */
        Sign GEOGRAM_API side4_3d_SOS(
            const double* p0, const double* p1, 
            const double* p2, const double* p3, const double* p4
        );
       
        /**
         * \brief Tests whether a 3d point is inside the circumscribed sphere of a
         *  3d tetrahedron.
         * \param[in] p0 first vertex of the tetrahedron
         * \param[in] p1 second vertex of the tetrahedron
         * \param[in] p2 third vertex of the tetrahedron
         * \param[in] p3 fourth vertex of the tetrahedron
         * \param[in] p4 the point to be tested
         * \retval POSITIVE whenever \p p4 is inside the circumscribed sphere
         *  of the tetrahedron \p p0, \p p1, \p p2, \p p3
         * \retval NEGATIVE whenever \p p4 is outside the circumscribed sphere
         *  of the tetrahedron \p p0, \p p1, \p p2, \p p3
         * \retval perturb() if \p p4 is exactly on the circumscribed sphere
         *  of the tetrahedron \p p0, \p p1, \p p2, \p p3, where \c perturb()
         *  denotes a globally consistent perturbation, that returns
         *  either POSITIVE or NEGATIVE
         * \pre orient_3d(p0,p1,p2,p3) > 0
         */
         Sign GEOGRAM_API in_sphere_3d_SOS(
            const double* p0, const double* p1, const double* p2, const double* p3,
            const double* p4
         );

        /**
         * \brief Computes the orientation predicate in 3d.
         * \details Computes the sign of the signed area of
         *  the ttriangle p0, p1, p2.
         * \param[in] p0 first vertex of the triangle
         * \param[in] p1 second vertex of the triangle
         * \param[in] p2 third vertex of the triangle
         * \retval POSITIVE if the triangle is oriented positively
         * \retval ZERO if the triangle is flat
         * \retval NEGATIVE if the triangle is oriented negatively
         * \todo check whether orientation is inverted as compared to 
         *   Shewchuk's version.
         */
        Sign GEOGRAM_API orient_2d(
            const double* p0, const double* p1, const double* p2
        );

        /**
         * \brief Computes the orientation predicate in 3d.
         * \details Computes the sign of the signed volume of
         *  the tetrahedron p0, p1, p2, p3.
         * \param[in] p0 first vertex of the tetrahedron
         * \param[in] p1 second vertex of the tetrahedron
         * \param[in] p2 third vertex of the tetrahedron
         * \param[in] p3 fourth vertex of the tetrahedron
         * \retval POSITIVE if the tetrahedron is oriented positively
         * \retval ZERO if the tetrahedron is flat
         * \retval NEGATIVE if the tetrahedron is oriented negatively
         * \todo check whether orientation is inverted as compared to 
         *   Shewchuk's version.
         */
        Sign GEOGRAM_API orient_3d(
            const double* p0, const double* p1,
            const double* p2, const double* p3
        );

        /**
         * \brief Computes the 4d orientation test.
         * \details Given four lifted points p0', p1', p2', and p3' in 
         * R^4, tests if the lifted point p4' in R^4 lies below or above 
         * the hyperplance passing through the four points p0', p1', p2', and p3'.
         *  This version does not apply symbolic perturbation.
         *  The first three coordinates and the
         *  fourth one are specified in separate arguments for each vertex.
         * \param[in] p0 first 3 coordinates of the first vertex of the 4-simplex
         * \param[in] p1 first 3 coordinates of the second vertex of the 4-simplex
         * \param[in] p2 first 3 coordinates of the third vertex of the 4-simplex
         * \param[in] p3 first 3 coordinates of the fourth vertex of the 4-simplex
         * \param[in] p4 first 3 coordinates of the fifth vertex of the 4-simplex
         * \param[in] h0 height of the first vertex of the 4-simplex
         * \param[in] h1 height of the second vertex of the 4-simplex
         * \param[in] h2 height of the third vertex of the 4-simplex
         * \param[in] h3 height of the fourth vertex of the 4-simplex
         * \param[in] h4 height of the fifth vertex of the 4-simplex
         * \retval POSITIVE if p4' lies below the hyperplane
         * \retval NEGATIVE if p4' lies above the hyperplane
         * \retval ZERO if p4' lies exactly on the hyperplane
         */
        Sign GEOGRAM_API orient_4d(
            const double* p0, const double* p1,
            const double* p2, const double* p3, const double* p4,
            double h0, double h1, double h2, double h3, double h4
        );


        /**
         * \brief Computes the 4d orientation test with symbolic perturbation.
         * \details Given four lifted points p0', p1', p2', and p3' in 
         * R^4, tests if the lifted point p4' in R^4 lies below or above 
         * the hyperplance passing through the four points p0', p1', p2', and p3'.
         *  Symbolic perturbation is applied whenever the 5 vertices are
         *  not linearly independent. The first three coordinates and the
         *  fourth one are specified in separate arguments for each vertex.
         * \param[in] p0 first 3 coordinates of the first vertex of the 4-simplex
         * \param[in] p1 first 3 coordinates of the second vertex of the 4-simplex
         * \param[in] p2 first 3 coordinates of the third vertex of the 4-simplex
         * \param[in] p3 first 3 coordinates of the fourth vertex of the 4-simplex
         * \param[in] p4 first 3 coordinates of the fifth vertex of the 4-simplex
         * \param[in] h0 height of the first vertex of the 4-simplex
         * \param[in] h1 height of the second vertex of the 4-simplex
         * \param[in] h2 height of the third vertex of the 4-simplex
         * \param[in] h3 height of the fourth vertex of the 4-simplex
         * \param[in] h4 height of the fifth vertex of the 4-simplex
         * \retval POSITIVE if p4' lies below the hyperplane
         * \retval NEGATIVE if p4' lies above the hyperplane
         * \retval perturb() if p4' lies exactly on the hyperplane
         *  where \c perturb() denotes a globally
         *  consistent perturbation, that returns either POSITIVE or NEGATIVE
         */
        Sign GEOGRAM_API orient_4d_SOS(
            const double* p0, const double* p1,
            const double* p2, const double* p3, const double* p4,
            double h0, double h1, double h2, double h3, double h4
        );

        /**
         * \brief Displays some statistics about predicates,
         *  including the number of calls, the number of exact arithmetics 
         *  calls, and the number of Simulation of Simplicity calls.
         */
        void GEOGRAM_API show_stats();
    }
}

#endif

