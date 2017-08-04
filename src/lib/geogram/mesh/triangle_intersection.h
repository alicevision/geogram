/*
 *  Copyright (c) 2012-2014, Bruno Levy
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
 *     Bruno.Levy@inria.fr
 *     http://www.loria.fr/~levy
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 */

#ifndef GEOGRAM_MESH_TRIANGLE_INTERSECTION
#define GEOGRAM_MESH_TRIANGLE_INTERSECTION

#include <geogram/basic/common.h>
#include <geogram/basic/geometry.h>
#include <geogram/basic/memory.h>
#include <utility>

/**
 * \file geogram/mesh/triangle_intersection.h
 * \brief Symbolic computation of triangle-triangle intersection
 */

namespace GEO {

    /**
     * \brief Encodes the location of a point within a triangle.
     * \details A point can be located in 6 different regions, that
     * correspond to the three vertices, three edges and interior
     * of a triangle.
     */
    enum TriangleRegion {
        T_RGN_P0 = 0, /**< The point is on p0 */
        T_RGN_P1 = 1, /**< The point is on p1 */
        T_RGN_P2 = 2, /**< The point is on p2 */
        T_RGN_E0 = 3, /**< The point is on edge E0 */
        T_RGN_E1 = 4, /**< The point is on edge E1 */
        T_RGN_E2 = 5, /**< The point is on edge E2 */
        T_RGN_T = 6   /**< The point is in the interior of the triangle */
    };

    /**
     * \brief Encodes the symbolic representation of a triangle intersection,
     *  as a pair of TriangleRegion.
     */
    typedef std::pair<TriangleRegion, TriangleRegion> TriangleIsect;

    /**
     * \brief Triangle-triangle intersection
     * \param[in] p0 , p1 , p2 first triangle
     * \param[in] q0 , q1 , q2 second triangle
     * \param[out] result the intersection in symbolic
     *  form, as TriangleRegion pairs. There can be
     *  between 0 and 6 intersection pairs in the result.
     * \param[in] degenerate if true, then input triangles can be
     *  degenerate (i.e. with three vertices aligned, or three vertices
     *  equal). If not, then a degenerate triangle may trigger an assertion
     *  failure.
     * \retval true if there is a non-degenerate intersection
     * \retval false otherwise. Degenerate intersection cases are:
     *  - one vertex in common
     *  - two vertices (an edge) in common
     *  - or duplicated triangles.
     */
    bool GEOGRAM_API triangles_intersections(
        const vec3& p0, const vec3& p1, const vec3& p2,
        const vec3& q0, const vec3& q1, const vec3& q2,
        vector<TriangleIsect>& result
    );
}

#endif

