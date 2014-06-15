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

#include <geogram/numerics/predicates.h>
#include <geogram/numerics/multi_precision.h>
#include <algorithm>

#define FPG_UNCERTAIN_VALUE 0

#include <geogram/numerics/predicates/side1.h>
#include <geogram/numerics/predicates/side2.h>
#include <geogram/numerics/predicates/side3.h>
#include <geogram/numerics/predicates/side4.h>
#include <geogram/numerics/predicates/side4h.h>
#include <geogram/numerics/predicates/orient2d.h>
#include <geogram/numerics/predicates/orient3d.h>

namespace {
    using namespace GEO;




    /**
     * \brief Arithmetic filter for the in_sphere_3d_SOS() predicate.
     * \details This filter was optimized by hand by Sylvain Pion
     *  (may be faster than FPG/PCK-generated filter).
     * Since it is used massively by Delaunay_3d, using the
     *   optimized version may be worth it.
     * \param[in] p first vertex of the tetrahedron
     * \param[in] q second vertex of the tetrahedron
     * \param[in] r third vertex of the tetrahedron
     * \param[in] s fourth vertex of the tetrahedron
     * \param[in] t point to be tested
     * \retval +1 if \p t was determined to be outside 
     *   the circumsphere of \p p,\p q,\p r,\p s
     * \retval -1 if \p t was determined to be inside 
     *   the circumsphere of  \p p,\p q,\p r,\p s
     * \retval 0 if the position of \p t could be be determined
     */
    inline int in_sphere_3d_filter_optim(
        const double* p, const double* q, 
        const double* r, const double* s, const double* t
    ) {

        double ptx = p[0] - t[0];
        double pty = p[1] - t[1];
        double ptz = p[2] - t[2];
        double pt2 = geo_sqr(ptx) + geo_sqr(pty) + geo_sqr(ptz);

        double qtx = q[0] - t[0];
        double qty = q[1] - t[1];
        double qtz = q[2] - t[2];
        double qt2 = geo_sqr(qtx) + geo_sqr(qty) + geo_sqr(qtz);

        double rtx = r[0] - t[0];
        double rty = r[1] - t[1];
        double rtz = r[2] - t[2];
        double rt2 = geo_sqr(rtx) + geo_sqr(rty) + geo_sqr(rtz);

        double stx = s[0] - t[0];
        double sty = s[1] - t[1];
        double stz = s[2] - t[2];
        double st2 = geo_sqr(stx) + geo_sqr(sty) + geo_sqr(stz);

        // Compute the semi-static bound.
        double maxx = ::fabs(ptx);
        double maxy = ::fabs(pty);
        double maxz = ::fabs(ptz);
        
        double aqtx = ::fabs(qtx);
        double artx = ::fabs(rtx);
        double astx = ::fabs(stx);
        
        double aqty = ::fabs(qty);
        double arty = ::fabs(rty);
        double asty = ::fabs(sty);
        
        double aqtz = ::fabs(qtz);
        double artz = ::fabs(rtz);
        double astz = ::fabs(stz);

        if (maxx < aqtx) maxx = aqtx;
        if (maxx < artx) maxx = artx;
        if (maxx < astx) maxx = astx;
        
        if (maxy < aqty) maxy = aqty;
        if (maxy < arty) maxy = arty;
        if (maxy < asty) maxy = asty;
        
        if (maxz < aqtz) maxz = aqtz;
        if (maxz < artz) maxz = artz;
        if (maxz < astz) maxz = astz;

        double eps = 1.2466136531027298e-13 * maxx * maxy * maxz;
  
        // Sort maxx < maxy < maxz.
        if (maxx > maxz)
            std::swap(maxx, maxz);
        if (maxy > maxz)
            std::swap(maxy, maxz);
        else if (maxy < maxx)
            std::swap(maxx, maxy);

        double det = det4x4(
                        ptx,pty,ptz,pt2,
                        rtx,rty,rtz,rt2,
                        qtx,qty,qtz,qt2,
                        stx,sty,stz,st2
                     );

        if (maxx < 1e-58)  { /* sqrt^5(min_double/eps) */
            // Protect against underflow in the computation of eps.
            return FPG_UNCERTAIN_VALUE;
        } else if (maxz < 1e61)  { /* sqrt^5(max_double/4 [hadamard]) */
            // Protect against overflow in the computation of det.
            eps *= (maxz * maxz);
            // Note: inverted as compared to CGAL
            if (det > eps)  return -1;
            if (det < -eps) return  1;
        }

        return FPG_UNCERTAIN_VALUE;
    }


    using namespace GEO;

    index_t cnt_side1_total = 0;
    index_t cnt_side1_exact = 0;
    index_t cnt_side1_SOS = 0;
    index_t len_side1 = 0;

    index_t cnt_side2_total = 0;
    index_t cnt_side2_exact = 0;
    index_t cnt_side2_SOS = 0;
    index_t len_side2_num = 0;
    index_t len_side2_denom = 0;
    index_t len_side2_SOS = 0;

    index_t cnt_side3_total = 0;
    index_t cnt_side3_exact = 0;
    index_t cnt_side3_SOS = 0;
    index_t len_side3_num = 0;
    index_t len_side3_denom = 0;
    index_t len_side3_SOS = 0;

    index_t cnt_side4_total = 0;
    index_t cnt_side4_exact = 0;
    index_t cnt_side4_SOS = 0;
    index_t len_side4_num = 0;
    index_t len_side4_denom = 0;
    index_t len_side4_SOS = 0;

    index_t cnt_orient2d_total = 0;
    index_t cnt_orient2d_exact = 0;
    index_t len_orient2d = 0;

    index_t cnt_orient3d_total = 0;
    index_t cnt_orient3d_exact = 0;
    index_t len_orient3d = 0;

    index_t cnt_orient4d_total = 0;
    index_t cnt_orient4d_exact = 0;
    index_t cnt_orient4d_SOS = 0;
    index_t len_orient4d_num = 0;
    index_t len_orient4d_denom = 0;
    index_t len_orient4d_SOS = 0;

    // ================= side1 =========================================

    /**
     * \brief Exact implementation of the side1() predicate using low-level
     *  exact arithmetics API (expansion class).
     */
    Sign side1_exact_SOS(
        const double* p0, const double* p1,
        const double* q0,
        coord_index_t dim
    ) {
        cnt_side1_exact++;
        expansion& l = expansion_sq_dist(p0, p1, dim);
        expansion& a = expansion_dot_at(p1, q0, p0, dim).scale_fast(2.0);
        expansion& r = expansion_diff(l, a);
        Sign r_sign = r.sign();
        // Symbolic perturbation, Simulation of Simplicity
        if(r_sign == ZERO) {
            cnt_side1_SOS++;
            return (p0 < p1) ? POSITIVE : NEGATIVE;
        }
        len_side1 = std::max(len_side1, r.length());
        return r_sign;
    }

    /**
     * \brief Implements side1() in 3d.
     */
    Sign side1_3d_SOS(
        const double* p0, const double* p1, const double* q0
    ) {
        Sign result = Sign(side1_3d_filter(p0, p1, q0));
        if(result == ZERO) {
            result = side1_exact_SOS(p0, p1, q0, 3);
        }
        return result;
    }

    /**
     * \brief Implements side1() in 4d.
     */
    Sign side1_4d_SOS(
        const double* p0, const double* p1, const double* q0
    ) {
        Sign result = Sign(side1_4d_filter(p0, p1, q0));
        if(result == ZERO) {
            result = side1_exact_SOS(p0, p1, q0, 4);
        }
        return result;
    }

    /**
     * \brief Implements side1() in 6d.
     */
    Sign side1_6d_SOS(
        const double* p0, const double* p1, const double* q0
    ) {
        Sign result = Sign(side1_6d_filter(p0, p1, q0));
        if(result == ZERO) {
            result = side1_exact_SOS(p0, p1, q0, 6);
        }
        return result;
    }

    /**
     * \brief Implements side1() in 7d.
     */
    Sign side1_7d_SOS(
        const double* p0, const double* p1, const double* q0
    ) {
        Sign result = Sign(side1_7d_filter(p0, p1, q0));
        if(result == ZERO) {
            result = side1_exact_SOS(p0, p1, q0, 7);
        }
        return result;
    }

    // ================= side2 =========================================

    /**
     * \brief Exact implementation of the side2() predicate using low-level
     *  exact arithmetics API (expansion class).
     */
    Sign side2_exact_SOS(
        const double* p0, const double* p1, const double* p2,
        const double* q0, const double* q1,
        coord_index_t dim
    ) {
        cnt_side2_exact++;

        const expansion& l1 = expansion_sq_dist(p1, p0, dim);
        const expansion& l2 = expansion_sq_dist(p2, p0, dim);

        const expansion& a10 = expansion_dot_at(p1, q0, p0, dim).scale_fast(2.0);
        const expansion& a11 = expansion_dot_at(p1, q1, p0, dim).scale_fast(2.0);
        const expansion& a20 = expansion_dot_at(p2, q0, p0, dim).scale_fast(2.0);
        const expansion& a21 = expansion_dot_at(p2, q1, p0, dim).scale_fast(2.0);

        const expansion& Delta = expansion_diff(a11, a10);

        Sign Delta_sign = Delta.sign();
        // Should not occur with symbolic
        // perturbation done at previous steps.
        geo_assert(Delta_sign != ZERO);

        //       [ Lambda0 ]   [ -1 ]        [  a11 ]
        // Delta [         ] = [    ] * l1 + [      ]
        //       [ Lambda1 ]   [  1 ]        [ -a10 ]

        const expansion& DeltaLambda0 = expansion_diff(a11, l1);
        const expansion& DeltaLambda1 = expansion_diff(l1, a10);

        // r = Delta*l2 - ( a20*DeltaLambda0 + a21*DeltaLambda1 )

        const expansion& r0 = expansion_product(Delta, l2);
        const expansion& r1 = expansion_product(a20, DeltaLambda0).negate();
        const expansion& r2 = expansion_product(a21, DeltaLambda1).negate();
        const expansion& r = expansion_sum3(r0, r1, r2);

        Sign r_sign = r.sign();

        // Statistics
        len_side2_num = std::max(len_side2_num, r.length());
        len_side2_denom = std::max(len_side2_denom, Delta.length());

        // Simulation of Simplicity (symbolic perturbation)
        if(r_sign == ZERO) {
            cnt_side2_SOS++;
            const double* p_sort[3];
            p_sort[0] = p0;
            p_sort[1] = p1;
            p_sort[2] = p2;
            std::sort(p_sort, p_sort + 3);
            for(index_t i = 0; i < 3; ++i) {
                if(p_sort[i] == p0) {
                    const expansion& z1 = expansion_diff(Delta, a21);
                    const expansion& z = expansion_sum(z1, a20);
                    Sign z_sign = z.sign();
                    len_side2_SOS = std::max(len_side2_SOS, z.length());
                    if(z_sign != ZERO) {
                        return Sign(Delta_sign * z_sign);
                    }
                }
                if(p_sort[i] == p1) {
                    const expansion& z = expansion_diff(a21, a20);
                    Sign z_sign = z.sign();
                    len_side2_SOS = std::max(len_side2_SOS, z.length());
                    if(z_sign != ZERO) {
                        return Sign(Delta_sign * z_sign);
                    }
                }
                if(p_sort[i] == p2) {
                    return NEGATIVE;
                }
            }
            geo_assert_not_reached;
        }

        return Sign(Delta_sign * r_sign);
    }

    /**
     * \brief Implements side2() in 3d.
     */
    Sign side2_3d_SOS(
        const double* p0, const double* p1, const double* p2,
        const double* q0, const double* q1
    ) {
        Sign result = Sign(side2_3d_filter(p0, p1, p2, q0, q1));
        if(result == ZERO) {
            result = side2_exact_SOS(p0, p1, p2, q0, q1, 3);
        }
        return result;
    }

    /**
     * \brief Implements side2() in 4d.
     */
    Sign side2_4d_SOS(
        const double* p0, const double* p1, const double* p2,
        const double* q0, const double* q1
    ) {
        Sign result = Sign(side2_4d_filter(p0, p1, p2, q0, q1));
        if(result == ZERO) {
            result = side2_exact_SOS(p0, p1, p2, q0, q1, 4);
        }
        return result;
    }

    /**
     * \brief Implements side2() in 6d.
     */
    Sign side2_6d_SOS(
        const double* p0, const double* p1, const double* p2,
        const double* q0, const double* q1
    ) {
        Sign result = Sign(side2_6d_filter(p0, p1, p2, q0, q1));
        if(result == ZERO) {
            result = side2_exact_SOS(p0, p1, p2, q0, q1, 6);
        }
        return result;
    }

    /**
     * \brief Implements side2() in 7d.
     */
    Sign side2_7d_SOS(
        const double* p0, const double* p1, const double* p2,
        const double* q0, const double* q1
    ) {
        Sign result = Sign(side2_7d_filter(p0, p1, p2, q0, q1));
        if(result == ZERO) {
            result = side2_exact_SOS(p0, p1, p2, q0, q1, 7);
        }
        return result;
    }

    // ================= side3 =========================================

    /**
     * \brief Exact implementation of the side3() predicate using low-level
     *  exact arithmetics API (expansion class).
     */
    Sign side3_exact_SOS(
        const double* p0, const double* p1, const double* p2, const double* p3,
        const double* q0, const double* q1, const double* q2,
        coord_index_t dim
    ) {
        cnt_side3_exact++;

        const expansion& l1 = expansion_sq_dist(p1, p0, dim);
        const expansion& l2 = expansion_sq_dist(p2, p0, dim);
        const expansion& l3 = expansion_sq_dist(p3, p0, dim);

        const expansion& a10 = expansion_dot_at(p1, q0, p0, dim).scale_fast(2.0);
        const expansion& a11 = expansion_dot_at(p1, q1, p0, dim).scale_fast(2.0);
        const expansion& a12 = expansion_dot_at(p1, q2, p0, dim).scale_fast(2.0);
        const expansion& a20 = expansion_dot_at(p2, q0, p0, dim).scale_fast(2.0);
        const expansion& a21 = expansion_dot_at(p2, q1, p0, dim).scale_fast(2.0);
        const expansion& a22 = expansion_dot_at(p2, q2, p0, dim).scale_fast(2.0);

        const expansion& a30 = expansion_dot_at(p3, q0, p0, dim).scale_fast(2.0);
        const expansion& a31 = expansion_dot_at(p3, q1, p0, dim).scale_fast(2.0);
        const expansion& a32 = expansion_dot_at(p3, q2, p0, dim).scale_fast(2.0);

        // [ b00 b01 b02 ]           [  1   1   1  ]-1
        // [ b10 b11 b12 ] = Delta * [ a10 a11 a12 ]
        // [ b20 b21 b22 ]           [ a20 a21 a22 ]

        const expansion& b00 = expansion_det2x2(a11, a12, a21, a22);
        const expansion& b01 = expansion_diff(a21, a22);
        const expansion& b02 = expansion_diff(a12, a11);
        const expansion& b10 = expansion_det2x2(a12, a10, a22, a20);
        const expansion& b11 = expansion_diff(a22, a20);
        const expansion& b12 = expansion_diff(a10, a12);
        const expansion& b20 = expansion_det2x2(a10, a11, a20, a21);
        const expansion& b21 = expansion_diff(a20, a21);
        const expansion& b22 = expansion_diff(a11, a10);

        const expansion& Delta = expansion_sum3(b00, b10, b20);
        Sign Delta_sign = Delta.sign();
        // Should not occur with symbolic
        // perturbation done at previous steps.
        geo_assert(Delta_sign != ZERO);

        //       [ Lambda0 ]   [ b01 b02 ]   [ l1 ]   [ b00 ]
        // Delta [ Lambda1 ] = [ b11 b12 ] * [    ] + [ b10 ]
        //       [ Lambda2 ]   [ b21 b22 ]   [ l2 ]   [ b20 ]

        const expansion& b01_l1 = expansion_product(b01, l1);
        const expansion& b02_l2 = expansion_product(b02, l2);
        const expansion& DeltaLambda0 = expansion_sum3(b01_l1, b02_l2, b00);

        const expansion& b11_l1 = expansion_product(b11, l1);
        const expansion& b12_l2 = expansion_product(b12, l2);
        const expansion& DeltaLambda1 = expansion_sum3(b11_l1, b12_l2, b10);

        const expansion& b21_l1 = expansion_product(b21, l1);
        const expansion& b22_l2 = expansion_product(b22, l2);
        const expansion& DeltaLambda2 = expansion_sum3(b21_l1, b22_l2, b20);

        // r = Delta*l3-(a30*DeltaLambda0+a31*DeltaLambda1+a32*DeltaLambda2)

        const expansion& r0 = expansion_product(Delta, l3);
        const expansion& r1 = expansion_product(a30, DeltaLambda0).negate();
        const expansion& r2 = expansion_product(a31, DeltaLambda1).negate();
        const expansion& r3 = expansion_product(a32, DeltaLambda2).negate();
        const expansion& r = expansion_sum4(r0, r1, r2, r3);
        Sign r_sign = r.sign();

        // Statistics
        len_side3_num = std::max(len_side3_num, r.length());
        len_side3_denom = std::max(len_side3_denom, Delta.length());

        // Simulation of Simplicity (symbolic perturbation)
        if(r_sign == ZERO) {
            cnt_side3_SOS++;
            const double* p_sort[4];
            p_sort[0] = p0;
            p_sort[1] = p1;
            p_sort[2] = p2;
            p_sort[3] = p3;
            std::sort(p_sort, p_sort + 4);
            for(index_t i = 0; i < 4; ++i) {
                if(p_sort[i] == p0) {
                    const expansion& z1_0 = expansion_sum(b01, b02);
                    const expansion& z1 = expansion_product(a30, z1_0).negate();
                    const expansion& z2_0 = expansion_sum(b11, b12);
                    const expansion& z2 = expansion_product(a31, z2_0).negate();
                    const expansion& z3_0 = expansion_sum(b21, b22);
                    const expansion& z3 = expansion_product(a32, z3_0).negate();
                    const expansion& z = expansion_sum4(Delta, z1, z2, z3);
                    Sign z_sign = z.sign();
                    len_side3_SOS = std::max(len_side3_SOS, z.length());
                    if(z_sign != ZERO) {
                        return Sign(Delta_sign * z_sign);
                    }
                } else if(p_sort[i] == p1) {
                    const expansion& z1 = expansion_product(a30, b01);
                    const expansion& z2 = expansion_product(a31, b11);
                    const expansion& z3 = expansion_product(a32, b21);
                    const expansion& z = expansion_sum3(z1, z2, z3);
                    Sign z_sign = z.sign();
                    len_side3_SOS = std::max(len_side3_SOS, z.length());
                    if(z_sign != ZERO) {
                        return Sign(Delta_sign * z_sign);
                    }
                } else if(p_sort[i] == p2) {
                    const expansion& z1 = expansion_product(a30, b02);
                    const expansion& z2 = expansion_product(a31, b12);
                    const expansion& z3 = expansion_product(a32, b22);
                    const expansion& z = expansion_sum3(z1, z2, z3);
                    Sign z_sign = z.sign();
                    len_side3_SOS = std::max(len_side3_SOS, z.length());
                    if(z_sign != ZERO) {
                        return Sign(Delta_sign * z_sign);
                    }
                } else if(p_sort[i] == p3) {
                    return NEGATIVE;
                }
            }
            geo_assert_not_reached;
        }
        return Sign(Delta_sign * r_sign);
    }

    /**
     * \brief Implements side3() in 3d.
     */
    Sign side3_3d_SOS(
        const double* p0, const double* p1, const double* p2, const double* p3,
        const double* q0, const double* q1, const double* q2
    ) {
        Sign result = Sign(side3_3d_filter(p0, p1, p2, p3, q0, q1, q2));
        if(result == ZERO) {
            result = side3_exact_SOS(p0, p1, p2, p3, q0, q1, q2, 3);
        }
        return result;
    }

    /**
     * \brief Implements side3() in 4d.
     */
    Sign side3_4d_SOS(
        const double* p0, const double* p1, const double* p2, const double* p3,
        const double* q0, const double* q1, const double* q2
    ) {
        Sign result = Sign(side3_4d_filter(p0, p1, p2, p3, q0, q1, q2));
        if(result == ZERO) {
            result = side3_exact_SOS(p0, p1, p2, p3, q0, q1, q2, 4);
        }
        return result;
    }

    /**
     * \brief Implements side3() in 6d.
     */
    Sign side3_6d_SOS(
        const double* p0, const double* p1, const double* p2, const double* p3,
        const double* q0, const double* q1, const double* q2
    ) {
        Sign result = Sign(side3_6d_filter(p0, p1, p2, p3, q0, q1, q2));
        if(result == ZERO) {
            result = side3_exact_SOS(p0, p1, p2, p3, q0, q1, q2, 6);
        }
        return result;
    }

    /**
     * \brief Implements side3() in 7d.
     */
    Sign side3_7d_SOS(
        const double* p0, const double* p1, const double* p2, const double* p3,
        const double* q0, const double* q1, const double* q2
    ) {
        Sign result = Sign(side3_7d_filter(p0, p1, p2, p3, q0, q1, q2));
        if(result == ZERO) {
            result = side3_exact_SOS(p0, p1, p2, p3, q0, q1, q2, 7);
        }
        return result;
    }

    // ================= side4 =========================================

    /**
     * \brief Exact implementation of the side4_3d_SOS() predicate using low-level
     *  exact arithmetics API (expansiont class).
     */
    Sign side4_3d_exact_SOS(
        const double* p0, const double* p1, const double* p2, const double* p3,
        const double* p4, bool sos = true
    ) {
        cnt_side4_exact++;

        const expansion& a11 = expansion_diff(p1[0], p0[0]);
        const expansion& a12 = expansion_diff(p1[1], p0[1]);
        const expansion& a13 = expansion_diff(p1[2], p0[2]);
        const expansion& a14 = expansion_sq_dist(p1, p0, 3).negate();

        const expansion& a21 = expansion_diff(p2[0], p0[0]);
        const expansion& a22 = expansion_diff(p2[1], p0[1]);
        const expansion& a23 = expansion_diff(p2[2], p0[2]);
        const expansion& a24 = expansion_sq_dist(p2, p0, 3).negate();

        const expansion& a31 = expansion_diff(p3[0], p0[0]);
        const expansion& a32 = expansion_diff(p3[1], p0[1]);
        const expansion& a33 = expansion_diff(p3[2], p0[2]);
        const expansion& a34 = expansion_sq_dist(p3, p0, 3).negate();

        const expansion& a41 = expansion_diff(p4[0], p0[0]);
        const expansion& a42 = expansion_diff(p4[1], p0[1]);
        const expansion& a43 = expansion_diff(p4[2], p0[2]);
        const expansion& a44 = expansion_sq_dist(p4, p0, 3).negate();

        // Note: we could probably reuse some of the 2x2 co-factors
        // (but for now I'd rather keep this form that is easier to
        //  read ... and to debug if need be !)
        const expansion& Delta1 = expansion_det3x3(
            a21, a22, a23,
            a31, a32, a33,
            a41, a42, a43
        );
        const expansion& Delta2 = expansion_det3x3(
            a11, a12, a13,
            a31, a32, a33,
            a41, a42, a43
        );
        const expansion& Delta3 = expansion_det3x3(
            a11, a12, a13,
            a21, a22, a23,
            a41, a42, a43
        );
        const expansion& Delta4 = expansion_det3x3(
            a11, a12, a13,
            a21, a22, a23,
            a31, a32, a33
        );

        Sign Delta4_sign = Delta4.sign();
        geo_assert(Delta4_sign != ZERO);

        const expansion& r_1 = expansion_product(Delta1, a14);
        const expansion& r_2 = expansion_product(Delta2, a24).negate();
        const expansion& r_3 = expansion_product(Delta3, a34);
        const expansion& r_4 = expansion_product(Delta4, a44).negate();
        const expansion& r = expansion_sum4(r_1, r_2, r_3, r_4);

        Sign r_sign = r.sign();

        // Statistics
        len_side4_num = std::max(len_side4_num, r.length());
        len_side4_denom = std::max(len_side4_denom, Delta1.length());

        // Simulation of Simplicity (symbolic perturbation)
        if(sos && r_sign == ZERO) {
            cnt_side4_SOS++;
            const double* p_sort[5];
            p_sort[0] = p0;
            p_sort[1] = p1;
            p_sort[2] = p2;
            p_sort[3] = p3;
            p_sort[4] = p4;
            std::sort(p_sort, p_sort + 5);
            for(index_t i = 0; i < 5; ++i) {
                if(p_sort[i] == p0) {
                    const expansion& z1 = expansion_diff(Delta2, Delta1);
                    const expansion& z2 = expansion_diff(Delta4, Delta3);
                    const expansion& z = expansion_sum(z1, z2);
                    Sign z_sign = z.sign();
                    len_side4_SOS = std::max(len_side4_SOS, z.length());
                    if(z_sign != ZERO) {
                        return Sign(Delta4_sign * z_sign);
                    }
                } else if(p_sort[i] == p1) {
                    Sign Delta1_sign = Delta1.sign();
                    if(Delta1_sign != ZERO) {
                        len_side4_SOS = std::max(len_side4_SOS, Delta1.length());
                        return Sign(Delta4_sign * Delta1_sign);
                    }
                } else if(p_sort[i] == p2) {
                    Sign Delta2_sign = Delta2.sign();
                    if(Delta2_sign != ZERO) {
                        len_side4_SOS = std::max(len_side4_SOS, Delta2.length());
                        return Sign(-Delta4_sign * Delta2_sign);
                    }
                } else if(p_sort[i] == p3) {
                    Sign Delta3_sign = Delta3.sign();
                    if(Delta3_sign != ZERO) {
                        len_side4_SOS = std::max(len_side4_SOS, Delta3.length());
                        return Sign(Delta4_sign * Delta3_sign);
                    }
                } else if(p_sort[i] == p4) {
                    return NEGATIVE;
                }
            }
        }
        return Sign(Delta4_sign * r_sign);
    }

    /**
     * \brief Exact implementation of the side4() predicate using low-level
     *  exact arithmetics API (expansion class).
     */
    Sign side4_exact_SOS(
        const double* p0, const double* p1, const double* p2, const double* p3,
        const double* p4,
        const double* q0, const double* q1, const double* q2, const double* q3,
        coord_index_t dim
    ) {
        cnt_side4_exact++;

        const expansion& l1 = expansion_sq_dist(p1, p0, dim);
        const expansion& l2 = expansion_sq_dist(p2, p0, dim);
        const expansion& l3 = expansion_sq_dist(p3, p0, dim);
        const expansion& l4 = expansion_sq_dist(p4, p0, dim);

        const expansion& a10 = expansion_dot_at(p1, q0, p0, dim).scale_fast(2.0);
        const expansion& a11 = expansion_dot_at(p1, q1, p0, dim).scale_fast(2.0);
        const expansion& a12 = expansion_dot_at(p1, q2, p0, dim).scale_fast(2.0);
        const expansion& a13 = expansion_dot_at(p1, q3, p0, dim).scale_fast(2.0);

        const expansion& a20 = expansion_dot_at(p2, q0, p0, dim).scale_fast(2.0);
        const expansion& a21 = expansion_dot_at(p2, q1, p0, dim).scale_fast(2.0);
        const expansion& a22 = expansion_dot_at(p2, q2, p0, dim).scale_fast(2.0);
        const expansion& a23 = expansion_dot_at(p2, q3, p0, dim).scale_fast(2.0);

        const expansion& a30 = expansion_dot_at(p3, q0, p0, dim).scale_fast(2.0);
        const expansion& a31 = expansion_dot_at(p3, q1, p0, dim).scale_fast(2.0);
        const expansion& a32 = expansion_dot_at(p3, q2, p0, dim).scale_fast(2.0);
        const expansion& a33 = expansion_dot_at(p3, q3, p0, dim).scale_fast(2.0);

        const expansion& a40 = expansion_dot_at(p4, q0, p0, dim).scale_fast(2.0);
        const expansion& a41 = expansion_dot_at(p4, q1, p0, dim).scale_fast(2.0);
        const expansion& a42 = expansion_dot_at(p4, q2, p0, dim).scale_fast(2.0);
        const expansion& a43 = expansion_dot_at(p4, q3, p0, dim).scale_fast(2.0);

        // [ b00 b01 b02 b03 ]           [  1   1   1   1  ]-1
        // [ b10 b11 b12 b13 ]           [ a10 a11 a12 a13 ]
        // [ b20 b21 b22 b23 ] = Delta * [ a20 a21 a22 a23 ]
        // [ b30 b31 b32 b33 ]           [ a30 a31 a32 a33 ]

        // Note: we could probably reuse some of the co-factors
        // (but for now I'd rather keep this form that is easier to
        //  read ... and to debug if need be !)

        const expansion& b00 = expansion_det3x3(a11, a12, a13, a21, a22, a23, a31, a32, a33);
        const expansion& b01 = expansion_det_111_2x3(a21, a22, a23, a31, a32, a33).negate();
        const expansion& b02 = expansion_det_111_2x3(a11, a12, a13, a31, a32, a33);
        const expansion& b03 = expansion_det_111_2x3(a11, a12, a13, a21, a22, a23).negate();

        const expansion& b10 = expansion_det3x3(a10, a12, a13, a20, a22, a23, a30, a32, a33).negate();
        const expansion& b11 = expansion_det_111_2x3(a20, a22, a23, a30, a32, a33);
        const expansion& b12 = expansion_det_111_2x3(a10, a12, a13, a30, a32, a33).negate();
        const expansion& b13 = expansion_det_111_2x3(a10, a12, a13, a20, a22, a23);

        const expansion& b20 = expansion_det3x3(a10, a11, a13, a20, a21, a23, a30, a31, a33);
        const expansion& b21 = expansion_det_111_2x3(a20, a21, a23, a30, a31, a33).negate();
        const expansion& b22 = expansion_det_111_2x3(a10, a11, a13, a30, a31, a33);
        const expansion& b23 = expansion_det_111_2x3(a10, a11, a13, a20, a21, a23).negate();

        const expansion& b30 = expansion_det3x3(a10, a11, a12, a20, a21, a22, a30, a31, a32).negate();
        const expansion& b31 = expansion_det_111_2x3(a20, a21, a22, a30, a31, a32);
        const expansion& b32 = expansion_det_111_2x3(a10, a11, a12, a30, a31, a32).negate();
        const expansion& b33 = expansion_det_111_2x3(a10, a11, a12, a20, a21, a22);

        const expansion& Delta = expansion_sum4(b00, b10, b20, b30);
        Sign Delta_sign = Delta.sign();
        geo_assert(Delta_sign != ZERO);

        //       [ Lambda0 ]   [ b01 b02 b03 ]   [ l1 ]   [ b00 ]
        //       [ Lambda1 ]   [ b11 b12 b13 ]   [ l2 ]   [ b10 ]
        // Delta [ Lambda2 ] = [ b21 b22 b23 ] * [ l3 ] + [ b20 ]
        //       [ Lambda3 ]   [ b31 b32 b33 ]   [ l4 ]   [ b30 ]

        const expansion& b01_l1 = expansion_product(b01, l1);
        const expansion& b02_l2 = expansion_product(b02, l2);
        const expansion& b03_l3 = expansion_product(b03, l3);
        const expansion& DeltaLambda0 = expansion_sum4(b01_l1, b02_l2, b03_l3, b00);

        const expansion& b11_l1 = expansion_product(b11, l1);
        const expansion& b12_l2 = expansion_product(b12, l2);
        const expansion& b13_l3 = expansion_product(b13, l3);
        const expansion& DeltaLambda1 = expansion_sum4(b11_l1, b12_l2, b13_l3, b10);

        const expansion& b21_l1 = expansion_product(b21, l1);
        const expansion& b22_l2 = expansion_product(b22, l2);
        const expansion& b23_l3 = expansion_product(b23, l3);
        const expansion& DeltaLambda2 = expansion_sum4(b21_l1, b22_l2, b23_l3, b20);

        const expansion& b31_l1 = expansion_product(b31, l1);
        const expansion& b32_l2 = expansion_product(b32, l2);
        const expansion& b33_l3 = expansion_product(b33, l3);
        const expansion& DeltaLambda3 = expansion_sum4(b31_l1, b32_l2, b33_l3, b30);

        // r = Delta*l4 - (
        //    a40*DeltaLambda0+
        //    a41*DeltaLambda1+
        //    a42*DeltaLambda2+
        //    a43*DeltaLambda3
        // )

        const expansion& r0 = expansion_product(Delta, l4);
        const expansion& r1 = expansion_product(a40, DeltaLambda0);
        const expansion& r2 = expansion_product(a41, DeltaLambda1);
        const expansion& r3 = expansion_product(a42, DeltaLambda2);
        const expansion& r4 = expansion_product(a43, DeltaLambda3);
        const expansion& r1234 = expansion_sum4(r1, r2, r3, r4);
        const expansion& r = expansion_diff(r0, r1234);
        Sign r_sign = r.sign();

        // Simulation of Simplicity (symbolic perturbation)
        if(r_sign == ZERO) {
            cnt_side4_SOS++;
            const double* p_sort[5];
            p_sort[0] = p0;
            p_sort[1] = p1;
            p_sort[2] = p2;
            p_sort[3] = p3;
            p_sort[4] = p4;
            std::sort(p_sort, p_sort + 5);
            for(index_t i = 0; i < 5; ++i) {
                if(p_sort[i] == p0) {
                    const expansion& z1_0 = expansion_sum3(b01, b02, b03);
                    const expansion& z1 = expansion_product(a30, z1_0);
                    const expansion& z2_0 = expansion_sum3(b11, b12, b13);
                    const expansion& z2 = expansion_product(a31, z2_0);
                    const expansion& z3_0 = expansion_sum3(b21, b22, b23);
                    const expansion& z3 = expansion_product(a32, z3_0);
                    const expansion& z4_0 = expansion_sum3(b31, b32, b33);
                    const expansion& z4 = expansion_product(a33, z4_0);
                    const expansion& z1234 = expansion_sum4(z1, z2, z3, z4);
                    const expansion& z = expansion_diff(Delta, z1234);
                    Sign z_sign = z.sign();
                    len_side4_SOS = std::max(len_side4_SOS, z.length());
                    if(z_sign != ZERO) {
                        return Sign(Delta_sign * z_sign);
                    }
                } else if(p_sort[i] == p1) {
                    const expansion& z1 = expansion_product(a30, b01);
                    const expansion& z2 = expansion_product(a31, b11);
                    const expansion& z3 = expansion_product(a32, b21);
                    const expansion& z4 = expansion_product(a33, b31);
                    const expansion& z = expansion_sum4(z1, z2, z3, z4);
                    Sign z_sign = z.sign();
                    len_side4_SOS = std::max(len_side4_SOS, z.length());
                    if(z_sign != ZERO) {
                        return Sign(Delta_sign * z_sign);
                    }
                } else if(p_sort[i] == p2) {
                    const expansion& z1 = expansion_product(a30, b02);
                    const expansion& z2 = expansion_product(a31, b12);
                    const expansion& z3 = expansion_product(a32, b22);
                    const expansion& z4 = expansion_product(a33, b32);
                    const expansion& z = expansion_sum4(z1, z2, z3, z4);
                    Sign z_sign = z.sign();
                    len_side4_SOS = std::max(len_side4_SOS, z.length());
                    if(z_sign != ZERO) {
                        return Sign(Delta_sign * z_sign);
                    }
                } else if(p_sort[i] == p3) {
                    const expansion& z1 = expansion_product(a30, b03);
                    const expansion& z2 = expansion_product(a31, b13);
                    const expansion& z3 = expansion_product(a32, b23);
                    const expansion& z4 = expansion_product(a33, b33);
                    const expansion& z = expansion_sum4(z1, z2, z3, z4);
                    Sign z_sign = z.sign();
                    len_side4_SOS = std::max(len_side4_SOS, z.length());
                    if(z_sign != ZERO) {
                        return Sign(Delta_sign * z_sign);
                    }
                } else if(p_sort[i] == p4) {
                    return NEGATIVE;
                }
            }
            geo_assert_not_reached;
        }
        return Sign(r_sign * Delta_sign);
    }

    /**
     * \brief Implements side4() in 4d.
     */
    Sign side4_4d_SOS(
        const double* p0,
        const double* p1, const double* p2, const double* p3, const double* p4,
        const double* q0, const double* q1, const double* q2, const double* q3
    ) {
        Sign result = Sign(side4_4d_filter(p0, p1, p2, p3, p4, q0, q1, q2, q3));
        if(result == ZERO) {
            result = side4_exact_SOS(p0, p1, p2, p3, p4, q0, q1, q2, q3, 4);
        }
        return result;
    }

    /**
     * \brief Implements side4() in 6d.
     */
    Sign side4_6d_SOS(
        const double* p0,
        const double* p1, const double* p2, const double* p3, const double* p4,
        const double* q0, const double* q1, const double* q2, const double* q3
    ) {
        Sign result = Sign(side4_6d_filter(p0, p1, p2, p3, p4, q0, q1, q2, q3));
        if(result == ZERO) {
            result = side4_exact_SOS(p0, p1, p2, p3, p4, q0, q1, q2, q3, 6);
        }
        return result;
    }

    /**
     * \brief Implements side4() in 7d.
     */
    Sign side4_7d_SOS(
        const double* p0,
        const double* p1, const double* p2, const double* p3, const double* p4,
        const double* q0, const double* q1, const double* q2, const double* q3
    ) {
        Sign result = Sign(side4_7d_filter(p0, p1, p2, p3, p4, q0, q1, q2, q3));
        if(result == ZERO) {
            result = side4_exact_SOS(p0, p1, p2, p3, p4, q0, q1, q2, q3, 7);
        }
        return result;
    }

    // ============ orient2d ==============================================

    Sign orient_2d_exact(
        const double* p0, const double* p1, const double* p2
    ) {
        cnt_orient2d_exact++;

        const expansion& a11 = expansion_diff(p1[0], p0[0]);
        const expansion& a12 = expansion_diff(p1[1], p0[1]);

        const expansion& a21 = expansion_diff(p2[0], p0[0]);
        const expansion& a22 = expansion_diff(p2[1], p0[1]);

        const expansion& Delta = expansion_det2x2(
            a11, a12, a21, a22
        );

        len_orient2d = std::max(len_orient2d, Delta.length());

        return Delta.sign();
    }


    // ============ orient3d ==============================================

    Sign orient_3d_exact(
        const double* p0, const double* p1,
        const double* p2, const double* p3
    ) {
        cnt_orient3d_exact++;

        const expansion& a11 = expansion_diff(p1[0], p0[0]);
        const expansion& a12 = expansion_diff(p1[1], p0[1]);
        const expansion& a13 = expansion_diff(p1[2], p0[2]);

        const expansion& a21 = expansion_diff(p2[0], p0[0]);
        const expansion& a22 = expansion_diff(p2[1], p0[1]);
        const expansion& a23 = expansion_diff(p2[2], p0[2]);

        const expansion& a31 = expansion_diff(p3[0], p0[0]);
        const expansion& a32 = expansion_diff(p3[1], p0[1]);
        const expansion& a33 = expansion_diff(p3[2], p0[2]);

        const expansion& Delta = expansion_det3x3(
            a11, a12, a13, a21, a22, a23, a31, a32, a33
        );

        len_orient3d = std::max(len_orient3d, Delta.length());

        return Delta.sign();
    }

    Sign side4h_3d_exact_SOS(
        const double* p0, const double* p1,
        const double* p2, const double* p3, const double* p4,
        double h0, double h1, double h2, double h3, double h4,
        bool sos = true
    ) {
        cnt_orient4d_exact++;

        const expansion& a11 = expansion_diff(p1[0], p0[0]);
        const expansion& a12 = expansion_diff(p1[1], p0[1]);
        const expansion& a13 = expansion_diff(p1[2], p0[2]);
        const expansion& a14 = expansion_diff(h0,h1);

        const expansion& a21 = expansion_diff(p2[0], p0[0]);
        const expansion& a22 = expansion_diff(p2[1], p0[1]);
        const expansion& a23 = expansion_diff(p2[2], p0[2]);
        const expansion& a24 = expansion_diff(h0,h2);

        const expansion& a31 = expansion_diff(p3[0], p0[0]);
        const expansion& a32 = expansion_diff(p3[1], p0[1]);
        const expansion& a33 = expansion_diff(p3[2], p0[2]);
        const expansion& a34 = expansion_diff(h0,h3);

        const expansion& a41 = expansion_diff(p4[0], p0[0]);
        const expansion& a42 = expansion_diff(p4[1], p0[1]);
        const expansion& a43 = expansion_diff(p4[2], p0[2]);
        const expansion& a44 = expansion_diff(h0,h4);

        // Note: we could probably reuse some of the 2x2 co-factors
        // (but for now I'd rather keep this form that is easier to
        //  read ... and to debug if need be !)
        const expansion& Delta1 = expansion_det3x3(
            a21, a22, a23,
            a31, a32, a33,
            a41, a42, a43
        );
        const expansion& Delta2 = expansion_det3x3(
            a11, a12, a13,
            a31, a32, a33,
            a41, a42, a43
        );
        const expansion& Delta3 = expansion_det3x3(
            a11, a12, a13,
            a21, a22, a23,
            a41, a42, a43
        );
        const expansion& Delta4 = expansion_det3x3(
            a11, a12, a13,
            a21, a22, a23,
            a31, a32, a33
        );

        Sign Delta4_sign = Delta4.sign();
        geo_assert(Delta4_sign != ZERO);

        const expansion& r_1 = expansion_product(Delta1, a14);
        const expansion& r_2 = expansion_product(Delta2, a24).negate();
        const expansion& r_3 = expansion_product(Delta3, a34);
        const expansion& r_4 = expansion_product(Delta4, a44).negate();
        const expansion& r = expansion_sum4(r_1, r_2, r_3, r_4);

        Sign r_sign = r.sign();

        // Statistics
        len_orient4d_num = std::max(len_orient4d_num, r.length());
        len_orient4d_denom = std::max(len_orient4d_denom, Delta1.length());

        // Simulation of Simplicity (symbolic perturbation)
        if(sos && r_sign == ZERO) {
            cnt_orient4d_SOS++;
            const double* p_sort[5];
            p_sort[0] = p0;
            p_sort[1] = p1;
            p_sort[2] = p2;
            p_sort[3] = p3;
            p_sort[4] = p4;
            std::sort(p_sort, p_sort + 5);
            for(index_t i = 0; i < 5; ++i) {
                if(p_sort[i] == p0) {
                    const expansion& z1 = expansion_diff(Delta2, Delta1);
                    const expansion& z2 = expansion_diff(Delta4, Delta3);
                    const expansion& z = expansion_sum(z1, z2);
                    Sign z_sign = z.sign();
                    len_orient4d_SOS = std::max(len_orient4d_SOS, z.length());
                    if(z_sign != ZERO) {
                        return Sign(Delta4_sign * z_sign);
                    }
                } else if(p_sort[i] == p1) {
                    Sign Delta1_sign = Delta1.sign();
                    if(Delta1_sign != ZERO) {
                        len_orient4d_SOS = std::max(len_orient4d_SOS, Delta1.length());
                        return Sign(Delta4_sign * Delta1_sign);
                    }
                } else if(p_sort[i] == p2) {
                    Sign Delta2_sign = Delta2.sign();
                    if(Delta2_sign != ZERO) {
                        len_orient4d_SOS = std::max(len_orient4d_SOS, Delta2.length());
                        return Sign(-Delta4_sign * Delta2_sign);
                    }
                } else if(p_sort[i] == p3) {
                    Sign Delta3_sign = Delta3.sign();
                    if(Delta3_sign != ZERO) {
                        len_orient4d_SOS = std::max(len_orient4d_SOS, Delta3.length());
                        return Sign(Delta4_sign * Delta3_sign);
                    }
                } else if(p_sort[i] == p4) {
                    return NEGATIVE;
                }
            }
        }
        return Sign(Delta4_sign * r_sign);
    }

    // ================================ statistics ========================

    /**
     * \brief Returns the percentage that a number represents
     *   relative to another one.
     */
    inline double percent(index_t a, index_t b) {
        if(a == 0 && b == 0) {
            return 0;
        }
        return double(a * 100) / b;
    }

    /**
     * \brief Displays statistic counters for exact predicates
     * \param[in] name name of the predicate
     * \param[in] cnt1 total number of invocations
     * \param[in] cnt2 number of exact invocations
     */
    void show_stats_plain(
        const std::string& name, index_t cnt1, index_t cnt2
    ) {
        std::cerr << name << ": " 
            << "Tot:" << cnt1
            << " Exact:" << cnt2
            << std::endl;
        std::cerr << name << ": "
            << " Exact: " << percent(cnt2, cnt1) << "% "
            << std::endl;
    }

    /**
     * \brief Displays statistic counters for exact predicates
     * \param[in] name name of the predicate
     * \param[in] cnt1 total number of invocations
     * \param[in] cnt2 number of exact invocations
     * \param[in] cnt3 number of SOS invocations
     */
    void show_stats_sos(
        const std::string& name, index_t cnt1, index_t cnt2, index_t cnt3
    ) {
        std::cerr << name << ": " 
            << "Tot:" << cnt1
            << " Exact:" << cnt2
            << " SOS:" << cnt3 << std::endl;
        std::cerr << name << ": "
            << " Exact: " << percent(cnt2, cnt1) << "% "
            << " SOS: " << percent(cnt3, cnt1) << "% "
            << std::endl;
    }

    /**
     * \brief Displays statistic counters for exact predicates
     * \param[in] name name of the predicate
     * \param[in] cnt1 total number of invocations
     * \param[in] cnt2 number of exact invocations
     * \param[in] cnt3 number of SOS invocations
     * \param[in] len maximum length of the expansions during exact computations
     */
    void show_stats_sos(
        const std::string& name, index_t cnt1, index_t cnt2, index_t cnt3,
        index_t len
    ) {
        show_stats_sos(name, cnt1, cnt2, cnt3);
        std::cerr << "name: " << " Len: " << len << std::endl;
    }

    /**
     * \brief Displays statistic counters for exact predicates
     * \param[in] name name of the predicate
     * \param[in] cnt1 total number of invocations
     * \param[in] cnt2 number of exact invocations
     * \param[in] len maximum length of the expansions during exact computations
     */
    void show_stats_plain(
        const std::string& name, index_t cnt1, index_t cnt2,
        index_t len
    ) {
        show_stats_plain(name, cnt1, cnt2);
        std::cerr << name << ":  Len: " << len << std::endl;
    }

    /**
     * \brief Displays statistic counters for exact predicates
     * \param[in] name name of the predicate
     * \param[in] cnt1 total number of invocations
     * \param[in] cnt2 number of exact invocations
     * \param[in] cnt3 number of SOS invocations
     * \param[in] num_len maximum length of the expansions
     *   that represent the numerator
     * \param[in] denom_len maximum length of the expansions
     *   that represent the denominator
     * \param[in] SOS_len maximum length of the expansions
     *   that represent the SOS terms
     */
    void show_stats_sos(
        const std::string& name, index_t cnt1, index_t cnt2, index_t cnt3,
        index_t num_len, index_t denom_len, index_t SOS_len
    ) {
        show_stats_sos(name, cnt1, cnt2, cnt3);
        std::cerr << name << ": "
            << " Num len: " << num_len
            << " Denom len: " << denom_len
            << " SOS len: " << SOS_len
            << std::endl;
    }
}

/****************************************************************************/

namespace GEO {

    namespace PCK {

        Sign side1_SOS(
            const double* p0, const double* p1,
            const double* q0,
            coord_index_t DIM
        ) {
            cnt_side1_total++;
            switch(DIM) {
            case 3:
                return side1_3d_SOS(p0, p1, q0);
            case 4:
                return side1_4d_SOS(p0, p1, q0);
            case 6:
                return side1_6d_SOS(p0, p1, q0);
            case 7:
                return side1_7d_SOS(p0, p1, q0);
            }
            geo_assert_not_reached;
            return ZERO;
        }

        Sign side2_SOS(
            const double* p0, const double* p1, const double* p2,
            const double* q0, const double* q1,
            coord_index_t DIM
        ) {
            cnt_side2_total++;
            switch(DIM) {
            case 3:
                return side2_3d_SOS(p0, p1, p2, q0, q1);
            case 4:
                return side2_4d_SOS(p0, p1, p2, q0, q1);
            case 6:
                return side2_6d_SOS(p0, p1, p2, q0, q1);
            case 7:
                return side2_7d_SOS(p0, p1, p2, q0, q1);
            }
            geo_assert_not_reached;
            return ZERO;
        }

        Sign side3_SOS(
            const double* p0, const double* p1, const double* p2, const double* p3,
            const double* q0, const double* q1, const double* q2,
            coord_index_t DIM
        ) {
            cnt_side3_total++;
            switch(DIM) {
            case 3:
                return side3_3d_SOS(p0, p1, p2, p3, q0, q1, q2);
            case 4:
                return side3_4d_SOS(p0, p1, p2, p3, q0, q1, q2);
            case 6:
                return side3_6d_SOS(p0, p1, p2, p3, q0, q1, q2);
            case 7:
                return side3_7d_SOS(p0, p1, p2, p3, q0, q1, q2);
            }
            geo_assert_not_reached;
            return ZERO;
        }

        Sign side4_SOS(
            const double* p0,
            const double* p1, const double* p2, const double* p3, const double* p4,
            const double* q0, const double* q1, const double* q2, const double* q3,
            coord_index_t DIM
        ) {
            switch(DIM) {
            case 3:
                // 3d is a special case for side4()
                //   (intrinsic dim == ambient dim)
                // therefore embedding tet q0,q1,q2,q3 is not needed.
                // WARNING: cnt_side4_total is not incremented here,
                // because it is
                // incremented in side4_3d_SOS().
                return side4_3d_SOS(p0, p1, p2, p3, p4);
            case 4:
                cnt_side4_total++;
                return side4_4d_SOS(p0, p1, p2, p3, p4, q0, q1, q2, q3);
            case 6:
                cnt_side4_total++;
                return side4_6d_SOS(p0, p1, p2, p3, p4, q0, q1, q2, q3);
            case 7:
                cnt_side4_total++;
                return side4_7d_SOS(p0, p1, p2, p3, p4, q0, q1, q2, q3);
            }
            geo_assert_not_reached;
            return ZERO;
        }


        Sign side4_3d(
            const double* p0, const double* p1, const double* p2, const double* p3,
            const double* p4
        ) {
            cnt_side4_total++;
            Sign result = Sign(side4_3d_filter(p0, p1, p2, p3, p4));
            if(result == 0) {
                // last argument is false: do not apply symbolic perturbation
                result = side4_3d_exact_SOS(p0, p1, p2, p3, p4, false);
            }
            return result;
        }

        Sign side4_3d_SOS(
            const double* p0, const double* p1, const double* p2, const double* p3,
            const double* p4
        ) {
            cnt_side4_total++;
            Sign result = Sign(side4_3d_filter(p0, p1, p2, p3, p4));
            if(result == 0) {
                result = side4_3d_exact_SOS(p0, p1, p2, p3, p4);
            }
            return result;
        }


        Sign in_sphere_3d_SOS(
            const double* p0, const double* p1, const double* p2, const double* p3,
            const double* p4
        ) {
            // in_sphere_3d is simply implemented using side4_3d.
            // Both predicates are equivalent through duality as can
            // be easily seen:
            // side4_3d(p0,p1,p2,p3,p4) returns POSITIVE if
            //    d(q,p0) < d(q,p4)
            //    where q denotes the circumcenter of (p0,p1,p2,p3)
            // Note that d(q,p0) = R  (radius of circumscribed sphere)
            // In other words, side4_3d(p0,p1,p2,p3,p4) returns POSITIVE if
            //   d(q,p4) > R which means whenever p4 is not in the
            //   circumscribed sphere of (p0,p1,p2,p3).
            // Therefore in_sphere_3d(p0,p1,p2,p3,p4) = -side4_3d(p0,p1,p2,p3,p4)

            cnt_side4_total++;
            
            // This specialized filter supposes that orient_3d(p0,p1,p2,p3) > 0
            Sign result = Sign(in_sphere_3d_filter_optim(p0, p1, p2, p3, p4));
            if(result == 0) {
                result = side4_3d_exact_SOS(p0, p1, p2, p3, p4);
            }
            return Sign(-result);
        }


        Sign orient_2d(
            const double* p0, const double* p1, const double* p2
        ) {
            cnt_orient2d_total++;
            Sign result = Sign(orient_2d_filter(p0, p1, p2));
            if(result == 0) {
                result = orient_2d_exact(p0, p1, p2);
            }
            return result;
        }

        Sign orient_3d(
            const double* p0, const double* p1,
            const double* p2, const double* p3
            ) {
            cnt_orient3d_total++;
            Sign result = Sign(orient_3d_filter(p0, p1, p2, p3));
            if(result == 0) {
                result = orient_3d_exact(p0, p1, p2, p3);
            }
            return result;
        }


        Sign orient_4d(
            const double* p0, const double* p1,
            const double* p2, const double* p3, const double* p4,
            double h0, double h1, double h2, double h3, double h4
        ) {
            cnt_orient4d_total++;
            Sign result = Sign(
                side4h_3d_filter(
                    p0, p1, p2, p3, p4, h0, h1, h2, h3, h4
                    )
                );
            if(result == 0) {
                // last argument is false -> do not perturb.
                result = side4h_3d_exact_SOS(
                    p0, p1, p2, p3, p4, h0, h1, h2, h3, h4, false
                );
            }
            // orient_4d() is opposite to side4h()
            // (like in_sphere() that is opposite to side4())
            return Sign(-result);
        }

        Sign orient_4d_SOS(
            const double* p0, const double* p1,
            const double* p2, const double* p3, const double* p4,
            double h0, double h1, double h2, double h3, double h4
        ) {
            cnt_orient4d_total++;
            Sign result = Sign(
                side4h_3d_filter(
                    p0, p1, p2, p3, p4, h0, h1, h2, h3, h4
                    )
                );
            if(result == 0) {
                result = side4h_3d_exact_SOS(
                    p0, p1, p2, p3, p4, h0, h1, h2, h3, h4
                );
            }
            // orient_4d() is opposite to side4h()
            // (like in_sphere() that is opposite to side4())
            return Sign(-result);
        }

        void show_stats() {
            show_stats_plain(
                "orient2d",
                cnt_orient2d_total, cnt_orient2d_exact,
                len_orient2d
            );
            show_stats_plain(
                "orient3d",
                cnt_orient3d_total, cnt_orient3d_exact,
                len_orient3d
            );
            show_stats_sos(
                "orient4d",
                cnt_orient4d_total, cnt_orient4d_exact, cnt_orient4d_SOS,
                len_orient4d_num, len_orient4d_denom, len_orient4d_SOS
            );
            show_stats_sos(
                "side1",
                cnt_side1_total, cnt_side1_exact, cnt_side1_SOS,
                len_side1
            );
            show_stats_sos(
                "side2",
                cnt_side2_total, cnt_side2_exact, cnt_side2_SOS,
                len_side2_num, len_side2_denom, len_side2_SOS
            );
            show_stats_sos(
                "side3",
                cnt_side3_total, cnt_side3_exact, cnt_side3_SOS,
                len_side3_num, len_side3_denom, len_side3_SOS
            );
            show_stats_sos(
                "side4/insph.",
                cnt_side4_total, cnt_side4_exact, cnt_side4_SOS,
                len_side4_num, len_side4_denom, len_side4_SOS
            );
        }
    }
}

