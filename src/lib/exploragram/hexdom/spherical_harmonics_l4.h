/*
 *  OGF/Graphite: Geometry and Graphics Programming Library + Utilities
 *  Copyright (C) 2000-2015 INRIA - Project ALICE
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact for Graphite: Bruno Levy - Bruno.Levy@inria.fr
 *  Contact for this Plugin: Nicolas Ray - nicolas.ray@inria.fr
 *
 *     Project ALICE
 *     LORIA, INRIA Lorraine,
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX
 *     FRANCE
 *
 *  Note that the GNU General Public License does not permit incorporating
 *  the Software into proprietary programs.
 *
 * As an exception to the GPL, Graphite can be linked with the following
 * (non-GPL) libraries:
 *     Qt, tetgen, SuperLU, WildMagic and CGAL
 */

#ifndef H_HEXDOM_ALGO_SPHERICALHARMONICSL4_H
#define H_HEXDOM_ALGO_SPHERICALHARMONICSL4_H

#include <exploragram/basic/common.h>
#include <geogram/mesh/mesh.h>

namespace GEO {
    
    struct EXPLORAGRAM_API SphericalHarmonicL4 {
        vecng<9, Numeric::float64> coeff;

        // easy manipulations as a 9D vector
        SphericalHarmonicL4(const SphericalHarmonicL4& ref);
        SphericalHarmonicL4(const vecng<9, Numeric::float64>& p_coeff);    
        SphericalHarmonicL4(double *fv);
        SphericalHarmonicL4(double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7, double x8);
        SphericalHarmonicL4();
        double& operator[](index_t i) {
	    geo_debug_assert(i<9);
	    return coeff[i];	    
	}
        double norm();
        double operator *(const SphericalHarmonicL4 &other);
        SphericalHarmonicL4 operator -(const SphericalHarmonicL4 &other);
        SphericalHarmonicL4 operator *(double s);
        SphericalHarmonicL4 operator /(double s);
        SphericalHarmonicL4 operator +(const SphericalHarmonicL4 &v);



        double value(vec3 v);
        static double basis(index_t id, vec3 v);
        void Rz(double alpha);
        void Ry(double alpha);
        void Rx(double alpha);
        void euler_rot(vec3 rot_vec);
        SphericalHarmonicL4 Ex();
        SphericalHarmonicL4 Ey();
        SphericalHarmonicL4 Ez();

        static SphericalHarmonicL4 rest_frame() {
            return SphericalHarmonicL4(0, 0, 0, 0, std::sqrt(7. / 12.), 0, 0, 0, std::sqrt(5. / 12.));
        }
        mat3 project_mat3(double grad_threshold = 1e-3, double dot_threshold = 1e-5, vec3* euler_prev = NULL);

    };

    inline std::istream& operator>> (std::istream& input, SphericalHarmonicL4 &gna) {
	return input >> gna.coeff;
    }

    inline std::ostream& operator<< (std::ostream& output, const SphericalHarmonicL4 &gna) {
	return output << gna.coeff;
    }
}


#endif //__SPHERICALHARMONICSL4_H__

	
