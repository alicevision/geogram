/*
 *  OGF/Graphite: Geometry and Graphics Programming Library + Utilities
 *  Copyright (C) 2000-2009 INRIA - Project ALICE
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
 *  Contact: Bruno Levy - levy@loria.fr
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
 * As an exception to the GPL, Graphite can be linked with 
 *     the following (non-GPL) libraries:
 *     Qt, SuperLU, WildMagic and CGAL
 */

#include <exploragram/optimal_transport/optimal_transport_3d.h>
#include <exploragram/optimal_transport/linear_least_squares.h>

#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_AABB.h>

#include <geogram/voronoi/CVT.h>
#include <geogram/voronoi/generic_RVD_vertex.h>
#include <geogram/voronoi/RVD_polyhedron_callback.h>
#include <geogram/voronoi/generic_RVD_cell.h>
#include <geogram/delaunay/delaunay_3d.h>

#include <geogram/points/nn_search.h>
#include <geogram/numerics/optimizer.h>
#include <geogram/numerics/lbfgs_optimizers.h>

#include <geogram/basic/stopwatch.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/process.h>
#include <geogram/basic/progress.h>
#include <geogram/basic/permutation.h>
#include <geogram/basic/command_line.h>

#include <geogram/NL/nl.h>

#include <geogram/bibliography/bibliography.h>

#include <stack>
#include <iterator>

namespace {
    using namespace GEO;

    /**
     * \brief Computes the contribution of a polyhedron
     *  to the objective function minimized by a semi-discrete
     *  optimal transport map.
     */
    class OTMPolyhedronCallback : public RVDPolyhedronCallback {
    public:
	/**
	 * \brief OTMPolyhedronCallback constructor.
	 * \param[in] OTM a pointer to the OptimalTransportMap3d
	 */
	OTMPolyhedronCallback(OptimalTransportMap3d* OTM) :
	    OTM_(OTM),
	    Newton_step_(false),
	    n_(0),
	    w_(nil),
	    g_(nil),
	    mg_(nil)
	{
	    weighted_ =
		OTM->mesh().vertices.attributes().is_defined("weight");
	}

        /**
         * \brief Sets where centroids should be output.
         * \details This computes mass times centroid. The mass can
         *  be retreived (and used to divide) from the gradient.
         * \param[in] mg a pointer to the 3*nb_points coordinates
         *  of the centroids times the mass of the Laguerre cells
         */
        void set_Laguerre_centroids(double* mg) {
            mg_ = mg;
        }

	/**
	 * \brief Tests whether Laguerre centroids should be computed.
	 * \retval true if Laguerre centroids should be computed.
	 * \retval false otherwise.
	 */
	bool has_Laguerre_centroids() const {
	    return (mg_ != nil);
	}

	double* Laguerre_centroids() {
	    return mg_;
	}
	
        /**
         * \brief Sets the weight vector
         * \param[in] w a const pointer to the weight vector.
         * \param[in] n the number of weights in the weight vector.
         */
        void set_w(const double* w, index_t n) { 
            w_ = w;
            n_ = n;
        }

	
	/**
	 * \copydoc RVDPolyhedronCallback::operator()
	 */
	virtual void operator() (
	    index_t v,
	    index_t t,
	    const GEOGen::ConvexCell& C
	) const {

	    geo_argused(t);

	    bool use_weighted = (weighted_ || mg_ != nil);
	    

	    Thread* thread = Thread::current();
	    index_t current_thread_id = (thread == nil) ? 0 : thread->id();
	    
	    // Remember that ConvexCell is represented in dual form !
	    //   - ConvexCell's vertices are facets
	    //   - ConvexCell's triangles are vertices

	    //   Find the first vertex, used in weighted mode to tesselate
	    // the cells. In non-weighted mode, cells are tesselated around
	    // the center vertex v.
	    const GEOGen::Vertex* V0 = nil;
	    if(use_weighted) {
		for(index_t ct=0; ct < C.max_t(); ++ct) {
		    if(C.triangle_is_used(ct)) {
			V0 = &C.triangle_dual(ct);
			break;
		    }
		}
		if(V0 == nil) {
		    return;
		}
	    }
	    
	    for(index_t cv = 0; cv < C.max_v(); ++cv) {
		signed_index_t ct = C.vertex_triangle(cv);
		if(ct == -1) {
		    continue;
		}
		geo_debug_assert(C.triangle_is_used(index_t(ct)));

		index_t v_adj = index_t(-1);
		signed_index_t adjacent = C.vertex_id(cv);
		if(adjacent > 0) {
		    // Positive adjacent indices correspond to
		    // Voronoi seed - Voronoi seed link
		    v_adj = index_t(adjacent - 1);
		}
		bool on_border = (adjacent == 0);
		// Zero adjacent indices corresponds to
		// tet facet on border.

		//   In non-weighted mode, we only traverse the facets on
		// the boundary of the Voronoi cells.
		if(!use_weighted && v_adj == index_t(-1) && !on_border) {
		    continue;
		}
		
		GEOGen::ConvexCell::Corner first(
		    index_t(ct), C.find_triangle_vertex(index_t(ct), cv)
		);
		const GEOGen::Vertex* V1 = &C.triangle_dual(first.t);
		const GEOGen::Vertex* V2 = nil;
		const GEOGen::Vertex* V3 = nil;		
		GEOGen::ConvexCell::Corner c = first;
		do {
		    V2 = V3;
		    V3 = &C.triangle_dual(c.t);
		    if(V2 != nil && V3 != V1) {
			double F = use_weighted ?
			    eval_simplex_weighted(v, V0, V1, V2, V3, index_t(v_adj)) :
			    eval_simplex(v, V1, V2, V3, index_t(v_adj)) ;			    
			const_cast<OTMPolyhedronCallback*>(this)->
			    funcval_[current_thread_id] += F;
		    }
		    C.move_to_next_around_vertex(c);
		} while(c != first);
	    }
	}

	void set_Newton_step(bool Newton) {
	    Newton_step_ = Newton;
	}

	bool is_Newton_step() const {
	    return Newton_step_;
	}

	void set_nb_threads(index_t nb) {
	    funcval_.assign(nb, 0.0);	    
	}

	void set_g(double* g) {
	    g_ = g;
	}
	
	/**
	 * \brief Gets the computed value of the objective function.
	 * \details This sums the contributions of all threads.
	 * \retval the value of the objective function.
	 */
	double funcval() const {
	    double result = 0.0;
	    FOR(i,funcval_.size()) {
		result += funcval_[i];
	    }
	    return result;
	}
	
    protected:

	/**
	 * \brief Evaluates the contribution of a simplex in non-weighted mode.
	 * \details The considered simplex connects one of the Delaunay vertices 
	 *  (\p center_vertex_index) and three other vertices of the current cell
	 *  (\p V1, \p V2, \p V3).
	 * \param[in] center_vertex_index the index of the current Delaunay vertex
	 * \param[in] V1 , V2 , V3 three vertices of the current cell
	 * \param[in] v_adj the index of the Delaunay vertex of the cell on the other
	 *  side of (\p V1, \p V2, \p V3) or index_t(-1) if there is no such a cell.
	 */
	double eval_simplex(
	    index_t center_vertex_index,
	    const GEOGen::Vertex* V1,
	    const GEOGen::Vertex* V2,
	    const GEOGen::Vertex* V3,
	    index_t v_adj
	) const {
	    const double* p0 = OTM_->point_ptr(center_vertex_index);
	    const double* p1 = V1->point();
	    const double* p2 = V2->point();
	    const double* p3 = V3->point();
	    double m = GEO::Geom::tetra_signed_volume(p0, p1, p2, p3);
	    double fT = 0.0;
	    for(coord_index_t c = 0; c < 3; ++c) {
		double Uc = p1[c] - p0[c];
		double Vc = p2[c] - p0[c];
		double Wc = p3[c] - p0[c];
		fT +=
		    Uc * Uc +
		    Vc * Vc +
		    Wc * Wc +
		    Uc * Vc +
		    Vc * Wc +
		    Wc * Uc;
	    }
                
	    fT = m * (fT / 10.0 - w_[center_vertex_index]);
	    
	    double hij = 0.0;
	    if(Newton_step_ && v_adj != index_t(-1)) {
		const double* pj = OTM_->point_ptr(v_adj);
		hij =
		    GEO::Geom::triangle_area_3d(p1,p2,p3) /
		    (2.0 * GEO::Geom::distance(p0,pj,3)) ;
	    }

	    //  Spinlocks are used in multithreading mode, to avoid
	    // that two threads update g_[center_vertex_index]
	    // simultaneously.
	    if(spinlocks_ != nil) {
		spinlocks_->acquire_spinlock(center_vertex_index);
	    }
                
	    // +m because we maximize F <=> minimize -F
	    g_[center_vertex_index] += m;
	    if(Newton_step_) {
		// ... but here -m because Newton step =
		//  solve H p = -g    (minus g in the RHS).
		OTM_->add_i_right_hand_side(
		    center_vertex_index,-m
                );

		// -hij because we maximize F <=> minimize -F
		if(::fabs(hij) > 1e-8) {
		    // Diagonal is positive, extra-diagonal
		    // coefficients are negative,
		    // this is a convex function.
		    OTM_->add_ij_coefficient(
			center_vertex_index, v_adj, -hij
                    );
		    OTM_->add_ij_coefficient(
			center_vertex_index, center_vertex_index, hij
                    );
		}
	    }
	    
	    if(spinlocks_ != nil) {
		spinlocks_->release_spinlock(center_vertex_index);
	    }
	    
            // -fT because we maximize F <=> minimize -F
            return -fT;
	}

	/**
	 * \brief Evaluates the contribution of a simplex in weighted mode.
	 * \details The considered simplex connects four vertices of the current cell
	 *  (\p V0, \p V1, \p V2, \p V3).
	 * \param[in] center_vertex_index the index of the current Delaunay vertex
	 * \param[in] V0 , V1 , V2 , V3 four vertices of the current cell
	 * \param[in] v_adj the index of the Delaunay vertex of the cell on the other
	 *  side of (\p V1, \p V2, \p V3) or index_t(-1) if there is no such a cell.
	 */
	double eval_simplex_weighted(
	    index_t center_vertex_index,
	    const GEOGen::Vertex* V0,
	    const GEOGen::Vertex* V1,
	    const GEOGen::Vertex* V2,
	    const GEOGen::Vertex* V3,
	    index_t v_adj
	) const {
	    const double* p0 = V0->point();
	    const double* p1 = V1->point();
	    const double* p2 = V2->point();
	    const double* p3 = V3->point();

	    //  Note: when mg_ needs to be computed, weighted mode
	    // is used (and then weights are set to 1.0)
	    double p0_mass = weighted_ ? V0->weight() : 1.0;
	    double p1_mass = weighted_ ? V1->weight() : 1.0;
	    double p2_mass = weighted_ ? V2->weight() : 1.0;
	    double p3_mass = weighted_ ? V3->weight() : 1.0;
            
	    const double* q = OTM_->point_ptr(center_vertex_index);

            double Tvol = GEO::Geom::tetra_volume<3>(p0,p1,p2,p3);
            double Sp = p0_mass + p1_mass + p2_mass + p3_mass;
            
            double m = (Tvol * Sp) / 4.0;
            double rho[4], alpha[4];
            
            rho[0] = p0_mass;
            rho[1] = p1_mass;
            rho[2] = p2_mass;
            rho[3] = p3_mass;
            
            alpha[0] = Sp + rho[0];
            alpha[1] = Sp + rho[1];
            alpha[2] = Sp + rho[2];
            alpha[3] = Sp + rho[3];
            
            double dotprod_00 = 0.0;
            double dotprod_10 = 0.0;
            double dotprod_11 = 0.0;
            double dotprod_20 = 0.0;
            double dotprod_21 = 0.0;
            double dotprod_22 = 0.0;
            double dotprod_30 = 0.0;
            double dotprod_31 = 0.0;
            double dotprod_32 = 0.0;
            double dotprod_33 = 0.0;
            
            for(coord_index_t c = 0; c < 3; c++) {
                double sp0 = q[c] - p0[c];
                double sp1 = q[c] - p1[c];
                double sp2 = q[c] - p2[c];
                double sp3 = q[c] - p3[c];                
                dotprod_00 += sp0 * sp0;
                dotprod_10 += sp1 * sp0;
                dotprod_11 += sp1 * sp1;
                dotprod_20 += sp2 * sp0;
                dotprod_21 += sp2 * sp1;
                dotprod_22 += sp2 * sp2;
                dotprod_30 += sp3 * sp0;
                dotprod_31 += sp3 * sp1;
                dotprod_32 += sp3 * sp2;
                dotprod_33 += sp3 * sp3;                
            }
            
            double fT = 0.0;
            fT += (alpha[0] + rho[0]) * dotprod_00;  
            fT += (alpha[1] + rho[0]) * dotprod_10;  
            fT += (alpha[1] + rho[1]) * dotprod_11;  
            fT += (alpha[2] + rho[0]) * dotprod_20;  
            fT += (alpha[2] + rho[1]) * dotprod_21;  
            fT += (alpha[2] + rho[2]) * dotprod_22;  
            fT += (alpha[3] + rho[0]) * dotprod_30;  
            fT += (alpha[3] + rho[1]) * dotprod_31;  
            fT += (alpha[3] + rho[2]) * dotprod_32;  
            fT += (alpha[3] + rho[3]) * dotprod_33;  

            fT = Tvol * fT / 60.0 - m * w_[center_vertex_index];

            double hij = 0.0;
            if(Newton_step_ && v_adj != index_t(-1)) {
                const double* pi = OTM_->point_ptr(center_vertex_index);                
                const double* pj = OTM_->point_ptr(v_adj);
                hij = (p1_mass + p2_mass + p3_mass) *
                    GEO::Geom::triangle_area_3d(p1,p2,p3) /
                    (3.0 * 2.0 * GEO::Geom::distance(pi,pj,3)) ;
            }

            //  Spinlocks are used in multithreading mode, to avoid
            // that two threads update g_[center_vertex_index]
            // simultaneously.
            if(spinlocks_ != nil) {
                spinlocks_->acquire_spinlock(center_vertex_index);
            }
            
            // +m because we maximize F <=> minimize -F                
            g_[center_vertex_index] += m;
            if(Newton_step_) {
                // ... but here -m because Newton step =
                //  solve H p = -g    (minus g in the RHS).
                OTM_->add_i_right_hand_side(center_vertex_index,-m);
                // -hij because we maximize F <=> minimize -F
                if(::fabs(hij) > 1e-8) {
                    // Diagonal is positive, extra-diagonal
                    // coefficients are negative,
                    // this is a convex function.
                    OTM_->add_ij_coefficient(center_vertex_index, v_adj, -hij);
                    OTM_->add_ij_coefficient(
                        center_vertex_index, center_vertex_index, hij
                    );
                }
            }

            if(mg_ != nil) {
                for(index_t c=0; c<3; ++c) {
                    mg_[3*center_vertex_index+c] +=
                        (m / 4.0) * (p0[c] + p1[c] + p2[c] + p3[c]);
                }
            }
            
            if(spinlocks_ != nil) {
                spinlocks_->release_spinlock(center_vertex_index);
            }

            // -fT because we maximize F <=> minimize -F
            return -fT;
	}

	
    private:
	OptimalTransportMap3d* OTM_;
	bool weighted_;
	bool Newton_step_;
	vector<double> funcval_;
	index_t n_;
	const double* w_;
	double* g_;
	double* mg_;
    };


/**********************************************************************/

    /**
     * \brief Gets the number of connected components
     *  of each tetrahedra regions in a mesh.
     * \param[in] RVD a const reference to the mesh
     * \param[out] nb_cnx_comps on exit, nb_cnx_comps[r]
     *  contains the number of connected components of
     *  region r.
     */
    void get_nb_connected_components(
        const Mesh& RVD, vector<index_t>& nb_cnx_comps
    ) {
        Attribute<index_t> tet_region(RVD.cells.attributes(),"region");
        vector<bool> marked(RVD.cells.nb(), false);
        std::stack<index_t> S;
        for(index_t t = 0; t < RVD.cells.nb(); ++t) {
            if(!marked[t]) {
                index_t cur_v = index_t(tet_region[t]);
                marked[t] = true;
                S.push(t);
                while(!S.empty()) {
                    index_t cur_t = S.top();
                    S.pop();
                    for(index_t lf = 0; lf < 4; ++lf) {
                        index_t neigh = RVD.cells.tet_adjacent(cur_t, lf);
                        if(neigh != NO_CELL) {
                            if(
                                tet_region[neigh] == cur_v &&
                                !marked[neigh]
                            ) {
                                marked[neigh] = true;
                                S.push(neigh);
                            }
                        }
                    }
                }
                if(cur_v >= nb_cnx_comps.size()) {
                    nb_cnx_comps.resize(cur_v,0);
                }
                ++nb_cnx_comps[cur_v];
            }
        }
    }


    /**
     * \brief Computes the original and final vertices of
     *  the optimal transport.
     * \param[in] CVT the Centroidal Voronoi Tesselation that
     *  samples the target mesh M2
     * \param[in] OTM the Optimal Transport Map that back-projects
     *  the samples of the target mesh M2 onto the source mesh M1
     * \param[out] M1_vertices the source vertices, computed from
     *  the centroids of the power cells restricted to M1
     * \param[out] M2_vertices the destination vertices, copied
     *  from the sampling of M2
     * \param[out] nb_cnx_comps gives for each vertex the number of
     *  connected components of the power cell restricted to M1
     */
    void get_vertices(
        CentroidalVoronoiTesselation& CVT,
        OptimalTransportMap3d& OTM,
        vector<vec3>& M1_vertices,
        vector<vec3>& M2_vertices,
        vector<index_t>& nb_cnx_comps
    ) {
        index_t nb_vertices = OTM.RVD()->delaunay()->nb_vertices();
        M1_vertices.resize(nb_vertices);
        M2_vertices.resize(nb_vertices);
        nb_cnx_comps.assign(nb_vertices, 0);
        {
            for(index_t v = 0; v < nb_vertices; ++v) {
                const double* p = CVT.RVD()->delaunay()->vertex_ptr(v);
                M2_vertices[v] = vec3(p[0], p[1], p[2]);
            }
        }

        {
            vector<vec3> mg(nb_vertices, vec3(0.0, 0.0, 0.0));
            vector<double> m(nb_vertices, 0.0);

            Mesh RVD;
            Attribute<index_t> tet_region(RVD.cells.attributes(),"region");
            OTM.RVD()->compute_RVD(
                RVD,
                0,     // dim (0 means use default)
                false, // cells_borders_only
                true   // integration_simplices
            );
            RVD.vertices.set_dimension(3);
            RVD.cells.connect();

	    /*
            if(CmdLine::get_arg_bool("RVD")) {
                MeshIOFlags flags;
                flags.set_element(MESH_CELLS);
                flags.set_attribute(MESH_CELL_REGION);
                mesh_save(RVD, "RVD.meshb", flags);
            }
	    */

            for(index_t t = 0; t < RVD.cells.nb(); ++t) {
                index_t v =  tet_region[t];
                index_t v0 = RVD.cells.tet_vertex(t, 0);
                index_t v1 = RVD.cells.tet_vertex(t, 1);
                index_t v2 = RVD.cells.tet_vertex(t, 2);
                index_t v3 = RVD.cells.tet_vertex(t, 3);
                vec3 p0(RVD.vertices.point_ptr(v0));
                vec3 p1(RVD.vertices.point_ptr(v1));
                vec3 p2(RVD.vertices.point_ptr(v2));
                vec3 p3(RVD.vertices.point_ptr(v3));
                double mt = GEO::Geom::tetra_signed_volume(p0, p1, p2, p3);
                mg[v] += (mt / 4.0) * (p0 + p1 + p2 + p3);
                m[v] += mt;
            }
            for(index_t v = 0; v < nb_vertices; ++v) {
                double s = ::fabs(m[v]);
                if(s != 0.0) {
                    s = 1.0 / s;
                }
                M1_vertices[v] = s * mg[v];
            }
            get_nb_connected_components(RVD, nb_cnx_comps);
        }
    }

    /**
     * \brief Tests whether a tetrahedron is included inside
     *  a given tetrahedral mesh.
     * \details Implemented by sampling the tetrahedron and
     *  testing whether each sample is inside the mesh.
     * \param[in] AABB a const reference to a MeshCellsAABB
     * \param[in] p1 a const reference to the first vertex 
     *  of the tetrahedron
     * \param[in] p2 a const reference to the second vertex 
     *  of the tetrahedron
     * \param[in] p3 a const reference to the third vertex 
     *  of the tetrahedron
     * \param[in] p4 a const reference to the fourth vertex 
     *  of the tetrahedron
     */
    bool mesh_contains_tet(
        const MeshCellsAABB& AABB,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4
    ) {
        const index_t NB = 5;
        for(index_t u1=0; u1<=NB; ++u1) {
            double s1 = double(u1)/double(NB);            
            for(index_t u2=0; u1+u2<=NB; ++u2) {
                double s2 = double(u2)/double(NB);                
                for(index_t u3=0; u1+u2+u3<=NB; ++u3) {
                    double s3 = double(u3)/double(NB);                    
                    index_t u4=NB-u1-u2-u3;
                    double s4 = double(u4)/double(NB);

                    // Skip the four vertices of the tetrahedron
                    if(
                        (u1 == NB) || (u2 == NB) ||
                        (u3 == NB) || (u4 == NB)
                    ) {
                        continue;
                    }
                    
                    vec3 g = s1*p1+s2*p2+s3*p3+s4*p4;
                    if(AABB.containing_tet(g) == MeshCellsAABB::NO_TET) {
                        return false;
                    }
                }
            }
        }
        return true;
    }
}


namespace GEO {

    OptimalTransportMap3d* OptimalTransportMap3d::instance_ = nil;
    
    OptimalTransportMap3d::OptimalTransportMap3d(
        Mesh* mesh, const std::string& delaunay_in, bool BRIO
    ) : mesh_(mesh) {

	geo_cite("DBLP:conf/compgeom/AurenhammerHA92");
	geo_cite("DBLP:journals/cgf/Merigot11");
	geo_cite("journals/M2AN/LevyNAL15");

	std::string delaunay = delaunay_in;
	if(delaunay == "default") {
	    delaunay = "PDEL";
	}
	
        epsilon_regularization_ = 0.0;

        // Mesh is supposed to be embedded in 4 dim (with
        // 4th dimension set to zero).
        geo_assert(mesh->vertices.dimension() == 4);

        // Note: we represent power diagrams as 4d Voronoi diagrams
        delaunay_ = Delaunay::create(4, delaunay);

        RVD_ = RestrictedVoronoiDiagram::create(delaunay_, mesh_);
        RVD_->set_volumetric(true);
        RVD_->set_check_SR(true);
        RVD_->create_threads();
        
        //   No need to reorder vertices if BRIO is activated since
        // vertices are then already reordered.
        if(BRIO) {
            RVD_->delaunay()->set_reorder(false);
        }

        newton_ = false;
        
        instance_ = nil;
        lambda_p_ = 0.0;
        total_mass_ = 0.0;
        current_call_iter_ = 0;
        epsilon_ = 0.01;
        level_ = 0;
        
        save_RVD_iter_ = false;
        save_RVD_last_iter_ = false;
        show_RVD_seed_ = false;
        current_iter_ = 0;
        
        pretty_log_ = false; // CmdLine::get_arg_bool("log:pretty");

        w_did_not_change_ = false;

        measure_of_smallest_cell_ = 0.0;

	polyhedron_callback_ = new OTMPolyhedronCallback(this);

	Laguerre_centroids_ = nil;
    }

    OptimalTransportMap3d::~OptimalTransportMap3d() {
	delete polyhedron_callback_;
    }
    
    void OptimalTransportMap3d::set_points(
        index_t nb_points, const double* points
    ) {
        // Note: we represent power diagrams as 4d Voronoi diagrams.
        // The target points are lifted to 4d.
        points_4d_.resize(nb_points * 4);
        for(index_t i = 0; i < nb_points; ++i) {
            points_4d_[i * 4] = points[i * 3];
            points_4d_[i * 4 + 1] = points[i * 3 + 1];
            points_4d_[i * 4 + 2] = points[i * 3 + 2];
            points_4d_[i * 4 + 3] = 0.0;
        }
        weights_.assign(nb_points, 0);
        total_mass_ = 0.0;
        
        //   This is terribly confusing, the parameters for
        // a power diagram are called "weights", and the
        // standard attribute name for vertices density is
        // also called "weight" (and is unrelated).
        //   In this program, what is called weight corresponds
        // to the parameters of the power diagram (except the
        // name of the attribute), and everything that corresponds
        // to mass/density is called mass.
        Attribute<double> vertex_mass;
        vertex_mass.bind_if_is_defined(
            mesh_->vertices.attributes(), "weight"
        );
        
        for(index_t t = 0; t < mesh_->cells.nb(); ++t) {
            double tet_mass = GEO::Geom::tetra_volume<3>(
                mesh_->vertices.point_ptr(mesh_->cells.tet_vertex(t, 0)),
                mesh_->vertices.point_ptr(mesh_->cells.tet_vertex(t, 1)),
                mesh_->vertices.point_ptr(mesh_->cells.tet_vertex(t, 2)),
                mesh_->vertices.point_ptr(mesh_->cells.tet_vertex(t, 3))
            );
            if(vertex_mass.is_bound()) {
                tet_mass *= (
                    vertex_mass[mesh_->cells.tet_vertex(t, 0)] +
                    vertex_mass[mesh_->cells.tet_vertex(t, 1)] +
                    vertex_mass[mesh_->cells.tet_vertex(t, 2)] +
                    vertex_mass[mesh_->cells.tet_vertex(t, 3)]
                ) / 4.0;
            }
            total_mass_ += tet_mass;
        }
        lambda_p_ = total_mass_ / double(nb_points);
    }


    // See http://arxiv.org/abs/1603.05579
    // Kitawaga, Merigot, Thibert,
    // A Newton Algorithm for semi-discrete OT
    
    void OptimalTransportMap3d::optimize_full_Newton(
        index_t max_iterations, index_t n
    ) {
        if(n == 0) {
            n = index_t(points_4d_.size() / 4);
        }

        vector<double> pk(n);
        vector<double> xk(n);
        vector<double> gk(n);
        double fk=0.0;
        bool converged = false;

        double epsilon0 = 0.0;
        w_did_not_change_ = false;
        
        for(index_t k=0; k<max_iterations; ++k) {
            std::cerr << "======= k = " << k << std::endl;
            xk=weights_;

            new_linear_system(n);
            eval_func_grad_Hessian(n,xk.data(),fk,gk.data());
            
            if(k == 0) {
                newiteration();
            }
            
            if(epsilon0 == 0.0) {
                epsilon0 = 0.5 * geo_min(measure_of_smallest_cell_, lambda_p_);
            }
            
            Logger::out("OTM") << "   Solving linear system" << std::endl;

            solve_linear_system(pk.data());

            std::cerr << "Line search ..." << std::endl;

            w_did_not_change_ = false;
            
            double alphak = 1.0;
            double gknorm = g_norm_;
            
            for(index_t inner_iter=0; inner_iter < 100; ++inner_iter) {
                std::cerr << "      inner iter = " << inner_iter << std::endl;

                // weights = xk + alphak pk
                for(index_t i=0; i<n; ++i) {
                    weights_[i] = xk[i] + alphak * pk[i];
                }

		OTMPolyhedronCallback* cb_otm =
		    dynamic_cast<OTMPolyhedronCallback*>(polyhedron_callback_);
		
                // Compute cell measures and nbZ.
		if(Laguerre_centroids_ != nil) {
		    cb_otm->set_Laguerre_centroids(Laguerre_centroids_);
		}
		
		funcgrad(n,weights_.data(),fk,gk.data());

		if(Laguerre_centroids_ != nil) {
		    cb_otm->set_Laguerre_centroids(nil);
		}
		
                std::cerr << "cell measure :"
                          << measure_of_smallest_cell_
                          << "(>=?)" << epsilon0 << std::endl;
                std::cerr << "gradient norm:"
                          << g_norm_ << "(<=?)"
                          << (1.0 - 0.5*alphak) * gknorm << std::endl;
                if(
                    (measure_of_smallest_cell_ >= epsilon0) &&
                    (g_norm_ <= (1.0 - 0.5*alphak) * gknorm) 
                ) {
                    if(g_norm_ < gradient_threshold(n)) {
                        converged = true;
                    }
                    break;
                } 
                // Else we halve the step.
                alphak /= 2.0;
            }
            newiteration();
            if(converged) {
                break;
            }
            // No need to update the power diagram at next iteration,
            // since we will evaluate the Hessian for the same weight
            // vector.
            w_did_not_change_ = true;
        }
        if(save_RVD_last_iter_) {
            save_RVD(current_iter_);
        }
    }
    
    void OptimalTransportMap3d::optimize(index_t max_iterations) {

        level_ = 0;
        
        if(newton_) {
            optimize_full_Newton(max_iterations);
            return;
        }
        
        index_t n = index_t(points_4d_.size() / 4);
        index_t m = 7;
        
        Optimizer_var optimizer = Optimizer::create("HLBFGS");
        
        optimizer->set_epsg(gradient_threshold(n));
        optimizer->set_epsf(0.0);
        optimizer->set_epsx(0.0);
        
        optimizer->set_newiteration_callback(newiteration_CB);
        optimizer->set_funcgrad_callback(funcgrad_CB);

        optimizer->set_N(n);
        optimizer->set_M(m);
        optimizer->set_max_iter(max_iterations);
        instance_ = this;
        current_call_iter_ = 0;
        current_iter_ = 0;
        optimizer->optimize(weights_.data());
        instance_ = nil;
        // To make sure everything is reset properly
        double dummy = 0;
        funcgrad(n, weights_.data(), dummy, nil);
        Logger::out("OTM")
            << "Used " << current_call_iter_ << " iterations" << std::endl;
        if(save_RVD_last_iter_) {
            save_RVD(current_iter_);
        }
    }


    void OptimalTransportMap3d::optimize_level(
        index_t b, index_t e, index_t max_iterations
    ) {

        // If this is not the first level, propagate the weights from
        // the lower levels.
        if(b != 0) {
            
            //   Create a nearest neighbor search data structure
            // and insert the [0..b) samples into it (they were
            // initialized at previous calls).
            NearestNeighborSearch_var NN = NearestNeighborSearch::create(3);
            NN->set_points(b, points_4d_.data(), 4);
            index_t degree = 2; // CmdLine::get_arg_uint("fitting_degree");
            
            // If degree \notin {1,2}, use weight of nearest sample
            if(degree < 1 || degree > 2) {
                for(index_t i = b; i < e; ++i) {
                    weights_[i] =
                        weights_[
                            NN->get_nearest_neighbor(&points_4d_[4 * i])
                        ];
                }
            } else {
                
                //   If degree \in {1,2} use linear least squares to
                // compute an estimate of the weight function.
                LinearLeastSquares LLS(degree);
                for(index_t i = b; i < e; ++i) {
                    const index_t nb = 10 * degree;
                    index_t neighbor[100];
                    double dist[100];
                    NN->get_nearest_neighbors(
                        nb, &points_4d_[4 * i], neighbor, dist
                    );
                    LLS.begin();
                    for(index_t jj = 0; jj < nb; ++jj) {
                        if(dist[jj] != 0.0) {
                            index_t j = neighbor[jj];
                            LLS.add_point(&points_4d_[4 * j], weights_[j]);
                        }
                    }
                    LLS.end();
                    weights_[i] = LLS.eval(&points_4d_[4 * i]);
                }
            }
        }
        
        // Optimize the weights associated with the sequence [0,e)
        index_t n = e;
        
        // Important! lambda_p_ (target measure of a cell) needs
        // to be updated, since it depends on the number of samples
        // (that varies at each level).
        lambda_p_ = total_mass_ / double(n);

        if(newton_) {
            optimize_full_Newton(max_iterations, n);
            return;
        }
        
        index_t m = 7;
        Optimizer_var optimizer = Optimizer::create("HLBFGS");

        optimizer->set_epsg(gradient_threshold(n));
        optimizer->set_epsf(0.0);
        optimizer->set_epsx(0.0);
        optimizer->set_newiteration_callback(newiteration_CB);
        optimizer->set_funcgrad_callback(funcgrad_CB);

        optimizer->set_N(n);
        optimizer->set_M(m);
        optimizer->set_max_iter(max_iterations);
        instance_ = this;
        current_call_iter_ = 0;
        optimizer->optimize(weights_.data());
        instance_ = nil;
        
        // To make sure everything is reset properly
        double dummy = 0;
        funcgrad(n, weights_.data(), dummy, nil);
    }

    void OptimalTransportMap3d::optimize_levels(
        const vector<index_t>& levels, index_t max_iterations
    ) {
        if(levels.size() > 2) {
            Logger::out("OTM") << "Using " << levels.size()-1
                               << " levels" << std::endl;
        } else {
            Logger::out("OTM") << "Using 1 level" << std::endl;
        }
        for(index_t l = 0; l + 1 < levels.size(); ++l) {
            level_ = l+1;
            index_t b = levels[l];
            index_t e = levels[l + 1];
            vector<index_t> brio_levels;
            for(index_t i=0; i<=l+1; ++i) {
                brio_levels.push_back(levels[i]);
            }
            RVD_->delaunay()->set_BRIO_levels(brio_levels);
            optimize_level(b, e, max_iterations);
        }
        if(save_RVD_last_iter_) {
            save_RVD(current_iter_);
        }
    }
    
    void OptimalTransportMap3d::funcgrad_CB(
        index_t n, double* x, double& f, double* g
    ) {
        instance_->funcgrad(n, x, f, g);
    }

    void OptimalTransportMap3d::newiteration_CB(
        index_t n, const double* x, double f, const double* g, double gnorm
    ) {
        geo_argused(n);
        geo_argused(x);
        geo_argused(f);
        geo_argused(g);
        geo_argused(gnorm);
        instance_->newiteration();
    }

    void OptimalTransportMap3d::newiteration() {
        //xxx std::cerr << "newiteration" << std::endl;
        if(save_RVD_iter_) {
            std::cerr << "  save iter" << std::endl;
            save_RVD(current_iter_);
        }
        ++current_iter_;
    }

    static void cell_shrink_to_animation(Mesh& mesh) {
        Attribute<index_t> cell_region;
        Attribute<double> W(mesh.vertices.attributes(),"w");
        cell_region.bind_if_is_defined(
            mesh.cells.attributes(), "region"
        );
        if(!cell_region.is_bound()) {
            Logger::err("OTM") << "region: no such cell attribute"
                               << std::endl;
            return;
        }
        index_t nb_RVD_cells = 0;
        for(index_t c=0; c<mesh.cells.nb(); ++c) {
            nb_RVD_cells = geo_max(nb_RVD_cells, cell_region[c]);
        }
        ++nb_RVD_cells;
        vector<vec3> center(nb_RVD_cells, vec3(0.0, 0.0, 0.0));
        vector<index_t> nb(nb_RVD_cells, 0);
        for(index_t c=0; c<mesh.cells.nb(); ++c) {
            index_t rvc = cell_region[c];
            for(index_t lv=0; lv<4; ++lv) {
                index_t v = mesh.cells.vertex(c,lv);
                center[rvc] += Geom::mesh_vertex(mesh,v);
                ++nb[rvc];
            }
        }
        for(index_t rvc=0; rvc<nb_RVD_cells; ++rvc) {
            if(nb[rvc] != 0) {
                center[rvc] /= double(nb[rvc]);
            }
        }
        for(index_t v=0; v<mesh.vertices.nb(); ++v) {
            W[v] = mesh.vertices.point_ptr(v)[3];
        }
        mesh.vertices.set_dimension(6);
        for(index_t c=0; c<mesh.cells.nb(); ++c) {
            index_t rvc = cell_region[c];
            for(index_t lv=0; lv<4; ++lv) {
                index_t v = mesh.cells.vertex(c,lv);
                const vec3& g = center[rvc];
                double* p = mesh.vertices.point_ptr(v);
                p[3] = g[0];
                p[4] = g[1];
                p[5] = g[2];
            }
        }
    }
    
    void OptimalTransportMap3d::get_RVD(Mesh& RVD_mesh) {
        RVD_mesh.clear();
        Attribute<index_t> tet_region(RVD_mesh.cells.attributes(),"region");
        RVD()->compute_RVD(
            RVD_mesh,
            0,             // dim (0 means use default)
            false,         // borders_only
            show_RVD_seed_ // integration_simplices
        );
        if(!show_RVD_seed_) {
            cell_shrink_to_animation(RVD_mesh);
        }
    }

    void OptimalTransportMap3d::compute_Laguerre_centroids(double* centroids) {
        vector<double> g(nb_points(), 0.0);
        Memory::clear(centroids, nb_points()*sizeof(double)*3);

	OTMPolyhedronCallback* cb_otm =
	    dynamic_cast<OTMPolyhedronCallback*>(polyhedron_callback_);
	
	cb_otm->set_Laguerre_centroids(centroids);
	cb_otm->set_g(g.data());
	{
	    Stopwatch* W = nil;
	    if(newton_) {
		W = new Stopwatch("RVD");
		Logger::out("OTM") << "In RVD (centroids)..." << std::endl;
	    }
	    RVD_->for_each_polyhedron(*cb_otm,false,false,true);
	    if(newton_) {
		delete W;
	    }
	}
	cb_otm->set_Laguerre_centroids(nil);	    
	
        for(index_t v=0; v<nb_points(); ++v) {
            centroids[3*v  ] /= g[v];
            centroids[3*v+1] /= g[v];
            centroids[3*v+2] /= g[v];            
        }
    }
    
    void OptimalTransportMap3d::save_RVD(index_t id) {
	/*
        if(scene_graph_ != nil) {
            MeshGrob* RVD_mesh =
                MeshGrob::find_or_create(
                    scene_graph_, "RVD_" + String::to_string(id)
                );
            get_RVD(*RVD_mesh);
            RVD_mesh->update();
        } else */ {
            Mesh RVD_mesh;
            get_RVD(RVD_mesh);
            MeshIOFlags flags;
            flags.set_attribute(MESH_CELL_REGION);
            flags.set_attribute(MESH_FACET_REGION);            
            mesh_save(
                RVD_mesh,
                "RVD_" + String::to_string(id) + ".geogram",
                flags
            );
        }
    }

    void OptimalTransportMap3d::funcgrad(
        index_t n, double* w, double& f, double* g
    ) {

	OTMPolyhedronCallback* cb_otm =
	    dynamic_cast<OTMPolyhedronCallback*>(polyhedron_callback_);
	
        bool is_Newton_step = cb_otm->is_Newton_step();	    
        
        // For now, always compute function and gradient
        bool update_fg = true;
        
        // Delaunay triangulation is only updated if function and
        // gradient is evaluated. If only Hessian needs to be evaluated,
        // then it is at the same point as the latest function and gradient
        // evaluation (see Yang Liu's CVT-Newton code).
        if(update_fg && !w_did_not_change_) {
            // Step 1: determine the 4d embedding from the weights
            double W = 0.0;
            for(index_t p = 0; p < n; ++p) {
                W = geo_max(W, w[p]);
            }
            for(index_t p = 0; p < n; ++p) {
                points_4d_[4 * p + 3] = ::sqrt(W - w[p]);
            }
        
            // Step 2: compute function and gradient
            {
		Stopwatch* W = nil;
		if(newton_) {
		    W = new Stopwatch("Power diagram");
		    Logger::out("OTM") << "In power diagram..." << std::endl;
		}
                delaunay_->set_vertices(n, points_4d_.data());
		if(newton_) {
		    delete W;
		}
            }
        }

        if(is_Newton_step) {
            update_sparsity_pattern();
        }
        
        if(g == nil) {
            if(pretty_log_) {
                CmdLine::ui_clear_line();
                CmdLine::ui_message(last_stats_ + "\n");
            }
            return;
        }
        
        if(update_fg) {
            f = 0.0;
            for(index_t p = 0; p < n; ++p) {
                g[p] = 0.0;
            }
        }

	cb_otm->set_w(w,n);
	cb_otm->set_g(g);
	cb_otm->set_nb_threads(Process::maximum_concurrent_threads());

	if(cb_otm->has_Laguerre_centroids()) {
	    Memory::clear(cb_otm->Laguerre_centroids(), nb_points()*sizeof(double)*3);	    
	}
	
	{
	    Stopwatch* W = nil;
	    if(newton_) {
		W = new Stopwatch("RVD");
		Logger::out("OTM") << "In RVD (funcgrad)..." << std::endl;
	    }
	    RVD_->for_each_polyhedron(*cb_otm,false,false,true);
	    if(newton_) {
		delete W;
	    }
	}
	f = cb_otm->funcval(); 

	if(cb_otm->has_Laguerre_centroids()) {
	    for(index_t v=0; v<nb_points(); ++v) {
		cb_otm->Laguerre_centroids()[3*v  ] /= g[v];
		cb_otm->Laguerre_centroids()[3*v+1] /= g[v];
		cb_otm->Laguerre_centroids()[3*v+2] /= g[v];            
	    }
	}
	
        if(update_fg) {
            measure_of_smallest_cell_ = Numeric::max_float64();
            for(index_t i=0; i<n; ++i) {
                measure_of_smallest_cell_ =
                    geo_min(measure_of_smallest_cell_, g[i]);
            }
        }
        
        if(update_fg) {
            for(index_t p = 0; p < n; ++p) {
                f += lambda_p_ * w[p];
                // Note: we minimize -f instead of maximizing f,
                // therefore, in the paper:
                //    g[p] = lambda_p - mesure(power cell associated with p)
                //
                // What is programmed:
                //    g[p] = mesure(power cell associated with p) - lambda_p
                g[p] -= lambda_p_;

                if(is_Newton_step) {
                    // Newton step: solve H deltax = -g
                    // (note the minus sign on the right hand side)
                    // g[p] -= lamnda_p_ -> RHS[p] += lambda_p_
                    add_i_right_hand_side(p, lambda_p_);
                }
            }
        }

        double max_diff = 0.0;
        double avg_diff = 0.0;
        index_t nb_empty_cells = 0;
        for(index_t p = 0; p < n; ++p) {
            double cur_diff = ::fabs(g[p]);
            max_diff = geo_max(max_diff, cur_diff);
            avg_diff += cur_diff / double(n);
            // At this step, g[p] = mu(Lag(p)) - lambda_p
            // We add lambda_p to retreive mu(Lag(p)) and to
            // count empty cells if mu(Lap(p)) is smaller than
            // a threshold.
            if(::fabs(g[p] + lambda_p_) < 1e-10) {
                nb_empty_cells++;
            }
        }

        
        // Regularisation: minimize the squared norm of the weight
        // vector to remove a translational degree of freedom.
        // It seems to make the overall convergence slower
        // (but this may be due to a wrong
        // scaling between the different levels, to be investigated...)
        if(epsilon_regularization_ != 0.0) {
            if(update_fg) {                                    
                for(index_t p = 0; p < n; ++p) {
                    f += 0.5 * epsilon_regularization_ * lambda_p_ * w[p]*w[p];
                    g[p] += epsilon_regularization_ * lambda_p_ * w[p];
                }
            }
            if(is_Newton_step) {
                for(index_t p = 0; p < n; ++p) {
                    add_ij_coefficient(
                        p,p,epsilon_regularization_*lambda_p_
                    );
                    add_i_right_hand_side(
                        p,-epsilon_regularization_*lambda_p_*w[p]
                    );
                }
            }
        }

        
        double gNorm = 0.0;
        for(index_t i = 0; i < n; ++i) {
            gNorm += geo_sqr(g[i]);
        }
        gNorm = ::sqrt(gNorm);

        nbZ_ = nb_empty_cells;
        
        std::ostringstream str;
        if(pretty_log_) {
            if(level_ == 0) {
                str << "o-[OTM         ] " ;
            } else {
                str << "o-[OTM Lvl." << level_ << "   ] " ;
            }
        } else {
            if(level_ == 0) {
                str << "   OTM      : " ;
            } else {
                str << "   OTM Lvl." << level_ << ": " ;
            }
        }
        
        str << "iter=" << current_call_iter_
            << " nbZ=" << nb_empty_cells
            //                << " f=" << f
            //                << " avg_diff=" << avg_diff
            //                << " max_diff=" << max_diff
            << " g=" << gNorm
            << " f=" << f 
            << " threshold=" << gradient_threshold(n);
        last_stats_ = str.str();

        g_norm_ = gNorm;
        
        // "custom task progress" (clears the standard message
        // and replaces it with another one).
        if(pretty_log_) {
            if(current_call_iter_ != 0) {
                CmdLine::ui_clear_line();
            }
            CmdLine::ui_message(str.str());
        } else {
            str << " f=" << f;
            CmdLine::ui_message(str.str() + "\n");
        }
        ++current_call_iter_;
    }

    void OptimalTransportMap3d::eval_func_grad_Hessian(
        index_t n, const double* w, double& f, double* g
    ) {
	OTMPolyhedronCallback* cb_otm =
	    dynamic_cast<OTMPolyhedronCallback*>(polyhedron_callback_);
	cb_otm->set_Newton_step(true);
	funcgrad(n,(double*)w,f,g);
	cb_otm->set_Newton_step(false);
    }
    
/************************************************************/

    void compute_morph(
        CentroidalVoronoiTesselation& CVT,
        OptimalTransportMap3d& OTM,
        Mesh& morph,
        bool filter_tets
    ) {
        geo_assert(CVT.volumetric());
        
        Logger::out("OTM")
            << "Computing coherent tet mesh"
            << std::endl;

        morph.clear();
        morph.vertices.set_dimension(6);
        
        // Step 1: Compute the candidate tets from the Delaunay triangulation
        // of the samples restricted to M2.
        vector<index_t> morph_tets;
        morph_tets.clear();
        {
            vector<double> embedding;
            CVT.RVD()->compute_RDT(
                morph_tets,
                embedding,
                RestrictedVoronoiDiagram::RDT_SEEDS_ALWAYS
            );
        }

        // Step 2: Compute the original vertices location (centroids
        //  of power cells restricted to M1) and get the final vertices
        //  locations (sampling of M2).
        vector<vec3> M1_vertices;
        vector<vec3> M2_vertices;
        vector<index_t> nb_cnx_comps;
        get_vertices(CVT, OTM, M1_vertices, M2_vertices, nb_cnx_comps);
        index_t nb_vertices = M1_vertices.size();

        //  Gather all the points coordinates in a vector of
        // 6d points.
        vector<double> morph_vertices(nb_vertices*6);
        for(index_t v=0; v<nb_vertices; ++v) {
            morph_vertices[v*6  ] = M2_vertices[v].x;
            morph_vertices[v*6+1] = M2_vertices[v].y;
            morph_vertices[v*6+2] = M2_vertices[v].z;
            morph_vertices[v*6+3] = M1_vertices[v].x;
            morph_vertices[v*6+4] = M1_vertices[v].y;
            morph_vertices[v*6+5] = M1_vertices[v].z;            
        }

        // Step 3: Filter-out the tets incident to a vertex
        // that splits during transport.
        index_t nb_tets = morph_tets.size()/4;
        vector<bool> tet_to_remove(nb_tets, false);
        for(index_t t=0; t<nb_tets; ++t) {
            for(index_t lv=0; lv<4; ++lv) {
                index_t v = morph_tets[4*t+lv];
                if(nb_cnx_comps[v] > 1) {
                    tet_to_remove[t] = true;
                    break;
                }
            }
        }

        // Step 4: Filter-out the tets that are not contained by
        // the initial mesh M1.
        if(filter_tets) {
            MeshCellsAABB AABB(*OTM.RVD()->mesh());
            try {
                ProgressTask progress("Classifying", 100);        
                for(index_t t=0; t<nb_tets; ++t) {
                    progress.progress(t * 100 / nb_tets);
                    if(!tet_to_remove[t]) {
                        vec3 p[4];
                        for(index_t lv=0; lv<4; ++lv) {
                            for(coord_index_t c=0; c<3; ++c) {
                                index_t v = morph_tets[4*t+lv];
                                p[lv][c] = morph_vertices[v*6+3+c];
                            }
                        }
                        if(!mesh_contains_tet(AABB, p[0], p[1], p[2], p[3])) {
                            tet_to_remove[t] = true;
                        }
                    }
                }
            } catch(...) {
            }
        }

        // Step 5: create the output mesh.
        vector<index_t> filtered_morph_tets;
        for(index_t t=0; t<nb_tets; ++t) {
            if(!tet_to_remove[t]) {
                filtered_morph_tets.push_back(morph_tets[4*t]);
                filtered_morph_tets.push_back(morph_tets[4*t+1]);
                filtered_morph_tets.push_back(morph_tets[4*t+2]);
                filtered_morph_tets.push_back(morph_tets[4*t+3]);
            }
        }
        
        morph.cells.assign_tet_mesh(
            6, morph_vertices, filtered_morph_tets, true
        );
        morph.cells.connect();
        morph.cells.compute_borders();
    }

    void compute_singular_surface(        
        CentroidalVoronoiTesselation& CVT,
        OptimalTransportMap3d& OTM,
        Mesh& singular
    ) {

        std::set<bindex> edges;
        {
            vector<index_t> simplices;
            vector<double> embedding;
            CVT.RVD()->compute_RDT(
                simplices, 
                embedding, 
                RestrictedVoronoiDiagram::RDT_SEEDS_ALWAYS
            );
            for(index_t t=0; t*4<simplices.size(); ++t) {
                index_t v1 = simplices[t*4];
                index_t v2 = simplices[t*4+1];
                index_t v3 = simplices[t*4+2];
                index_t v4 = simplices[t*4+3];
                edges.insert(bindex(v1,v2));
                edges.insert(bindex(v1,v3));
                edges.insert(bindex(v1,v4));
                edges.insert(bindex(v2,v3));
                edges.insert(bindex(v2,v4));
                edges.insert(bindex(v3,v4));
            }
        }

        Mesh RVD;
        MeshIOFlags flags;
        flags.set_element(MESH_CELLS);
        flags.set_attribute(MESH_CELL_REGION);
        OTM.RVD()->compute_RVD(
            RVD,
            0,
            false, // cells_borders_only
            true   // integration_simplices
        );
        RVD.vertices.set_dimension(3);
        RVD.cells.connect();

        singular.clear();
        singular.vertices.set_dimension(3);
        
        vector<index_t> triangles;

        Attribute<index_t> tet_region(RVD.cells.attributes(),"region");
        
        for(index_t t=0; t<RVD.cells.nb(); ++t) {
            index_t v1 = tet_region[t];
            for(index_t f=0; f<4; ++f) {
                index_t nt = RVD.cells.tet_adjacent(t,f);
                if(nt != NO_CELL) {
                    index_t v2 = tet_region[nt];
                    if(v1 != v2 && edges.find(bindex(v1,v2)) == edges.end()) {
                        for(index_t i=0; i<3; ++i) {
                            index_t lv =
                                RVD.cells.local_tet_facet_vertex_index(f,i);
                            index_t v = RVD.cells.tet_vertex(t,lv);
                            triangles.push_back(v);
                        }
                    }
                }
            }
        }

        singular.vertices.assign_points(
            RVD.vertices.point_ptr(0),
            RVD.vertices.dimension(), RVD.vertices.nb()
        );
        singular.facets.assign_triangle_mesh(
            triangles, true 
        );
    }


    /*************************************************************************/

    // TODO: use OpenNL buffers to avoid data copy and allocation.
    // TODO: in the Euler code, see if we do not have duplicated computations, i.e.
    //    - Power diagrams when leaving and entering iteration ?
    //    - Centroids: do we restart a RVD computation ?
    
    void OptimalTransportMap3d::update_sparsity_pattern() {
        // Does nothing for now,
	// (we let OpenNL discover the sparsity pattern)
	// Tryed smarter things, but was not faster...
    }

    void OptimalTransportMap3d::new_linear_system(index_t n) {
        nlNewContext();
            
        bool use_SUPERLU = false;
        
        if(use_SUPERLU) {
            nlInitExtension("SUPERLU");
        }

	nlEnable(NL_VERBOSE);
	
        nlSolverParameteri(NL_NB_VARIABLES, NLint(n));
        if(use_SUPERLU) {
            nlSolverParameteri(NL_SOLVER, NL_PERM_SUPERLU_EXT);
        } else {
            nlSolverParameteri(NL_SOLVER, NL_CG);
            nlSolverParameteri(NL_PRECONDITIONER, NL_PRECOND_JACOBI);
            nlSolverParameteri(NL_SYMMETRIC, NL_TRUE);
            nlSolverParameterd(NL_THRESHOLD, 0.001);
            nlSolverParameteri(NL_MAX_ITERATIONS, 1000);                
        }
        nlBegin(NL_SYSTEM);
        nlBegin(NL_MATRIX);
    }

    void OptimalTransportMap3d::solve_linear_system(double* x) {
        nlEnd(NL_MATRIX);
        nlEnd(NL_SYSTEM);
        nlSolve();
        // Query and display OpenNL stats.
        NLint n;
        nlGetIntegerv(NL_NB_VARIABLES, &n);
        {
            int used_iters;
            double elapsed_time;
            double gflops;
            double error;
            nlGetIntegerv(NL_USED_ITERATIONS, &used_iters);
            nlGetDoublev(NL_ELAPSED_TIME, &elapsed_time);
            nlGetDoublev(NL_GFLOPS, &gflops);
            nlGetDoublev(NL_ERROR, &error);                
            std::cerr << "   "
                      << used_iters << " iters in "
                      << elapsed_time << " seconds "
                      << gflops << " GFlop/s"
                      << "  ||Ax-b||/||b||="
                      << error
                      << std::endl;
        }
        for(NLint i=0; i<n; ++i) {
            x[i] = nlGetVariable(NLuint(i));
        }
        nlDeleteContext(nlGetCurrent());
    }
    
    /**********************************************************************/
    
    void compute_Laguerre_centroids_3d(
        Mesh* omega,
        index_t nb_points,
        const double* points,
        double* centroids,
	bool parallel_pow
    ) {
        omega->vertices.set_dimension(4);

        // false = no BRIO
        // (OTM does not use multilevel and lets Delaunay
        //  reorder the vertices)
        OptimalTransportMap3d OTM(
	    omega,
	    std::string(parallel_pow ? "PDEL" : "BPOW"),
	    false
	);

        OTM.set_regularization(1e-3);
        OTM.set_Newton(true);
        OTM.set_points(nb_points, points);
        OTM.set_epsilon(0.01);
	OTM.set_Laguerre_centroids(centroids);
        OTM.optimize(1000);
        omega->vertices.set_dimension(3);        
    }

    
}
