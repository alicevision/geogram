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

#include <exploragram/optimal_transport/optimal_transport_2d.h>
#include <geogram/voronoi/generic_RVD_vertex.h>
#include <geogram/voronoi/generic_RVD_polygon.h>
#include <geogram/voronoi/RVD_callback.h>
#include <geogram/basic/geometry.h>
#include <geogram/basic/geometry_nd.h>
#include <geogram/basic/stopwatch.h>

// TODO: rewrite OTMPolygonCallback and SurfaceOTMPolygonCallback
// as OTMPolyhedronCallback in optimal_transport_3d.cpp (much clearer...)

namespace {
    using namespace GEO;

    /**
     * \brief A RVDPolygonCallback that stores the Restricted Voronoi
     *  Diagram in a Mesh.
     */
    class ComputeRVDPolygonCallback : public RVDPolygonCallback {
    public:
	ComputeRVDPolygonCallback(OptimalTransportMap* OTM, Mesh* target) :
	    OTM_(OTM), target_(target) {
	    target_->clear();
	    target_->vertices.set_dimension(3);
	    chart_.bind(target_->facets.attributes(), "chart");
	}

	~ComputeRVDPolygonCallback() {
	    chart_.unbind();
	}

	virtual void operator() (
	    index_t v,
	    index_t t,
	    const GEOGen::Polygon& P
	) const {
	    geo_argused(v);
	    geo_argused(t);
	    index_t voffset = target_->vertices.nb();
	    FOR(i,P.nb_vertices()) {
		const double* p = P.vertex(i).point();
		if(OTM_->dimension() == 2) {
		    target_->vertices.create_vertex(
			vec3(p[0], p[1], 0.0).data()
		    );
		} else {
		    target_->vertices.create_vertex(
			vec3(p[0], p[1], p[2]).data()
		    );
		}
	    }
	    index_t f = target_->facets.create_polygon(P.nb_vertices());
	    FOR(i,P.nb_vertices()) {
		target_->facets.set_vertex(f,i,voffset+i);
	    }
	    const_cast<Attribute<index_t>&>(chart_)[f] = v;
	}
	
    private:
	OptimalTransportMap* OTM_;
	Mesh* target_;
	Attribute<index_t> chart_;
    };
    
    /**
     * \brief Computes the contribution of a polygon
     *  to the objective function minimized by a semi-discrete
     *  optimal transport map in 2D.
     */
    class OTMPolygonCallback :
	public OptimalTransportMap::Callback,
	public RVDPolygonCallback {
    public:
	
	/**
	 * \brief OTMPolygonCallback constructor.
	 * \param[in] OTM a pointer to the OptimalTransportMap2d
	 */
	OTMPolygonCallback(OptimalTransportMap2d* OTM) :
	    OptimalTransportMap::Callback(OTM) {
	}

	/**
	 * \copydoc RVDPolygonCallback::operator()
	 */
	virtual void operator() (
	    index_t v,
	    index_t t,
	    const GEOGen::Polygon& P
	) const {
	    geo_argused(t);
	    double F = 0.0;
	    Thread* thread = Thread::current();
	    index_t current_thread_id = (thread == nil) ? 0 : thread->id();

	    for(index_t i1=0; i1<P.nb_vertices(); ++i1) {
		index_t i2 = i1 + 1;
		if(i2 == P.nb_vertices()) {
		    i2 = 0;
		}
		// Note:
		// It is P.vertex(i2).adjacent_seed(),
		// not P.vertex(i1).adjacent_seed() !!!
		F += eval_simplex(
		    v, &P.vertex(i1), &P.vertex(i2),
		    index_t(P.vertex(i2).adjacent_seed())
		);
	    }
	    const_cast<OTMPolygonCallback*>(this)->
		funcval_[current_thread_id] += F;

	    if(mg_ != nil) {

		double mgx = 0.0;
		double mgy = 0.0;
		const GEOGen::Vertex& V0 = P.vertex(0);
		const double* p0 = V0.point();
		for(index_t i=1; i+1<P.nb_vertices(); ++i) {
		    const GEOGen::Vertex& V1 = P.vertex(i);
		    const double* p1 = V1.point();
		    const GEOGen::Vertex& V2 = P.vertex(i+1);
		    const double* p2 = V2.point();
		    double m = Geom::triangle_area(vec2(p0),vec2(p1),vec2(p2));
		    mgx += m * (p0[0] + p1[0] + p2[0]) / 3.0;
		    mgy += m * (p0[1] + p1[1] + p2[1]) / 3.0;
		}
		
		if(spinlocks_ != nil) {
		    spinlocks_->acquire_spinlock(v);
		}

		mg_[2*v] += mgx;  
		mg_[2*v+1] += mgy;
		
		if(spinlocks_ != nil) {
		    spinlocks_->release_spinlock(v);
		}
	    }
	}

    protected:

	/**
	 * \brief Evaluates the contribution of a simplex in non-weighted mode.
	 * \details The considered simplex connects one of the Delaunay vertices
	 *  (\p center_vertex_index) and two other vertices of the 
	 *  current cell (\p V1, \p V2).
	 * \param[in] center_vertex_index the index of the current 
	 *  Delaunay vertex
	 * \param[in] V1 , V2 three vertices of the current cell
	 * \param[in] v_adj the index of the Delaunay vertex of the cell 
	 *  on the other side of (\p V1, \p V2) or index_t(-1) if 
	 *  there is no such a cell.
	 */
	double eval_simplex(
	    index_t center_vertex_index,
	    const GEOGen::Vertex* V1,
	    const GEOGen::Vertex* V2,
	    index_t v_adj
	) const {
	    const double* q0 = OTM_->point_ptr(center_vertex_index);
	    const double* p1 = V1->point();
	    const double* p2 = V2->point();
	    double m = -GEO::Geom::triangle_signed_area(
		vec2(q0), vec2(p1), vec2(p2)
	    );
	    double fT = 0.0;
	    for(coord_index_t c = 0; c < 2; ++c) {
		double Uc = p1[c] - q0[c];
		double Vc = p2[c] - q0[c];
		fT +=
		    Uc * Uc +
		    Vc * Vc +
		    Uc * Vc ;
	    }
                
	    fT = m * (fT / 6.0 - w_[center_vertex_index]);
	    
	    double hij = 0.0;
	    if(Newton_step_ && v_adj != index_t(-1)) {
		const double* pj = OTM_->point_ptr(v_adj);
		hij =
		    GEO::Geom::distance(p1,p2,2) /
		    (2.0 * GEO::Geom::distance(q0,pj,2)) ;
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
		if(hij != 0.0) {
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
    };

    /********************************************************************/
}

namespace GEO {

    OptimalTransportMap2d::OptimalTransportMap2d(
        Mesh* mesh, const std::string& delaunay, bool BRIO
    ) :
	OptimalTransportMap(
	    2, 
	    mesh,
	    (delaunay == "default") ? "BPOW2d" : delaunay,
	    BRIO
        ) {
	callback_ = new OTMPolygonCallback(this);
	total_mass_ = total_mesh_mass();
    }

    OptimalTransportMap2d::~OptimalTransportMap2d() {
    }
    
    void OptimalTransportMap2d::get_RVD(Mesh& RVD_mesh) {
	ComputeRVDPolygonCallback callback(this, &RVD_mesh);
	RVD()->for_each_polygon(callback, false, false, false);
	/*
	  // NOTE: Does not work, TODO: determine why
        Attribute<index_t> tet_region(RVD_mesh.cells.attributes(),"region");
        RVD()->compute_RVD(
            RVD_mesh,
            0,             // dim (0 means use default)
            false,         // borders_only
            show_RVD_seed_ // integration_simplices
        );
	*/
    }

    void OptimalTransportMap2d::compute_Laguerre_centroids(double* centroids) {
        vector<double> g(nb_points(), 0.0);
        Memory::clear(centroids, nb_points()*sizeof(double)*2);

	callback_->set_Laguerre_centroids(centroids);
	callback_->set_g(g.data());
	{
	    Stopwatch* W = nil;
	    if(newton_ && verbose_) {
		W = new Stopwatch("RVD");
		Logger::out("OTM") << "In RVD (centroids)..." << std::endl;
	    }
	    RVD_->for_each_polygon(
		*dynamic_cast<RVDPolygonCallback*>(callback_), false, false, true
	    );
	    if(newton_ && verbose_) {
		delete W;
	    }
	}
	
	callback_->set_Laguerre_centroids(nil);	    
	
        for(index_t v=0; v<nb_points(); ++v) {
            centroids[2*v  ] /= g[v];
            centroids[2*v+1] /= g[v];
        }
    }

    double OptimalTransportMap2d::total_mesh_mass() const {
	double result = 0.0;
	
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
        
        for(index_t t = 0; t < mesh_->facets.nb(); ++t) {
            double tri_mass = GEO::Geom::triangle_area(
                vec2(mesh_->vertices.point_ptr(mesh_->facets.vertex(t, 0))),
		vec2(mesh_->vertices.point_ptr(mesh_->facets.vertex(t, 1))),
		vec2(mesh_->vertices.point_ptr(mesh_->facets.vertex(t, 2)))
            );
            if(vertex_mass.is_bound()) {
                tri_mass *= (
                    vertex_mass[mesh_->facets.vertex(t, 0)] +
                    vertex_mass[mesh_->facets.vertex(t, 1)] +
                    vertex_mass[mesh_->facets.vertex(t, 2)] 
                ) / 3.0;
            }
            result += tri_mass;
        }
	return result;
    }

    void OptimalTransportMap2d::call_callback_on_RVD() {
	RVD_->for_each_polygon(
	    *dynamic_cast<RVDPolygonCallback*>(callback_),false,false,true
	);
    }
    
    /**********************************************************************/
    
    void compute_Laguerre_centroids_2d(
        Mesh* omega,
        index_t nb_points,
        const double* points,
        double* centroids,
	bool parallel_pow,
	Mesh* RVD,
	bool verbose
    ) {
	geo_argused(parallel_pow); // Not implemented yet.

	// Omega can be either 2d or 3d with third coordinate set to
	// zero.
	index_t omega_dim_backup = omega->vertices.dimension();
        omega->vertices.set_dimension(3);

        // false = no BRIO
        // (OTM does not use multilevel and lets Delaunay
        //  reorder the vertices)
        OptimalTransportMap2d OTM(
	    omega,
	    std::string("BPOW2d"),
	    false
	);

        OTM.set_regularization(1e-3);
        OTM.set_Newton(true);
        OTM.set_points(nb_points, points);
        OTM.set_epsilon(0.01);
	OTM.set_Laguerre_centroids(centroids);
	OTM.set_verbose(verbose);
        OTM.optimize(1000);

	if(RVD != nil) {
	    OTM.get_RVD(*RVD);
	}
	
        omega->vertices.set_dimension(omega_dim_backup);
    }

    /**********************************************************************/    

    /**
     * \brief Computes the contribution of a polygon
     *  to the objective function minimized by a semi-discrete
     *  optimal transport map computed on a 3D surface.
     */
    class SurfaceOTMPolygonCallback :
	public OptimalTransportMap::Callback,
	public RVDPolygonCallback {
    public:
	
	/**
	 * \brief SurfaceOTMPolygonCallback constructor.
	 * \param[in] OTM a pointer to the OptimalTransportMap2d
	 */
	SurfaceOTMPolygonCallback(OptimalTransportMapOnSurface* OTM) :
	    OptimalTransportMap::Callback(OTM) {
	}

	/**
	 * \copydoc RVDPolygonCallback::operator()
	 */
	virtual void operator() (
	    index_t v,
	    index_t t,
	    const GEOGen::Polygon& P
	) const {
	    geo_argused(t);
	    double F = 0.0;
	    Thread* thread = Thread::current();
	    index_t current_thread_id = (thread == nil) ? 0 : thread->id();

	    index_t i0 = 0;
	    for(index_t i1=1; i1+1<P.nb_vertices(); ++i1) {
		index_t i2 = i1+1;
		F += eval_simplex(
		    v,
		    &P.vertex(i0), &P.vertex(i1), &P.vertex(i2)
		);
	    }

	    if(Newton_step_) {
		if(spinlocks_ != nil) {
		    spinlocks_->acquire_spinlock(v);
		}
	    
		for(index_t i1=0; i1<P.nb_vertices(); ++i1) {
		    index_t i2 = i1+1;
		    if(i2 == P.nb_vertices()) {
			i2 = 0;
		    }
		    
		    // Note:
		    // It is P.vertex(i2).adjacent_seed(),
		    // not P.vertex(i1).adjacent_seed() !!!
		    eval_border(
			v,
			&P.vertex(i1), &P.vertex(i2),
			index_t(P.vertex(i2).adjacent_seed())			
		    );
		}
		if(spinlocks_ != nil) {
		    spinlocks_->release_spinlock(v);
		}
	    }

	    
	    const_cast<SurfaceOTMPolygonCallback*>(this)->
		funcval_[current_thread_id] += F;

	    if(mg_ != nil) {

		double mgx = 0.0;
		double mgy = 0.0;
		double mgz = 0.0;
		
		const GEOGen::Vertex& V0 = P.vertex(0);
		const double* p0 = V0.point();
		for(index_t i=1; i+1<P.nb_vertices(); ++i) {
		    const GEOGen::Vertex& V1 = P.vertex(i);
		    const double* p1 = V1.point();
		    const GEOGen::Vertex& V2 = P.vertex(i+1);
		    const double* p2 = V2.point();
		    double m = triangle_mass(&V0,&V1,&V2);
		    if(weighted_) {
			double w0 = V0.weight();
			double w1 = V1.weight();
			double w2 = V2.weight();
			double s = m/(w0+w1+w2);
			mgx += s * (w0*p0[0] + w1*p1[0] + w2*p2[0]);
			mgy += s * (w0*p0[1] + w1*p1[1] + w2*p2[1]);
			mgz += s * (w0*p0[2] + w1*p1[2] + w2*p2[2]);
		    } else {
			mgx += m * (p0[0] + p1[0] + p2[0]) / 3.0;
			mgy += m * (p0[1] + p1[1] + p2[1]) / 3.0;
			mgz += m * (p0[2] + p1[2] + p2[2]) / 3.0;
		    }
		}
		
		if(spinlocks_ != nil) {
		    spinlocks_->acquire_spinlock(v);
		}

		mg_[3*v]   += mgx;  
		mg_[3*v+1] += mgy;
		mg_[3*v+2] += mgz;
		
		if(spinlocks_ != nil) {
		    spinlocks_->release_spinlock(v);
		}
	    }
	}

    protected:

	/**
	 * \brief Computes the mass of a triangle.
	 * \details Weights are taken into account in weighted_ mode.
	 * \param[in] V0 , V1 , V2 the three vertices of the triangle,
	 *  given as RVD vertices.
	 * \return the area of the triangle in 3D times the average value of the
	 *  three weights.
	 */
	double triangle_mass(
	    const GEOGen::Vertex* V0,
	    const GEOGen::Vertex* V1,
	    const GEOGen::Vertex* V2
	) const {
	    double m = Geom::triangle_area_3d(V0->point(), V1->point(), V2->point());
	    if(weighted_) {
		m *= ((V0->weight() + V1->weight() + V2->weight())/3.0);
	    } 
	    return m;
	}

	/**
	 * \brief Computes the mass of an edge.
	 * \details Weights are taken into account in weighted_ mode.
	 * \param[in] V0 , V1 the two vertices of the edge,
	 *  given as RVD vertices.
	 * \return the length of the edge in 3D times the average value of the
	 *  two weights.
	 */
	double edge_mass(
	    const GEOGen::Vertex* V0,
	    const GEOGen::Vertex* V1
	) const {
	    double m = Geom::distance(V0->point(), V1->point(), 3);
	    if(weighted_) {
		m *= ((V0->weight() + V1->weight())/2.0);
	    } 
	    return m;
	}
	
	
	/**
	 * \brief Evaluates the contribution of a simplex in non-weighted mode.
	 * \details The considered simplex connects three vertices of the 
	 *  current cell (\p V0, \p V1, \p V2).
	 * \param[in] center_vertex_index the index of the current 
	 *  Delaunay vertex
	 * \param[in] V0 , V1 , V2 three vertices of the current cell
	 * \param[in] v_adj the index of the Delaunay vertex of the cell 
	 *  on the other side of (\p V1, \p V2) or index_t(-1) if 
	 *  there is no such a cell.
	 */
	double eval_simplex(
	    index_t center_vertex_index,
	    const GEOGen::Vertex* V0,
	    const GEOGen::Vertex* V1,
	    const GEOGen::Vertex* V2
	) const {
	    
	    const double* q0 = OTM_->point_ptr(center_vertex_index);
	    const double* p0 = V0->point();
	    const double* p1 = V1->point();
	    const double* p2 = V2->point();
	    
	    double m = triangle_mass(V0,V1,V2);

	    // TODO: check formula for fT...
	    // TODO: formula is not correct in weighted mode.
	    // (note: it is * not used * in Newton mode, the  value
	    //  of the objective function is only used in multilevel mode).
	    
	    double fT = 0.0;
	    for(coord_index_t c = 0; c < 3; ++c) {
		double Uc = p0[c] - q0[c];		
		double Vc = p1[c] - q0[c];
		double Wc = p2[c] - q0[c];
		fT +=
		    Uc * Uc +
		    Vc * Vc +
		    Wc * Wc +
		    Uc * Vc +
		    Uc * Wc +
		    Vc * Wc ;
	    }
	    fT = m * (fT / 30.0 - w_[center_vertex_index]);
	    
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
	    }

	    if(spinlocks_ != nil) {
		spinlocks_->release_spinlock(center_vertex_index);
	    }
	    
            // -fT because we maximize F <=> minimize -F
            return -fT;
	}

	void eval_border(
	    index_t center_vertex_index,
	    const GEOGen::Vertex* V0,
	    const GEOGen::Vertex* V1,
	    index_t v_adj
	) const {
	    if(v_adj == index_t(-1)) {
		return;
	    }
	    const double* q0 = OTM_->point_ptr(center_vertex_index);	    
	    const double* pj = OTM_->point_ptr(v_adj);
	    
	    double hij =
		edge_mass(V0,V1) /
		(2.0 * GEO::Geom::distance(q0,pj,3));
	
	    // -hij because we maximize F <=> minimize -F
	    if(hij != 0.0) {
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
    };
    
    /**********************************************************************/
    
    OptimalTransportMapOnSurface::OptimalTransportMapOnSurface(
        Mesh* mesh, const std::string& delaunay, bool BRIO
    ) :
	OptimalTransportMap(
	    3, 
	    mesh,
	    (delaunay == "default") ? "BPOW" : delaunay,
	    BRIO
        ) {
	callback_ = new SurfaceOTMPolygonCallback(this);
	total_mass_ = total_mesh_mass();
    }

    OptimalTransportMapOnSurface::~OptimalTransportMapOnSurface() {
    }
    
    void OptimalTransportMapOnSurface::get_RVD(Mesh& RVD_mesh) {
	ComputeRVDPolygonCallback callback(this, &RVD_mesh);
	RVD()->for_each_polygon(callback, false, false, false);
	/*
	  // NOTE: Does not work, TODO: determine why
        Attribute<index_t> tet_region(RVD_mesh.cells.attributes(),"region");
        RVD()->compute_RVD(
            RVD_mesh,
            0,             // dim (0 means use default)
            false,         // borders_only
            show_RVD_seed_ // integration_simplices
        );
	*/
    }

    void OptimalTransportMapOnSurface::compute_Laguerre_centroids(double* centroids) {
        vector<double> g(nb_points(), 0.0);
        Memory::clear(centroids, nb_points()*sizeof(double)*3);

	callback_->set_Laguerre_centroids(centroids);
	callback_->set_g(g.data());
	{
	    Stopwatch* W = nil;
	    if(newton_ && verbose_) {
		W = new Stopwatch("RVD");
		Logger::out("OTM") << "In RVD (centroids)..." << std::endl;
	    }
	    RVD_->for_each_polygon(
		*dynamic_cast<RVDPolygonCallback*>(callback_), false, false, true
	    );
	    if(newton_ && verbose_) {
		delete W;
	    }
	}
	
	callback_->set_Laguerre_centroids(nil);	    
	
        for(index_t v=0; v<nb_points(); ++v) {
            centroids[3*v  ] /= g[v];
            centroids[3*v+1] /= g[v];
            centroids[3*v+2] /= g[v];	    
        }
    }

    double OptimalTransportMapOnSurface::total_mesh_mass() const {
	double result = 0.0;
	
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
        
        for(index_t t = 0; t < mesh_->facets.nb(); ++t) {
            double tri_mass = GEO::Geom::triangle_area(
                vec3(mesh_->vertices.point_ptr(mesh_->facets.vertex(t, 0))),
		vec3(mesh_->vertices.point_ptr(mesh_->facets.vertex(t, 1))),
		vec3(mesh_->vertices.point_ptr(mesh_->facets.vertex(t, 2)))
            );
            if(vertex_mass.is_bound()) {
                tri_mass *= (
                    vertex_mass[mesh_->facets.vertex(t, 0)] +
                    vertex_mass[mesh_->facets.vertex(t, 1)] +
                    vertex_mass[mesh_->facets.vertex(t, 2)] 
                ) / 3.0;
            }
            result += tri_mass;
        }
	return result;
    }

    void OptimalTransportMapOnSurface::call_callback_on_RVD() {
	RVD_->for_each_polygon(
	    *dynamic_cast<RVDPolygonCallback*>(callback_),false,false,true
	);
    }

    /**********************************************************************/    
    
}

