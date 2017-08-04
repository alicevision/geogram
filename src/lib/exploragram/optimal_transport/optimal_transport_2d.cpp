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
#include <exploragram/optimal_transport/linear_least_squares.h>

#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/mesh/mesh_geometry.h>

#include <geogram/voronoi/CVT.h>
#include <geogram/voronoi/generic_RVD_vertex.h>
#include <geogram/delaunay/delaunay_2d.h>

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

namespace {
    using namespace GEO;

    /**
     * \brief Computes the contribution of a simplex.
     *  to the objective function minimized by a semi-discrete
     *  optimal transport map.
     */
    class OTMIntegrationSimplex : public IntegrationSimplex {
    public:
    public:
	/**
	 * \brief OTMIntegrationSimplex constructor.
	 * \param[in] OTM a pointer to the OptimalTransportMap3d
	 */
	OTMIntegrationSimplex(OptimalTransportMap2d* OTM) :
	    IntegrationSimplex(OTM->mesh(), false, 0, 0, nil),            
	    OTM_(OTM),
	    Newton_step_(false),
	    w_(nil),
	    mg_(nil)
	{
	    weighted_ =
		OTM->mesh().vertices.attributes().is_defined("weight");
	}


	/**
	 * \copydoc IntegrationSimplex::eval()
	 */
	virtual double eval(
	    index_t center_vertex_index,
	    const GEOGen::Vertex& v0,
	    const GEOGen::Vertex& v1,
	    const GEOGen::Vertex& v2,
	    index_t t,
	    index_t t_adj = index_t(-1),
	    index_t v_adj = index_t(-1)
	) {
	    geo_argused(center_vertex_index);
	    geo_argused(v0);
	    geo_argused(v1);
	    geo_argused(v2);
	    geo_argused(t);
	    geo_argused(t_adj);
	    geo_argused(v_adj);
	    return 0.0;
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
	    geo_argused(n);
        }

	void set_Newton_step(bool Newton) {
	    Newton_step_ = Newton;
	}

	bool is_Newton_step() const {
	    return Newton_step_;
	}

    protected:

    private:
	OptimalTransportMap2d* OTM_;
	bool weighted_;
	bool Newton_step_;
	const double* w_;
	double* mg_;
    };
}

namespace GEO {

    OptimalTransportMap2d* OptimalTransportMap2d::instance_ = nil;
    
    OptimalTransportMap2d::OptimalTransportMap2d(
        Mesh* mesh, const std::string& delaunay_in, bool BRIO
    ) : mesh_(mesh) {

	geo_cite("DBLP:conf/compgeom/AurenhammerHA92");
	geo_cite("DBLP:journals/cgf/Merigot11");
	geo_cite("journals/M2AN/LevyNAL15");

	std::string delaunay = delaunay_in;
	if(delaunay == "default") {
	    delaunay = "BPOW2d";
	}
	
        epsilon_regularization_ = 0.0;

        // Mesh is supposed to be embedded in 3 dim (with
        // 3rd dimension set to zero).
        geo_assert(mesh->vertices.dimension() == 0);

        // Note: we represent power diagrams as 3d Voronoi diagrams
        delaunay_ = Delaunay::create(3, delaunay);

        RVD_ = RestrictedVoronoiDiagram::create(delaunay_, mesh_);
        RVD_->set_volumetric(true);
        RVD_->set_check_SR(true);
        RVD_->create_threads();
        
        //   No need to reorder vertices if BRIO is activated since
        // vertices are then already reordered.
        if(BRIO) {
            RVD_->delaunay()->set_reorder(false);
        }

        newton_ = true;
        
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

	integration_simplex_ = new OTMIntegrationSimplex(this);

	Laguerre_centroids_ = nil;
    }

    OptimalTransportMap2d::~OptimalTransportMap2d() {
    }
    
    void OptimalTransportMap2d::set_points(
        index_t nb_points, const double* points
    ) {
        // Note: we represent power diagrams as 3d Voronoi diagrams.
        // The target points are lifted to 3d.
        points_3d_.resize(nb_points * 3);
        for(index_t i = 0; i < nb_points; ++i) {
            points_3d_[i * 3] = points[i * 3];
            points_3d_[i * 3 + 1] = points[i * 3 + 1];
            points_3d_[i * 3 + 2] = 0.0;
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
            total_mass_ += tri_mass;
        }
        lambda_p_ = total_mass_ / double(nb_points);
    }


    // See http://arxiv.org/abs/1603.05579
    // Kitawaga, Merigot, Thibert,
    // A Newton Algorithm for semi-discrete OT
    
    void OptimalTransportMap2d::optimize_full_Newton(
        index_t max_iterations, index_t n
    ) {
        if(n == 0) {
            n = index_t(points_3d_.size() / 3);
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

		OTMIntegrationSimplex* smplx_otm =
		    dynamic_cast<OTMIntegrationSimplex*>(integration_simplex_.get());
                // Compute cell measures and nbZ.
		if(Laguerre_centroids_ != nil) {
		    smplx_otm->set_Laguerre_centroids(Laguerre_centroids_);
		}
		
		funcgrad(n,weights_.data(),fk,gk.data());

		if(Laguerre_centroids_ != nil) {
		    smplx_otm->set_Laguerre_centroids(nil);
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
    
    void OptimalTransportMap2d::optimize(index_t max_iterations) {
        level_ = 0;
        if(newton_) {
            optimize_full_Newton(max_iterations);
            return;
        }
    }

    void OptimalTransportMap2d::get_RVD(Mesh& RVD_mesh) {
        RVD_mesh.clear();
        Attribute<index_t> tet_region(RVD_mesh.cells.attributes(),"region");
        RVD()->compute_RVD(
            RVD_mesh,
            0,             // dim (0 means use default)
            false,         // borders_only
            show_RVD_seed_ // integration_simplices
        );
    }

    void OptimalTransportMap2d::compute_Laguerre_centroids(double* centroids) {
        vector<double> g(nb_points(), 0.0);
        Memory::clear(centroids, nb_points()*sizeof(double)*2);

	OTMIntegrationSimplex* smplx_otm =
	    dynamic_cast<OTMIntegrationSimplex*>(integration_simplex_.get());
	
	smplx_otm->set_Laguerre_centroids(centroids);
	{
	    double f = 0.0;
	    Stopwatch* W = nil;
	    if(newton_) {
		W = new Stopwatch("RVD");
		Logger::out("OTM") << "In RVD (centroids)..." << std::endl;
	    }
	    RVD_->compute_integration_simplex_func_grad(f, g.data(), smplx_otm);
	    if(newton_) {
		delete W;
	    }
	}
	
	smplx_otm->set_Laguerre_centroids(nil);	    
	
        for(index_t v=0; v<nb_points(); ++v) {
            centroids[2*v  ] /= g[v];
            centroids[2*v+1] /= g[v];
        }
    }
    
    void OptimalTransportMap2d::save_RVD(index_t id) {
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

    void OptimalTransportMap2d::funcgrad(
        index_t n, double* w, double& f, double* g
    ) {

	OTMIntegrationSimplex* smplx_otm =
	    dynamic_cast<OTMIntegrationSimplex*>(integration_simplex_.get());
	
        bool is_Newton_step = smplx_otm->is_Newton_step();	    
        
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
                points_3d_[3 * p + 2] = ::sqrt(W - w[p]);
            }
        
            // Step 2: compute function and gradient
            {
		Stopwatch* W = nil;
		if(newton_) {
		    W = new Stopwatch("Power diagram");
		    Logger::out("OTM") << "In power diagram..." << std::endl;
		}
                delaunay_->set_vertices(n, points_3d_.data());
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

	smplx_otm->set_w(w,n);

	if(smplx_otm->has_Laguerre_centroids()) {
	    Memory::clear(smplx_otm->Laguerre_centroids(), nb_points()*sizeof(double)*3);	    
	}
	
	{
	    Stopwatch* W = nil;
	    if(newton_) {
		W = new Stopwatch("RVD");
		Logger::out("OTM") << "In RVD (funcgrad)..." << std::endl;
	    }
	    RVD_->compute_integration_simplex_func_grad(
                f, g, smplx_otm
            );
	    if(newton_) {
		delete W;
	    }
	}

	if(smplx_otm->has_Laguerre_centroids()) {
	    for(index_t v=0; v<nb_points(); ++v) {
		smplx_otm->Laguerre_centroids()[2*v  ] /= g[v];
		smplx_otm->Laguerre_centroids()[2*v+1] /= g[v];
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

    void OptimalTransportMap2d::eval_func_grad_Hessian(
        index_t n, const double* w, double& f, double* g
    ) {
	OTMIntegrationSimplex* smplx_otm =
	    dynamic_cast<OTMIntegrationSimplex*>(integration_simplex_.get());
	smplx_otm->set_Newton_step(true);
	funcgrad(n,(double*)w,f,g);
	smplx_otm->set_Newton_step(false);
    }
    
/************************************************************/

    // TODO: use OpenNL buffers to avoid data copy and allocation.
    // TODO: in the Euler code, see if we do not have duplicated computations, i.e.
    //    - Power diagrams when leaving and entering iteration ?
    //    - Centroids: do we restart a RVD computation ?
    
    void OptimalTransportMap2d::update_sparsity_pattern() {
        // Does nothing for now,
	// (we let OpenNL discover the sparsity pattern)
	// Tryed smarter things, but was not faster...
    }

    void OptimalTransportMap2d::new_linear_system(index_t n) {
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

    void OptimalTransportMap2d::solve_linear_system(double* x) {
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

    void OptimalTransportMap2d::newiteration() {
        //xxx std::cerr << "newiteration" << std::endl;
        if(save_RVD_iter_) {
            std::cerr << "  save iter" << std::endl;
            save_RVD(current_iter_);
        }
        ++current_iter_;
    }
    
    /**********************************************************************/
    
    void compute_Laguerre_centroids_2d(
        Mesh* omega,
        index_t nb_points,
        const double* points,
        double* centroids,
	bool parallel_pow
    ) {
	geo_argused(parallel_pow); // Not implemented yet.
	
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
        OTM.optimize(1000);
        omega->vertices.set_dimension(3);        
    }

    
}
