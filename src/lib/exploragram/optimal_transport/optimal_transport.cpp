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

#include <exploragram/optimal_transport/optimal_transport.h>
#include <exploragram/optimal_transport/linear_least_squares.h>

#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_AABB.h>

#include <geogram/voronoi/CVT.h>
#include <geogram/voronoi/generic_RVD_vertex.h>
#include <geogram/voronoi/RVD_callback.h>
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


namespace GEO {

    OptimalTransportMap* OptimalTransportMap::instance_ = nil;

    OptimalTransportMap::Callback::~Callback() {
    }
    
    OptimalTransportMap::OptimalTransportMap(
	index_t dimension,
        Mesh* mesh, const std::string& delaunay, bool BRIO
    ) : mesh_(mesh) {

	geo_cite("DBLP:conf/compgeom/AurenhammerHA92");
	geo_cite("DBLP:journals/cgf/Merigot11");
	geo_cite("journals/M2AN/LevyNAL15");

	dimension_ = dimension;
	dimp1_ = dimension_+1;
	
        epsilon_regularization_ = 0.0;

        // Mesh is supposed to be embedded in d+1 dim (with
        // (d+1-th dimension set to zero).
        geo_assert(mesh->vertices.dimension() == dimp1_);

        // Note: we represent power diagrams as d+1 Voronoi diagrams
        delaunay_ = Delaunay::create(coord_index_t(dimp1_), delaunay);

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
	verbose_ = false;
        
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
	callback_ = nil;
	Laguerre_centroids_ = nil;
    }

    OptimalTransportMap::~OptimalTransportMap() {
	delete callback_;
	callback_ = nil;
    }
    
    void OptimalTransportMap::set_points(
        index_t nb_points, const double* points, index_t stride
    ) {
	if(stride == 0) {
	    stride = dimension_;
	}
	
        // Note: we represent power diagrams as (d+1)-dim Voronoi diagrams.
        // The target points are lifted to (d+1)-dim.
        points_dimp1_.resize(nb_points * dimp1_);
        for(index_t i = 0; i < nb_points; ++i) {
	    for(index_t c=0; c<dimension_; ++c) {
		points_dimp1_[i*dimp1_+c] = points[i*stride+c];
	    }
            points_dimp1_[i*dimp1_ + dimension()] = 0.0;
        }
        weights_.assign(nb_points, 0);
        lambda_p_ = total_mass_ / double(nb_points);
    }

    // See http://arxiv.org/abs/1603.05579
    // Kitawaga, Merigot, Thibert,
    // A Newton Algorithm for semi-discrete OT
    
    void OptimalTransportMap::optimize_full_Newton(
        index_t max_iterations, index_t n
    ) {
        if(n == 0) {
            n = index_t(points_dimp1_.size() / dimp1_);
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

                // Compute cell measures and nbZ.
		if(Laguerre_centroids_ != nil) {
		    callback_->set_Laguerre_centroids(Laguerre_centroids_);
		}
		
		funcgrad(n,weights_.data(),fk,gk.data());

		if(Laguerre_centroids_ != nil) {
		    callback_->set_Laguerre_centroids(nil);
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
    
    void OptimalTransportMap::optimize(index_t max_iterations) {

        level_ = 0;
        
        if(newton_) {
            optimize_full_Newton(max_iterations);
            return;
        }
        
        index_t n = index_t(points_dimp1_.size() / dimp1_);
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


    void OptimalTransportMap::optimize_level(
        index_t b, index_t e, index_t max_iterations
    ) {

        // If this is not the first level, propagate the weights from
        // the lower levels.
        if(b != 0) {
            
            //   Create a nearest neighbor search data structure
            // and insert the [0..b) samples into it (they were
            // initialized at previous calls).
            NearestNeighborSearch_var NN =
		NearestNeighborSearch::create(coord_index_t(dimension()));
	    
            NN->set_points(b, points_dimp1_.data(), dimp1_);
            index_t degree = 2; // CmdLine::get_arg_uint("fitting_degree");
            
            // If degree \notin {1,2}, use weight of nearest sample
            if(degree < 1 || degree > 2) {
                for(index_t i = b; i < e; ++i) {
                    weights_[i] =
                        weights_[
                            NN->get_nearest_neighbor(&points_dimp1_[dimp1_ * i])
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
                        nb, &points_dimp1_[dimp1_ * i], neighbor, dist
                    );
                    LLS.begin();
                    for(index_t jj = 0; jj < nb; ++jj) {
                        if(dist[jj] != 0.0) {
                            index_t j = neighbor[jj];
                            LLS.add_point(&points_dimp1_[dimp1_ * j], weights_[j]);
                        }
                    }
                    LLS.end();
                    weights_[i] = LLS.eval(&points_dimp1_[dimp1_ * i]);
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

    void OptimalTransportMap::optimize_levels(
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
    
    void OptimalTransportMap::funcgrad_CB(
        index_t n, double* x, double& f, double* g
    ) {
        instance_->funcgrad(n, x, f, g);
    }

    void OptimalTransportMap::newiteration_CB(
        index_t n, const double* x, double f, const double* g, double gnorm
    ) {
        geo_argused(n);
        geo_argused(x);
        geo_argused(f);
        geo_argused(g);
        geo_argused(gnorm);
        instance_->newiteration();
    }

    void OptimalTransportMap::newiteration() {
        //xxx std::cerr << "newiteration" << std::endl;
        if(save_RVD_iter_) {
            std::cerr << "  save iter" << std::endl;
            save_RVD(current_iter_);
        }
        ++current_iter_;
    }

    void OptimalTransportMap::save_RVD(index_t id) {
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

    void OptimalTransportMap::funcgrad(
        index_t n, double* w, double& f, double* g
    ) {

        bool is_Newton_step = callback_->is_Newton_step();	    
        
        // For now, always compute function and gradient
        bool update_fg = true;
        
        // Delaunay triangulation is only updated if function and
        // gradient is evaluated. If only Hessian needs to be evaluated,
        // then it is at the same point as the latest function and gradient
        // evaluation (see Yang Liu's CVT-Newton code).
        if(update_fg && !w_did_not_change_) {
            // Step 1: determine the (dim+1)d embedding from the weights
            double W = 0.0;
            for(index_t p = 0; p < n; ++p) {
                W = geo_max(W, w[p]);
            }
            for(index_t p = 0; p < n; ++p) {
                points_dimp1_[dimp1_ * p + dimension_] = ::sqrt(W - w[p]);
            }
        
            // Step 2: compute function and gradient
            {
		Stopwatch* SW = nil;
		if(newton_) {
		    SW = new Stopwatch("Power diagram");
		    Logger::out("OTM") << "In power diagram..." << std::endl;
		}
                delaunay_->set_vertices(n, points_dimp1_.data());
		if(newton_) {
		    delete SW;
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

	callback_->set_w(w,n);
	callback_->set_g(g);
	callback_->set_nb_threads(Process::maximum_concurrent_threads());

	if(callback_->has_Laguerre_centroids()) {
	    Memory::clear(callback_->Laguerre_centroids(), nb_points()*sizeof(double)*dimension());	    
	}
	
	{
	    Stopwatch* W = nil;
	    if(newton_) {
		W = new Stopwatch("RVD");
		Logger::out("OTM") << "In RVD (funcgrad)..." << std::endl;
	    }
	    call_callback_on_RVD();
	    if(newton_) {
		delete W;
	    }
	}
	f = callback_->funcval(); 

	if(callback_->has_Laguerre_centroids()) {
	    for(index_t v=0; v<nb_points(); ++v) {
		for(index_t c=0; c<dimension_; ++c) {
		    callback_->Laguerre_centroids()[dimension_*v+c] /= g[v];
		}
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

    void OptimalTransportMap::eval_func_grad_Hessian(
        index_t n, const double* w, double& f, double* g
    ) {
	callback_->set_Newton_step(true);
	funcgrad(n,(double*)w,f,g);
	callback_->set_Newton_step(false);
    }
    
/************************************************************/

    // TODO: use OpenNL buffers to avoid data copy and allocation.
    // TODO: in the Euler code, see if we do not have duplicated computations, i.e.
    //    - Power diagrams when leaving and entering iteration ?
    //    - Centroids: do we restart a RVD computation ?
    
    void OptimalTransportMap::update_sparsity_pattern() {
        // Does nothing for now,
	// (we let OpenNL discover the sparsity pattern)
	// Tryed smarter things, but was not faster...
    }

    void OptimalTransportMap::new_linear_system(index_t n) {
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

    void OptimalTransportMap::solve_linear_system(double* x) {
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
}
