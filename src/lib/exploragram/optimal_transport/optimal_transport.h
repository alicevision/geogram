
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

#ifndef H_EXPLORAGRAM_OPTIMAL_TRANSPORT_OPTIMAL_TRANSPORT_H
#define H_EXPLORAGRAM_OPTIMAL_TRANSPORT_OPTIMAL_TRANSPORT_H

#include <exploragram/basic/common.h>
#include <geogram/mesh/mesh.h>
#include <geogram/voronoi/RVD.h>
#include <geogram/voronoi/integration_simplex.h>
#include <geogram/delaunay/delaunay.h>
#include <geogram/third_party/HLBFGS/HLBFGS.h>

/**
 * \file exploragram/optimal_transport/optimal_transport.h
 * \brief Solver for semi-discrete optimal transport (
 *  multilevel and Newton).
 */

namespace OGF {
    class SceneGraph;
}

namespace GEO {
    class CentroidalVoronoiTesselation;
    class SparseMatrix;

    /**
     * \brief Computes the centroids of the Laguerre cells that
     *  correspond to optimal transport.
     * \param[in] omega a reference to the mesh that represents the
     *  domain
     * \param[in] nb_points number of points
     * \param[in] points a pointer to the coordinates of the points
     * \param[out] centroids a pointer to the computed centroids of 
     *  the Laguerre cells that correspond to the optimal transport of
     *  the uniform measure to the points
     * \param[in] parallel_pow if true, use parallel power diagram algorithm
     */
    void EXPLORAGRAM_API compute_Laguerre_centroids(
        Mesh* omega,
        index_t nb_points,
        const double* points,
        double* centroids,
	bool parallel_pow=true
    );
    
    /**
     * \brief Computes semi-discrete optimal transport maps.
     * \details Computes an optimal transport map between two
     *  distributions in 3D. The first distribution is represented
     *  by a 3D tetrahedral mesh. The second distribution is a sum
     *  of Diracs.
     *  The algorithm is described in the following references:
     *   - 3D algorithm: http://arxiv.org/abs/1409.1279
     *   - Earlier 2D version by Quentin M\'erigot: 
     *    Q. Merigot. A multiscale approach to optimal transport.
     *    Computer Graphics Forum 30 (5) 1583--1592, 2011 (Proc SGP 2011).
     *   - Earlier article on OT and power diagrams: 
     *    F. Aurenhammer, F. Hoffmann, and B. Aronov. Minkowski-type theorems 
     *    and least-squares clustering. Algorithmica, 20:61-76, 1998.
     */
    class EXPLORAGRAM_API OptimalTransportMap {
    public:
        /**
         * \brief Initializes a new OptimalTransportMap.
         * \param[in] mesh the source distribution, represented as a 3d mesh
         * \param[in] delaunay factory name of the Delaunay triangulation, one
	 *  of "PDEL" (parallel), "BPOW" (sequential)
         * \param[in] BRIO true if vertices are already ordered using BRIO
         */
        OptimalTransportMap(
            Mesh* mesh,
            const std::string& delaunay = "PDEL",
            bool BRIO = false
        );

        /**
         * \brief Gets the mesh.
         * \return a reference to the mesh
         */
        Mesh& mesh() {
            return *mesh_;
        }

        /**
         * \brief Sets whether Newton algorithm should be used.
         * \details It is (for now) incompatible with multilevel.
         * \param[in] x if set, Newton algorithm is used instead
         *  of BFGS. 
         */
        void set_Newton(bool x) {
            newton_ = x;
        }
        
        /**
         * \brief Sets the points that define the target distribution.
         * \param[in] nb_points number of points in the target distribution
         * \param[in] points coordinates of the Diracs centers in the target
         *  distribution.
         */
        void set_points(index_t nb_points, const double* points);

        /**
         * \brief Sets the maximum error.
         * \param eps acceptable relative deviation for the measure of a
         *   Voronoi cell.
         */
        void set_epsilon(double eps) {
            epsilon_ = eps;
        }

        /**
         * \brief Sets the weight of the regularization term.
         * \details The regularization term (norm of the weight vector) cancels
         *  the translational degree of freedom of the weights.
         * \param[in] eps_reg the weight of the regularization term. Use 0.0 for
         *  no regularization.
         */
        void set_regularization(double eps_reg) {
            epsilon_regularization_ = eps_reg;
        }

        /**
         * \brief Computes the weights that realize the optimal
         *  transport map between the source mesh and the target
         *  pointset.
         * \param[in] max_iterations maximum number of solver iterations.
         */
        void optimize(index_t max_iterations);

        /**
         * \brief Computes the weights that realize the optimal
         *  transport map between the source mesh and the target
         *  pointset.
         * \param[in] max_iterations maximum number of solver iterations.
         * \param[in] n number of weights to optimize, used in hierarchical
         *  mode. If zero, optimizes all the weights.
         */
        void optimize_full_Newton(index_t max_iterations, index_t n=0);
        
        /**
         * \brief Optimizes one level of the multilevel algorithm.
         * \details The function supposes that the sequence [0,b)
         *  has been previously optimized. It is used to initialize
         *  the sequence [b,e). The whole sequence [0,e) is then
         *  optimized.
         * \param[in] b index fo the first point in the level
         * \param[in] e one position past the last index of the level
         * \param[in] max_iterations maximum number of iterations
         */
        void optimize_level(index_t b, index_t e, index_t max_iterations);

        /**
         * \brief Multi-level optimization.
         * \details The points specified by set_points() need to have
         *   a hierarchical structure. They can be constructed by
         *   compute_hierarchical_sampling().
         * \param[in] levels sample indices that correspond to level l are
         *   in the range levels[l] (included) ... levels[l+1] (excluded)
         * \param[in] max_iterations maximum number of iterations
         * \see compute_hierarchical_sampling()
         */
        void optimize_levels(
            const vector<index_t>& levels, index_t max_iterations
        );

        /**
         * \brief Gets the number of points.
         * \return The number of points, that was previously defined
         *  by set_points()
         */
        index_t nb_points() const {
            return weights_.size();
        }

        /**
         * \brief Gets a point.
         * \param[in] i index of the point
         * \return a const pointer to the coordinates of the 4d point \p i
         */
        const double* point_ptr(index_t i) const {
            geo_debug_assert(i < nb_points());
            return &(points_4d_[4 * i]);
        }

        /**
         * \brief Gets weight of a point.
         * \param[in] i index of the point
         * \return the weight that was computed for point \p i
         */
        double weight(index_t i) const {
            return weights_[i];
        }

        /**
         * \brief Sets a weight of a point.
         * \param[in] i index of the point
         * \param[in] val new value of the weight
         */
        void set_weight(index_t i, double val) {
            weights_[i] = val;
        }
        
        /**
         * \brief Gets the d+1-th coordinate of the embedding for a point.
         * \param[in] i index of the point
         * \return the d+1-th coordinate that was computed for point \p i
         */
        double potential(index_t i) const {
            return points_4d_[4 * i + 3];
        }

        /**
         * \brief Callback for the numerical solver.
         * \details Evaluates the objective function and its gradient.
         * \param[in] n number of variables
         * \param[in] x current value of the variables
         * \param[out] f current value of the objective function
         * \param[out] g gradient of the objective function
         */
        static void funcgrad_CB(
            index_t n, double* x, double& f, double* g
        );

        /**
         * \brief Callback for the numerical solver.
         * \param[in] n number of variables
         * \param[in] x current value of the variables
         * \param[in] f current value of the objective function
         * \param[in] g gradient of the objective function
         * \param[in] gnorm norm of the gradient of the objective function
         */
        static void newiteration_CB(
            index_t n, const double* x, double f, const double* g, double gnorm
        );

        /**
         * \brief Gets the restricted Voronoi diagram.
         * \return a pointer to the restricted Voronoi diagram
         */
        RestrictedVoronoiDiagram* RVD() {
            return RVD_;
        }

        /**
         * \brief Sets whether the restricted Voronoi diagram at
         *  each iteration should be saved.
         * \details If flag is set, then each iteration is saved 
         *  in file "RVD_nnn.geogram".
         * \param[in] x true if each iteration should be saved, 
         *  false otherwise.
         * \param[in] scene_graph if specified, for each iteration
         *  an object is created in the specified scene graph,
         *  else iterations are saved in files.
         * \param[in] show_RVD_seed if true, the seed associated
         *  with each restricted Voronoi cell is connected to it
         * \param[in] last_iter_only if true, only the last iteration
         *  is saved
         */
        void set_save_RVD_iter(
            bool x, OGF::SceneGraph* scene_graph = nil,
            bool show_RVD_seed = false,
            bool last_iter_only = false
        ) {
            if(last_iter_only) {
                save_RVD_iter_ = false;
                save_RVD_last_iter_ = true;
            } else {
                save_RVD_iter_ = x;
            }
            scene_graph_ = scene_graph;
            show_RVD_seed_ = show_RVD_seed;
        }

        void get_RVD(Mesh& M);

        /**
         * \brief Computes the centroids of the Laguerre cells.
         * \param[out] centroids a pointer to the 3*nb_points coordinates
         *  of the centroids.
         */
        void compute_Laguerre_centroids(double* centroids);

        /**
         * \brief Updates the sparsity pattern of the Hessian right after
         *  a new Laguerre diagram was computed.
         */
        void update_sparsity_pattern();
        
        /**
         * \brief Starts a new linear system.
         * \param[in] n the dimension of the system
         */
        void new_linear_system(index_t n);

        /**
         * \brief Adds a coefficient to the matrix of the system.
         * \param[in] i , j the indices of the coefficient
         * \param[in] a the value to be added to the coefficient
         */
        void add_ij_coefficient(index_t i, index_t j, double a);

        /**
         * \brief Adds a coefficient to the right hand side.
         * \param[in] i the index of the coefficient
         * \param[in] a the value to be added to the coefficient
         */
        void add_i_right_hand_side(index_t i, double a);
        
        /**
         * \brief Solves a linear system.
         * \param[out] x a pointer to the solution of the linear system.
         */
        void solve_linear_system(double* x);

    protected:

        /**
         * \brief Callback for the numerical solver.
         */
        void newiteration();
        
        /**
         * \brief Saves the RVD at each iteration if
         *   specified on command line (just for debugging/
         *   explaining the algorithm).
         * \param[in] id index to be used for the file, that
         *   will be named RVD_id.meshb
         */
        void save_RVD(index_t id);

        /**
         * \brief Computes the objective function and its gradient.
         * \param[in] n number of variables
         * \param[in] w current value of the variables
         * \param[out] f current value of the objective function
         * \param[out] g gradient of the objective function
         */
        void funcgrad(index_t n, double* w, double& f, double* g);

        /**
         * \brief Computes the objective function, its gradient and its Hessian.
         * \details Gradient and Hessian are used to solve a Newton
         *  step H p = -g
         * \param[in] n number of variables
         * \param[in] w current value of the variables
         * \param[out] f current value of the objective function
         * \param[out] g gradient of the objective function
         */
        void eval_func_grad_Hessian(
            index_t n, const double* w,
            double& f, double* g
        );

        /**
         * \brief Computes the stopping criterion of the solver.
         * \details The stopping criterion is determined from
         *  the user-specified epsilon, number of samples and
         *  target measure of a cell (lambda_p_).
         * \param n number of samples
         * \return the gradient threshold
         * \see set_epsilon()
         */
        double gradient_threshold(index_t n) const {
            return ::sqrt(double(n) * geo_sqr(epsilon_ * lambda_p_));
        }

    private:
        static OptimalTransportMap* instance_;
        Mesh* mesh_;
        Delaunay_var delaunay_;
        RestrictedVoronoiDiagram_var RVD_;
        vector<double> points_4d_;
        vector<double> weights_;
        double total_mass_;
        double lambda_p_; /**< \brief Value of one of the Diracs */
        double epsilon_;
        /**< \brief Acceptable relative deviation for the measure of a cell */
        index_t current_call_iter_;
        IntegrationSimplex_var simplex_func_;
        std::string last_stats_;
        bool pretty_log_;
        index_t level_;

        bool save_RVD_iter_;
        bool save_RVD_last_iter_;
        bool show_RVD_seed_;
        index_t current_iter_;
	OGF::SceneGraph* scene_graph_;

        bool newton_;

        /**
         * \brief Add a regularization term to remove 
         *  translational degree of freedom for the
         *  weights.
         */
        double epsilon_regularization_;

        /**
         * \brief Number of empty cells in last iteration.
         */
        index_t nbZ_;

        /**
         * \brief Norm of the gradient in last iteration.
         */
        double g_norm_;

        /**
         * \brief Measure of the smallest Laguerre cell.
         */
        double measure_of_smallest_cell_;

        /**
         * \brief True if w did not change, thus there is 
         *  no need to recompute the power diagram.
         */
        bool w_did_not_change_;

        /** \brief vertex to tet corner mapping */
        vector<index_t> v_to_c_;

        /** \brief linked corners incident to same vrtx */
        vector<index_t> nxt_c_around_v_;

        /** \brief CRS matrix - row pointer */
        vector<index_t> rowptr_;
        
        /** \brief CRS matrix - column index */
        vector<index_t> colind_;

        /** \brief CRS matrix - value array */
        vector<double>  val_;

        /** \brief CRS matrix - CRS storage */
        bool symmetric_storage_;
        
        /** \brief right-hand side of the linear system */
        vector<double>  rhs_;
    };


    /**
     * \brief Computes a shape that interpolates the two input tet
     *  meshes.
     * \details The shape is composed of tetrahedra
     * \param [in] CVT the Centroidal Voronoi Tesselation
     *   used to sample the second shape
     * \param [in] OTM the Optimal Transport Map
     * \param [out] morph mesh where to store the morphing shape. It uses
     *   6d coordinates (original location + final location).
     * \param[in] filter_tets if true, remove the tetrahedra that are outside
     *   the source mesh.
     */
    void EXPLORAGRAM_API compute_morph(
        CentroidalVoronoiTesselation& CVT,
        OptimalTransportMap& OTM,
        Mesh& morph,
        bool filter_tets=true
    );


    /**
     * \brief Computes the surface that corresponds to discontinuities
     *  in the optimal transport map.
     * \details The surface is determined as the facets of Voronoi cells
     *  that are adjacent in Pow(X)|M1 but not in Vor(X)|M2 
     * \param [in] CVT the Centroidal Voronoi Tesselation
     *   used to sample the second shape M2
     * \param [in] OTM the Optimal Transport Map with the
     *   power diagram that samples the first shape M1
     * \param [out] singular_set where to store the singular surface
     */
    void EXPLORAGRAM_API compute_singular_surface(        
        CentroidalVoronoiTesselation& CVT,
        OptimalTransportMap& OTM,
        Mesh& singular_set
    );

}

#endif
