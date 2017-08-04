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

#include <geogram/basic/common.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/stopwatch.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/process.h>
#include <geogram/basic/progress.h>
#include <geogram/basic/matrix.h>
#include <geogram/basic/permutation.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/mesh/mesh_tetrahedralize.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/voronoi/CVT.h>
#include <geogram/voronoi/RVD.h>
#include <geogram/voronoi/generic_RVD_vertex.h>
#include <geogram/delaunay/delaunay.h>
#include <geogram/delaunay/delaunay_3d.h>
#include <geogram/points/nn_search.h>
#include <geogram/numerics/optimizer.h>
#include <stack>
#include <iterator>

#ifdef GEOGRAM_WITH_VORPALINE
#include <vorpalib/voronoi/LpCVT.h>
#define CentroidalVoronoiTesselation LpCentroidalVoronoiTesselation
#endif

/*
 * CLANG complains about some functions that are unused.
 */
#ifdef __clang__
#pragma GCC diagnostic ignored "-Wunused-member-function"
#endif

namespace {
    using namespace GEO;

    const char* banner[] = {
        " __        __                   _      _           \n",
        " \\ \\      / /_ _ _ __ _ __   __| |_ __(_)_   _____ \n",
        "  \\ \\ /\\ / / _` | '__| '_ \\ / _` | '__| \\ \\ / / _ \\\n", 
        "   \\ V  V / (_| | |  | |_) | (_| | |  | |\\ V /  __/\n",
        "    \\_/\\_/ \\__,_|_|  | .__/ \\__,_|_|  |_| \\_/ \\___|\n",
        "                     |_|                           \n",
        "\n",
        nil
    };

    /**
     * \brief Computes the linear least squares
     * regression of a function evaluated
     * in 3d.
     */

    // TODO: have a linear solve function that does
    // not require a template argument...

    class LinearLeastSquares {
    public:
        /**
         * \brief Constructs a new LinearLeastSquares
         * \param[in] degree one of 1 (linear), 2 (quadratic)
         */
        LinearLeastSquares(
            index_t degree
        ) :
            degree_(degree)
        {
            switch(degree_) {
                case 1:
                    dim_ = 4;
                    break;
                case 2:
                    dim_ = 10;
                    break;
                default:
                    geo_assert_not_reached;
            }
        }

        /**
         * \brief Starts a new computation.
         */
        void begin() {
            AtA_4_.load_zero();
            AtA_10_.load_zero();
            for(index_t i = 0; i < MAX_DIM; ++i) {
                Atb_[i] = 0.0;
            }
        }

        /**
         * \brief Ends the current computation.
         * \details Computes the current equation
         *  from the set of samples declared with
         *  add_point().
         */
        void end() {
            switch(degree_) {
                case 1:
                {
                    Matrix<4,double> M = AtA_4_.inverse();
                    mult(M, Atb_, eqn_);
                } break;
                case 2:
                {
                    Matrix<10,double> M = AtA_10_.inverse();
                    mult(M, Atb_, eqn_);
                } break;
                default:
                    geo_assert_not_reached;
            }
        }

        /**
         * \brief Adds a sample to the current computation.
         * \details This function needs to be called between
         *  a begin() / end() pair.
         * \param[in] p 3d coordinates of the point
         * \param[in] v function value associated with \p p_in
         */
        void add_point(const double* p, double v) {
            double b[MAX_DIM];
            eval_basis(p, b);

            for(index_t i = 0; i < dim(); ++i) {
                for(index_t j = 0; j < dim(); ++j) {
                    switch(degree_) {
                    case 1: {
                        AtA_4_(i, j) += b[i] * b[j];
                    } break;
                    case 2: {
                        AtA_10_(i, j) += b[i] * b[j];
                    } break;
                    default: 
                        geo_assert_not_reached;
                    }
                }
                Atb_[i] += b[i] * v;
            }
        }

        /**
         * \brief Evaluates the least-squares linear estimate
         *  at a given point.
         * \details This function beeds to be called after end().
         * \param[in] p 3d coordinates of the point
         * \return the linear estimate at \p p
         */
        double eval(const double* p) const {
            double b[MAX_DIM];
            for(index_t i = 0; i < MAX_DIM; ++i) {
                b[i] = 0.0;
            }
            eval_basis(p, b);
            double result = 0;
            for(index_t i = 0; i < dim(); ++i) {
                result += eqn_[i] * b[i];
            }
            return result;
        }

    protected:
        /**
         * \brief Gets the dimension of the function basis.
         */
        index_t dim() const {
            return dim_;
        }

        /**
         * \brief Evaluates the function basis at a given
         *  point.
         * \param[in] p 3d coordinates of the point
         * \param[out] b array of size dim(), value of the
         *  function basis at \p p
         */
        void eval_basis(const double* p, double* b) const {
            double x = p[0];
            double y = p[1];
            double z = p[2];
            b[0] = 1.0;
            b[1] = x;
            b[2] = y;
            b[3] = z;
            if(degree_ >= 2) {
                b[4] = x * x;
                b[5] = y * y;
                b[6] = z * z;
                b[7] = x * y;
                b[8] = y * z;
                b[9] = z * x;
            }
        }

        /**
         * \brief Maximum dimension of the function basis
         */
        static const int MAX_DIM = 10;

    private:
        index_t degree_;
        index_t dim_;
        Matrix<4,double> AtA_4_;
        Matrix<10,double> AtA_10_;
        double Atb_[MAX_DIM];
        double eqn_[MAX_DIM];
    };

    /*************************************************************************/

    /**
     * \brief Computes the contribution of an integration simplex
     *  to the objective function minimized by a semi-discrete
     *  optimal transport map.
     */
    class OTMIntegrationSimplex : public IntegrationSimplex {
    public:

        /**
         * \brief Constructs a new OTMIntegrationSimplex.
         * \param[in] M the input mesh
         */
        OTMIntegrationSimplex(
            const Mesh& M
        ) : IntegrationSimplex(
            M, true, 0, 0, nil
        ), w_(nil) {
            weighted_ = M.vertices.attributes().is_defined("weight");
            varying_background_ = weighted_;
        }
        
        /**
         * \brief Sets the weight vector
         * \param[in] w a const pointer to the weight vector.
         */
        void set_w(const double* w) { 
            w_ = w ; 
        }

        /**
         * \copydoc IntegrationSimplex::eval()
         */
        virtual double eval(
            index_t center_vertex_index,
            const GEOGen::Vertex& v1,
            const GEOGen::Vertex& v2,
            const GEOGen::Vertex& v3,
            index_t t,
            index_t t_adj,
            index_t v_adj
        ) {
            geo_argused(t_adj);
            geo_argused(v_adj);
            double fT = 0.0;
            if(weighted_) {
                if(
                    current_tet_ != t ||
                    current_seed_ != center_vertex_index
                ) {
                    current_tet_ = t;
                    current_seed_ = center_vertex_index;
                    p0_[0] = v1.point()[0];
                    p0_[1] = v1.point()[1];
                    p0_[2] = v1.point()[2];
                    p0_mass_ = v1.weight();
                } else {
                    fT = eval_with_density(
                        center_vertex_index,
                        p0_, p0_mass_,
                        v1.point(), v1.weight(),
                        v2.point(), v2.weight(),
                        v3.point(), v3.weight()
                    );
                }
            } else {
                const double* p0 = point(center_vertex_index);
                const double* p1 = v1.point();
                const double* p2 = v2.point();
                const double* p3 = v3.point();
                double m = Geom::tetra_signed_volume(p0, p1, p2, p3);

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
                
                //  Spinlocks are used in multithreading mode, to avoid
                // that two threads update g_[center_vertex_index]
                // simultaneously.
                if(spinlocks_ != nil) {
                    spinlocks_->acquire_spinlock(center_vertex_index);
                }
                // -m because we maximize F <=> minimize -F                
                g_[center_vertex_index] -= m;
                if(spinlocks_ != nil) {
                    spinlocks_->release_spinlock(center_vertex_index);
                }
            }
            // -fT because we maximize F <=> minimize -F
            return -fT;
        }

        /**
         * \brief Evaluates the objective function maximized by
         *  OTM on an integration simplex when the background 
         *  mesh has a varying density interpolated on the vertices.
         * \param[in] center_vertex_index index of the first vertex
         *  of the integration simplex, that corresponds to one of 
         *  the points to be optimized.
         * \param[in] p0 const pointer to the three coordinates of
         *   the first vertex
         * \param[in] p0_mass the mass of the first vertex
         * \param[in] p1 const pointer to the three coordinates of
         *   the second vertex
         * \param[in] p1_mass the mass of the second vertex
         * \param[in] p2 const pointer to the three coordinates of
         *   the third vertex
         * \param[in] p2_mass the mass of the third vertex
         * \param[in] p3 const pointer to the three coordinates of
         *   the fourth vertex
         * \param[in] p3_mass the mass of the fourth vertex
         */
        double eval_with_density(
            index_t center_vertex_index,
            const double* p0, double p0_mass,
            const double* p1, double p1_mass,
            const double* p2, double p2_mass,
            const double* p3, double p3_mass            
        ) {
            const double* q = point(center_vertex_index);

            double Tvol = Geom::tetra_volume<3>(p0,p1,p2,p3);
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
            fT += (alpha[0] + rho[0]) * dotprod_00;  // 0 0
            fT += (alpha[1] + rho[0]) * dotprod_10;  // 1 0
            fT += (alpha[1] + rho[1]) * dotprod_11;  // 1 1
            fT += (alpha[2] + rho[0]) * dotprod_20;  // 2 0
            fT += (alpha[2] + rho[1]) * dotprod_21;  // 2 1
            fT += (alpha[2] + rho[2]) * dotprod_22;  // 2 2
            fT += (alpha[3] + rho[0]) * dotprod_30;  // 3 0
            fT += (alpha[3] + rho[1]) * dotprod_31;  // 3 1
            fT += (alpha[3] + rho[2]) * dotprod_32;  // 3 2
            fT += (alpha[3] + rho[3]) * dotprod_33;  // 3 3            

            fT = Tvol * fT / 60.0 - m * w_[center_vertex_index];
            
            //  Spinlocks are used in multithreading mode, to avoid
            // that two threads update g_[center_vertex_index]
            // simultaneously.
            if(spinlocks_ != nil) {
                spinlocks_->acquire_spinlock(center_vertex_index);
            }
            // -m because we maximize F <=> minimize -F                
            g_[center_vertex_index] -= m;
            if(spinlocks_ != nil) {
                spinlocks_->release_spinlock(center_vertex_index);
            }

            // Note Ft is negated (minimize -F) by caller
            return fT;
        }
        
        virtual void reset_thread_local_storage() {
            current_tet_ = index_t(-1);
            current_seed_ = index_t(-1);
        }
        
    private:
        const double* w_;
        bool weighted_;

        //   For each Voronoi Cell / background tet intersection,
        // we need to tessellate the corresponding polyhedron.
        //   The callback traverses the border of the polyhedron,
        // one triangle at a time. To tesselate the polyhedron, we
        // need to keep track of the first vertex. This needs to
        // be a thread-local variable.

        static GEO_THREAD_LOCAL double p0_[3];
        static GEO_THREAD_LOCAL double p0_mass_;
        static GEO_THREAD_LOCAL index_t current_seed_;
        static GEO_THREAD_LOCAL index_t current_tet_;
    } ;

    GEO_THREAD_LOCAL double OTMIntegrationSimplex::p0_[3];
    GEO_THREAD_LOCAL double OTMIntegrationSimplex::p0_mass_;    
    GEO_THREAD_LOCAL index_t OTMIntegrationSimplex::current_tet_;
    GEO_THREAD_LOCAL index_t OTMIntegrationSimplex::current_seed_;        
    
    /*************************************************************************/

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
    class OptimalTransportMap {
    public:
        /**
         * \brief Initializes a new OptimalTransportMap.
         * \param[in] mesh the source distribution, represented as a 3d mesh
         * \param[in] delaunay factory name of the Delaunay triangulation
         */
        OptimalTransportMap(
            Mesh* mesh,
            const std::string& delaunay = "default"
        ) :
            mesh_(mesh)
        {
            // Note: we represent power diagrams as 4d Voronoi diagrams
            delaunay_ = Delaunay::create(4, delaunay);
            RVD_ = RestrictedVoronoiDiagram::create(delaunay_, mesh_);
            RVD_->set_volumetric(true);
            RVD_->set_check_SR(true);
            RVD_->create_threads();

            //   No need to reorder vertices if BRIO is activated since
            // vertices are then already reordered.
            if(CmdLine::get_arg_bool("BRIO")) {
                RVD_->delaunay()->set_reorder(false);
            }

            instance_ = nil;
            lambda_p_ = 0.0;
            total_mass_ = 0.0;
            current_call_iter_ = 0;
            epsilon_ = 0.01;
            level_ = 0;
            
            simplex_func_ = new OTMIntegrationSimplex(*mesh);

            save_RVD_iter_ = CmdLine::get_arg_bool("RVD_iter");
            current_iter_ = 0;

            pretty_log_ = CmdLine::get_arg_bool("log:pretty");
        }

        /**
         * \brief Sets the points that define the target distribution.
         * \param[in] nb_points number of points in the target distribution
         * \param[in] points coordinates of the Diracs centers in the target
         *  distribution.
         */
        void set_points(index_t nb_points, const double* points) {
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
                double tet_mass = Geom::tetra_volume<3>(
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

        /**
         * \brief Sets the maximum error.
         * \param eps acceptable relative deviation for the measure of a
         *   Voronoi cell.
         */
        void set_epsilon(double eps) {
            epsilon_ = eps;
        }

        /**
         * \brief Computes the weights that realize the optimal
         *  transport map between the source mesh and the target
         *  pointset.
         * \param[in] max_iterations maximum number of solver iterations.
         */
        void optimize(index_t max_iterations) {
            level_ = 0;
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
            optimizer->optimize(&weights_[0]);
            instance_ = nil;
            // To make sure everything is reset properly
            double dummy = 0;
            funcgrad(n, &weights_[0], dummy, nil);
            Logger::out("OTM")
                << "Used " << current_call_iter_ << " iterations" << std::endl;

        }

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
        void optimize_level(index_t b, index_t e, index_t max_iterations) {

            // If this is not the first level, propagate the weights from
            // the lower levels.
            if(b != 0) {

                //   Create a nearest neighbor search data structure
                // and insert the [0..b) samples into it (they were
                // initialized at previous calls).
                NearestNeighborSearch_var NN = NearestNeighborSearch::create(3);
                NN->set_points(b, &points_4d_[0], 4);
                index_t degree = CmdLine::get_arg_uint("fitting_degree");

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
            optimizer->optimize(&weights_[0]);
            instance_ = nil;

            // To make sure everything is reset properly
            double dummy = 0;
            funcgrad(n, &weights_[0], dummy, nil);
        }

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
        }

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
         * \brief Gets a weight.
         * \param[in] i index of the point
         * \return the weight that was computed for point \p i
         */
        double weight(index_t i) const {
            return weights_[i];
        }

        /**
         * \brief Gets the value of Kantorowich potential 
         *  at a given point.
         * \param[in] i index of the point
         * \return the potential that was computed for point \p i
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
        ) {
            instance_->funcgrad(n, x, f, g);
        }

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
        ) {
            geo_argused(n);
            geo_argused(x);
            geo_argused(f);
            geo_argused(g);
            geo_argused(gnorm);
            instance_->newiteration();
        }

        /**
         * \brief Gets the restricted Voronoi diagram.
         * \return a pointer to the restricted Voronoi diagram
         */
        RestrictedVoronoiDiagram* RVD() {
            return RVD_;
        }

    protected:

        /**
         * \brief Callback for the numerical solver.
         */
        void newiteration() {
            if(save_RVD_iter_) {
                save_RVD(current_iter_);
            }
            ++current_iter_;
        }
        
        /**
         * \brief Saves the RVD at each iteration if
         *   specified on command line (just for debugging/
         *   explaining the algorithm).
         * \param[in] id index to be used for the file, that
         *   will be named RVD_id.meshb
         */
        void save_RVD(index_t id) {
            
            if(!save_RVD_iter_) {
                return;
            }
            
            Mesh RVD_mesh;
            Attribute<index_t> tet_region(RVD_mesh.cells.attributes(),"region");
            RVD()->compute_RVD(
                RVD_mesh,
                0,     // dim (0 means use default)
                CmdLine::get_arg_bool("RVD:borders_only"),
                CmdLine::get_arg_bool("RVD:integration_simplices")
            );
            RVD_mesh.vertices.set_dimension(3);
            RVD_mesh.cells.connect();
            MeshIOFlags flags;
            flags.set_attribute(MESH_CELL_REGION);
            flags.set_attribute(MESH_FACET_REGION);            
            mesh_save(
                RVD_mesh,
                "RVD_" + String::to_string(id) + ".meshb",
                flags
            );
        }

        /**
         * \brief Computes the objective function and its gradient.
         * \param[in] n number of variables
         * \param[in] w current value of the variables
         * \param[out] f current value of the objective function
         * \param[out] g gradient of the objective function
         */
        void funcgrad(index_t n, double* w, double& f, double* g) {

            // Step 1: determine the 4d embedding from the weights
            double W = 0.0;
            for(index_t p = 0; p < n; ++p) {
                W = geo_max(W, w[p]);
            }
            for(index_t p = 0; p < n; ++p) {
                points_4d_[4 * p + 3] = ::sqrt(W - w[p]);
            }

            // Step 2: compute function and gradient 
            delaunay_->set_vertices(n, &points_4d_[0]);

            if(g == nil) {
                if(pretty_log_) {
                    CmdLine::ui_clear_line();
                    CmdLine::ui_message(last_stats_ + "\n");
                }
                return;
            }

            index_t nb_empty_cells = 0;

            f = 0.0;
            for(index_t p = 0; p < n; ++p) {
                g[p] = 0.0;
            }

            OTMIntegrationSimplex* simplex_func_otm = 
                dynamic_cast<OTMIntegrationSimplex*>(simplex_func_.get());
            simplex_func_otm->set_w(w);
            RVD_->compute_integration_simplex_func_grad(f, g, simplex_func_otm);

            double max_diff = 0.0;
            double avg_diff = 0.0;
            for(index_t p = 0; p < n; ++p) {
                f += lambda_p_ * w[p];

                // Note: we minimize -f instead of maximizing f,
                // therefore, in the paper:
                //    g[p] = lambda_p - mesure(power cell associated with p)
                //
                // What is programmed:
                //    g[p] = mesure(power cell associated with p) - lambda_p

                if(::fabs(g[p]) < 1e-10) {
                    nb_empty_cells++;
                }

                g[p] = -g[p] - lambda_p_;

                double cur_diff = ::fabs(g[p]);
                max_diff = geo_max(max_diff, cur_diff);
                avg_diff += cur_diff / double(n);
            }

            double gNorm = 0.0;
            for(index_t i = 0; i < n; ++i) {
                gNorm += geo_sqr(g[i]);
            }
            gNorm = ::sqrt(gNorm);

            std::ostringstream str;
            if(level_ == 0) {
                str << "o-[OTM         ] " ;
            } else {
                str << "o-[OTM Lvl." << level_ << "   ] " ;
            }
            str << "iter=" << current_call_iter_
                << " nbZ=" << nb_empty_cells
                //                << " f=" << f
                //                << " avg_diff=" << avg_diff
                //                << " max_diff=" << max_diff
                << " g=" << gNorm
                << " threshold=" << gradient_threshold(n);
            last_stats_ = str.str();

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
        index_t current_iter_;
    };

    OptimalTransportMap* OptimalTransportMap::instance_ = nil;
    
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
        OptimalTransportMap& OTM,
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
            if(CmdLine::get_arg_bool("RVD")) {
                MeshIOFlags flags;
                flags.set_element(MESH_CELLS);
                flags.set_attribute(MESH_CELL_REGION);
                mesh_save(RVD, "RVD.meshb", flags);
            }

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
                double mt = Geom::tetra_signed_volume(p0, p1, p2, p3);
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
    
    /**
     * \brief Computes a shape that interpolates the two input tet
     *  meshes.
     * \details The shape is composed of tetrahedra
     * \param [in] CVT the Centroidal Voronoi Tesselation
     *   used to sample the second shape
     * \param [in] OTM the Optimal Transport Map
     * \param [in] filename where to store the morphing shape
     *  (.tet6 tet mesh 6d coordinates)
     */
    void compute_morph(
        CentroidalVoronoiTesselation& CVT,
        OptimalTransportMap& OTM,
        const std::string& filename
    ) {
        Logger::out("OTM")
            << "Computing coherent tet mesh and saving result to:"
            << filename << std::endl;

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
        Mesh morph;

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
        
        std::string extension = FileSystem::extension(filename);
        if(extension == "obj6") {
            mesh_save(morph, filename);            
        } else if(extension == "tet6") {
            MeshIOFlags flags;
            flags.set_element(MESH_CELLS);
            mesh_save(morph, filename, flags);
        } else if(extension == "eobj") {
            mesh_repair(morph);
            std::ofstream out(filename.c_str());
            if(!out) {
                Logger::err("OTM")
                    << filename
                    << ": cannot create file"
                    << std::endl;
                return;
            }
            out << "# attribute geom2 vertex vec3" << std::endl;
            out << "# attribute potential vertex real" << std::endl;
            for(index_t v=0; v<morph.vertices.nb(); ++v) {
                out << "v " << M2_vertices[v] << std::endl;
                out << "# attrs v " << v+1 << " " << M1_vertices[v] << " "
                    << OTM.potential(v)
                    << std::endl;
            }
            for(index_t f=0; f<morph.facets.nb(); ++f) {
                out << "f ";
                for(
                    index_t c=morph.facets.corners_begin(f);
                    c<morph.facets.corners_end(f); ++c) {
                    out << morph.facet_corners.vertex(c)+1 << " ";
                }
                out << std::endl;
            }
        } else {
            Logger::warn("OTM")
                << filename
                << "Unknown extension for transport"
                << std::endl;
            MeshIOFlags flags;
            flags.set_element(MESH_CELLS);
            mesh_save(morph, filename, flags);
        }
        
    }

    /**
     * \brief Computes the surface that corresponds to discontinuities
     *  in the optimal transport map.
     * \details The surface is determined as the facets of Voronoi cells
     *  that are adjacent in Pow(X)|M1 but not in Vor(X)|M2 
     * \param [in] CVT the Centroidal Voronoi Tesselation
     *   used to sample the second shape M2
     * \param [in] OTM the Optimal Transport Map with the
     *   power diagram that samples the first shape M1
     * \param [in] filename where to store the singular surface
     *  (Graphite .obj file format)
     */
    void compute_singular_surface(        
        CentroidalVoronoiTesselation& CVT,
        OptimalTransportMap& OTM,
        const std::string& filename
    ) {
        std::ofstream out(filename.c_str());
        if(!out) {
            Logger::err("Singular") 
                << filename << ":could not create file" 
                << std::endl;
            return;
        } else {
            Logger::out("Singular") 
                << "saving singular surface to:" << filename 
                << std::endl;
        }

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

        Mesh singular;
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

        mesh_repair(singular);

        mesh_save(singular,filename);
        
    }

    /**
     * \brief Internal implementation function for
     *   compute_hierarchical_sampling().
     * \param[in,out] CVT the CentroidalVoronoiTesselation, initialized
     *  with the volume to be sampled. On output, it stores the samples
     * \param[in] nb_samples total number of samples to generate
     * \param[out] levels sample indices that correspond to level l are
     *   in the range levels[l] (included) ... levels[l+1] (excluded)
     * \param[in] ratio number of samples ratio between two consecutive
     *   levels
     * \param[in] threshold minimum number of samples in a level
     * \param[in] b first element of the level to be generated
     * \param[in] e one position past the last element of the
     *  level to be generated
     * \param[in,out] points work vector allocated by caller,
     *  of size 3*nb_samples
     */
    void compute_hierarchical_sampling_recursive(
        CentroidalVoronoiTesselation& CVT,
        index_t nb_samples,
        vector<index_t>& levels,
        double ratio,
        index_t threshold,
        index_t b, index_t e,
        vector<double>& points
    ) {
        index_t m = b;

        // Recurse in [b...m) range
        if(e - b > threshold) {
            m = b + index_t(double(e - b) * ratio);
            compute_hierarchical_sampling_recursive(
                CVT, nb_samples, levels, ratio, threshold, b, m, points
            );
        }

        // Initialize random points in [m...e) range
        CVT.RVD()->compute_initial_sampling(&points[3 * m], e - m);

        //  Set the points in [b...e) range
        CVT.set_points(e - b, &points[0]);

        // Lock [b...m) range
        if(CmdLine::get_arg_bool("lock")) {
            for(index_t i = b; i < e; ++i) {
                if(i < m) {
                    CVT.lock_point(i);
                } else {
                    CVT.unlock_point(i);
                }
            }
        }

        Logger::div(
            std::string("Generating level ") +
            String::to_string(levels.size())
        );

        Logger::out("Sample") << " generating a level with " << e - m
            << " samples" << std::endl;

        try {
            ProgressTask progress("Lloyd", 100);
            CVT.set_progress_logger(&progress);
            CVT.Lloyd_iterations(CmdLine::get_arg_uint("opt:nb_Lloyd_iter"));
        }
        catch(const TaskCanceled&) {
        }

        try {
            ProgressTask progress("Newton", 100);
            CVT.set_progress_logger(&progress);
            CVT.Newton_iterations(CmdLine::get_arg_uint("opt:nb_Newton_iter"));
        }
        catch(const TaskCanceled&) {
        }

        levels.push_back(e);
    }

    /**
     * \brief Computes a hierarchical sampling of a volume.
     * \param[in,out] CVT the CentroidalVoronoiTesselation, initialized
     *  with the volume to be sampled. On output, it stores the samples
     * \param[in] nb_samples total number of samples to generate
     * \param[out] levels sample indices that correspond to level l are
     *   in the range levels[l] (included) ... levels[l+1] (excluded)
     * \param[in] ratio number of samples ratio between two consecutive
     *   levels
     * \param[in] threshold minimum number of samples in a level
     */
    void compute_hierarchical_sampling(
        CentroidalVoronoiTesselation& CVT,
        index_t nb_samples,
        vector<index_t>& levels,
        double ratio = 0.125,
        index_t threshold = 300
    ) {
        levels.push_back(0);
        vector<double> points(nb_samples * 3);
        compute_hierarchical_sampling_recursive(
            CVT, nb_samples, levels, ratio, threshold,
            0, nb_samples,
            points
        );
        CVT.unlock_all_points();
    }

    /**
     * \brief Computes a sampling of a volume.
     * \param[in,out] CVT the CentroidalVoronoiTesselation, initialized
     *  with the volume to be sampled. On output, it stores the samples
     * \param[in] nb_samples total number of samples to generate
     */
    void compute_single_level_sampling(
        CentroidalVoronoiTesselation& CVT,
        index_t nb_samples
    ) {

        CVT.compute_initial_sampling(nb_samples);

        try {
            ProgressTask progress("Lloyd", 100);
            CVT.set_progress_logger(&progress);
            CVT.Lloyd_iterations(CmdLine::get_arg_uint("opt:nb_Lloyd_iter"));
        }
        catch(const TaskCanceled&) {
        }

        try {
            ProgressTask progress("Newton", 100);
            CVT.set_progress_logger(&progress);
            CVT.Newton_iterations(CmdLine::get_arg_uint("opt:nb_Newton_iter"));
        }
        catch(const TaskCanceled&) {
        }
    }

    /**
     * \brief Projects the points of a volumetric sampling
     *  onto the border of the volume.
     */
    void project_sampling_on_border(
        CentroidalVoronoiTesselation& CVT
    ) {
        try {
            ProgressTask progress("Surf. Lloyd", 100);
            CVT.set_progress_logger(&progress);
            CVT.set_volumetric(false);
            CVT.Lloyd_iterations(
                CmdLine::get_arg_uint("opt:nb_Lloyd_iter") * 2
            );
        }
        catch(const TaskCanceled&) {
        }

        if(CmdLine::get_arg_bool("feature_sensitive")) {

#ifdef GEOGRAM_WITH_VORPALINE            
            try {
                ProgressTask progress("LpCVT", 100);
                CVT.set_progress_logger(&progress);
                CVT.set_normal_anisotropy(5.0);
                CVT.Newton_iterations(30, 7);
            }
            catch(const TaskCanceled&) {
            }
            CVT.set_normal_anisotropy(1.0);
#endif
        }

        vector<double> mg(3 * CVT.nb_points());
        vector<double> m(CVT.nb_points());
        CVT.RVD()->compute_centroids(&mg[0], &m[0]);
        for(index_t i = 0; i < CVT.nb_points(); ++i) {
            if(m[i] == 0.0) {
                CVT.unlock_point(i);
            } else {
                CVT.lock_point(i);
            }
        }

        CVT.set_volumetric(true);

        try {
            ProgressTask progress("Relax. vol.", 100);
            CVT.set_progress_logger(&progress);
            CVT.Lloyd_iterations(
                CmdLine::get_arg_uint("opt:nb_Lloyd_iter") * 2
            );
        }
        catch(const TaskCanceled&) {
        }
    }

    /**
     * \brief Reorders the points in a Centroidal Voronoi Tesselation
     *  in such a way that continguous index ranges correspond to
     *  multiple resolutions.
     * \param[in,out] CVT the CentroidalVoronoiTesselation
     * \param[out] levels sample indices that correspond to level l are
     *   in the range levels[l] (included) ... levels[l+1] (excluded)
     * \param[in] ratio number of samples ratio between two consecutive
     *   levels
     * \param[in] threshold minimum number of samples in a level
     */
    void BRIO_reorder(
        CentroidalVoronoiTesselation& CVT,
        vector<index_t>& levels,
        double ratio,
        index_t threshold
    ) {
        vector<index_t> sorted_indices;
        compute_BRIO_order(
            CVT.nb_points(), CVT.embedding(0), sorted_indices,
            CVT.dimension(), threshold, ratio, &levels
        );
        Permutation::apply(
            CVT.embedding(0), sorted_indices, 
            index_t(CVT.dimension() * sizeof(double))
        );
    }

    /**
     * \brief Translates a mesh in such a way that its center matches
     *  the center of another mesh.
     * \param[in] M1 a const reference to the reference mesh
     * \param[in,out] M2 a reference to the mesh that will be recentered
     */
    void recenter_mesh(const Mesh& M1, Mesh& M2) {
        double xyzmin1[3];
        double xyzmax1[3];
        double xyzmin2[3];
        double xyzmax2[3];
        double xlat[3];
        get_bbox(M1, xyzmin1, xyzmax1);
        get_bbox(M2, xyzmin2, xyzmax2);
        for(coord_index_t c=0; c<3; ++c) {
            xlat[c] = 0.5*
                ((xyzmin1[c] + xyzmax1[c]) - (xyzmin2[c] + xyzmax2[c]));
        }
        for(index_t v=0; v<M2.vertices.nb(); ++v) {
            for(coord_index_t c=0; c<3; ++c) {
                M2.vertices.point_ptr(v)[c] += xlat[c];
            }
        }
    }


    /**
     * \brief Computes the volume of a tetrahedral mesh.
     * \param[in] M a const reference to the mesh
     * \return the volume of the tetrahedra of M
     */
    double mesh_tets_volume(const Mesh& M) {
        double result = 0.0;
        for(index_t t = 0; t < M.cells.nb(); ++t) {
            result += Geom::tetra_volume<3>(
                M.vertices.point_ptr(M.cells.tet_vertex(t, 0)),
                M.vertices.point_ptr(M.cells.tet_vertex(t, 1)),
                M.vertices.point_ptr(M.cells.tet_vertex(t, 2)),
                M.vertices.point_ptr(M.cells.tet_vertex(t, 3))
            );
        }
        return result;
    }

    /**
     * \brief Rescales a mesh in such a way that its total volume
     *  matches the volume of a reference mesh.
     * \param[in] M1 a const reference to the reference mesh
     * \param[in,out] M2 a reference to the mesh that will be rescaled
     */
    void rescale_mesh(const Mesh& M1, Mesh& M2) {
        double xyzmin[3];
        double xyzmax[3];
        get_bbox(M2, xyzmin, xyzmax);
        double s = pow(mesh_tets_volume(M1)/mesh_tets_volume(M2), 1.0/3.0);
        for(unsigned int v=0; v<M2.vertices.nb(); ++v) {
            for(index_t c=0; c<3; ++c) {
                double gc = 0.5*(xyzmin[c]+xyzmax[c]);
                M2.vertices.point_ptr(v)[c] =
                    gc + s * (M2.vertices.point_ptr(v)[c] - gc);
            }
        }
    }


    /**
     * \brief Loads a volumetric mesh.
     * \details If the specified file contains a surface, try to
     *  tesselate it. If the surface has self-intersections, try to
     *  remove them.
     * \param[in] filename the name of the file
     * \param[out] M the mesh
     * \retval true if the file was successfully loaded
     * \retval false otherwise
     */
    bool load_volume_mesh(const std::string& filename, Mesh& M) {
        MeshIOFlags flags;
        flags.set_element(MESH_CELLS);
        flags.set_attribute(MESH_CELL_REGION);

        if(!mesh_load(filename, M, flags)) {
            return 1;
        }
        if(!M.cells.are_simplices()) {
            Logger::err("I/O") << "File "
                               << filename
                               << " should only have tetrahedra" << std::endl;
            return false;
        }
        if(M.cells.nb() == 0) {
            Logger::out("I/O") << "File "
                               << filename
                               << " does not contain a volume" << std::endl;
            Logger::out("I/O") << "Trying to tetrahedralize..." << std::endl;
            if(!mesh_tetrahedralize(M,true,false)) {
                return false;
            }
        }
        return true;
    }

    enum DensityFunction {
        DENSITY_X=0,
        DENSITY_Y=1,
        DENSITY_Z=2,
        DENSITY_R,
        DENSITY_SIN,
        DENSITY_DIST
    };
    
    void set_density(Mesh& M, double mass1, double mass2) {

        if(mass1 == mass2) {
            return;
        }

        std::string function_str =
            CmdLine::get_arg("density_function");

        bool minus = false;
        if(function_str.length() > 1 && function_str[0] == '-') {
            minus = true;
            function_str = function_str.substr(1,function_str.length()-1);
        }
        double density_pow = 1.0;
        {
            std::size_t found = function_str.find('^');
            if(found != std::string::npos) {
                std::string pow_str =
                    function_str.substr(found+1, function_str.length()-found-1);
                density_pow = String::to_double(pow_str);
                function_str = function_str.substr(0,found);
            }
        }

        Logger::out("OTM")
            << "Using density: "
            << (minus ? "-" : "+")
            << function_str << "^"
            << density_pow
            << " rescaled to ("
            << mass1 << "," << mass2
            << ")"
            << std::endl;
        
        DensityFunction function;
        if(function_str == "X") {
            function = DENSITY_X;
        } else if(function_str == "Y") {
            function = DENSITY_Y;            
        } else if(function_str == "Z") {
            function = DENSITY_Z;            
        } else if(function_str == "R") {
            function = DENSITY_R;            
        } else if(function_str == "sin") {
            function = DENSITY_SIN;
        } else if(function_str == "dist") {
            function = DENSITY_DIST;
        } else {
            Logger::err("OTM") << function_str << ": no such density function"
                               << std::endl;
            return;
        }
        
        Attribute<double> mass(M.vertices.attributes(),"weight");

        switch(function) {
        case DENSITY_X:
        case DENSITY_Y:
        case DENSITY_Z: {
            for(index_t v=0; v<M.vertices.nb(); ++v) {
                mass[v] = M.vertices.point_ptr(v)[index_t(function)];
            }
        } break;
        case DENSITY_R: {
            double xyz_min[3];
            double xyz_max[3];
            get_bbox(M, xyz_min, xyz_max);
            for(index_t v=0; v<M.vertices.nb(); ++v) {
                double r=0;
                const double* p = M.vertices.point_ptr(v);
                for(coord_index_t c=0; c<3; ++c) {
                    r += geo_sqr(p[c] - 0.5*(xyz_min[c] + xyz_max[c]));
                }
                r = ::sqrt(r);
                mass[v] = r;
            }
        } break;
        case DENSITY_SIN: {
            double xyz_min[3];
            double xyz_max[3];
            get_bbox(M, xyz_min, xyz_max);
            for(index_t v=0; v<M.vertices.nb(); ++v) {
                double f = 1.0;
                const double* p = M.vertices.point_ptr(v);
                for(coord_index_t c=0; c<3; ++c) {
                    double coord = (p[c] - xyz_min[c]) / (xyz_max[c] - xyz_min[c]);
                    f *= sin(coord *  M_PI * 2.0 * 2.0);
                }
                mass[v] = f;                
            }
        } break;
        case DENSITY_DIST: {
            std::string ref_filename = CmdLine::get_arg("density_distance_reference");
            if(ref_filename != "")  {
                Mesh reference;
                MeshIOFlags flags;
                flags.reset_element(MESH_CELLS);
                if(!mesh_load(ref_filename, reference, flags)) {
                    exit(1);
                }
                MeshFacetsAABB AABB(reference);
                for(index_t v=0; v<M.vertices.nb(); ++v) {
                    mass[v] = ::sqrt(AABB.squared_distance(vec3(M.vertices.point_ptr(v))));
                }
            } else {
                MeshFacetsAABB AABB(M);
                for(index_t v=0; v<M.vertices.nb(); ++v) {
                    mass[v] = ::sqrt(AABB.squared_distance(vec3(M.vertices.point_ptr(v))));
                }
            }
        } break;
        }

        // Compute min and max mass
        double mass_min = Numeric::max_float64();
        double mass_max = Numeric::min_float64();
        for(index_t v=0; v<M.vertices.nb(); ++v) {
            mass_min = geo_min(mass_min, mass[v]);
            mass_max = geo_max(mass_max, mass[v]);
        }

        // Normalize mass, apply power, and rescale to (mass1 - mass2)
        for(index_t v=0; v<M.vertices.nb(); ++v) {
            double f = (mass[v] - mass_min) / (mass_max - mass_min);
            if(minus) {
                f = 1.0 - f;
            }
            f = ::pow(f,density_pow);
            mass[v] = mass1 + f*(mass2 - mass1);
        }
    }
}

int main(int argc, char** argv) {
    using namespace GEO;

    GEO::initialize();

    try {
        
        std::vector<std::string> filenames;

        CmdLine::import_arg_group("standard");
        CmdLine::import_arg_group("algo");
        CmdLine::import_arg_group("opt");
        CmdLine::declare_arg("nb_pts", 1000, "number of points");
        CmdLine::declare_arg("nb_iter", 1000, "number of iterations for OTM");
        CmdLine::declare_arg("RDT", false, "save regular triangulation");
        CmdLine::declare_arg_group(
            "RVD", "RVD output options", CmdLine::ARG_ADVANCED
        );
        CmdLine::declare_arg("RVD", false, "save restricted Voronoi diagram");
        CmdLine::declare_arg(
            "RVD_iter", false, "save restricted Voronoi diagram at each iteration"
        );
        CmdLine::declare_arg(
            "RVD:borders_only", false, "save only border of RVD"
        );        
        CmdLine::declare_arg(
            "RVD:integration_simplices", true, "export RVD as integration simplices"
        );        
        
        CmdLine::declare_arg("multilevel", true, "use multilevel algorithm");
        CmdLine::declare_arg("BRIO", true, 
                             "use BRIO reordering to compute the levels"
        );
        CmdLine::declare_arg("ratio", 0.125, "ratio between levels");
        CmdLine::declare_arg(
            "epsilon", 0.01, "relative measure error in a cell"
        );
        CmdLine::declare_arg(
            "lock", true, "Lock lower levels when sampling shape"
        );
        CmdLine::declare_arg(
            "fitting_degree", 2, "degree for interpolating weights"
        );
        CmdLine::declare_arg(
            "project", true, "project sampling on border"
        );
        CmdLine::declare_arg(
            "feature_sensitive", true, "attempt to recover hard edges"
        );
        CmdLine::declare_arg(
            "singular", false, "compute and save singular surface"
        );
        CmdLine::set_arg("algo:delaunay", "BPOW");
        CmdLine::declare_arg(
            "recenter", true, "recenter target onto source mesh"
        );
        CmdLine::declare_arg(
            "rescale", true, "rescale target to match source volume"
        );
        CmdLine::declare_arg(
            "density_min", 1.0, "min density in first mesh"
        );
        CmdLine::declare_arg(
            "density_max", 1.0, "max density in first mesh"
        );
        CmdLine::declare_arg(
            "density_function", "x", "used function for density"
        );
        CmdLine::declare_arg(
            "density_distance_reference", "",
            "filename of the reference surface"
        );
        CmdLine::declare_arg(
            "out", "morph.tet6", "output filename"
        );
        
        Logger::div("Warpdrive - Optimal Transport");
        const char** banner_line = banner;
        while(*banner_line) {
            CmdLine::ui_message(*banner_line);
            banner_line++;
        }
        
        if(
            !CmdLine::parse(
                argc, argv, filenames, "mesh1 mesh2"
            )
        ) {
            return 1;
        }

        std::string mesh1_filename = filenames[0];
        std::string mesh2_filename = filenames[1];
        std::string output_filename = CmdLine::get_arg("out");
        if(filenames.size() == 3) {
            output_filename = filenames[2];
        }
        
        Logger::div("Loading data");

        Mesh M1;
        Mesh M2;
        Mesh M2_samples;
        
        if(!load_volume_mesh(mesh1_filename, M1)) {
            return 1;
        }
        
        if(!load_volume_mesh(mesh2_filename, M2)) {
            return 1;
        }

        set_density(
            M1,
            CmdLine::get_arg_double("density_min"),
            CmdLine::get_arg_double("density_max")
        );
        
        if(CmdLine::get_arg_bool("recenter")) {
            recenter_mesh(M1,M2);
        }

        if(CmdLine::get_arg_bool("rescale")) {
            rescale_mesh(M1,M2);
        }

        if(M1.cells.nb() == 0) {
            Logger::err("Mesh") << "M1 does not have any tetrahedron, exiting"
                << std::endl;
            return 1;
        }

        if(M2.cells.nb() == 0) {
            Logger::err("Mesh") << "M2 does not have any tetrahedron, exiting"
                << std::endl;
            return 1;
        }

        
        Logger::div("Sampling target shape");

        CentroidalVoronoiTesselation CVT(&M2, 0, "NN");
        vector<index_t> levels;
        CVT.set_volumetric(true);

        bool multilevel =
            CmdLine::get_arg_bool("multilevel") || 
            CmdLine::get_arg_bool("BRIO");

        if(CmdLine::get_arg_bool("RVD_iter") && multilevel) {
            Logger::warn("OTM") << "Deactivating multilevel mode" << std::endl;
            Logger::warn("OTM") << "(because RVD_iter is set)" << std::endl;            
            multilevel = false;
        }
        
        if(multilevel) {
            if(CmdLine::get_arg_bool("BRIO")) {
                compute_single_level_sampling(
                    CVT,
                    CmdLine::get_arg_uint("nb_pts")
                );
                BRIO_reorder(
                    CVT, levels, CmdLine::get_arg_double("ratio"), 300
                );
            } else {
                compute_hierarchical_sampling(
                    CVT,
                    CmdLine::get_arg_uint("nb_pts"),
                    levels,
                    CmdLine::get_arg_double("ratio")
                );
            }
        } else {
            compute_single_level_sampling(
                CVT,
                CmdLine::get_arg_uint("nb_pts")
            );
        }

        if(CmdLine::get_arg_bool("project")) {
            project_sampling_on_border(CVT);
        }

        M2_samples.vertices.assign_points(
            CVT.embedding(0), CVT.dimension(), CVT.nb_points()
        );
        
        Logger::div("Optimal transport");
        // Everything happens in dimension 4 (power diagram is seen
        // as Voronoi diagram in dimension 4), therefore the dimension
        // of M1 needs to be changed as well (even if it is not used).
        M1.vertices.set_dimension(4);
        OptimalTransportMap OTM(&M1);
        OTM.set_points(
            M2_samples.vertices.nb(), M2_samples.vertices.point_ptr(0)
        );
        OTM.set_epsilon(CmdLine::get_arg_double("epsilon"));
        index_t nb_iter = CmdLine::get_arg_uint("nb_iter");

        {
            Stopwatch W("OTM Total");
            if(multilevel) {
                OTM.optimize_levels(levels, nb_iter);
            } else {
                OTM.optimize(nb_iter);
            }
        }

        Logger::div("Morphing");
        Logger::out("OTM") <<  "Time-coherent triangulation." << std::endl;

        compute_morph(
            CVT, OTM, output_filename
        );
        

        if(CmdLine::get_arg_bool("singular")) {
            Logger::out("OTM") << "Computing singular set." << std::endl;
            compute_singular_surface(CVT,OTM,"singular.obj");
        }
    }
    catch(const std::exception& e) {
        std::cerr << "Received an exception: " << e.what() << std::endl;
        return 1;
    }

    Logger::out("") << "Everything OK, Returning status 0" << std::endl;
    return 0;
}

