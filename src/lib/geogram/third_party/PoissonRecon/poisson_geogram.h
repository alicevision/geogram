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
 

#ifndef H__OGF_SHAPENET_ALGOS_POISSON__H
#define H__OGF_SHAPENET_ALGOS_POISSON__H

#include <geogram/basic/common.h>
#include <geogram/basic/numeric.h>

// Note: This file is not part of Kahzdan's original Poisson reconstruction
// code.

namespace GEO {

    class Mesh;
    
    /**
     * \brief A wrapper to interface Kahzdan's Poisson reconstruction
     *  with Geogram datastructures.
     */
    class GEOGRAM_API PoissonReconstruction {
    public:
        
        /**
         * \brief PoissonReconstruction constructor.
         */
        PoissonReconstruction();

        /**
         * \brief Reconstructs a surface.
         * \param[in] points a pointer to a Mesh with the points.
         *  It needs to have a vector attribute of dimension 3 
         *  called "normal" and attached to the vertices.
         * \param[out] surface the reconstructed surface
         */
        void reconstruct(Mesh* points, Mesh* surface);

        /**
         * \brief Sets the depth of the octree.
         * \param[in] x the new depth of the octree
         * \details Default value is 8. Use 10 or 11 for highly 
         *  detailed models.
         */
        void set_depth(index_t x) {
            depth_ = int(x);
        }

        /**
         * \brief Gets the depth of the octree.
         * \return the depth of the octree
         */
        index_t get_depth() const {
            return index_t(depth_);
        }
        
        
    private:

        bool performance_;
        bool complete_;
        bool no_comments_;
        bool polygon_mesh_;
        bool confidence_;
        bool normal_weights_;
        bool non_manifold_;
        bool dirichlet_;
        bool ascii_;
        bool density_;
        bool linear_fit_;
        bool primal_voxel_;
        bool verbose_;
        
        int degree_;
        int depth_;
        int cg_depth_;
        int kernel_depth_;
        int adaptive_exponent_;
        int iters_;
        int voxel_depth_;
        int full_depth_;
        int min_depth_;
        int max_solve_depth_;
        int threads_;

        float color_;
        float samples_per_node_;
        float scale_;
        float cg_accuracy_;
        float point_weight_;
    };
    
}

#endif
