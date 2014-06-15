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

#ifndef __GEOGRAM_POINTS_KD_TREE__
#define __GEOGRAM_POINTS_KD_TREE__

#include <geogram/basic.h>
#include <algorithm>

namespace GEO {

    /**
     * \brief A balanced kd-tree for nearest neighbor search.
     */
    class GEOGRAM_API KdTree {
    public:
        /**
         * \brief Creates a new KdTree.
         * \param[in] dim dimension of the points
         */
        KdTree(coord_index_t dim);

        /**
         * \brief KdTree destructor
         */
        ~KdTree();

        /**
         * \brief Sets the points and create the search data structure.
         * \param[in] nb_points number of points
         * \param[in] points an array of nb_points * dimension()
         */
        void set_points(index_t nb_points, const double* points);

        /**
         * \brief Sets the points and create the search data structure.
         * \details This variant has a stride parameter. 
         * \param[in] nb_points number of points
         * \param[in] points an array of nb_points * dimension()
         * \param[in] stride number of doubles between two consecutive
         *  points (stride=dimension() by default).
         */
        void set_points(
            index_t nb_points, const double* points, index_t stride
        );

        /**
         * \brief Finds the nearest neighbors of a point given by
         *  coordinates.
         * \param[in] nb_neighbors number of neighbors to be searched.
         *  Should be smaller or equal to nb_points() (else it triggers
         *  an assertion)
         * \param[in] query_point as an array of dimension() doubles
         * \param[out] neighbors array of nb_neighbors index_t
         * \param[out] neighbors_sq_dist array of nb_neighbors doubles
         */
        void get_nearest_neighbors(
            index_t nb_neighbors,
            const double* query_point,
            index_t* neighbors,
            double* neighbors_sq_dist
        ) const;

        /**
         * \brief Finds the nearest neighbors of a point given by
         *  its index.
         * \details For some implementation, may be faster than
         *  nearest neighbor search by point coordinates.
         * \param[in] nb_neighbors number of neighbors to be searched.
         *  Should be smaller or equal to nb_points() (else it triggers
         *  an assertion)
         * \param[in] query_point as the index of one of the points that
         *  was inserted in this NearestNeighborSearch
         * \param[out] neighbors array of nb_neighbors index_t
         * \param[out] neighbors_sq_dist array of nb_neighbors doubles
         */
        void get_nearest_neighbors(
            index_t nb_neighbors,
            index_t query_point,
            index_t* neighbors,
            double* neighbors_sq_dist
        ) const;

        /**
         * \brief Nearest neighbor search.
         * \param[in] query_point array of dimension() doubles
         * \return the index of the nearest neighbor from \p query_point
         */
        index_t get_nearest_neighbor(
            const double* query_point
        ) const {
            index_t result;
            double sq_dist;
            get_nearest_neighbors(1, query_point, &result, &sq_dist);
            geo_assert(signed_index_t(result) >= 0);
            return index_t(result);
        }

        
        /**
         * \brief Gets the dimension of the points.
         * \return the dimension
         */
        coord_index_t dimension() const {
            return dimension_;
        }

        /**
         * \brief Gets the number of points.
         * \return the number of points
         */
        index_t nb_points() const {
            return nb_points_;
        }

        /**
         * \brief Gets a point by its index
         * \param[in] i index of the point
         * \return a const pointer to the coordinates of the point
         */
        const double* point_ptr(index_t i) const {
            geo_debug_assert(i < nb_points());
            return points_ + i * stride_;
        }


    protected:

        /**
         * \brief Number of points stored in the leafs of the tree.
         */
        static const index_t MAX_LEAF_SIZE = 16;

        /**
         * \brief Returns the maximum node index in subtree.
         * \param[in] node_id node index of the subtree
         * \param[in] b first index of the points sequence in the subtree
         * \param[in] e one position past the last index of the point
         *  sequence in the subtree
         */
        static index_t max_node_index(
            index_t node_id, index_t b, index_t e
        ) {
            if(e - b <= MAX_LEAF_SIZE) {
                return node_id;
            }
            index_t m = b + (e - b) / 2;
            return std::max(
                max_node_index(2 * node_id, b, m),
                max_node_index(2 * node_id + 1, m, e)
            );
        }

        /**
         * \brief The context for traversing a KdTree.
         * \details Stores a sorted sequence of (point,distance)
         *  couples.
         */
        struct NearestNeighbors {

            /**
             * \brief Creates a new NearestNeighbors
             * \details Storage is provided
             * and managed by the caller.
             * Initializes neighbors_sq_dist[0..nb_neigh-1]
             * to Numeric::max_float64() and neighbors[0..nb_neigh-1]
             * to index_t(-1).
             * \param[in] nb_neighbors_in number of neighbors to retreive
             * \param[in] neighbors_in storage for the neighbors, allocated
             *  and managed by caller
             * \param[in] neighbors_sq_dist_in storage for neighbors squared
             *  distance, allocated and managed by caller
             */
            NearestNeighbors(
                index_t nb_neighbors_in,
                index_t* neighbors_in,
                double* neighbors_sq_dist_in
            ) :
                nb_neighbors(nb_neighbors_in),
                neighbors(neighbors_in),
                neighbors_sq_dist(neighbors_sq_dist_in) {
                for(index_t i = 0; i < nb_neighbors; ++i) {
                    neighbors[i] = index_t(-1);
                    neighbors_sq_dist[i] = Numeric::max_float64();
                }
            }

            /**
             * \brief Gets the squared distance to the furthest
             *  neighbor.
             */
            double furthest_neighbor_sq_dist() const {
                return neighbors_sq_dist[nb_neighbors - 1];
            }

            /**
             * \brief Inserts a new neighbor.
             * \details Only the nb_neighbor nearest points are kept.
             * \param[in] neighbor the index of the point
             * \param[in] sq_dist the squared distance between the point
             *  and the query point.
             */
            void insert(
                index_t neighbor, double sq_dist
            ) {
                if(sq_dist >= furthest_neighbor_sq_dist()) {
                    return;
                }
                index_t i = nb_neighbors;
                while(i != 0 && neighbors_sq_dist[i - 1] > sq_dist) {
                    if(i < nb_neighbors) {
                        neighbors[i] = neighbors[i - 1];
                        neighbors_sq_dist[i] = neighbors_sq_dist[i - 1];
                    }
                    --i;
                }
                geo_debug_assert(i < nb_neighbors);
                neighbors[i] = neighbor;
                neighbors_sq_dist[i] = sq_dist;
            }

            index_t nb_neighbors;
            index_t* neighbors;
            double* neighbors_sq_dist;
        };

        /**
         * \brief Computes the coordinate along which a point
         *   sequence will be splitted.
         * \param[in] b first index of the point sequence
         * \param[in] e one position past the last index of the point sequence
         */
        coord_index_t best_splitting_coord(index_t b, index_t e);

        /**
         * \brief Computes the extent of a point sequence along a coordinate.
         * \param[in] b first index of the point sequence
         * \param[in] e one position past the last index of the point sequence
         * \param[in] coord coordinate along which the extent is measured
         */
        double spread(
            index_t b, index_t e, coord_index_t coord
        ) {
            double minval = Numeric::max_float64();
            double maxval = Numeric::min_float64();
            for(index_t i = b; i < e; ++i) {
                double val = point_ptr(point_index_[i])[coord];
                minval = std::min(minval, val);
                maxval = std::max(maxval, val);
            }
            return maxval - minval;
        }

        /**
         * \brief Creates the subtree under a node.
         * \param[in] node_index index of the node that represents
         *  the subtree to create
         * \param[in] b first index of the point sequence in the subtree
         * \param[in] e one position past the last index of the point
         *  index in the subtree
         */
        void create_kd_tree_recursive(
            index_t node_index, index_t b, index_t e
        ) {
            if(e - b <= MAX_LEAF_SIZE) {
                return;
            }
            index_t m = split_kd_node(node_index, b, e);
            create_kd_tree_recursive(2 * node_index, b, m);
            create_kd_tree_recursive(2 * node_index + 1, m, e);
        }

        /**
         * \brief Computes and stores the splitting coordinate
         *  and splitting value of the node node_index, that
         *  corresponds to the [b,e) points sequence.
         *
         * \return a node index m. The point sequences
         *  [b,m) and [m,e) correspond to the left
         *  child (2*node_index) and right child (2*node_index+1)
         *  of node_index.
         */
        index_t split_kd_node(
            index_t node_index, index_t b, index_t e
        );

        /**
         * \brief The recursive function to implement KdTree traversal and
         *  nearest neighbors computation.
         * \details Traverses the subtree under the
         *  node_index node that corresponds to the
         *  [b,e) point sequence. Nearest neighbors
         *  are inserted into neighbors during
         *  traversal.
         * \param[in] node_index index of the current node in the Kd tree
         * \param[in] b index of the first point in the subtree under
         *  node \p node_index
         * \param[in] e one position past the index of the last point in the
         *  subtree under node \p node_index
         * \param[in,out] bbox_min coordinates of the lower
         *  corner of the bounding box.
         *  Allocated and managed by caller.
         *  Modified by the function and restored on exit.
         * \param[in,out] bbox_max coordinates of the
         *  upper corner of the bounding box.
         *  Allocated and managed by caller.
         *  Modified by the function and restored on exit.
         * \param[in] bbox_dist squared distance between
         *  the query point and a bounding box of the
         *  [b,e) point sequence. It is used to early
         *  prune traversals that do not generate nearest
         *  neighbors.
         * \param[in] query_point the query point
         * \param[in,out] neighbors the computed nearest neighbors
         */
        void get_nearest_neighbors_recursive(
            index_t node_index, index_t b, index_t e,
            double* bbox_min, double* bbox_max,
            double bbox_dist, const double* query_point,
            NearestNeighbors& neighbors
        ) const;

    protected:
        coord_index_t dimension_;
        index_t nb_points_;
        index_t stride_;
        const double* points_;

        vector<index_t> point_index_;
        vector<coord_index_t> splitting_coord_;
        vector<double> splitting_val_;
        vector<double> bbox_min_;
        vector<double> bbox_max_;
    };
}

#endif

