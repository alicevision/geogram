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

#include <geogram/points/kd_tree.h>
#include <geogram/basic/geometry_nd.h>
#include <geogram/basic/process.h>
#include <geogram/basic/algorithm.h>

namespace {

    using namespace GEO;

    /**
     * \brief Comparison functor used to
     * sort the point indices.
     */
    class ComparePointCoord {
    public:
        /**
         * \brief Creates a new ComparePointCoord
         * \param[in] nb_points number of points
         * \param[in] points pointer to first point
         * \param[in] stride number of doubles between two
         *  consecutive points in array (=dimension if point
         *  array is compact).
         * \param[in] splitting_coord the coordinate to compare
         */
        ComparePointCoord(
            index_t nb_points,
            const double* points,
            index_t stride,
            coord_index_t splitting_coord
        ) :
            nb_points_(nb_points),
            points_(points),
            stride_(stride),
            splitting_coord_(splitting_coord) {
        }

        /**
         * \brief Compares to point indices (does the
         * indirection and coordinate lookup).
         * \param[in] i index of first point to compare
         * \param[in] j index of second point to compare
         * \return true if point \p i is before point \p j, false otherwise
         */
        bool operator() (index_t i, index_t j) const {
            geo_debug_assert(i < nb_points_);
            geo_debug_assert(j < nb_points_);
            return
                (points_ + i * stride_)[splitting_coord_] <
                (points_ + j * stride_)[splitting_coord_]
            ;
        }

    private:
        index_t nb_points_;
        const double* points_;
        index_t stride_;
        coord_index_t splitting_coord_;
    };
}

/****************************************************************************/

namespace GEO {

    KdTree::KdTree(coord_index_t dim) :
        NearestNeighborSearch(dim),
        bbox_min_(dim),
        bbox_max_(dim),
        m0_(max_index_t()),
        m1_(max_index_t()),
        m2_(max_index_t()),
        m3_(max_index_t()),
        m4_(max_index_t()),
        m5_(max_index_t()),
        m6_(max_index_t()),
        m7_(max_index_t()),
        m8_(max_index_t()) {
    }

    KdTree::~KdTree() {
    }

    bool KdTree::stride_supported() const {
        return true;
    }

    void KdTree::set_points(
        index_t nb_points, const double* points, index_t stride
    ) {
        nb_points_ = nb_points;
        points_ = points;
        stride_ = stride;

        index_t sz = max_node_index(1, 0, nb_points) + 1;

        point_index_.resize(nb_points);
        splitting_coord_.resize(sz);
        splitting_val_.resize(sz);

        for(index_t i = 0; i < nb_points; i++) {
            point_index_[i] = i;
        }

        // If there are more than 16*MAX_LEAF_SIZE (=256) points,
        // create the tree in parallel
        if(
            nb_points >= (16 * MAX_LEAF_SIZE) &&
            Process::maximum_concurrent_threads() > 1
        ) {
            m0_ = 0;
            m8_ = nb_points;
            // Create the first level of the tree
            m4_ = split_kd_node(1, m0_, m8_);
            // Create the second level of the tree
            //  (using two threads)
            parallel_for(*this, 2, 4);
            // Create the third level of the tree
            //  (using four threads)
            parallel_for(*this, 4, 8);
            // Create the fourth level of the tree
            //  (using eight threads)
            parallel_for(*this, 8, 16);
        } else {
            create_kd_tree_recursive(1, 0, nb_points);
        }

        // Compute the bounding box.
        for(coord_index_t c = 0; c < dimension(); ++c) {
            bbox_min_[c] = 1e30;
            bbox_max_[c] = -1e30;
        }
        for(index_t i = 0; i < nb_points; ++i) {
            const double* p = point_ptr(i);
            for(coord_index_t c = 0; c < dimension(); ++c) {
                bbox_min_[c] = geo_min(bbox_min_[c], p[c]);
                bbox_max_[c] = geo_max(bbox_max_[c], p[c]);
            }
        }
    }

    void KdTree::set_points(
        index_t nb_points, const double* points
    ) {
        set_points(nb_points, points, dimension());
    }

    void KdTree::operator() (index_t i) {
        switch(i) {
            // Second level of the tree: create two nodes in
            //  parallel.
            case 2:
                m2_ = split_kd_node(2, m0_, m4_);
                break;
            case 3:
                m6_ = split_kd_node(3, m4_, m8_);
                break;

            // Third level of the tree: create four nodes in
            // parallel.
            case 4:
                m1_ = split_kd_node(4, m0_, m2_);
                break;
            case 5:
                m3_ = split_kd_node(5, m2_, m4_);
                break;
            case 6:
                m5_ = split_kd_node(6, m4_, m6_);
                break;
            case 7:
                m7_ = split_kd_node(7, m6_, m8_);
                break;

            // Fourth level of the tree: create eight subtrees
            // in parallel.
            case 8:
                create_kd_tree_recursive(8, m0_, m1_);
                break;
            case 9:
                create_kd_tree_recursive(9, m1_, m2_);
                break;
            case 10:
                create_kd_tree_recursive(10, m2_, m3_);
                break;
            case 11:
                create_kd_tree_recursive(11, m3_, m4_);
                break;
            case 12:
                create_kd_tree_recursive(12, m4_, m5_);
                break;
            case 13:
                create_kd_tree_recursive(13, m5_, m6_);
                break;
            case 14:
                create_kd_tree_recursive(14, m6_, m7_);
                break;
            case 15:
                create_kd_tree_recursive(15, m7_, m8_);
                break;
        }
    }

    index_t KdTree::split_kd_node(
        index_t node_index, index_t b, index_t e
    ) {

        geo_debug_assert(e > b);
        // Do not split leafs
        if(b + 1 == e) {
            return b;
        }

        coord_index_t splitting_coord = best_splitting_coord(b, e);
        index_t m = b + (e - b) / 2;
        geo_debug_assert(m < e);

        // sorts the indices in such a way that points's
        // coordinates splitting_coord in [b,m) are smaller
        // than m's and points in [m,e) are
        // greater or equal to m's
        std::nth_element(
            point_index_.begin() + std::ptrdiff_t(b),
            point_index_.begin() + std::ptrdiff_t(m),
            point_index_.begin() + std::ptrdiff_t(e),
            ComparePointCoord(
                nb_points_, points_, stride_, splitting_coord
            )
        );

        // Initialize node's variables (splitting coord and
        // splitting value)
        splitting_coord_[node_index] = splitting_coord;
        splitting_val_[node_index] =
            point_ptr(point_index_[m])[splitting_coord];
        return m;
    }

    coord_index_t KdTree::best_splitting_coord(
        index_t b, index_t e
    ) {
        // Returns the coordinates that maximizes
        // point's spread. We should probably
        // use a tradeoff between spread and
        // bbox shape ratio, as done in ANN, but
        // this simple method seems to give good
        // results in our case.
        coord_index_t result = 0;
        double max_spread = spread(b, e, 0);
        for(coord_index_t c = 1; c < dimension(); ++c) {
            double coord_spread = spread(b, e, c);
            if(coord_spread > max_spread) {
                result = c;
                max_spread = coord_spread;
            }
        }
        return result;
    }

    void KdTree::get_nearest_neighbors(
        index_t nb_neighbors,
        const double* query_point,
        index_t* neighbors,
        double* neighbors_sq_dist
    ) const {

        geo_debug_assert(nb_neighbors <= nb_points());

        // Compute distance between query point and global bounding box
        // and copy global bounding box to local variables (bbox_min, bbox_max),
        // allocated on the stack. bbox_min and bbox_max are updated during the
        // traversal of the KdTree (see get_nearest_neighbors_recursive()). They
        // are necessary to compute the distance between the query point and the
        // bbox of the current node.
        double box_dist = 0.0;
        double* bbox_min = (double*) (alloca(dimension() * sizeof(double)));
        double* bbox_max = (double*) (alloca(dimension() * sizeof(double)));
        for(coord_index_t c = 0; c < dimension(); ++c) {
            bbox_min[c] = bbox_min_[c];
            bbox_max[c] = bbox_max_[c];
            if(query_point[c] < bbox_min_[c]) {
                box_dist += geo_sqr(bbox_min_[c] - query_point[c]);
            } else if(query_point[c] > bbox_max_[c]) {
                box_dist += geo_sqr(bbox_max_[c] - query_point[c]);
            }
        }
        NearestNeighbors NN(
            nb_neighbors, neighbors, neighbors_sq_dist
        );
        get_nearest_neighbors_recursive(
            1, 0, nb_points(), bbox_min, bbox_max, box_dist, query_point, NN
        );
    }


    void KdTree::get_nearest_neighbors(
        index_t nb_neighbors,
        index_t q_index,
        index_t* neighbors,
        double* neighbors_sq_dist
    ) const {
        // TODO: optimized version that uses the fact that
        // we know that query_point is in the search data
        // structure already.
        // (I tryed something already, see in the Attic, 
        //  but it did not give any significant speedup).
        get_nearest_neighbors(
            nb_neighbors, point_ptr(q_index), 
            neighbors, neighbors_sq_dist
        );
    }

    void KdTree::get_nearest_neighbors_recursive(
        index_t node_index, index_t b, index_t e,
        double* bbox_min, double* bbox_max, double box_dist,
        const double* query_point, NearestNeighbors& NN
    ) const {
        geo_debug_assert(e > b);

        // Simple case (node is a leaf)
        if((e - b) <= MAX_LEAF_SIZE) {
            for(index_t i = b; i < e; ++i) {
                index_t p = point_index_[i];
                double d2 = Geom::distance2(
                    query_point, point_ptr(p), dimension()
                );
                NN.insert(p, d2);
            }
            return;
        }

        coord_index_t coord = splitting_coord_[node_index];
        double val = splitting_val_[node_index];
        double cut_diff = query_point[coord] - val;
        index_t m = b + (e - b) / 2;

        // If the query point is on the left side
        if(cut_diff < 0.0) {

            // Traverse left subtree
            {
                double bbox_max_save = bbox_max[coord];
                bbox_max[coord] = val;
                get_nearest_neighbors_recursive(
                    2 * node_index, b, m, 
                    bbox_min, bbox_max, box_dist, query_point, NN
                );
                bbox_max[coord] = bbox_max_save;
            }

            // Update bbox distance (now measures the
            // distance to the bbox of the right subtree)
            double box_diff = bbox_min[coord] - query_point[coord];
            if(box_diff > 0.0) {
                box_dist -= geo_sqr(box_diff);
            }
            box_dist += geo_sqr(cut_diff);

            // Traverse the right subtree, only if bbox
            // distance is nearer than furthest neighbor,
            // else there is no chance that the right
            // subtree contains points that will change
            // anything in the nearest neighbors NN.
            if(box_dist <= NN.furthest_neighbor_sq_dist()) {
                double bbox_min_save = bbox_min[coord];
                bbox_min[coord] = val;
                get_nearest_neighbors_recursive(
                    2 * node_index + 1, m, e, 
                    bbox_min, bbox_max, box_dist, query_point, NN
                );
                bbox_min[coord] = bbox_min_save;
            }
        } else {
            // else the query point is on the right side
            // (then do the same with left and right subtree
            //  permutted).
            {
                double bbox_min_save = bbox_min[coord];
                bbox_min[coord] = val;
                get_nearest_neighbors_recursive(
                    2 * node_index + 1, m, e, 
                    bbox_min, bbox_max, box_dist, query_point, NN
                );
                bbox_min[coord] = bbox_min_save;
            }

            // Update bbox distance (now measures the
            // distance to the bbox of the left subtree)
            double box_diff = query_point[coord] - bbox_max[coord];
            if(box_diff > 0.0) {
                box_dist -= geo_sqr(box_diff);
            }
            box_dist += geo_sqr(cut_diff);

            if(box_dist <= NN.furthest_neighbor_sq_dist()) {
                double bbox_max_save = bbox_max[coord];
                bbox_max[coord] = val;
                get_nearest_neighbors_recursive(
                    2 * node_index, b, m, 
                    bbox_min, bbox_max, box_dist, query_point, NN
                );
                bbox_max[coord] = bbox_max_save;
            }
        }
    }
}

