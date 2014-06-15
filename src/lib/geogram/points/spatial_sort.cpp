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

#include <geogram/points/spatial_sort.h>
#include <algorithm>

namespace {

    using namespace GEO;

    /**
     * \brief Splits a sequence into two ordered halves.
     * \details The algorithm shuffles the sequence and
     *  partitions its into two halves with the same number of elements
     *  and such that the elements of the first half are smaller
     *  than the elements of the second half.
     * \param[in] begin an iterator to the first element
     * \param[in] end an iterator one position past the last element
     * \param[in] cmp the comparator object
     * \return an iterator to the middle of the sequence that separates
     *  the two halves
     */
    template <class IT, class CMP>
    inline IT split(
        IT begin, IT end, CMP cmp
    ) {
        if(begin >= end) {
            return begin;
        }
        IT middle = begin + (end - begin) / 2;
        std::nth_element(begin, middle, end, cmp);
        return middle;
    }

    /************************************************************************/

    /**
     * \brief Exposes an interface compatible with the requirement
     * of Hilbert sort templates for a raw array of vertices.
     */
    class VertexArray {
    public:
        /**
         * \brief Constructs a new VertexArray.
         * \param[in] base address of the points
         * \param[in] stride number of doubles between
         *  two consecutive points
         */
        VertexArray(
            index_t nb_vertices,
            const double* base, index_t stride
        ) :
            base_(base),
            stride_(stride) {
            nb_vertices_ = nb_vertices;
        }

        /**
         * \brief Gets a vertex by its index.
         * \param[in] i the index of the point
         * \return a const pointer to the coordinates of the vertex
         */
        const double* vertex_ptr(index_t i) const {
            geo_debug_assert(i < nb_vertices_);
            return base_ + i * stride_;
        }

    private:
        const double* base_;
        index_t stride_;
        index_t nb_vertices_;
    };

    /************************************************************************/

    /**
     * \brief The generic comparator class for Hilbert vertex
     *  ordering.
     * \tparam COORD the coordinate to compare
     * \tparam UP    if true, use direct order, else use reverse order
     * \tparam MESH  the class that represents meshes
     */
    template <int COORD, bool UP, class MESH>
    struct Hilbert_vcmp {
    };

    /**
     * \brief Specialization (UP=true) of the generic comparator class
     *  for Hilbert vertex ordering.
     * \see Hilbert_vcmp
     * \tparam COORD the coordinate to compare
     * \tparam MESH  the class that represents meshes
     */
    template <int COORD, class MESH>
    struct Hilbert_vcmp<COORD, true, MESH> {

        /**
         * \brief Constructs a new Hilbert_vcmp.
         * \param[in] mesh the mesh in which the compared
         *  points reside.
         */
        Hilbert_vcmp(const MESH& mesh) :
            mesh_(mesh) {
        }

        /**
         * \brief Compares two points.
         * \param[in] i1 index of the first point to compare
         * \param[in] i2 index of the second point to compare
         * \return true if point \p i1 is before point \p i2,
         *  false otherwise.
         */
        bool operator() (index_t i1, index_t i2) {
            return
                mesh_.vertex_ptr(i1)[COORD] <
                mesh_.vertex_ptr(i2)[COORD];
        }

        const MESH& mesh_;
    };

    /**
     * \brief Specialization (UP=false) of the generic comparator class
     *  for Hilbert vertex ordering.
     * \see Hilbert_vcmp
     * \tparam COORD the coordinate to compare
     * \tparam MESH  the class that represents meshes
     */
    template <int COORD, class MESH>
    struct Hilbert_vcmp<COORD, false, MESH> {

        /**
         * \brief Constructs a new Hilbert_vcmp.
         * \param[in] mesh the mesh in which the compared
         *  points reside.
         */
        Hilbert_vcmp(const MESH& mesh) :
            mesh_(mesh) {
        }

        /**
         * \brief Compares two points.
         * \param[in] i1 index of the first point to compare
         * \param[in] i2 index of the second point to compare
         * \return true if point \p i1 is before point \p i2,
         *  false otherwise.
         */
        bool operator() (index_t i1, index_t i2) {
            return
                mesh_.vertex_ptr(i1)[COORD] >
                mesh_.vertex_ptr(i2)[COORD];
        }

        const MESH& mesh_;
    };

    /************************************************************************/

    /**
     * \brief Generic class for sorting arbitrary elements in
     *  Hilbert and Morton orders.
     * \details The implementation is inspired by:
     *  Christophe Delage and Olivier Devillers. 
     *  Spatial Sorting. In CGAL User and Reference Manual. 
     *  CGAL Editorial Board, 3.9 edition, 2011
     * \tparam CMP the comparator class for ordering the elements. CMP
     *  is itself a template parameterized by~:
     *    - COORD the coordinate along which elements should be
     *      sorted
     *    - UP a boolean that indicates whether direct or reverse
     *      order should be used
     *    - MESH the class that represents meshes
     * \tparam MESH  the class that represents meshes
     */
    template <template <int COORD, bool UP, class MESH> class CMP, class MESH>
    struct HilbertSort {

        /**
         * \brief Low-level recursive spatial sorting function
         * \details This function is recursive
         * \param[in] M the mesh in which the elements reside
         * \param[in] begin an iterator that points to the
         *  first element of the sequence
         * \param[in] end an interator that points one position past the
         *  last element of the sequence
         * \param[in] limit subsequences smaller than limit are left unsorted
         * \tparam COORDX the first coordinate, can be 0,1 or 2. The second
         *  and third coordinates are COORDX+1 modulo 3 and COORDX+2 modulo 3
         *  respectively
         * \tparam UPX whether ordering along the first coordinate
         *  is direct or inverse
         * \tparam UPY whether ordering along the second coordinate
         *  is direct or inverse
         * \tparam UPZ whether ordering along the third coordinate
         *  is direct or inverse
         */
        template <int COORDX, bool UPX, bool UPY, bool UPZ, class IT>
        static void sort(
            const MESH& M, IT begin, IT end, index_t limit = 1
        ) {
            const int COORDY = (COORDX + 1) % 3, COORDZ = (COORDY + 1) % 3;
            if(end - begin <= signed_index_t(limit)) {
                return;
            }
            IT m0 = begin, m8 = end;
            IT m4 = split(m0, m8, CMP<COORDX, UPX, MESH>(M));
            IT m2 = split(m0, m4, CMP<COORDY, UPY, MESH>(M));
            IT m1 = split(m0, m2, CMP<COORDZ, UPZ, MESH>(M));
            IT m3 = split(m2, m4, CMP<COORDZ, !UPZ, MESH>(M));
            IT m6 = split(m4, m8, CMP<COORDY, !UPY, MESH>(M));
            IT m5 = split(m4, m6, CMP<COORDZ, UPZ, MESH>(M));
            IT m7 = split(m6, m8, CMP<COORDZ, !UPZ, MESH>(M));
            sort<COORDZ, UPZ, UPX, UPY>(M, m0, m1);
            sort<COORDY, UPY, UPZ, UPX>(M, m1, m2);
            sort<COORDY, UPY, UPZ, UPX>(M, m2, m3);
            sort<COORDX, UPX, !UPY, !UPZ>(M, m3, m4);
            sort<COORDX, UPX, !UPY, !UPZ>(M, m4, m5);
            sort<COORDY, !UPY, UPZ, !UPX>(M, m5, m6);
            sort<COORDY, !UPY, UPZ, !UPX>(M, m6, m7);
            sort<COORDZ, !UPZ, !UPX, UPY>(M, m7, m8);
        }

        /**
         * \brief Sorts a sequence of elements spatially.
         * \details This function does an indirect sort, in the sense that a sequence 
         *  of indices that refer to the elements is sorted. This function uses a 
         *  multithreaded implementation.
         * \param[in] M the mesh in which the elements to sort reside
         * \param[in] b an interator to the first index to be sorted
         * \param[in] e an interator one position past the last index to be sorted
         * \param[in] limit subsequences smaller than limit are left unsorted
         */
        HilbertSort(
            const MESH& M,
            vector<index_t>::iterator b,
            vector<index_t>::iterator e,
            index_t limit = 1
        ) {
            geo_debug_assert(e > b);
            sort<0,false,false,false>(M, b, e, limit);
        }
    };

    /************************************************************************/

    /**
     * \brief Computes the BRIO order for a set of 3D points.
     * \details It is used to accelerate incremental insertion
     *  in Delaunay triangulation
     * \param[in] nb_vertices number of vertices to sort
     * \param[in] vertices pointer to the coordinates of the vertices
     * \param[in] stride number of doubles between two consecutive vertices
     * \param[in,out] sorted_indices indices to sort
     * \param[in] b iterator to the first index to sort
     * \param[in] e iterator one position past the last index to sort
     * \param[in] threshold minimum size of interval to be sorted
     * \param[in] ratio splitting ratio between current interval and
     *  the rest to be sorted
     * \param[in,out] depth iteration depth
     * \param[out] levels if non-null, bounds of each level
     */
    void compute_BRIO_order_recursive(
        index_t nb_vertices, const double* vertices,
        index_t stride,
        vector<index_t>& sorted_indices,
        vector<index_t>::iterator b,
        vector<index_t>::iterator e,
        index_t threshold,
        double ratio,
        index_t& depth,
        vector<index_t>* levels
    ) {
        geo_debug_assert(e > b);

        vector<index_t>::iterator m = b;
        if(index_t(e - b) > threshold) {
            ++depth;
            m = b + int(double(e - b) * ratio);
            compute_BRIO_order_recursive(
                nb_vertices, vertices, stride,
                sorted_indices, b, m,
                threshold, ratio, depth,
                levels
            );
        }

        VertexArray M(nb_vertices, vertices, stride);
        HilbertSort<Hilbert_vcmp, VertexArray>(
            M, m, e
        );

        if(levels != nil) {
            levels->push_back(index_t(e - sorted_indices.begin()));
        }
    }
}

/****************************************************************************/

namespace GEO {

    void compute_Hilbert_order(
        index_t nb_vertices, const double* vertices,
        vector<index_t>& sorted_indices,
        index_t stride
    ) {
        sorted_indices.resize(nb_vertices);
        for(index_t i = 0; i < nb_vertices; ++i) {
            sorted_indices[i] = i;
        }
        VertexArray M(nb_vertices, vertices, stride);
        HilbertSort<Hilbert_vcmp, VertexArray>(
            M, sorted_indices.begin(), sorted_indices.end()
        );
    }

    void compute_BRIO_order(
        index_t nb_vertices, const double* vertices,
        vector<index_t>& sorted_indices,
        index_t stride,
        index_t threshold,
        double ratio,
        vector<index_t>* levels
    ) {
        if(levels != nil) {
            levels->clear();
            levels->push_back(0);
        }
        index_t depth = 0;
        sorted_indices.resize(nb_vertices);
        for(index_t i = 0; i < nb_vertices; ++i) {
            sorted_indices[i] = i;
        }
        std::random_shuffle(sorted_indices.begin(), sorted_indices.end());
        compute_BRIO_order_recursive(
            nb_vertices, vertices, stride,
            sorted_indices,
            sorted_indices.begin(), sorted_indices.end(),
            threshold, ratio, depth, levels
        );
    }
}

