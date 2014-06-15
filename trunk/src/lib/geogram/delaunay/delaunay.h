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

#ifndef __GEOGRAM_DELAUNAY_DELAUNAY__
#define __GEOGRAM_DELAUNAY_DELAUNAY__

#include <geogram/basic.h>


namespace GEO {

    /************************************************************************/

    /**
     * \brief Abstract interface for Delaunay triangulation in Nd.
     */
    class GEOGRAM_API Delaunay {
    public:

        /**
         * \brief Gets the dimension of this Delaunay.
         * \return the dimension of this Delauna
         */
        coord_index_t dimension() const {
            return dimension_;
        }

        /**
         * \brief Gets the number of vertices in each cell
         * \details Cell_size =  dimension + 1
         * \return the number of vertices in each cell
         */
        index_t cell_size() const {
            return cell_size_;
        }

        /**
         * \brief Sets the vertices of this Delaunay, and recomputes the cells.
         * \param[in] nb_vertices number of vertices
         * \param[in] vertices a pointer to the coordinates of the vertices, as
         *  a contiguous array of doubles
         */
        virtual void set_vertices(
            index_t nb_vertices, const double* vertices
        );

        /**
         * \brief Gets a pointer to the array of vertices.
         * \return A const pointer to the array of vertices.
         */
        const double* vertices_ptr() const {
            return vertices_;
        }

        /**
         * \brief Gets a pointer to a vertex by its global index.
         * \param[in] i global index of the vertex
         * \return a pointer to vertex \p i
         */
        const double* vertex_ptr(index_t i) const {
            geo_debug_assert(i < nb_vertices());
            return vertices_ + vertex_stride_ * i;
        }

        /**
         * \brief Gets the number of vertices.
         * \return the number of vertices in this Delaunay
         */
        index_t nb_vertices() const {
            return nb_vertices_;
        }

        /**
         * \brief Gets the number of cells.
         * \return the number of cells in this Delaunay
         */
        index_t nb_cells() const {
            return nb_cells_;
        }

        /**
         * \brief Gets a pointer to the cell-to-vertex incidence array.
         * \return a const pointer to the cell-to-vertex incidence array
         */
        const signed_index_t* cell_to_v() const {
            return cell_to_v_;
        }

        /**
         * \brief Gets a pointer to the cell-to-cell adjacency array.
         * \return a const pointer to the cell-to-cell adjacency array
         */
        const signed_index_t* cell_to_cell() const {
            return cell_to_cell_;
        }

        /**
         * \brief Computes the nearest vertex from a query point.
         * \param[in] p query point
         * \return the index of the nearest vertex
         */
        virtual index_t nearest_vertex(const double* p) const;

        /**
         * \brief Gets a vertex index by cell index and local vertex index.
         * \param[in] c cell index
         * \param[in] lv local vertex index in cell \p c
         * \return the index of the lv-th vertex of cell c.
         */
        signed_index_t cell_vertex(index_t c, index_t lv) const {
            geo_debug_assert(c < nb_cells());
            geo_debug_assert(lv < cell_size());
            return cell_to_v_[c * cell_v_stride_ + lv];
        }

        /**
         * \brief Gets an adjacent cell index by cell index and
         *  local facet index.
         * \param[in] c cell index
         * \param[in] lf local facet index
         * \return the index of the cell adjacent to \p c accross
         *  facet \p lf if it exists, or -1 if on border
         */
        signed_index_t cell_adjacent(index_t c, index_t lf) const {
            geo_debug_assert(c < nb_cells());
            geo_debug_assert(lf < cell_size());
            return cell_to_cell_[c * cell_neigh_stride_ + lf];
        }

        /**
         * \brief Retreives a local vertex index from cell index
         *  and global vertex index.
         * \param[in] c cell index
         * \param[in] v global vertex index
         * \return the local index of vertex \p v in cell \p c
         * \pre cell \p c is incident to vertex \p v
         */
        index_t index(index_t c, signed_index_t v) const {
            geo_debug_assert(c < nb_cells());
            geo_debug_assert(v < (signed_index_t) nb_vertices());
            for(index_t iv = 0; iv < cell_size(); iv++) {
                if(cell_vertex(c, iv) == v) {
                    return iv;
                }
            }
            geo_assert_not_reached;
            return cell_size();
        }

        /**
         * \brief Retreives a local facet index from two adacent
         *  cell global indices.
         * \param[in] c1 global index of first cell
         * \param[in] c2 global index of second cell
         * \return the local index of the face accross which
         *  \p c2 is adjacent to \p c1
         * \pre cell \p c1 and cell \p c2 are adjacent
         */
        index_t adjacent_index(index_t c1, index_t c2) const {
            geo_debug_assert(c1 < nb_cells());
            geo_debug_assert(c2 < nb_cells());
            for(index_t f = 0; f < cell_size(); f++) {
                if(cell_adjacent(c1, f) == signed_index_t(c2)) {
                    return f;
                }
            }
            geo_assert_not_reached;
            return cell_size();
        }

    protected:
        /**
         * \brief Creates a new Delaunay triangulation
         * \details This creates a new Delaunay triangulation for the
         * specified \p dimension. 
         */
        Delaunay(coord_index_t dimension);

        /**
         * \brief Delaunay destructor.
         */
        virtual ~Delaunay();

        /**
         * \brief Sets the arrays that represent the combinatorics
         *  of this Delaunay.
         * \param[in] nb_cells number of cells
         * \param[in] cell_to_v the cell-to-vertex incidence array
         * \param[in] cell_to_cell the cell-to-cell adjacency array
         */
        virtual void set_arrays(
            index_t nb_cells,
            const signed_index_t* cell_to_v, const signed_index_t* cell_to_cell
        );

        /**
         * \brief Sets the dimension of this Delaunay.
         * \details Updates all the parameters related with
         *  the dimension. This includes vertex_stride (number
         *  of doubles between two consecutive vertices),
         *  cell size (number of vertices in a cell),
         *  cell_v_stride (number of integers between two
         *  consecutive cell vertex arrays),
         *  cell_neigh_stride (number of integers
         *  between two consecutive cell adjacency arrays).
         * \param[in] dim the dimension
         */
        void set_dimension(coord_index_t dim) {
            dimension_ = dim;
            vertex_stride_ = dim;
            cell_size_ = index_t(dim) + 1;
            cell_v_stride_ = cell_size_;
            cell_neigh_stride_ = cell_size_;
        }

        coord_index_t dimension_;
        index_t vertex_stride_;
        index_t cell_size_;
        index_t cell_v_stride_;
        index_t cell_neigh_stride_;
        const double* vertices_;
        index_t nb_vertices_;
        index_t nb_cells_;
        const signed_index_t* cell_to_v_;
        const signed_index_t* cell_to_cell_;
    };

    /************************************************************************/

    
}

#endif

