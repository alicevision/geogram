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
#include <geogram/numerics/predicates.h>
#include <string>

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

    /**
     * \brief Implementation of Delaunay triangulation and 
     *  regular triangulation in 3d.
     * \details This package uses concepts inspired by 
     *  two triangulation softwares, CGAL and tetgen,
     *  described in the following references. This package follows the
     *  idea used in CGAL of traversing the cavity from inside, since
     *  it traverses less tetrahedra than when traversing from outside.
     *  - Jean-Daniel Boissonnat, Olivier Devillers, Monique Teillaud, 
     *   and Mariette Yvinec. Triangulations in CGAL. 
     *   In Proc. 16th Annu. ACM Sympos. Comput. Geom., pages 11–18, 2000.
     *  - Hang Si, Constrained Delaunay tetrahedral mesh generation and 
     *   refinement. Finite elements in Analysis and Design, 
     *   46 (1-2):33--46, 2010.
     *
     *  Note that the algorithm here does not support vertex deletion nor
     *  degenerate input with all coplanar or all colinear points (use CGAL
     *  instead if you have these requirements).
     *
     *  The core algorithm used in this code, CGAL and tetgen was
     *  independently and simultaneously discovered by Bowyer and Watson:
     *  - Adrian Bowyer, "Computing Dirichlet tessellations", 
     *   Comput. J., vol. 24, no 2,‎ 1981, p. 162-166 
     *  - David F. Watson, "Computing the n-dimensional Delaunay tessellation 
     *   with application to Voronoi polytopes", Comput. J., vol. 24, 
     *   no 2, 1981, p. 167-172
     *
     *  The spatial reordering method, that dramatically increases the 
     *  performances, also used in this code, CGAL and tetgen was introduced
     *  in the following references. The second one is a smart implementation
     *  based on the std::nth_element() function of the STL, that inspired
     *  the compute_BRIO_ordering() function of this package.
     *  - Nina Amenta, Sunghee Choi and Gunter Rote, "Incremental constructions
     *   con brio", ACM Symposium on Computational Geometry 2003.
     *  - Christophe Delage and Olivier Devillers. Spatial Sorting. 
     *   In CGAL User and Reference Manual. CGAL Editorial Board, 
     *   3.9 edition, 2011
     *
     *  The locate() function is based on the following two references. 
     *  The first one randomizes the choice of the next tetrahedron.
     *  The second one uses an inexact locate() function to initialize 
     *  the exact one (it is called "structural filtering"). The first
     *  idea is used in both CGAL and tetgen, and the second one is used
     *  in CGAL.
     *  - Walking in a triangulation, O Devillers, S Pion, M Teillaud
     *   17th Annual Symposium on Computational geometry, 106-114
     *  - Stefan Funke , Kurt Mehlhorn and Stefan Naher, "Structural filtering,
     *  a paradigm for efficient and exact geometric programs", 1999
     *
     */
    class GEOGRAM_API Delaunay3d : public Delaunay {
    public:
        /**
         * \brief Constructs a new Delaunay3d.
         * \param[in] weighted if true, this creates a regular triangulation
         *  (dual of a power diagram). In this case:
         *  - the input points are 4d points, were the fourth coordinate
         *   of point \f$ i \f$ is \f$ \sqrt{W - w_i} \f$ where \f$ W \f$ is
         *   the maximum of the  weights of all the points and \d$ w_i \$ is
         *   the weight associated with vertex \f$ i \f$.
         *  - the constructed combinatorics is a tetrahedralized volume (3d and
         *   not 4d although dimension() returns 4). This tetrahedralized volume
         *   corresponds to the regular triangulation of the weighted points.
         */
        Delaunay3d(bool weighted = false);

        /**
         * \brief Delaunay3d destructor
         */
        virtual ~Delaunay3d();

        virtual void set_vertices(
            index_t nb_vertices, const double* vertices
        );

        virtual index_t nearest_vertex(const double* p) const;

        /**
         * \brief Specifies whether vertices should be reordered.
         * \details Reordering is activated by default. Some special
         *  usages of Delaunay3d may require to deactivate it (for
         *  instance if vertices are already known to be ordered).
         * \param[in] x if true, then vertices are reordered using
         *  BRIO-Hilbert ordering. This improves speed significantly
         *  (enabled by default).
         */
        void set_reorder(bool x) {
            do_reorder_ = x;
        }

    protected:

        /**
         * \brief Symbolic constant for uninitialized hint.
         * \details Locate functions can be accelerated by
         *  specifying a hint. This constant indicates that
         *  no hint is given.
         */
        static const index_t NO_TETRAHEDRON = index_t(-1);

        /**
         * \brief Finds in the pointset a set of four non-coplanar
         *  points.
         * \details This function is used to initiate the incremental
         *  Delaunay construction.
         * \param[out] iv0 index of the first vertex
         * \param[out] iv1 index of the second vertex
         * \param[out] iv2 index of the third vertex
         * \param[out] iv3 index of the fourth vertex
         * \retval true if a set of four non-coplanar points was found
         * \retval false if all the points are coplanar
         */
        bool create_first_tetrahedron(
            index_t& iv0, index_t& iv1, index_t& iv2, index_t& iv3
        );

        /**
         * \brief Finds the tetrahedron that contains a point.
         * \details If the point is on a face, edge or vertex,
         *  the function returns one of the tetrahedra incident
         *  to that face, edge or vertex.
         * \param[in] p a pointer to the coordinates of the point
         * \param[out] orient a pointer to an array of four Sign%s
         *  or nil. If non-nil, returns the orientation with respect
         *  to the four facets of the tetrahedron that contains \p p.
         * \return the index of a tetrahedron that contains \p p.
         *  If the point is outside the convex hull of
         *  the inserted so-far points, then the returned tetrahedron
         *  is a virtual one (first vertex is the "vertex at infinity"
         *  of index -1) or NO_TETRAHEDRON if the virtual tetrahedra 
         *  were previously removed.
         */
         index_t locate(
            const double* p, index_t hint = NO_TETRAHEDRON,
            Sign* orient = nil
         ) const;
         
        /**
         * \brief Finds the tetrahedron that (approximately) 
         *  contains a point using inexact predicates.
         * \details The result of this function can be used as a hint
         *  for locate(). It accelerates locate as compared to calling
         *  it directly. This technique is referred to as "structural
         *  filtering".
         * \param[in] p a pointer to the coordinates of the point
         * \param[in] max_iter maximum number of traversed tets
         * \return the index of a tetrahedron that (approximately) 
         *  contains \p p.
         *  If the point is outside the convex hull of
         *  the inserted so-far points, then the returned tetrahedron
         *  is a virtual one (first vertex is the "vertex at infinity"
         *  of index -1) or NO_TETRAHEDRON if the virtual tetrahedra 
         *  were previously removed.
         */
         index_t locate_inexact(
             const double* p, index_t hint, index_t max_iter
         ) const;

        /**
         * \brief Inserts a point in the triangulation.
         * \param[in] v the index of the point to be inserted
         * \param[in] hint the index of a tetrahedron as near as
         *  possible to \p v, or -1 if unspecified
         * \return the index of one of the tetrahedra incident to
         *  point \p v
         */
         index_t insert(index_t v, index_t hint = NO_TETRAHEDRON);

        /**
         * \brief Determines the list of tetrahedra in conflict
         *  with a given point.
         * \param[in] v the index of the point to be inserted
         * \param[in] t the index of a tetrahedron that contains
         *  \p p, as returned by locate()
         * \param[in] orient an array of four signs indicating
         *  the orientation of \p p with respect to the four
         *  faces of \p t, as returned by locate()
         * \param[out] t_bndry a tetrahedron adjacent to the
         *  boundary of the conflict zone
         * \param[out] f_boundary the facet along which t_bndry is
         *  adjacent to the boundary of the conflict zone
         * \param[out] first the index of the first tetrahedron in conflict
         * \param[out] last the index of the last tetrahedron in conflict
         *  The other tetrahedra are linked, and can be traversed 
         *  from \p first by using tet_next() until \p last or END_OF_LIST 
         *  is reached.
         *  The conflict zone can be empty under two circumstances:
         *  - the vertex \p v already exists in the triangulation
         *  - the triangulation is weighted and \p v is not visible
         *  in either cases, both \p first and \p last contain END_OF_LIST
         */
         void find_conflict_zone(
             index_t v, 
             index_t t, const Sign* orient,
             index_t& t_bndry, index_t& f_bndry,
             index_t& first, index_t& last
         );
         
         
         /**
          * \brief This function is used to implement find_conflict_zone.
          * \details This function detects the neighbors of \p t that are
          *  in the conflict zone and calls itself recursively on them.
          * \param[in] p the point to be inserted
          * \param[in] t index of a tetrahedron in the fonflict zone
          * \param[out] t_bndry a tetrahedron adjacent to the
          *  boundary of the conflict zone
          * \param[out] f_boundary the facet along which t_bndry is
          *  adjacent to the boundary of the conflict zone
          * \param[out] first the index of the first tetrahedron in conflict
          * \param[out] last the index of the last tetrahedron in conflict
          * \pre The tetrahedron \p t was alredy marked as 
          *  conflict (tet_is_in_list(t))
          */
         void find_conflict_zone_recursive(
             const double* p, index_t t,
             index_t& t_bndry, index_t& f_bndry,
             index_t& first, index_t& last
         );

         /**
          * \brief Creates a star of tetrahedra filling the conflict
          *  zone.
          * \details For each tetrahedron facet on the border of the
          *  conflict zone, a new tetrahedron is created, resting on
          *  the facet and incident to vertex \p v. The function is 
          *  called recursively until the entire conflict zone is filled.
          * \param[in] v the index of the point to be inserted
          * \param[in] t_bndry index of a tetrahedron on the border
          *  of the conflict zone.
          * \param[in] f_bndry index of the facet along which \p t_bndry
          *  is incident to the border of the conflict zone
          * \param[in] prev_f the facet of \p t_bndry connected to the
          *  tetrahedron that \p t_bndry was reached from
          * \return the index of one the newly created tetrahedron
          */
         index_t stellate_conflict_zone(
             index_t v, 
             index_t t_bndry, index_t f_bndry, 
             index_t prev_f=index_t(~0)
         );

         // _________ Combinatorics - new and delete _________________________

         /**
         * \brief Maximum valid index for a tetrahedron.
         * \details This includes not only real tetrahedra,
         *  but also the virtual ones on the border, the conflict
         *  list and the free list.
         * \return the maximum valid index for a tetrahedron
         */
        index_t max_t() const {
            return index_t(cell_to_v_store_.size() / 4);
        }


        /**
         * \brief Default symbolic value of the cell_next_ field
         *  that indicates that a tetrahedron is not
         *  in a linked list.
         * \details This is the default value. Note that it suffices
         *  that NOT_IN_LIST_BIT is set for a tetrahedron
         *  to be not in any list.
         * A tetrahedron can be:
         *  - in a list (cell_next_[t] & NOT_IN_LIST_BIT == 0)
         *  - not in a list and not marked 
         *    (cell_next_[t] & NOT_IN_LIST_BIT != 0) && 
         *    (cell_next_[t] != cur_stamp_)
         *  - not in a list and marked
         *    (cell_next_[t] == cur_stamp_)
         */
        static const index_t NOT_IN_LIST  = index_t(~0);

        /**
         * \brief If cell_next_[t] & NOT_IN_LIST_BIT != 0,
         *  then t is not in a linked list.
         * \details The other bits of cell_next_[t] are used
         *  to store the stamp (i.e. index of the current point
         *  being inserted). The stamp is used for marking tetrahedra
         *  that were detected as non-conflict when inserting a point.
         * A tetrahedron can be:
         *  - in a list (cell_next_[t] & NOT_IN_LIST_BIT == 0)
         *  - not in a list and not marked 
         *    (cell_next_[t] & NOT_IN_LIST_BIT != 0) && 
         *    (cell_next_[t] != cur_stamp_)
         *  - not in a list and marked
         *    (cell_next_[t] == cur_stamp_)
         */
        static const index_t NOT_IN_LIST_BIT = index_t(1 << 31);

        /**
         * \brief Symbolic value of the cell_next_ field
         *  that indicates the end of list in a linked
         *  list of tetrahedra.
         */
        static const index_t END_OF_LIST = ~(NOT_IN_LIST_BIT);


        /**
         * \brief Tests whether a tetrahedron belongs to a linked
         *  list.
         * \details Tetrahedra can be linked, it is used to manage
         *  both the free list that recycles deleted tetrahedra,
         *  the conflict region and the list of newly created
         *  tetrahedra. In addition, a tetrahedron that is not
         *  in a list can be marked. The same space is used for
         *  marking and chaining tetrahedra in lists.
         *  A tetrahedron can be in the following states:
         *  - in list
         *  - not in list and marked
         *  - not in list and not marked
         * \param[in] t the index of the tetrahedron
         * \retval true if tetrahedron \p t belongs to a linked list
         * \retval false otherwise
         */
        bool tet_is_in_list(index_t t) const {
            geo_debug_assert(t < max_t());
            return (cell_next_[t] & NOT_IN_LIST_BIT) == 0;
        }

        /**
         * \brief Gets the index of a successor of a tetrahedron.
         * \details Tetrahedra can be linked, it is used to manage
         *  both the free list that recycles deleted tetrahedra.
         * \param[in] t the index of the tetrahedron
         * \retval END_OF_LIST if the end of the list is reached
         * \retval the index of the successor of
         *   tetrahedron \t otherwise
         * \pre tet_is_in_list(t)
         */
        index_t tet_next(index_t t) const {
            geo_debug_assert(t < max_t());
            geo_debug_assert(tet_is_in_list(t));
            return cell_next_[t];
        }

        /**
         * \brief Adds a tetrahedron to a linked list.
         * \details Tetrahedra can be linked, it is used to manage
         *  the free list that recycles deleted tetrahedra.
         * \param[in] t the index of the tetrahedron
         * \param[in,out] first first item of the list or END_OF_LIST if
         *  the list is empty
         * \param[in,out] last last item of the list or END_OF_LIST if
         *  the list is empty
         */
        void add_tet_to_list(index_t t, index_t& first, index_t& last) {
            geo_debug_assert(t < max_t());
            geo_debug_assert(!tet_is_in_list(t));
            if(last == END_OF_LIST) {
                geo_debug_assert(first == END_OF_LIST);
                first = last = t;
                cell_next_[t] = END_OF_LIST;
            } else {
                cell_next_[t] = first;
                first = t;
            }
        }

        /**
         * \brief Removes a tetrahedron from the linked list it
         *  belongs to.
         * \details Tetrahedra can be linked, it is used to manage
         *  both the free list that recycles deleted tetrahedra.
         * \param[in] t the index of the tetrahedron
         */
        void remove_tet_from_list(index_t t) {
            geo_debug_assert(t < max_t());
            geo_debug_assert(tet_is_in_list(t));
            cell_next_[t] = NOT_IN_LIST;
        }

        /**
         * \brief Symbolic value for a vertex of a
         *  tetrahedron that indicates a virtual tetrahedron.
         * \details The three other vertices then correspond to a
         *  facet on the convex hull of the points.
         */
        static const signed_index_t VERTEX_AT_INFINITY = -1;

        /**
         * \brief Tests whether a tetrahedron is
         *  a real one.
         * \details Real tetrahedra are incident to
         *  four user-specified vertices (there are also
         *  virtual tetrahedra that are incident to the
         *  vertex at infinity, with index -1)
         * \param[in] t index of the tetrahedron
         * \retval true if tetrahedron \p t is a real one
         * \retval false otherwise
         */
        bool tet_is_real(index_t t) const {
            return
                !tet_is_free(t) &&
                cell_to_v_store_[4 * t] >= 0 &&
                cell_to_v_store_[4 * t + 1] >= 0 &&
                cell_to_v_store_[4 * t + 2] >= 0 &&
                cell_to_v_store_[4 * t + 3] >= 0;
        }

        /**
         * \brief Tests whether a tetrahedron is
         *  a virtual one.
         * \details Virtual tetrahedra are tetrahedra
         *  incident to the vertex at infinity.
         * \param[in] t index of the tetrahedron
         * \retval true if tetrahedron \p t is virtual
         * \retval false otherwise
         */
        bool tet_is_virtual(index_t t) const {
            return
                !tet_is_free(t) && (
                cell_to_v_store_[4 * t] == VERTEX_AT_INFINITY ||
                cell_to_v_store_[4 * t + 1] == VERTEX_AT_INFINITY ||
                cell_to_v_store_[4 * t + 2] == VERTEX_AT_INFINITY ||
                cell_to_v_store_[4 * t + 3] == VERTEX_AT_INFINITY) ;
        }

        /**
         * \brief Tests whether a tetrahedron is
         *  in the free list.
         * \details Deleted tetrahedra are recycled
         *  in a free list.
         * \param[in] t index of the tetrahedron
         * \retval true if tetrahedron \p t is in
         * the free list
         * \retval false otherwise
         */
        bool tet_is_free(index_t t) const {
            return tet_is_in_list(t);
        }

        /**
         * \brief Creates a new tetrahedron.
         * \details Uses either a tetrahedron recycled
         *  from the free list, or creates a new one by
         *  expanding the two indices arrays.
         * \return the index of the newly created tetrahedron
         */
        index_t new_tetrahedron() {
            index_t result;
            if(first_free_ == END_OF_LIST) {
                cell_to_v_store_.resize(cell_to_v_store_.size() + 4, -1);
                cell_to_cell_store_.resize(cell_to_cell_store_.size() + 4, -1);
                // index_t(NOT_IN_LIST) is necessary else with
                // NOT_IN_LIST alone the compiler tries to generate a
                // reference to NOT_IN_LIST resulting in a link error.
                cell_next_.push_back(index_t(NOT_IN_LIST));
                result = max_t() - 1;
            } else {
                result = first_free_;
                first_free_ = tet_next(first_free_);
                remove_tet_from_list(result);
            }

            cell_to_cell_store_[4 * result] = -1;
            cell_to_cell_store_[4 * result + 1] = -1;
            cell_to_cell_store_[4 * result + 2] = -1;
            cell_to_cell_store_[4 * result + 3] = -1;

            return result;
        }

        /**
         * \brief Creates a new tetrahedron.
         * \details Sets the vertices. Adjacent tetrahedra index are
         *  left uninitialized. Uses either a tetrahedron recycled
         *  from the free list, or creates a new one by
         *  expanding the two indices arrays.
         * \param[in] v1 index of the first vertex
         * \param[in] v2 index of the second vertex
         * \param[in] v3 index of the third vertex
         * \param[in] v4 index of the fourth vertex
         * \return the index of the newly created tetrahedron
         */
        index_t new_tetrahedron(
            signed_index_t v1, signed_index_t v2, 
            signed_index_t v3, signed_index_t v4
        ) {
            index_t result = new_tetrahedron();
            cell_to_v_store_[4 * result] = v1;
            cell_to_v_store_[4 * result + 1] = v2;
            cell_to_v_store_[4 * result + 2] = v3;
            cell_to_v_store_[4 * result + 3] = v4;
            return result;
        }

        /**
         * \brief Generates a unique stamp for marking tets.
         * \details Storage is shared for list-chaining and stamp-marking 
         * both are mutually exclusive), therefore the stamp has
         * the NOT_IN_LIST_BIT set.
         * \param[in] stamp the unique stamp for marking tets
         */
        void set_tet_mark_stamp(index_t stamp) {
            cur_stamp_ = (stamp | NOT_IN_LIST_BIT);
        }

        /**
         * \brief Tests whether a tetrahedron is marked.
         * \details A tetrahedron is marked whenever it is
         *  detected as non-conflict. The index of the
         *  point being inserted is used as a time-stamp
         *  for marking tetrahedra. The same space is used
         *  for marking and for chaining the conflict list.
         *  A tetrahedron can be in the following states:
         *  - in list
         *  - not in list and marked
         *  - not in list and not marked
         * \param[in] t index of the tetrahedron
         * \retval true if tetrahedron \p t is marked
         * \retval false otherwise
         */
        bool tet_is_marked(index_t t) const {
            return cell_next_[t] == cur_stamp_;
        }

        /**
         * \brief Marks a tetrahedron.
         * \details A tetrahedron is marked whenever it is
         *  detected as non-conflict. The same space is used
         *  for marking and for chaining the conflict list.
         *  The index of the point being inserted is used as a 
         *  time-stamp for marking tetrahedra.
         *  A tetrahedron can be in the following states:
         *  - in list
         *  - not in list and marked
         *  - not in list and not marked
         * \param[in] t index of the tetrahedron to be marked
         */
        void mark_tet(index_t t) {
            cell_next_[t] = cur_stamp_;
        }

        // _________ Combinatorics ___________________________________

        /**
         * \brief Returns the local index of a vertex by 
         *   facet and by local vertex index in the facet.
         * \details
         * tet facet vertex is such that the tetrahedron
         * formed with:
         * - vertex lv
         * - tet_facet_vertex(lv,0)
         * - tet_facet_vertex(lv,1)
         * - tet_facet_vertex(lv,2)
         * has the same orientation as the original tetrahedron for
         * any vertex lv.
         * \param[in] f local facet index, in (0,1,2,3)
         * \param[in] v local vertex index, in (0,1,2)
         * \return the local tetrahedron vertex index of 
         *  vertex \p v in facet \p f
         */
        static index_t tet_facet_vertex(index_t f, index_t v) {
            geo_debug_assert(f < 4);
            geo_debug_assert(v < 3);
            return index_t(tet_facet_vertex_[f][v]);
        }

        /**
         * \brief Gets the index of a vertex of a tetrahedron
         * \param[in] t index of the tetrahedron
         * \param[in] lv local vertex (0,1,2 or 3) index in \p t
         * \return the global index of the \p lv%th vertex of tetrahedron \p t
         *  or -1 if the vertex is at infinity
         */
        signed_index_t tet_vertex(index_t t, index_t lv) const {
            geo_debug_assert(t < max_t());
            geo_debug_assert(lv < 4);
            return cell_to_v_store_[4 * t + lv];
        }

        /**
         * \brief Finds the index of the vertex in a tetrahedron.
         * \param[in] t the tetrahedron
         * \param[in] v the vertex
         * \return iv such that tet_vertex(t,v)==iv
         * \pre \p t is incident to \p v
         */
        index_t find_tet_vertex(index_t t, signed_index_t v) const {
            geo_debug_assert(t < max_t());
            //   Find local index of v in tetrahedron t vertices.
            const signed_index_t* T = &(cell_to_v_store_[4 * t]);
            return find_4(T,v);
        }


        /**
         * \brief Gets the index of a vertex of a tetrahedron
         * \param[in] t index of the tetrahedron
         * \param[in] lv local vertex (0,1,2 or 3) index in \p t
         * \return the global index of the \p lv%th vertex of tetrahedron \p t
         * \pre Vertex \p lv of tetrahedron \p t is not at infinity
         */
         index_t finite_tet_vertex(index_t t, index_t lv) const {
            geo_debug_assert(t < max_t());
            geo_debug_assert(lv < 4);
            geo_debug_assert(cell_to_v_store_[4 * t + lv] != -1);
            return index_t(cell_to_v_store_[4 * t + lv]);
        }

        /**
         * \brief Sets a tetrahedron-to-vertex adjacency.
         * \param[in] t index of the tetrahedron
         * \param[in] lv local vertex index (0,1,2 or 3) in \p t
         * \param[in] v global index of the vertex
         */
        void set_tet_vertex(index_t t, index_t lv, signed_index_t v) {
            geo_debug_assert(t < max_t());
            geo_debug_assert(lv < 4);
            cell_to_v_store_[4 * t + lv] = v;
        }

        /**
         * \brief Gets the index of a tetrahedron adjacent to another one.
         * \param[in] t index of the tetrahedron
         * \param[in] lf local facet (0,1,2 or 3) index in \p t
         * \return the tetrahedron adjacent to \p t accorss facet \p lf
         */
        signed_index_t tet_adjacent(index_t t, index_t lf) const {
            geo_debug_assert(t < max_t());
            geo_debug_assert(lf < 4);
            signed_index_t result = cell_to_cell_store_[4 * t + lf];
            return result;
        }

        /**
         * \brief Sets a tetrahedron-to-tetrahedron adjacency.
         * \param[in] t1 index of the first tetrahedron
         * \param[in] lf1 local facet index (0,1,2 or 3) in t1
         * \param[in] t2 index of the tetrahedron
         *  adjacent to \p t1 accross \p lf1
         */
        void set_tet_adjacent(index_t t1, index_t lf1, index_t t2) {
            geo_debug_assert(t1 < max_t());
            geo_debug_assert(t2 < max_t());
            geo_debug_assert(lf1 < 4);
            cell_to_cell_store_[4 * t1 + lf1] = signed_index_t(t2);
        }
        
        /**
         * \brief Finds the index of the facet accross which t1 is 
         *  adjacent to t2.
         * \param[in] t1 first tetrahedron
         * \param[in] t2 second tetrahedron
         * \return f such that tet_adjacent(t1,f)==t2
         * \pre \p t1 and \p t2 are adjacent
         */
        index_t find_tet_adjacent(
            index_t t1, index_t t2_in
        ) const {
            geo_debug_assert(t1 < max_t());
            geo_debug_assert(t2_in < max_t());
            geo_debug_assert(t1 != t2_in);

            signed_index_t t2 = signed_index_t(t2_in);

            // Find local index of t2 in tetrahedron t1 adajcent tets.
            const signed_index_t* T = &(cell_to_cell_store_[4 * t1]);
            index_t result = find_4(T,t2);

            // Sanity check: make sure that t1 is adjacent to t2
            // only once!
            geo_debug_assert(tet_adjacent(t1,(result+1)%4) != t2);
            geo_debug_assert(tet_adjacent(t1,(result+2)%4) != t2);
            geo_debug_assert(tet_adjacent(t1,(result+3)%4) != t2);
            return result;
        }



        /**
         * \brief Sets the vertices and adjacent tetrahedra of
         *  a tetrahedron.
         * \param[in] t index of the tetrahedron
         * \param[in] v0 index of the first vertex
         * \param[in] v1 index of the second vertex
         * \param[in] v2 index of the third vertex
         * \param[in] v3 index of the fourth vertex
         * \param[in] a0 index of the adjacent tetrahedron opposite to \p v0
         * \param[in] a1 index of the adjacent tetrahedron opposite to \p v1
         * \param[in] a2 index of the adjacent tetrahedron opposite to \p v2
         * \param[in] a3 index of the adjacent tetrahedron opposite to \p v3
         */
        void set_tet(
            index_t t,
            signed_index_t v0, signed_index_t v1,
            signed_index_t v2, signed_index_t v3,
            index_t a0, index_t a1, index_t a2, index_t a3
        ) {
            geo_debug_assert(t < max_t());
            cell_to_v_store_[4 * t] = v0;
            cell_to_v_store_[4 * t + 1] = v1;
            cell_to_v_store_[4 * t + 2] = v2;
            cell_to_v_store_[4 * t + 3] = v3;
            cell_to_cell_store_[4 * t] = signed_index_t(a0);
            cell_to_cell_store_[4 * t + 1] = signed_index_t(a1);
            cell_to_cell_store_[4 * t + 2] = signed_index_t(a2);
            cell_to_cell_store_[4 * t + 3] = signed_index_t(a3);
        }

        // _________ Combinatorics - traversals ______________________________

        /**
         *  Gets the local facet index incident to an
         * oriented halfedge.
         * \param[in] t index of the tetrahedron
         * \param[in] v1 global index of the first extremity
         * \param[in] v2 global index of the second extremity
         * \return the local index of the facet incident to
         *  the oriented edge \p v1, \p v2.
         */
        index_t get_facet_by_halfedge(
            index_t t, signed_index_t v1, signed_index_t v2
        ) const {
            geo_debug_assert(t < max_t());
            geo_debug_assert(v1 != v2);
            //   Find local index of v1 and v2 in tetrahedron t
            const signed_index_t* T = &(cell_to_v_store_[4 * t]);
            index_t lv1 = find_4(T,v1);
            index_t lv2 = find_4(T,v2);
            geo_debug_assert(lv1 != lv2);
            return index_t(halfedge_facet_[lv1][lv2]);
        }


        /**
         *  Gets the local facet indices incident to an
         * oriented halfedge.
         * \param[in] t index of the tetrahedron
         * \param[in] v1 global index of the first extremity
         * \param[in] v2 global index of the second extremity
         * \param[out] f12 the local index of the facet 
         *  indicent to the halfedge [v1,v2]
         * \param[out] f21 the local index of the facet 
         *  indicent to the halfedge [v2,v1]
         */
        void get_facets_by_halfedge(
            index_t t, signed_index_t v1, signed_index_t v2,
            index_t& f12, index_t& f21
        ) const {
            geo_debug_assert(t < max_t());
            geo_debug_assert(v1 != v2);

            //   Find local index of v1 and v2 in tetrahedron t
            // The following expression is 10% faster than using
            // if() statements (multiply by boolean result of test).
            // Thank to Laurent Alonso for this idea.
            const signed_index_t* T = &(cell_to_v_store_[4 * t]);

            signed_index_t lv1 = 
                (T[1] == v1) | ((T[2] == v1) * 2) | ((T[3] == v1) * 3);

            signed_index_t lv2 = 
                (T[1] == v2) | ((T[2] == v2) * 2) | ((T[3] == v2) * 3);

            geo_debug_assert(lv1 != 0 || T[0] == v1);
            geo_debug_assert(lv2 != 0 || T[0] == v2);
            geo_debug_assert(lv1 >= 0);
            geo_debug_assert(lv2 >= 0);
            geo_debug_assert(lv1 != lv2);

            f12 = index_t(halfedge_facet_[lv1][lv2]);
            f21 = index_t(halfedge_facet_[lv2][lv1]);
        }


        /**
         * \brief Gets the next tetrahedron around an oriented edge of
         *  a tetrahedron.
         * \param[in,out] t the tetrahedron
         * \param[in] v1 global index of the first extremity of the edge
         * \param[in] v2 global index of the second extremity of the edge
         * \return the next tetrahedron from \p t around the oriented edge
         *   (\p v1 \p v2).
         */
        index_t next_around_halfedge(
            index_t& t, signed_index_t v1, signed_index_t v2
        ) const {
            return (index_t)tet_adjacent(
                t, get_facet_by_halfedge(t, v1, v2)
            );
        }

        // _________ Predicates _____________________________________________

        /**
         * \brief Tests whether a given tetrahedron is in conflict with
         *  a given 3d point.
         * \details A real tetrahedron is in conflict with a point whenever
         *  the point is contained by its circumscribed sphere, and a
         *  virtual tetrahedron is in conflict with a point whenever the
         *  tetrahedron formed by its real face and with the point has
         *  positive orientation.
         * \param[in] t the index of the tetrahedron
         * \param[in] p a pointer to the coordinates of the point
         * \retval true if point \p p is in conflict with tetrahedron \p t
         * \retval false otherwise
         */
        bool tet_is_conflict(index_t t, const double* p) const {

            // Lookup tetrahedron vertices
            const double* pv[4];
            for(index_t i=0; i<4; ++i) {
                signed_index_t v = tet_vertex(t,i);
                pv[i] = (v == -1) ? nil : vertex_ptr(index_t(v));
            }

            // Check for virtual tetrahedra (then in_sphere()
            // is replaced with orient3d())
            for(index_t lf = 0; lf < 4; ++lf) {

                if(pv[lf] == nil) {

                    // Facet of a virtual tetrahedron opposite to
                    // infinite vertex corresponds to
                    // the triangle on the convex hull of the points.
                    // Orientation is obtained by replacing vertex lf
                    // with p.
                    pv[lf] = p;
                    Sign sign = PCK::orient_3d(pv[0],pv[1],pv[2],pv[3]);

                    if(sign > 0) {
                        return true;
                    }

                    if(sign < 0) {
                        return false;
                    }

                    // If sign is zero, we check the real tetrahedron
                    // adjacent to the facet on the convex hull.
                    geo_debug_assert(tet_adjacent(t, lf) >= 0);
                    index_t t2 = index_t(tet_adjacent(t, lf));
                    geo_debug_assert(!tet_is_virtual(t2));

                    //  If t2 is already chained in the conflict list,
                    // then it is conflict
                    if(tet_is_in_list(t2)) {
                        return true;
                    }

                    //  If t2 is marked, then it is not in conflict.
                    if(tet_is_marked(t2)) {
                        return false;
                    }

                    return tet_is_conflict(t2, p);
                }
            }

            //   If the tetrahedron is a finite one, it is in conflict
            // if its circumscribed sphere contains the point (this is
            // the standard case).

            if(weighted_) {
                double h0 = heights_[finite_tet_vertex(t, 0)];
                double h1 = heights_[finite_tet_vertex(t, 1)];
                double h2 = heights_[finite_tet_vertex(t, 2)];
                double h3 = heights_[finite_tet_vertex(t, 3)];
                index_t pindex = index_t(
                    (p - vertex_ptr(0)) / int(vertex_stride_)
                );
                double h = heights_[pindex];
                return (PCK::orient_4d_SOS(
                            pv[0],pv[1],pv[2],pv[3],p,h0,h1,h2,h3,h
                       ) > 0) ;
            }

            return (PCK::in_sphere_3d_SOS(pv[0], pv[1], pv[2], pv[3], p) > 0);
        }

    protected:

        /**
         * \brief Finds the index of an integer in an array of four integers.
         * \param[in] T a const pointer to an array of four integers
         * \param[in] v the integer to retreive in \p T
         * \return the index (0,1,2 or 3) of \p v in \p T
         * \pre The four entries of \p T are different and one of them is
         *  equal to \p v.
         */
        static index_t find_4(const signed_index_t* T, signed_index_t v) {
            // The following expression is 10% faster than using
            // if() statements. This uses the C++ norm, that 
            // ensures that the 'true' boolean value converted to 
            // an int is always 1. With most compilers, this avoids 
            // generating branching instructions.
            // Thank to Laurent Alonso for this idea.
            // Note: Laurent also has this version:
            //    (T[0] != v)+(T[2]==v)+2*(T[3]==v)
            // that avoids a *3 multiply, but it is not faster in
            // practice.
            index_t result = index_t(
                (T[1] == v) | ((T[2] == v) * 2) | ((T[3] == v) * 3)
            );
            // Sanity check, important if it was T[0], not explicitely
            // tested (detects input that does not meet the precondition).
            geo_debug_assert(T[result] == v);
            return result; 
        }



    private:
        vector<signed_index_t> cell_to_v_store_;
        vector<signed_index_t> cell_to_cell_store_;
        vector<index_t> cell_next_;
        vector<index_t> reorder_;
        index_t cur_stamp_; // used for marking
        index_t first_free_;
        bool weighted_;
        vector<double> heights_; // only used in weighted mode
        bool do_reorder_; // If true, uses BRIO reordering

        /**
         * \brief Gives the indexing of tetrahedron facet
         *  vertices.
         * \details tet_facet_vertex[lf][lv] gives the
         *  local vertex index (in 0,1,2,3) from a
         *  local facet index lf (in 0,1,2,3) and a
         *  local vertex index within the facet (in 0,1,2).
         */
        static char tet_facet_vertex_[4][3];

        /**
         * \brief Gives a local facet index by
         *  halfedge extremities local indices.
         */
        static char halfedge_facet_[4][4];
    };

    /************************************************************************/
    
}

#endif

