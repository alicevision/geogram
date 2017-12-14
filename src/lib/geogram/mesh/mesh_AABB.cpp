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

#include <geogram/mesh/mesh_AABB.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/numerics/predicates.h>
#include <geogram/basic/geometry_nd.h>

namespace {

    using namespace GEO;

    /**
     * \brief Computes the axis-aligned bounding box of a mesh facet.
     * \param[in] M the mesh
     * \param[out] B the bounding box of the facet
     * \param[in] f the index of the facet in mesh \p M
     */
    void get_facet_bbox(
        const Mesh& M, Box& B, index_t f
    ) {
        index_t c = M.facets.corners_begin(f);
        const double* p = M.vertices.point_ptr(M.facet_corners.vertex(c));
        for(coord_index_t coord = 0; coord < 3; ++coord) {
            B.xyz_min[coord] = p[coord];
            B.xyz_max[coord] = p[coord];
        }
        for(++c; c < M.facets.corners_end(f); ++c) {
            p = M.vertices.point_ptr(M.facet_corners.vertex(c));
            for(coord_index_t coord = 0; coord < 3; ++coord) {
                B.xyz_min[coord] = geo_min(B.xyz_min[coord], p[coord]);
                B.xyz_max[coord] = geo_max(B.xyz_max[coord], p[coord]);
            }
        }
    }


    /**
     * \brief Computes the axis-aligned bounding box of a mesh tetrahedron.
     * \param[in] M the mesh
     * \param[out] B the bounding box of the facet
     * \param[in] t the index of the tetrahedron in mesh \p M
     */
    void get_tet_bbox(
        const Mesh& M, Box& B, index_t t
    ) {
        const double* p = M.vertices.point_ptr(M.cells.vertex(t,0));
        for(coord_index_t coord = 0; coord < 3; ++coord) {
            B.xyz_min[coord] = p[coord];
            B.xyz_max[coord] = p[coord];
        }
        for(index_t lv=1; lv<4; ++lv) {
            p = M.vertices.point_ptr(M.cells.vertex(t,lv));
            for(coord_index_t coord = 0; coord < 3; ++coord) {
                B.xyz_min[coord] = geo_min(B.xyz_min[coord], p[coord]);
                B.xyz_max[coord] = geo_max(B.xyz_max[coord], p[coord]);
            }
        }
    }

    /**
     * \brief Computes the axis-aligned bounding box of a mesh cell
     * \param[in] M the mesh
     * \param[out] B the bounding box of the facet
     * \param[in] c the index of the cell in mesh \p M
     */
    void get_cell_bbox(
        const Mesh& M, Box& B, index_t c
    ) {
        const double* p = M.vertices.point_ptr(M.cells.vertex(c,0));
        for(coord_index_t coord = 0; coord < 3; ++coord) {
            B.xyz_min[coord] = p[coord];
            B.xyz_max[coord] = p[coord];
        }
        for(index_t lv=1; lv<M.cells.nb_vertices(c); ++lv) {
            p = M.vertices.point_ptr(M.cells.vertex(c,lv));
            for(coord_index_t coord = 0; coord < 3; ++coord) {
                B.xyz_min[coord] = geo_min(B.xyz_min[coord], p[coord]);
                B.xyz_max[coord] = geo_max(B.xyz_max[coord], p[coord]);
            }
        }
    }
    
    /**
     * \brief Computes the maximum node index in a subtree
     * \param[in] node_index node index of the root of the subtree
     * \param[in] b first facet index in the subtree
     * \param[in] e one position past the last facet index in the subtree
     * \return the maximum node index in the subtree rooted at \p node_index
     */
    index_t max_node_index(index_t node_index, index_t b, index_t e) {
        geo_debug_assert(e > b);
        if(b + 1 == e) {
            return node_index;
        }
        index_t m = b + (e - b) / 2;
        index_t childl = 2 * node_index;
        index_t childr = 2 * node_index + 1;
        return geo_max(
            max_node_index(childl, b, m),
            max_node_index(childr, m, e)
        );
    }

    /**
     * \brief Computes the hiearchy of bounding boxes recursively.
     * \details This function is generic and can be used to compute
     *  a bbox hierarchy of arbitrary elements.
     * \param[in] M the mesh
     * \param[in] bboxes the array of bounding boxes
     * \param[in] node_index the index of the root of the subtree
     * \param[in] b first element index in the subtree
     * \param[in] e one position past the last element index in the subtree
     * \param[in] get_bbox a function that computes the bbox of an element
     * \tparam GET_BBOX a function (or a functor) with the following arguments:
     *  - mesh: a const reference to the mesh
     *  - box: a reference where the computed bounding box of the element 
     *   will be stored
     *  - element: the index of the element
     */
    template <class GET_BBOX>
    void init_bboxes_recursive(
        const Mesh& M, vector<Box>& bboxes,
        index_t node_index,
        index_t b, index_t e,
        const GET_BBOX& get_bbox
    ) {
        geo_debug_assert(node_index < bboxes.size());
        geo_debug_assert(b != e);
        if(b + 1 == e) {
            get_bbox(M, bboxes[node_index], b);
            return;
        }
        index_t m = b + (e - b) / 2;
        index_t childl = 2 * node_index;
        index_t childr = 2 * node_index + 1;
        geo_debug_assert(childl < bboxes.size());
        geo_debug_assert(childr < bboxes.size());
        init_bboxes_recursive(M, bboxes, childl, b, m, get_bbox);
        init_bboxes_recursive(M, bboxes, childr, m, e, get_bbox);
        geo_debug_assert(childl < bboxes.size());
        geo_debug_assert(childr < bboxes.size());
        bbox_union(bboxes[node_index], bboxes[childl], bboxes[childr]);
    }

    /**
     * \brief Finds the nearest point in a mesh facet from a query point.
     * \param[in] M the mesh
     * \param[in] p the query point
     * \param[in] f index of the facet in \p M
     * \param[out] nearest_p the point of facet \p f nearest to \p p
     * \param[out] squared_dist the squared distance between
     *  \p p and \p nearest_p
     * \pre the mesh \p M is triangulated
     */
    void get_point_facet_nearest_point(
        const Mesh& M,
        const vec3& p,
        index_t f,
        vec3& nearest_p,
        double& squared_dist
    ) {
        geo_debug_assert(M.facets.nb_vertices(f) == 3);
        index_t c = M.facets.corners_begin(f);
        const vec3& p1 = Geom::mesh_vertex(M, M.facet_corners.vertex(c));
        ++c;
        const vec3& p2 = Geom::mesh_vertex(M, M.facet_corners.vertex(c));
        ++c;
        const vec3& p3 = Geom::mesh_vertex(M, M.facet_corners.vertex(c));
        double lambda1, lambda2, lambda3;  // barycentric coords, not used.
        squared_dist = Geom::point_triangle_squared_distance(
            p, p1, p2, p3, nearest_p, lambda1, lambda2, lambda3
        );
    }

    /**
     * \brief Computes the squared distance between a point and a Box.
     * \param[in] p the point
     * \param[in] B the box
     * \return the squared distance between \p p and \p B
     * \pre p is inside B
     */
    double inner_point_box_squared_distance(
        const vec3& p,
        const Box& B
    ) {
        geo_debug_assert(B.contains(p));
        double result = geo_sqr(p[0] - B.xyz_min[0]);
        result = geo_min(result, geo_sqr(p[0] - B.xyz_max[0]));
        for(coord_index_t c = 1; c < 3; ++c) {
            result = geo_min(result, geo_sqr(p[c] - B.xyz_min[c]));
            result = geo_min(result, geo_sqr(p[c] - B.xyz_max[c]));
        }
        return result;
    }

    /**
     * \brief Computes the squared distance between a point and a Box
     *  with negative sign if the point is inside the Box.
     * \param[in] p the point
     * \param[in] B the box
     * \return the signed squared distance between \p p and \p B
     */
    double point_box_signed_squared_distance(
        const vec3& p,
        const Box& B
    ) {
        bool inside = true;
        double result = 0.0;
        for(coord_index_t c = 0; c < 3; c++) {
            if(p[c] < B.xyz_min[c]) {
                inside = false;
                result += geo_sqr(p[c] - B.xyz_min[c]);
            } else if(p[c] > B.xyz_max[c]) {
                inside = false;
                result += geo_sqr(p[c] - B.xyz_max[c]);
            }
        }
        if(inside) {
            result = -inner_point_box_squared_distance(p, B);
        }
        return result;
    }

    /**
     * \brief Computes the squared distance between a point and the
     *  center of a box.
     * \param[in] p the point
     * \param[in] B the box
     * \return the squared distance between \p p and the center of \p B
     */
    double point_box_center_squared_distance(
        const vec3& p, const Box& B
    ) {
        double result = 0.0;
        for(coord_index_t c = 0; c < 3; ++c) {
            double d = p[c] - 0.5 * (B.xyz_min[c] + B.xyz_max[c]);
            result += geo_sqr(d);
        }
        return result;
    }

    /**
     * \brief Tests whether a mesh tetrahedron contains a given point
     * \param[in] M a const reference to the mesh
     * \param[in] t the index of the tetrahedron in \p M
     * \param[in] p a const reference to the point
     * \param[in] exact specifies whether exact predicates should be used
     * \retval true if the tetrahedron \p t or its boundary contains 
     *  the point \p p
     * \retval false otherwise
     */
    bool mesh_tet_contains_point(
        const Mesh& M, index_t t, const vec3& p, bool exact = true
    ) {
        // Inexact mode is not implemented yet.
        geo_argused(exact);

        const vec3& p0 = Geom::mesh_vertex(M, M.cells.vertex(t,0));
        const vec3& p1 = Geom::mesh_vertex(M, M.cells.vertex(t,1));
        const vec3& p2 = Geom::mesh_vertex(M, M.cells.vertex(t,2));
        const vec3& p3 = Geom::mesh_vertex(M, M.cells.vertex(t,3));

        Sign s[4];
        s[0] = PCK::orient_3d(p, p1, p2, p3);
        s[1] = PCK::orient_3d(p0, p, p2, p3);
        s[2] = PCK::orient_3d(p0, p1, p, p3);
        s[3] = PCK::orient_3d(p0, p1, p2, p);

        return (
            (s[0] >= 0 && s[1] >= 0 && s[2] >= 0 && s[3] >= 0) ||
            (s[0] <= 0 && s[1] <= 0 && s[2] <= 0 && s[3] <= 0)
        );
    }
}

/****************************************************************************/

namespace GEO {

    MeshFacetsAABB::MeshFacetsAABB(
        Mesh& M, bool reorder
    ) :
        mesh_(M) {
        if(!M.facets.are_simplices()) {
            mesh_repair(
		M,
		MeshRepairMode(
		    MESH_REPAIR_TRIANGULATE | MESH_REPAIR_QUIET
		 )
	    );
        }
        if(reorder) {
            mesh_reorder(mesh_, MESH_ORDER_MORTON);
        }
        bboxes_.resize(
            max_node_index(
                1, 0, mesh_.facets.nb()
            ) + 1 // <-- this is because size == max_index + 1 !!!
        );
        init_bboxes_recursive(
            mesh_, bboxes_, 1, 0, mesh_.facets.nb(), get_facet_bbox
        );
    }

    void MeshFacetsAABB::get_nearest_facet_hint(
        const vec3& p,
        index_t& nearest_f, vec3& nearest_point, double& sq_dist
    ) const {

        // Find a good initial value for nearest_f by traversing
        // the boxes and selecting the child such that the center
        // of its bounding box is nearer to the query point.
        // For a large mesh (20M facets) this gains up to 10%
        // performance as compared to picking nearest_f randomly.
        index_t b = 0;
        index_t e = mesh_.facets.nb();
        index_t n = 1;
        while(e != b + 1) {
            index_t m = b + (e - b) / 2;
            index_t childl = 2 * n;
            index_t childr = 2 * n + 1;
            if(
                point_box_center_squared_distance(p, bboxes_[childl]) <
                point_box_center_squared_distance(p, bboxes_[childr])
            ) {
                e = m;
                n = childl;
            } else {
                b = m;
                n = childr;
            }
        }
        nearest_f = b;

        index_t v = mesh_.facet_corners.vertex(
            mesh_.facets.corners_begin(nearest_f)
        );
        nearest_point = Geom::mesh_vertex(mesh_, v);
        sq_dist = Geom::distance2(p, nearest_point);
    }

    void MeshFacetsAABB::nearest_facet_recursive(
        const vec3& p,
        index_t& nearest_f, vec3& nearest_point, double& sq_dist,
        index_t n, index_t b, index_t e
    ) const {
        geo_debug_assert(e > b);

        // If node is a leaf: compute point-facet distance
        // and replace current if nearer
        if(b + 1 == e) {
            vec3 cur_nearest_point;
            double cur_sq_dist;
            get_point_facet_nearest_point(
                mesh_, p, b, cur_nearest_point, cur_sq_dist
            );
            if(cur_sq_dist < sq_dist) {
                nearest_f = b;
                nearest_point = cur_nearest_point;
                sq_dist = cur_sq_dist;
            }
            return;
        }
        index_t m = b + (e - b) / 2;
        index_t childl = 2 * n;
        index_t childr = 2 * n + 1;

        double dl = point_box_signed_squared_distance(p, bboxes_[childl]);
        double dr = point_box_signed_squared_distance(p, bboxes_[childr]);

        // Traverse the "nearest" child first, so that it has more chances
        // to prune the traversal of the other child.
        if(dl < dr) {
            if(dl < sq_dist) {
                nearest_facet_recursive(
                    p,
                    nearest_f, nearest_point, sq_dist,
                    childl, b, m
                );
            }
            if(dr < sq_dist) {
                nearest_facet_recursive(
                    p,
                    nearest_f, nearest_point, sq_dist,
                    childr, m, e
                );
            }
        } else {
            if(dr < sq_dist) {
                nearest_facet_recursive(
                    p,
                    nearest_f, nearest_point, sq_dist,
                    childr, m, e
                );
            }
            if(dl < sq_dist) {
                nearest_facet_recursive(
                    p,
                    nearest_f, nearest_point, sq_dist,
                    childl, b, m
                );
            }
        }
    }

/****************************************************************************/

    MeshCellsAABB::MeshCellsAABB(Mesh& M, bool reorder) : mesh_(M) {
        if(reorder) {
            mesh_reorder(mesh_, MESH_ORDER_MORTON);
        }
        bboxes_.resize(
            max_node_index(
                1, 0, mesh_.cells.nb()
            ) + 1 // <-- this is because size == max_index + 1 !!!
        );
        if(mesh_.cells.are_simplices()) {
            init_bboxes_recursive(
                mesh_, bboxes_, 1, 0, mesh_.cells.nb(), get_tet_bbox
            );
        } else {
            init_bboxes_recursive(
                mesh_, bboxes_, 1, 0, mesh_.cells.nb(), get_cell_bbox
            );
        }
    }

    index_t MeshCellsAABB::containing_tet_recursive(
        const vec3& p, bool exact,
        index_t n, index_t b, index_t e        
    ) const {

        if(!bboxes_[n].contains(p)) {
            return NO_TET;
        }
        
        if(e==b+1) {
            if(mesh_tet_contains_point(mesh_, b, p, exact)) {
                return b;
            } else {
                return NO_TET;
            }
        }
        
        index_t m = b + (e - b) / 2;
        index_t childl = 2 * n;
        index_t childr = 2 * n + 1;

        index_t result = containing_tet_recursive(
            p, exact, childl, b, m
        );
        if(result == NO_TET) {
            result = containing_tet_recursive(p, exact, childr, m, e);
        }
        return result;
    }
    
/****************************************************************************/
        
}

