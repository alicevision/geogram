/*
 *  Copyright (c) 2010-2017, ALICE project, Inria
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

#include <geogram/voronoi/RVD_callback.h>
#include <geogram/voronoi/RVD_mesh_builder.h>
#include <geogram/voronoi/generic_RVD_cell.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/basic/argused.h>

namespace {
    using namespace GEO;

    /**
     * \brief Gets the maximum region index in the mesh, not counting
     *   index_t(-1) triangles.
     * \return The maximum region index, or index_t(-1) if all region
     *   indices are index_t(-1)
     */
    index_t max_region(Mesh& mesh, Attribute<index_t>& facet_region) {
	index_t result = index_t(-1);
	for(index_t f=0; f<mesh.facets.nb(); ++f) {
	    if(facet_region[f] != index_t(-1)) {
		if(result == index_t(-1)) {
		    result = facet_region[f];
		} else {
		    result = geo_max(result, facet_region[f]);
		}
	    }
	}
	return result;
    }

    /**
     * \brief Splits the regions along hard edges.
     * \param[in,out] mesh a reference to the surface mesh.
     * \param[in,out] facet_region the attribute that defines the regions. 
     *  The regions with facet_region == index_t(-1) will be split.
     * \param[in] threshold (in degrees). Edges with two adjacent facets 
     *  with normals that form an  angle larger than \p threshold degrees are 
     *  considered as hard edges.
     */
    void split_regions_along_hard_edges(
	Mesh& mesh, Attribute<index_t>& facet_region, double threshold
    ) {

	threshold *= (M_PI / 180.0);
	
	index_t R = max_region(mesh, facet_region);
	if(R == index_t(-1)) {
	    R = 0;
	} else {
	    ++R;
	}
	vector<bool> is_crease(mesh.facet_corners.nb(), false);
	for(index_t f=0; f<mesh.facets.nb(); ++f) {
	    vec3 n = Geom::mesh_facet_normal(mesh,f);
	    for(
		index_t c=mesh.facets.corners_begin(f);
		c<mesh.facets.corners_end(f); ++c
	    ) {
		index_t f2 = mesh.facet_corners.adjacent_facet(c);
		if(f2 != index_t(-1)) {
		    vec3 n2 = Geom::mesh_facet_normal(mesh,f2);
		    double alpha = Geom::angle(n,n2);
		    if(alpha > threshold) {
			is_crease[c] = true;
		    }
		}
	    }
	}
	for(index_t f=0; f<mesh.facets.nb(); ++f) {
	    if(facet_region[f] == index_t(-1)) {
		std::stack<index_t> S;
		facet_region[f] = R;
		S.push(f);
		while(!S.empty()) {
		    index_t f2 = S.top();
		    S.pop();
		    for(
			index_t c=mesh.facets.corners_begin(f2);
			c<mesh.facets.corners_end(f2); ++c
		    ) {
			if(is_crease[c]) {
			    continue;
			}
			index_t f3 = mesh.facet_corners.adjacent_facet(c);
			if(
			    f3 != index_t(-1) &&
			    facet_region[f3] == index_t(-1)
			) {
			    facet_region[f3] = R;
			    S.push(f3);
			}
		    }
		}
		++R;
	    }
	}
    }
    
    /**
     * \brief Simplifies the facets of a surface mesh based on an attribute.
     * \details Groups of connected facets with the same attribute value are
     *  replaced with a single facet.
     * \param[in,out] mesh a reference to the surface mesh to be simplified
     * \param[in] facet_region a reference to the attribute with the facet 
     *  region
     * \param[in] angle_threshold (in degrees). In the outer region 
     *  (i.e. facet_region == index_t(-1)), an edge shared by two adjacent 
     *  facets is suppressed if the angle between the 
     *  two facet normals is smaller than \p angle_threshold.
     */
    bool simplify(
	Mesh& mesh,
	Attribute<index_t>& facet_region,
	double angle_threshold
    ) {

	bool keep_outer_region = (angle_threshold == 0.0);

	index_t max_r = 0;
	
	if(!keep_outer_region) {
	    max_r = max_region(mesh, facet_region);
	    split_regions_along_hard_edges(mesh, facet_region, angle_threshold);
	}
	
	vector<bool> is_corner(mesh.vertices.nb(),false);

	{
	    vector<index_t> rgn1(mesh.vertices.nb(), index_t(-2));
	    vector<index_t> rgn2(mesh.vertices.nb(), index_t(-2));
	    
	    // Keep all the vertices adjacent to 3 regions or more
	    for(index_t f=0; f<mesh.facets.nb(); ++f) {
		index_t r = facet_region[f];
		for(index_t lv=0; lv<mesh.facets.nb_vertices(f); ++lv) {
		    index_t v = mesh.facets.vertex(f,lv);
		    if(rgn1[v] == r || rgn2[v] == r) {
			continue;
		    }
		    if(rgn1[v] == index_t(-2)) {
			rgn1[v] = r;
		    } else if(rgn2[v] == index_t(-2)) {
			rgn2[v] = r;
		    } else {
			is_corner[v] = true;
		    }
		}
	    }

	    // Keep also the vertices adjacent to region -1 and
	    // to another region
	    
	    if(keep_outer_region) {
		for(index_t v=0; v<mesh.vertices.nb(); ++v) {
		    if(
			(rgn1[v] == index_t(-1) && rgn2[v] != index_t(-1)) ||
			(rgn2[v] == index_t(-1) && rgn1[v] != index_t(-1)) 
		    ) {
			is_corner[v] = true;
		    }
		}
	    }
	}

	// -1: not visited, 0: keep, 1: delete
	vector<index_t> facet_status(mesh.facets.nb(),index_t(-1));

	// Needs to be backed-up, we are modifying the mesh !!
	index_t nf = mesh.facets.nb();
	
	for(index_t f=0; f<nf; ++f) {
	    if(facet_status[f] == index_t(-1)) {
		index_t r = facet_region[f];
		if(keep_outer_region && (r == index_t(-1))) {
		    facet_status[f] = 0;
		    continue;
		} else {
		    std::stack<index_t> S;
		    std::map<index_t, index_t> border_next;
		    S.push(f);
		    facet_status[f] = 1;
		    while(!S.empty()) {
		    	index_t f2 = S.top();
			S.pop();
			for(index_t c1=mesh.facets.corners_begin(f2);
			    c1<mesh.facets.corners_end(f2); ++c1) {
			    index_t f3 = mesh.facet_corners.adjacent_facet(c1);
			    if(f3 == index_t(-1) || facet_region[f3] != r) {
				index_t c2 =
				    mesh.facets.next_corner_around_facet(f2,c1);
				index_t v1 = mesh.facet_corners.vertex(c1);
				index_t v2 = mesh.facet_corners.vertex(c2);
				if(border_next.find(v1) != border_next.end()) {
				    Logger::warn("Simplify")
					<< "Region has non-manifold border"
					<< std::endl;
			    // Yes, goto!, why not, what's wrong with goto ?
			    // Why would it be ok to throw() and not goto ?
				    goto rollback;
				}
				border_next[v1] = v2;
			    } else {
				if(facet_status[f3] == index_t(-1)) {
				    facet_status[f3] = 1;
				    S.push(f3);
				}
			    }
			}
		    }
		    index_t nb_border_visited = 0;
		    vector<index_t> new_facet;
		    index_t v = border_next.begin()->first;
		    do {
			if(is_corner[v]) {
			    new_facet.push_back(v);
			}
			++nb_border_visited;
			v = border_next[v];
			if(nb_border_visited > mesh.vertices.nb()) {
			    Logger::warn("Simplify")
				<< "Region has singular border topology"
				<< std::endl;
			    goto rollback;
			}
		    } while(v != border_next.begin()->first);

		    if(nb_border_visited != border_next.size()) {
			Logger::warn("Simplify")
			    << "Region has multiple borders"
			    << std::endl;
			goto rollback;
		    }
		    
		    if(new_facet.size() < 3) {
			Logger::warn("Simplify")
			    << "Region has border with less than 3 corners"
			    << std::endl;
			goto rollback;
		    }

		    index_t new_f = mesh.facets.nb();
		    mesh.facets.create_polygon(new_facet.size());
		    for(index_t i=0; i<new_facet.size(); ++i) {
			mesh.facets.set_vertex(new_f,i,new_facet[i]);
		    }
		    facet_region[new_f] = r;
		}
	    }
	}
	
	facet_status.resize(mesh.facets.nb(), 0);
	mesh.facets.delete_elements(facet_status);
	if(!keep_outer_region) {
	    for(index_t f=0; f<mesh.facets.nb(); ++f) {
		if(facet_region[f] > max_r) {
		    facet_region[f] = index_t(-1);
		}
	    }
	}
	return true;

      rollback:
	Logger::out("Simplify") << "...Rolling back." << std::endl;
	// delete all facets...
	facet_status.resize(mesh.facets.nb(), 1);
	// ... except the initial ones !
	for(index_t f=0; f<nf; ++f) {
	    facet_status[f] = 0;
	}
	mesh.facets.delete_elements(facet_status);
	if(!keep_outer_region) {
	    for(index_t f=0; f<mesh.facets.nb(); ++f) {
		if(facet_region[f] > max_r) {
		    facet_region[f] = index_t(-1);
		}
	    }
	}
	return false;
    }
}


namespace GEO {

    RVDCallback::RVDCallback() :
	seed_(index_t(-1)),
	simplex_(index_t(-1)),
	spinlocks_(nil) {
    }

    RVDCallback::~RVDCallback() {
    }

    void RVDCallback::begin() {
	seed_ = index_t(-1);
	simplex_ = index_t(-1);
    }

    void RVDCallback::end() {
    }

    /*********************************************************************/

    RVDPolygonCallback::RVDPolygonCallback() {
    }
    
    RVDPolygonCallback::~RVDPolygonCallback() {
    }

    void RVDPolygonCallback::operator() (
	index_t v,
	index_t t,
	const GEOGen::Polygon& C
    ) const {
	const_cast<RVDPolygonCallback*>(this)->seed_ = v;
	const_cast<RVDPolygonCallback*>(this)->simplex_ = t;
	geo_argused(C);
    }
    
    void RVDPolygonCallback::begin() {
    }
    
    void RVDPolygonCallback::end() {
    }
    
    /*********************************************************************/
    
    RVDPolyhedronCallback::RVDPolyhedronCallback() :
	facet_seed_(index_t(-1)),
	facet_tet_(index_t(-1)),
	last_seed_(index_t(-1)),
	simplify_internal_tet_facets_(false),
	simplify_voronoi_facets_(false),
	simplify_boundary_facets_(false),
	simplify_boundary_facets_angle_threshold_(0.0),
	use_mesh_(false),	
	facet_is_skipped_(false),
	vertex_map_(nil)
    {
    }

    RVDPolyhedronCallback::~RVDPolyhedronCallback() {
	if(mesh_vertex_sym_.is_bound()) {
	    mesh_vertex_sym_.unbind();
	}
	if(mesh_facet_seed_.is_bound()) {
	    mesh_facet_seed_.unbind();
	}
	if(mesh_facet_tet_.is_bound()) {
	    mesh_facet_tet_.unbind();
	}
    }

    void RVDPolyhedronCallback::set_use_mesh(bool x) {
	use_mesh_ = x;
	if(!mesh_vertex_sym_.is_bound()) {
	    mesh_vertex_sym_.bind(mesh_.vertices.attributes(), "sym");
	}
	if(!mesh_facet_seed_.is_bound()) {
	    mesh_facet_seed_.bind(mesh_.facets.attributes(),"seed");
	}
	if(!mesh_facet_tet_.is_bound()) {
	    mesh_facet_tet_.bind(mesh_.facets.attributes(),"tet");
	}
    }

    /********************************************************************/    
    
    void RVDPolyhedronCallback::begin_polyhedron(
	index_t seed, index_t tetrahedron
    ) {
	geo_argused(seed);
	geo_argused(tetrahedron);
    }

    void RVDPolyhedronCallback::begin_facet(
	index_t facet_seed, index_t facet_tet
    ) {
	geo_argused(facet_seed);
	geo_argused(facet_tet);
    }

    void RVDPolyhedronCallback::vertex(
	const double* geometry, const GEOGen::SymbolicVertex& symb
    ) {
	geo_argused(geometry);
	geo_argused(symb);
    }

    void RVDPolyhedronCallback::end_facet() {
    }

    void RVDPolyhedronCallback::end_polyhedron() {
    }

    /********************************************************************/

    void RVDPolyhedronCallback::begin_polyhedron_internal(
	index_t seed, index_t tetrahedron
    ) {
	last_seed_ = seed;
	seed_ = seed;
	simplex_ = tetrahedron;
	if(use_mesh_) {
	    vertex_map_ = new RVDVertexMap;
	} else {
	    begin_polyhedron(seed, tetrahedron);
	}
    }

    void RVDPolyhedronCallback::begin_facet_internal(
	index_t facet_seed, index_t facet_tet
    ) {
	facet_seed_ = facet_seed;
	facet_tet_ = facet_tet;
	facet_is_skipped_ = (
	    simplify_internal_tet_facets_ && (facet_tet != index_t(-1))
	);
	if(!facet_is_skipped_) {
	    if(use_mesh_) {
	    } else {
		begin_facet(facet_seed, facet_tet);
	    }
	}
    }

    void RVDPolyhedronCallback::vertex_internal(
	const double* geometry, const GEOGen::SymbolicVertex& symb
    ) {
	if(!facet_is_skipped_) {	
	    if(use_mesh_) {
		index_t v = vertex_map_->find_or_create_vertex(seed(),symb);
		if(v >= mesh_.vertices.nb()) {
		    mesh_.vertices.create_vertex(geometry);
		    mesh_vertex_sym_[v] = symb;
		}
		current_facet_.push_back(v);
	    } else {
		vertex(geometry, symb);
	    }
	}
    }

    void RVDPolyhedronCallback::end_facet_internal() {
	if(!facet_is_skipped_) {	
	    if(use_mesh_) {
		index_t f = mesh_.facets.nb();
		mesh_.facets.create_polygon(current_facet_.size());
		for(index_t i=0; i<current_facet_.size(); ++i) {
		    mesh_.facets.set_vertex(f, i, current_facet_[i]);
		}
		mesh_facet_seed_[f] = facet_seed();
		mesh_facet_tet_[f] = facet_tet();
		current_facet_.resize(0);
	    } else {
		end_facet();
	    }
	}
	facet_seed_ = index_t(-1);
	facet_tet_ = index_t(-1);
    }

    void RVDPolyhedronCallback::end_polyhedron_internal() {
	if(use_mesh_) {
	    mesh_.facets.connect();
	    process_polyhedron_mesh();
	    mesh_.clear(true,true);
	    delete vertex_map_;
	    vertex_map_ = nil;
	} else {
	    end_polyhedron();
	}
	seed_ = index_t(-1);
	simplex_ = index_t(-1);
    }
    
    /********************************************************************/

    void RVDPolyhedronCallback::process_polyhedron_mesh() {
	if(simplify_voronoi_facets_) {
	    simplify(
		mesh_,
		mesh_facet_seed_,
		simplify_boundary_facets_angle_threshold_
	    );
	}
	begin_polyhedron(seed(), tet());
	for(index_t f=0; f<mesh_.facets.nb(); ++f) {
	    facet_seed_ = mesh_facet_seed_[f];
	    facet_tet_ = mesh_facet_tet_[f];
	    begin_facet(facet_seed_, facet_tet_);
	    for(index_t lv=0; lv<mesh_.facets.nb_vertices(f); ++lv) {
		index_t v = mesh_.facets.vertex(f,lv);
		vertex(mesh_.vertices.point_ptr(v), mesh_vertex_sym_[v]);
	    }
	    end_facet();
	}
	end_polyhedron();
    }
    
    /********************************************************************/

    void RVDPolyhedronCallback::begin() {
    }

    void RVDPolyhedronCallback::end() {
	
	GEO::RVDPolyhedronCallback& callbacks =
	    const_cast<GEO::RVDPolyhedronCallback&>(*this);
	
	if(simplify_internal_tet_facets_ && seed_ != index_t(-1)) {
	    callbacks.end_polyhedron_internal();	    
	}
    }
    
    void RVDPolyhedronCallback::operator() (
	index_t v,
	index_t t,
	const GEOGen::ConvexCell& C
    ) const {

	GEO::RVDPolyhedronCallback& callbacks =
	    const_cast<GEO::RVDPolyhedronCallback&>(*this);

	if(simplify_internal_tet_facets_) {
	    if(v != last_seed_) {
		if(last_seed_ != index_t(-1)) {
		    callbacks.end_polyhedron_internal();
		}
		callbacks.begin_polyhedron_internal(v,t);
	    }
	} else {
	    callbacks.begin_polyhedron_internal(v,t);	    
	}

	// Remember that ConvexCell is represented in dual form !
	//   - ConvexCell's vertices are facets
	//   - ConvexCell's triangles are vertices
	
	for(index_t cv = 0; cv < C.max_v(); ++cv) {
	    signed_index_t ct = C.vertex_triangle(cv);
	    if(ct == -1) {
		continue;
	    }
	    geo_debug_assert(C.triangle_is_used(index_t(ct)));
	    
	    signed_index_t adjacent = C.vertex_id(cv);
	    signed_index_t v_adj = -1;
	    signed_index_t t_adj = -1;

	    if(adjacent < 0) {
		// Negative adjacent indices correspond to
		// tet-tet links
		t_adj = -adjacent - 1;
	    } else if(adjacent > 0) {
		// Positive adjacent indices correspond to
		// Voronoi seed - Voronoi seed link
		v_adj = adjacent - 1;
	    } // Zero adjacent indices corresponds to
 	      // tet facet on border.
  	    
	    callbacks.begin_facet_internal(index_t(v_adj), index_t(t_adj));
	    
	    GEOGen::ConvexCell::Corner first(
		index_t(ct), C.find_triangle_vertex(index_t(ct), cv)
	    );
	    
	    GEOGen::ConvexCell::Corner c = first;
	    do {
		const GEOGen::Vertex& vx = C.triangle_dual(c.t);
		callbacks.vertex_internal(vx.point(), vx.sym());
		C.move_to_next_around_vertex(c);
	    } while(c != first);
	    callbacks.end_facet_internal();
	}

	if(!simplify_internal_tet_facets_) {
	    callbacks.end_polyhedron_internal();
	}
    }
    
}

