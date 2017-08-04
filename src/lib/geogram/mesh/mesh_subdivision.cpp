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

#include <geogram/mesh/mesh_subdivision.h>
#include <geogram/mesh/mesh.h>

namespace GEO {

    void mesh_split_triangles(
	Mesh& M, index_t facets_begin, index_t facets_end
    ) {
	
	geo_assert(M.facets.are_simplices());

	if(facets_end == index_t(-1)) {
	    facets_end = M.facets.nb();
	}
	
	index_t nv0 = M.vertices.nb();
	index_t nf0 = M.facets.nb();

	// Compute corner to new vertex mapping
	vector<index_t> ctov(M.facet_corners.nb(), NO_VERTEX);
	index_t nbnewv=0;
	for(index_t f=facets_begin; f<facets_end; ++f) {
	    for(index_t c=M.facets.corners_begin(f);
		c<M.facets.corners_end(f); ++c
	    ) {
		if(ctov[c] == index_t(-1)) {
		    ctov[c] = nbnewv;
		    index_t f2 = M.facet_corners.adjacent_facet(c);
		    if(f2 != NO_FACET) {
			for(index_t c2=M.facets.corners_begin(f2);
			    c2!=M.facets.corners_end(f2); ++c2) {
			    if(M.facet_corners.adjacent_facet(c2) == f) {
				ctov[c2] = nbnewv;
				break;
			    }
			}
		    }
		    ++nbnewv;
		}
	    }	    
	}

	// Create vertices
	M.vertices.create_vertices(nbnewv);
	for(index_t f=facets_begin; f<facets_end; ++f) {
	    for(index_t c1=M.facets.corners_begin(f);
		c1<M.facets.corners_end(f); ++c1
	    ) {
		index_t c2 = M.facets.next_corner_around_facet(f,c1);
		index_t v1 = M.facet_corners.vertex(c1);
		index_t v2 = M.facet_corners.vertex(c2);
		index_t v12 = ctov[c1] + nv0;
		const double* p1 = M.vertices.point_ptr(v1);
		const double* p2 = M.vertices.point_ptr(v2);
		double* p12 = M.vertices.point_ptr(v12);
		for(index_t coord=0; coord<M.vertices.dimension(); ++coord) {
		    p12[coord] = 0.5*(p1[coord]+p2[coord]);
		}
	    }
	}

	// Create facets
	M.facets.create_triangles(3*(facets_end - facets_begin));
	for(index_t f=facets_begin; f<facets_end; ++f) {
	    index_t v1 = M.facets.vertex(f,0);
	    index_t v2 = M.facets.vertex(f,1);	
	    index_t v3 = M.facets.vertex(f,2);
	    index_t v12 = ctov[M.facets.corners_begin(f)    ] + nv0;
	    index_t v23 = ctov[M.facets.corners_begin(f) + 1] + nv0;
	    index_t v31 = ctov[M.facets.corners_begin(f) + 2] + nv0;
	    M.facets.set_vertex(f,0,v31);
	    M.facets.set_vertex(f,1,v12);
	    M.facets.set_vertex(f,2,v23);
	    M.facets.set_vertex(nf0+3*(f-facets_begin),0,v1);
	    M.facets.set_vertex(nf0+3*(f-facets_begin),1,v12);
	    M.facets.set_vertex(nf0+3*(f-facets_begin),2,v31);
	    M.facets.set_vertex(nf0+3*(f-facets_begin)+1,0,v12);
	    M.facets.set_vertex(nf0+3*(f-facets_begin)+1,1,v2);
	    M.facets.set_vertex(nf0+3*(f-facets_begin)+1,2,v23); 
	    M.facets.set_vertex(nf0+3*(f-facets_begin)+2,0,v31);
	    M.facets.set_vertex(nf0+3*(f-facets_begin)+2,1,v23);
	    M.facets.set_vertex(nf0+3*(f-facets_begin)+2,2,v3);	    	    
	}
	M.facets.connect();
    }
    
}
