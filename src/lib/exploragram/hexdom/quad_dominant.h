#ifndef H_QUAD_DOMINANT_H
#define H_QUAD_DOMINANT_H

#include <exploragram/basic/common.h>
#include <geogram/mesh/mesh.h>
namespace GEO {


	/**
	* All informations required to produce a hex-dom mesh are exported
	* into a 2D non-manifold mesh with uv coordinates
	*
	* This new mesh is the only input of the following steps in the
	* hex dom generation pipeline
	*/

	void export_boundary_with_uv(Mesh* m, Mesh* hex, const char* uv_name, const char* singtri_name);

	void split_edges_by_iso_uvs(Mesh* hex, const char * uv_name, const char *singular_name);

	void facets_split(Mesh* m, const char *uv_name, const char *singular_name);

	void mark_charts(Mesh* m, const char *uv_name, const char *charts_name, const char *singular_name);

	void simplify_quad_charts(Mesh* m);

	void export_quadtri_from_charts(Mesh* hex);

}
#endif
