#ifndef H_HEX_H
#define H_HEX_H

#include <geogram/mesh/mesh.h>
namespace GEO {

	void  hex_set_2_hex_mesh(Mesh* hex, Mesh* quadtri = NULL);
	void  kill_intersecting_hexes(Mesh* hex);
}
#endif
