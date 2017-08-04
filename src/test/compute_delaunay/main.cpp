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

#include <geogram/basic/common.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/stopwatch.h>
#include <geogram/basic/file_system.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/delaunay/delaunay.h>
#include <algorithm>

namespace {
    using namespace GEO;
    // local user functions to be put here
}

int main(int argc, char** argv) {
    using namespace GEO;

    GEO::initialize();

    try {

        Stopwatch Wtot("Total time");

        std::vector<std::string> filenames;

        CmdLine::import_arg_group("standard");
        CmdLine::import_arg_group("algo");

        CmdLine::declare_arg(
            "convex_hull", false,
            "compute just the convex hull of the points"
        );

        CmdLine::declare_arg(
            "dimension", 3, "3 for 3D, 2 for 2D"
        );

        CmdLine::set_arg("algo:delaunay","default");
        
        if(
            !CmdLine::parse(
                argc, argv, filenames, "pointsfile <outputfile|none>"
            )
        ) {
            return 1;
        }


        std::string points_filename = filenames[0];

        std::string output_filename =
            filenames.size() >= 2 ? filenames[1] : std::string("out.meshb");

        bool output = (output_filename != "none");

        Logger::div("Data I/O");

        Logger::out("I/O") << "Output = " << output_filename << std::endl;

        Mesh M_in;
        Mesh M_out;

        {
            MeshIOFlags flags;
            flags.reset_element(MESH_FACETS);
            flags.reset_element(MESH_CELLS);
            if(!mesh_load(points_filename, M_in, flags)) {
                return 1;
            }
        }

        bool convex_hull_only = CmdLine::get_arg_bool("convex_hull");
        index_t dimension = index_t(CmdLine::get_arg_int("dimension"));

        std::string del = CmdLine::get_arg("algo:delaunay");
        if(del == "default") {
            if(dimension == 3) {
                if(DelaunayFactory::has_creator("PDEL")) {
                    CmdLine::set_arg("algo:delaunay", "PDEL");
                } else {
                    CmdLine::set_arg("algo:delaunay", "BDEL");
                }
            } else if(dimension == 2) {
                CmdLine::set_arg("algo:delaunay", "triangle");                
            }
        }

        if(convex_hull_only && dimension != 3) {
            Logger::err("Delaunay")
                << "Convex hull only implemented in dimension 3"
                << std::endl;
            return 1;
        }
        
        Logger::div("Computing 3D Delaunay triangulation");

        Logger::out("Delaunay")
            << "Using " << CmdLine::get_arg("algo:delaunay") << std::endl;

        Delaunay_var delaunay = Delaunay::create(coord_index_t(dimension));

        M_in.vertices.set_dimension(dimension);
        
        {
            Stopwatch Wdel("Delaunay");

            //   If we want the convex hull, we keep the infinite facets,
            // because the convex hull can be retreived as the finite facets
            // of the infinite cells (note: it would be also possible to
            // throw away the infinite cells and get the convex hull as
            // the facets adjacent to no cell).
            if(convex_hull_only) {
                delaunay->set_keeps_infinite(true);
            }
            
            delaunay->set_vertices(
                M_in.vertices.nb(), M_in.vertices.point_ptr(0)
            );
        }

        Logger::out("Delaunay") << delaunay->nb_cells() << " tetrahedra"
            << std::endl;

        if(output) {
            vector<double> pts(delaunay->nb_vertices() * 3);
            for(index_t v = 0; v < delaunay->nb_vertices(); ++v) {
                pts[3 * v] = delaunay->vertex_ptr(v)[0];
                pts[3 * v + 1] = delaunay->vertex_ptr(v)[1];
                pts[3 * v + 2] = (dimension == 3) ? delaunay->vertex_ptr(v)[2] : 0.0;
            }

            if(convex_hull_only) {
                
                // The convex hull can be retrieved as the finite facets
                // of the infinite cells (note: it would be also possible to
                // throw away the infinite cells and get the convex hull as
                // the facets adjacent to no cell). Here we use the infinite
                // cells to show an example with them.


                // This block is just a sanity check
                {
                    for(index_t t=0; t < delaunay->nb_finite_cells(); ++t) {
                        geo_debug_assert(delaunay->cell_is_finite(t));
                    }
                
                    for(index_t t=delaunay->nb_finite_cells();
                        t < delaunay->nb_cells(); ++t) {
                        geo_debug_assert(delaunay->cell_is_infinite(t));
                    }
                }
                
                vector<index_t> tri2v;

                // This iterates on the infinite cells
                for(
                    index_t t = delaunay->nb_finite_cells();
                    t < delaunay->nb_cells(); ++t
                ) {
                    for(index_t lv=0; lv<4; ++lv) {
                        signed_index_t v = delaunay->cell_vertex(t,lv);
                        if(v != -1) {
                            tri2v.push_back(index_t(v));
                        }
                    }
                }
                M_out.facets.assign_triangle_mesh(3, pts, tri2v, true);
                
            } else if(dimension == 3) {
                vector<index_t> tet2v(delaunay->nb_cells() * 4);
                for(index_t t = 0; t < delaunay->nb_cells(); ++t) {
                    tet2v[4 * t] = index_t(delaunay->cell_vertex(t, 0));
                    tet2v[4 * t + 1] = index_t(delaunay->cell_vertex(t, 1));
                    tet2v[4 * t + 2] = index_t(delaunay->cell_vertex(t, 2));
                    tet2v[4 * t + 3] = index_t(delaunay->cell_vertex(t, 3));
                }
                M_out.cells.assign_tet_mesh(3, pts, tet2v, true);
            } else if(dimension == 2) {
                vector<index_t> tri2v(delaunay->nb_cells() * 3);
                for(index_t t = 0; t < delaunay->nb_cells(); ++t) {
                    tri2v[3 * t] = index_t(delaunay->cell_vertex(t, 0));
                    tri2v[3 * t + 1] = index_t(delaunay->cell_vertex(t, 1));
                    tri2v[3 * t + 2] = index_t(delaunay->cell_vertex(t, 2));
                }
                M_out.facets.assign_triangle_mesh(3, pts, tri2v, true);
            }
            M_out.show_stats();
                
            Logger::div("Saving the result");
            MeshIOFlags flags;
            flags.set_element(MESH_FACETS);            
            flags.set_element(MESH_CELLS);
            mesh_save(M_out, output_filename, flags);
        }
    }
    catch(const std::exception& e) {
        std::cerr << "Received an exception: " << e.what() << std::endl;
        return 1;
    }

    Logger::out("") << "Everything OK, Returning status 0" << std::endl;
    return 0;
}

