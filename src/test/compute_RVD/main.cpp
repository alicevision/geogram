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
#include <geogram/basic/process.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_topology.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/mesh_fill_holes.h>
#include <geogram/mesh/mesh_preprocessing.h>
#include <geogram/mesh/mesh_degree3_vertices.h>
#include <geogram/delaunay/delaunay.h>
#include <geogram/voronoi/RVD.h>
#include <geogram/numerics/predicates.h>

namespace {

    using namespace GEO;

    /**
     * \brief Removes zero area facets in a mesh
     * \param[in] M the input mesh
     */
    void check_for_zero_area_facets(Mesh& M) {
        vector<index_t> remove_f;
        vec3 q1(0, 0, 0);
        vec3 q2(0, 0, 1);
        vec3 q3(0, 1, 0);
        vec3 q4(1, 0, 0);
        for(index_t f = 0; f < M.facets.nb(); ++f) {
            index_t c = M.facets.corners_begin(f);
            index_t v1 = M.facet_corners.vertex(c);
            index_t v2 = M.facet_corners.vertex(c + 1);
            index_t v3 = M.facet_corners.vertex(c + 2);
            const vec3& p1 = Geom::mesh_vertex(M, v1);
            const vec3& p2 = Geom::mesh_vertex(M, v2);
            const vec3& p3 = Geom::mesh_vertex(M, v3);

            // Colinearity is tested by using four coplanarity
            // tests with points q1,q2,q3,q4 that are
            // not coplanar.
            if(
                PCK::orient_3d(p1, p2, p3, q1) == 0.0 &&
                PCK::orient_3d(p1, p2, p3, q2) == 0.0 &&
                PCK::orient_3d(p1, p2, p3, q3) == 0.0 &&
                PCK::orient_3d(p1, p2, p3, q4) == 0.0
            ) {
                Logger::warn("Validate") << "Found a zero-area facet"
                    << std::endl;
                remove_f.resize(M.facets.nb(), 0);
                remove_f[f] = 1;
            }
        }
        if(remove_f.size() != 0) {
            Logger::warn("Validate") << "Removing zero-area facet(s)"
                << std::endl;
            M.facets.delete_elements(remove_f);
        }
    }
}

int main(int argc, char** argv) {
    using namespace GEO;

    GEO::initialize();

    try {

        Stopwatch W("Total time");

        std::vector<std::string> filenames;

        CmdLine::import_arg_group("standard");
        CmdLine::import_arg_group("algo");
        CmdLine::declare_arg("volumetric", false, "volumetric or surfacic RVD");
        CmdLine::declare_arg(
            "cell_borders", false, "generate only cell borders"
        );
        CmdLine::declare_arg(
            "integration_smplx", false,
            "in volumetric mode, generate integration simplices"
        );
        CmdLine::declare_arg("RDT", false, "save RDT");
        CmdLine::declare_arg("RVD", true, "save RVD");       
        CmdLine::declare_arg_percent(
            "epsilon",0.001,
            "Tolerance for merging vertices relative to bbox diagonal"
        );
        CmdLine::declare_arg("constrained", false, "constrained Delaunay");       
        CmdLine::declare_arg("prefer_seeds", false, "in constrained mode, use seeds whenever possible");

        if(
            !CmdLine::parse(
                argc, argv, filenames, "meshfile <pointsfile> <outputfile>"
            )
        ) {
            return 1;
        }

        std::string mesh_filename = filenames[0];
        std::string points_filename = filenames[0];
        if(filenames.size() >= 2) {
            points_filename = filenames[1];
        }

        bool volumetric = CmdLine::get_arg_bool("volumetric");
        bool cell_borders = CmdLine::get_arg_bool("cell_borders");
        bool integ_smplx = CmdLine::get_arg_bool("integration_smplx");
        std::string output_filename;
        if(filenames.size() >= 3) {
            output_filename = filenames[2];
        } else {
            if(volumetric || CmdLine::get_arg_bool("constrained")) {
                output_filename = "out.meshb";
            } else {
                output_filename = "out.eobj";
            }
        }

        Logger::out("I/O") << "Output = " << output_filename << std::endl;

        Logger::div("Loading data");

        Mesh M_in, points_in;
        Mesh M_out;

        {
            MeshIOFlags flags;
            if(volumetric) {
                flags.set_element(MESH_CELLS);
            }
            if(!mesh_load(mesh_filename, M_in, flags)) {
                return 1;
            }
        }

        if(!volumetric) {
            mesh_repair(M_in);
            check_for_zero_area_facets(M_in);
        }

        if(!mesh_load(points_filename, points_in)) {
            return 1;
        }
        
        double epsilon = CmdLine::get_arg_percent(
            "epsilon",bbox_diagonal(points_in)
        );
        points_in.facets.clear();
        mesh_repair(points_in, MESH_REPAIR_COLOCATE, epsilon);

        geo_assert(points_in.vertices.dimension() == 3);

        if(CmdLine::get_arg_bool("constrained")) {

            Mesh surface;
            vector<double> inner_points;


            {
                Logger::div("Computing the surface");
                Delaunay_var delaunay = Delaunay::create(3);
                RestrictedVoronoiDiagram_var RVD = 
                    RestrictedVoronoiDiagram::create(delaunay,&M_in);
                delaunay->set_vertices(
                    points_in.vertices.nb(), points_in.vertices.point_ptr(0)
                );


                RestrictedVoronoiDiagram::RDTMode mode =                 
                    RestrictedVoronoiDiagram::RDTMode(
                        RestrictedVoronoiDiagram::RDT_MULTINERVE |
                        RestrictedVoronoiDiagram::RDT_RVC_CENTROIDS 
                    );

                if(CmdLine::get_arg_bool("prefer_seeds")) {
                    mode = RestrictedVoronoiDiagram::RDTMode(
                        mode | RestrictedVoronoiDiagram::RDT_PREFER_SEEDS
                    );
                }


                RVD->compute_RDT(surface, mode);

                mesh_repair(surface);
                remove_small_connected_components(surface,0.0,100);
                fill_holes(surface, 1e30);
                double radius = bbox_diagonal(surface);
                remove_degree3_vertices(surface, 0.01*radius);
                mesh_save(surface,"surface.meshb");

                vector<double> m(points_in.vertices.nb());
                vector<double> mg(points_in.vertices.nb()*3);
                RVD->compute_centroids_on_surface(&mg[0], &m[0]);
                for(index_t v=0; v<points_in.vertices.nb(); ++v) {
                    if(m[v] == 0.0) {
                        inner_points.push_back(
                            points_in.vertices.point_ptr(v)[0]
                        );
                        inner_points.push_back(
                            points_in.vertices.point_ptr(v)[1]
                        );
                        inner_points.push_back(
                            points_in.vertices.point_ptr(v)[2]
                        );
                    }
                }
            }

            Logger::div("Calling tetgen");
            Delaunay_var delaunay = Delaunay::create(3,"tetgen");
            delaunay->set_constraints(&surface);
            delaunay->set_vertices(inner_points.size()/3, &inner_points[0]);

            vector<double> pts(delaunay->nb_vertices() * 3);
            vector<index_t> tet2v(delaunay->nb_cells() * 4);
            for(index_t v = 0; v < delaunay->nb_vertices(); ++v) {
                pts[3 * v] = delaunay->vertex_ptr(v)[0];
                pts[3 * v + 1] = delaunay->vertex_ptr(v)[1];
                pts[3 * v + 2] = delaunay->vertex_ptr(v)[2];
            }
            for(index_t t = 0; t < delaunay->nb_cells(); ++t) {
                tet2v[4 * t] = index_t(delaunay->cell_vertex(t, 0));
                tet2v[4 * t + 1] = index_t(delaunay->cell_vertex(t, 1));
                tet2v[4 * t + 2] = index_t(delaunay->cell_vertex(t, 2));
                tet2v[4 * t + 3] = index_t(delaunay->cell_vertex(t, 3));
            }
            M_out.cells.assign_tet_mesh(3, pts, tet2v, true);
            M_out.show_stats();

            Logger::div("Saving the result");
            MeshIOFlags flags;
            flags.set_element(MESH_CELLS);
            mesh_save(M_out, output_filename, flags);
        } else {
            Delaunay_var delaunay = Delaunay::create(3);
            RestrictedVoronoiDiagram_var RVD = RestrictedVoronoiDiagram::create(
                delaunay, &M_in
            );
            delaunay->set_vertices(
                points_in.vertices.nb(), points_in.vertices.point_ptr(0)
            );

            RVD->set_volumetric(volumetric);
       
       
            if(CmdLine::get_arg_bool("RVD")) {
                Logger::div("Restricted Voronoi Diagram");
                
                RVD->compute_RVD(M_out, 0, cell_borders, integ_smplx);
                if(integ_smplx && volumetric) {
                    M_out.cells.connect();
                    M_out.cells.compute_borders();
                }
                Logger::div("Result");

                MeshIOFlags flags;
                flags.set_attribute(MESH_FACET_REGION);
                flags.set_attribute(MESH_CELL_REGION);
                flags.set_element(MESH_CELLS);
                mesh_save(M_out, output_filename, flags);
            }
            
            if(CmdLine::get_arg_bool("RDT")) {
                Logger::out("RDT") << "Computing RDT..." << std::endl;
                Mesh RDT ;
                RVD->compute_RDT(RDT);
                MeshIOFlags flags;
                if(volumetric) {
                    flags.set_elements(MESH_CELLS);
                }
                mesh_save(RDT, "RDT.meshb", flags);
            }


            if(!meshes_have_same_topology(M_in, M_out, true)) {
                Logger::out("") << "Returning error code (2)" << std::endl;
                return 2;
            }
        }
    }
    catch(const std::exception& e) {
        std::cerr << "Received an exception: " << e.what() << std::endl;
        return 1;
    }

    Logger::out("") << "Everything OK, Returning status 0" << std::endl;
    return 0;
}

