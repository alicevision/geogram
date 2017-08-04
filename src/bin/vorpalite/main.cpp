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
#include <geogram/basic/progress.h>
#include <geogram/basic/stopwatch.h>
#include <geogram/basic/process.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/geometry_nd.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_preprocessing.h>
#include <geogram/mesh/mesh_intersection.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/mesh/mesh_frame_field.h>
#include <geogram/mesh/mesh_tetrahedralize.h>
#include <geogram/mesh/mesh_remesh.h>
#include <geogram/delaunay/LFS.h>
#include <geogram/points/co3ne.h>
#include <typeinfo>

namespace {

    using namespace GEO;

    /**
     * \brief Removes the small connected components from a mesh.
     * \param[in,out] M the mesh
     * \param[in] min_comp_area connected components smaller than
     *   this threshold are discarded
     */
    void remove_small_components(Mesh& M, double min_comp_area) {
        if(min_comp_area == 0.0) {
            return;
        }
        double nb_f_removed = M.facets.nb();
        remove_small_connected_components(M, min_comp_area);
        nb_f_removed -= M.facets.nb();
        if(nb_f_removed != 0) {
            double radius = bbox_diagonal(M);
            double epsilon = CmdLine::get_arg_percent(
                "pre:epsilon", radius
            );
            mesh_repair(M, MESH_REPAIR_DEFAULT, epsilon);
        }
    }

    /**
     * \brief Reconstructs a mesh from a point set.
     * \param[in,out] M_in the input point set and the reconstructed mesh
     */
    void reconstruct(Mesh& M_in) {
        Logger::div("reconstruction");

        Logger::out("Co3Ne") << "Preparing data" << std::endl;
        // Remove all facets
        M_in.facets.clear();
        double bbox_diag = bbox_diagonal(M_in);
        double epsilon = CmdLine::get_arg_percent(
            "pre:epsilon", bbox_diag
        );
        mesh_repair(M_in, MESH_REPAIR_COLOCATE, epsilon);
        index_t nb_neigh = CmdLine::get_arg_uint("co3ne:nb_neighbors");
        index_t Psmooth_iter = CmdLine::get_arg_uint("co3ne:Psmooth_iter");
        double radius = CmdLine::get_arg_percent(
            "co3ne:radius", bbox_diag
        );
        Co3Ne_smooth_and_reconstruct(M_in, nb_neigh, Psmooth_iter, radius);
    }
    
    /**
     * \brief Applies pre-processing to a mesh.
     * \details Pre-processing operations and their parameters are
     *  obtained from the command line.
     * \param[in,out] M_in the mesh to pre-process
     * \return true if resulting mesh is valid, false otherwise
     */
    bool preprocess(Mesh& M_in) {

        Logger::div("preprocessing");
        Stopwatch W("Pre");
        bool pre = CmdLine::get_arg_bool("pre");

        double radius = bbox_diagonal(M_in);
        double area = Geom::mesh_area(M_in, 3);

        int nb_kills = CmdLine::get_arg_int("pre:brutal_kill_borders");
        for(int k=0; k<nb_kills; ++k) {
            vector<index_t> to_kill(M_in.facets.nb(), 0);
            for(index_t f=0; f<M_in.facets.nb(); ++f) {
                for(index_t c=M_in.facets.corners_begin(f); c<M_in.facets.corners_end(f); ++c) {
                    if(M_in.facet_corners.adjacent_facet(c) == NO_FACET) {
                        to_kill[f] = 1;
                    }
                }
            }
            index_t nb_facet_kill=0;
            for(index_t i=0; i<to_kill.size(); ++i) {
                if(to_kill[i] != 0) {
                    ++nb_facet_kill;
                }
            }
            if(nb_facet_kill != 0) {
                Logger::out("Pre")
                    << "Killed " << nb_facet_kill << " facet(s) on border"
                    << std::endl;
                M_in.facets.delete_elements(to_kill);
                mesh_repair(M_in);
            } else {
                Logger::out("Pre") << "No facet on border (good)" << std::endl;
                break;
            }
        }
        
        index_t nb_bins = CmdLine::get_arg_uint("pre:vcluster_bins");
        if(pre && nb_bins != 0) {
            Logger::err("Remesh")
                << "Mesh clustering only available in Vorpaline"
                << std::endl;
            exit(-1);
        } else if(pre && CmdLine::get_arg_bool("pre:repair")) {
            MeshRepairMode mode = MESH_REPAIR_DEFAULT;
            double epsilon = CmdLine::get_arg_percent(
                "pre:epsilon", radius
            );
            mesh_repair(M_in, mode, epsilon);
        }

        if(pre) {
            remove_small_components(
                M_in, CmdLine::get_arg_percent(
                    "pre:min_comp_area", area
                )
            );
        }

        if(pre) {
            double max_area = CmdLine::get_arg_percent(
                "pre:max_hole_area", area
            );
            index_t max_edges = CmdLine::get_arg_uint(
                "pre:max_hole_edges"
            );
            if(max_area != 0.0 && max_edges != 0) {
                fill_holes(M_in, max_area, max_edges);
            }
        }

        double anisotropy = 0.02 * CmdLine::get_arg_double("remesh:anisotropy");
        if(anisotropy != 0.0) {
            compute_normals(M_in);
            unsigned int nb_normal_smooth =
                CmdLine::get_arg_uint("pre:Nsmooth_iter");
            if(nb_normal_smooth != 0) {
                Logger::out("Nsmooth") << "Smoothing normals, "
                    << nb_normal_smooth
                    << " iteration(s)" << std::endl;
                simple_Laplacian_smooth(M_in, nb_normal_smooth, true);
            }
            set_anisotropy(M_in, anisotropy);
        }

        if(CmdLine::get_arg_bool("remesh")) {
            unsigned int nb_removed = M_in.facets.nb();
            remove_small_facets(M_in, 1e-30);
            nb_removed -= M_in.facets.nb();
            if(nb_removed == 0) {
                Logger::out("Validate")
                    << "Mesh does not have 0-area facets (good)" << std::endl;
            } else {
                Logger::out("Validate")
                    << "Removed " << nb_removed
                    << " 0-area facets" << std::endl;
            }
        }

        double margin = CmdLine::get_arg_percent(
            "pre:margin", radius
        );
        if(pre && margin != 0.0) {
            expand_border(M_in, margin);
        }

        if(M_in.facets.nb() == 0) {
            Logger::warn("Preprocessing")
                << "After pre-processing, got an empty mesh"
                << std::endl;
            // return false;
        }

        return true;
    }

    /**
     * \brief Applies post-processing to a mesh
     * \details Post-processing operations and their parameters are
     *  obtained from the command line.
     * \param[in,out] M_out the mesh to pre-process
     * \return true if resulting mesh is valid, false otherwise
     */
    bool postprocess(Mesh& M_out) {
        Logger::div("postprocessing");
        {
            Stopwatch W("Post");
            if(CmdLine::get_arg_bool("post")) {
                double radius = bbox_diagonal(M_out);
                double area = Geom::mesh_area(M_out, 3);
                if(CmdLine::get_arg_bool("post:repair")) {
                    double epsilon = CmdLine::get_arg_percent(
                        "pre:epsilon", radius
                    );
                    mesh_repair(M_out, MESH_REPAIR_DEFAULT, epsilon);
                }
                remove_small_components(
                    M_out, CmdLine::get_arg_percent(
                        "post:min_comp_area", area
                    )
                );
                double max_area = CmdLine::get_arg_percent(
                    "post:max_hole_area", area
                );
                index_t max_edges = CmdLine::get_arg_uint(
                    "post:max_hole_edges"
                );
                if(max_area != 0.0 && max_edges != 0) {
                    fill_holes(M_out, max_area, max_edges);
                }
                double deg3_dist = CmdLine::get_arg_percent(
                    "post:max_deg3_dist", radius
                );
                while(remove_degree3_vertices(M_out, deg3_dist) != 0) {}
                if(CmdLine::get_arg_bool("post:isect")) {
                    mesh_remove_intersections(M_out);
                }
            }
            orient_normals(M_out);
            if(CmdLine::get_arg_bool("post:compute_normals")) {
                Attribute<double> normal;
                normal.bind_if_is_defined(
                    M_out.vertices.attributes(),
                    "normal"
                );
                if(!normal.is_bound()) {
                    normal.create_vector_attribute(
                        M_out.vertices.attributes(),
                        "normal",
                        3
                    );
                }
                for(index_t f=0; f<M_out.facets.nb(); ++f) {
                    vec3 N = Geom::mesh_facet_normal(M_out,f);
                    N = normalize(N);
                    for(index_t lv=0; lv<M_out.facets.nb_vertices(f); ++lv) {
                        index_t v = M_out.facets.vertex(f,lv);
                        normal[3*v  ] = N.x;
                        normal[3*v+1] = N.y;
                        normal[3*v+2] = N.z;                        
                    }
                }
            }
        }

        Logger::div("result");
        M_out.show_stats("FinalMesh");
        if(M_out.facets.nb() == 0) {
            Logger::warn("Postprocessing")
                << "After post-processing, got an empty mesh"
                << std::endl;
            // return false;
        }

        return true;
    }

    /**
     * \brief Generates a tetrahedral mesh.
     * \param[in] input_filename name of the input file, can be
     *   either a closed surfacic mesh or a tetrahedral mesh
     * \param[in] output_filename name of the output file
     * \retval 0 on success
     * \retval non-zero value otherwise
     */
    int tetrahedral_mesher(
        const std::string& input_filename, const std::string& output_filename
    ) {
        MeshIOFlags flags;
        flags.set_element(MESH_CELLS);

        Mesh M_in;
        if(!mesh_load(input_filename, M_in, flags)) {
            return 1;
        }
        mesh_tetrahedralize(
            M_in,
            CmdLine::get_arg_bool("tet:preprocess"),
            CmdLine::get_arg_bool("tet:refine"),
            CmdLine::get_arg_double("tet:quality")
        );
        M_in.cells.compute_borders();
        if(!mesh_save(M_in, output_filename, flags)) {
            return 1;
        }
        return 0;
    }
    
}

int main(int argc, char** argv) {
    using namespace GEO;

    GEO::initialize();

    try {

        Stopwatch total("Total time");

        CmdLine::import_arg_group("standard");
        CmdLine::import_arg_group("pre");
        CmdLine::import_arg_group("remesh");
        CmdLine::import_arg_group("algo");
        CmdLine::import_arg_group("post");
        CmdLine::import_arg_group("opt");
        CmdLine::import_arg_group("co3ne");        
        CmdLine::import_arg_group("tet");
        
        std::vector<std::string> filenames;

#ifdef GEO_OS_ANDROID
        Logger::out("Math") << "Checking math functions: sin 45 = "
                            << sin(M_PI/4.0) << std::endl;
        Logger::out("Math") << "Checking double conversion: 1.2345 = "
                            << atof("1.2345")
                            << std::endl;
        Logger::out("Math") << "Checking math functions: atan(1.0) = "
                            << atan(1.0) << std::endl;
        
#endif        

        
        if(!CmdLine::parse(argc, argv, filenames, "inputfile <outputfile>")) {
            return 1;
        }
        std::string input_filename = filenames[0];
        std::string output_filename =
            filenames.size() >= 2 ? filenames[1] : std::string("out.meshb");
        Logger::out("I/O") << "Output = " << output_filename << std::endl;
        CmdLine::set_arg("input", input_filename);
        CmdLine::set_arg("output", output_filename);

        if(CmdLine::get_arg_bool("tet")) {
            return tetrahedral_mesher(input_filename, output_filename);
        }
        
        Mesh M_in, M_out;
        {
            Stopwatch W("Load");
            if(!mesh_load(input_filename, M_in)) {
                return 1;
            }
        }

        if(CmdLine::get_arg_bool("co3ne")) {
            reconstruct(M_in);
        }
        
        if(!preprocess(M_in)) {
            return 1;
        }

        if(!CmdLine::get_arg_bool("remesh")) {
            if(!postprocess(M_in)) {
                return 1;
            }
            if(!mesh_save(M_in, output_filename)) {
                return 1;
            }
            return 0;
        }

        Logger::div("remeshing");
        if(CmdLine::get_arg_bool("remesh:by_parts")) {
            Logger::err("Remesh")
                << "By parts remeshing only available in Vorpaline"
                << std::endl;
            exit(-1);
        } else if(CmdLine::get_arg_bool("remesh:sharp_edges")) {
            Logger::err("Remesh")
                << "Feature-sensitive remeshing only available in Vorpaline"
                << std::endl;
            exit(-1);
        } else {
            double gradation = CmdLine::get_arg_double("remesh:gradation");
            if(gradation != 0.0) {
                compute_sizing_field(
                    M_in, gradation, CmdLine::get_arg_uint("remesh:lfs_samples")
                );
            }
            coord_index_t dim = 3;
            double anisotropy =
                0.02 * CmdLine::get_arg_double("remesh:anisotropy");
            if(anisotropy != 0.0 && M_in.vertices.dimension() >= 6) {
                dim = 6;
            }
            remesh_smooth(
                M_in, M_out,
                CmdLine::get_arg_uint("remesh:nb_pts"),
                dim,
                CmdLine::get_arg_uint("opt:nb_Lloyd_iter"),
                CmdLine::get_arg_uint("opt:nb_Newton_iter"),
                CmdLine::get_arg_uint("opt:Newton_m")
            );
        }

        if(M_out.facets.nb() == 0) {
            Logger::err("Remesh") << "After remesh, got an empty mesh"
                << std::endl;
            return 1;
        }

        if(!postprocess(M_out)) {
            return 1;
        }

        if(!mesh_save(M_out, output_filename)) {
            return 1;
        }

    }
    catch(const std::exception& e) {
        std::cerr << "Received an exception: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

