/* Vorpaline - geogram demo program */



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
#include <geogram/mesh/mesh_decimate.h>
#include <geogram/mesh/mesh_remesh.h>

#include <geogram/delaunay/LFS.h>
#include <geogram/voronoi/CVT.h>
#include <geogram/voronoi/RVD_mesh_builder.h>
#include <geogram/points/co3ne.h>

#include <geogram/third_party/PoissonRecon/poisson_geogram.h>

#include <typeinfo>
#include <algorithm>

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

        index_t Psmooth_iter = CmdLine::get_arg_uint("co3ne:Psmooth_iter");
        index_t nb_neigh = CmdLine::get_arg_uint("co3ne:nb_neighbors");
        
        if(CmdLine::get_arg("algo:reconstruct") == "Poisson") {
            if(Psmooth_iter != 0) {
                Co3Ne_smooth(M_in, nb_neigh, Psmooth_iter);
            }
            bool has_normals = false;
            {
                Attribute<double> normal;
                normal.bind_if_is_defined(M_in.vertices.attributes(), "normal");
                has_normals = (normal.is_bound() && normal.dimension() == 3);
            }
            
            if(!has_normals) {
                // TODO: add a way of making normals orientation coherent
                // in Co3Ne_compute_normals...

                if(M_in.facets.nb() != 0) {
                    Attribute<double> normal;
                    normal.bind_if_is_defined(
                        M_in.vertices.attributes(),"normal"
                    );
                    if(!normal.is_bound()) {
                        normal.create_vector_attribute(
                            M_in.vertices.attributes(), "normal", 3
                        );
                    }
                    for(index_t i=0; i<M_in.vertices.nb()*3; ++i) {
                        normal[i]=0.0;
                    }
                    for(index_t f=0; f<M_in.facets.nb(); ++f) {
                        vec3 N = Geom::mesh_facet_normal(M_in, f);
                        for(index_t lv=0; lv<M_in.facets.nb_vertices(f); ++lv) {
                            index_t v = M_in.facets.vertex(f,lv);
                            normal[3*v  ] += N.x;
                            normal[3*v+1] += N.y;
                            normal[3*v+2] += N.z;
                        }
                    }
                    for(index_t v=0; v<M_in.vertices.nb(); ++v) {
                        vec3 N(normal[3*v],normal[3*v+1],normal[3*v+2]);
                        N = normalize(N);
                        normal[3*v  ]=N.x;
                        normal[3*v+1]=N.y;
                        normal[3*v+2]=N.z;
                    }
                } else {
                    Logger::out("Poisson")
                        << "Dataset has no normals, estimating them"
                        << std::endl;
                    Logger::out("Poisson")
                    << "(result may be not so good, normals may be incoherent)"
                    << std::endl;
                    Co3Ne_compute_normals(M_in, nb_neigh, true);
                }
            }
            
            index_t depth = CmdLine::get_arg_uint("poisson:octree_depth");
            Mesh M_out;
            PoissonReconstruction recons;
            recons.set_depth(depth);
            Logger::out("Reconstruct")
                << "Starting Poisson reconstruction..."
                << std::endl;
            recons.reconstruct(&M_in, &M_out);
            Logger::out("Reconstruct")
                << "Poisson reconstruction done."
                << std::endl;
            MeshElementsFlags what = MeshElementsFlags(
                MESH_VERTICES | MESH_FACETS 
            );
            M_in.copy(M_out, true, what);
        } else {
            // Remove all facets
            M_in.facets.clear();
            
            double bbox_diag = bbox_diagonal(M_in);
            double epsilon = CmdLine::get_arg_percent(
                "pre:epsilon", bbox_diag
            );
            mesh_repair(M_in, MESH_REPAIR_COLOCATE, epsilon);
            
            double radius = CmdLine::get_arg_percent(
                "co3ne:radius", bbox_diag
            );
            Co3Ne_smooth_and_reconstruct(M_in, nb_neigh, Psmooth_iter, radius);
        }
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
                for(index_t c=M_in.facets.corners_begin(f);
                    c<M_in.facets.corners_end(f); ++c
                ) {
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
            mesh_decimate_vertex_clustering(M_in, nb_bins);
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
            index_t nb_normal_smooth =
                CmdLine::get_arg_uint("pre:Nsmooth_iter");
            if(nb_normal_smooth != 0) {
                Logger::out("Nsmooth") << "Smoothing normals, "
                    << nb_normal_smooth
                    << " iteration(s)" << std::endl;
                simple_Laplacian_smooth(M_in, index_t(nb_normal_smooth), true);
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
     * \brief The callback called for each RVD polyhedron. Constructs a 
     *  mesh with the boundary of all cells.
     * \details Its member functions are called for each RVD polyhedron, 
     *  i.e. the intersections between the volumetric mesh tetrahedra and
     *  the Voronoi cells. Based on set_simplify_xxx(), a smaller number of
     *  polyhedra can be generated.
     */
    class SaveRVDCells : public RVDPolyhedronCallback {
    public:

	/**
	 * \brief SaveRVDCells constructor.
	 * \param[out] output_mesh a reference to the generated mesh 
	 */
	SaveRVDCells(Mesh& output_mesh) :
	    output_mesh_(output_mesh), shrink_(0.0) {
	    cell_vertex_map_ = nil;
	    global_vertex_map_ = nil;
	    current_cell_id_ = 0;
	    generate_ids_ = false;
	}

	/**
	 * \brief SaveRVDCells destructor.
	 */
	~SaveRVDCells() {
	    if(generate_ids_) {
		cell_id_.unbind();
		seed_id_.unbind();
		vertex_id_.unbind();
		facet_seed_id_.unbind();
		delete global_vertex_map_;
		global_vertex_map_ = nil;
	    }
	    delete cell_vertex_map_;
	    cell_vertex_map_ = nil;
	}

	/**
	 * \brief Specifies whether ids should be generated.
	 * \details If enabled, unique vertex ids, seed ids and cell ids are
	 *  generated and attached to the mesh vertices ("vertex_id" attribute)
	 *  and mesh facets ("seed_id" and "cell_id" attributes) respectively.
	 *  There is a cell_id per generated polyhedron, and seed_id refers to
	 *  the Voronoi seed (the point that the Voronoi cell is associated 
	 *  with).
	 * \param[in] x true if ids should be generated, false
	 *  otherwise (default)
	 */
	void set_generate_ids(bool x) {
	    if(x == generate_ids_) {
		return;
	    }
	    generate_ids_ = x;
	    if(generate_ids_) {
		cell_id_.bind(
		    output_mesh_.facets.attributes(), "cell_id"
		);
		seed_id_.bind(
		    output_mesh_.facets.attributes(), "seed_id"		    
		);
		vertex_id_.bind(
		    output_mesh_.vertices.attributes(), "vertex_id"
		);
		facet_seed_id_.bind(
		    output_mesh_.facets.attributes(), "facet_seed_id"
		);
		global_vertex_map_ = new RVDVertexMap;
	    } else {
		cell_id_.unbind();
		seed_id_.unbind();
		vertex_id_.unbind();
		facet_seed_id_.unbind();
		delete global_vertex_map_;
		global_vertex_map_ = nil;
	    }
	}

	
	/**
	 * \brief Defines the optional shrink factor for cells.
	 * \param[in] x shrink factor, 0.0 means no shrink, 1.0 means
	 *  maximum shrink (cell reduced to a point).
	 */
	void set_shrink(double x) {
	    shrink_ = x;
	    if(shrink_ != 0.0) {
		set_use_mesh(true);
	    }
	}
	
	/**
	 * \brief Called at the beginning of RVD traversal.
	 */
	virtual void begin() {
	    RVDPolyhedronCallback::begin();
	    output_mesh_.clear();
	    output_mesh_.vertices.set_dimension(3);
	}

	/**
	 * \brief Called at the end of RVD traversal.
	 */
	virtual void end() {
	    RVDPolyhedronCallback::end();
	    output_mesh_.facets.connect();
	}

	/**
	 * \brief Called at the beginning of each RVD polyhedron.
	 * \param[in] seed , tetrahedron the (seed,tetrahedron) pair that
	 *  defines the RVD polyhedron, as the intersection between the Voronoi
	 *  cell of the seed and the tetrahedron.
	 */
	virtual void begin_polyhedron(index_t seed, index_t tetrahedron) {
	    geo_argused(tetrahedron);
	    geo_argused(seed);
	    delete cell_vertex_map_;
	    cell_vertex_map_ = new RVDVertexMap;
	    cell_vertex_map_->set_first_vertex_index(
		output_mesh_.vertices.nb()
	    );
	}

	virtual void begin_facet(index_t facet_seed, index_t facet_tet_facet) {
	    geo_argused(facet_seed);
	    geo_argused(facet_tet_facet);
	    current_facet_.resize(0);
	}

	virtual void vertex(
	    const double* geometry, const GEOGen::SymbolicVertex& symb
	) {
	    index_t v = cell_vertex_map_->find_or_create_vertex(seed(), symb);
	    if(v >= output_mesh_.vertices.nb()) {
		output_mesh_.vertices.create_vertex(geometry);
		if(generate_ids_) {
		    vertex_id_[v] = int(
			global_vertex_map_->find_or_create_vertex(seed(), symb)
		    );
		}
	    }
	    current_facet_.push_back(v);
	}

	virtual void end_facet() {
	    index_t f = output_mesh_.facets.nb();
	    output_mesh_.facets.create_polygon(current_facet_.size());
	    for(index_t i=0; i<current_facet_.size(); ++i) {
		output_mesh_.facets.set_vertex(f,i,current_facet_[i]);
	    }
	    if(generate_ids_) {
		seed_id_[f] = int(seed());
		cell_id_[f] = int(current_cell_id_);
		facet_seed_id_[f] = int(facet_seed());
	    }
	}

	virtual void end_polyhedron() {
	    ++current_cell_id_;
	}

	virtual void process_polyhedron_mesh() {
	    if(shrink_ != 0.0 && mesh_.vertices.nb() != 0) {
		vec3 center(0.0, 0.0, 0.0);
		for(index_t v=0; v<mesh_.vertices.nb(); ++v) {
		    center += vec3(mesh_.vertices.point_ptr(v));
		}
		center = (1.0 / double(mesh_.vertices.nb())) * center;
		for(index_t v=0; v<mesh_.vertices.nb(); ++v) {
		    vec3 p(mesh_.vertices.point_ptr(v));
		    p = shrink_ * center + (1.0 - shrink_) * p;
		    mesh_.vertices.point_ptr(v)[0] = p.x;
		    mesh_.vertices.point_ptr(v)[1] = p.y;
		    mesh_.vertices.point_ptr(v)[2] = p.z;		    
		}
	    }
	    RVDPolyhedronCallback::process_polyhedron_mesh();
	}
	
    private:
	vector<index_t> current_facet_;
	Mesh& output_mesh_;
	RVDVertexMap* global_vertex_map_;
	RVDVertexMap* cell_vertex_map_;
	double shrink_;
	bool generate_ids_;
	Attribute<int> cell_id_;
	Attribute<int> seed_id_;
	Attribute<int> vertex_id_;
	Attribute<int> facet_seed_id_;
	index_t current_cell_id_;
    };


    /**
     * \brief Generates a polyhedral mesh.
     * \param[in] input_filename name of the input file, can be
     *   either a closed surfacic mesh or a tetrahedral mesh
     * \param[in] output_filename name of the output file
     * \retval 0 on success
     * \retval non-zero value otherwise
     */
    int polyhedral_mesher(
        const std::string& input_filename, std::string output_filename
    ) {
        Mesh M_in;
        Mesh M_out;
	Mesh M_points;

        Logger::div("Polyhedral meshing");

        if(!mesh_load(input_filename, M_in)) {
            return 1;
        }
        
        if(M_in.cells.nb() == 0) {
            Logger::out("Poly") << "Mesh is not a volume" << std::endl;
            Logger::out("Poly") << "Trying to tetrahedralize" << std::endl;
            if(!mesh_tetrahedralize(M_in)) {
                return 1;
            }
            M_in.cells.compute_borders();
        }

	index_t dim = M_in.vertices.dimension();
	index_t spec_dim = CmdLine::get_arg_uint("poly:embedding_dim");
	if(spec_dim != 0 && spec_dim <= dim) {
	    dim = spec_dim;
	}
	
        CentroidalVoronoiTesselation CVT(&M_in, coord_index_t(dim));
        CVT.set_volumetric(true);

	if(CmdLine::get_arg("poly:points_file") == "") {

	    Logger::div("Generate random samples");
	    
	    CVT.compute_initial_sampling(
		CmdLine::get_arg_uint("remesh:nb_pts")
	    );

	    Logger::div("Optimize sampling");

	    try {
		index_t nb_iter = CmdLine::get_arg_uint("opt:nb_Lloyd_iter");
		ProgressTask progress("Lloyd", nb_iter);
		CVT.set_progress_logger(&progress);
		CVT.Lloyd_iterations(nb_iter);
	    }
	    catch(const TaskCanceled&) {
	    }

	    try {
		index_t nb_iter = CmdLine::get_arg_uint("opt:nb_Newton_iter");
		    ProgressTask progress("Newton", nb_iter);
		CVT.set_progress_logger(&progress);
		CVT.Newton_iterations(nb_iter);
	    }
	    catch(const TaskCanceled&) {
	    }
        
	    CVT.set_progress_logger(nil);
	} else {
	    if(!mesh_load(CmdLine::get_arg("poly:points_file"), M_points)) {
		return 1;
	    }
	    CVT.delaunay()->set_vertices(
		M_points.vertices.nb(), M_points.vertices.point_ptr(0)
	    );
	}

	CVT.RVD()->set_exact_predicates(true);
	{
	    SaveRVDCells callback(M_out);
	    std::string simplify = CmdLine::get_arg("poly:simplify");
	    if(simplify == "tets_voronoi_boundary") {
		double angle_threshold =
		    CmdLine::get_arg_double("poly:normal_angle_threshold");
		callback.set_simplify_boundary_facets(true, angle_threshold);
	    } else if(simplify == "tets_voronoi") {
		callback.set_simplify_voronoi_facets(true);
	    } else if(simplify == "tets") {
		callback.set_simplify_internal_tet_facets(true);		
	    } else if(simplify == "none") {
		callback.set_simplify_internal_tet_facets(false);
	    } else {
		Logger::err("Poly")
		    << simplify << " invalid cells simplification mode"
		    << std::endl;
	    }
	    callback.set_shrink(CmdLine::get_arg_double("poly:cells_shrink"));
	    callback.set_generate_ids(
		CmdLine::get_arg_bool("poly:generate_ids")
	    );
	    CVT.RVD()->for_each_polyhedron(callback);				
	}

	if(
	    FileSystem::extension(output_filename) == "mesh" ||
	    FileSystem::extension(output_filename) == "meshb"
	) {
	    Logger::warn("Poly")
		<< "Specified file format does not handle polygons"
		<< " (falling back to .obj)"
		<< std::endl;
	    output_filename =
		FileSystem::dir_name(output_filename) + "/" +
		FileSystem::base_name(output_filename) + ".obj";
	}

	if(
	    CmdLine::get_arg_bool("poly:generate_ids") &&
	    FileSystem::extension(output_filename) != "geogram" &&
	    FileSystem::extension(output_filename) != "geogram_ascii"
	) {
	    Logger::warn("Poly") << "Speficied file format does not handle ids"
				 << " (use .geogram or .geogram_ascii instead)"
				 << std::endl;
	}
	
	{
	    MeshIOFlags flags;
	    flags.set_attributes(MESH_ALL_ATTRIBUTES);
	    mesh_save(M_out, output_filename, flags);
	}
        
        return 0;
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
        CmdLine::import_arg_group("poly");
        
        std::vector<std::string> filenames;

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

	
        if(CmdLine::get_arg_bool("poly")) {
            return polyhedral_mesher(input_filename, output_filename);
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
	{
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

