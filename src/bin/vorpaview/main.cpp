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

#include "commands.h"
#include <geogram_gfx/basic/GLSL.h>
#include <geogram_gfx/glup_viewer/glup_viewer.h>
#include <geogram_gfx/mesh/mesh_gfx.h>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>

#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/stopwatch.h>

#ifdef GEO_OS_EMSCRIPTEN
#include <emscripten.h>
#endif

namespace {

    GEO::Mesh M(3);
    GEO::MeshGfx M_gfx;

    bool show_colors   = true;
    bool white_bg      = true;
    bool lighting      = true;
    
    bool show_vertices = false;
    GLfloat vertices_size = 1.0f;

    bool show_surface  = true;
    bool show_mesh = true;
    bool show_surface_borders = false;
    
    bool show_volume   = false;
    GLfloat cells_shrink = 0.0f;    
    bool show_hexes = true;    
    bool show_colored_cells = false;

    GLfloat OTM_time = 0.0f;
    GLfloat OTM_speed = 1.0f;

    /**
     * \brief The path to the files displayed in the file
     *  browser (Load menu).
     */
    std::string viewer_path = "./";
    
    /**
     * \brief Increments the time of the Optimal Transport animation.
     */
    void increment_time() {
        OTM_time += 0.05f;
        if(OTM_time > 1.0f) {
            OTM_time = 1.0f;
        }
    }

    /**
     * \brief Decrements the time of the Optimal Transport animation.
     */
    void decrement_time() {
        OTM_time -= 0.05f;
        if(OTM_time < 0.0f) {
            OTM_time = 0.0f;
        }
    }


    /**
     * \brief Increments cells shrinkage (shrinks more).
     */
    void increment_shrink() {
        cells_shrink += 0.05f;
        if(cells_shrink > 1.0f) {
            cells_shrink = 1.0f;
        }
    }

    /**
     * \brief Decrements cells shrinkage (shrinks less).
     */
    void decrement_shrink() {
        cells_shrink -= 0.05f;
        if(cells_shrink < 0.0f) {
            cells_shrink = 0.0f;
        }
    }

    /**
     * \brief Zooms in.
     * \details Zooming factor is 1.1x.
     */
    void zoom_in() {
        *glup_viewer_float_ptr(GLUP_VIEWER_ZOOM) *= 1.1f;
    }

    /**
     * \brief Zooms out.
     * \details De-zooming factor is (1/1.1)x.
     */
    void zoom_out() {
        *glup_viewer_float_ptr(GLUP_VIEWER_ZOOM) /= 1.1f;
    }


    void load_mesh(const std::string& filename);
    bool can_load(const std::string& filename);
    
    /**
     * \brief Initializes OpenGL objects.
     * \details Specifed as glup_viewer_set_init_func() callback.
     */
    void init() {
        GEO::Graphics::initialize();

        glup_viewer_disable(GLUP_VIEWER_BACKGROUND);        
        glup_viewer_set_background_color(1.0, 1.0, 1.0);
        glup_viewer_add_toggle('b', &white_bg, "white background");
        glup_viewer_add_toggle('c', &show_colors, "colors");

        glup_viewer_add_key_func('r', decrement_time, "Decrement time");
        glup_viewer_add_key_func('t', increment_time, "Increment time");

        glup_viewer_add_key_func('x', decrement_shrink, "Decrement shrink");
        glup_viewer_add_key_func('w', increment_shrink, "Increment shrink");
        
        glup_viewer_add_toggle('p', &show_vertices, "vertices");
        
        glup_viewer_add_toggle('S', &show_surface, "surface");
        glup_viewer_add_toggle('B', &show_surface_borders, "borders");
        glup_viewer_add_toggle('m', &show_mesh, "mesh");

        glup_viewer_add_toggle('V', &show_volume, "volume");
        glup_viewer_add_toggle('j',  &show_hexes, "hexes");
        glup_viewer_add_toggle('C', &show_colored_cells, "colored cells");

#ifdef GEO_OS_EMSCRIPTEN
        {
            std::vector<std::string> all_files;
            GEO::FileSystem::get_directory_entries("/",all_files);
            if(
                all_files.size() > 0 && 
                can_load(all_files[0])
            ) {
                 load_mesh(all_files[0]);
            }
        }
#endif        
    }

    /**
     * \brief Draws the current mesh according to current rendering mode.
     * \details Specifed as glup_viewer_set_display_func() callback.
     */
    void display() {

        if(glup_viewer_is_enabled(GLUP_VIEWER_IDLE_REDRAW)) {
            OTM_time = float(
                         sin(double(OTM_speed) * GEO::SystemStopwatch::now())
                       );
            OTM_time = 0.5f * (OTM_time + 1.0f);
        }
        
        M_gfx.set_lighting(lighting);
        M_gfx.set_time(double(OTM_time));
        
        if(M_gfx.mesh() != &M) {
            M_gfx.set_mesh(&M);
        }
        
        if(show_vertices) {
            M_gfx.set_points_size(vertices_size);
            M_gfx.draw_vertices();
        }

        if(white_bg) {
            glup_viewer_set_background_color(1.0, 1.0, 1.0);
            M_gfx.set_mesh_color(0.0, 0.0, 0.0);
        } else {
            glup_viewer_set_background_color(0.0, 0.0, 0.0);
            M_gfx.set_mesh_color(1.0, 1.0, 1.0);
        }
        
        if(show_colors) {
            if(M.cells.nb() == 0) {
                M_gfx.set_surface_color(0.5f, 0.75f, 1.0f);
                M_gfx.set_backface_surface_color(1.0f, 0.0f, 0.0f);
            } else {
                M_gfx.set_surface_color(0.7f, 0.0f, 0.0f);
                M_gfx.set_backface_surface_color(1.0f, 1.0f, 0.0f);
            }
        } else {
            if(white_bg) {
                M_gfx.set_surface_color(0.9f, 0.9f, 0.9f);
            } else {
                M_gfx.set_surface_color(0.1f, 0.1f, 0.1f);
            }
        }

        M_gfx.set_show_mesh(show_mesh);

        if(show_surface) {
            M_gfx.draw_surface();
        }
        
        if(show_surface_borders) {
            M_gfx.draw_surface_borders();
        }

        if(show_mesh) {
            M_gfx.draw_edges();
        }

        if(show_volume) {
            M_gfx.set_shrink(double(cells_shrink));
            M_gfx.set_draw_cells(GEO::MESH_HEX, show_hexes);
            if(show_colored_cells) {
                M_gfx.set_cells_colors_by_type();
            } else {
                M_gfx.set_cells_color(0.9f, 0.9f, 0.9f);
            }
            M_gfx.draw_volume();
        }
    }

    /**
     * \brief Gets the bounding box of a mesh animation.
     * \details In animated mode, the mesh animation is stored as a mesh with
     *  6d coordinates, that correspond to the geometric location
     *  at the vertices at time 0 and at time 1.
     * \param[in] M_in the mesh
     * \param[out] xyzmin a pointer to the three minimum coordinates
     * \param[out] xyzmax a pointer to the three maximum coordinates
     * \param[in] animate true if displaying a mesh animation
     */
    void get_bbox(
        const GEO::Mesh& M_in, double* xyzmin, double* xyzmax,
        bool animate
    ) {
        geo_assert(M_in.vertices.dimension() >= GEO::index_t(animate ? 6 : 3));
        for(GEO::index_t c = 0; c < 3; c++) {
            xyzmin[c] = GEO::Numeric::max_float64();
            xyzmax[c] = GEO::Numeric::min_float64();
        }

        for(GEO::index_t v = 0; v < M_in.vertices.nb(); ++v) {
            if(M_in.vertices.single_precision()) {
                const float* p = M_in.vertices.single_precision_point_ptr(v);
                for(GEO::coord_index_t c = 0; c < 3; ++c) {
                    xyzmin[c] = GEO::geo_min(xyzmin[c], double(p[c]));
                    xyzmax[c] = GEO::geo_max(xyzmax[c], double(p[c]));
                    if(animate) {
                        xyzmin[c] = GEO::geo_min(xyzmin[c], double(p[c+3]));
                        xyzmax[c] = GEO::geo_max(xyzmax[c], double(p[c+3]));
                    }
                }
            } else {
                const double* p = M_in.vertices.point_ptr(v);
                for(GEO::coord_index_t c = 0; c < 3; ++c) {
                    xyzmin[c] = GEO::geo_min(xyzmin[c], p[c]);
                    xyzmax[c] = GEO::geo_max(xyzmax[c], p[c]);
                    if(animate) {
                        xyzmin[c] = GEO::geo_min(xyzmin[c], p[c+3]);
                        xyzmax[c] = GEO::geo_max(xyzmax[c], p[c+3]);
                    }
                }
            }
        }
    }

    /**
     * \brief Inverts the normals of a mesh.
     * \details In color mode, this swaps the red and the blue sides.
     */
    void invert_normals() {
        for(GEO::index_t f = 0; f < M.facets.nb(); ++f) {
            M.facets.flip(f);
        }
        M_gfx.set_mesh(&M);
    }

    /**
     * \brief Tests whether a given file can be loaded as a mesh.
     * \param[in] filename a const reference to the complete path 
     *  to the file
     * \retval true if the file can be loaded as a mesh
     * \retval false otherwise
     */
    bool can_load(const std::string& filename) {
        std::string ext = GEO::FileSystem::extension(filename);
        return GEO::MeshIOHandlerFactory::has_creator(ext);
    }
    
    /**
     * \brief Loads a mesh from a file.
     */
    void load_mesh(const std::string& filename) {
        if(!GEO::FileSystem::is_file(filename)) {
            GEO::Logger::out("I/O") << "is not a file" << std::endl;
        }
        M_gfx.set_mesh(nil);
        
        GEO::MeshIOFlags flags;
        if(GEO::CmdLine::get_arg_bool("attributes")) {
            flags.set_attribute(GEO::MESH_FACET_REGION);
            flags.set_attribute(GEO::MESH_CELL_REGION);            
        } 
        if(!GEO::mesh_load(filename, M, flags)) {
            return;
        }

        if(
            GEO::FileSystem::extension(filename) == "obj6" ||
            GEO::FileSystem::extension(filename) == "tet6"
        ) {
            GEO::Logger::out("Vorpaview")
                << "Displaying optimal transport" << std::endl;

            glup_viewer_enable(GLUP_VIEWER_IDLE_REDRAW);
            
            M_gfx.set_animate(true);
            double xyzmin[3];
            double xyzmax[3];
            get_bbox(M, xyzmin, xyzmax, true);
            glup_viewer_set_region_of_interest(
                float(xyzmin[0]), float(xyzmin[1]), float(xyzmin[2]),
                float(xyzmax[0]), float(xyzmax[1]), float(xyzmax[2])
            );
        } else {
            M_gfx.set_animate(false);            
            M.vertices.set_dimension(3);
            double xyzmin[3];
            double xyzmax[3];
            get_bbox(M, xyzmin, xyzmax, false);
            glup_viewer_set_region_of_interest(
                float(xyzmin[0]), float(xyzmin[1]), float(xyzmin[2]),
                float(xyzmax[0]), float(xyzmax[1]), float(xyzmax[2])
            );
        }

        show_vertices = (M.facets.nb() == 0);
    }

    /**
     * \brief Loads a mesh from a file icon dropped into the window.
     * \details Specifed as glup_viewer_set_drag_drop_func() callback.
     */
    void dropped_file_cb(char* filename) {
        load_mesh(std::string(filename));
    }
}

/*************************************************************************/
            
/**
 * \brief Tests whether a directory should be skipped in the
 *  file browser.
 * \details Emscripten has some default directories in its
 *  root file system that are there just for compatibility 
 *  reasons and that are not likely to contain any 3D file.
 */
static bool skip_directory(const std::string& dirname) {
    GEO::geo_argused(dirname);
#ifdef GEO_OS_EMSCRIPTEN            
    if(
        dirname == "proc" ||
        dirname == "dev" ||
        dirname == "home" ||
        dirname == "tmp"
    ) {
        return true;
    }
#endif
    return false;
}

/**
 * \brief Converts a complete path to a file to a label
 *  displayed in the file browser.
 * \details Strips viewer_path from the input path.
 * \param[in] path the complete path, can be either a directory or
 *  a file
 * \return the label to be displayed in the menu
 */
static std::string path_to_label(const std::string& path) {
    if(GEO::String::string_starts_with(path, viewer_path)) {
        return path.substr(
            viewer_path.length(), path.length()-viewer_path.length()
        );
    }
    return path;
}

/**
 * \brief Recursively browes a directory and generates
 *  menu items.
 * \param[in] path the path to be browsed
 */
static void browse(const std::string& path) {
    std::vector<std::string> files;
    GEO::FileSystem::get_directory_entries(path,files);
    for(GEO::index_t i=0; i<files.size(); ++i) {
        if(GEO::FileSystem::is_directory(files[i])) {
            if(skip_directory(files[i])) {
                continue;
            }
            if(ImGui::BeginMenu(path_to_label(files[i]).c_str())) {
                browse(files[i]);
                ImGui::EndMenu();
            }
        } else {
            if(can_load(files[i])) {
                if(ImGui::MenuItem(path_to_label(files[i]).c_str())) {
                    load_mesh(files[i]);
                }
            }
        }
    }
}

/**
 * \brief Saves the current mesh to a file.
 * \details Under Emscripten, this opens the webbrowser's file
 *  dialog to save the file on the disk.
 * \param[in] filename the name of the file. Format is
 *  determined from the extension.
 */
static void save(const std::string& filename) {
    GEO::mesh_save(M, filename);
#ifdef GEO_OS_EMSCRIPTEN
    std::string command =
        "saveFileFromMemoryFSToDisk(\'" +
        filename +
        "\');"
        ;
    emscripten_run_script(command.c_str());
#endif    
}

/**
 * \brief Handles the "save as..." menu.
 * \details This generates a submenu for each supported file
 *  extension.
 */
static void save_as() {
    std::vector<std::string> extensions;
    GEO::MeshIOHandlerFactory::list_creators(extensions);
    for(GEO::index_t i=0; i<extensions.size(); ++i) {
        std::string filename = "out." + extensions[i];
        if(ImGui::MenuItem(filename.c_str())) {
            save(filename);
        }
    }
}

static bool win_object_properties=true;
static bool win_viewer_properties=true;
static bool win_console=false;

static const int MENU_HEIGHT = 20;
static const int PANE_WIDTH = 140;
static const int CONSOLE_HEIGHT = 200;

/**
 * \brief Displays and handles the object properties GUI pane.
 */
static void object_properties() {
    int w,h;
    glup_viewer_get_screen_size(&w, &h);
    
    ImGui::SetNextWindowPos(
        ImVec2(float(w-PANE_WIDTH), float(MENU_HEIGHT)),
        ImGuiSetCond_Always
    );
    ImGui::SetNextWindowSize(
        ImVec2(
            float(PANE_WIDTH),
            float(h-MENU_HEIGHT- (win_console ? CONSOLE_HEIGHT : 0))
        ),
        ImGuiSetCond_Always
    );

    ImGui::Begin(
        "Object", NULL,
        ImGuiWindowFlags_NoResize |
        ImGuiWindowFlags_NoMove |
        ImGuiWindowFlags_NoCollapse 
    );
    
    if(M.vertices.dimension() >= 6) {
        ImGui::Separator();
        ImGui::Checkbox(
            "Animate [a]",
            (bool*)glup_viewer_is_enabled_ptr(GLUP_VIEWER_IDLE_REDRAW)
        );
        ImGui::SliderFloat("spd.", &OTM_speed, 1.0f, 10.0f, "%.1f");
        ImGui::SliderFloat("t.", &OTM_time, 0.0f, 1.0f, "%.2f");
    }
    
    ImGui::Separator();    
    ImGui::Checkbox("Vertices [p]", &show_vertices);
    ImGui::SliderFloat("sz.", &vertices_size, 0.1f, 5.0f, "%.1f");

    if(M.facets.nb() != 0) {
        ImGui::Separator();
        ImGui::Checkbox("Surface [S]", &show_surface);
        if(show_surface) {
            ImGui::Checkbox("mesh [m]", &show_mesh);
            ImGui::Checkbox("borders [B]", &show_surface_borders);
        }
    }

    if(M.cells.nb() != 0) {
        ImGui::Separator();
        ImGui::Checkbox("Volume [V]", &show_volume);
        if(show_volume) {
            ImGui::SliderFloat(
                "shrk.", &cells_shrink, 0.0f, 1.0f, "%.2f"
            );        
            if(!M.cells.are_simplices()) {
                ImGui::Checkbox("colored cells [C]", &show_colored_cells);  
                ImGui::Checkbox("hexes [j]", &show_hexes);
            }
        }
    }
    
    ImGui::End();
}

/**
 * \brief Displays and handles the viewer properties GUI pane.
 */

static void viewer_properties() {

    int w,h;
    glup_viewer_get_screen_size(&w, &h);
    
    ImGui::SetNextWindowPos(
        ImVec2(0.0f, float(MENU_HEIGHT)),
        ImGuiSetCond_Always
    );

    float height = float(h - MENU_HEIGHT);
    if(win_console) {
        height -= float(CONSOLE_HEIGHT);
    }
    if(GEO::Command::current() != nil) {
        height /= 2.0f;
    }
    
    ImGui::SetNextWindowSize(
        ImVec2(float(PANE_WIDTH), height),
        ImGuiSetCond_Always
    );

    ImGui::Begin(
        "Viewer", NULL,
        ImGuiWindowFlags_NoResize |
        ImGuiWindowFlags_NoMove |
        ImGuiWindowFlags_NoCollapse 
    );

    ImGui::Checkbox(
        "Lighting [L]", &lighting
    );
    if(lighting) {
        ImGui::Checkbox(
            "edit [l]",
            (bool*)glup_viewer_is_enabled_ptr(GLUP_VIEWER_ROTATE_LIGHT)
        );
    }
    
    ImGui::Separator();
    ImGui::Checkbox(
        "Clipping [F1]", (bool*)glup_viewer_is_enabled_ptr(GLUP_VIEWER_CLIP)
    );
    if(glup_viewer_is_enabled(GLUP_VIEWER_CLIP)) {
        ImGui::Checkbox(
            "edit [F2]",
            (bool*)glup_viewer_is_enabled_ptr(GLUP_VIEWER_EDIT_CLIP)
        );
        ImGui::Checkbox(
            "fixed [F3]",
            (bool*)glup_viewer_is_enabled_ptr(GLUP_VIEWER_FIXED_CLIP)
        );
    }
    
    ImGui::Separator();
    ImGui::Text("Colors");
    ImGui::Checkbox("colors [c]", &show_colors);
    ImGui::Checkbox("white bkgnd [b]", &white_bg);

    ImGui::End();
}

/**
 * \brief Displays and handles the console.
 */
static void console() {
    int w,h;
    glup_viewer_get_screen_size(&w, &h);

    ImGui::SetNextWindowPos(
        ImVec2(0.0f, float(h-CONSOLE_HEIGHT+1)),
        ImGuiSetCond_Always
    );
    ImGui::SetNextWindowSize(
        ImVec2(float(w),float(CONSOLE_HEIGHT-1)),
        ImGuiSetCond_Always
    );
    glup_viewer_draw_console();
}

/*****************************************************************************/

/**
 * \brief Displays and handles the main menu.
 */
static void main_menu() {
    if(ImGui::BeginMainMenuBar()) {
        if(ImGui::BeginMenu("File")) {

#ifdef GEO_OS_EMSCRIPTEN
            ImGui::Text("To load a file,");
            ImGui::Text("use the \"Browse\"");
            ImGui::Text("button on the top");
            ImGui::Text("(or \"recent files\"");
            ImGui::Text("below)");
            ImGui::Separator();
            if(ImGui::BeginMenu("recent files...")) {
                browse(viewer_path);
                ImGui::EndMenu(); 
            }
#else
            if(ImGui::BeginMenu("load...")) {
                browse(viewer_path);
                ImGui::EndMenu();
            }
#endif            
            if(ImGui::BeginMenu("save as...")) {
                save_as();
                ImGui::EndMenu();
            }

            ImGui::Separator();
            if(ImGui::BeginMenu("about...")) {
                ImGui::Text(
                    "   This is VorpaView, the GEOGRAM demonstration program.\n"
                    "\n"
                    "\n"
                    "It demonstrates some geometric algorithms stemming from:\n"
                    "\n"
                    "  * ERC StG GOODSHAPE (ERC-StG-205693) \n"
                    "  * ERC PoC VORPALINE (ERC-PoC-334829) \n"
                    "\n"
                    "  ...as well as completely new algorithms.\n"
                    "\n"
#ifdef GEO_OS_EMSCRIPTEN
                    "  This version runs in your webbrowser using Emscripten.\n"
                    "To get a (much faster!) native executable and the sources,"
#endif
                    "\n"
                    "            see project's homepage:\n"
                    "  http://alice.loria.fr/software/geogram/doc/html/\n"
                    "\n"
                    "      (C)opyright 2006-2016, The ALICE project, Inria\n"
                );
                ImGui::EndMenu();
            }
            
#ifndef GEO_OS_EMSCRIPTEN                        
            ImGui::Separator();
            if(ImGui::MenuItem("quit [q]")) {
                glup_viewer_exit_main_loop();
            }
#endif

            ImGui::EndMenu();
        }

        if(ImGui::BeginMenu("Windows")) {
            ImGui::MenuItem("object properties", NULL, &win_object_properties);
            ImGui::MenuItem("viewer properties", NULL, &win_viewer_properties);
            ImGui::MenuItem("console", NULL, &win_console);
            if(ImGui::MenuItem(
                "show/hide GUI [T]", NULL,
                (bool*)glup_viewer_is_enabled_ptr(GLUP_VIEWER_TWEAKBARS)
            )) {
                glup_viewer_post_redisplay();
            }
            ImGui::EndMenu();
        }

        vorpaview_commands_menus();
        
        ImGui::EndMainMenuBar();            
    }
}


/**
 * \brief Drawns and manages the graphic user interface.
 */
static void overlay() {
    main_menu();
    if(win_object_properties) {
        object_properties();
    }
    if(win_viewer_properties) {
        viewer_properties();
        vorpaview_commands_gui();
    }
    if(win_console) {
        console();
    }
}

int main(int argc, char** argv) {

    GEO::initialize();
    GEO::Logger::instance()->set_quiet(false);

    // Import the arg groups needed by graphics.
    GEO::CmdLine::import_arg_group("standard");
    GEO::CmdLine::import_arg_group("algo");
    GEO::CmdLine::import_arg_group("gfx");

    // Import all the other arg groups needed by
    // the geogram algorithms used by the commands,
    // and indicate to the commands where the mesh is
    // as well as some flags.
    vorpaview_commands_init(
        M, M_gfx,
        show_vertices, show_surface, show_volume,
        win_console
    );
    
    // Default value for activating GLSL is determined by the
    // name of the executable, so that Windows users can
    // determine defaut behavior simply by changing the name
    // of the executable.
    const std::string& program_name = GEO::FileSystem::base_name(argv[0]);
    if(program_name == "vorpaview0") {
        GEO::CmdLine::set_arg("gfx:GLUP_profile","VanillaGL");
    }
    
    GEO::CmdLine::declare_arg(
        "attributes", true, "load mesh attributes"
    );

    GEO::CmdLine::declare_arg(
        "single_precision", true, "use single precision vertices (FP32)"
    );
    
    std::vector<std::string> filenames;
    if(!GEO::CmdLine::parse(argc, argv, filenames, "<filename>")) {
        return 1;
    }

    if(GEO::CmdLine::get_arg_bool("single_precision")) {
        M.vertices.set_single_precision();
    }
    
    if(filenames.size() == 1) {
        if(GEO::FileSystem::is_directory(filenames[0])) {
            viewer_path = filenames[0];
        } else if(can_load(filenames[0])) {
            load_mesh(filenames[0]);
        } else {
            GEO::Logger::err("Vorpaview")
                << filenames[0]
                << " is neither a directory nor a 3D file."
                << std::endl;
        }
    } 


    glup_viewer_set_window_title(
        (char*) "||||||(G)||E||(O)|(G)||R||/A\\|M|||||||"
    );
    glup_viewer_set_init_func(init);
    glup_viewer_set_display_func(display);
    glup_viewer_set_overlay_func(overlay);
    glup_viewer_set_drag_drop_func(dropped_file_cb);
    
    glup_viewer_add_toggle('L', &lighting, "toggle lighting");
    glup_viewer_add_toggle(
        'a', glup_viewer_is_enabled_ptr(GLUP_VIEWER_IDLE_REDRAW), "animate"
    );
    glup_viewer_add_key_func('n', invert_normals, "invert normals");
    
    glup_viewer_add_toggle(
        'T', glup_viewer_is_enabled_ptr(GLUP_VIEWER_TWEAKBARS), "tweakbars"
    );

    glup_viewer_add_key_func('z', zoom_in, "Zoom in");
    glup_viewer_add_key_func('Z', zoom_out, "Zoom out");

    if(GEO::CmdLine::get_arg_bool("gfx:full_screen")) {
       glup_viewer_enable(GLUP_VIEWER_FULL_SCREEN);
    }

    M_gfx.set_points_color(0.0, 1.0, 0.0);
    
    glup_viewer_main_loop(argc, argv);

    return 0;
}
