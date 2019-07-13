/*
 *  Copyright (c) 2012-2016, Bruno Levy
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

#include <geogram_gfx/gui/simple_application.h>
#include <geogram_gfx/gui/geogram_logo_256.xpm>
#include <geogram/basic/file_system.h>
#include <geogram/basic/string.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/command_line.h>
#include <geogram/bibliography/bibliography.h>

#include <geogram_gfx/full_screen_effects/ambient_occlusion.h>
#include <geogram_gfx/full_screen_effects/unsharp_masking.h>


#include <geogram_gfx/gui/colormaps/french.xpm>
#include <geogram_gfx/gui/colormaps/black_white.xpm>
#include <geogram_gfx/gui/colormaps/viridis.xpm>
#include <geogram_gfx/gui/colormaps/rainbow.xpm>
#include <geogram_gfx/gui/colormaps/cei_60757.xpm>
#include <geogram_gfx/gui/colormaps/inferno.xpm>
#include <geogram_gfx/gui/colormaps/magma.xpm>
#include <geogram_gfx/gui/colormaps/parula.xpm>
#include <geogram_gfx/gui/colormaps/plasma.xpm>
#include <geogram_gfx/gui/colormaps/blue_red.xpm>

#ifdef GEOGRAM_WITH_LUA
#include <geogram_gfx/lua/lua_glup.h>
#include <geogram_gfx/lua/lua_simple_application.h>
#include <geogram_gfx/lua/lua_imgui.h>
#include <geogram/lua/lua_io.h>

extern "C" {
#include <geogram/third_party/lua/lua.h>    
#include <geogram/third_party/lua/lauxlib.h>
#include <geogram/third_party/lua/lualib.h>
}

#endif

#ifdef GEO_OS_EMSCRIPTEN
#include <emscripten.h>
#endif

namespace {
#include <geogram_gfx/gui/gui_state_bigfont_v.h>
#include <geogram_gfx/gui/gui_state_bigfont_h.h>
#include <geogram_gfx/gui/gui_state.h>
}

/******************************************************************************/

namespace {

   /**
    * \brief Converts a complete path to a file to a label
    *  displayed in the file browser.
    * \details Strips viewer_path from the input path.
    * \param[in] path the complete path, can be either a directory or
    *  a file
    * \return the label to be displayed in the menu
    */
    std::string path_to_label(
        const std::string& viewer_path, const std::string& path
    ) {
        std::string result = path;
        if(GEO::String::string_starts_with(result, viewer_path)) {
            result = result.substr(
                viewer_path.length(), result.length()-viewer_path.length()
            );
        }
        return result;
    }
}

/******************************************************************************/

namespace GEO {

    SimpleApplication::SimpleApplication(
	const std::string& name
    ) :
	Application(name),
	text_editor_(&text_editor_visible_)
    {
	lighting_ = true;
	edit_light_ = false;
	clipping_ = false;
	clip_mode_ = GLUP_CLIP_WHOLE_CELLS;
	edit_clip_ = false;
	fixed_clip_ = false;
	background_color_ = vec4f(0.0f, 0.0f, 0.0f, 1.0f);	
	effect_ = GLenum(0);
	
	filename_[0] = '\0';
	geogram_logo_texture_ = 0;
	viewer_properties_visible_ = true;
	object_properties_visible_ = true;
#ifdef GEO_OS_ANDROID
	console_visible_           = false;	
#else	
	console_visible_           = true;
#endif
	use_text_editor_           = false;
	text_editor_visible_       = false;
	menubar_visible_           = true;

        console_ = new Console(&console_visible_);
	console_->hide_command_prompt();
	text_editor_.set_fixed_layout(false);
        status_bar_ = new StatusBar;

	add_key_func("q", [this]() { stop(); });
	add_key_func("z", [this]() { zoom_in(); });
	add_key_func("Z", [this]() { zoom_out(); });
	add_key_func("H", [this]() { home(); });		
	add_key_toggle("L",   &lighting_);
	add_key_toggle("l",   &edit_light_);
	add_key_toggle("F1",  &clipping_);
	add_key_toggle("F2",  &edit_clip_);
	add_key_toggle("F3",  &fixed_clip_);
	add_key_toggle("a",   animate_ptr());
	add_key_toggle("F6",  &text_editor_visible_);				
	add_key_toggle("F7",  &viewer_properties_visible_);
	add_key_toggle("F8",  &object_properties_visible_);
	add_key_toggle("F9",  &console_visible_);
	add_key_toggle("F12", &menubar_visible_);
	set_region_of_interest(
	    0.0, 0.0, 0.0, 1.0, 1.0, 1.0
	);

	object_translation_ = vec3(0.0, 0.0, 0.0);

	mouse_op_ = MOUSE_NOOP;
	mouse_target_ = MOUSE_NOTARGET;

	three_D_ = true;
	zoom_ = 1.0;
	zoom_down_ = 1.0;

#ifdef GEOGRAM_WITH_LUA	
	lua_error_occured_ = false;
	lua_state_ = luaL_newstate();
	luaL_openlibs(lua_state_);
	init_lua_io(lua_state_);
	init_lua_glup(lua_state_);
	init_lua_simple_application(lua_state_);		
	init_lua_imgui(lua_state_);
#else
	lua_error_occured_ = false;
	lua_state_ = nullptr;
#endif
	geo_cite_with_info(
	    "WEB:ImGUI",
	    "Used to create the GUI of GEOGRAM utilities "
	    "(vorpaview, geobox, geocod)."
	);
    }

    SimpleApplication::~SimpleApplication() {
#ifdef GEOGRAM_WITH_LUA
	if(lua_state_ != nullptr) {
	    lua_close(lua_state_);
	    lua_state_ = nullptr;
	}
#endif	
    }
    
    void SimpleApplication::home() {
	zoom_ = 1.0;
	object_translation_ = vec3(0.0, 0.0, 0.0);
	object_rotation_.reset();
	light_rotation_.reset();
	clip_rotation_.reset();
	clip_translation_ = vec3(0.0, 0.0, 0.0);
	clipping_ = false;
	edit_clip_ = false;
	fixed_clip_ = false;
	edit_light_ = false;
    }

    void SimpleApplication::add_key_func(
	const std::string& key, std::function<void()> cb
    ) {
	key_funcs_[key] = cb;
    }
    
    void SimpleApplication::add_key_toggle(
	const std::string& key, bool* p_val
    ) {
	add_key_func(
	    key, [p_val]() { *p_val = !*p_val; }
	);
    }

    void SimpleApplication::char_callback(unsigned int c){
	Application::char_callback(c);
	if(text_editor_visible_) {
	    return;
	}
	std::string k=" ";
	k[0] = char(c);
	auto F = key_funcs_.find(k);
	if(F != key_funcs_.end()) {
	    F->second();
	}
    }

    void SimpleApplication::key_callback(
	int key, int scancode, int action, int mods
    ) {
	Application::key_callback(key, scancode, action, mods);
	if(action == 0) {
	    auto F = key_funcs_.find(std::string(key_to_string(key)));
	    if(F != key_funcs_.end()) {
		F->second();
	    }
	}
    }
    
    void SimpleApplication::set_style(const std::string& style) {
	Application::set_style(style);
	if(String::string_starts_with(style, "Light")) {
	    background_color_ = vec4f(1.0f, 1.0f, 1.0f, 1.0f);
	} else {
	    background_color_ = vec4f(0.0f, 0.0f, 0.0f, 1.0f);	    
	}
    }

    void SimpleApplication::set_region_of_interest(
	double xmin, double ymin, double zmin,
	double xmax, double ymax, double zmax
    ) {
	roi_.xyz_min[0] = xmin;
	roi_.xyz_min[1] = ymin;
	roi_.xyz_min[2] = zmin;
	roi_.xyz_max[0] = xmax;
	roi_.xyz_max[1] = ymax;
	roi_.xyz_max[2] = zmax;
	roi_radius_ = sqrt(
	    0.25 * (xmax - xmin) * (xmax - xmin) +
	    0.25 * (ymax - ymin) * (ymax - ymin) +
	    0.25 * (zmax - zmin) * (zmax - zmin)
	);
    }

    void SimpleApplication::get_region_of_interest(
	double& xmin, double& ymin, double& zmin,
	double& xmax, double& ymax, double& zmax
    ) const {
	xmin = roi_.xyz_min[0];
	ymin = roi_.xyz_min[1];
	zmin = roi_.xyz_min[2];
	xmax = roi_.xyz_max[0];
	ymax = roi_.xyz_max[1];
	zmax = roi_.xyz_max[2];
    }
    
    void SimpleApplication::draw_gui() {
	draw_menu_bar();
	draw_dock_space();
	draw_viewer_properties_window();
	draw_object_properties_window();
	draw_console();
	draw_command_window();
	if(text_editor_visible_) {
	    text_editor_.draw();
	}
	if(
	    ImGui::FileDialog("##load_dlg", filename_, geo_imgui_string_length)
	) {
	    load(filename_);
	}
	if(
	    ImGui::FileDialog("##save_dlg", filename_, geo_imgui_string_length)
	) {
	    save(filename_);
	}
	if(status_bar_->active()) {
	    float w = float(get_width());
	    float h = float(get_height());
	    float STATUS_HEIGHT = 50.0f * float(scaling());
	    ImGui::SetNextWindowPos(
		ImVec2(0.0f, h-STATUS_HEIGHT),
		ImGuiCond_Always
	    );
	    ImGui::SetNextWindowSize(
		ImVec2(w,STATUS_HEIGHT-1.0f),
		ImGuiCond_Always
	    );
	    status_bar_->draw();
	}
    }

    void SimpleApplication::draw_scene_begin() {
	glClearColor(
	    background_color_.x,
	    background_color_.y,
	    background_color_.z,
	    background_color_.w
	);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);

	glViewport(0, 0, int(get_width()), int(get_height()));

	double zScreen = 5.0; // screen projection plane 
	{
	    glupMatrixMode(GLUP_PROJECTION_MATRIX);
	    glupLoadIdentity();
	    
	    double aspect = double(get_width()) / double(get_height());
	    double zNear = 1.0;   // near clipping plane.
	    double zFar = 10.0;   // far clipping plane.
	    if(three_D_) {
		double camera_aperture = 9.0; // field of view in degrees.
		double view_max_size = zScreen * tan(
		    (camera_aperture * M_PI / 180.0) / 2.0
		);
		double right;
		double top;
		if(aspect < 1.0) {
		    top = view_max_size;
		    right = top * aspect;
		} else {
		    right = view_max_size;
		    top = right / aspect;
		}
		right /= zoom_;
		top   /= zoom_;
		glupFrustum(-right, right, -top, top, zNear, zFar);
	    } else {
		double x = 1.0 / zoom_;
		double y = 1.0 / zoom_;
		if(aspect > 1.0) {
		    x *= aspect;
		} else {
		    y /= aspect;
		}
		glupOrtho(-x, x, -y, y, zNear, zFar);
	    }
	}

	static vec3 light0 = normalize(vec3(1.0, 1.0, 1.0));
	vec3 light = transform_vector(
	    light0, light_rotation_.get_value()
	);
	glupLightVector3f(float(light.x), float(light.y), float(light.z));
    
	glupMatrixMode(GLUP_MODELVIEW_MATRIX);
	glupLoadIdentity();
	glupTranslated(0.0, 0.0, -zScreen);

	// Clipping
	{
	    static double clip_eqn[4];
	    clip_eqn[0] = 0.0;
	    clip_eqn[1] = 0.0;
	    clip_eqn[2] = 0.0;
	    clip_eqn[3] = 0.0;
	    glupClipPlane(clip_eqn);
	    glupDisable(GLUP_CLIPPING);
	    
	    if(clipping_) {
		glupPushMatrix();
		glupTranslated(
		    clip_translation_.x,
		    clip_translation_.y,
		    clip_translation_.z
		);
		glupMultMatrix(clip_rotation_.get_value());

		{
		    glupSetColor3f(
			GLUP_FRONT_AND_BACK_COLOR,
			1.0f - background_color_.x,
			1.0f - background_color_.y,
			1.0f - background_color_.z
		    );
		    float sq_w = 1.25f / float(zoom_);

		    GLboolean vertex_colors_save =
			glupIsEnabled(GLUP_VERTEX_COLORS);
		    GLboolean texturing_save = glupIsEnabled(GLUP_TEXTURING);
		    glupDisable(GLUP_VERTEX_COLORS);
		    glupDisable(GLUP_TEXTURING);

		    // Draw the cross 
		    glupBegin(GLUP_LINES);
		    glupVertex3f(-sq_w, 0.0f, 0.0f);
		    glupVertex3f(sq_w, 0.0f, 0.0f);
		    glupVertex3f(0.0f, -sq_w, 0.0f);
		    glupVertex3f(0.0f, sq_w, 0.0f);

		    // Draw the square around the cross 
		    for(index_t i=0; i<3; ++i) {
			glupVertex3f(sq_w, -sq_w, 0.0f);
			glupVertex3f(sq_w, sq_w, 0.0f);
			glupVertex3f(sq_w, sq_w, 0.0f);            
			glupVertex3f(-sq_w, sq_w, 0.0f);
			glupVertex3f(-sq_w, sq_w, 0.0f);            
			glupVertex3f(-sq_w, -sq_w, 0.0f);
			glupVertex3f(-sq_w, -sq_w, 0.0f);
			glupVertex3f(sq_w, -sq_w, 0.0f);
			sq_w = sq_w * 1.01f;
		    }
		    glupEnd();
		    if(vertex_colors_save) {
			glupEnable(GLUP_VERTEX_COLORS);
		    }
		    if(texturing_save) {
			glupEnable(GLUP_TEXTURING);
		    }
		}
		clip_eqn[0] = 0.0;
		clip_eqn[1] = 0.0;
		clip_eqn[2] = -1.0;
		clip_eqn[3] = 0.0;
		glupEnable(GLUP_CLIPPING); 
		glupClipPlane(clip_eqn);
		glupClipMode(clip_mode_);
		glupPopMatrix();
	    } 
	}
	
	glupTranslate(object_translation_);
	glupMultMatrix(object_rotation_.get_value());
	
	glupScaled(
	    1.5 / roi_radius_, 1.5 / roi_radius_, 1.5 / roi_radius_
	);
	glupTranslated(
	    -0.5 * (roi_.xyz_min[0] + roi_.xyz_max[0]),
	    -0.5 * (roi_.xyz_min[1] + roi_.xyz_max[1]),
	    -0.5 * (roi_.xyz_min[2] + roi_.xyz_max[2])
	);

	if(lighting_) {
	    glupEnable(GLUP_LIGHTING);
	} else {
	    glupDisable(GLUP_LIGHTING);	    
	}

	// Save transform for picking
	{
	    glGetIntegerv(GL_VIEWPORT, viewport_);
	    // Note: OpenGL uses column-major order for matrices
	    // (thus what we get is the transpose of each matrix)
	    glupGetMatrixdv(GLUP_MODELVIEW_MATRIX, modelview_transpose_.data());
	    glupGetMatrixdv(GLUP_PROJECTION_MATRIX, project_transpose_.data());
	}
    }

    vec3 SimpleApplication::project(const vec3& p) {
	vec3 result;
	glupProject(
	    p.x, p.y, p.z,
	    modelview_transpose_.data(), project_transpose_.data(), viewport_,
	    &result.x, &result.y, &result.z
	);
	return result;
    }

    vec3 SimpleApplication::unproject(const vec3& p) {
	vec3 result;
	glupUnProject(
	    p.x, p.y, p.z,
	    modelview_transpose_.data(), project_transpose_.data(), viewport_,
	    &result.x, &result.y, &result.z
	);
	return result;
    }

    vec2 SimpleApplication::unproject_2d(const vec2& p) {
	double z = project(vec3(0.0, 0.0, 0.0)).z;
	vec3 result3d = unproject(vec3(p.x, p.y, z));
	return vec2(result3d.x, result3d.y);
    }
    
    void SimpleApplication::draw_scene_end() {
    }

    void SimpleApplication::draw_scene() {
    }
    
    void SimpleApplication::draw_graphics() {
	if(!full_screen_effect_.is_null()) {
	    full_screen_effect_->pre_render(get_width(), get_height());
	}
	draw_scene_begin();
	draw_scene();
	draw_scene_end();
	if(!full_screen_effect_.is_null()) {
	    full_screen_effect_->post_render();
	}
    }

    void SimpleApplication::draw_viewer_properties_window() {
	if(!viewer_properties_visible_) {
	    return;
	}
	if(ImGui::Begin(
	       (icon_UTF8("camera")+" Viewer").c_str(),
	       &viewer_properties_visible_)
	) {
	    draw_viewer_properties();
	}
	ImGui::End();
    }

    void SimpleApplication::draw_viewer_properties() {
	if(ImGui::Button(
	       (icon_UTF8("home") + " Home [H]").c_str(), ImVec2(-1.0, 0.0))
	) {
	    home();
	}
	ImGui::Separator();
	ImGui::Checkbox("Animate [a]", animate_ptr());
	if(three_D_) {
	    ImGui::Checkbox("Lighting [L]", &lighting_);
	    if(lighting_) {
		ImGui::Checkbox("Edit light [l]", &edit_light_);	    
	    }
	    ImGui::Separator();
	    ImGui::Checkbox("Clipping [F1]", &clipping_);
	    if(clipping_) {
		ImGui::Combo(
		    "mode", (int*)&clip_mode_,                
		    "std. GL\0cells\0straddle\0slice\0\0"
		);
		ImGui::Checkbox(
		    "edit clip [F2]", &edit_clip_
		);
		ImGui::Checkbox(
		    "fixed clip [F3]", &fixed_clip_
		);
	    }
	    ImGui::Separator();
	}
	ImGui::ColorEdit3WithPalette("Background", background_color_.data());
#ifndef GEO_OS_ANDROID
	if(three_D_) {
	    if(ImGui::Combo("sfx", (int*)&effect_, "none\0SSAO\0cartoon\0\0")) {
		switch(effect_) {
		    case 0:
			full_screen_effect_.reset();
			break;
		    case 1:
			full_screen_effect_ = new AmbientOcclusionImpl();
			break;
		    case 2:
			full_screen_effect_ = new UnsharpMaskingImpl();
			break;
		}
	    }
	}
#endif	
    }

    void SimpleApplication::draw_object_properties_window() {
	if(!object_properties_visible_) {
	    return;
	}
	if(ImGui::Begin(
	       (icon_UTF8("edit")+" Object").c_str(),
	       &object_properties_visible_)
	) {
	    draw_object_properties();
	}
	ImGui::End();
    }

    void SimpleApplication::draw_object_properties() {
    }

    void SimpleApplication::draw_command_window() {
	if(Command::current() == nullptr) {
	    return;
	}
	if(!Command::current()->is_visible()) {
	    Command::reset_current();
	    return;
	}
        if(ImGui::Begin(
	    "Command",
            Command::current()->is_visible_ptr()
        )) {
	    Command::current()->draw();
	}
        ImGui::End();
    }

    void SimpleApplication::draw_console() {
	console_->draw(&console_visible_);
    }
    
    void SimpleApplication::draw_menu_bar() {
	if(!menubar_visible_) {
	    return;
	}
        if(ImGui::BeginMainMenuBar()) {
            if(ImGui::BeginMenu("File")) {
                if(supported_read_file_extensions() != "") {
                    draw_load_menu();
                }
#ifndef GEO_OS_EMSCRIPTEN		
		if(current_file_ != "") {
		    if(ImGui::MenuItem(icon_UTF8("save") + " Save")) {
			if(save(current_file_)) {
			    Logger::out("I/O") << "Saved "
					       << current_file_ << std::endl;
			} else {
			    Logger::out("I/O") << "Could not save "
					       << current_file_ << std::endl;
			}
		    }
		}
#endif		
                if(supported_write_file_extensions() != "") {
                    draw_save_menu();
                }
		draw_fileops_menu();
#ifndef GEO_OS_EMSCRIPTEN                        
                ImGui::Separator();
                if(ImGui::MenuItem(icon_UTF8("door-open") + " quit",
				   "[q]", false, true)
		) {
                    this->stop();
                }
#endif
                draw_about();
                ImGui::EndMenu();
            }
	    if(ImGui::BeginMenu("Windows")) {
		draw_windows_menu();
		ImGui::EndMenu();
	    }
            draw_application_menus();
            
            ImGui::EndMainMenuBar();            
	}
    }

    void SimpleApplication::draw_load_menu() {
#ifdef GEO_OS_EMSCRIPTEN
            ImGui::Text("To load a file,");
            ImGui::Text("use the \"Browse\"");
            ImGui::Text("button on the top");
            ImGui::Text("(or \"recent files\"");
            ImGui::Text("below)");
            ImGui::Separator();
            if(ImGui::BeginMenu("Recent files...")) {
                browse(path_);
                ImGui::EndMenu(); 
            }
#else
	    if(ImGui::MenuItem(icon_UTF8("folder-open") + " Load...")) {
		ImGui::OpenFileDialog(
		    "##load_dlg",
		    supported_read_file_extensions().c_str(),
		    filename_,
		    ImGuiExtFileDialogFlags_Load		
		);
	    }
#endif        
    }

    void SimpleApplication::draw_save_menu() {
#ifdef GEO_OS_EMSCRIPTEN
        if(ImGui::BeginMenu(icon_UTF8("save") + " Save as...")) {
	    ImGui::MenuItem("Supported extensions:", nullptr, false, false);
            std::vector<std::string> extensions;
            String::split_string(
                supported_write_file_extensions(), ';', extensions
            );
            for(index_t i=0; i<extensions.size(); ++i) {
		ImGui::MenuItem(
		    " ." + extensions[i], nullptr, false, false
		);	    		
	    }
	    ImGui::Separator();
	    static char buff[geo_imgui_string_length];
	    if(current_file_ != "") {
		strcpy(buff, current_file_.c_str());
	    } else if (extensions.size() != 0) {
		strcpy(buff, ("out." + extensions[0]).c_str());		
	    }

	    if(ImGui::InputText(
		   "##MenuFileName",buff,geo_imgui_string_length,
		   ImGuiInputTextFlags_EnterReturnsTrue)
	    ) {
		current_file_ = buff;
		if(String::string_starts_with(current_file_, "/")) {
		    current_file_ = current_file_.substr(
			1,current_file_.length()-1
		    );
		}
		if(save(current_file_)) {
		    std::string command =
			"saveFileFromMemoryFSToDisk(\'" +
			current_file_ +
			"\');" ;
		    emscripten_run_script(command.c_str());
		}
	    }
            ImGui::EndMenu();
        }
#else        
        if(ImGui::MenuItem(icon_UTF8("save") + " Save as...")) {
	    ImGui::OpenFileDialog(
		"##save_dlg",
		supported_write_file_extensions().c_str(),
		filename_,
		ImGuiExtFileDialogFlags_Save
	    );
        }
#endif        
    }

    void SimpleApplication::draw_fileops_menu() {
    }

    void SimpleApplication::draw_about() {
        ImGui::Separator();
        if(ImGui::BeginMenu(icon_UTF8("info") + " About...")) {
            ImGui::Text("%s : a GEOGRAM application", name().c_str());
	    float sz = float(280.0 * std::min(scaling(), 2.0));
            ImGui::Image(
                convert_to_ImTextureID(geogram_logo_texture_),
                ImVec2(sz, sz)
            );
            ImGui::Text("\n");            
            ImGui::Separator();
            ImGui::Text("\n");
            ImGui::Text("GEOGRAM website: ");
            ImGui::Text("http://alice.loria.fr/software/geogram");

            ImGui::Text("\n");
            ImGui::Separator();
            ImGui::Text(
                "%s",
                (
                    "GEOGRAM version:" +
                    Environment::instance()->get_value("version")
                ).c_str()
            );
            ImGui::EndMenu();
        }
    }

    void SimpleApplication::draw_windows_menu() {
	if(use_text_editor_) {
	    ImGui::MenuItem(
		"text editor", "[F6]", &text_editor_visible_
	    );
	}
	ImGui::MenuItem(
	    "viewer properties", "[F7]", &viewer_properties_visible_
	);
	ImGui::MenuItem(
	    "object properties", "[F8]", &object_properties_visible_
	);
	ImGui::MenuItem(
	    "console", "[F9]", &console_visible_
	);
#ifndef GEO_OS_ANDROID	
	ImGui::MenuItem(
	    "menubar", "[F12]", &menubar_visible_
	);
#endif	
	ImGui::Separator();
	if(ImGui::BeginMenu("Font size")) {
	    static index_t font_sizes[] = {
		10, 12, 18, 22, 24, 30, 32, 36, 40, 46, 50, 56
	    };
	    for(index_t i=0; i<sizeof(font_sizes)/sizeof(int); ++i) {
		bool selected = (get_font_size() == font_sizes[i]);
		if(ImGui::MenuItem(
		       String::to_string(
			   font_sizes[i]).c_str(), nullptr, &selected)
		    ) {
		    set_font_size(font_sizes[i]);
		}
	    }
	    ImGui::EndMenu();
	}
	if(ImGui::BeginMenu("Style")) {
	    std::vector<std::string> styles;
	    String::split_string(get_styles(), ';', styles);
	    for(index_t i=0; i<styles.size(); ++i) {
		bool selected = (get_style() == styles[i]);
		if(ImGui::MenuItem(styles[i].c_str(), nullptr, &selected)) {
		    set_style(styles[i]);
		}
	    }
	    ImGui::EndMenu();
	}
	ImGui::Separator();
	if(ImGui::MenuItem("Restore default layout")) {
	    set_default_layout();
	}
	if(CmdLine::get_arg_bool("gui:expert")) {
	    if(ImGui::MenuItem("Test android vertical layout")) {
		set_font_size(56);
		set_gui_state(default_layout_android_vertical());
	    }
	    if(ImGui::MenuItem("Test android horizontal layout")) {
		set_font_size(56);
		set_gui_state(default_layout_android_horizontal());		
	    }
	    if(ImGui::MenuItem("Export gui state to C++")) {
		std::string filename =
		    FileSystem::get_current_working_directory() +
		    "/gui_state.h";
		Logger::out("GUI")
		    << "Exporting current GUI state to C++ file: "
		    << filename
		    << std::endl;
		std::string state = get_gui_state();
		std::ofstream out("gui_state.h");
		out << "// Serialized ImGui windows docking configuration"
		    << std::endl;
		out << "// generated using <geogram_program> gui:expert=true"
		    << std::endl;
		out << "// then Windows->Export gui state to C++"
		    << std::endl;
		out << "const char gui_state[] = u8\""
		    << state
		    << "\";" << std::endl;
	    }
	}
    }
    
    void SimpleApplication::set_default_layout() {
	set_gui_state(default_layout());
    }

    const char* SimpleApplication::default_layout_android_vertical() const {
	return gui_state_v;
    }

    const char* SimpleApplication::default_layout_android_horizontal() const {
	return gui_state_h;
    }
    
    const char* SimpleApplication::default_layout() const {
#ifdef GEO_OS_ANDROID
	if(get_height() >= get_width()) {
	    return default_layout_android_vertical();
	} else {
	    return default_layout_android_horizontal();
	}
#else
	return gui_state;
#endif	
    }

    
    void SimpleApplication::resize(index_t w, index_t h) {
	Application::resize(w,h);
#ifdef GEO_OS_ANDROID	
	set_default_layout();
#endif	
    }
    
    void SimpleApplication::draw_application_menus() {
    }

    void SimpleApplication::post_draw() {
	Command::flush_queue();
    }

    void SimpleApplication::mouse_button_callback(
	int button, int action, int mods, int source
    ) {
	geo_argused(mods);

	// Swap "buttons" if using fingers (it is more
	// natural to do the rotation with one finger,
	// zoom with two fingers and then translation)
	// Same thing if editing light and event source
	// is stylus.
	if(
	    source == EVENT_SOURCE_FINGER ||
	    (lighting_ && source == EVENT_SOURCE_STYLUS && edit_light_)
	) {
	    if(button == 0) {
		button = 1;
	    } else if(button == 1) {
		button = 0;
	    }
	}
	
	if(action == EVENT_ACTION_DOWN) {
	    mouse_down_xy_ = mouse_xy_;
	    if(button == 1) {
		if(three_D_) {
		    mouse_op_ = MOUSE_ROTATE;
		} else {
		    mouse_op_ = MOUSE_TRANSLATE;
		}
	    } else if(button == 0) {
		mouse_op_ = MOUSE_TRANSLATE;
	    } else if(button == 2) {
		mouse_op_ = MOUSE_ZOOM;
		zoom_down_ = zoom_;
	    }
	    if(clipping_ && edit_clip_) {
		mouse_target_ = MOUSE_CLIP;
	    } else if(lighting_ && edit_light_) {
		mouse_target_ = MOUSE_LIGHT;
	    } else {
		mouse_target_ = MOUSE_OBJECT;
	    }
	    if(three_D_ && mouse_op_ == MOUSE_ROTATE) {
		switch(mouse_target_) {
		    case MOUSE_NOTARGET:
			break;
		    case MOUSE_OBJECT:
			object_rotation_.grab(mouse_xy_);
			if(fixed_clip_) {
			    clip_rotation_.grab(mouse_xy_); 
			}
			break;
		    case MOUSE_LIGHT:
			light_rotation_.grab(mouse_xy_);		    
			break;
		    case MOUSE_CLIP:
			clip_rotation_.grab(mouse_xy_);		    
			break;
		}
	    }
	} else if(action == EVENT_ACTION_UP) {
	    if(three_D_ && mouse_op_ == MOUSE_ROTATE) {
		switch(mouse_target_) {
		    case MOUSE_NOTARGET:
			break;
		    case MOUSE_OBJECT:
			object_rotation_.release(mouse_xy_);
			if(fixed_clip_) {
			    clip_rotation_.release(mouse_xy_); 
			}
			break;
		    case MOUSE_LIGHT:
			light_rotation_.release(mouse_xy_);		    
			break;
		    case MOUSE_CLIP:
			clip_rotation_.release(mouse_xy_);		    
			break;
		}
	    }
	    mouse_op_ = MOUSE_NOOP;
	    mouse_target_ = MOUSE_NOTARGET;
	}
    }

    void SimpleApplication::cursor_pos_callback(
	double x, double y, int source
    ) {
	geo_argused(source);
	mouse_xy_ = vec2(
	    double(x) / double(get_width()),
	    double(y) / double(get_height())
	);
	mouse_xy_ *= 2.0;
	mouse_xy_ -= vec2(1.0, 1.0);
	if(three_D_ && mouse_op_ == MOUSE_ROTATE) {
	    switch(mouse_target_) {
		case MOUSE_NOTARGET:
		    break;
		case MOUSE_OBJECT:
		    object_rotation_.drag(mouse_xy_);
		    if(fixed_clip_) {
			clip_rotation_.drag(mouse_xy_);			
		    }
		    break;
		case MOUSE_LIGHT:
		    light_rotation_.drag(mouse_xy_);		    
		    break;
		case MOUSE_CLIP:
		    clip_rotation_.drag(mouse_xy_);		    
		    break;
	    }
	} else if(mouse_op_ == MOUSE_ZOOM && mouse_target_ == MOUSE_OBJECT) {
	    double R = mouse_xy_.y - mouse_down_xy_.y;
	    double fact = (1.0 + R);
	    fact = std::min(fact, 2.0);
	    fact = std::max(fact, 0.1);
	    zoom_ = zoom_down_ * fact;
	} else if( mouse_op_ == MOUSE_TRANSLATE) {
	    double dx = mouse_xy_.x - mouse_down_xy_.x;	    	    
	    double dy = mouse_xy_.y - mouse_down_xy_.y;
	    if(mouse_target_ == MOUSE_OBJECT) {
		object_translation_.x += 2.0 * dx / zoom_;
		object_translation_.y -= 2.0 * dy / zoom_;
		if(fixed_clip_) {
		    clip_translation_.x += 2.0 * dx / zoom_;
		    clip_translation_.y -= 2.0 * dy / zoom_;
		}
	    } else if(mouse_target_ == MOUSE_CLIP) {
		clip_translation_.x += 2.0 * dx / zoom_;
		clip_translation_.y -= 2.0 * dy / zoom_;
	    }
	    mouse_down_xy_ = mouse_xy_;
	}
    }

    void SimpleApplication::scroll_callback(double xoffset, double yoffset) {
	geo_argused(xoffset);
	double dy = -40.0*double(yoffset) / double(get_height());
	zoom_ *= (1.0 + dy);
    }

    bool SimpleApplication::save(const std::string& filename) {
        Logger::warn("GLUP")
	    << "Could not save " << filename << std::endl;
        Logger::warn("GLUP")
	    << "SimpleApplication::save() needs to be overloaded"
	    << std::endl;
        return false;
    }
    
    bool SimpleApplication::load(const std::string& filename) {
        Logger::warn("GLUP")
	    << "Could not load " << filename << std::endl;
        Logger::warn("GLUP")
	    << "SimpleApplication::load() needs to be overloaded"
	    << std::endl;
        return false;
    }
    
    bool SimpleApplication::can_load(const std::string& filename) {
        std::string extensions_str = supported_read_file_extensions();
        if(extensions_str == "") {
            return false;
        }
        if(extensions_str == "*") {
            return true;
        }
        std::string extension = FileSystem::extension(filename);
        std::vector<std::string> extensions;
        String::split_string(extensions_str, ';', extensions);
        for(index_t i=0; i<extensions.size(); ++i) {
            if(extensions[i] == extension) {
                return true;
            }
        }
        return false;
    }
    
    std::string SimpleApplication::supported_read_file_extensions() {
        return "";
    }

    std::string SimpleApplication::supported_write_file_extensions() {
        return "";
    }

    ImTextureID SimpleApplication::convert_to_ImTextureID(
	GLuint gl_texture_id_in
    ) {
        // It is not correct to directly cast a GLuint into a void*
        // (generates warnings), therefore I'm using a union.
        union {
            GLuint gl_texture_id;
            ImTextureID imgui_texture_id;
        };
        imgui_texture_id = nullptr;
        gl_texture_id = gl_texture_id_in;
        return imgui_texture_id;
    }
    
    void SimpleApplication::GL_initialize() {
	Application::GL_initialize();
        glGenTextures(1, &geogram_logo_texture_);
        glActiveTexture(GL_TEXTURE0 + GLUP_TEXTURE_2D_UNIT);
        glBindTexture(GL_TEXTURE_2D, geogram_logo_texture_);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexImage2Dxpm(geogram_logo_256_xpm);
	init_colormaps();
	std::string keys = CmdLine::get_arg("gfx:keypress");
	for(index_t i=0; i<index_t(keys.length()); ++i) {
	    char_callback((unsigned int)(keys[i]));
	}
    }
    
    void SimpleApplication::GL_terminate() {
        for(index_t i=0; i<colormaps_.size(); ++i) {
            if(colormaps_[i].texture != 0) {
                glDeleteTextures(1, &colormaps_[i].texture);
		colormaps_[i].texture = 0;
            }
        }
        if(geogram_logo_texture_ != 0) {
            glDeleteTextures(1, &geogram_logo_texture_);
	    geogram_logo_texture_ = 0;
        }
	Application::GL_terminate();
    }


    void SimpleApplication::browse(const std::string& path, bool subdirs) {
        std::vector<std::string> files;
        FileSystem::get_directory_entries(path,files);
        
        for(index_t i=0; i<files.size(); ++i) {
            if(FileSystem::is_directory(files[i]) && subdirs) {
                if(ImGui::BeginMenu(path_to_label(path_,files[i]).c_str())) {
                    browse(files[i]);
                    ImGui::EndMenu();
                }
            } else {
                if(can_load(files[i])) {
                    if(ImGui::MenuItem(path_to_label(path_,files[i]).c_str())) {
                        load(files[i]);
                    }
                }
            }
        }
    }
    
    void SimpleApplication::geogram_initialize(int argc, char** argv) {
	GEO::Application::geogram_initialize(argc, argv);
        if(filenames().size() == 1 &&
	   FileSystem::is_directory(filenames()[0])
	) {
            path_ = filenames()[0];
        } else if(filenames().size() > 0) {
            for(index_t i=0; i<filenames().size(); ++i) {
                load(filenames()[i]);
            }
            if(filenames().size() > 0) {
                path_ = FileSystem::dir_name(filenames()[filenames().size()-1]);
            }
        } else {
            path_ = FileSystem::documents_directory();
        }
        Logger::instance()->register_client(console_);
        Progress::set_client(status_bar_);
	set_default_layout();
    }

    void SimpleApplication::init_colormap(
        const std::string& name, const char** xpm_data
    ) {
        colormaps_.push_back(ColormapInfo());
        colormaps_.rbegin()->name = name;
        glGenTextures(1, &colormaps_.rbegin()->texture);
        glBindTexture(GL_TEXTURE_2D, colormaps_.rbegin()->texture);
        glTexImage2Dxpm(xpm_data);
        glGenerateMipmap(GL_TEXTURE_2D);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(
            GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR
        );
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glBindTexture(GL_TEXTURE_2D, 0);
    }

    void SimpleApplication::init_colormaps() {
        init_colormap("french", french_xpm);
        init_colormap("black_white", black_white_xpm);
        init_colormap("viridis", viridis_xpm);
        init_colormap("rainbow", rainbow_xpm);
        init_colormap("cei_60757", cei_60757_xpm);
        init_colormap("inferno", inferno_xpm);
        init_colormap("magma", magma_xpm);
        init_colormap("parula", parula_xpm);
        init_colormap("plasma", plasma_xpm);
        init_colormap("blue_red", blue_red_xpm);
    }

    bool SimpleApplication::exec_command(const char* command) {
#ifdef GEOGRAM_WITH_LUA	
	if(luaL_dostring(lua_state_,command)) {
	    adjust_lua_glup_state(lua_state_);
	    const char* msg = lua_tostring(lua_state_,-1);
	    const char* msg2 = strchr(msg,']');
	    if(msg2 != nullptr) {
		msg = msg2+2;
	    }
	    Logger::err("LUA") << "line " << msg << std::endl;
	    lua_error_occured_ = true;
	} else {
	    lua_error_occured_ = false;
	}
	return !lua_error_occured_;
#else
	geo_argused(command);
	Logger::err("LUA") << "Compiled without LUA support"
			   << std::endl;
	return false;
#endif	
    }

    void SimpleApplication::ImGui_initialize() {
	CmdLine::set_arg("gui:style","Light");
	Application::ImGui_initialize();
// Tooltips do not play well with touch screens.	
#ifdef GEO_OS_ANDROID	
	ImGui::DisableTooltips();
#endif	
    }
}
