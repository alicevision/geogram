/*
 *    _____   _       _   _     ____   
 *   /  ___| | |     | | | |   /  _ \  
 *   | |     | |     | | | |   | |_\ \ 
 *   | |  _  | |     | | | |   |  __ / 
 *   | |_| | | |___  | |_| |   | |     
 *   \_____/ |_____| \_____/   |_|     
 *
 *    _     _   _   _____   _          __  _____   _____
 *   | |   / / | | | ____| | |        / / | ____| |  _  \
 *   | |  / /  | | | |__   | |  __   / /  | |__   | |_| |
 *   | | / /   | | |  __|  | | /  | / /   |  __|  |  _  /
 *   | |/ /    | | | |___  | |/   |/ /    | |___  | | \ \
 *   |___/     |_| |_____| |___/|___/     |_____| |_|  \_\
 *
 *  Version 1.0
 *  Bruno Levy, April 2016
 *  INRIA, Project ALICE
 *
 *  Used internally by GLUP Viewer for interfacing with ImGUI
 *
 */

#include <geogram_gfx/basic/GL.h>
#include <geogram_gfx/glup_viewer/glup_viewer_gui_private.h>
#include <geogram_gfx/glup_viewer/glup_viewer_gui.h>
#include <geogram_gfx/glup_viewer/glup_viewer.h>
#include <geogram_gfx/third_party/ImGui/imgui.h>
#include <geogram_gfx/third_party/ImGui/imgui_impl_glfw_gl3.h>
#include <geogram_gfx/third_party/ImGui/imgui_impl_glfw.h>
#include <geogram_gfx/third_party/quicktext/glQuickText.h>

#include <geogram/basic/logger.h>
#include <geogram/basic/progress.h>
#include <geogram/basic/command_line.h>

/* 
 *Lots of documentation tags in GLFW that are
 * not understood by CLANG.
 */
#ifdef __clang__
#pragma GCC diagnostic ignored "-Wdocumentation"
#endif

#ifdef GEO_OS_EMSCRIPTEN
#include <GLFW/glfw3.h>
#include <emscripten.h>
#include <geogram/basic/file_system.h>
#else
#include <third_party/glfw/include/GLFW/glfw3.h>
#endif

#include <geogram_gfx/GLUP/GLUP.h>

#include <string.h>
#include <iostream>

#ifdef GEO_GL_LEGACY
static bool vanillaGL = false;
#endif

/***************************************************************************/

// Commands may need to update the GUI (when using
// the console or the progressbar). The Command class
// does that outside of the ImGui handler, but if client
// code does that directly, then we may have two nested
// invocations of the ImGui handler, which is not correct.
// This variable avoids to have two nested invocations of
// the ImGui handler.
static bool glup_viewer_gui_locked = false;

extern "C" {
    void glup_viewer_one_frame(void);
}

void glup_viewer_gui_update() {
    // It's a pity, but under Emscripten, only the browser can have the
    // control of the rendering loop, therefore it is not possible (or I
    // don't know how) to update the graphics during computations.
#ifndef GEO_OS_EMSCRIPTEN
    glup_viewer_post_redisplay();    
    if(!glup_viewer_gui_locked) {
        glup_viewer_one_frame();
    }
#endif    
}

/***************************************************************************/

void glup_viewer_gui_init(GLFWwindow* w) {
#ifdef GEO_GL_LEGACY        
    vanillaGL = (strcmp(glupCurrentProfileName(), "VanillaGL") == 0);
    if(vanillaGL) {
        GEO::Logger::out("ImGUI") << "Viewer GUI init (Vanilla)"
                                  << std::endl;        
        ImGui_ImplGlfw_Init(w, false);        
    } else
#endif
    {
        GEO::Logger::out("ImGUI") << "Viewer GUI init (GL3)"
                                  << std::endl;        
        ImGui_ImplGlfwGL3_Init(w, false);        
    }

    ImGuiStyle& style = ImGui::GetStyle();
    style.WindowRounding = 10.0f;
    style.FrameRounding = 10.0f;
    style.GrabRounding = 10.0f;
    ImGuiIO& io = ImGui::GetIO();
    io.IniFilename = NULL;
}

void glup_viewer_gui_cleanup() {
#ifdef GEO_GL_LEGACY        
    if(vanillaGL) {    
        ImGui_ImplGlfw_Shutdown();
    } else
#endif
    {
        ImGui_ImplGlfwGL3_Shutdown();        
    }
}


void glup_viewer_gui_begin_frame() {
    glup_viewer_gui_locked = true;
#ifdef GEO_GL_LEGACY            
    if(vanillaGL) {
        ImGui_ImplGlfw_NewFrame();        
    } else
#endif        
    {
        ImGui_ImplGlfwGL3_NewFrame();
    }
}

void glup_viewer_gui_end_frame() {
    GlupViewerDisplayFunc overlay_func = glup_viewer_get_overlay_func();
    if(overlay_func != nil) {
        overlay_func();
        ImGui::Render();
    }
    glup_viewer_gui_locked = false;
    // We flush the queued command here, once ImGui::Render() was
    // called, so that if it triggers a frame rendering (e.g. through
    // the Logger or the ProgressLogger), then ImGui calls are correctly
    // nested.
    GEO::Command::flush_queue();
    if(glup_viewer_gui_takes_input()) {
        glup_viewer_post_redisplay();
    }
}

int glup_viewer_gui_takes_input() {
    if(!glup_viewer_is_enabled(GLUP_VIEWER_TWEAKBARS)) {
        return 0;
    }
    return (
        ImGui::GetIO().WantCaptureMouse ||
        ImGui::GetIO().WantCaptureKeyboard
    ) ? 1 : 0;
}

void glup_viewer_gui_mouse_button_callback(
    GLFWwindow* window, int button, int action, int mods
) {
    glup_viewer_post_redisplay();
#ifdef GEO_GL_LEGACY            
    if(vanillaGL) {
        ImGui_ImplGlfw_MouseButtonCallback(window, button, action, mods);
    } else
#endif        
    {
        ImGui_ImplGlfwGL3_MouseButtonCallback(window, button, action, mods);
    }
}

void glup_viewer_gui_scroll_callback(
    GLFWwindow* window, double xoffset, double yoffset
) {
    glup_viewer_post_redisplay();    
#ifdef GEO_GL_LEGACY            
    if(vanillaGL) {
        ImGui_ImplGlfw_ScrollCallback(window, xoffset, yoffset);
    } else
#endif        
    {
        ImGui_ImplGlfwGL3_ScrollCallback(window, xoffset, yoffset);
    }
}

void glup_viewer_gui_key_callback(
    GLFWwindow* window, int key, int scancode, int action, int mods
) {
    glup_viewer_post_redisplay();    
#ifdef GEO_GL_LEGACY            
    if(vanillaGL) {
        ImGui_ImplGlFw_KeyCallback(window, key, scancode, action, mods);
    } else
#endif        
    {
        ImGui_ImplGlfwGL3_KeyCallback(window, key, scancode, action, mods);
    }
}

void glup_viewer_gui_char_callback(GLFWwindow* window, unsigned int c) {
    glup_viewer_post_redisplay();    
#ifdef GEO_GL_LEGACY            
    if(vanillaGL) {
        ImGui_ImplGlfw_CharCallback(window, c);
    } else
#endif        
    {
        ImGui_ImplGlfwGL3_CharCallback(window, c);
    }
}

void glup_viewer_gui_resize(int width, int height) {
    glup_viewer_post_redisplay();    
    ImGui::GetIO().DisplaySize.x = float(width);
    ImGui::GetIO().DisplaySize.y = float(height);
}

GLboolean glup_viewer_get_arg_bool(const char* param) {
    return GEO::CmdLine::get_arg_bool(param) ? GL_TRUE : GL_FALSE;
}

GLboolean glup_viewer_test_arg_string(const char* param, const char* arg) {
    return (GEO::CmdLine::get_arg(param) == arg) ? GL_TRUE : GL_FALSE;
}

void glup_viewer_set_screen_size_from_args() {
    std::string geometry = GEO::CmdLine::get_arg("gfx:geometry");
    int w,h;
    sscanf(geometry.c_str(),"%dx%d",&w,&h);
    glup_viewer_set_screen_size(w,h);
}


#ifdef GEO_OS_EMSCRIPTEN

extern "C" {
    
    void drag_drop(GLFWwindow* w, int nb, const char** p);
    void file_system_changed_callback();

    /**
     * \brief This function is called by the HTML shell each
     *  time a file is loaded.
     * \details For now, it uses the fact that the latest loaded
     *  file appears first in the list (to be improved).
     */
    EMSCRIPTEN_KEEPALIVE void file_system_changed_callback() {
        std::vector<std::string> all_files;
        GEO::FileSystem::get_directory_entries("/",all_files);
        if(
            all_files.size() > 0 && 
            GEO::FileSystem::is_file(all_files[0])
        ) {
            const char* pname = all_files[0].c_str();
            drag_drop(nil, 1, &pname);
        }
    }
}
#endif
