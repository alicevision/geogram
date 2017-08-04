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
#include <geogram_gfx/glup_viewer/glup_viewer_gui.h>
#include <geogram_gfx/glup_viewer/glup_viewer.h>
#include <geogram_gfx/third_party/ImGui/imgui.h>
#include <geogram_gfx/third_party/ImGui/imgui_impl_glfw_gl3.h>
#include <geogram_gfx/third_party/ImGui/imgui_impl_glfw.h>
#include <geogram_gfx/third_party/quicktext/glQuickText.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/progress.h>

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

/**
 * \brief ApplicationLogger, grabbed from ImGui example.
 * \details Example of code:
 * \code
 *  static ExampleAppLog my_log;
 *  my_log.AddLog("Hello %d world\n", 123);
 *  my_log.Draw("title");
 * \endcode
 */

struct ExampleAppLog {
    ImGuiTextBuffer Buf;
    ImGuiTextFilter Filter;
    /** \brief Index to lines offset */
    ImVector<int>   LineOffsets;   
    bool            ScrollToBottom;

    void Clear() {
        Buf.clear();
        LineOffsets.clear();
    }

    void AddLog(const char* fmt, ...) IM_PRINTFARGS(2) {
        int old_size = Buf.size();
        va_list args;
        va_start(args, fmt);
        Buf.appendv(fmt, args);
        va_end(args);
        for (int new_size = Buf.size(); old_size < new_size; old_size++) {
            if (Buf[old_size] == '\n') {
                LineOffsets.push_back(old_size);
            }
        }
        ScrollToBottom = true;
    }

    void Draw(const char* title, bool* p_open = NULL) {
        ImGui::Begin(
            title, p_open,
            ImGuiWindowFlags_NoResize |
            ImGuiWindowFlags_NoMove |
            ImGuiWindowFlags_NoCollapse 
        );
        if (ImGui::Button("Clear")) Clear();
        ImGui::SameLine();
        bool copy = ImGui::Button("Copy");
        ImGui::SameLine();
        Filter.Draw("Filter", -100.0f);
        ImGui::Separator();
        ImGui::BeginChild(
            "scrolling", ImVec2(0,0), false,
            ImGuiWindowFlags_HorizontalScrollbar
        );
        if (copy) {
            ImGui::LogToClipboard();
        }

        if (Filter.IsActive()) {
            const char* buf_begin = Buf.begin();
            const char* line = buf_begin;
            for (int line_no = 0; line != NULL; line_no++) {
                const char* line_end =
                    (line_no < LineOffsets.Size) ?
                    buf_begin + LineOffsets[line_no] : NULL;
                
                if (Filter.PassFilter(line, line_end)) {
                    ImGui::TextUnformatted(line, line_end);
                }
                line = line_end && line_end[1] ? line_end + 1 : NULL;
            }
        } else {
            ImGui::TextUnformatted(Buf.begin());
        }

        if (ScrollToBottom) {
            ImGui::SetScrollHere(1.0f);
        }
        ScrollToBottom = false;
        ImGui::EndChild();
        ImGui::End();
    }
};

static ExampleAppLog* logger = nil;

void glup_viewer_gui_draw_console() {
    logger->Draw("Console");
}

namespace {
    class GlupViewerLoggerClient : public GEO::LoggerClient {
    public:
        /**
         * \copydoc GEO::LoggerClient::div()
         */
        virtual void div(const std::string& value) {
            if(logger != nil) {
                logger->AddLog("========== %s", value.c_str());
            }
        }
        
        /**
         * \copydoc GEO::LoggerClient::out()
         */
        virtual void out(const std::string& value) {
            if(logger != nil) {
                logger->AddLog("    %s", value.c_str());
            }
        }

        /**
         * \copydoc GEO::LoggerClient::warn()
         */
        virtual void warn(const std::string& value) {
            if(logger != nil) {
                logger->AddLog("[W] %s", value.c_str());
            }
        }
        
        /**
         * \copydoc GEO::LoggerClient::err()
         */
        virtual void err(const std::string& value) {
            if(logger != nil) {
                logger->AddLog("[E] %s", value.c_str());
            }
        }
        
        /**
         * \copydoc GEO::LoggerClient::status()
         */
        virtual void status(const std::string& value) {
            if(logger != nil) {
                logger->AddLog("[status] %s", value.c_str());
            }
        }
        
        virtual ~GlupViewerLoggerClient() {
        }
    };
    
    class GlupViewerProgressClient : public GEO::ProgressClient {
    public:
        /**
         * \copydoc GEO::ProgressClient::begin()
         */
        virtual void begin() {
        }
        
        /**
         * \copydoc GEO::ProgressClient::progress()
         */
        virtual void progress(GEO::index_t step, GEO::index_t percent) {
            GEO::geo_argused(step);
            GEO::geo_argused(percent);
        }
        
        /**
         * \copydoc GEO::ProgressClient::end()
         */
        virtual void end(bool canceled) {
            GEO::geo_argused(canceled);
        }
        
        virtual ~GlupViewerProgressClient() {
        }
    };
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
    style.WindowRounding = 0.0f;
    style.FrameRounding = 10.0f;
    style.GrabRounding = 10.0f;
    
    logger = new ExampleAppLog;
    GlupViewerLoggerClient* logger_client = new GlupViewerLoggerClient;
    GEO::Logger::instance()->register_client(logger_client);
}

void glup_viewer_gui_cleanup() {
    delete logger;
    logger = nil;
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

void glup_viewer_gui_mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
#ifdef GEO_GL_LEGACY            
    if(vanillaGL) {
        ImGui_ImplGlfw_MouseButtonCallback(window, button, action, mods);
    } else
#endif        
    {
        ImGui_ImplGlfwGL3_MouseButtonCallback(window, button, action, mods);
    }
}

void glup_viewer_gui_scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
#ifdef GEO_GL_LEGACY            
    if(vanillaGL) {
        ImGui_ImplGlfw_ScrollCallback(window, xoffset, yoffset);
    } else
#endif        
    {
        ImGui_ImplGlfwGL3_ScrollCallback(window, xoffset, yoffset);
    }
}

void glup_viewer_gui_key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
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
    ImGui::GetIO().DisplaySize.x = float(width);
    ImGui::GetIO().DisplaySize.y = float(height);
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
