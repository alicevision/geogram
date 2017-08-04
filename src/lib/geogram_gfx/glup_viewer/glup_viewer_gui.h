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


#ifndef GLUP_VIEWER_GUI
#define GLUP_VIEWER_GUI

#ifdef __cplusplus
extern "C" {
#endif

    /**
     * \file geogram_gfx/glup_viewer/glup_viewer_gui.h
     * \short These functions are used internally by
     *  glup_viewer for interfacing with ImGUI.
     */
    
    struct GLFWwindow;
    
    void glup_viewer_gui_init(GLFWwindow* w);
    void glup_viewer_gui_cleanup();
    void glup_viewer_gui_begin_frame();
    void glup_viewer_gui_end_frame();
    int  glup_viewer_gui_takes_input();

    void glup_viewer_gui_mouse_button_callback(
        GLFWwindow* window, int button, int action, int mods
    );
    
    void glup_viewer_gui_scroll_callback(
        GLFWwindow* window, double xoffset, double yoffset
    );

    void glup_viewer_gui_key_callback(
        GLFWwindow* window, int key, int scancode, int action, int mods
    );

    void glup_viewer_gui_char_callback(GLFWwindow* window, unsigned int c);

    void glup_viewer_gui_resize(int width, int height);

    void glup_viewer_gui_draw_console(void);
    
#ifdef __cplusplus
}
#endif

#endif
