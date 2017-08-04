//import <GLUP/current_profile/fragment_shader_preamble.h>
//import <GLUPES/fragment_shader_state.h>
//import <GLUP/stdglup.h>
//import <GLUP/current_profile/toggles.h>
//import <GLUP/current_profile/primitive.h>
//import <GLUPES/fragment_shader_utils.h>

#ifdef GL_ES                                                  
   varying vec4 color;                                        
   varying vec4 tex_coord;                                    
   varying float clip_dist;                                   
   varying highp float primitive_id;
#define glup_FragColor gl_FragColor
#else                                                         
   in vec4 color;                                             
   in vec4 tex_coord;                                         
   in float clip_dist;                                        
   flat in highp int primitive_id;
   out vec4 glup_FragColor;
#endif                                                        


void main() {
    
    if(glupIsEnabled(GLUP_CLIPPING) && (clip_dist < 0.0)) {              
        discard;                                                
    }                                                          

    if(glupIsEnabled(GLUP_PICKING)) {
        glup_FragColor = glup_picking(int(primitive_id));        
        return;
    }

    vec4 result;
    if(glupIsEnabled(GLUP_VERTEX_COLORS)) {
        result = color;
    } else {
        result = GLUP.front_color;
    }
    if(glupIsEnabled(GLUP_TEXTURING)) {
        result = glup_texturing(result, tex_coord);
    }
    glup_FragColor = result;
}                                                             
