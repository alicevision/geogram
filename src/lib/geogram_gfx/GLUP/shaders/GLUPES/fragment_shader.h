//import <GLUP/current_profile/fragment_shader_preamble.h>
//import <GLUPES/fragment_shader_state.h>
//import <GLUP/stdglup.h>
//import <GLUP/current_profile/toggles.h>
//import <GLUP/current_profile/primitive.h>
//import <GLUPES/fragment_shader_utils.h>

#ifdef GL_ES
   varying vec3 vertex_clip_space;                            
   varying float clip_dist;                                   
   varying vec4 color;                                        
   varying vec4 tex_coord;                                    
   varying vec4 mesh_tex_coord;                               
   varying highp float primitive_id;                          
#define glup_FragColor gl_FragColor
#else
   in vec3 vertex_clip_space;                                 
   in float clip_dist;                                 
   in vec4 color;                                     
   in vec4 tex_coord;                                 
   in vec4 mesh_tex_coord;                            
   flat in highp int primitive_id;
   out vec4 glup_FragColor;
#endif                                                        


void main() {

    if(glupIsEnabled(GLUP_CLIPPING)) {
        if(glup_primitive_dimension == 2) {
            if(clip_dist < 0.0) {                                       
                discard;                               
            }                                         
        } else if(glup_primitive_dimension == 3) {
            if(
                clip_dist < 0.0 &&
                GLUP.clipping_mode == GLUP_CLIP_STANDARD
            ) {
                discard;
            }
        }
    }

    vec3 N;
    if(glupIsEnabled(GLUP_LIGHTING)) {
        vec3 U = dFdx(vertex_clip_space);                     
        vec3 V = dFdy(vertex_clip_space);                 
        N = normalize(cross(U,V));                       
    }
    
    glup_FragColor = glup_shading(
        color, tex_coord, N, int(primitive_id), mesh_tex_coord
    );
}                                                             
