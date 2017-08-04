//import <GLUP/current_profile/vertex_shader_preamble.h>
//import <GLUPES/vertex_shader_state.h>
//import <GLUP/stdglup.h>
//import <GLUP/current_profile/toggles.h>
//import <GLUP/current_profile/primitive.h>

#ifdef GL_ES                                          
   attribute vec4 vertex_in;                          
   attribute vec4 color_in;                           
   attribute vec4 tex_coord_in;                       
   attribute highp float vertex_id_in;                
   varying vec4 color;                                
   varying vec4 tex_coord;                            
   varying float clip_dist;                           
   varying highp float primitive_id;                  
#else                                                 
   in vec4 vertex_in;                                 
   in vec4 color_in;                                  
   in vec4 tex_coord_in;                              
   in highp float vertex_id_in;                       
   out vec4 color;                                    
   out vec4 tex_coord;                                
   out float clip_dist;                               
   flat out highp int primitive_id;                           
#endif                                                
                                                              
void main() {                                         
    if(glupIsEnabled(GLUP_CLIPPING)) {                     
        clip_dist = dot(vertex_in, GLUP_VS.world_clip_plane);
    }                                                  
    if(glupIsEnabled(GLUP_PICKING)) {                      
#ifdef GL_ES                                          
        primitive_id = float(int(vertex_id_in + 0.5)) + 0.5;
#else                                                 
        primitive_id = int(vertex_id_in + 0.5);              
#endif                                                
    }
    if(glupIsEnabled(GLUP_VERTEX_COLORS)) {                
        color = color_in;                              
    }                                                  
    if(glupIsEnabled(GLUP_TEXTURING)) {                    
        tex_coord = GLUP_VS.texture_matrix * tex_coord_in; 
    }
    if(glup_primitive == GLUP_POINTS) {
        gl_PointSize = GLUP_VS.point_size;
    }
    gl_Position = GLUP_VS.modelviewprojection_matrix * vertex_in;
}


