//import <GLUP/current_profile/vertex_shader_preamble.h>
//import <GLUP/stdglup.h>
//import <GLUP/current_profile/toggles.h>
//import <GLUP/current_profile/primitive.h>
//import <GLUPES/vertex_shader_state.h>

#ifdef GL_ES                                          
   attribute vec4 vertex_in;                          
   attribute vec4 color_in;                           
   attribute vec4 tex_coord_in;                       
   attribute float vertex_id_in;                      
   varying vec3 vertex_clip_space;                    
   varying float clip_dist;                                   
   varying vec4 color;                                
   varying vec4 tex_coord;                            
   varying vec4 mesh_tex_coord;                       
   varying highp float primitive_id;                  
#else                                                 
   in vec4 vertex_in;                                 
   in vec4 color_in;                                  
   in vec4 tex_coord_in;                              
   in highp float vertex_id_in;                       
   out vec3 vertex_clip_space;                        
   out float clip_dist;                               
   out vec4 color;                                    
   out vec4 tex_coord;                                
   out vec4 mesh_tex_coord;                           
   flat out highp int primitive_id;                   
#endif                                                

        
void main() {                                         
    
    if(glupIsEnabled(GLUP_CLIPPING)) {                     
        clip_dist = dot(                                
            vertex_in, GLUP_VS.world_clip_plane         
        );                                              
    }                                                  

    if(glupIsEnabled(GLUP_PICKING)) {                        
#ifdef GL_ES                                            
        // Note: we need to add 0.5, else there are some precision
        // issues, and the integer mod() operation creates random
        // values...
        primitive_id = float(
            int(vertex_id_in+0.5)/glup_primitive_nb_vertices
        )+0.5;
#else                                                   
        primitive_id = int(vertex_id_in + 0.5)/glup_primitive_nb_vertices; 
#endif                                                  
    }                                                    
        
    if(glupIsEnabled(GLUP_LIGHTING)) {                             
        vertex_clip_space = (GLUP_VS.modelview_matrix * vertex_in).xyz;  
    }                                                  
    
    if(glupIsEnabled(GLUP_VERTEX_COLORS)) {                
        color = color_in;                               
    }                                                  
        
    if(glupIsEnabled(GLUP_TEXTURING)) {                    
        tex_coord = GLUP_VS.texture_matrix * tex_coord_in; 
    }                                                  

    if(glupIsEnabled(GLUP_DRAW_MESH)) {
        mesh_tex_coord = get_mesh_tex_coord(int(vertex_id_in + 0.5));
    }                                                  
    
    gl_Position = GLUP_VS.modelviewprojection_matrix * vertex_in;  
}                                                     
