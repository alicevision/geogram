//import <GLUP/current_profile/fragment_shader_preamble.h>
//import <GLUPGLSL/state.h>
//import <GLUP/stdglup.h>
//import <GLUP/current_profile/toggles.h>
//import <GLUP/current_profile/primitive.h>
//import <GLUP/fragment_shader_utils.h>

in float gl_ClipDistance[];                                

in VertexData {                            
    vec4 color;                             
    vec4 tex_coord;                         
} FragmentIn;                              

void main() {

#ifdef GL_ES    
    if(glupIsEnabled(GLUP_CLIPPING) && (gl_ClipDistance[0] < 0.0)) {
        discard;                                                
    }                                                          
#endif
    
    if(glupIsEnabled(GLUP_PICKING)) {
        glup_FragColor = glup_picking(gl_PrimitiveID);        
        return;
    }
    
    vec4 result;
    if(glupIsEnabled(GLUP_VERTEX_COLORS)) {
        result = FragmentIn.color;
    } else {
        result = GLUP.front_color;
    }
    if(glupIsEnabled(GLUP_TEXTURING)) {
        result = glup_texturing(result, FragmentIn.tex_coord);
    }
    
    glup_FragColor = result;
    glup_alpha_discard();    
}                                                                  

