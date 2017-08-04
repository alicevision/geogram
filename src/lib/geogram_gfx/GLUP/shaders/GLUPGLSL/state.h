layout(shared)                              
uniform GLUPStateBlock {                    
    
    bool vertex_colors_enabled;                      
    
    vec4  front_color;                       
    vec4  back_color;                        
    
    bool draw_mesh_enabled;                  
    vec4  mesh_color;                        
    float mesh_width;                        
                                                              
    bool lighting_enabled;                   
    vec3 light_vector;                       
    vec3 light_half_vector;                  
                                                              
    bool texturing_enabled;                  
    bool indirect_texturing_enabled;                      
    int  texture_mode;                       
    int  texture_type;                               
                                                              
    float cells_shrink;                      
                                                              
    bool picking_enabled;                             
    int   picking_mode;                      
    int   picking_id;                         
    int   base_picking_id;                    
                                                      
    bool clipping_enabled;                   
    int   clipping_mode;                     
    vec4  clip_plane;                        
    vec4  world_clip_plane;                  
                                         
    mat4 modelviewprojection_matrix;         
    mat4 modelview_matrix;
    mat4 projection_matrix;    
    mat3 normal_matrix;                      
    mat4 texture_matrix;
    
    float point_size;                        
} GLUP;                                     

// Note: the 1D colormap is stored in a 2D texture, because
// 1D textures are not supported by all OpenGL implementations
uniform sampler2D texture1Dsampler;         
uniform sampler2D texture2Dsampler;         
uniform sampler3D texture3Dsampler;         
