struct VSUniformState {              
    mat3  normal_matrix;             
    mat4 modelviewprojection_matrix; 
    mat4 modelview_matrix;           
    mat4 texture_matrix;
    vec4  world_clip_plane;          
    float point_size;                
};

uniform VSUniformState GLUP_VS;     

