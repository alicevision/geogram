
#ifdef GL_ES
int glup_mod(in int x, in int y) {
    return x - (x/y)*y;
}
#else
int glup_mod(in int x, in int y) {
    return x % y;
}
#endif


// Definitions for GLUPES shaders
// These macro are used to have portable
// declarations with 
// OpenGLES, GLSL 1.3 and GLSL 1.5,

#ifdef GL_ES
#ifdef GLUP_VERTEX_SHADER
#define glup_in attribute
#define glup_out varying
#elif defined GLUP_FRAGMENT_SHADER
#define glup_in varying
#endif
#define glup_flat
#define glup_id highp float
#else
#define glup_in in
#define glup_out out
#define glup_flat flat
#define glup_id highp int
#endif

#ifdef GL_ES
vec4 glup_texture(in sampler2D samp, in vec2 uv) {
    return texture2D(samp, uv);                 
}                                               
#else
vec4 glup_texture(in sampler2D samp, in vec2 uv) {
    return texture(samp, uv);
}                                              
#endif

#ifdef GLUP_FRAGMENT_SHADER
#ifdef GL_ES        
#define glup_FragColor gl_FragColor
#ifdef GL_EXT_frag_depth
#define glup_FragDepth gl_FragDepthEXT
#else
float glup_FragDepth; // depth updates will be ignored
#endif
#else
out vec4 glup_FragColor;
#define glup_FragDepth gl_FragDepth 
#endif
#endif

