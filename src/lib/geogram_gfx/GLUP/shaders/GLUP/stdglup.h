//import <GLUP/constants.h>

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

