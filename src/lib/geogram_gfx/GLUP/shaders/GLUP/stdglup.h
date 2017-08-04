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

