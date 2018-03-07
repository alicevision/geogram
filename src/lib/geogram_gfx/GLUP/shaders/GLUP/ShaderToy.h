//import <GLUP/portable_fragment_shader.h>

uniform float iTime;
const vec2 iResolution=vec2(1024.0, 1024.0);
const vec4 iMouse=vec4(0.0, 0.0, 0.0, 0.0);
uniform sampler2D iChannel0;
uniform sampler2D iChannel1;
uniform sampler2D iChannel2;
uniform sampler2D iChannel3;

void mainImage(out vec4 fragColor, in vec2 fragCoord);

void main(void) {
     mainImage(glup_FragColor, 1024.0*tex_coord.xy);
}

