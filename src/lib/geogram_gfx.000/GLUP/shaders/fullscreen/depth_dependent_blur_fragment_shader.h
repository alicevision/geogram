//import <fullscreen/current_profile/fragment_shader_preamble.h>
//import <GLUP/stdglup.h>

glup_in vec2 tex_coord;  

uniform sampler2D texture2D;
uniform sampler2D depth_texture;

uniform float blur_width;
uniform bool vertical;

const float PI = 3.14159265;
const float threshold = 0.005;

float width  = float(textureSize(texture2D,0).x);
float height = float(textureSize(texture2D,0).y);

// 1D Gaussian distribution, s is standard deviation
float gaussian(in float x, in float s) {
    return exp(-x * x / (2 * s * s)) / (s * sqrt(2 * PI));
}

float get_z_coeff(in vec2 pos) {
    float zCoef = texture(depth_texture, pos).r;
    zCoef = 3.0 * (zCoef - 0.1) ;
    return zCoef;
}

float get_z_dist(in vec2 center_pos, in vec2 other_pos) {
    return abs(get_z_coeff(center_pos) - get_z_coeff(other_pos));
}

void compute_blur() {
    int n = int(floor(3.0f * blur_width) - 1);
    float sum = 0;
    int i;
    
    vec2 cur_pix_coords;
    vec4 cur_pix_tex;
    vec4 final_pix_tex = vec4(0.0);
    
    // Calculate the sum of weights for the blur
    for (i = -n; i <= n; i++) {
        float x_offset, y_offset;
        if (vertical) {
            x_offset = 0;
            y_offset = i;
        } else {
            x_offset = i;
            y_offset = 0;
        }
        
        x_offset = x_offset / width;
        y_offset = y_offset / height;
        
        cur_pix_coords = vec2(x_offset, y_offset) + tex_coord;
        
        if(get_z_dist(tex_coord, cur_pix_coords) <= threshold) {
            float weight = gaussian(i, blur_width);
            sum += weight;
        }
    }
    
    // Calculate the blurred color
    for (i = -n; i <= n; i++) {
        float x_offset, y_offset;
        if (vertical) {
            x_offset = 0;
            y_offset = i;
        } else {
            x_offset = i;
            y_offset = 0;
        }
        
        x_offset = x_offset / width;
        y_offset = y_offset / height;
        
        cur_pix_coords = vec2(x_offset, y_offset) + tex_coord;
        
        if(get_z_dist(tex_coord, cur_pix_coords) <= threshold) {
            cur_pix_tex = texture(texture2D, cur_pix_coords);
            float weight = gaussian(i, blur_width) / sum;
            final_pix_tex += cur_pix_tex * weight;
        }
    }
    glup_FragColor.rgb = final_pix_tex.rgb;
}

void main() {
    compute_blur();
}

