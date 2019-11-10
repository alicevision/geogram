const float FARAWAY=1e30;
const float EPSILON=1e-6;

struct Camera {
    vec3 Obs;
    vec3 View;
    vec3 Up;
    vec3 Horiz;
    float H;
    float W;
    float z;
};

struct Ray {
    vec3 Origin;
    vec3 Dir;
};

Camera camera(in vec3 Obs, in vec3 LookAt, in float aperture) {
   Camera C;
   C.Obs = Obs;
   C.View = normalize(LookAt - Obs);
   C.Horiz = normalize(cross(vec3(0.0, 0.0, 1.0), C.View));
   C.Up = cross(C.View, C.Horiz);
   C.W = float(iResolution.x);
   C.H = float(iResolution.y);
   C.z = (C.H/2.0) / tan((aperture * 3.1415 / 180.0) / 2.0);
   return C;
}

Ray launch(in Camera C, in vec2 XY) {
   return Ray(
      C.Obs,
      C.z*C.View+(XY.x-C.W/2.0)*C.Horiz+(XY.y-C.H/2.0)*C.Up 
   );
}

 
struct Sphere {
   vec3 Center;
   float R;
};

bool intersect_sphere(in Ray R, in Sphere S, out float t) {
   vec3 CO = R.Origin - S.Center;
   float a = dot(R.Dir, R.Dir);
   float b = 2.0*dot(R.Dir, CO);
   float c = dot(CO, CO) - S.R*S.R;
   float delta = b*b - 4.0*a*c;
   if(delta < 0.0) {
      return false;
   }
   t = (-b-sqrt(delta)) / (2.0*a);
   return true;
}
 
void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
 
   Camera C = camera(
       vec3(2.0, 2.0, 1.5),
       vec3(0.5, 0.5, 0.5),
       50.0       
   );
   Ray R = launch(C, fragCoord);
   Sphere S = Sphere(vec3(0.0, 0.0, 0.0), 0.5);
   
   fragColor = vec4(0.5, 0.5, 1.0, 1.0);
   
   float t;
   if(intersect_sphere(R,S,t)) {
      fragColor = vec4(1.0, 1.0, 1.0, 1.0);
   }

}




