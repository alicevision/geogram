const float FARAWAY=1e30;
const float EPSILON=1e-3;

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

struct Material {
    vec3 Kd; // diffuse color
    vec3 Ke; // emissive color
};

Material diffuse(in vec3 Kd) {
   return Material(Kd, vec3(0.0, 0.0, 0.0));
}

Material light(in vec3 Ke) {
   return Material(vec3(0.0, 0.0, 0.0), Ke);
}

struct Object {
   Sphere sphere;
   Material material;
};

Object scene[3];

void init_scene() {
   float beta = float(iFrame)/30.0;
   float s = sin(beta);
   float c = cos(beta); 

   scene[0] = Object(
      Sphere(vec3(0.0, 0.0, 0.0),0.5), 
      diffuse(vec3(1.0, 1.0, 1.0))
   );

   scene[1] = Object(
      Sphere(vec3(0.7*s, 0.7*c, 0.0),0.1), 
      diffuse(vec3(1.0, 0.0, 0.0))
   );

   scene[2] = Object(
      Sphere(vec3(5.0, 0.0, 3.0),0.02),
      light(vec3(1.0, 1.0, 1.0)) 
   );

//   scene[3] = Object(
//      Sphere(vec3(1.5, 0.0, 1.5),0.02),
//      light(vec3(1.0, 1.0, 1.0)) 
//   );


}

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
 
bool shadow(in Ray R) {
   for(int i=0; i<scene.length(); ++i) {
        float t;
        if(
          scene[i].material.Ke == vec3(0.0, 0.0, 0.0) &&
          intersect_sphere(R, scene[i].sphere, t) &&
          t > EPSILON && t < 1.0
        ) {
          return true;
        }
    }
    return false;
}

vec3 lighting(in vec3 P, in vec3 N, in Material material) {
   if(material.Ke != vec3(0.0, 0.0, 0.0)) {
      return material.Ke;
   }  

   vec3 result = vec3(0.0, 0.0, 0.0);

   for(int i=0; i<scene.length(); ++i) {
      if(scene[i].material.Ke != vec3(0.0, 0.0, 0.0)) {
         Ray R2 = Ray(P, scene[i].sphere.Center);
         if(!shadow(R2)) {
           vec3 E = scene[i].sphere.Center - P;
           float lamb = max(0.0, dot(E,N) / length(E));
           result += lamb * material.Kd * scene[i].material.Ke;
         }
      }
   }

   return result;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {

   init_scene();

   Camera C = camera(
       vec3(2.0, 2.0, 1.5),
       vec3(0.5, 0.5, 0.5),
       50.0       
   );
   Ray R = launch(C, fragCoord);
   
   
   fragColor = vec4(0.5, 0.5, 1.0, 1.0);

   const float FARAWAY=1e30; 
   float t = FARAWAY;

   vec3 P;  // Point courant
   vec3 N;  // Normale
   Material material; // Couleur

   for(int i=0; i<scene.length(); ++i) {
       float cur_t;
       if(
          intersect_sphere(R, scene[i].sphere, cur_t) 
          && cur_t < t
       ) {
           t = cur_t;
           P = R.Origin + t*R.Dir;
           N = normalize(P - scene[i].sphere.Center);
           material = scene[i].material;
       } 
   }

   if(t != FARAWAY) {
      fragColor.rgb = lighting(P,N,material);
   }

}






