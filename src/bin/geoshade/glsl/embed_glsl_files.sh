#!/bin/sh

glsl_files="
  course/raytrace_step1.glsl
  course/raytrace_step2.glsl
  course/raytrace_step3.glsl
  course/raytrace_step4.glsl
  course/raytrace_step5.glsl
  course/raytrace_step6.glsl  
  course/raytrace_step7.glsl
ShaderToy/AlloyPlatedVoronoi.glsl
ShaderToy/AndromedaJewel.glsl
ShaderToy/Circuits.glsl
ShaderToy/ContouredLayers.glsl
ShaderToy/FractalLand.glsl
ShaderToy/GeodesicTiling.glsl
ShaderToy/Geomechanical.glsl
ShaderToy/HexFlow.glsl
ShaderToy/JellyTubes.glsl
ShaderToy/MengerTunnel.glsl
ShaderToy/QuadTreeTruchet.glsl
ShaderToy/rabbit.glsl
ShaderToy/RayMarchingPrimitives.glsl
ShaderToy/RoundedVoronoiEdges.glsl
ShaderToy/RounderVoronoi.glsl
ShaderToy/SiggraphLogo.glsl
ShaderToy/TentacleObject.glsl
ShaderToy/ThePopularShader.glsl
ShaderToy/TracedMinkwskiTube.glsl
ShaderToy/VoxelPacMan.glsl
ShaderToy/Voxels.glsl
"

cat <<EOF
/*
 * This file was automatically generated, do not edit.
 */

#include <string>

void register_embedded_glsl_file(
      const char* filename, const char* body
);
void register_embedded_glsl_files(void);
   
void register_embedded_glsl_files() {
EOF

for f in $glsl_files
do
    echo "     register_embedded_glsl_file(\"$f\","
    cat $f | sed -e 's|\\|\\\\|' \
                 -e 's|"|\\\"|g' \
	         -e 's|^|        \"|' \
	         -e 's| *$|\\n\"|' 
    echo "     );"
    echo
done

echo "}"
