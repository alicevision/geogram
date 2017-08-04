/*
 *  Copyright (c) 2012-2016, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     Bruno.Levy@inria.fr
 *     http://www.loria.fr/~levy
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 */

#include <geogram_gfx/GLUP/GLUP_context_GLSL.h>
#include <geogram_gfx/basic/GLSL.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/progress.h>
#include <geogram/basic/logger.h>

#ifdef GEO_GL_150

namespace GLUP {
    using namespace GEO;
    
    /***********************************************************************/
    /***********************************************************************/
    /**** GLUP implementation using GLSL 1.5                             ***/
    /***********************************************************************/
    /***********************************************************************/

    void Context_GLSL150::setup() {
        GL_ARB_conservative_depth_ = extension_is_supported(
            "GL_ARB_conservative_depth"
        );
        Logger::out("GLUP") << "GL_ARB_conservative_depth = "
                            << GL_ARB_conservative_depth_
                            << std::endl;
        Context::setup();
        marching_tet_.create_UBO();        
        marching_hex_.create_UBO();        
        marching_prism_.create_UBO();
        marching_pyramid_.create_UBO();
        marching_connector_.create_UBO();
    }
    
    const char* Context_GLSL150::profile_name() const {
        return "GLUP150";
    }

#ifdef GEO_OS_APPLE
    static const char* GLUP150_shader_source_header =
        "#version 150                               \n"
        ;
#else
    static const char* GLUP150_shader_source_header =
        "#version 150 core                          \n"
        ;
#endif
    
    static const char* GLUP150_vshader_in_out_declaration =
        "in vec4 vertex_in;                         \n"
        "in vec4 color_in;                          \n"
        "in vec4 tex_coord_in;                      \n"
        "out VertexData {                           \n"
        "   vec4 transformed;                       \n"
        "   vec4 color;                             \n"
        "   vec4 tex_coord;                         \n"
        "} VertexOut;                               \n";


    // The geometry shader gets the input geometry and
    // attributes (color, tex_coords) through the functions
    // vertex_in(), color_in() and tex_coord_in(),
    // and the predicate prim_is_discarded().
    // There can be other ways of implementing these
    // functions, see vertex_gather_mode and gather
    // tesselation evaluation shader.
    // The function prim_is_discarded() returns true if
    // the primitive should be ignored by the geometry
    // shader (this happens with the gather tesselation
    // evaluation shader).
    
    static const char* GLUP150_gshader_in_out_declaration =
        "in VertexData {                            \n"
        "    vec4 transformed;                      \n"
        "    vec4 color;                            \n"
        "    vec4 tex_coord;                        \n"
        "} VertexIn[];                              \n"
        "                                           \n"
        "vec4 vertex_in(in int i) {                 \n"
        "    return gl_in[i].gl_Position;           \n"
        "}                                          \n"
        "                                           \n"
        "vec4 transformed_in(in int i) {            \n"
        "    return VertexIn[i].transformed;        \n"
        "}                                          \n"
        "                                           \n"
        "vec4 color_in(in int i) {                  \n"
        "    return VertexIn[i].color;              \n"
        "}                                          \n"
        "                                           \n"
        "vec4 tex_coord_in(in int i) {              \n"
        "    return VertexIn[i].tex_coord;          \n"
        "}                                          \n"
        "                                           \n"        
        "out FragmentData {                         \n"
        "   vec4 color;                             \n"
        "   vec4 tex_coord;                         \n"
        "   flat float diffuse;                     \n"
        "   flat float specular;                    \n"
        "   vec2 bary;                              \n"
        "   flat vec3 edge_mask;                    \n"
        "} FragmentOut;                             \n"
        "                                           \n"
        "bool prim_is_discarded() {                 \n"
        "   return false;                           \n"
        "}                                          \n";

   
    static const char* GLUP150_fshader_in_out_declaration =
        "out vec4 frag_color ;                      \n"
        "in float gl_ClipDistance[];                \n"        
        "                                           \n"
        "in FragmentData {                          \n"
        "   vec4 color;                             \n"
        "   vec4 tex_coord;                         \n"
        "   flat float diffuse;                     \n"
        "   flat float specular;                    \n"
        "   vec2 bary;                              \n"
        "   flat vec3 edge_mask;                    \n"
        "} FragmentIn;                              \n"
        "                                           \n"
        ;

    /**
     * \brief Declaration of input/output for simple fragment shaders
     *   (for points and lines).
     * \details This one is used when the fragment shader is directly
     *  plugged to the vertex shader (i.e. without geometry shader).
     */
    static const char* GLUP150_simple_fshader_in_out_declaration =
        "out vec4 frag_color ;                      \n"
        "in float gl_ClipDistance[];                \n"                
        "                                           \n"
        "in VertexData {                            \n"
        "   vec4 transformed;                       \n"
        "   vec4 color;                             \n"
        "   vec4 tex_coord;                         \n"
        "} FragmentIn;                              \n"
        "                                           \n"
        ;

    // Note There is packUnorm4x8() and unpackUnorm4x8() that does what
    // we want, but it is only supported in GLSL 4.1...
    static const char* GLUP150_fshader_utils =
        "vec4 int_to_vec4(in int x) {                                        \n"
        "  return vec4(                                                      \n"
        "     float(x         & 255)/255.0,                                  \n"
        "     float((x >>  8) & 255)/255.0,                                  \n"
        "     float((x >> 16) & 255)/255.0,                                  \n"
        "     float((x >> 24) & 255)/255.0                                   \n"
        "  );                                                                \n"
        "}                                                                   \n"
        "                                                                    \n"
        "void output_picking_id() {                                          \n"
        "   if(GLUP.picking_mode == GLUP_PICK_PRIMITIVE) {                   \n"
        "      frag_color = int_to_vec4(gl_PrimitiveID+GLUP.base_picking_id);\n"
        "   } else {                                                         \n"
        "      frag_color = int_to_vec4(GLUP.picking_id);                    \n"
        "   }                                                                \n"
        "}                                                                   \n"
        "                                                                    \n"
        "void get_color() {                                                  \n"
        "   if(vertex_colors_enabled()) {                                    \n"
        "        frag_color = FragmentIn.color;                              \n"
        "   } else {                                                         \n"
        "        frag_color = gl_FrontFacing ?                               \n"
        "                       GLUP.front_color : GLUP.back_color;          \n"
        "   }                                                                \n"
        "   if(texturing_enabled()) {                                        \n"
        "       vec4 tex_color;                                              \n"
        "       if(GLUP.texture_type == GLUP_TEXTURE_1D) {                   \n"
        "           tex_color = texture(                                     \n"
        "               texture1Dsampler, FragmentIn.tex_coord.xy            \n"
        "           );                                                       \n"
        "       } else if(GLUP.texture_type == GLUP_TEXTURE_2D) {            \n"
        "           tex_color = texture(                                     \n"
        "                texture2Dsampler, FragmentIn.tex_coord.xy           \n"
        "           );                                                       \n"
        "       } else if(GLUP.texture_type == GLUP_TEXTURE_3D) {            \n"
        "           tex_color = texture(                                     \n"
        "               texture3Dsampler, FragmentIn.tex_coord.xyz           \n"
        "           );                                                       \n"
        "       }                                                            \n"
        "       if(indirect_texturing_enabled()) {                           \n"
        "           tex_color = GLUP.texture_matrix * tex_color;             \n"
        "           tex_color = texture(texture1Dsampler, tex_color.xy);     \n"
        "       }                                                            \n"
        "       if(GLUP.texture_mode == GLUP_TEXTURE_REPLACE) {              \n"
        "             frag_color = tex_color;                                \n"
        "       } else if(GLUP.texture_mode == GLUP_TEXTURE_MODULATE) {      \n"
        "             frag_color *= tex_color;                               \n"
        "       } else if(GLUP.texture_mode == GLUP_TEXTURE_ADD) {           \n"
        "             frag_color += tex_color;                               \n"
        "       }                                                            \n"
        "   }                                                                \n"
        "}                                                                   \n"
        "                                                                    \n"
        "void output_lighting(in float diff, in float spec) {                \n"
        "   float s = gl_FrontFacing ? 1.0 : -1.0 ;                          \n"
        "   float sdiffuse = s * diff;                                       \n"
        "   if((glup_primitive_dimension == 3) &&                            \n"
        "      clipping_enabled() &&                                         \n"
        "      GLUP.clipping_mode == GLUP_CLIP_SLICE_CELLS) {                \n"
        "       sdiffuse = abs(sdiffuse);                                    \n"
        "   }                                                                \n"
        "   if(sdiffuse > 0.0) {                                             \n"
        "       vec3 vspec = spec*vec3(1.0,1.0,1.0);                         \n"
        "       frag_color = sdiffuse*frag_color + vec4(vspec,1.0);          \n"
        "       frag_color.rgb += vec3(0.2, 0.2, 0.2);                       \n"
        "   } else {                                                         \n"
        "       frag_color = vec4(0.2, 0.2, 0.2, 1.0);                       \n"
        "   }                                                                \n"
        "}                                                                   \n"
        "                                                                    \n"
        "void clip_fragment() {                                              \n"
        "   if(ES_profile &&                                                 \n"
        "      clipping_enabled() &&                                         \n"
        "      gl_ClipDistance[0] < 0.0                                      \n"
        "   ) {                                                              \n"
        "     discard;                                                       \n"
        "   }                                                                \n"
        "}                                                                   \n"
        ;
    
#define GLUP150_std(prim)                    \
        GLUP150_shader_source_header,        \
        profile_dependent_declarations(),    \
        uniform_state_declaration(),         \
        primitive_declaration(prim).c_str(), \
        toggles_declaration()

    // GLUP_POINTS ********************************************************

    static const char* GLUP150_points_and_lines_vshader_source =
        "void main() {                                              \n"
        "    if(clipping_enabled()) {                               \n"
        "       gl_ClipDistance[0] = dot(                           \n"
        "            vertex_in, GLUP.world_clip_plane               \n"
        "       );                                                  \n"
        "   } else {                                                \n"
        "      gl_ClipDistance[0] = 0.0;                            \n"
        "   }                                                       \n"
        "   if(vertex_colors_enabled()) {                           \n"
        "      VertexOut.color = color_in;                          \n"
        "   }                                                       \n"
        "   if(texturing_enabled()) {                               \n"
        "      if(indirect_texturing_enabled()) {                   \n"
        "          VertexOut.tex_coord = tex_coord_in;              \n"
        "      } else {                                             \n"
        "          VertexOut.tex_coord =                            \n"
        "                       GLUP.texture_matrix * tex_coord_in; \n"
        "      }                                                    \n"
        "   }                                                       \n"
        "   gl_PointSize = GLUP.point_size;                         \n"
        "   gl_Position =                                           \n"
        "                GLUP.modelviewprojection_matrix*vertex_in; \n"
        "}                                                          \n";
    
    // Note: depth update is not correct, it should be something like:
    // (to be checked...)
    // gl_FragDepth = gl_FragCoord.z +
    //   (pt_size*0.0001)/3.0 * gl_ProjectionMatrix[2].z * sqrt(1.0 - r2);

    static const char* GLUP150_points_fshader_source =
        "void main() {                                                      \n"
        "   clip_fragment();                                                \n"
        "   vec2 V = 2.0*(gl_PointCoord - vec2(0.5, 0.5));                  \n"
        "   float one_minus_r2 = 1.0 - dot(V,V);                            \n"
        "   if(one_minus_r2 < 0.0) {                                        \n"
        "      discard;                                                     \n"
        "   }                                                               \n"
        "   vec3 W = vec3(V.x, -V.y, sqrt(one_minus_r2));                   \n"
        "   update_depth(gl_FragCoord.z - 0.001 * W.z);                     \n"
        "   if(picking_enabled()) {                                         \n"
        "        output_picking_id();                                       \n"
        "   } else {                                                        \n"
        "        get_color();                                               \n"
        "        if(lighting_enabled()) {                                   \n"
        "            float diff = dot(W,GLUP.light_vector);                 \n"
        "            float spec = dot(W,GLUP.light_half_vector);            \n"
        "            spec = pow(spec,30.0);                                 \n"
        "            output_lighting(diff,spec);                            \n"
        "        }                                                          \n"
        "    }                                                              \n"
        "}                                                                  \n";

    void Context_GLSL150::setup_GLUP_POINTS() {

        if(!use_core_profile_) {
            glEnable(GL_POINT_SPRITE);
            // Not needed anymore it seems.
            // glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
        }
        
        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            GLUP150_std(GLUP_POINTS),
            GLUP150_vshader_in_out_declaration,
            GLUP150_points_and_lines_vshader_source,
            0
        );

        GLuint fshader = 0;

        if(GL_ARB_conservative_depth_) {
            fshader = GLSL::compile_shader(
                GL_FRAGMENT_SHADER,
                GLUP150_std(GLUP_POINTS),
                GLUP150_simple_fshader_in_out_declaration,
                GLUP150_fshader_utils,
                "#extension GL_ARB_conservative_depth : enable \n",
                "layout (depth_less) out float gl_FragDepth;   \n",
                "void update_depth(in float f) { gl_FragDepth = f; } \n",
                GLUP150_points_fshader_source,
                0
             );
        } else {
            fshader = GLSL::compile_shader(
                GL_FRAGMENT_SHADER,
                GLUP150_std(GLUP_POINTS),
                GLUP150_simple_fshader_in_out_declaration,
                GLUP150_fshader_utils,
                "void update_depth(in float f) { } \n",
                GLUP150_points_fshader_source,
                0
             );
        }

        GLuint program = GLSL::create_program_from_shaders_no_link(
            vshader, fshader, 0 
        );
        
        set_primitive_info(GLUP_POINTS, GL_POINTS, program);

        glDeleteShader(vshader);
        glDeleteShader(fshader);
    }

    // GLUP_LINES *******************************************************

    static const char* GLUP150_lines_fshader_source =
        "void main() {                                              \n"
        "   clip_fragment();                                        \n"
        "   if(picking_enabled()) {                                 \n"
        "      output_picking_id();                                 \n"
        "   } else {                                                \n"
        "      get_color();                                         \n"
        "   }                                                       \n"
        "}                                                          \n";
    
    void Context_GLSL150::setup_GLUP_LINES() {
        
        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            GLUP150_std(GLUP_LINES),
            GLUP150_vshader_in_out_declaration,
            GLUP150_points_and_lines_vshader_source,
            0
        );
        
        GLuint fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER,
            GLUP150_std(GLUP_LINES),
            GLUP150_simple_fshader_in_out_declaration,
            GLUP150_fshader_utils,                        
            GLUP150_lines_fshader_source,
            0
        );

        GLuint program = GLSL::create_program_from_shaders_no_link(
            vshader, fshader, 0 
        );
        
        set_primitive_info(GLUP_LINES, GL_LINES, program);

        glDeleteShader(vshader);
        glDeleteShader(fshader);
    }

    // GLUP_TRIANGLES *******************************************************

    /**
     * \brief The fragment shader for polygons if GLSL version is 1.50.
     */
    static const char* GLUP150_triangle_fshader_source = 
        "float edge_factor() {                                              \n"
        "    vec3 bary3 = vec3(                                             \n"
        "       FragmentIn.bary.x,                                          \n"
        "       FragmentIn.bary.y,                                          \n"
        "       1.0-FragmentIn.bary.x-FragmentIn.bary.y                     \n"
        "    ) ;                                                            \n"
        "    vec3 d = fwidth(bary3);                                        \n"
        "    vec3 a3 = smoothstep(                                          \n"
        "                  vec3(0.0,0.0,0.0), d*GLUP.mesh_width, bary3      \n"
        "    );                                                             \n"
        "    a3 = vec3(1.0, 1.0, 1.0)                                       \n"
        "           - FragmentIn.edge_mask + FragmentIn.edge_mask*a3;       \n"
        "    return min(min(a3.x, a3.y), a3.z);                             \n"
        "}                                                                  \n"
        "                                                                   \n"
        "void main() {                                                      \n"
        "    clip_fragment();                                               \n"
        "    if(picking_enabled()) {                                        \n"
        "        output_picking_id();                                       \n"
        "    } else {                                                       \n"
        "        get_color();                                               \n"
        "        if(lighting_enabled()) {                                   \n"
        "          output_lighting(FragmentIn.diffuse, FragmentIn.specular);\n"
        "        }                                                          \n"
        "        if(draw_mesh_enabled()) {                                  \n"
        "            frag_color = mix(                                      \n"
        "                  GLUP.mesh_color,frag_color,edge_factor()         \n"
        "            );                                                     \n"
        "        }                                                          \n"
        "    }                                                              \n"
        "}                                                                  \n";
    
    /** 
     * \brief Some utility functions for the geometry shaders.
     * \details Provides functions for clipping, projection, and
     *  for generating shaded polygons.
     *  - flat_shaded_triangle(p1,p2,p3,pp1,pp2,pp3,do_clip) where
     *   (p1,p2,p3) are the coordinates in world space, (pp1,pp2,pp3) the
     *   transformed coordinates in clip space and do_clip specifies whether
     *   the triangle should be clipped.
     *  - flat_shaded_quad(p1,p2,p3,p4,pp1,pp2,pp3,pp4,do_clip,edges)
     */
    static const char* GLUP150_gshader_utils_source =
        "out float gl_ClipDistance[];                                       \n"
        "vec4 projected[nb_vertices];                                       \n"
        "                                                                   \n"
        "void project_vertices() {                                          \n"
        "   for(int i=0; i<nb_vertices; ++i) {                              \n"
        "     if((glup_primitive_dimension == 3) && clipping_enabled()) {   \n"
        "        projected[i] =                                             \n"
        "                  GLUP.modelviewprojection_matrix * vertex_in(i);  \n"
        "     } else {                                                      \n"
        "        projected[i] = transformed_in(i);                          \n"
        "     }                                                             \n"
        "   }                                                               \n"
        "   if(GLUP.cells_shrink != 0.0) {                                  \n"
        "       vec4 g = vec4(0.0, 0.0, 0.0, 0.0);                          \n"
        "       for(int i=0; i<nb_vertices; ++i) {                          \n"
        "            g += projected[i];                                     \n"
        "       }                                                           \n"
        "       g /= float(nb_vertices);                                    \n"
        "       float s = GLUP.cells_shrink;                                \n"
        "       for(int i=0; i<nb_vertices; ++i) {                          \n"
        "            projected[i] = mix(projected[i], g, s);                \n"
        "       }                                                           \n"
        "   }                                                               \n"
        "}                                                                  \n"
        "                                                                   \n"
        "float clip(in vec4 V, in bool do_clip) {                           \n"
        "    return do_clip?dot(V,GLUP.world_clip_plane):1.0;               \n"
        "}                                                                  \n"
        "                                                                   \n"
        "bool cell_is_clipped() {                                           \n"
        "  if(prim_is_discarded()) {                                        \n"
        "     return true;                                                  \n"
        "  }                                                                \n"
        "  if(                                                              \n"
        "     (glup_primitive_dimension != 3) ||                            \n"
        "     !clipping_enabled() ||                                        \n"
        "     GLUP.clipping_mode==GLUP_CLIP_STANDARD  ||                    \n"
        "     GLUP.clipping_mode==GLUP_CLIP_SLICE_CELLS                     \n" 
        "  ) {                                                              \n"
        "     return false;                                                 \n"
        "  }                                                                \n"
        "  int count = 0;                                                   \n"
        "  for(int i=0; i<nb_vertices; ++i) {                               \n"
        "       count += int(clip(vertex_in(i),true) >= 0.0);               \n"
        "  }                                                                \n"
        "  if(GLUP.clipping_mode==GLUP_CLIP_WHOLE_CELLS && count == 0) {    \n"
        "    return true;                                                   \n"
        "  }                                                                \n"
        "  if(                                                              \n"
        "      GLUP.clipping_mode==GLUP_CLIP_STRADDLING_CELLS &&            \n"
        "      (count==0 || count==nb_vertices)                             \n"
        "  ) {                                                              \n"
        "    return true;                                                   \n"
        "  }                                                                \n"
        "  return false;                                                    \n"
        "}                                                                  \n"
        "                                                                   \n"
        " /* L is supposed to be normalized */                              \n"
        "float cosangle(in vec3 N, in vec3 L) {                             \n"
        "   float s = inversesqrt(dot(N,N)) ;                               \n"
        "   return s*dot(N,L) ;                                             \n"
        "}                                                                  \n"
        "                                                                   \n"
        "void compute_lighting(in vec3 N) {                                 \n"
        "     FragmentOut.diffuse = cosangle(N,GLUP.light_vector) ;         \n"
        "     FragmentOut.specular = abs(cosangle(N,GLUP.light_half_vector));\n"
        "     FragmentOut.specular = pow(FragmentOut.specular,30.0);        \n"
        "}                                                                  \n"
        "                                                                   \n"
        "void emit_vertex(in int i, in bool do_clip) {                      \n"
        "   gl_ClipDistance[0] = clip(vertex_in(i),do_clip);                \n"
        "   gl_Position = projected[i];                                     \n"
        "   if(vertex_colors_enabled()) {                                   \n"
        "      FragmentOut.color = color_in(i);                             \n"
        "   }                                                               \n"
        "   if(texturing_enabled()) {                                       \n"
        "      FragmentOut.tex_coord = tex_coord_in(i);                     \n"
        "   }                                                               \n"
        "   EmitVertex();                                                   \n"
        "}                                                                  \n"
        "                                                                   \n"
        "void flat_shaded_triangle(                                         \n"
        "     in int i1, in int i2, in int i3,                              \n"
        "     in bool do_clip                                               \n"
        ") {                                                                \n"
        "   if(lighting_enabled() && !picking_enabled()) {                  \n"
        "      vec4 p1 = vertex_in(i1);                                     \n"
        "      vec4 p2 = vertex_in(i2);                                     \n"
        "      vec4 p3 = vertex_in(i3);                                     \n"
        "      vec3 N = GLUP.normal_matrix *                                \n"
        "            cross((p2-p1).xyz,(p3-p1).xyz) ;                       \n"
        "      compute_lighting(N);                                         \n"
        "   }                                                               \n"
        "   FragmentOut.edge_mask = vec3(1.0,1.0,1.0);                      \n"
        "   FragmentOut.bary = vec2(0.0,0.0);                               \n"
        "   emit_vertex(i1,do_clip);                                        \n"
        "   FragmentOut.bary = vec2(1.0,0.0);                               \n"
        "   emit_vertex(i2,do_clip);                                        \n"
        "   FragmentOut.bary = vec2(0.0,1.0);                               \n"
        "   emit_vertex(i3,do_clip);                                        \n"
        "   EndPrimitive();                                                 \n"
        "}                                                                  \n"
        "                                                                   \n"
        "void flat_shaded_quad(                                             \n"
        "     in int i1, in int i2, in int i3, in int i4,                   \n"
        "     in bool do_clip                                               \n"
        "  ) {                                                              \n"
        "   if(lighting_enabled() && !picking_enabled()) {                  \n"
        "      vec4 p1 = vertex_in(i1);                                     \n"
        "      vec4 p2 = vertex_in(i2);                                     \n"
        "      vec4 p3 = vertex_in(i3);                                     \n"
        "      vec4 p4 = vertex_in(i4);                                     \n"
        "      vec3 N = GLUP.normal_matrix * (                              \n" 
        "           cross((p2-p1).xyz,(p4-p1).xyz) -                        \n"
        "           cross((p4-p3).xyz,(p2-p3).xyz)                          \n"
        "      );                                                           \n"
        "      compute_lighting(N);                                         \n"
        "   }                                                               \n"
        "   FragmentOut.edge_mask = vec3(0.0, 1.0, 1.0);                    \n"
        "   FragmentOut.bary=vec2(1.0,0.0);                                 \n"
        "   emit_vertex(i1,do_clip);                                        \n"
        "   FragmentOut.bary=vec2(0.0,1.0);                                 \n"
        "   emit_vertex(i2,do_clip);                                        \n"
        "   FragmentOut.bary=vec2(0.0,0.0);                                 \n"
        "   emit_vertex(i3,do_clip);                                        \n"
        "   FragmentOut.edge_mask = vec3(0.0, 1.0, 1.0);                    \n"
        "   FragmentOut.bary=vec2(1.0,0.0);                                 \n"
        "   emit_vertex(i4,do_clip);                                        \n"
        "   EndPrimitive();                                                 \n"
        "}                                                                  \n"
        ;

    /**
     * \brief The vertex shader.
     * \details Used by points, quads, tets, prisms
     */
    static const char* GLUP150_vshader_transform_source =
        " void main(void) {                                                 \n"
        "     if(vertex_colors_enabled()) {                                 \n"
        "        VertexOut.color = color_in;                                \n"
        "     }                                                             \n"
        "     if(texturing_enabled()) {                                     \n"
        "         if(indirect_texturing_enabled()) {                        \n"
        "             VertexOut.tex_coord = tex_coord_in;                   \n"
        "         } else {                                                  \n"
        "             VertexOut.tex_coord =                                 \n"
        "                          GLUP.texture_matrix * tex_coord_in;      \n"
        "         }                                                         \n"
        "     }                                                             \n"
        "     if(                                                           \n"
        "         glup_primitive_dimension != 3 ||                          \n"
        "         !clipping_enabled()                                       \n"
        "     ) {                                                           \n"
        "        VertexOut.transformed =                                    \n"
        "                      GLUP.modelviewprojection_matrix * vertex_in; \n"
        "     }                                                             \n"
        "     gl_Position = vertex_in;                                      \n"
        " }                                                                 \n";

    /**
     * \brief The geometry shader for triangles.
     * \details Uses vshader_transform and gshader_utils.
     */
    static const char* GLUP150_gshader_tri_source =
        "void main() {                                                      \n"
        "    gl_PrimitiveID = gl_PrimitiveIDIn;                             \n"
        "    project_vertices();                                            \n"
        "    flat_shaded_triangle(0,1,2,true);                              \n"
        "}                                                                  \n";
    
    void Context_GLSL150::setup_GLUP_TRIANGLES() {

        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            GLUP150_std(GLUP_TRIANGLES),
            GLUP150_vshader_in_out_declaration,
            GLUP150_vshader_transform_source,
            0
        );
        
        GLuint gshader = GLSL::compile_shader(
            GL_GEOMETRY_SHADER,
            GLUP150_std(GLUP_TRIANGLES),
            "layout(triangles) in;                         \n",
            "layout(triangle_strip, max_vertices = 3) out; \n",
            GLUP150_gshader_in_out_declaration,
            "const int nb_vertices = 3;",
            GLUP150_gshader_utils_source,
            GLUP150_gshader_tri_source,
            0
        );

        GLuint fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER,
            GLUP150_std(GLUP_TRIANGLES),
            GLUP150_fshader_in_out_declaration,
            GLUP150_fshader_utils,
            GLUP150_triangle_fshader_source,
            0
        );
                
        GLuint program = 
            GLSL::create_program_from_shaders_no_link(
                vshader, gshader, fshader, 0 
            );

        set_primitive_info(GLUP_TRIANGLES, GL_TRIANGLES, program);
        
        glDeleteShader(vshader);
        glDeleteShader(gshader);
        glDeleteShader(fshader);        
    }

    // GLUP_QUADS ***********************************************************

    /**
     * \brief The geometry shader for quads.
     * \details Uses vshader_transform and gshader_utils.
     */
    static const char* GLUP150_gshader_quad_source =
        "void main() {                                                      \n"
        "    gl_PrimitiveID = gl_PrimitiveIDIn;                             \n"
        "    project_vertices();                                            \n"
        "    flat_shaded_quad(0,1,3,2,true);                                \n"
        "}                                                                  \n";
    
    void Context_GLSL150::setup_GLUP_QUADS() {

        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            GLUP150_std(GLUP_QUADS),
            GLUP150_vshader_in_out_declaration,
            GLUP150_vshader_transform_source,
            0
        );

        
        GLuint gshader = GLSL::compile_shader(
            GL_GEOMETRY_SHADER,
            GLUP150_std(GLUP_QUADS),
            "layout(lines_adjacency) in;                   \n",
            "layout(triangle_strip, max_vertices = 4) out; \n",
            GLUP150_gshader_in_out_declaration,
            "const int nb_vertices = 4;",
            GLUP150_gshader_utils_source,
            GLUP150_gshader_quad_source,
            0
        );

        GLuint fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER,
            GLUP150_std(GLUP_QUADS),
            GLUP150_fshader_in_out_declaration,
            GLUP150_fshader_utils,                        
            GLUP150_triangle_fshader_source,
            0
        );

        GLuint program = 
            GLSL::create_program_from_shaders_no_link(
                vshader, gshader, fshader, 0
            );
        
        set_primitive_info(GLUP_QUADS, GL_LINES_ADJACENCY, program);
        
        glDeleteShader(vshader);
        glDeleteShader(gshader);
        glDeleteShader(fshader);        
    }

    // GLUP_TETRAHEDRA ******************************************************

    static const char* GLUP150_marching_cells_utils =
       "int compute_config() {                               \n"
       "  int result = 0;                                    \n"
       "  for(int v=0; v<cell_nb_vertices; ++v) {            \n"
       "    if(                                              \n"
       "      dot(vertex_in(v),GLUP.world_clip_plane) > 0.0  \n"
       "    ) {                                              \n"
       "       result = result | (1 << v);                   \n"
       "    }                                                \n"
       "  }                                                  \n" 
       "  return result;                                     \n" 
       "}                                                    \n"
       "void emit_isect_vertex(in int i) {                   \n"
       "   gl_ClipDistance[0] = 1.0;                         \n"
       "   gl_Position = isect_point[i];                     \n"
       "   if(vertex_colors_enabled()) {                     \n"
       "      FragmentOut.color = isect_color[i];            \n"
       "   }                                                 \n"
       "   if(texturing_enabled()) {                         \n"
       "      FragmentOut.tex_coord = isect_tex_coord[i];    \n"
       "   }                                                 \n"
       "   EmitVertex();                                     \n"
       "}                                                    \n"
       "void isect_triangle(int i, int j, int k) {           \n"
       "   FragmentOut.bary = vec2(0.0,0.0);                 \n"
       "   emit_isect_vertex(i);                             \n"
       "   FragmentOut.bary = vec2(0.0,1.0);                 \n"        
       "   emit_isect_vertex(j);                             \n"
       "   FragmentOut.bary = vec2(1.0,0.0);                 \n"                
       "   emit_isect_vertex(k);                             \n"
       "   EndPrimitive();                                   \n"
       "}                                                    \n"
       "void draw_marching_cell() {                          \n"
       "   int config = compute_config();                    \n"
       "   if(config_size(config) == 0) { return; }          \n"
       "   compute_intersections();                          \n"
       "   if(lighting_enabled() && !picking_enabled()) {    \n"
       "      compute_lighting(GLUP.world_clip_plane.xyz);   \n"
       "   }                                                 \n"
       "   int size = config_size(config);                   \n"
       "   if(draw_mesh_enabled()) {                         \n"
       "      if(size == 3) {                                \n"
       "        FragmentOut.edge_mask = vec3(1.0,1.0,1.0);   \n"
       "        isect_triangle(                              \n"
       "           config_edge(config,0),                    \n"
       "           config_edge(config,1),                    \n"
       "           config_edge(config,2)                     \n"
       "        );                                           \n"
       "      } else {                                       \n"
       "        FragmentOut.edge_mask = vec3(1.0,0.0,1.0);   \n"
       "        isect_triangle(                              \n"
       "           config_edge(config,0),                    \n"
       "           config_edge(config,1),                    \n"
       "           config_edge(config,2)                     \n"
       "        );                                           \n"
       "        FragmentOut.edge_mask = vec3(0.0,0.0,1.0);   \n"        
       "        for(int i=1; i+2<size; ++i) {                \n"
       "          isect_triangle(                            \n"
       "             config_edge(config,0),                  \n"
       "             config_edge(config,i),                  \n"
       "             config_edge(config,i+1)                 \n"
       "          );                                         \n"
       "        }                                            \n"
       "        FragmentOut.edge_mask = vec3(0.0,1.0,1.0);   \n"
       "        isect_triangle(                              \n"
       "           config_edge(config,0),                    \n"
       "           config_edge(config,size-2),               \n"
       "           config_edge(config,size-1)                \n"
       "        );                                           \n"
       "      }                                              \n" 
       "   } else {                                          \n"
       "      FragmentOut.edge_mask = vec3(1.0,1.0,1.0);     \n"        
       "      for(int i=1; i+1<size; ++i) {                  \n"
       "        isect_triangle(                              \n"
       "           config_edge(config,0),                    \n"
       "           config_edge(config,i),                    \n"
       "           config_edge(config,i+1)                   \n"
       "        );                                           \n"
       "      }                                              \n"
       "   }                                                 \n" 
       "}                                                    \n"
       ;
    
    /**
     * \brief The geometry shader for tetrahedra.
     * \details Uses v_shader_transform and gshader_utils.
     */
    static const char* GLUP150_gshader_tet_source =
        "void main() {                                                      \n"
        "    if(cell_is_clipped()) {                                        \n"
        "        return;                                                    \n"
        "    }                                                              \n"
        "    gl_PrimitiveID = gl_PrimitiveIDIn;                             \n"
        "    if(clipping_enabled() &&                                       \n"
        "        GLUP.clipping_mode == GLUP_CLIP_SLICE_CELLS) {             \n"
        "       draw_marching_cell();                                       \n"
        "       return;                                                     \n"
        "    }                                                              \n"
        "    project_vertices();                                            \n"
        "    bool do_clip = (clipping_enabled() &&                          \n"
        "                    GLUP.clipping_mode==GLUP_CLIP_STANDARD);       \n"
        "    flat_shaded_triangle(0,1,2,do_clip);                           \n"
        "    flat_shaded_triangle(1,0,3,do_clip);                           \n"
        "    flat_shaded_triangle(0,2,3,do_clip);                           \n"
        "    flat_shaded_triangle(2,1,3,do_clip);                           \n"
        "}                                                                  \n";
    
    
    void Context_GLSL150::setup_GLUP_TETRAHEDRA() {

        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            GLUP150_std(GLUP_TETRAHEDRA),
            GLUP150_vshader_in_out_declaration,
            GLUP150_vshader_transform_source,
            0
        );

        GLuint gshader = GLSL::compile_shader(
            GL_GEOMETRY_SHADER,
            GLUP150_std(GLUP_TETRAHEDRA),
            marching_tet_.GLSL_uniform_state_declaration(),
            "layout(lines_adjacency) in;                    \n",
            "layout(triangle_strip, max_vertices = 12) out; \n",
            GLUP150_gshader_in_out_declaration,
            "const int nb_vertices = 4; ",
            GLUP150_gshader_utils_source,
            marching_tet_.GLSL_compute_intersections(),
            GLUP150_marching_cells_utils,
            GLUP150_gshader_tet_source,
            0
        );

        GLuint fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER,
            GLUP150_std(GLUP_TETRAHEDRA),
            GLUP150_fshader_in_out_declaration,
            GLUP150_fshader_utils,                        
            GLUP150_triangle_fshader_source,
            0
        );

        GLuint program = 
            GLSL::create_program_from_shaders_no_link(
                vshader, gshader, fshader, 0 
            );

        set_primitive_info(GLUP_TETRAHEDRA, GL_LINES_ADJACENCY, program);
        marching_tet_.bind_uniform_state(program);
        
        glDeleteShader(vshader);
        glDeleteShader(gshader);
        glDeleteShader(fshader);        
    }


    /**
     * \brief The geometry shader for connectors.
     */
    static const char* GLUP150_gshader_connector_source =
        "void main() {                                                      \n"
        "    if(cell_is_clipped()) {                                        \n"
        "        return;                                                    \n"
        "    }                                                              \n"
        "    gl_PrimitiveID = gl_PrimitiveIDIn;                             \n"
        "    if(clipping_enabled() &&                                       \n"
        "        GLUP.clipping_mode == GLUP_CLIP_SLICE_CELLS) {             \n"
        "       draw_marching_cell();                                       \n"
        "       return;                                                     \n"
        "    }                                                              \n"
        "    project_vertices();                                            \n"
        "    bool do_clip = (clipping_enabled() &&                          \n"
        "                    GLUP.clipping_mode==GLUP_CLIP_STANDARD);       \n"
        "    flat_shaded_quad(0,1,3,2,do_clip);                             \n"
        "    flat_shaded_triangle(2,1,0,do_clip);                           \n"
        "    flat_shaded_triangle(3,2,0,do_clip);                           \n"
        "}                                                                  \n";
    
    void Context_GLSL150::setup_GLUP_CONNECTORS() {

        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            GLUP150_std(GLUP_CONNECTORS),
            GLUP150_vshader_in_out_declaration,
            GLUP150_vshader_transform_source,
            0
        );

        GLuint gshader = GLSL::compile_shader(
            GL_GEOMETRY_SHADER,
            GLUP150_std(GLUP_CONNECTORS),
            marching_connector_.GLSL_uniform_state_declaration(),
            "layout(lines_adjacency) in;                    \n",
            "layout(triangle_strip, max_vertices = 12) out; \n",
            GLUP150_gshader_in_out_declaration,
            "const int nb_vertices = 4; ",
            GLUP150_gshader_utils_source,
            marching_connector_.GLSL_compute_intersections(),
            GLUP150_marching_cells_utils,
            GLUP150_gshader_connector_source,
            0
        );

        GLuint fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER,
            GLUP150_std(GLUP_CONNECTORS),
            GLUP150_fshader_in_out_declaration,
            GLUP150_fshader_utils,                        
            GLUP150_triangle_fshader_source,
            0
        );

        GLuint program = 
            GLSL::create_program_from_shaders_no_link(
                vshader, gshader, fshader, 0 
            );
        
        set_primitive_info(GLUP_CONNECTORS, GL_LINES_ADJACENCY, program);
        marching_connector_.bind_uniform_state(program);
        
        glDeleteShader(vshader);
        glDeleteShader(gshader);
        glDeleteShader(fshader);        
    }

    // GLUP_PRISMS **********************************************************

    /**
     * \brief The geometry shader for prisms
     * \details Uses v_shader_transform and gshader_utils.
     */
    static const char* GLUP150_gshader_prism_source =
        "void main() {                                                      \n"
        "    if(cell_is_clipped()) {                                        \n"
        "        return;                                                    \n"
        "    }                                                              \n"
        "    gl_PrimitiveID = gl_PrimitiveIDIn;                             \n"
        "    if(clipping_enabled() &&                                       \n"
        "        GLUP.clipping_mode == GLUP_CLIP_SLICE_CELLS) {             \n"
        "       draw_marching_cell();                                       \n"
        "       return;                                                     \n"
        "    }                                                              \n"
        "    project_vertices();                                            \n"
        "    bool do_clip = (clipping_enabled() &&                          \n"
        "                    GLUP.clipping_mode==GLUP_CLIP_STANDARD);       \n"
        "    flat_shaded_triangle(0,1,2,do_clip);                           \n"
        "    flat_shaded_triangle(5,4,3,do_clip);                           \n"
        "    flat_shaded_quad(0,3,1,4,do_clip);                             \n"
        "    flat_shaded_quad(0,2,3,5,do_clip);                             \n"
        "    flat_shaded_quad(1,4,2,5,do_clip);                             \n" 
        "}                                                                  \n";

    
    void Context_GLSL150::setup_GLUP_PRISMS() {
        
        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            GLUP150_std(GLUP_PRISMS),
            GLUP150_vshader_in_out_declaration,
            GLUP150_vshader_transform_source,
            0
        );

        GLuint gshader = GLSL::compile_shader(
            GL_GEOMETRY_SHADER,
            GLUP150_std(GLUP_PRISMS),
            marching_prism_.GLSL_uniform_state_declaration(),
            "layout(triangles_adjacency) in;",
            "layout(triangle_strip, max_vertices = 18) out;",
            GLUP150_gshader_in_out_declaration,
            "const int nb_vertices = 6; ",
            GLUP150_gshader_utils_source,
            marching_prism_.GLSL_compute_intersections(),
            GLUP150_marching_cells_utils,
            GLUP150_gshader_prism_source,
            0
        );

        GLuint fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER,
            GLUP150_std(GLUP_PRISMS),
            GLUP150_fshader_in_out_declaration,
            GLUP150_fshader_utils,                        
            GLUP150_triangle_fshader_source,
            0
        );

        GLuint program = 
            GLSL::create_program_from_shaders_no_link(
                vshader, gshader, fshader, 0 
            );
        
        set_primitive_info(GLUP_PRISMS, GL_TRIANGLES_ADJACENCY, program);
        marching_prism_.bind_uniform_state(program);
        
        glDeleteShader(vshader);
        glDeleteShader(gshader);
        glDeleteShader(fshader);        
    }

    // GLUP_HEXAHEDRA *******************************************************

    static const char* GLUP150_vshader_gather_source =
        " const int nb_vertices_per_GL =                      \n"
        "                      nb_vertices / nb_vertices_GL;  \n"
        " in vec4 vertex_in[nb_vertices_per_GL];              \n"
        " in vec4 color_in[nb_vertices_per_GL];               \n"
        " in vec4 tex_coord_in[nb_vertices_per_GL];           \n"
        "                                                     \n"
        "out GVertexData {                                    \n"
        "    vec4 other_vertex[nb_vertices_per_GL-1];         \n"
        "    vec4 transformed[nb_vertices_per_GL];            \n"
        "    vec4 color[nb_vertices_per_GL];                  \n"
        "    vec4 tex_coord[nb_vertices_per_GL];              \n"
        "} VertexOut;                                         \n"
        "                                                     \n"
        "void main() {                                        \n"
        "   for(int i=1; i<nb_vertices_per_GL; ++i) {         \n"
        "       VertexOut.other_vertex[i-1] = vertex_in[i];   \n"
        "   }                                                 \n"
        "   for(int i=0; i<nb_vertices_per_GL; ++i) {         \n"
        "     VertexOut.transformed[i] =                      \n"
        "      GLUP.modelviewprojection_matrix * vertex_in[i];\n"
        "   }                                                 \n"
        "   if(texturing_enabled()) {                         \n"
        "       for(int i=0; i<nb_vertices_per_GL; ++i) {     \n"
        "           if(indirect_texturing_enabled()) {        \n"
        "              VertexOut.tex_coord[i] =               \n" 
        "                                    tex_coord_in[i]; \n"
        "           } else {                                  \n"
        "              VertexOut.tex_coord[i] =               \n"
        "              GLUP.texture_matrix * tex_coord_in[i]; \n"
        "           }                                         \n"
        "       }                                             \n"
        "   }                                                 \n"
        "   if(vertex_colors_enabled()) {                     \n"
        "       for(int i=0; i<nb_vertices_per_GL; ++i) {     \n"
        "           VertexOut.color[i] = color_in[i];         \n"
        "       }                                             \n"
        "   }                                                 \n"
        "   gl_Position = vertex_in[0];                       \n"
        "}                                                    \n"
        ;

    // To be used for primitives that have a number of vertices
    // that does not match existing OpenGL primitives. In that
    // case, the other vertices are gathered by the vertex shader,
    // and all vertices are merged into the attributes of a single
    // vertex. The geometry shader then expands this single vertex
    // into the primitive.
    static const char* GLUP150_gshader_gather_in_out_declaration =
        " const int nb_vertices_per_GL =            \n"
        "            nb_vertices / nb_vertices_GL;  \n"
        "in GVertexData {                           \n"
        "   vec4 other_vertex[nb_vertices_per_GL-1];\n"
        "   vec4 transformed[nb_vertices_per_GL];   \n"
        "   vec4 color[nb_vertices_per_GL];         \n"
        "   vec4 tex_coord[nb_vertices_per_GL];     \n"
        "} VertexIn[];                              \n"
        "                                           \n"
        "vec4 vertex_in(in int i) {                 \n"
        "   int i0 = i / nb_vertices_per_GL;        \n"
        "   int i1 = i % nb_vertices_per_GL;        \n"
        "   return (i1==0) ? gl_in[i0].gl_Position :\n"
        "          VertexIn[i0].other_vertex[i1-1]; \n"
        "}                                          \n"
        "                                           \n"
        "vec4 transformed_in(in int i) {            \n"
        "    int i0 = i / nb_vertices_per_GL;       \n"
        "    int i1 = i % nb_vertices_per_GL;       \n"
        "    return VertexIn[i0].transformed[i1];   \n"
        "}                                          \n"
        "                                           \n"
        "vec4 color_in(in int i) {                  \n"
        "    int i0 = i / nb_vertices_per_GL;       \n"
        "    int i1 = i % nb_vertices_per_GL;       \n"
        "    return VertexIn[i0].color[i1];         \n"
        "}                                          \n"
        "                                           \n"
        "vec4 tex_coord_in(in int i) {              \n"  
        "    int i0 = i / nb_vertices_per_GL;       \n"
        "    int i1 = i % nb_vertices_per_GL;       \n"
        "    return VertexIn[i0].tex_coord[i1];     \n"
        "}                                          \n"
        "                                           \n"        
        "out FragmentData {                         \n"
        "   vec4 color;                             \n"
        "   vec4 tex_coord;                         \n"
        "   flat float diffuse;                     \n"
        "   flat float specular;                    \n"
        "   vec2 bary;                              \n"
        "   flat vec3 edge_mask;                    \n"
        "} FragmentOut;                             \n"
        "                                           \n"
        "bool prim_is_discarded() {                 \n"
        "   return false;                           \n"
        "}                                          \n";

    /**
     * \brief The geometry shader for hexahedra
     * \details Uses v_shader_gather and gshader_utils.
     */
    static const char* GLUP150_gshader_hex_source =
        "void main() {                                                      \n"
        "    if(cell_is_clipped()) {                                        \n"
        "        return;                                                    \n"
        "    }                                                              \n"
        "    gl_PrimitiveID = gl_PrimitiveIDIn;                             \n"
        "    if(clipping_enabled() &&                                       \n"
        "        GLUP.clipping_mode == GLUP_CLIP_SLICE_CELLS) {             \n"
        "       draw_marching_cell();                                       \n"
        "       return;                                                     \n"
        "    }                                                              \n"
        "    project_vertices();                                            \n"
        "    bool do_clip = (clipping_enabled() &&                          \n"
        "                    GLUP.clipping_mode==GLUP_CLIP_STANDARD);       \n"
        "    flat_shaded_quad(0,2,4,6,do_clip);                             \n"
        "    flat_shaded_quad(3,1,7,5,do_clip);                             \n"
        "    flat_shaded_quad(1,0,5,4,do_clip);                             \n"
        "    flat_shaded_quad(2,3,6,7,do_clip);                             \n"
        "    flat_shaded_quad(1,3,0,2,do_clip);                             \n"
        "    flat_shaded_quad(4,6,5,7,do_clip);                             \n" 
        "}                                                                  \n";
    
    void Context_GLSL150::setup_GLUP_HEXAHEDRA() {

        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            GLUP150_std(GLUP_HEXAHEDRA),
            "const int nb_vertices = 8;",
            "const int nb_vertices_GL = 4;",
            GLUP150_vshader_gather_source,
            0
        );

        GLuint gshader = GLSL::compile_shader(
            GL_GEOMETRY_SHADER,
            GLUP150_std(GLUP_HEXAHEDRA),
            marching_hex_.GLSL_uniform_state_declaration(),            
            "const int nb_vertices = 8;",
            "const int nb_vertices_GL = 4;",
            "layout(lines_adjacency) in;",
            "layout(triangle_strip, max_vertices = 24) out;",
            GLUP150_gshader_gather_in_out_declaration,
            GLUP150_gshader_utils_source,
            marching_hex_.GLSL_compute_intersections(),
            GLUP150_marching_cells_utils,
            GLUP150_gshader_hex_source,
            0
        );

        GLuint fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER,
            GLUP150_std(GLUP_HEXAHEDRA),
            GLUP150_fshader_in_out_declaration,
            GLUP150_fshader_utils,                        
            GLUP150_triangle_fshader_source,
            0
        );

        GLuint program = 
            GLSL::create_program_from_shaders_no_link(
                vshader, gshader, fshader, 0 
            );

        set_primitive_info_vertex_gather_mode(
            GLUP_HEXAHEDRA, GL_LINES_ADJACENCY, program
        );
        
        marching_hex_.bind_uniform_state(program);
        
        glDeleteShader(vshader);
        glDeleteShader(gshader);
        glDeleteShader(fshader);                
    }
    
    // GLUP_PYRAMIDS ********************************************************

    /**
     * \brief The geometry shader for pyramids
     * \details Uses v_shader_gather and gshader_utils.
     */
    static const char* GLUP150_gshader_pyramid_source =
        "void main() {                                                      \n"
        "    if(cell_is_clipped()) {                                        \n"
        "        return;                                                    \n"
        "    }                                                              \n"
        "    gl_PrimitiveID = gl_PrimitiveIDIn;                             \n"
        "    if(clipping_enabled() &&                                       \n"
        "        GLUP.clipping_mode == GLUP_CLIP_SLICE_CELLS) {             \n"
        "       draw_marching_cell();                                       \n"
        "       return;                                                     \n"
        "    }                                                              \n"
        "    project_vertices();                                            \n"
        "    bool do_clip = (clipping_enabled() &&                          \n"
        "                    GLUP.clipping_mode==GLUP_CLIP_STANDARD);       \n"
        "    flat_shaded_quad(0,1,3,2,do_clip);                             \n"
        "    flat_shaded_triangle(0,4,1,do_clip);                           \n"
        "    flat_shaded_triangle(0,3,4,do_clip);                           \n"
        "    flat_shaded_triangle(2,4,3,do_clip);                           \n"
        "    flat_shaded_triangle(2,1,4,do_clip);                           \n"
        "}                                                                  \n";
    
    void Context_GLSL150::setup_GLUP_PYRAMIDS() {

        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            GLUP150_std(GLUP_PYRAMIDS),
            "const int nb_vertices = 5;",
            "const int nb_vertices_GL = 1;",            
            GLUP150_vshader_gather_source,
            0
        );

        GLuint gshader = GLSL::compile_shader(
            GL_GEOMETRY_SHADER,
            GLUP150_std(GLUP_PYRAMIDS),
            marching_pyramid_.GLSL_uniform_state_declaration(),            
            "const int nb_vertices = 5;",
            "const int nb_vertices_GL = 1;",
            "layout(points) in;",
            "layout(triangle_strip, max_vertices = 28) out;",
            GLUP150_gshader_gather_in_out_declaration,
            GLUP150_gshader_utils_source,
            marching_pyramid_.GLSL_compute_intersections(),
            GLUP150_marching_cells_utils,
            GLUP150_gshader_pyramid_source,
            0
        );

        GLuint fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER,
            GLUP150_std(GLUP_PYRAMIDS),
            GLUP150_fshader_in_out_declaration,
            GLUP150_fshader_utils,                        
            GLUP150_triangle_fshader_source,
            0
        );

        GLuint program = 
            GLSL::create_program_from_shaders_no_link(
                vshader, gshader, fshader, 0 
            );
        
        set_primitive_info_vertex_gather_mode(
            GLUP_PYRAMIDS, GL_POINTS, program
        );
        
        marching_pyramid_.bind_uniform_state(program);
        
        glDeleteShader(vshader);
        glDeleteShader(gshader);
        glDeleteShader(fshader);                
    }

    /***********************************************************************/
    /***********************************************************************/
    /**** GLUP implementation using GLSL 4.40                            ***/
    /***********************************************************************/
    /***********************************************************************/

    const char* Context_GLSL440::profile_name() const {
        return "GLUP440";
    }

#ifdef GEO_OS_APPLE
    static const char* GLUP440_shader_source_header =
        "#version 440 \n"
        ;
#else    
    static const char* GLUP440_shader_source_header =
        "#version 440 core \n"
        ;
#endif
    
#define GLUP440_std(prim)                    \
        GLUP440_shader_source_header,        \
        profile_dependent_declarations(),    \
        uniform_state_declaration(),         \
        primitive_declaration(prim).c_str(), \
        toggles_declaration()

     // This version of the tesselation shader gathers all input vertices into
     // a single vertex. It is used for pyramids.
     static const char* GLUP440_teshader_gather_single_vertex_source =
        "layout(isolines, point_mode) in;                                   \n"
        "                                                                   \n"
        "in VertexData {                                                    \n"
        "   vec4 transformed;                                               \n"
        "   vec4 color;                                                     \n"
        "   vec4 tex_coord;                                                 \n"
        "} VertexIn[];                                                      \n"
        "                                                                   \n"
        "out GTVertexData {                                                 \n"
        "    vec4 vertex[nb_vertices];                                      \n"
        "    vec4 transformed[nb_vertices];                                 \n"
        "    vec4 color[nb_vertices];                                       \n"
        "    vec4 tex_coord[nb_vertices];                                   \n"
        "    bool discard_me;                                               \n"
        "} VertexOut;                                                       \n"
        "                                                                   \n"
        "void main() {                                                      \n"
        "   VertexOut.discard_me = (gl_TessCoord.x > 0.5);                  \n"
        "   if(VertexOut.discard_me) {                                      \n"
        "        return;                                                    \n"
        "   }                                                               \n"
        "   for(int i=0; i<nb_vertices; ++i) {                              \n"
        "        VertexOut.vertex[i] = gl_in[i].gl_Position;                \n"
        "        VertexOut.transformed[i] = VertexIn[i].transformed;        \n"
        "   }                                                               \n"
        "   if(vertex_colors_enabled()) {                                   \n"
        "     for(int i=0; i<nb_vertices; ++i) {                            \n"
        "        VertexOut.color[i] = VertexIn[i].color;                    \n"
        "     }                                                             \n"
        "   }                                                               \n"
        "   if(texturing_enabled()) {                                       \n"
        "     for(int i=0; i<nb_vertices; ++i) {                            \n"
        "        VertexOut.tex_coord[i] = VertexIn[i].tex_coord;            \n"
        "     }                                                             \n"
        "   }                                                               \n"
        "}                                                                  \n"
        ;

     // This version of the tesselation shader gathers the input vertices into
     // several vertices. It is used for hexahedra.
     static const char* GLUP440_teshader_gather_multi_vertices_source =
        "layout(isolines) in;                                               \n"
        "                                                                   \n"
        " const int nb_vertices_per_GL = nb_vertices / nb_vertices_GL;      \n"
        "                                                                   \n" 
        "in VertexData {                                                    \n"
        "   vec4 transformed;                                               \n"
        "   vec4 color;                                                     \n"
        "   vec4 tex_coord;                                                 \n"
        "} VertexIn[];                                                      \n"
        "                                                                   \n"
        "out GTVertexData {                                                 \n"
        "    vec4 vertex[nb_vertices_per_GL];                               \n"
        "    vec4 transformed[nb_vertices_per_GL];                          \n"
        "    vec4 color[nb_vertices_per_GL];                                \n"
        "    vec4 tex_coord[nb_vertices_per_GL];                            \n"
        "    bool discard_me;                                               \n"
        "} VertexOut;                                                       \n"
        "                                                                   \n"
        "void main() {                                                      \n"
        "   int i0 = int(gl_TessCoord.x);                                   \n"
        "   for(int i1=0; i1<nb_vertices_per_GL; ++i1) {                    \n"
        "        int i = i0*nb_vertices_per_GL + i1;                        \n"
        "        VertexOut.vertex[i1] = gl_in[i].gl_Position;               \n"
        "        VertexOut.transformed[i1] = VertexIn[i].transformed;       \n"
        "   }                                                               \n"
        "   if(vertex_colors_enabled()) {                                   \n"
        "     for(int i1=0; i1<nb_vertices_per_GL; ++i1) {                  \n"
        "        int i = i0*nb_vertices_per_GL + i1;                        \n"
        "        VertexOut.color[i1] = VertexIn[i].color;                   \n"
        "     }                                                             \n"
        "   }                                                               \n"
        "   if(texturing_enabled()) {                                       \n"
        "     for(int i1=0; i1<nb_vertices_per_GL; ++i1) {                  \n"
        "        int i = i0*nb_vertices_per_GL + i1;                        \n"
        "        VertexOut.tex_coord[i1] = VertexIn[i].tex_coord;           \n"
        "     }                                                             \n"
        "   }                                                               \n"
        "}                                                                  \n"
        ;

    
    // The previous two shaders are to be used for primitives
    // that have a number of vertices that does not match
    // existing OpenGL primitives. In that
    // case, a tessellation shader fetches the vertices, using
    // GL_PATCHES primitive and generates two vertices (with a
    // lot of attributes). The second one needs to be discarded
    // (discard_me = true).
    
    static const char* GLUP440_gshader_tegather_in_out_declaration =
        " const int nb_vertices_per_GL =            \n"
        "       nb_vertices / nb_vertices_GL;       \n"
        "                                           \n"
        "in GTVertexData {                          \n"
        "    vec4 vertex[nb_vertices_per_GL];       \n"
        "    vec4 transformed[nb_vertices_per_GL];  \n"
        "    vec4 color[nb_vertices_per_GL];        \n"
        "    vec4 tex_coord[nb_vertices_per_GL];    \n"
        "    bool discard_me;                       \n"
        "} VertexIn[];                              \n"
        "                                           \n"
        "vec4 vertex_in(in int i) {                 \n"
        "   int i0 = i / nb_vertices_per_GL;        \n"
        "   int i1 = i % nb_vertices_per_GL;        \n"
        "   return VertexIn[i0].vertex[i1];         \n"
        "}                                          \n"
        "                                           \n"
        "vec4 transformed_in(in int i) {            \n"
        "   int i0 = i / nb_vertices_per_GL;        \n"
        "   int i1 = i % nb_vertices_per_GL;        \n"
        "   return VertexIn[i0].transformed[i1];    \n"
        "}                                          \n"
        "                                           \n"
        "vec4 color_in(in int i) {                  \n"
        "   int i0 = i / nb_vertices_per_GL;        \n"
        "   int i1 = i % nb_vertices_per_GL;        \n"
        "   return VertexIn[i0].color[i1];          \n"
        "}                                          \n"
        "                                           \n"
        "vec4 tex_coord_in(in int i) {              \n"
        "   int i0 = i / nb_vertices_per_GL;        \n"
        "   int i1 = i % nb_vertices_per_GL;        \n"
        "   return VertexIn[i0].tex_coord[i1];      \n"
        "}                                          \n"
        "                                           \n"        
        "out FragmentData {                         \n"
        "   vec4 color;                             \n"
        "   vec4 tex_coord;                         \n"        
        "   flat float diffuse;                     \n"
        "   flat float specular;                    \n"
        "   vec2 bary;                              \n"
        "   flat vec3 edge_mask;                    \n"
        "} FragmentOut;                             \n"
        "                                           \n"        
        "bool prim_is_discarded() {                 \n"
        "   if(nb_vertices_per_GL == nb_vertices) { \n"
        "     return VertexIn[0].discard_me;        \n"
        "   } else {                                \n"
        "     return false;                         \n"
        "   }                                       \n"
        "}                                          \n";
    
    void Context_GLSL440::setup_GLUP_HEXAHEDRA() {
        
        if(!GEO::CmdLine::get_arg_bool("gfx:GLSL_tesselation")) {
            Context_GLSL150::setup_GLUP_HEXAHEDRA();
            return;
        }

        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            GLUP150_std(GLUP_HEXAHEDRA),
            GLUP150_vshader_in_out_declaration,
            GLUP150_vshader_transform_source,
            0
        );

        GLuint teshader = GLSL::compile_shader(
            GL_TESS_EVALUATION_SHADER,
            GLUP440_std(GLUP_HEXAHEDRA),
            "const int nb_vertices = 8;",
            "const int nb_vertices_GL = 2;",
            GLUP440_teshader_gather_multi_vertices_source,
            0
        );
        
        GLuint gshader = GLSL::compile_shader(
            GL_GEOMETRY_SHADER,
            GLUP440_std(GLUP_HEXAHEDRA),
            marching_hex_.GLSL_uniform_state_declaration(),              
            "const int nb_vertices = 8;",
            "const int nb_vertices_GL = 2;",
            "layout(lines) in;",
            "layout(triangle_strip, max_vertices = 24) out;",
            GLUP440_gshader_tegather_in_out_declaration,            
            GLUP150_gshader_utils_source,
            marching_hex_.GLSL_compute_intersections(),
            GLUP150_marching_cells_utils,
            GLUP150_gshader_hex_source,
            0
        );

        GLuint fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER,
            GLUP150_std(GLUP_HEXAHEDRA),
            GLUP150_fshader_in_out_declaration,
            GLUP150_fshader_utils,                        
            GLUP150_triangle_fshader_source,
            0
        );

        GLuint program = 
            GLSL::create_program_from_shaders_no_link(
                vshader, teshader, gshader, fshader, 0 
            );

        set_primitive_info(
            GLUP_HEXAHEDRA, GL_PATCHES, program
        );
        marching_hex_.bind_uniform_state(program);
        
        glDeleteShader(teshader);
        glDeleteShader(vshader);        
        glDeleteShader(gshader);
        glDeleteShader(fshader);                
    }
    
    void Context_GLSL440::setup_GLUP_PYRAMIDS() {

        if(!GEO::CmdLine::get_arg_bool("gfx:GLSL_tesselation")) {
            Context_GLSL150::setup_GLUP_PYRAMIDS();
            return;
        }

        GLuint teshader = GLSL::compile_shader(
            GL_TESS_EVALUATION_SHADER,
            GLUP440_std(GLUP_PYRAMIDS),
            "const int nb_vertices = 5;",
            "const int nb_vertices_GL = 1;",            
            GLUP440_teshader_gather_single_vertex_source,
            0
        );

        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            GLUP150_std(GLUP_PYRAMIDS),
            GLUP150_vshader_in_out_declaration,
            GLUP150_vshader_transform_source,
            0
        );
        
        GLuint gshader = GLSL::compile_shader(
            GL_GEOMETRY_SHADER,
            GLUP440_std(GLUP_PYRAMIDS),
            marching_pyramid_.GLSL_uniform_state_declaration(),              
            "const int nb_vertices = 5;",
            "const int nb_vertices_GL = 1;",
            "layout(points) in;",            
            "layout(triangle_strip, max_vertices = 28) out;",            
            GLUP440_gshader_tegather_in_out_declaration,            
            GLUP150_gshader_utils_source,
            marching_pyramid_.GLSL_compute_intersections(),
            GLUP150_marching_cells_utils,
            GLUP150_gshader_pyramid_source,
            0
        );

        GLuint fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER,
            GLUP150_std(GLUP_PYRAMIDS),
            GLUP150_fshader_in_out_declaration,
            GLUP150_fshader_utils,                        
            GLUP150_triangle_fshader_source,
            0
        );

        GLuint program = 
            GLSL::create_program_from_shaders_no_link(
                vshader, teshader, gshader, fshader, 0 
            );

        set_primitive_info(
            GLUP_PYRAMIDS, GL_PATCHES, program
        );
        marching_pyramid_.bind_uniform_state(program);
        
        glDeleteShader(teshader);
        glDeleteShader(vshader);        
        glDeleteShader(gshader);
        glDeleteShader(fshader);                
    }

    /***********************************************************************/
}

#endif
