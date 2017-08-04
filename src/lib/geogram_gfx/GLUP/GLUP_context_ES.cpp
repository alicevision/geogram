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

#include <geogram_gfx/GLUP/GLUP_context_ES.h>
#include <geogram_gfx/basic/GLSL.h>
#include <geogram/basic/string.h>

#ifdef GEO_OS_EMSCRIPTEN
#include <emscripten.h>
#endif

#ifdef GEO_GL_ES2

/*
 * Notes:
 *
 * Tex coords pour le mesh, et picking:
 * ====================================
 * On pourrait passer de maniere optionelle un attribut de sommets
 * supplementaire, avec le numero du sommet (ou utiliser gl_VertexID
 * si OpenGL ES 3), ceci permettrait les choses suivantes:
 * 1) tex coords implicites (avec un petit tableau constant, et en
 *  utilisant gl_VertexID % nb_vertices_per_glup_primitive comme indice.
 * 2) picking, en utilisant 
 *  glup_primitiveID = gl_VertexID / nb_vertices_per_glup_primitive.
 * Si on fait \c{c}a, attention au code d'emulation de VAO, la fonction
 *  qui recuperer l'etat d'un attribut sous WebGL est en general buggee,
 *  il faudra memoriser l'etat soi-meme...
 *
 * Slice cells:
 * ============
 * Pour implanter le mode "slice cells", je me demande s'il vaut mieux
 *  - avoir un tout petit buffer, et faire un flush par cellule ?
 *  - tout assembler dans des gros buffers, et faire un flush tous les
 *    65K sommets ?
 *  -> essayer le premier (plus simple a programmer), et si \c{c}a rame
 *    faire le deuxi\`eme:
 *       Ca a l'air de marcher pas trop mal.
 *       Rem: le ELEMENT_ARRAY_BUFFER est en entiers 32 bits, certains archis
 *    (telephones etc...) peuvent avoir besoin d'un ELEMENT_ARRAY_BUFFER 16 bits
 *    (ce qu'on pourrait assez facilement avoir...)
 */

namespace GLUP {

    extern bool vertex_array_emulate;
    
    const char* Context_ES2::profile_name() const {
        return "GLUPES2";
    }

    void Context_ES2::prepare_to_draw(GLUPprimitive primitive) {
        Context::prepare_to_draw(primitive);
#ifdef GEO_GL_150        
        glEnable(GL_POINT_SPRITE);
#endif
    }

    void Context_ES2::done_draw(GLUPprimitive primitive) {
        Context::done_draw(primitive);
    }
    
    bool Context_ES2::primitive_supports_array_mode(GLUPprimitive prim) const {
        // Crashes if first frame is drawn like that, it
        // seems that it makes it skip an initialization, to
        // be checked... Or it is probably connected with this
        // mesh enabled / mesh disabled thingy...
        /*
        if(prim == GLUP_TRIANGLES) {
            return true;
        }
        */
        return Context::primitive_supports_array_mode(prim);
    }


    static const char* uniform_state_vs_ES2 = 
        " struct VSUniformState {              \n"
        "     mat3  normal_matrix;             \n"
        "     mat4 modelviewprojection_matrix; \n"
        "     mat4 modelview_matrix;           \n"
        "     mat4 texture_matrix;             \n"
        "     vec4  world_clip_plane;          \n"
        "     float point_size;                \n"
        " };                                   \n"
        "  uniform VSUniformState GLUP_VS;     \n"        
        ;
        
    static const char* uniform_state_fs_ES2 =
        " struct UniformState {                \n"        
        "     bool vertex_colors_enabled;      \n"        
        "                                      \n"        
        "     vec4  front_color;               \n"
        "     vec4  back_color;                \n"
        "                                      \n"        
        "     bool  draw_mesh_enabled;         \n"
        "     vec4  mesh_color;                \n"
        "     float mesh_width;                \n"

        "     bool  lighting_enabled;          \n"
        "     vec3  light_vector;              \n"
        "     vec3  light_half_vector;         \n"

        "     bool texturing_enabled;          \n"
        "     bool indirect_texturing_enabled; \n"             
        "     int  texture_mode;               \n"
        "     int  texture_type;               \n"        
        "                                      \n"        
        "     float cells_shrink;              \n"

        "     bool  picking_enabled;           \n"         
        "     int   picking_mode;              \n"
        "     int   picking_id;                \n" 
        "     int   base_picking_id;           \n" 

        "     bool  clipping_enabled;          \n"
        "     int   clipping_mode;             \n"
        "     vec4  clip_plane;                \n"

        "  };                                  \n"
        
        "  uniform UniformState GLUP;          \n"
        "  uniform sampler2D texture1Dsampler; \n"                
        "  uniform sampler2D texture2Dsampler; \n"        
        ;

    static const char* glup_constants = 
        "  const int GLUP_CLIP_STANDARD         = 0;   \n"
        "  const int GLUP_CLIP_WHOLE_CELLS      = 1;   \n"
        "  const int GLUP_CLIP_STRADDLING_CELLS = 2;   \n"
        "  const int GLUP_CLIP_SLICE_CELLS      = 3;   \n"

        "  const int GLUP_TEXTURE_1D = 1;              \n"
        "  const int GLUP_TEXTURE_2D = 2;              \n"
        "  const int GLUP_TEXTURE_3D = 3;              \n"

        "  const int GLUP_TEXTURE_REPLACE  = 0;        \n"
        "  const int GLUP_TEXTURE_MODULATE = 1;        \n"
        "  const int GLUP_TEXTURE_ADD      = 2;        \n"

        "  const int GLUP_PICK_PRIMITIVE   = 1;        \n"
        "  const int GLUP_PICK_CONSTANT    = 2;        \n"

        "  const int GLUP_POINTS     =0;               \n"
        "  const int GLUP_LINES      =1;               \n"
        "  const int GLUP_TRIANGLES  =2;               \n"
        "  const int GLUP_QUADS      =3;               \n"
        "  const int GLUP_TETRAHEDRA =4;               \n"
        "  const int GLUP_HEXAHEDRA  =5;               \n"
        "  const int GLUP_PRISMS     =6;               \n"
        "  const int GLUP_PYRAMIDS   =7;               \n"
        "  const int GLUP_CONNECTORS =8;               \n"                     
        "  const int GLUP_NB_PRIMITIVES = 9;           \n"
        ;

    static const char* fshader_utils = 
        "void output_lighting(                                              \n"
        "   in float sdiffuse, in float spec                                \n"
        ") {                                                                \n"
        "   if(sdiffuse > 0.0) {                                            \n"
        "       vec3 vspec = spec*vec3(1.0,1.0,1.0);                        \n"
        "       glup_FragColor =                                            \n"
        "              (sdiffuse*glup_FragColor) + vec4(vspec,1.0);         \n"
        "       glup_FragColor.rgb += vec3(0.2, 0.2, 0.2);                  \n"
        "   } else {                                                        \n"
        "       glup_FragColor = vec4(0.2, 0.2, 0.2, 1.0);                  \n"
        "   }                                                               \n"
        "}                                                                  \n"

        "float min3(float x, float y, float z) {                            \n"
        "   return min(min(x,y),z);                                         \n"
        "}                                                                  \n"

        "float min4(float x, float y, float z, float w) {                   \n"
        "   return min(min(x,y),min(z,w));                                  \n"
        "}                                                                  \n"
        
        "float edge_factor1(float bary) {                                   \n"
        "   float d = fwidth(bary);                                         \n"
        "   float a = smoothstep(0.0, d*GLUP.mesh_width, bary);             \n"
        "   return a;                                                       \n"
        "}                                                                  \n"

        "float edge_factor3(vec3 bary) {                                    \n"
        "   vec3 d = fwidth(bary);                                          \n"
        "   vec3 a = smoothstep(                                            \n"
        "      vec3(0.0, 0.0, 0.0), d*GLUP.mesh_width, bary                 \n"
        "   );                                                              \n"
        "   return min3(a.x, a.y ,a.z);                                     \n"
        "}                                                                  \n"

        "float edge_factor4(vec4 bary) {                                    \n"
        "   vec4 d = fwidth(bary);                                          \n"
        "   vec4 a = smoothstep(                                            \n"
        "      vec4(0.0, 0.0, 0.0, 0.0), d*GLUP.mesh_width, bary            \n"
        "   );                                                              \n"
        "   return min4(a.x, a.y, a.z, a.w);                                \n"
        "}                                                                  \n"

        "float cell_edge_factor(vec2 bary) {                                \n"
        "   return edge_factor1(1.0-(1.0 - bary.x)*(1.0 - bary.y));         \n"
        "}                                                                  \n"
        
        ;

    static const char* fshader_header =
#ifdef GEO_OS_EMSCRIPTEN
    "#extension GL_OES_standard_derivatives : enable \n"        
    "precision mediump float;                        \n"
    "vec4 glup_FragColor;                            \n"
    "void glup_write_fragment() {                    \n"
    "    gl_FragColor = glup_FragColor;              \n"
    "}                                               \n"
    "vec4 glup_texture(                              \n"
    "    in sampler2D samp, in vec2 uv               \n"
    ") {                                             \n"
    "    return texture2D(samp, uv);                 \n"
    "}                                               \n"
#else
    "#version 330                                    \n"
    "out vec4 glup_FragColor;                        \n"
    "vec4 glup_texture(                              \n"
    "    in sampler2D samp, in vec2 uv               \n"
    ") {                                             \n"
    "    return texture(samp, uv);                   \n"
    "}                                               \n"
    "void glup_write_fragment() { }                  \n"
#endif
    ;

    static const char* vshader_header =
#ifndef GEO_OS_EMSCRIPTEN
    "#version 330 \n"    
#endif        
    "  "
    ;


    /**
     * \brief Computes the number of elements in the buffer for
     *  a given primitive.
     * \param[in] nb_vertices_per_glup_primitive number of vertices
     *  per primitive
     * \param[in] nb_elements_per_glup_primitive number of elements
     *  per primitive
     * \return total number of elements required to draw all the
     *  primitives in the buffer.
     */
    inline index_t nb_elements(
        index_t nb_vertices_per_glup_primitive,
        index_t nb_elements_per_glup_primitive
    ) {
        index_t nb_glup_primitives = 
            IMMEDIATE_BUFFER_SIZE / nb_vertices_per_glup_primitive;
        return (
            nb_glup_primitives * nb_elements_per_glup_primitive
        );
    }

    
    void Context_ES2::setup() {

#ifdef GEO_OS_EMSCRIPTEN

/**
 * \brief Tests whether a WebGL extension is supported.       
 * \details WebGL requires getExtension() to be called for initializing
 *  the extension (if supported, this returns an extension object,
 *  else this returns null).
 * \param[in] ext the name of the extension surrounded by single quotes
 * \retval true if the extension is supported
 * \retval false otherwise
 */
#define WEBGL_EXTENSION_SUPPORTED(ext)                      \
        EM_ASM_INT({return (                                \
             Module.ctx.getExtension(ext) == null ? 0 : 1   \
        );},0) != 0 

        
        GL_OES_standard_derivatives_ = WEBGL_EXTENSION_SUPPORTED(
            'OES_standard_derivatives'
        );

        GL_OES_vertex_array_object_ = WEBGL_EXTENSION_SUPPORTED(
            'OES_vertex_array_object'
        );

        GL_EXT_frag_depth_ = WEBGL_EXTENSION_SUPPORTED(
            'EXT_frag_depth'
        );
#else
        GL_OES_standard_derivatives_ =
            extension_is_supported("GL_OES_standard_derivatives");

        GL_OES_vertex_array_object_ =
            extension_is_supported("GL_OES_vertex_array_object");

        GL_EXT_frag_depth_ =
            extension_is_supported("GL_EXT_frag_depth");
#endif

#ifdef GEO_OS_APPLE
        GL_OES_vertex_array_object_ = true;
        GL_OES_standard_derivatives_ = true;
#endif
        
        //   Switch-on GLUP's emulation for Vertex Array
        // Objects if they are not supported.
        if(!GL_OES_vertex_array_object_) {
            GLUP::vertex_array_emulate = true;
        }

#ifdef GLUP_DEBUG        
        GLUP::vertex_array_emulate = true;
#endif

        Logger::out("GLUP") 
            << "GL_OES_standard_derivatives: "
            << GL_OES_standard_derivatives_
            << std::endl;

        Logger::out("GLUP")
            << "GL_OES_vertex_array_object: "
            << GL_OES_vertex_array_object_
            << std::endl;

        Logger::out("GLUP")
            << "GL_EXT_frag_depth:"
            << GL_EXT_frag_depth_
            << std::endl;
            
        create_CPU_side_uniform_buffer();


        /************** VAOs and VBOs for whole cells clipping ******/
        
        //   Compute the maximum number of elements required to draw
        // the vertices in a buffer, for all types of volumetric
        // primitives.
        index_t max_nb_elements = 0;
        max_nb_elements = geo_max(max_nb_elements, nb_elements(4,12));
        max_nb_elements = geo_max(max_nb_elements, nb_elements(8,36));
        max_nb_elements = geo_max(max_nb_elements, nb_elements(6,24));
        max_nb_elements = geo_max(max_nb_elements, nb_elements(5,18));        
        max_nb_elements = geo_max(max_nb_elements, nb_elements(4,12));        
        
        nb_clip_cells_elements_ = max_nb_elements;        
        clip_cells_elements_ = new Numeric::uint16[nb_clip_cells_elements_];

        glupGenVertexArrays(1,&clip_cells_VAO_);
        glupBindVertexArray(clip_cells_VAO_);
        
        bind_immediate_state_buffers_to_VAO();

        glGenBuffers(1, &clip_cells_elements_VBO_);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, clip_cells_elements_VBO_);

        update_buffer_object(
            clip_cells_elements_VBO_,
            GL_ELEMENT_ARRAY_BUFFER,
            max_nb_elements * sizeof(Numeric::uint32),
            nil // no need to copy the buffer, it will be overwritten after.
        );
        
        glupBindVertexArray(0);
        
        /************** VAOs and VBOs for sliced cells clipping ******/

        max_nb_elements = 6;
        index_t max_nb_vertices = 12;
        
        glupGenVertexArrays(1,&sliced_cells_VAO_);
        glupBindVertexArray(sliced_cells_VAO_);

        glGenBuffers(1, &sliced_cells_elements_VBO_);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, sliced_cells_elements_VBO_);

        update_buffer_object(
            sliced_cells_elements_VBO_,
            GL_ELEMENT_ARRAY_BUFFER,
            max_nb_elements * sizeof(Numeric::uint16),
            nil // no need to copy the buffer, it will be overwritten after.
        );
        
        for(index_t i=0; i<3; ++i) {
            glGenBuffers(1, &sliced_cells_vertex_attrib_VBO_[i]);
            glBindBuffer(GL_ARRAY_BUFFER, sliced_cells_vertex_attrib_VBO_[i]);

            update_buffer_object(
                sliced_cells_vertex_attrib_VBO_[i],
                GL_ARRAY_BUFFER,
                max_nb_vertices * sizeof(Numeric::float32) * 4,
                nil // no need to copy the buffer, it will be overwritten after.
            );
            
            glVertexAttribPointer(
                i,
                4,
                GL_FLOAT, GL_FALSE,
                0,
                0
            );
        }
        glEnableVertexAttribArray(0);
        glupBindVertexArray(0);
    }

    Context_ES2::Context_ES2() :
        nb_clip_cells_elements_(0),
        clip_cells_elements_(nil),        
        clip_cells_elements_VBO_(0),
        clip_cells_VAO_(0),
        sliced_cells_elements_VBO_(0),
        sliced_cells_VAO_(0) {
        GL_OES_standard_derivatives_ = false;        
        GL_OES_vertex_array_object_ = false;
        GL_EXT_frag_depth_ = false;
        for(index_t i=0; i<3; ++i) {
            sliced_cells_vertex_attrib_VBO_[i] = 0;
        }
        uses_mesh_tex_coord_ = true;        
    }
    
    Context_ES2::~Context_ES2() {
        glDeleteBuffers(1,&clip_cells_elements_VBO_);
        glupDeleteVertexArrays(1, &clip_cells_VAO_);
        for(index_t i=0; i<3; ++i) {
            glDeleteBuffers(1, &sliced_cells_vertex_attrib_VBO_[i]);
        }
        glupDeleteVertexArrays(1, &sliced_cells_VAO_);
        delete[] clip_cells_elements_;
    }
    
    Memory::pointer Context_ES2::get_state_variable_address(
        const char* name
    ) {
        geo_assert(variable_to_offset_.find(name) != variable_to_offset_.end());
        return uniform_buffer_data_ + variable_to_offset_[name];
    }

    void Context_ES2::do_update_uniform_buffer() {
        GEO_CHECK_GLUP();
        
        if(!uniform_buffer_dirty_) {
            return;
        }

        GEO_CHECK_GLUP();
        
        if(matrices_dirty_) {
            update_matrices();
        }

        GEO_CHECK_GLUP();
        
        if(lighting_dirty_) {
            update_lighting();
        }

        GEO_CHECK_GLUP();
        
        glUseProgram(latest_program_);

        GEO_CHECK_GLUP();

        if(latest_program_ != 0) {
            copy_uniform_state_to_current_program();
        }
        
        GEO_CHECK_GLUP();
        
        uniform_buffer_dirty_ = false;
    }

    void Context_ES2::copy_uniform_state_to_current_program() {
        GLint loc;

        // Matrices (in vertex shader).
        {
            loc = glGetUniformLocation(
                latest_program_, "GLUP_VS.modelviewprojection_matrix"
            );
            glUniformMatrix4fv(
                loc, 1, GL_FALSE,
                uniform_state_.modelviewprojection_matrix.get_pointer()
            );
            loc = glGetUniformLocation(
                latest_program_, "GLUP_VS.modelview_matrix"
            );
            glUniformMatrix4fv(
                loc, 1, GL_FALSE,
                uniform_state_.modelview_matrix.get_pointer()
            );
            loc = glGetUniformLocation(
                latest_program_, "GLUP_VS.texture_matrix"
            );
            glUniformMatrix4fv(
                loc, 1, GL_FALSE,
                uniform_state_.texture_matrix.get_pointer()
            );
            loc = glGetUniformLocation(
                latest_program_, "GLUP_VS.normal_matrix"
            );

            // Normal matrix is stored in state buffer with padding.
            float normal_matrix[9];
            for(index_t i=0; i<3; ++i) {
                for(index_t j=0; j<3; ++j) {
                    normal_matrix[i*3+j] =
                        uniform_state_.normal_matrix.get_pointer()[i*4+j];
                }
            }
            
            glUniformMatrix3fv(loc, 1, GL_FALSE, normal_matrix);
        }

        // Points (in vertex shader).
        {
            loc = glGetUniformLocation(
                latest_program_, "GLUP_VS.point_size"
            );
            glUniform1f(loc, uniform_state_.point_size.get());
        }
        
        
        // Lighting.
        {
            loc = glGetUniformLocation(
                latest_program_, "GLUP.light_vector"
            );
            glUniform3fv(
                loc, 1, uniform_state_.light_vector.get_pointer()
            );
            loc = glGetUniformLocation(
                latest_program_, "GLUP.light_half_vector"
            );
            glUniform3fv(
                loc, 1, uniform_state_.light_half_vector.get_pointer()
            );
        }

        // All the toggles.
        for(index_t i=0; i<uniform_state_.toggle.size(); ++i) {
            loc = glGetUniformLocation(
                latest_program_,
                ("GLUP." + uniform_state_.toggle[i].name()).c_str()
            );
            glUniform1i(
                loc, GLint(uniform_state_.toggle[i].get())
            );
        }

        // Colors.
        if(!uniform_state_.toggle[GLUP_VERTEX_COLORS].get()) {
            for(index_t i=0; i<uniform_state_.color.size(); ++i) {
                loc = glGetUniformLocation(
                    latest_program_,
                    ("GLUP." + uniform_state_.color[i].name()).c_str()
                    );
                glUniform4fv(
                    loc, 1, uniform_state_.color[i].get_pointer()
                );
            }
        }

        // Mesh.
        if(uniform_state_.toggle[GLUP_DRAW_MESH].get()) {
            loc = glGetUniformLocation(latest_program_, "GLUP.mesh_width");
            glUniform1f(loc, uniform_state_.mesh_width.get());
        }

        // Texturing.
        if(uniform_state_.toggle[GLUP_TEXTURING].get()) {
            loc = glGetUniformLocation(
                latest_program_, "GLUP.texture_mode"
            );
            glUniform1i(loc, uniform_state_.texture_mode.get());
            loc = glGetUniformLocation(
                latest_program_, "GLUP.texture_type"
            );
            glUniform1i(loc, uniform_state_.texture_type.get());
        }

        // Cell shrink
        loc = glGetUniformLocation(latest_program_, "GLUP.cells_shrink");
        glUniform1f(loc, uniform_state_.cells_shrink.get());

        // Picking
        if(uniform_state_.toggle[GLUP_PICKING].get()) {
            loc = glGetUniformLocation(
                latest_program_, "GLUP.picking_mode"
            );
            glUniform1i(loc, uniform_state_.picking_mode.get());
            loc = glGetUniformLocation(
                latest_program_, "GLUP.picking_id"
            );
            glUniform1i(loc, uniform_state_.picking_id.get());
            loc = glGetUniformLocation(
                latest_program_, "GLUP.base_picking_id"
            );
            glUniform1i(loc, uniform_state_.base_picking_id.get());
        }

        // Clipping.
        if(uniform_state_.toggle[GLUP_CLIPPING].get()) {
            loc = glGetUniformLocation(
                latest_program_, "GLUP.clipping_mode"
            );
            glUniform1i(loc, uniform_state_.clipping_mode.get());
            loc = glGetUniformLocation(
                latest_program_, "GLUP.clip_plane"
            );
            glUniform4fv(loc, 1, uniform_state_.clip_plane.get_pointer());
            loc = glGetUniformLocation(
                latest_program_, "GLUP_VS.world_clip_plane"
            );
            glUniform4fv(loc, 1, uniform_state_.world_clip_plane.get_pointer());
        }
    }


/*****************************************************************************/

    static const char* points_vshader_source =
#ifdef GEO_OS_EMSCRIPTEN        
        "attribute vec4 vertex_in;                             \n"
        "attribute vec4 color_in;                              \n"
        "attribute vec4 tex_coord_in;                          \n"
        "varying vec4 color;                                   \n"
        "varying vec4 tex_coord;                               \n"
        "varying float clip_dist;                              \n"
#else
        "in vec4 vertex_in;                                    \n"
        "in vec4 color_in;                                     \n"
        "in vec4 tex_coord_in;                                 \n"
        "out vec4 color;                                       \n"
        "out vec4 tex_coord;                                   \n"
        "out float clip_dist;                                  \n"
#endif
        
        "void main() {                                         \n"
        "   if(maybe_clipping_enabled()) {                     \n"
        "      clip_dist = dot(                                \n"
        "          vertex_in, GLUP_VS.world_clip_plane         \n"
        "      );                                              \n"
        "   }                                                  \n"
        "   if(maybe_vertex_colors_enabled()) {                \n"
        "       color = color_in;                              \n"
        "   }                                                  \n"
        "   if(maybe_texturing_enabled()) {                    \n"
        "   tex_coord = GLUP_VS.texture_matrix * tex_coord_in; \n"
        "   }                                                  \n"
        "   gl_PointSize = GLUP_VS.point_size;                 \n"
        "   gl_Position =                                      \n"
        "     GLUP_VS.modelviewprojection_matrix * vertex_in;  \n"
        "}                                                     \n"
        ;


    static const char* update_depth =
        "   void update_depth(in vec3 N) {                             \n"
        "        gl_FragDepthEXT = gl_FragCoord.z - 0.001 * N.z;       \n"
        "   }                                                          \n"
        ;

    static const char* no_update_depth =
        "   void update_depth(in vec3 N) {                             \n"
        "   }                                                          \n"
        ;
                
    static const char* points_fshader_source =
#ifdef GEO_OS_EMSCRIPTEN        
        "varying vec4 color;                                           \n"
        "varying vec4 tex_coord;                                       \n"
        "varying float clip_dist;                                      \n"
#else
        "in vec4 color;                                                \n"
        "in vec4 tex_coord;                                            \n"
        "in float clip_dist;                                           \n"
#endif        
        "void main() {                                                 \n"
        "   if(clipping_enabled() && (clip_dist < 0.0)) {              \n"
        "      discard;                                                \n"
        "   }                                                          \n"
        "   if(vertex_colors_enabled()) {                              \n"
        "      glup_FragColor = color;                                 \n"
        "   } else {                                                   \n"
        "      glup_FragColor = GLUP.front_color;                      \n"
        "   }                                                          \n"
        "   if(texturing_enabled()) {                                  \n"
        "      vec4 tex_color;                                         \n"
        "      if(GLUP.texture_type == GLUP_TEXTURE_1D) {              \n"
        "           tex_color = glup_texture(                          \n"
        "               texture1Dsampler, tex_coord.xy                 \n"
        "           );                                                 \n"
        "      } else {                                                \n"
        "           tex_color = glup_texture(                          \n"
        "                texture2Dsampler, tex_coord.xy                \n"
        "           );                                                 \n"
        "      }                                                       \n"
        "      if(GLUP.texture_mode == GLUP_TEXTURE_REPLACE) {         \n"
        "             glup_FragColor = tex_color;                      \n"
        "      } else if(GLUP.texture_mode == GLUP_TEXTURE_MODULATE) { \n"
        "             glup_FragColor *= tex_color;                     \n"
        "      } else {                                                \n"
        "             glup_FragColor += tex_color;                     \n"
        "      }                                                       \n"
        "   }                                                          \n"
        "   vec2 V = 2.0*(gl_PointCoord - vec2(0.5, 0.5));             \n"
        "   float one_minus_r2 = 1.0 - dot(V,V);                       \n"
        "   if(one_minus_r2 < 0.0) {                                   \n"
        "      discard;                                                \n"
        "   }                                                          \n"
        "   vec3 N = vec3(V.x, -V.y, sqrt(one_minus_r2));              \n"
        "   update_depth(N);                                           \n"
        "   float diff = dot(N,GLUP.light_vector);                     \n"
        "   float spec = dot(N,GLUP.light_half_vector);                \n"
        "   spec = pow(spec,30.0);                                     \n"
        "   output_lighting(diff,spec);                                \n"
        "   glup_write_fragment();                                     \n"
        "}                                                             \n"
        ;
    
    void Context_ES2::setup_GLUP_POINTS() {
        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            vshader_header,            
            uniform_state_vs_ES2,
            glup_constants,
            potential_toggles_declaration(),            
            points_vshader_source,
            0
        );

        GLuint fshader = 0;
        if(GL_EXT_frag_depth_) {
            fshader = GLSL::compile_shader(
                GL_FRAGMENT_SHADER,
                "#extension GL_EXT_frag_depth : enable\n",            
                fshader_header,
                uniform_state_fs_ES2,
                glup_constants,
                toggles_declaration(),
                fshader_utils,
                update_depth,
                points_fshader_source,
                0
            );
        } else {
            fshader = GLSL::compile_shader(
                GL_FRAGMENT_SHADER,
                fshader_header,
                uniform_state_fs_ES2,
                glup_constants,
                toggles_declaration(),
                fshader_utils,
                no_update_depth,                
                points_fshader_source,
                0
            );
        }

        GLuint program = GLSL::create_program_from_shaders_no_link(
            vshader, fshader, 0 
        );

        glDeleteShader(vshader);
        glDeleteShader(fshader);
        
        set_primitive_info(GLUP_POINTS, GL_POINTS, program);
    }

/*****************************************************************************/

    static const char* lines_vshader_source =
#ifdef GEO_OS_EMSCRIPTEN        
        "attribute vec4 vertex_in;                         \n"
        "attribute vec4 color_in;                          \n"
        "attribute vec4 tex_coord_in;                      \n"
        "varying vec4 color;                               \n"
        "varying vec4 tex_coord;                           \n"
        "varying float clip_dist;                          \n"
#else
        "in vec4 vertex_in;                                \n"
        "in vec4 color_in;                                 \n"
        "in vec4 tex_coord_in;                             \n"
        "out vec4 color;                                   \n"
        "out vec4 tex_coord;                               \n"
        "out float clip_dist;                              \n"
#endif        
        "void main() {                                           \n"
        "   if(maybe_clipping_enabled()) {                       \n"
        "      clip_dist = dot(                                  \n"
        "          vertex_in, GLUP_VS.world_clip_plane           \n"
        "      );                                                \n"
        "   }                                                    \n"
        "   if(maybe_vertex_colors_enabled()) {                  \n"
        "       color = color_in;                                \n"
        "   }                                                    \n"
        "   if(maybe_texturing_enabled()) {                      \n"
        "     tex_coord = GLUP_VS.texture_matrix * tex_coord_in; \n"
        "   }                                                    \n"
        "   gl_Position =                                        \n"
        "     GLUP_VS.modelviewprojection_matrix * vertex_in;    \n"
        "}                                                       \n"
        ;

    static const char* lines_fshader_source =
#ifdef GEO_OS_EMSCRIPTEN        
        "varying vec4 color;                                           \n"
        "varying vec4 tex_coord;                                       \n"
        "varying float clip_dist;                                      \n"
#else
        "in vec4 color;                                                \n"
        "in vec4 tex_coord;                                            \n"
        "in float clip_dist;                                           \n"
#endif        
        "void main() {                                                 \n"
        "   if(clipping_enabled() && (clip_dist < 0.0)) {              \n"
        "      discard;                                                \n"
        "   }                                                          \n"
        "   if(vertex_colors_enabled()) {                              \n"
        "      glup_FragColor = color;                                 \n"
        "   } else {                                                   \n"
        "      glup_FragColor = GLUP.front_color;                      \n"
        "   }                                                          \n"
        "   if(texturing_enabled()) {                                  \n"
        "      vec4 tex_color;                                         \n"
        "      if(GLUP.texture_type == GLUP_TEXTURE_1D) {              \n"
        "           tex_color = glup_texture(                          \n"
        "               texture1Dsampler, tex_coord.xy                 \n"
        "           );                                                 \n"
        "      } else {                                                \n"
        "           tex_color = glup_texture(                          \n"
        "                texture2Dsampler, tex_coord.xy                \n"
        "           );                                                 \n"
        "      }                                                       \n"
        "      if(GLUP.texture_mode == GLUP_TEXTURE_REPLACE) {         \n"
        "             glup_FragColor = tex_color;                      \n"
        "      } else if(GLUP.texture_mode == GLUP_TEXTURE_MODULATE) { \n"
        "             glup_FragColor *= tex_color;                     \n"
        "      } else {                                                \n"
        "             glup_FragColor += tex_color;                     \n"
        "      }                                                       \n"
        "   }                                                          \n"
        "   glup_write_fragment();                                     \n"
        "}                                                             \n"
        ;
    
    void Context_ES2::setup_GLUP_LINES() {
        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            vshader_header,            
            uniform_state_vs_ES2,
            glup_constants,
            potential_toggles_declaration(),            
            lines_vshader_source,
            0
        );

        GLuint fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER,
            fshader_header,
            uniform_state_fs_ES2,
            glup_constants,
            toggles_declaration(),
            fshader_utils,
            lines_fshader_source,
            0
        );

        GLuint program = GLSL::create_program_from_shaders_no_link(
            vshader, fshader, 0 
        );

        glDeleteShader(vshader);
        glDeleteShader(fshader);
        
        set_primitive_info(GLUP_LINES, GL_LINES, program);
    }

/*****************************************************************************/

    static const char* polygons_vshader_source =
#ifdef GEO_OS_EMSCRIPTEN        
        "attribute vec4 vertex_in;                             \n"
        "attribute vec4 color_in;                              \n"
        "attribute vec4 tex_coord_in;                          \n"
        "attribute vec4 mesh_tex_coord_in;                     \n"
        "varying vec3 vertex_clip_space;                       \n"
        "varying float clip_dist;                              \n"        
        "varying vec4 color;                                   \n"
        "varying vec4 tex_coord;                               \n"
        "varying vec4 mesh_tex_coord;                          \n"
#else
        "in vec4 vertex_in;                                    \n"
        "in vec4 color_in;                                     \n"
        "in vec4 tex_coord_in;                                 \n"
        "in vec4 mesh_tex_coord_in;                            \n"
        "out vec3 vertex_clip_space;                           \n"
        "out float clip_dist;                                  \n"
        "out vec4 color;                                       \n"
        "out vec4 tex_coord;                                   \n"
        "out vec4 mesh_tex_coord;                              \n"
#endif        

        "void main() {                                         \n"
        "   if(maybe_clipping_enabled()) {                     \n"
        "      clip_dist = dot(                                \n"
        "          vertex_in, GLUP_VS.world_clip_plane         \n"
        "      );                                              \n"
        "   }                                                  \n"
        "   if(maybe_lighting_enabled()) {                     \n"        
        "      vertex_clip_space =                             \n"
        "         (GLUP_VS.modelview_matrix * vertex_in).xyz;  \n"
        "   }                                                  \n"
        "   if(maybe_vertex_colors_enabled()) {                \n"
        "      color = color_in;                               \n"
        "   }                                                  \n"
        "   if(maybe_texturing_enabled()) {                    \n"
        "   tex_coord = GLUP_VS.texture_matrix * tex_coord_in; \n"
        "   }                                                  \n"
        "   if(maybe_draw_mesh_enabled()) {                    \n"
        "      mesh_tex_coord=mesh_tex_coord_in;               \n"
        "   }                                                  \n"
        "   gl_Position =                                      \n"
        "     GLUP_VS.modelviewprojection_matrix * vertex_in;  \n"
        "}                                                     \n"
        ;

    static const char* polygons_fragment_clipping =
        "void clip_fragment(in float dist) {          \n"
        "   if(                                       \n"
        "     clipping_enabled() &&                   \n"
        "       dist < 0.0                            \n"
        "   ) {                                       \n"
        "      discard;                               \n"
        "   }                                         \n"
        "}                                            \n"
        ;

    static const char* polyhedra_fragment_clipping =
        "void clip_fragment(in float dist) {          \n"
        "   if(                                       \n"
        "     clipping_enabled() &&                   \n"
        "     GLUP.clipping_mode ==                   \n"
        "                       GLUP_CLIP_STANDARD && \n"
        "          dist < 0.0                         \n"
        "   ) {                                       \n"
        "      discard;                               \n"
        "   }                                         \n"
        "}                                            \n"        
        ;
    

    
    // Clipping and lighting are done in the fragment shader
    // it is not classical, but:
    //  - I do not want to compute normals on the CPU side.
    //  - OpenGLES2 does not have geometry shaders.
    // It is probably not optimum, but it works reasonably
    // well, and satisfies my goal of directly playing back
    // a vertex buffer object.

    
    static const char* polygons_fshader_source =
#ifdef GEO_OS_EMSCRIPTEN        
        "varying vec3 vertex_clip_space;              \n"
        "varying float clip_dist;                     \n"
        "varying vec4 color;                          \n"
        "varying vec4 tex_coord;                      \n"
        "varying vec4 mesh_tex_coord;                 \n"
#else
        "in vec3 vertex_clip_space;                   \n"
        "in float clip_dist;                          \n"
        "in vec4 color;                               \n"
        "in vec4 tex_coord;                           \n"
        "in vec4 mesh_tex_coord;                      \n"
#endif        
        "void main() {                                \n"
        "   clip_fragment(clip_dist);                 \n"
        "   if(vertex_colors_enabled()) {             \n"
        "      glup_FragColor = color;                \n"
        "   } else {                                  \n"
        "      glup_FragColor = gl_FrontFacing ?      \n"
        "         GLUP.front_color : GLUP.back_color; \n"
        "   }                                         \n"
        "   if(texturing_enabled()) {                 \n"
        "      vec4 tex_color;                                         \n"
        "      if(GLUP.texture_type == GLUP_TEXTURE_1D) {              \n"
        "           tex_color = glup_texture(                          \n"
        "               texture1Dsampler, tex_coord.xy                 \n"
        "           );                                                 \n"
        "      } else {                                                \n"
        "           tex_color = glup_texture(                          \n"
        "                texture2Dsampler, tex_coord.xy                \n"
        "           );                                                 \n"
        "      }                                                       \n"
        "      if(GLUP.texture_mode == GLUP_TEXTURE_REPLACE) {         \n"
        "             glup_FragColor = tex_color;                      \n"
        "      } else if(GLUP.texture_mode == GLUP_TEXTURE_MODULATE) { \n"
        "             glup_FragColor *= tex_color;                     \n"
        "      } else {                                                \n"
        "             glup_FragColor += tex_color;                     \n"
        "      }                                                       \n"
        "   }                                                          \n"
        "   if(lighting_enabled()) {                                   \n"
        "      vec3 U = dFdx(vertex_clip_space);                       \n"
        "      vec3 V = dFdy(vertex_clip_space);                       \n"
        "      vec3 N = normalize(cross(U,V));                         \n"
        "      float diff = dot(N,GLUP.light_vector);                  \n"
        "      float spec = dot(N,GLUP.light_half_vector);             \n"
        "      spec = pow(spec,30.0);                                  \n"
        "      output_lighting(diff,spec);                             \n"
        "   }                                                          \n"
        "   if(                                                        \n"
        "      draw_mesh_enabled() &&                                  \n"
        "        (!clipping_enabled() ||                               \n"
        "         GLUP.clipping_mode != GLUP_CLIP_SLICE_CELLS)         \n"
        "    ) {                                                       \n"
        "       glup_FragColor = mix(                                  \n"
        "            GLUP.mesh_color,                                  \n"
        "            glup_FragColor,                                   \n"
        "            edge_factor(mesh_tex_coord)                       \n"
        "       );                                                     \n"
        "   }                                                          \n"
        "   glup_write_fragment();                                     \n"
        "}                                                             \n"
        ;

/*****************************************************************************/

    static const char* triangle_edge_factor_source =
        " float edge_factor(in vec4 bary) {    \n"
        "    return edge_factor3(bary.xyz);    \n"
        " }                                    \n"
        ;
    
    void Context_ES2::setup_GLUP_TRIANGLES() {
        static index_t element_indices[3]  = {
            0, 1, 2
        };

        static GLUPfloat tex_coords[12] = {
            1.0f, 0.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 1.0f, 0.0f            
        };

        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            vshader_header,
            uniform_state_vs_ES2,
            glup_constants,
            potential_toggles_declaration(),            
            polygons_vshader_source,
            0
        );

        GLuint fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER,
            fshader_header,
            uniform_state_fs_ES2,
            glup_constants,
            toggles_declaration(),
            fshader_utils,
            triangle_edge_factor_source,
            polygons_fragment_clipping,
            polygons_fshader_source,
            0
        );

        GLuint program = GLSL::create_program_from_shaders_no_link(
            vshader, fshader, 0 
        );

        glDeleteShader(vshader);
        glDeleteShader(fshader);

        set_primitive_info_immediate_index_mode(
            GLUP_TRIANGLES, GL_TRIANGLES, program,
            3, element_indices, tex_coords
        );
    }

/*****************************************************************************/

    static const char* quad_edge_factor_source =
        " float edge_factor(in vec4 bary) {                             \n"
        "    return edge_factor4(bary);                                 \n"
        " }                                                             \n"
        ;
    
    void Context_ES2::setup_GLUP_QUADS() {
        
        static index_t element_indices[6]  = {
            0, 1, 2,
            0, 2, 3
        };

        static GLUPfloat tex_coords[16] = {
            0.0f, 0.0f, 1.0f, 1.0f,
            1.0f, 0.0f, 0.0f, 1.0f,
            1.0f, 1.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 1.0f, 0.0f,            
        };

        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            vshader_header,            
            uniform_state_vs_ES2,
            glup_constants,
            potential_toggles_declaration(),            
            polygons_vshader_source,
            0
        );

        GLuint fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER,
            fshader_header,
            uniform_state_fs_ES2,
            glup_constants,
            toggles_declaration(),
            fshader_utils,
            quad_edge_factor_source,
            polygons_fragment_clipping,            
            polygons_fshader_source,
            0
        );

        GLuint program = GLSL::create_program_from_shaders_no_link(
            vshader, fshader, 0 
        );

        glDeleteShader(vshader);
        glDeleteShader(fshader);
        
        set_primitive_info_immediate_index_mode(
            GLUP_QUADS, GL_TRIANGLES, program,
            6, element_indices, tex_coords
        );        
    }

/*****************************************************************************/

    static const char* tet_edge_factor_source =
        " float edge_factor(in vec4 bary) {                             \n"
        "   float e1 = cell_edge_factor(bary.xy);                       \n"
        "   float e2 = cell_edge_factor(bary.xz);                       \n"
        "   float e3 = cell_edge_factor(bary.xw);                       \n"
        "   float e4 = cell_edge_factor(bary.yz);                       \n"
        "   float e5 = cell_edge_factor(bary.yw);                       \n"
        "   float e6 = cell_edge_factor(bary.zw);                       \n"
        "   return min3(min(e1,e2),min(e3,e4),min(e5,e6));              \n"
        " }                                                             \n"
        ;

    void Context_ES2::setup_GLUP_TETRAHEDRA() {
        static index_t element_indices[12] = {
            1,3,2,
            0,2,3,
            3,1,0,
            0,1,2
        };
        
        static GLUPfloat tex_coords[16] = {
            1.0f, 0.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 1.0f, 0.0f,
            0.0f, 0.0f, 0.0f, 1.0f            
        };

        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            vshader_header,            
            uniform_state_vs_ES2,
            glup_constants,
            potential_toggles_declaration(),            
            polygons_vshader_source,
            0
        );

        GLuint fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER,
            fshader_header,
            uniform_state_fs_ES2,
            glup_constants,
            toggles_declaration(),
            fshader_utils,
            tet_edge_factor_source,
            polyhedra_fragment_clipping,                        
            polygons_fshader_source,
            0
        );

        GLuint program = GLSL::create_program_from_shaders_no_link(
            vshader, fshader, 0 
        );

        glDeleteShader(vshader);
        glDeleteShader(fshader);
        
        set_primitive_info_immediate_index_mode(
            GLUP_TETRAHEDRA, GL_TRIANGLES, program,
            12, element_indices, tex_coords
        );        
    }

/*****************************************************************************/

    static const char* hex_edge_factor_source =
        " float edge_factor(in vec4 bary) {               \n"
        "     vec3 u = bary.xyz;                          \n"
        "     vec3 U = vec3(1.0, 1.0, 1.0) - u;           \n"
        "     float e1 = cell_edge_factor(vec2(u.x,u.y)); \n"
        "     float e2 = cell_edge_factor(vec2(u.x,u.z)); \n"
        "     float e3 = cell_edge_factor(vec2(u.x,U.y)); \n"
        "     float e4 = cell_edge_factor(vec2(u.x,U.z)); \n"

        "     float e5 = cell_edge_factor(vec2(U.x,u.y)); \n"
        "     float e6 = cell_edge_factor(vec2(U.x,u.z)); \n"
        "     float e7 = cell_edge_factor(vec2(U.x,U.y)); \n"
        "     float e8 = cell_edge_factor(vec2(U.x,U.z)); \n"

        "     float e9  = cell_edge_factor(vec2(u.y,u.z)); \n"
        "     float e10 = cell_edge_factor(vec2(u.y,U.z)); \n"
        "     float e11 = cell_edge_factor(vec2(U.y,u.z)); \n"
        "     float e12 = cell_edge_factor(vec2(U.y,U.z)); \n"
        
        "     float r1 = min4(e1,e2,e3,e4);                \n"
        "     float r2 = min4(e5,e6,e7,e8);                \n"
        "     float r3 = min4(e9,e10,e11,e12);             \n"
        
        "     return min3(r1,r2,r3);                       \n"
        " }                                                \n"
        ;
    
    void Context_ES2::setup_GLUP_HEXAHEDRA() {

        static index_t element_indices[36] = {
            0,2,6, 0,6,4,
            3,1,5, 3,5,7,
            1,0,4, 1,4,5,
            2,3,7, 2,7,6,
            1,3,2, 1,2,0,
            4,6,7, 4,7,5
        };

        static GLUPfloat tex_coords[32] = {
            0.0f, 0.0f, 0.0f, 0.0f,            
            0.0f, 0.0f, 1.0f, 0.0f,
            0.0f, 1.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 1.0f, 0.0f,
            1.0f, 0.0f, 0.0f, 0.0f,            
            1.0f, 0.0f, 1.0f, 0.0f,
            1.0f, 1.0f, 0.0f, 0.0f,
            1.0f, 1.0f, 1.0f, 0.0f,
        };

        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            vshader_header,            
            uniform_state_vs_ES2,
            glup_constants,
            potential_toggles_declaration(),            
            polygons_vshader_source,
            0
        );

        GLuint fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER,
            fshader_header,
            uniform_state_fs_ES2,
            glup_constants,
            toggles_declaration(),
            fshader_utils,
            hex_edge_factor_source,
            polyhedra_fragment_clipping,                        
            polygons_fshader_source,
            0
        );

        GLuint program = GLSL::create_program_from_shaders_no_link(
            vshader, fshader, 0 
        );

        glDeleteShader(vshader);
        glDeleteShader(fshader);
        
        set_primitive_info_immediate_index_mode(
            GLUP_HEXAHEDRA, GL_TRIANGLES, program,
            36, element_indices, tex_coords
        );        
    }

/*****************************************************************************/

    static const char* prism_edge_factor_source =
        " float edge_factor(in vec4 bary) {                             \n"
        "   vec4 bary2 = vec4(bary.x, bary.y, bary.z, 1.0 - bary.w);    \n"
        "   float e1 = cell_edge_factor(bary.xw);                       \n"
        "   float e2 = cell_edge_factor(bary.yw);                       \n"
        "   float e3 = cell_edge_factor(bary.zw);                       \n"
        "   float e4 = cell_edge_factor(bary2.xw);                      \n"
        "   float e5 = cell_edge_factor(bary2.yw);                      \n"
        "   float e6 = cell_edge_factor(bary2.zw);                      \n"
        "   float e7 = cell_edge_factor(bary.xy);                       \n"
        "   float e8 = cell_edge_factor(bary.yz);                       \n"
        "   float e9 = cell_edge_factor(bary.zx);                       \n"
        "   return min(min3(e7,e8,e9),                                  \n"
        "              min3(min(e1,e2),min(e3,e4),min(e5,e6))           \n"
        "          );                                                   \n"
        " }                                                             \n"
        ;
    
    void Context_ES2::setup_GLUP_PRISMS() {
        static index_t element_indices[24] = {
            0,1,2,
            3,5,4,
            0,3,4, 0,4,1,
            0,2,5, 0,5,3,
            1,4,5, 1,5,2
        };

        static GLUPfloat tex_coords[32] = {
            1.0f, 0.0f, 0.0f, 0.0f,            
            0.0f, 1.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 1.0f, 0.0f,
            1.0f, 0.0f, 0.0f, 1.0f,            
            0.0f, 1.0f, 0.0f, 1.0f,
            0.0f, 0.0f, 1.0f, 1.0f
        };
        
        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            vshader_header,            
            uniform_state_vs_ES2,
            glup_constants,
            potential_toggles_declaration(),            
            polygons_vshader_source,
            0
        );

        GLuint fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER,
            fshader_header,
            uniform_state_fs_ES2,
            glup_constants,
            toggles_declaration(),
            fshader_utils,
            prism_edge_factor_source,
            polyhedra_fragment_clipping,                                    
            polygons_fshader_source,
            0
        );

        GLuint program = GLSL::create_program_from_shaders_no_link(
            vshader, fshader, 0 
        );

        glDeleteShader(vshader);
        glDeleteShader(fshader);
        
        set_primitive_info_immediate_index_mode(
            GLUP_PRISMS, GL_TRIANGLES, program,
            24, element_indices, tex_coords
        );        
    }

    /************************************************************************/

    // Draw mesh not implemented for pyramids because
    // I did not find yet a nice way of defining barycentric
    // coordinates that I could interpret in the shader.
    
    static const char* no_edge_factor_source =
        " float edge_factor(in vec4 bary) {    \n"
        "    return 1.0;                       \n"
        " }                                    \n"
        ;

    void Context_ES2::setup_GLUP_PYRAMIDS() {
        static index_t element_indices[18] = {
            0,1,2, 0,2,3,
            0,4,1,
            0,3,4,
            2,4,3,
            2,1,4
        };

        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            vshader_header,            
            uniform_state_vs_ES2,
            glup_constants,
            potential_toggles_declaration(),            
            polygons_vshader_source,
            0
        );

        GLuint fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER,
            fshader_header,
            uniform_state_fs_ES2,
            glup_constants,
            toggles_declaration(),
            fshader_utils,
            no_edge_factor_source,
            polyhedra_fragment_clipping,                                    
            polygons_fshader_source,
            0
        );

        GLuint program = GLSL::create_program_from_shaders_no_link(
            vshader, fshader, 0 
        );

        glDeleteShader(vshader);
        glDeleteShader(fshader);
        
        set_primitive_info_immediate_index_mode(
            GLUP_PYRAMIDS, GL_TRIANGLES, program,
            18, element_indices
        );        
    }

    void Context_ES2::setup_GLUP_CONNECTORS() {
        static index_t element_indices[12] = {
            0,1,2, 0,2,3,
            2,1,0,
            3,2,0
        };

        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            vshader_header,            
            uniform_state_vs_ES2,
            glup_constants,
            potential_toggles_declaration(),            
            polygons_vshader_source,
            0
        );

        GLuint fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER,
            fshader_header,
            uniform_state_fs_ES2,
            glup_constants,
            toggles_declaration(),
            fshader_utils,
            no_edge_factor_source,
            polyhedra_fragment_clipping,                                    
            polygons_fshader_source,
            0
        );

        GLuint program = GLSL::create_program_from_shaders_no_link(
            vshader, fshader, 0 
        );

        glDeleteShader(vshader);
        glDeleteShader(fshader);
        
        set_primitive_info_immediate_index_mode(
            GLUP_CONNECTORS, GL_TRIANGLES, program,
            12, element_indices
        );        
    }

    bool Context_ES2::cell_by_cell_clipping() const {
        if(!uniform_state_.toggle[GLUP_CLIPPING].get()) {
            return false;
        }

        if(
            immediate_state_.primitive() != GLUP_TETRAHEDRA &&
            immediate_state_.primitive() != GLUP_HEXAHEDRA &&
            immediate_state_.primitive() != GLUP_PRISMS &&
            immediate_state_.primitive() != GLUP_PYRAMIDS &&
            immediate_state_.primitive() != GLUP_CONNECTORS 
        ) {
            return false;
        }
        
        return (
            uniform_state_.clipping_mode.get() == GLUP_CLIP_WHOLE_CELLS ||
            uniform_state_.clipping_mode.get() == GLUP_CLIP_STRADDLING_CELLS
        );
    }

    bool Context_ES2::sliced_cells_clipping() const {
        if(!uniform_state_.toggle[GLUP_CLIPPING].get()) {
            return false;
        }

        if(
            immediate_state_.primitive() != GLUP_TETRAHEDRA &&
            immediate_state_.primitive() != GLUP_HEXAHEDRA &&
            immediate_state_.primitive() != GLUP_PRISMS &&
            immediate_state_.primitive() != GLUP_PYRAMIDS &&
            immediate_state_.primitive() != GLUP_CONNECTORS 
        ) {
            return false;
        }
        
        return (
            uniform_state_.clipping_mode.get() == GLUP_CLIP_SLICE_CELLS 
        );
    }
    
    void Context_ES2::flush_immediate_buffers_with_cell_by_cell_clipping() {
        index_t cur_vertex = 0;
        index_t cur_element_out = 0;
        index_t nb_vertices_per_cell =
            nb_vertices_per_primitive[immediate_state_.primitive()];
        index_t nb_elements_per_cell = primitive_info_[
            immediate_state_.primitive()].nb_elements_per_primitive;

        // Assemble the list of indices (elements) that correspond
        // to the vertices of the visible cells.
        while(cur_vertex < immediate_state_.nb_vertices()) {
            if(!cell_is_clipped(cur_vertex)) {
                for(index_t lei=0; lei < nb_elements_per_cell; ++lei) {
                    geo_debug_assert(
                        cur_element_out < nb_clip_cells_elements_
                    );
                    clip_cells_elements_[cur_element_out] =
                        Numeric::uint16(
                            cur_vertex + primitive_info_[
                                immediate_state_.primitive()
                            ].primitive_elements[lei]
                        );
                    ++cur_element_out;
                }
            }
            cur_vertex += nb_vertices_per_cell;
        }

        glupBindVertexArray(clip_cells_VAO_);

        // Stream the values in all bound buffers (and attach
        // them to the VAO).
        for(index_t i=0; i<immediate_state_.buffer.size(); ++i) {
            if(immediate_state_.buffer[i].is_enabled()) {
                glEnableVertexAttribArray(i);
                if(immediate_state_.buffer[i].VBO() != 0) {
                    stream_buffer_object(
                        immediate_state_.buffer[i].VBO(),
                        GL_ARRAY_BUFFER,
                        immediate_state_.buffer[i].size_in_bytes(),
                        immediate_state_.buffer[i].data()
                    );
                }
            }
        }

        // If mesh should be drawn, then bind the VBO that has the
        // texture coordinates used by mesh drawing.
        if(
            uniform_state_.toggle[GLUP_DRAW_MESH].get() &&
            primitive_info_[immediate_state_.primitive()].tex_coords_VBO
            != 0
        ) {
            glEnableVertexAttribArray(3);
            glBindBuffer(
                GL_ARRAY_BUFFER,
                primitive_info_[immediate_state_.primitive()].tex_coords_VBO
            );
            glVertexAttribPointer(
                immediate_state_.buffer.size(), 
                4,
                GL_FLOAT,
                GL_FALSE,
                0,  // stride
                0   // pointer (relative to bound VBO beginning)
            );
        } else {
            glDisableVertexAttribArray(3);                
        }

        // Stream the indices into the elements VBO.
        stream_buffer_object(
            clip_cells_elements_VBO_,
            GL_ELEMENT_ARRAY_BUFFER,
            cur_element_out * sizeof(Numeric::uint16),
            clip_cells_elements_
        );
        
        glDrawElements(
            primitive_info_[immediate_state_.primitive()].GL_primitive,
            GLsizei(cur_element_out),
            GL_UNSIGNED_SHORT,
            0
        );
        
        glupBindVertexArray(0);
        glDisableVertexAttribArray(3);
        immediate_state_.reset();
    }

    void Context_ES2::flush_immediate_buffers_with_sliced_cells_clipping() {
        
        MarchingCell* marching_cell = nil;
        switch(immediate_state_.primitive()) {
        case GLUP_TETRAHEDRA:
            marching_cell = &marching_tet_;
            break;
        case GLUP_HEXAHEDRA:
            marching_cell = &marching_hex_;
            break;
        case GLUP_PRISMS:
            marching_cell = &marching_prism_;
            break;
        case GLUP_PYRAMIDS:
            marching_cell = &marching_pyramid_;
            break;
        case GLUP_CONNECTORS:
            marching_cell = &marching_connector_;
            break;
        case GLUP_POINTS:
        case GLUP_LINES:
        case GLUP_TRIANGLES:
        case GLUP_QUADS:
        case GLUP_NB_PRIMITIVES:
            geo_assert_not_reached;
        }

        glupBindVertexArray(sliced_cells_VAO_);

        index_t v0=0;
        while(v0 < immediate_state_.nb_vertices()) {
            index_t config = get_config(v0, marching_cell->nb_vertices());

            //   Compute all the intersection vertices (plus their
            // attributes if enabled).
            for(index_t i=0; i<marching_cell->config_size(config); ++i) {
                index_t e = marching_cell->config_edges(config)[i];
                compute_intersection(
                    v0+marching_cell->edge_vertex(e,0),
                    v0+marching_cell->edge_vertex(e,1),
                    e
                );
            }

            if(marching_cell->config_size(config) != 0) {
                
                stream_buffer_object(
                    sliced_cells_elements_VBO_,
                    GL_ELEMENT_ARRAY_BUFFER,
                    marching_cell->max_config_size() * sizeof(Numeric::uint32),
                    marching_cell->config_edges(config)
                );

                stream_buffer_object(
                    sliced_cells_vertex_attrib_VBO_[0],
                    GL_ARRAY_BUFFER,
                    12 * sizeof(Numeric::float32) * 4,
                    &isect_vertex_attribute_[0][0]
                );

                for(index_t i=1; i<3; ++i) {
                    if(immediate_state_.buffer[i].is_enabled()) {
                        glEnableVertexAttribArray(i);
                        stream_buffer_object(
                            sliced_cells_vertex_attrib_VBO_[i],
                            GL_ARRAY_BUFFER,
                            12 * sizeof(Numeric::float32) * 4,
                            &isect_vertex_attribute_[i][0]
                        );
                    } else {
                        glDisableVertexAttribArray(i);
                    }
                }

                glDrawElements(
                    GL_TRIANGLE_FAN,
                    GLsizei(marching_cell->config_size(config)), 
                    GL_UNSIGNED_INT,
                    0
                );
            }
            
            v0 += marching_cell->nb_vertices();
        }
        
        glupBindVertexArray(0);        
        immediate_state_.reset();
    }
    
    void Context_ES2::flush_immediate_buffers() {
        classify_vertices_in_immediate_buffers();        
        shrink_cells_in_immediate_buffers();
        if(cell_by_cell_clipping()) {
            flush_immediate_buffers_with_cell_by_cell_clipping();            
        } else if(sliced_cells_clipping()) {
            flush_immediate_buffers_with_sliced_cells_clipping();
        } {
            Context::flush_immediate_buffers();            
        }
    }
}

#endif

