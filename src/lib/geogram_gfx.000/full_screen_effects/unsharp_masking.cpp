/*
 *  Copyright (c) 2012-2014, Bruno Levy
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


#include <geogram_gfx/full_screen_effects/unsharp_masking.h>
#include <geogram_gfx/basic/GLSL.h>
#include <geogram/basic/logger.h>
#include <geogram/bibliography/bibliography.h>

namespace GEO {

    UnsharpMaskingImpl::UnsharpMaskingImpl() {
        unsharp_masking_program_ = 0;
        blur_program_ = 0;
        positive_shadows_ = true;
        contrast_ = 50;
        intensity_ = 50;
        halos_ = true;
        blur_width_ = 1;
	geo_cite("DBLP:journals/tog/LuftCD06");
	ES_profile_ = true;
    }
    
    UnsharpMaskingImpl::~UnsharpMaskingImpl() {
        glDeleteTextures(1, &depth_tex_);    
        depth_tex_ = 0;
        if (unsharp_masking_program_ != 0) { 
            glDeleteProgram(unsharp_masking_program_);       
        }
        if (blur_program_ != 0) { 
            glDeleteProgram(blur_program_);       
        }
    } 

    double UnsharpMaskingImpl::required_GLSL_version() const {
#ifdef GEO_OS_EMSCRIPTEN
	return 1.0;
#else	
        return 1.3;
#endif	
    }
    
    void UnsharpMaskingImpl::initialize(index_t w, index_t h) {
	
        FullScreenEffectImpl::initialize(w,h);
	
        if(!OK()) {
            return;
        }

        glGenTextures(1, &depth_tex_);
        glBindTexture(GL_TEXTURE_2D, depth_tex_);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

        // previously: const GLint internal_format = GL_LUMINANCE_FLOAT32_ATI;
        // Now deprecated, we use "red-only" texture format with floats
#ifdef GEO_OS_EMSCRIPTEN
        const GLint internal_format = GL_RGBA;
#else
        const GLint internal_format = GL_R32F;
#endif       
        
        if(!blur_1_.initialize(
	       width(), height(), false, internal_format
	   )
	) {
            Logger::err("UnsharpM")
                << "blur_1_ FBO is not initialized" << std::endl;
        }
        if(!blur_2_.initialize(
	       width(), height(), false, internal_format
	   )
	) {
            Logger::err("UnsharpM")
                << "blur_2_ FBO is not initialized" << std::endl;
        }


	Logger::out("SFX") << "Initialized FBOs" << std::endl;

	Logger::out("SFX") << "Compiling blur shader..." << std::endl;
	Logger::out("SFX") << "ES_profile_ = " << ES_profile_ << std::endl;
	
	// Shader sources are embedded in source code,
	// Initial sourcecode is in:
	// geogram_gfx/GLUP/shaders/fullscreen
	blur_program_ = GLSL::compile_program_with_includes_no_link(
	    this,
	    "//stage GL_VERTEX_SHADER\n"
	    "//import <fullscreen/vertex_shader.h>\n",
	    "//stage GL_FRAGMENT_SHADER\n"
	    "//import <fullscreen/blur_fragment_shader.h>\n"                
	);

        glBindAttribLocation(blur_program_, 0, "vertex_in");
        glBindAttribLocation(blur_program_, 1, "tex_coord_in");

        GLSL::link_program(blur_program_);
        
        GLSL::set_program_uniform_by_name(
            blur_program_, "source_tex", 0
        );
        GLSL::set_program_uniform_by_name(
            blur_program_, "blur_width", 2.0f
        );
	if(ES_profile_) {
	    GLSL::set_program_uniform_by_name(
		blur_program_, "source_tex_size",
		float(blur_1_.width), float(blur_1_.height)
	    );
	}
	
	Logger::out("SFX") << "Compiling unsharp masking shader..." << std::endl;
	Logger::out("SFX") << "ES_profile_ = " << ES_profile_ << std::endl;
	
	// Shader sources are embedded in source code,
	// Initial sourcecode is in:
	// geogram_gfx/GLUP/shaders/fullscreen
	unsharp_masking_program_ = GLSL::compile_program_with_includes_no_link(
	    this,
	    "//stage GL_VERTEX_SHADER\n"
	    "//import <fullscreen/vertex_shader.h>\n",
	    "//stage GL_FRAGMENT_SHADER\n"
	    "//import <fullscreen/unsharp_masking_fragment_shader.h>\n"                
	);

	Logger::out("SFX") << "Compiled shaders" << std::endl;	
        
        glBindAttribLocation(unsharp_masking_program_, 0, "vertex_in");
        glBindAttribLocation(unsharp_masking_program_, 1, "tex_coord_in");

        GLSL::link_program(unsharp_masking_program_);

	Logger::out("SFX") << "Linked shaders" << std::endl;
	
        GLSL::set_program_uniform_by_name(
            unsharp_masking_program_, "blur_texture", 0
        );
        GLSL::set_program_uniform_by_name(
            unsharp_masking_program_, "depth_texture", 1
        );

        update();
    }

    void UnsharpMaskingImpl::pre_render(index_t w, index_t h) {
        FullScreenEffectImpl::pre_render(w,h);
    }

    void UnsharpMaskingImpl::resize(index_t width, index_t height) {
        blur_1_.resize(width, height);
        blur_2_.resize(width, height);
                
        glBindTexture(GL_TEXTURE_2D, depth_tex_);
        glTexImage2D(
            GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT,
            GLsizei(width), GLsizei(height), 0,
            GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE, 0
        );
        glBindTexture(GL_TEXTURE_2D, 0);
        FullScreenEffectImpl::resize(width, height);
    }

    void UnsharpMaskingImpl::copy_depth_buffer_to_texture() {
        glBindTexture(GL_TEXTURE_2D, depth_tex_);
        glCopyTexImage2D(
            GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, 0, 0,
            GLsizei(width()), GLsizei(height()), 0
        );
    }

    void UnsharpMaskingImpl::copy_depth_texture_to_FBO() {
        //   Default Graphite viewport does not cover the [-1,1]x[-1,1]
        // normalized device coordinates space.
        glViewport(0, 0, GLsizei(width()), GLsizei(height()));
        blur_1_.bind_as_framebuffer();
        glBindTexture(GL_TEXTURE_2D, depth_tex_);
        draw_unit_textured_quad();
        blur_1_.unbind();
    }
    
    void UnsharpMaskingImpl::blur() {
        // Horizontal blur: blur_1_ -> blur_2_

        blur_2_.bind_as_framebuffer();
	
        GLSL::set_program_uniform_by_name(blur_program_, "vertical", false);
        glUseProgram(blur_program_);

        blur_1_.bind_as_texture();
	
        draw_unit_textured_quad();
        blur_2_.unbind();
	
        // Vertical blur: blur_2_ -> blur_1_
        blur_1_.bind_as_framebuffer();
        GLSL::set_program_uniform_by_name(blur_program_, "vertical", true);
        glUseProgram(blur_program_);

        blur_2_.bind_as_texture();
        draw_unit_textured_quad();	
        blur_1_.unbind();
	
        glUseProgram(0);
    }

    void UnsharpMaskingImpl::display_final_texture() {
        glDisable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        glUseProgram(unsharp_masking_program_);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, depth_tex_);
        glActiveTexture(GL_TEXTURE0);

        blur_1_.bind_as_texture();
        draw_unit_textured_quad();
        glUseProgram(0);

        glEnable(GL_DEPTH_TEST);
        glDisable(GL_BLEND);
        blur_1_.unbind();
    }

    void UnsharpMaskingImpl::post_render() {
#ifndef GEO_OS_EMSCRIPTEN       
	if(!core_profile_) {
	    glDisable(GL_CLIP_PLANE0);
	    glDisable(GL_LIGHTING);
	}
#endif       
        copy_depth_buffer_to_texture();
        copy_depth_texture_to_FBO();
        blur();
        display_final_texture();
	reset_alpha();
        FullScreenEffectImpl::post_render();
    }

    void UnsharpMaskingImpl::update() {
        if(!OK()) {
            return;
        }
        if(unsharp_masking_program_ != 0) {
            float fx = 1.0f - float(contrast_)/100.0f;
            geo_clamp(fx, 0.0f, 1.0f);
            GLSL::set_program_uniform_by_name(
                unsharp_masking_program_,
                "shadows_gamma", fx
            );
            GLSL::set_program_uniform_by_name(
                unsharp_masking_program_,
                "shadows_intensity", float(intensity_) / 400.0f
            );
            GLSL::set_program_uniform_by_name(
                unsharp_masking_program_, "shadows_halo", halos_
            );
            GLSL::set_program_uniform_by_name(
                unsharp_masking_program_, "do_positive_shadows",
                positive_shadows_
            );
        }
        if(blur_program_ != 0) {
            GLSL::set_program_uniform_by_name(
                blur_program_, "blur_width", float(blur_width_)
            );
        }
    }
}
