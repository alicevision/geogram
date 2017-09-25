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


#include <geogram_gfx/full_screen_effects/full_screen_effect.h>
#include <geogram_gfx/basic/GLSL.h>
#include <geogram_gfx/GLUP/shaders/embedded_shaders.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/logger.h>

namespace GEO {

    FullScreenEffectImpl::FullScreenEffectImpl() :
        initialized_(false),
        OK_(true),
        width_(0),
        height_(0),
        core_profile_(false),
        ES_profile_(false) {
    }

    double FullScreenEffectImpl::required_GLSL_version() const {
        return 1.0;
    }

    void FullScreenEffectImpl::get_vertex_shader_preamble_pseudo_file(
        std::vector<GLSL::Source>& sources
    ) {
        if(ES_profile_) {
            sources.push_back(
                "#version 100               \n"
                "#define GLUP_VERTEX_SHADER \n"
            );
        } else {
            sources.push_back(
#ifdef GEO_OS_APPLE		
		"#version 150 core          \n"
#else		
                "#version 130               \n"
#endif		
                "#define GLUP_VERTEX_SHADER \n"
            );
        }
    }
    

    void FullScreenEffectImpl::get_fragment_shader_preamble_pseudo_file(
        std::vector<GLSL::Source>& sources
    ) {
        if(ES_profile_) {
            sources.push_back(
                "#version 100                                    \n"
                "#define GLUP_FRAGMENT_SHADER                    \n"
                "precision mediump float;                        \n"
                "#ifdef GL_FRAGMENT_PRECISION_HIGH               \n"
                "   precision highp int;                         \n"
                "#endif                                          \n"
            );
        } else {
            sources.push_back(
#ifdef GEO_OS_APPLE		
		"#version 150 core                               \n"
#else		
                "#version 130                                    \n"
#endif		
                "#define GLUP_FRAGMENT_SHADER                    \n"
            );
        }
	sources.push_back(
	    "#ifdef GL_ES                            \n"
	    "#define glup_FragColor gl_FragColor     \n"
	    "#else                                   \n"
	    "out vec4 glup_FragColor;                \n"
	    "#endif	                             \n"
	);
    }
    
    void FullScreenEffectImpl::pre_render(index_t w, index_t h) {
        if(!initialized_) {
            initialize(w,h) ;
        }
        if(!OK_) {
            return;
        }
        if(w != width() || h != height()) {
	    std::cerr << "resize" << std::endl;
            resize(w,h) ;
        }
	glGetIntegerv(GL_VIEWPORT, viewport_);	
	
    }

    void FullScreenEffectImpl::post_render() {
        if(OK_) {
	    glViewport(viewport_[0], viewport_[1], viewport_[2], viewport_[3]);	    
        }
    }

    void FullScreenEffectImpl::update() {
    }
    
    void FullScreenEffectImpl::initialize(index_t width, index_t height) {
        width_       = width;
        height_      = height;
        initialized_ = true ;
	core_profile_ = CmdLine::get_arg("gfx:GL_profile") == "core";
#ifdef GEO_OS_EMSCRIPTEN
	ES_profile_ = true;
#endif


	static bool first_time = true;

	if(first_time) {
	    GLSL::register_GLSL_include_file(
		"fullscreen/current_profile/vertex_shader_preamble.h",
		static_cast<GLSL::PseudoFile>(
		    &GEO::FullScreenEffectImpl::get_vertex_shader_preamble_pseudo_file
		)
	    );
            
	    GLSL::register_GLSL_include_file(
		"fullscreen/current_profile/fragment_shader_preamble.h",
		static_cast<GLSL::PseudoFile>(
		    &GEO::FullScreenEffectImpl::get_fragment_shader_preamble_pseudo_file
		)
	    );

	    GLUP::register_embedded_shaders_fullscreen();	    
	    
	    first_time = false;
	}
	
        double supp_ver = GLSL::supported_language_version();
        double req_ver = required_GLSL_version();
        OK_ = (supp_ver >= req_ver);
        if(!OK_) {
            Logger::err("Renderer")
                << "This FullScreenEffect requires GLSL ver. >= "
                << req_ver
                << std::endl;
            Logger::err("Renderer")
                << "This GPU has GLSL ver. " << supp_ver
                << std::endl;
            return;
        }
    }

    void FullScreenEffectImpl::resize(index_t width, index_t height) {
        width_ = width;
        height_ = height;
    }

    void FullScreenEffectImpl::reset_alpha() {
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_TRUE);
	glClear(GL_COLOR_BUFFER_BIT);
	glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
    }
}

