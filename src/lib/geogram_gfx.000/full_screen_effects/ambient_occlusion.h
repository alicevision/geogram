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

#ifndef H_OGF_RENDERER_CONTEXT_AMBIENT_OCCLUSION_H
#define H_OGF_RENDERER_CONTEXT_AMBIENT_OCCLUSION_H

#include <geogram_gfx/basic/common.h>
#include <geogram_gfx/full_screen_effects/full_screen_effect.h>
#include <geogram_gfx/basic/frame_buffer_object.h>

/**
 * \file geogram_gfx/full_screen_effects/ambient_occlusion.h
 * \brief Implementation of AmbientOcclusion full screen effect.
 */

namespace GEO {

    /**
     * \brief Implementation of AmbientOcclusion full screen effect.
     */
    class GEOGRAM_GFX_API AmbientOcclusionImpl : public FullScreenEffectImpl {
    public:

        /**
         * \brief AmbientOcclusionImpl constructor.
         */
        AmbientOcclusionImpl();

        /**
         * \brief AmbientOcclusionImpl destructor.
         */
        ~AmbientOcclusionImpl();

        /**
         * \copydoc FullScreenEffectImpl::required_GLSL_version()
         */
        virtual double required_GLSL_version() const;
        
        /**
         * \copydoc FullScreenEffectImpl::pre_render()
         */
        virtual void pre_render(index_t w, index_t h);

        /**
         * \copydoc FullScreenEffectImpl::post_render()
         */
        virtual void post_render();

        /**
         * \copydoc FullScreenEffectImpl::update()
         */
        virtual void update();

        /**
         * \brief Gets the lightness.
         * \return the lightness, as an integer.
         * \details A value of 10 (default) corresponds to average lightness.
         */
        index_t get_lightness() const {
            return lightness_;
        }

        /**
         * \brief Sets the lightness.
         * \param[in] x the lightness, as an integer.
         * \details A value of 10 (default) corresponds to average lightness.
         */
        void set_lightness(index_t x) {
            lightness_ = x;
        }

        /**
         * \brief Gets the contrast.
         * \return the contrast, as an integer.
         * \details A value of 10 (default) corresponds to average contrast.
         */
        index_t get_contrast() const {
            return contrast_;
        }

        /**
         * \brief Sets the contrast.
         * \details A value of 10 (default) corresponds to average contrast.
         */
        void set_contrast(index_t x) {
            contrast_ = x;
        }

        /**
         * \brief Gets the size of the blurring kernel.
         * \return the size of the blurring kernel, in pixels.
         */
        index_t get_blur_width() const {
            return blur_width_;
        }

        /**
         * \brief Sets the size of the blurring kernel.
         * \param[in] x the size of the blurring kernel, in pixels.
         */
        void set_blur_width(index_t x) {
            blur_width_ = x;
        }

        /**
         * \brief Gets the number of directions.
         * \details These directions are used to compute the
         *  visibility integrals (the higher, the better),
         *  typical value is 5 to 10.
         * \return x the number of directions.
         */
        index_t get_nb_directions() const {
            return nb_directions_;
        }

        /**
         * \brief Sets the number of directions.
         * \details These directions are used to compute the
         *  visibility integrals (the higher, the better),
         *  typical value is 5 to 10.
         * \param[in] x the number of directions.
         */
        void set_nb_directions(index_t x) {
            nb_directions_ = x;
        }

    protected:
        /**
         * \copydoc FullScreenEffectImpl::initialize()
         */
        virtual void initialize(index_t w, index_t h);
        
        /**
         * \copydoc FullScreenEffectImpl::resize()
         */
        virtual void resize(index_t w, index_t h);

        /**
         * \brief Creates a texture with random values
         *  in it.
         * \details Used to randomly rotate the direction
         *   sampling pattern.
         */
        void create_random_tex();

        /**
         * \brief Copies the depth buffer into a texture.
         * \details Creates a copy of the depth buffer,
         *   with the same resolution as the depth buffer
         *   in depth_tex_.
         */
        void copy_depth_buffer_to_texture();

        /**
         * \brief Copies the depth buffer into a texture.
         * \details Creates a scaled copy of the depth buffer,
         *   with 1/REDUCTION resolution in reduced_depth_FBO_.
         */
        void copy_depth_texture_to_FBO();

        /**
         * \brief Displays the final result.
         * \details It composites the final result with the
         *   previously rendered frame.
         */
        void display_final_texture();

        /**
         * \brief Computes the ambient occlusion in the blur_1_
         *   FrameBufferObject, from the depth textures.
         */
        void compute_SSAO();

        /**
         * \brief Applies a Gaussian blur to the (raw) ambient occlusion 
         *  computed by apply_shader().
         * \details The input and the result are both in blur_1_. 
         *  The function uses blur_2_ as a work variable. It does two passes 
         *  of 1D blurring (horizontal and vertical, one from blur_1_ to 
         *  blur_2_ and the other from blur_2_ to blur_1_.
         */
        void blur();

        /**
         * \brief Gets the inverse of the projection transform.
         * \details It is used by the SSAO shader, to inverse-map 
         *   screen-space coordinates into world space.
         */
        void get_proj_inv();
        
    private:
        index_t lightness_;
        index_t contrast_;
        index_t blur_width_;
        index_t nb_directions_;

        GLint proj_inv_loc_;
        GLfloat proj_inv_[16];
        
        GLuint depth_tex_; // 24 bits depth
        FrameBufferObject reduced_depth_FBO_;
                   // transformed into 32 bits floating points
        
        GLuint random_tex_;
        GLuint SSAO_program_;

        // Blur
        FrameBufferObject blur_1_;
        FrameBufferObject blur_2_;
        GLuint blur_program_;
    };

    /**
     * \brief An automatic reference-counted pointer to an AmbientOcclusionImpl.
     */
    typedef SmartPointer<AmbientOcclusionImpl> AmbientOcclusionImpl_var ;
    
}

#endif
