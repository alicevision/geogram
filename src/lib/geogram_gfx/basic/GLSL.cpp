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

#include <geogram_gfx/basic/GLSL.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/command_line.h>
#include <cstdarg>
#include <cstdio>

namespace {

    using namespace GEO;
    
    /**
     * \brief Loads the content of an ASCII file in a buffer.
     * \details Memory ownership is transfered
     *   to the caller. Memory should be deallocated with
     *   delete[].
     * \param[in] filename the name of the file
     * \return a pointer to a buffer that contains the
     *   contents of the file. 
     */
     char* load_ASCII_file(const char* filename) {
        FILE* f = fopen(filename, "rt") ;
        if(!f) {
            Logger::err("GLSL")
                << "Could not open file: \'"
                << filename << "\'" << std::endl;
            return nil ;
        }
        /* 
         * An easy way of determining the length of a file:
         * Go to the end of the file, and ask where we are.
         */
        fseek(f, 0, SEEK_END) ;
        size_t size = size_t(ftell(f)) ;

        /* Let's go back to the beginning of the file */
        fseek(f, 0, SEEK_SET) ;
        
        char* result = new char[size+1] ;
        size_t read_size = fread(result, 1, size, f);
        if(read_size != size) {
            Logger::warn("GLSL")
                << "Could not read completely file \'"
                << filename << "\'" << std::endl;
        }
        result[size] = '\0' ;
        fclose(f) ;
        return result ;
    }

    /**
     * \brief Links a GLSL program and displays errors if any.
     * \details If errors where encountered, program is deleted
     *  and reset to zero.
     * \param[in,out] program the handle to the GLSL program
     */
    void link_program_and_check_status(GLuint& program) {
        glLinkProgram(program);
        GLint link_status;
        glGetProgramiv(program, GL_LINK_STATUS, &link_status);
        if(!link_status) {
            GLchar linker_message[4096];
            glGetProgramInfoLog(
                program, sizeof(linker_message), 0, linker_message
                );
            Logger::err("GLSL") << "linker status :"
                                << link_status << std::endl;
            Logger::err("GLSL") << "linker message:"
                                << linker_message << std::endl;
            glDeleteProgram(program);
            program = 0;
        }
        if(CmdLine::get_arg_bool("dbg:gfx")) {
            GLSL::introspect_program(program);
        }
    }
}

namespace GEO {

    /***********************************************************************/
    
    namespace GLSL {
            
        void initialize() {
        }
            
        void terminate() {
        }
            
        /*************************************************************/

        const char* GLSLCompileError::what() const GEO_NOEXCEPT {
            return "GLSL Compile Error";
        }
            
        /*************************************************************/

        

        /**
         * \brief Parses in a string what looks like a version number.
         * \param[in] version_string a const reference to the string
         * \return the parsed version number, as a double precision floating
         *  point number.
         */
        static double find_version_number(const std::string& version_string) {
            // The way the driver exposes the version of GLSL may differ,
            // in some drivers the number comes in first position, in some
            // others it comes in last position, therefore we take the first
            // word that contains a valid number.
            std::vector<std::string> version_words;
            String::split_string(
                version_string, ' ', version_words
            );
            double version = 0.0;
            for(index_t i=0; i<version_words.size(); ++i) {
                version = atof(version_words[i].c_str());
                if(version != 0.0) {
                    break;
                }
            }
            // Some drivers expose version 4.4 as 4.4 and some others
            // as 440 !!
            if(version > 100.0) {
                version /= 100.0;
            }
            return version;
        }

        // If some drivers do not implement glGetString(GLSL_VERSION),
        // then we can determine the GLSL version from the OpenGL version,
        // may be more reliable...
        //
        //GLSL Version      OpenGL Version
        //1.10              2.0
        //1.20              2.1
        //1.30              3.0
        //1.40              3.1
        //1.50              3.2
        //3.30              3.3
        //4.00              4.0
        //4.10              4.1
        //4.20              4.2
        //4.30              4.3
        //4.40              4.4
        //4.50              4.5
        
        static double GLSL_version_from_OpenGL_version() {
            const char* opengl_ver_str = (const char*)glGetString(GL_VERSION);
            if(opengl_ver_str == NULL) {
                Logger::warn("GLSL")
                    << "glGetString(GL_VERSION)"
                    << " did not answer, falling back to VanillaGL"
                    << std::endl;
               return 0.0;
            }
            double OpenGL_version = find_version_number(opengl_ver_str);
            
            Logger::out("GLSL")
                << "Determining GLSL version from OpenGL version"
                << std::endl;
            
            Logger::out("GLSL")
                << "OpenGL version = " << OpenGL_version
                << std::endl;

            double GLSL_version = 0.0;
            if(OpenGL_version >= 3.3) {
                GLSL_version = OpenGL_version;
            } else if(OpenGL_version == 2.0) {
                GLSL_version = 1.1;
            } else if(OpenGL_version == 2.1) {
                GLSL_version = 1.2;                
            } else if(OpenGL_version == 3.0) {
                GLSL_version = 1.3;                                
            } else if(OpenGL_version == 3.1) {
                GLSL_version = 1.4;
            } else if(OpenGL_version == 3.2) {
                GLSL_version = 1.5;                
            }

            if(GLSL_version == 0.0) {
                Logger::warn("GLSL") << "Could not determine GLSL version"
                                     << std::endl;
            } else {
                Logger::out("GLSL") << "GLSL version = "
                                     << GLSL_version
                                     << std::endl;
            }
            return GLSL_version;
        }

        
        double supported_language_version() {

            double GLSL_version = CmdLine::get_arg_double("gfx:GLSL_version");
            
            if(GLSL_version != 0.0) {
                Logger::out("GLSL") << "forced to version "
                                    << GLSL_version 
                                    << " (gfx:GLSL_version)" << std::endl;
                return GLSL_version;
            }

            // New OpenGL API: one should use glGetStringi()

            const char* shading_language_ver_str = nil;

#ifdef GEO_GL_150
            if(glGetStringi != nil) {
                shading_language_ver_str = (const char*)glGetStringi(
                    GL_SHADING_LANGUAGE_VERSION, 0
                );
            }
#endif            
            if(shading_language_ver_str == nil) {
                // Some buggy drivers do not implement glGetStringi(),
                // so I try also glGetString() (without the "i")
                shading_language_ver_str =
                    (const char*)glGetString(GL_SHADING_LANGUAGE_VERSION);
            }

            // If the driver does not implement glGetString neither
            // glGetStringi with GL_SHADING_LANGUAGE_VERSION, then try
            // to deduce it from OpenGL version.
            if(shading_language_ver_str == nil) {
                return GLSL_version_from_OpenGL_version();
            }
           
            const char* vendor = (const char*)glGetString(GL_VENDOR);
            
            Logger::out("GLSL") << "vendor = " << vendor << std::endl;
            Logger::out("GLSL") << "version string = "
                                << shading_language_ver_str << std::endl;


            GLSL_version = find_version_number(shading_language_ver_str);
            Logger::out("GLSL") << "version = " << GLSL_version
                                << std::endl;
            
            if(!CmdLine::get_arg_bool("gfx:GLSL")) {
                Logger::out("GLSL")
                    << "OpenGL shaders deactivated (gfx:GLSL=false)"
                    << std::endl;
                
                GLSL_version = 0.0;
            }
            return GLSL_version;
        }
            

        
        GLuint compile_shader(
            GLenum target, const char** sources, index_t nb_sources
        ) {
            GLuint s_handle = glCreateShader(target);
            if(s_handle == 0) {
                Logger::err("GLSL") << "Could not create shader for target"
                                    << std::endl;
                switch(target) {
                case GL_VERTEX_SHADER:
                    Logger::err("GLSL") << " (target = GL_VERTEX_SHADER)"
                                        << std::endl;
                    break;
                case GL_FRAGMENT_SHADER:
                    Logger::err("GLSL")
                        << " (target = GL_FRAGMENT_SHADER)"
                        << std::endl;
                    break;
#ifdef GEO_GL_150
                case GL_COMPUTE_SHADER:
                    Logger::err("GLSL") << " (target = GL_COMPUTE_SHADER)"
                                        << std::endl;
                    break;
                case GL_TESS_CONTROL_SHADER:
                    Logger::err("GLSL") << " (target = GL_TESS_CONTROL_SHADER)"
                                        << std::endl;
                    break;
                case GL_TESS_EVALUATION_SHADER:
                    Logger::err("GLSL")
                        << " (target = GL_TESS_EVALUATION_SHADER)"
                        << std::endl;
                    break;
                case GL_GEOMETRY_SHADER:
                    Logger::err("GLSL")
                        << " (target = GL_GEOMETRY_SHADER)"
                        << std::endl;
                    break;
#endif                    
                default:
                    Logger::err("GLSL")
                        << " (unknown target)"
                        << std::endl;
                    break;
                }
                throw GLSL::GLSLCompileError();
            }
            glShaderSource(s_handle, (GLsizei)nb_sources, sources, 0);
            glCompileShader(s_handle);
            GLint compile_status;
            glGetShaderiv(s_handle, GL_COMPILE_STATUS, &compile_status);
            if(!compile_status) {
                GLchar compiler_message[4096];
                glGetShaderInfoLog(
                    s_handle, sizeof(compiler_message), 0, compiler_message
                );
                Logger::err("GLSL")
                    << "compiler status :"
                    << compile_status << std::endl;
                Logger::err("GLSL")
                    << "compiler message:"
                    << compiler_message << std::endl;

                Logger::out("GLSL") << "Erroneous program source:"
                                    << std::endl;
                for(index_t i=0; i<nb_sources; ++i) {
                    Logger::out("GLSL") << sources[i] << std::endl;
                }

                
                glDeleteShader(s_handle);
                s_handle = 0;
                throw GLSLCompileError();
            }
            return s_handle;
        }

        GLuint compile_shader(
            GLenum target,
            const char* source1,
            const char* source2,
            const char* source3,
            const char* source4,
            const char* source5,
            const char* source6,
            const char* source7,
            const char* source8,
            const char* source9,
            const char* source10,
            const char* source11,
            const char* source12,
            const char* source13,
            const char* source14,
            const char* source15,
            const char* source16,
            const char* source17,
            const char* source18,
            const char* source19,
            const char* source20
        ) {
            vector<const char*> sources;
            geo_assert(source1 != nil);
            if(source1 != nil) {
                sources.push_back(source1);
            }
            if(source2 != nil) {
                sources.push_back(source2);
            }
            if(source3 != nil) {
                sources.push_back(source3);
            }
            if(source4 != nil) {
                sources.push_back(source4);
            }
            if(source5 != nil) {
                sources.push_back(source5);
            }
            if(source6 != nil) {
                sources.push_back(source6);
            }
            if(source7 != nil) {
                sources.push_back(source7);
            }
            if(source8 != nil) {
                sources.push_back(source8);
            }
            if(source9 != nil) {
                sources.push_back(source9);
            }
            if(source10 != nil) {
                sources.push_back(source10);
            }
            if(source11 != nil) {
                sources.push_back(source11);
            }
            if(source12 != nil) {
                sources.push_back(source12);
            }
            if(source13 != nil) {
                sources.push_back(source13);
            }
            if(source14 != nil) {
                sources.push_back(source14);
            }
            if(source15 != nil) {
                sources.push_back(source15);
            }
            if(source16 != nil) {
                sources.push_back(source16);
            }
            if(source17 != nil) {
                sources.push_back(source17);
            }
            if(source18 != nil) {
                sources.push_back(source18);
            }
            if(source19 != nil) {
                sources.push_back(source19);
            }
            if(source20 != nil) {
                sources.push_back(source20);
            }

            if(CmdLine::get_arg_bool("dbg:gfx")) {
                std::ofstream out("last_shader.glsl");
                
                for(index_t i=0; i<sources.size(); ++i) {
                    out << sources[i];
                }
                
                Logger::out("GLSL") << "===== Shader source ===="
                                    << std::endl;
                
                for(index_t i=0; i<sources.size(); ++i) {
                    Logger::out("GLSL") << sources[i];
                }
            }
            
            return compile_shader(target, &sources[0], sources.size());
        }


        void link_program(GLuint program) {
            link_program_and_check_status(program);
            if(program == 0) {
                throw GLSL::GLSLCompileError();                
            }
        }

        GLuint create_program_from_shaders_no_link(GLuint shader1, ...) {
            va_list args;            
            GLuint program = glCreateProgram();
            va_start(args,shader1);
            GLuint shader = shader1;
            while(shader != 0) {
                glAttachShader(program, shader);
                shader = va_arg(args, GLuint);
            }
            va_end(args);
            return program;
        }
        
        GLuint create_program_from_shaders(GLuint shader1, ...) {
            va_list args;            
            GLuint program = glCreateProgram();
            va_start(args,shader1);
            GLuint shader = shader1;
            while(shader != 0) {
                glAttachShader(program, shader);
                shader = va_arg(args, GLuint);
            }
            va_end(args);
            link_program_and_check_status(program);
            return program;
        }

        /*****************************************************************/

        GLuint create_program_from_string_no_link(
            const char* string_in, bool copy_string
        ) {
            GLuint program = glCreateProgram();
            
            // string will be temporarily modified (to insert '\0' markers)
            // but will be restored to its original state right after.
            char* string = const_cast<char*>(string_in);
            if(copy_string) {
                string = strdup(string_in);
            }
            
            char* src = string;
            
            bool err_flag = false;
            
            for(;;) {
                char* begin = strstr(src, "#BEGIN(");
                char* end = strstr(src, "#END(");
                
                if(begin == nil && end == nil) {
                    break;
                }
                
                if(begin == nil) {
                    Logger::err("GLSL") << "missing #BEGIN() statement"
                                        << std::endl;
                    err_flag = true;
                    break;
                }
                
                if(end == nil) {
                    Logger::err("GLSL") << "missing #END() statement"
                                        << std::endl;
                    err_flag = true;
                    break;
                }

                    
                if(begin > end) {
                    Logger::err("GLSL") << "#END() before #BEGIN()"
                                        << std::endl;
                    err_flag = true;
                    break;
                }

                char* begin_opening_brace = begin + strlen("#BEGIN");
                char* end_opening_brace = end + strlen("#END");
                
                char* begin_closing_brace = strchr(begin_opening_brace,')');
                char* end_closing_brace = strchr(end_opening_brace, ')');
                
                if(begin_closing_brace == nil) {
                    Logger::err("GLSL") << "#BEGIN: missing closing brace"
                                        << std::endl;
                    err_flag = true;
                    break;
                }

                if(end_closing_brace == nil) {
                    Logger::err("GLSL") << "#END: missing closing brace"
                                        << std::endl;
                    err_flag = true;
                    break;
                }
                    
                std::string begin_kw(
                    begin_opening_brace+1,
                    size_t((begin_closing_brace - begin_opening_brace) - 1)
                    );
                std::string end_kw(
                    end_opening_brace+1,
                    size_t((end_closing_brace - end_opening_brace) - 1)
                    );
                if(end_kw != begin_kw) {
                    Logger::err("GLSL")
                        << "Mismatch: #BEGIN(" << begin_kw
                        << ") / #END(" << end_kw << ")"
                        << std::endl;
                    err_flag = true;
                    break;
                }
                
                // Replace '#END(...)' with string end marker
                *end = '\0';
                
                GLenum shader_type = GLenum(0);
                if(begin_kw == "GL_VERTEX_SHADER") {
                    shader_type = GL_VERTEX_SHADER;
                } else if(begin_kw == "GL_FRAGMENT_SHADER") {
                    shader_type = GL_FRAGMENT_SHADER;
                }
#ifdef GEO_GL_150
                  else if(begin_kw == "GL_GEOMETRY_SHADER") {
                    shader_type = GL_GEOMETRY_SHADER;
                } else if(begin_kw == "GL_TESS_CONTROL_SHADER") {
                    shader_type = GL_TESS_CONTROL_SHADER;
                } else if(begin_kw == "GL_TESS_EVALUATION_SHADER") {
                    shader_type = GL_TESS_EVALUATION_SHADER;
                }
#endif
                  else {
                    Logger::err("GLSL") << begin_kw
                                        << ": No such shader type"
                                        << std::endl;
                    err_flag = true;
                    break;
                }
                
                src = begin_closing_brace+1;
                GLuint shader = 0;
                try {
                    shader = compile_shader(shader_type, src, 0);
                } catch(...) {
                    err_flag = true;
                    break;
                }
                glAttachShader(program, shader);
                
                // Restore '#END(...)' statement
                // ('#' was replaced with string end marker).
                *end = '#';
                
                src = end_closing_brace + 1;
            }
            
            if(copy_string) {
                free(string);
            }
            
            if(err_flag) {
                glDeleteProgram(program);
                return 0;
            }
            return program;
        }

        /*****************************************************************/

        GLuint create_program_from_file_no_link(const std::string& filename) {
            char* buffer = load_ASCII_file(filename.c_str());
            if(buffer == nil) {
                return 0;
            }
            GLuint result = 0;
            try {
                // last argument to false:
                //  no need to copy the buffer, we know it
                // is not a string litteral.
                result = create_program_from_string_no_link(buffer,false);
            } catch(...) {
                delete[] buffer;
                throw;
            }
            return result;
        }

        /*****************************************************************/

        GLint GEOGRAM_GFX_API get_uniform_variable_offset(
            GLuint program, const char* varname
        ) {
#ifndef GEO_GL_150
            geo_argused(program);
            geo_argused(varname);
            return -1;
#else
            GLuint index = GL_INVALID_INDEX;
            glGetUniformIndices(program, 1, &varname, &index);
            if(index == GL_INVALID_INDEX) {
                Logger::err("GLUP")
                    << varname 
                    << ":did not find uniform state variable"
                    << std::endl;
                throw GLSL::GLSLCompileError();
            }
            geo_assert(index != GL_INVALID_INDEX);
            GLint offset = -1;
            glGetActiveUniformsiv(
                program, 1, &index, GL_UNIFORM_OFFSET, &offset
            );
            geo_assert(offset != -1);
            return offset;
#endif            
        }


        void introspect_program(GLuint program) {
            Logger::out("GLSL") << "Program " << program << " introspection:"
                                << std::endl;
            if(!glIsProgram(program)) {
                Logger::out("GLSL") << "  not a program !"
                                    << std::endl;
                return;
            }

            {
                GLint link_status;
                glGetProgramiv(program, GL_LINK_STATUS, &link_status);
                Logger::out("GLSL") << "  link status=" << link_status
                                    << std::endl;
            }

            {
                GLint active_attributes;
                glGetProgramiv(
                    program, GL_ACTIVE_ATTRIBUTES, &active_attributes
                );
                Logger::out("GLSL")
                    << "  active attributes=" << active_attributes
                    << std::endl;
                for(GLuint i=0; i<GLuint(active_attributes); ++i) {
                    GLsizei length;
                    GLint size;
                    GLenum type;
                    GLchar name[1024];
                    glGetActiveAttrib(
                        program, i, GLsizei(1024), &length, &size, &type, name
                    );
                    Logger::out("GLSL") << "    Attribute " << i << " : "
                                        << name
                                        << std::endl;
                }
            }

            {
                GLint active_uniforms;
                glGetProgramiv(program, GL_ACTIVE_UNIFORMS, &active_uniforms);
                Logger::out("GLSL") << "  active uniforms=" << active_uniforms
                                    << std::endl;
                for(GLuint i=0; i<GLuint(active_uniforms); ++i) {
                    GLsizei length;
                    GLint size;
                    GLenum type;
                    GLchar name[1024];
                    glGetActiveUniform(
                        program, i, GLsizei(1024), &length, &size, &type, name
                    );
                    Logger::out("GLSL") << "    Uniform " << i << " : "
                                        << name
                                        << std::endl;
                }
            }

#ifdef GEO_GL_150            
            {
                GLint active_uniform_blocks;
                glGetProgramiv(
                    program, GL_ACTIVE_UNIFORM_BLOCKS, &active_uniform_blocks
                );
                Logger::out("GLSL") << "  active uniform blocks="
                                    << active_uniform_blocks
                                    << std::endl;
            }
#endif
            
        }
        
    }

}

