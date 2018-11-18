/**
 * Copyright (c) 2018 Eric Bruneton
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holders nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

/*<h2>atmosphere/demo/webgl/precompute.cc</h2>

<p>This file precomputes the atmosphere textures and saves them to disk. It also
saves to disk the shaders necessary for the demo. For this a C++
<a href="../demo.h.html">Demo</a> instance is created (which precomputes the
textures and creates the shaders), its shaders and textures are read back using
the OpenGL API, and are saved to disk:
*/

#include <glad/glad.h>
#include <GL/freeglut.h>

#include <memory>
#include <fstream>

#include "atmosphere/demo/demo.h"
#include "atmosphere/constants.h"

using atmosphere::demo::Demo;

void SaveShader(const GLuint shader, const std::string& filename) {
  GLint sourceLength;
  glGetShaderiv(shader, GL_SHADER_SOURCE_LENGTH, &sourceLength);

  GLsizei actualLength;
  std::unique_ptr<GLchar[]> buffer(new GLchar[sourceLength]);
  glGetShaderSource(shader, sourceLength, &actualLength, buffer.get());

  std::ofstream output_stream(filename, std::ofstream::out);
  output_stream << std::string(buffer.get());
  output_stream.close();
}

void SaveTexture(const GLenum texture_unit, const GLenum texture_target,
    const int texture_size, const std::string& filename) {
  std::unique_ptr<float[]> pixels(new float[texture_size * 4]);
  glActiveTexture(texture_unit);
  glGetTexImage(texture_target, 0, GL_RGBA, GL_FLOAT, pixels.get());

  std::ofstream output_stream(
      filename, std::ofstream::out | std::ofstream::binary);
  output_stream.write((const char*) pixels.get(), texture_size * 16);
  output_stream.close();
}

int main(int argc, char** argv) {
  glutInitContextVersion(3, 3);
  glutInitContextProfile(GLUT_CORE_PROFILE);
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);

  std::unique_ptr<Demo> demo(new Demo(0, 0));
  demo->model().SetProgramUniforms(demo->program(), 0, 1, 2);

  const std::string output_dir(argv[1]);
  SaveShader(demo->model().shader(), output_dir + "atmosphere_shader.txt");
  SaveShader(demo->vertex_shader(), output_dir + "vertex_shader.txt");
  SaveShader(demo->fragment_shader(), output_dir + "fragment_shader.txt");
  SaveTexture(
      GL_TEXTURE0,
      GL_TEXTURE_2D,
      atmosphere::TRANSMITTANCE_TEXTURE_WIDTH *
          atmosphere::TRANSMITTANCE_TEXTURE_HEIGHT,
      output_dir + "transmittance.dat");
  SaveTexture(
      GL_TEXTURE1,
      GL_TEXTURE_3D,
      atmosphere::SCATTERING_TEXTURE_WIDTH *
          atmosphere::SCATTERING_TEXTURE_HEIGHT *
          atmosphere::SCATTERING_TEXTURE_DEPTH,
      output_dir + "scattering.dat");
  SaveTexture(
      GL_TEXTURE2,
      GL_TEXTURE_2D,
      atmosphere::IRRADIANCE_TEXTURE_WIDTH *
          atmosphere::IRRADIANCE_TEXTURE_HEIGHT,
      output_dir + "irradiance.dat");

  return 0;
}
