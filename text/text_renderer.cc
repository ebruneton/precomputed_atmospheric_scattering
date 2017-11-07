/**
 * Copyright (c) 2017 Ruslan Kabatsayev
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

#include "text/text_renderer.h"

#include "font.inc"

namespace {

const char kVertexShader[] = R"(
    #version 330
    uniform mat4 clip_from_model;
    uniform mat3x2 texture_coord_from_model;
    layout(location = 0) in vec4 vertex;
    out vec2 texture_coord;
    void main() {
      gl_Position = clip_from_model * vertex;
      texture_coord = texture_coord_from_model * vertex.xyw;
    })";

const char kFragmentShader[] = R"(
    #version 330
    uniform vec3 text_color;
    uniform sampler2D font_sampler;
    in vec2 texture_coord;
    layout(location = 0) out vec4 color;
    void main() {
      color = vec4(text_color, 1) * texture(font_sampler, texture_coord).rrrr;
    })";

}  // anonymous namespace

void TextRenderer::SetupTexture() {
  glGenTextures(1, &font_texture_);

  // Avoid interfering with caller's assumptions.
  glActiveTexture(GL_TEXTURE0);
  GLint old_texture;
  glGetIntegerv(GL_TEXTURE_BINDING_2D, &old_texture);

  glBindTexture(GL_TEXTURE_2D, font_texture_);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, font.atlas_width, font.atlas_height,
               0, GL_RED, GL_UNSIGNED_BYTE, font.data.data());
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

  glBindTexture(GL_TEXTURE_2D, old_texture);
}

void TextRenderer::SetupBuffers() {
  glGenVertexArrays(1, &char_vao_);
  glBindVertexArray(char_vao_);
  glGenBuffers(1, &char_vbo_);
  glBindBuffer(GL_ARRAY_BUFFER, char_vbo_);
  const GLfloat vertices[] = {
    0, 0,
    1, 0,
    0, 1,
    1, 1,
  };
  glBufferData(GL_ARRAY_BUFFER, sizeof vertices, vertices, GL_STATIC_DRAW);
  constexpr GLuint kAttribIndex = 0;
  constexpr int kCoordsPerVertex = 2;
  glVertexAttribPointer(kAttribIndex, kCoordsPerVertex, GL_FLOAT, false, 0, 0);
  glEnableVertexAttribArray(kAttribIndex);
  glBindVertexArray(0);
}

void TextRenderer::SetupProgram() {
  const auto vertex_shader = glCreateShader(GL_VERTEX_SHADER);
  const char* const vertex_shader_source = kVertexShader;
  glShaderSource(vertex_shader, 1, &vertex_shader_source, nullptr);
  glCompileShader(vertex_shader);
  const auto fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
  const char* const fragment_shader_source = kFragmentShader;
  glShaderSource(fragment_shader, 1, &fragment_shader_source, nullptr);
  glCompileShader(fragment_shader);
  program_ = glCreateProgram();
  glAttachShader(program_, vertex_shader);
  glAttachShader(program_, fragment_shader);
  glLinkProgram(program_);
  glDetachShader(program_, fragment_shader);
  glDeleteShader(fragment_shader);
  glDetachShader(program_, vertex_shader);
  glDeleteShader(vertex_shader);
}

void TextRenderer::DrawChar(char c, int x, int y,
                            int viewport_width, int viewport_height) {
  if (c < 0x20 || c > 0x7e) {
    c = '?';
  }

  const GLfloat char_width = font.char_width;
  const GLfloat char_height = font.char_height;
  {
    const GLfloat scale_x = char_width / font.atlas_width;
    const GLfloat scale_y = -char_height / font.atlas_height;
    const int characters_per_line = font.atlas_width / font.char_width;
    const GLfloat translate_x =
        (c % characters_per_line) * char_width / font.atlas_width;
    const GLfloat translate_y =
        (c / characters_per_line - 1) * char_height / font.atlas_height;
    const GLfloat texture_coord_from_model[] = {
      scale_x, 0, translate_x,
      0, scale_y, translate_y
    };
    glUniformMatrix3x2fv(
        glGetUniformLocation(program_, "texture_coord_from_model"),
        1, true, texture_coord_from_model);
  }
  {
    const GLfloat scale_x = 2 * char_width / viewport_width;
    const GLfloat scale_y = 2 * char_height / viewport_height;
    const GLfloat translate_x = (2.f * x) / viewport_width - 1;
    const GLfloat translate_y = (2.f * y) / viewport_height - 1;
    const GLfloat clip_from_model[] = {
      scale_x, 0, 0, translate_x,
      0, scale_y, 0, translate_y,
      0, 0, 1, 0,
      0, 0, 0, 1
    };
    glUniformMatrix4fv(glGetUniformLocation(program_, "clip_from_model"),
                       1, true, clip_from_model);
  }

  glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
}

TextRenderer::TextRenderer() {
  SetupTexture();
  SetupBuffers();
  SetupProgram();

  SetColor(1, 1, 1);
}

TextRenderer::~TextRenderer() {
  glDeleteProgram(program_);
  glDeleteBuffers(1, &char_vbo_);
  glDeleteVertexArrays(1, &char_vao_);
  glDeleteTextures(1, &font_texture_);
}

void TextRenderer::SetColor(float r, float g, float b) {
  color_[0] = r;
  color_[1] = g;
  color_[2] = b;
  GLint old_program;
  glGetIntegerv(GL_CURRENT_PROGRAM, &old_program);
  glUseProgram(program_);
  glUniform3fv(glGetUniformLocation(program_, "text_color"), 1, color_);
  glUseProgram(old_program);
}

void TextRenderer::DrawText(const std::string& text, int left, int top) {
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);

  // Avoid interfering with caller's assumptions.
  GLint old_program;
  glGetIntegerv(GL_CURRENT_PROGRAM, &old_program);
  glActiveTexture(GL_TEXTURE0);
  GLint old_texture;
  glGetIntegerv(GL_TEXTURE_BINDING_2D, &old_texture);

  glBindVertexArray(char_vao_);
  glUseProgram(program_);

  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, font_texture_);
  glUniform1i(glGetUniformLocation(program_, "font_sampler"), 0);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  int x = left;
  int y = viewport[3] - top - font.char_height;
  for (char c : text) {
    switch (c) {
      case ' ':
        break;
      case '\n':
        x = left;
        y -= font.char_height + 1;
        continue;
      default:
        DrawChar(c, x, y, viewport[2], viewport[3]);
        break;
    }
    x += font.char_width;
  }

  glDisable(GL_BLEND);
  glBindVertexArray(0);

  glBindTexture(GL_TEXTURE_2D, old_texture);
  glUseProgram(old_program);
}
