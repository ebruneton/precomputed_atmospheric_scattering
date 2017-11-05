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

void TextRenderer::SetupTexture() {
  glGenTextures(1, &fontTexture_);

  // Avoid interfering with caller's assumptions
  glActiveTexture(GL_TEXTURE0);
  GLint oldTexture;
  glGetIntegerv(GL_TEXTURE_BINDING_2D, &oldTexture);

  glBindTexture(GL_TEXTURE_2D, fontTexture_);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, font.atlasWidth, font.atlasHeight,
               0, GL_RED, GL_UNSIGNED_BYTE, font.data.data());
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

  glBindTexture(GL_TEXTURE_2D, oldTexture);
}

void TextRenderer::SetupBuffers() {
  glGenVertexArrays(1, &vao_);
  glBindVertexArray(vao_);
  glGenBuffers(1, &vbo_);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_);
  const GLfloat vertices[] = {
    0, 0,
    1, 0,
    0, 1,
    1, 1,
  };
  glBufferData(GL_ARRAY_BUFFER, sizeof vertices, vertices, GL_STATIC_DRAW);
  constexpr GLuint kAttribIndex = 0;
  constexpr int kCoordsPerVertex = 2;
  glVertexAttribPointer(kAttribIndex, kCoordsPerVertex,
                        GL_FLOAT, false, 0, 0);
  glEnableVertexAttribArray(kAttribIndex);
  glBindVertexArray(0);
}

void TextRenderer::SetupProgram() {
  const auto vertexShader = glCreateShader(GL_VERTEX_SHADER);
  const char*const kVertexShaderSrc = R"(#version 330
uniform mat4 mvp;
uniform vec2 texCoordIn;
uniform vec2 charSizeInTexture;
layout(location = 0) in vec4 vertex;
out vec2 texCoord;
void main()
{
  gl_Position = mvp*vertex;
  texCoord = texCoordIn+charSizeInTexture*(vec2(vertex.x, 1-vertex.y));
}
)";
  glShaderSource(vertexShader, 1, &kVertexShaderSrc, nullptr);
  glCompileShader(vertexShader);
  const auto fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
  const char*const kFragmentShaderSrc = R"(#version 330
uniform vec3 textColor;
uniform sampler2D font;
in vec2 texCoord;
layout(location = 0) out vec4 color;
void main()
{
  color = vec4(textColor.rgb, 1)*texture(font, texCoord).rrrr;
}
)";
  glShaderSource(fragmentShader, 1, &kFragmentShaderSrc, nullptr);
  glCompileShader(fragmentShader);
  program_ = glCreateProgram();
  glAttachShader(program_, vertexShader);
  glAttachShader(program_, fragmentShader);
  glLinkProgram(program_);
  glDetachShader(program_, fragmentShader);
  glDeleteShader(fragmentShader);
  glDetachShader(program_, vertexShader);
  glDeleteShader(vertexShader);
}

void TextRenderer::DrawChar(char c, int x, int y,
                            int viewportWidth, int viewportHeight) {
  if (c < 0x20 || c > 0x7e) c = '?';
  glBindVertexArray(vao_);
  glUseProgram(program_);

  const GLfloat charW = font.charWidth, charH = font.charHeight;
  glUniform2f(glGetUniformLocation(program_, "charSizeInTexture"),
              charW / font.atlasWidth,
              charH / font.atlasHeight);
  glUniform2f(glGetUniformLocation(program_, "texCoordIn"),
              (c & 0xf) * charW / font.atlasWidth,
              ((c >> 4) - 2) * charH / font.atlasHeight);
  {
    const GLfloat dx = 2.f * x / viewportWidth  - 1,
                  dy = 2.f * y / viewportHeight - 1;
    const GLfloat sx = 2 * charW / viewportWidth,
                  sy = 2 * charH / viewportHeight;
    const GLfloat mvp[] = {sx, 0,  0, 0,
                           0,  sy, 0, 0,
                           0,  0,  1, 0,
                           dx, dy, 0, 1};
    glUniformMatrix4fv(glGetUniformLocation(program_, "mvp"),
                       1, false, mvp);
  }

  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, fontTexture_);
  glUniform1i(glGetUniformLocation(program_, "font"), 0);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

  glDisable(GL_BLEND);
  glBindVertexArray(0);
}

TextRenderer::TextRenderer() {
  SetupTexture();
  SetupBuffers();
  SetupProgram();

  SetColor(1, 1, 1);
}

TextRenderer::~TextRenderer() {
  glDeleteProgram(program_);
  glDeleteBuffers(1, &vbo_);
  glDeleteVertexArrays(1, &vao_);
  glDeleteTextures(1, &fontTexture_);
}

void TextRenderer::SetColor(float r, float g, float b) {
  color_[0] = r;
  color_[1] = g;
  color_[2] = b;
  GLint oldProgram;
  glGetIntegerv(GL_CURRENT_PROGRAM, &oldProgram);
  glUseProgram(program_);
  glUniform3fv(glGetUniformLocation(program_, "textColor"), 1, color_);
  glUseProgram(oldProgram);
}

void TextRenderer::DrawText(const char* text, int left, int top) {
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);

  // Avoid interfering with caller's assumptions
  GLint oldProgram;
  glGetIntegerv(GL_CURRENT_PROGRAM, &oldProgram);
  glActiveTexture(GL_TEXTURE0);
  GLint oldTexture;
  glGetIntegerv(GL_TEXTURE_BINDING_2D, &oldTexture);

  int x = left, y = viewport[3] - top - font.charHeight;
  for (int i = 0; text[i] != 0; ++i) {
    const char c = text[i];
    switch (c) {
    case ' ':
      break;
    case '\n':
      y -= font.charHeight + 1;
      x = left;
      continue;
    default:
      DrawChar(c, x, y, viewport[2], viewport[3]);
      break;
    }
    x += font.charWidth;
  }

  glBindTexture(GL_TEXTURE_2D, oldTexture);
  glUseProgram(oldProgram);
}
