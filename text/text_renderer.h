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

#ifndef TEXT_TEXT_RENDERER_H_
#define TEXT_TEXT_RENDERER_H_

#include <glad/glad.h>

#include <string>

class TextRenderer {
 public:
  TextRenderer();
  TextRenderer(TextRenderer const&) = delete;
  TextRenderer(TextRenderer&&) = delete;
  ~TextRenderer();

  void SetColor(float r, float g, float b);
  void DrawText(const std::string& text, int left, int top);

 private:
  void SetupTexture();
  void SetupBuffers();
  void SetupProgram();
  void DrawChar(char c, int x, int y, int viewport_width, int viewport_height);

  GLuint font_texture_;
  GLuint char_vao_;
  GLuint char_vbo_;
  GLuint program_;
  GLfloat color_[3];
};

#endif  // TEXT_TEXT_RENDERER_H_
