/**
 * Copyright (c) 2017 Eric Bruneton
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

#include <fstream>
#include <iostream>
#include <regex>
#include <string>

namespace {

std::string GenerateHtml(const std::string& source,
    const std::string& html_template) {

  std::string src = source;
  // Put the copyright notice in a comment.
  src = std::regex_replace(src, std::regex("\\/\\*\\*"), "<!--",
      std::regex_constants::format_first_only);
  src = std::regex_replace(src, std::regex(" \\*\\/"), "-->",
      std::regex_constants::format_first_only);

  // Put the code parts in <pre> blocks.
  const std::string kPreBegin = "<pre class=\"prettyprint\">";
  const std::string kPreEnd = "</pre>";
  src = std::regex_replace(src, std::regex("\\n\\/\\*"), kPreEnd);
  src = std::regex_replace(src, std::regex("\\*\\/\\n"), kPreBegin);
  if (src.find(kPreBegin) == std::string::npos) {
    src = std::regex_replace(src, std::regex("-->\\n"), "-->\n" + kPreBegin,
        std::regex_constants::format_first_only);
  }
  src += kPreEnd;

  // Escape the < and > characters in <pre> blocks.
  std::stringstream body;
  const std::string kBodyPlaceHolder = "BODY";
  size_t start_pos = html_template.find(kBodyPlaceHolder);
  body << html_template.substr(0, start_pos);
  body << "<body>\n";
  bool in_pre_block = false;
  for (unsigned int i = 0; i < src.length(); ++i) {
    if (src[i] == '>') {
      if (in_pre_block) {
        body << "&gt;";
      } else {
        if (i + 1 >= kPreBegin.length() &&
            src.substr(i + 1 - kPreBegin.length(), kPreBegin.length()) ==
                kPreBegin) {
          in_pre_block = true;
        }
        body << ">";
      }
    } else if (src[i] == '<') {
      if (in_pre_block &&
          i + kPreEnd.length() <= src.length() &&
          src.substr(i, kPreEnd.length()) == kPreEnd) {
        in_pre_block = false;
      }
      body << (in_pre_block ? "&lt;" : "<");
    } else {
      body << src[i];
    }
  }
  body << "\n</body>";
  body << html_template.substr(start_pos + kBodyPlaceHolder.length());
  return body.str();
}

}  // anonymous namespace

int main(int argc, char** argv) {
  if (argc != 4) {
    std::cout << "Usage: " << argv[0]
              << " <source_file> <template_file> <output_file>" << std::endl;
    return -1;
  }

  std::ifstream source_stream(argv[1]);
  std::string source(
      (std::istreambuf_iterator<char>(source_stream)),
      std::istreambuf_iterator<char>());

  std::ifstream html_template_stream(argv[2]);
  std::string html_template(
      (std::istreambuf_iterator<char>(html_template_stream)),
      std::istreambuf_iterator<char>());

  std::ofstream output_file(argv[3]);
  output_file << GenerateHtml(source, html_template);
  output_file.close();

  return 0;
}
