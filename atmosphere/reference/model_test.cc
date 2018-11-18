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

/*<h2>atmosphere/reference/model_test.cc</h2>

<p>This file provides tests that compare the results of our GPU atmosphere model
with reference images computed with the same code, but executed on CPU in full
spectral mode and with double precision floats. The goal is to make sure that
the approximations made in the GPU version (half precision floats, single Mie
scattering extrapolated from a single wavelength, and luminance values computed
from 3 or 15 wavelengths instead of 47 wavelengths on CPU) do not significantly
reduce the image quality.

<p>The code is organized as follows:
<ul>
<li><a href="#constants">Constants</a></li>
<li><a href="#fixture">Test fixture</a>
<ul>
<li><a href="#setup">Setup methods</a>
<li><a href="#rendering">Rendering methods</a>
<li><a href="#comparison">Comparison methods</a>
</ul>
</li>
<li><a href="#cases">Test cases</a></li>
</ul>

<p>The test results can be seen <a href="test_report.html">here</a>.

<h3 id="constants">Constants</h3>

<p>We start by including the files we need. Since we want to render images with
the GPU and CPU versions of our algorithm, we need to include the GPU and CPU
model definitions:
*/

#include "atmosphere/reference/model.h"

#include <glad/glad.h>
#include <GL/freeglut.h>

#include <array>
#include <fstream>
#include <memory>

#include "atmosphere/model.h"
#include "atmosphere/reference/definitions.h"
#include "minpng/minpng.h"
#include "test/test_case.h"
#include "util/progress_bar.h"

/*
<p>Our test scene is a sphere on a purely spherical planet. Its position and
size are specified by the following constants (note that we use a large sphere
so that it can produce visible light shafts, in order to test them):
*/

namespace atmosphere {
namespace reference {

namespace {

constexpr Length kSphereRadius = 1.0 * km;
constexpr Position kSphereCenter = Position(0.0 * km, 0.0 * km, kSphereRadius);

/*
<p>Our tests use length values expressed in kilometers and, for the tests based
on radiance values (as opposed to luminance values), make use of the following
3 wavelengths:
*/

constexpr Length kLengthUnit = 1.0 * km;
constexpr Wavelength kLambdaR = atmosphere::Model::kLambdaR * nm;
constexpr Wavelength kLambdaG = atmosphere::Model::kLambdaG * nm;
constexpr Wavelength kLambdaB = atmosphere::Model::kLambdaB * nm;

/*
<p>The test scene is rendered on GPU by the following shaders. The vertex shader
simply renders a full screen quad, and outputs the view ray direction in model
space:
*/

const char kVertexShader[] = R"(
    #version 330
    uniform mat3 model_from_clip;
    layout(location = 0) in vec4 vertex;
    out vec3 view_ray;
    void main() {
      view_ray = model_from_clip * vertex.xyw;
      gl_Position = vertex;
    })";

/*
<p>The fragment shader computes the radiance (or luminance, if the USE_LUMINANCE
preprocessor macro is defined) corresponding to this view ray and uses a simple
tone mapping function to convert it to a final color. This shader takes as input
some uniforms describing the camera and the scene:
*/

const char kFragmentShader[] = R"(
    #define OUT(x) out x
    uniform vec3 camera_;
    uniform float exposure_;
    uniform vec3 earth_center_;
    uniform vec3 sun_direction_;
    uniform vec2 sun_size_;
    uniform vec3 ground_albedo_;
    uniform vec3 sphere_albedo_;
    in vec3 view_ray;
    layout(location = 0) out vec3 color;

    #ifdef USE_LUMINANCE
    #define GetSolarRadiance GetSolarLuminance
    #define GetSkyRadiance GetSkyLuminance
    #define GetSkyRadianceToPoint GetSkyLuminanceToPoint
    #define GetSunAndSkyIrradiance GetSunAndSkyIlluminance
    #endif

    vec3 GetSolarRadiance();
    vec3 GetSkyRadiance(vec3 camera, vec3 view_ray, float shadow_length,
        vec3 sun_direction, out vec3 transmittance);
    vec3 GetSkyRadianceToPoint(vec3 camera, vec3 point, float shadow_length,
        vec3 sun_direction, out vec3 transmittance);
    vec3 GetSunAndSkyIrradiance(
        vec3 p, vec3 normal, vec3 sun_direction, out vec3 sky_irradiance);
    vec3 GetViewRayRadiance(vec3 view_ray, vec3 view_ray_diff);

    void main() {
      color = GetViewRayRadiance(view_ray, dFdx(view_ray) + dFdy(view_ray));
      color = pow(vec3(1.0) - exp(-color * exposure_), vec3(1.0 / 2.2));
    })";

/*
<p>The main function <code>GetViewRayRadiance</code> is defined in
<a href="model_test.glsl.html">model_test.glsl</a>, which must be appended to
the above code (together with
<a href="../definitions.glsl.html">definitions.glsl</a>) to get a complete
shader. These GLSL files are provided as C++ string literals by the generated
files included below:
*/

#include "atmosphere/definitions.glsl.inc"
#include "atmosphere/reference/model_test.glsl.inc"

/*
<p>Each test case produces two images, using two different methods, and checks
that the difference between the two is small enough. These images are stored on
disk with the following function, which is a simple wrapper around the
<a href="https://github.com/jrmuizel/minpng">minpng</a> library:
*/

typedef std::unique_ptr<unsigned int[]> Image;

const char kOutputDir[] = "output/Doc/atmosphere/reference/";
constexpr unsigned int kWidth = 640;
constexpr unsigned int kHeight = 360;

void WritePngArgb(const std::string& name, void* pixels) {
  write_png((std::string(kOutputDir) + name).c_str(), pixels, kWidth, kHeight);
}

}  // anonymous namespace

/*
<h3 id="fixture">Test fixture</h3>

<p>The test fixture provides a shared class and shared methods for all the
test cases. It extends the <code>TestCase</code> class provided by the <a href=
"https://github.com/ebruneton/dimensional_types">dimensional_types</a> library:
*/

using std::max;
using std::min;

class ModelTest : public dimensional::TestCase {
 public:
  template<typename T>
  ModelTest(const std::string& name, T test)
      : TestCase("ModelTest " + name, static_cast<Test>(test)), name_(name) {}

/*
<h4 id="setup">Setup methods</h4>

<p>The <code>SetUp</code> method is called before each test case. We put here
the initialization code which must be executed before any test case, i.e. the
initialization of the atmosphere parameters (all our tests use the same
atmosphere parameters) and of the constant scene parameters (earth center, sun
size and radiance, surface albedos):
*/

  void SetUp() override {
    // Values from "Reference Solar Spectral Irradiance: ASTM G-173", ETR column
    // (see http://rredc.nrel.gov/solar/spectra/am1.5/ASTMG173/ASTMG173.html),
    // summed and averaged in each bin (e.g. the value for 360nm is the average
    // of the ASTM G-173 values for all wavelengths between 360 and 370nm).
    // Values in W.m^-2.
    constexpr int kLambdaMin = 360;
    constexpr int kLambdaMax = 830;
    constexpr double kSolarIrradiance[48] = {
      1.11776, 1.14259, 1.01249, 1.14716, 1.72765, 1.73054, 1.6887, 1.61253,
      1.91198, 2.03474, 2.02042, 2.02212, 1.93377, 1.95809, 1.91686, 1.8298,
      1.8685, 1.8931, 1.85149, 1.8504, 1.8341, 1.8345, 1.8147, 1.78158, 1.7533,
      1.6965, 1.68194, 1.64654, 1.6048, 1.52143, 1.55622, 1.5113, 1.474, 1.4482,
      1.41018, 1.36775, 1.34188, 1.31429, 1.28303, 1.26758, 1.2367, 1.2082,
      1.18737, 1.14683, 1.12362, 1.1058, 1.07124, 1.04992
    };
    constexpr ScatteringCoefficient kRayleigh = 1.24062e-6 / m;
    constexpr Length kRayleighScaleHeight = 8000.0 * m;
    constexpr Length kMieScaleHeight = 1200.0 * m;
    constexpr double kMieAngstromAlpha = 0.0;
    constexpr double kMieAngstromBeta = 5.328e-3;
    constexpr double kMieSingleScatteringAlbedo = 0.9;
    constexpr double kMiePhaseFunctionG = 0.8;
    // Values from http://www.iup.uni-bremen.de/gruppen/molspec/databases/
    // referencespectra/o3spectra2011/index.html for 233K, summed and averaged
    // in each bin (e.g. the value for 360nm is the average of the original
    // values for all wavelengths between 360 and 370nm). Values in m^2.
    constexpr double kOzoneCrossSection[48] = {
      1.18e-27, 2.182e-28, 2.818e-28, 6.636e-28, 1.527e-27, 2.763e-27, 5.52e-27,
      8.451e-27, 1.582e-26, 2.316e-26, 3.669e-26, 4.924e-26, 7.752e-26,
      9.016e-26, 1.48e-25, 1.602e-25, 2.139e-25, 2.755e-25, 3.091e-25, 3.5e-25,
      4.266e-25, 4.672e-25, 4.398e-25, 4.701e-25, 5.019e-25, 4.305e-25,
      3.74e-25, 3.215e-25, 2.662e-25, 2.238e-25, 1.852e-25, 1.473e-25,
      1.209e-25, 9.423e-26, 7.455e-26, 6.566e-26, 5.105e-26, 4.15e-26,
      4.228e-26, 3.237e-26, 2.451e-26, 2.801e-26, 2.534e-26, 1.624e-26,
      1.465e-26, 2.078e-26, 1.383e-26, 7.105e-27
    };
    // From https://en.wikipedia.org/wiki/Dobson_unit, in molecules.m^-2.
    constexpr dimensional::Scalar<-2, 0, 0, 0, 0> kDobsonUnit = 2.687e20 / m2;
    // Maximum number density of ozone molecules, in m^-3 (computed so at to get
    // 300 Dobson units of ozone - for this we divide 300 DU by the integral of
    // the ozone density profile defined below, which is equal to 15km).
    constexpr NumberDensity kMaxOzoneNumberDensity =
        300.0 * kDobsonUnit / (15.0 * km);

    std::vector<SpectralIrradiance> solar_irradiance;
    std::vector<ScatteringCoefficient> rayleigh_scattering;
    std::vector<ScatteringCoefficient> mie_scattering;
    std::vector<ScatteringCoefficient> mie_extinction;
    std::vector<ScatteringCoefficient> absorption_extinction;
    for (int l = kLambdaMin; l <= kLambdaMax; l += 10) {
      double lambda = static_cast<double>(l) * 1e-3;  // micro-meters
      SpectralIrradiance solar = kSolarIrradiance[(l - kLambdaMin) / 10] *
          watt_per_square_meter_per_nm;
      ScatteringCoefficient rayleigh = kRayleigh * pow(lambda, -4);
      ScatteringCoefficient mie = kMieAngstromBeta / kMieScaleHeight *
          pow(lambda, -kMieAngstromAlpha);
      solar_irradiance.push_back(solar);
      rayleigh_scattering.push_back(rayleigh);
      mie_scattering.push_back(mie * kMieSingleScatteringAlbedo);
      mie_extinction.push_back(mie);
      absorption_extinction.push_back(kMaxOzoneNumberDensity *
          kOzoneCrossSection[(l - kLambdaMin) / 10] * m2);
    }

    atmosphere_parameters_.solar_irradiance = IrradianceSpectrum(
        kLambdaMin * nm, kLambdaMax * nm, solar_irradiance);
    atmosphere_parameters_.sun_angular_radius = 0.2678 * deg;
    atmosphere_parameters_.bottom_radius = 6360.0 * km;
    atmosphere_parameters_.top_radius = 6420.0 * km;
    atmosphere_parameters_.rayleigh_density.layers[1] = DensityProfileLayer(
        0.0 * m, 1.0, -1.0 / kRayleighScaleHeight, 0.0 / m, 0.0);
    atmosphere_parameters_.rayleigh_scattering = ScatteringSpectrum(
        kLambdaMin * nm, kLambdaMax * nm, rayleigh_scattering);
    atmosphere_parameters_.mie_density.layers[1] = DensityProfileLayer(
        0.0 * m, 1.0, -1.0 / kMieScaleHeight, 0.0 / m, 0.0);
    atmosphere_parameters_.mie_scattering = ScatteringSpectrum(
        kLambdaMin * nm, kLambdaMax * nm, mie_scattering);
    atmosphere_parameters_.mie_extinction = ScatteringSpectrum(
        kLambdaMin * nm, kLambdaMax * nm, mie_extinction);
    atmosphere_parameters_.mie_phase_function_g = kMiePhaseFunctionG;
    // Density profile increasing linearly from 0 to 1 between 10 and 25km, and
    // decreasing linearly from 1 to 0 between 25 and 40km. Approximate profile
    // from http://www.kln.ac.lk/science/Chemistry/Teaching_Resources/Documents/
    // Introduction%20to%20atmospheric%20chemistry.pdf (page 10).
    atmosphere_parameters_.absorption_density.layers[0] = DensityProfileLayer(
        25.0 * km, 0.0, 0.0 / km, 1.0 / (15.0 * km), -2.0 / 3.0);
    atmosphere_parameters_.absorption_density.layers[1] = DensityProfileLayer(
        0.0 * km, 0.0, 0.0 / km, -1.0 / (15.0 * km), 8.0 / 3.0);
    atmosphere_parameters_.absorption_extinction = ScatteringSpectrum(
        kLambdaMin * nm, kLambdaMax * nm, absorption_extinction);
    atmosphere_parameters_.ground_albedo = DimensionlessSpectrum(0.1);
    atmosphere_parameters_.mu_s_min = cos(102.0 * deg);

    earth_center_ =
        Position(0.0 * m, 0.0 * m, -atmosphere_parameters_.bottom_radius);

    sun_size_ = dimensional::vec2(
        tan(atmosphere_parameters_.sun_angular_radius),
        cos(atmosphere_parameters_.sun_angular_radius));

    ground_albedo_ = GetGrassAlbedo();
    sphere_albedo_ = GetSnowAlbedo();
    program_ = 0;
  }

/*
<p>where the ground and sphere albedos are provided by the following methods:
*/

  DimensionlessSpectrum GetGrassAlbedo() {
    // Grass spectral albedo from Uwe Feister and Rolf Grewe, "Spectral albedo
    // measurements in the UV and visible region over different types of
    // surfaces", Photochemistry and Photobiology, 62, 736-744, 1995.
    constexpr double kGrassAlbedo[45] = {
      0.018, 0.019, 0.019, 0.020, 0.022, 0.024, 0.027, 0.029, 0.030, 0.031,
      0.032, 0.032, 0.032, 0.033, 0.035, 0.040, 0.055, 0.073, 0.084, 0.089,
      0.089, 0.079, 0.069, 0.063, 0.061, 0.057, 0.052, 0.051, 0.048, 0.042,
      0.039, 0.035, 0.035, 0.043, 0.087, 0.156, 0.234, 0.334, 0.437, 0.513,
      0.553, 0.571, 0.579, 0.581, 0.587
    };
    std::vector<Number> grass_albedo_samples;
    for (int i = 0; i < 45; ++i) {
      grass_albedo_samples.push_back(kGrassAlbedo[i]);
    }
    return DimensionlessSpectrum(360.0 * nm, 800.0 * nm, grass_albedo_samples);
  }

  DimensionlessSpectrum GetSnowAlbedo() {
    // Snow 5cm spectral albedo from Uwe Feister and Rolf Grewe, "Spectral
    // albedo measurements in the UV and visible region over different types of
    // surfaces", Photochemistry and Photobiology, 62, 736-744, 1995.
    constexpr double kSnowAlbedo[7] = {
      0.796, 0.802, 0.807, 0.810, 0.818, 0.825, 0.826
    };
    std::vector<Number> snow_albedo_samples;
    for (int i = 0; i < 7; ++i) {
      snow_albedo_samples.push_back(kSnowAlbedo[i]);
    }
    return DimensionlessSpectrum(360.0 * nm, 420.0 * nm, snow_albedo_samples);
  }

/*
<p>The GPU model is initialized differently depending on the test case, so we
provide a separate method to initialize it:
*/

  void InitGpuModel(bool combine_textures, bool precomputed_luminance) {
    if (!glutGet(GLUT_INIT_STATE)) {
      int argc = 0;
      char** argv = nullptr;
      glutInitContextVersion(3, 3);
      glutInitContextProfile(GLUT_CORE_PROFILE);
      glutInit(&argc, argv);
      glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
      glutInitWindowSize(kWidth, kHeight);
      glutCreateWindow("ModelTest");
      glutHideWindow();
      if (!gladLoadGL()) {
        throw std::runtime_error("GLAD initialization failed");
      }
      if (!GLAD_GL_VERSION_3_3) {
        throw std::runtime_error("OpenGL 3.3 or higher is required");
      }
    }

    std::vector<double> wavelengths;
    const auto& spectrum = atmosphere_parameters_.solar_irradiance;
    for (unsigned int i = 0; i < spectrum.size(); ++i) {
      wavelengths.push_back(spectrum.GetSample(i).to(nm));
    }
    auto profile = [](DensityProfileLayer layer) {
      return atmosphere::DensityProfileLayer(layer.width.to(m),
          layer.exp_term(), layer.exp_scale.to(1.0 / m),
          layer.linear_term.to(1.0 / m), layer.constant_term());
    };
    model_.reset(new atmosphere::Model(
        wavelengths,
        atmosphere_parameters_.solar_irradiance.to(
            watt_per_square_meter_per_nm),
        atmosphere_parameters_.sun_angular_radius.to(rad),
        atmosphere_parameters_.bottom_radius.to(m),
        atmosphere_parameters_.top_radius.to(m),
        {profile(atmosphere_parameters_.rayleigh_density.layers[1])},
        atmosphere_parameters_.rayleigh_scattering.to(1.0 / m),
        {profile(atmosphere_parameters_.mie_density.layers[1])},
        atmosphere_parameters_.mie_scattering.to(1.0 / m),
        atmosphere_parameters_.mie_extinction.to(1.0 / m),
        atmosphere_parameters_.mie_phase_function_g(),
        {profile(atmosphere_parameters_.absorption_density.layers[0]),
         profile(atmosphere_parameters_.absorption_density.layers[1])},
        atmosphere_parameters_.absorption_extinction.to(1.0 / m),
        atmosphere_parameters_.ground_albedo.to(Number::Unit()),
        acos(atmosphere_parameters_.mu_s_min()),
        kLengthUnit.to(m),
        precomputed_luminance ? 15 : 3 /* num_computed_wavelengths */,
        combine_textures,
        true /* half_precision */));
    model_->Init();
    glutSwapBuffers();
  }

/*
<p>Likewise, the CPU model might not be needed by all test cases, so we provide
a separate method to initialize it:
*/

  void InitCpuModel() {
    reference_model_.reset(
        new reference::Model(atmosphere_parameters_, "output/"));
    reference_model_->Init();
  }

/*
<p>Finally, before rendering an image with the GPU or CPU model, we must
initialize the camera (position, transform matrix, exposure) and the sun
direction and choose the rendering output (radiance or luminance). For this we
provide the following method:
*/

  void SetViewParameters(Angle sun_theta, Angle sun_phi, bool use_luminance) {
    // Transform matrix from camera frame to world space (i.e. the inverse of a
    // GL_MODELVIEW matrix).
    const float kCameraPos[3] = { 2000.0, -8000.0, 500.0 };
    constexpr float kPitch = PI / 30.0;
    const float model_from_view[16] = {
      1.0, 0.0, 0.0, kCameraPos[0],
      0.0, -sinf(kPitch), -cosf(kPitch), kCameraPos[1],
      0.0, cosf(kPitch), -sinf(kPitch), kCameraPos[2],
      0.0, 0.0, 0.0, 1.0
    };

    // Transform matrix from clip space to camera space (i.e. the inverse of a
    // GL_PROJECTION matrix).
    constexpr float kFovY = 50.0 / 180.0 * PI;
    const float kTanFovY = std::tan(kFovY / 2.0);
    const float view_from_clip[16] = {
      kTanFovY * static_cast<float>(kWidth) / kHeight, 0.0, 0.0, 0.0,
      0.0, kTanFovY, 0.0, 0.0,
      0.0, 0.0, 0.0, -1.0,
      0.0, 0.0, 1.0, 1.0
    };

    // Transform matrix from clip space to world space.
    for (int row = 0; row < 3; ++row) {
      for (int col = 0; col < 3; ++col) {
        int col2 = col < 2 ? col : 3;
        model_from_clip_[col + 3 * row] =
            model_from_view[0 + 4 * row] * view_from_clip[col2 + 0] +
            model_from_view[1 + 4 * row] * view_from_clip[col2 + 4] +
            model_from_view[2 + 4 * row] * view_from_clip[col2 + 8];
      }
    }

    camera_ = Position(kCameraPos[0] * m, kCameraPos[1] * m, kCameraPos[2] * m);
    exposure_ = use_luminance ? 1e-4 : 10.0;
    use_luminance_ = use_luminance;
    sun_direction_ = Direction(
        cos(sun_phi) * sin(sun_theta),
        sin(sun_phi) * sin(sun_theta),
        cos(sun_theta));
  }

/*
<p>The last "setup" method is the <code>TearDown</code> method, which is called
after each test case. It releases all the resources that might have been
allocated during the test:
*/

  void TearDown() override {
    model_ = nullptr;
    reference_model_ = nullptr;
    if (program_) {
     glDeleteProgram(program_);
    }
  }

/*
<h4 id="rendering">Rendering methods</h4>

<p>Once the test fixture is initialized with the above methods, we can render
one or more images, either with the GPU or with the CPU model. For the GPU model
we must first initialize the GPU shader, which can be done with the following
method:
*/

  void InitShader() {
    GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    const char* const vertex_shader_source = kVertexShader;
    glShaderSource(vertex_shader, 1, &vertex_shader_source, NULL);
    glCompileShader(vertex_shader);

    const std::string fragment_shader_str =
        "#version 330\n" +
        std::string(use_luminance_ ? "#define USE_LUMINANCE\n" : "") +
        std::string(kFragmentShader) + definitions_glsl +
        "const vec3 kSphereCenter = vec3(0.0, 0.0, " +
            std::to_string(kSphereRadius.to(kLengthUnit)) + ");\n" +
        "const float kSphereRadius = " +
            std::to_string(kSphereRadius.to(kLengthUnit)) + ";\n" +
        model_test_glsl;
    const char* fragment_shader_source = fragment_shader_str.c_str();
    GLuint fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment_shader, 1, &fragment_shader_source, NULL);
    glCompileShader(fragment_shader);

    if (program_) {
      glDeleteProgram(program_);
    }
    program_ = glCreateProgram();
    glAttachShader(program_, vertex_shader);
    glAttachShader(program_, fragment_shader);
    glAttachShader(program_, model_->shader());
    glLinkProgram(program_);
    glDetachShader(program_, vertex_shader);
    glDetachShader(program_, fragment_shader);
    glDetachShader(program_, model_->shader());
    glDeleteShader(vertex_shader);
    glDeleteShader(fragment_shader);

    glUseProgram(program_);
    model_->SetProgramUniforms(program_, 0, 1, 2, 3);
    glUniformMatrix3fv(glGetUniformLocation(program_, "model_from_clip"),
        1, true, model_from_clip_.data());
    glUniform3f(glGetUniformLocation(program_, "camera_"),
        camera_.x.to(kLengthUnit),
        camera_.y.to(kLengthUnit),
        camera_.z.to(kLengthUnit));
    glUniform1f(glGetUniformLocation(program_, "exposure_"),
        exposure_());
    glUniform3f(glGetUniformLocation(program_, "earth_center_"),
        earth_center_.x.to(kLengthUnit),
        earth_center_.y.to(kLengthUnit),
        earth_center_.z.to(kLengthUnit));
    glUniform3f(glGetUniformLocation(program_, "sun_direction_"),
        sun_direction_.x(),
        sun_direction_.y(),
        sun_direction_.z());
    glUniform2f(glGetUniformLocation(program_, "sun_size_"),
        sun_size_.x(), sun_size_.y());
    glUniform3f(glGetUniformLocation(program_, "ground_albedo_"),
        ground_albedo_(kLambdaR)(),
        ground_albedo_(kLambdaG)(),
        ground_albedo_(kLambdaB)());
    glUniform3f(glGetUniformLocation(program_, "sphere_albedo_"),
        sphere_albedo_(kLambdaR)(),
        sphere_albedo_(kLambdaG)(),
        sphere_albedo_(kLambdaB)());
  }

/*
<p>With the help of this method, we can now implement a method to render an
image with the GPU model. For this we just need to render a full screen quad
with the GPU program, and then read back the framebuffer pixels.
*/

  Image RenderGpuImage() {
    InitShader();

    glViewport(0, 0, kWidth, kHeight);
    {
      GLuint full_screen_quad_vao;
      glGenVertexArrays(1, &full_screen_quad_vao);
      glBindVertexArray(full_screen_quad_vao);
      GLuint full_screen_quad_vbo;
      glGenBuffers(1, &full_screen_quad_vbo);
      glBindBuffer(GL_ARRAY_BUFFER, full_screen_quad_vbo);
      const GLfloat vertices[] = {
        -1.0, -1.0, 0.0, 1.0,
        +1.0, -1.0, 0.0, 1.0,
        -1.0, +1.0, 0.0, 1.0,
        +1.0, +1.0, 0.0, 1.0,
      };
      glBufferData(GL_ARRAY_BUFFER, sizeof vertices, vertices, GL_STATIC_DRAW);
      constexpr GLuint kAttribIndex = 0;
      constexpr int kCoordsPerVertex = 4;
      glVertexAttribPointer(kAttribIndex, kCoordsPerVertex, GL_FLOAT,
                            false, 0, 0);
      glEnableVertexAttribArray(kAttribIndex);
      glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
      glBindBuffer(GL_ARRAY_BUFFER, 0);
      glDeleteBuffers(1, &full_screen_quad_vbo);
      glBindVertexArray(0);
      glDeleteVertexArrays(1, &full_screen_quad_vao);
    }
    glutSwapBuffers();

    std::unique_ptr<unsigned char[]> gl_pixels(
        new unsigned char[4 * kWidth * kHeight]);
    glReadPixels(
        0, 0, kWidth, kHeight, GL_RGBA, GL_UNSIGNED_BYTE, gl_pixels.get());

    Image pixels(new unsigned int[kWidth * kHeight]);
    for (unsigned int j = 0; j < kHeight; ++j) {
      for (unsigned int i = 0; i < kWidth; ++i) {
        int gl_offset = 4 * (i + (kHeight - 1 - j) * kWidth);
        int offset = i + j * kWidth;
        pixels[offset] = (255 << 24) |
            (gl_pixels[gl_offset] << 16) |
            (gl_pixels[gl_offset + 1] << 8) |
            gl_pixels[gl_offset + 2];
      }
    }
    return pixels;
  }

/*
<p>In order to render an image with the CPU model, we must first provide its
CPU implementation. For this, as for the unit tests of the GPU model, we simply
view the GLSL shader <a href="model_test.glsl.html">model_test.glsl</a> as C++
code (see the <a href="../index.html">Introduction</a>), that we include here
directly (after the definitions of the functions and macros it requires - the
"uniforms" are provided by the fields of the test fixture class, defined at the
end of this file):
*/

  RadianceSpectrum GetSolarRadiance() {
    return reference_model_->GetSolarRadiance();
  }

  RadianceSpectrum GetSkyRadiance(Position camera, Direction view_ray,
      Length shadow_length, Direction sun_direction,
      DimensionlessSpectrum& transmittance) {
    return reference_model_->GetSkyRadiance(
        camera, view_ray, shadow_length, sun_direction, &transmittance);
  }

  RadianceSpectrum GetSkyRadianceToPoint(Position camera, Position point,
      Length shadow_length, Direction sun_direction,
      DimensionlessSpectrum& transmittance) {
    return reference_model_->GetSkyRadianceToPoint(
        camera, point, shadow_length, sun_direction, &transmittance);
  }

  IrradianceSpectrum GetSunAndSkyIrradiance(Position point, Direction normal,
      Direction sun_direction, IrradianceSpectrum& sky_irradiance) {
    return reference_model_->GetSunAndSkyIrradiance(
        point, normal, sun_direction, &sky_irradiance);
  }

#define OUT(x) x&
#include "atmosphere/reference/model_test.glsl"

/*
<p>With this CPU implementation, we can render an image with a simple loop over
all the pixels, calling <code>GetViewRayRadiance</code> for each pixel, and
using the same tone mapping function as in the GPU version to convert the result
to a final color. The main difference with the GPU model is the conversion from
a radiance spectrum to an sRGB value, which must be done explicitely if a
luminance output is desired (otherwise, for radiance outputs, we simply need to
sample the radiance spectrum at the 3 predefined wavelengths):
*/

  Image RenderCpuImage() {
    constexpr auto kMaxLuminousEfficacy = MAX_LUMINOUS_EFFICACY * lm / watt;
    std::vector<Wavelength> wavelengths;
    std::vector<Number> x_values;
    std::vector<Number> y_values;
    std::vector<Number> z_values;
    for (unsigned int i = 0; i < 95 * 4; i += 4) {
      wavelengths.push_back(CIE_2_DEG_COLOR_MATCHING_FUNCTIONS[i] * nm);
      x_values.push_back(CIE_2_DEG_COLOR_MATCHING_FUNCTIONS[i + 1]);
      y_values.push_back(CIE_2_DEG_COLOR_MATCHING_FUNCTIONS[i + 2]);
      z_values.push_back(CIE_2_DEG_COLOR_MATCHING_FUNCTIONS[i + 3]);
    }
    const auto cie_x_bar = DimensionlessSpectrum(wavelengths, x_values);
    const auto cie_y_bar = DimensionlessSpectrum(wavelengths, y_values);
    const auto cie_z_bar = DimensionlessSpectrum(wavelengths, z_values);

    Image pixels(new unsigned int[kWidth * kHeight]);
    ProgressBar progress_bar(kWidth * kHeight);
    RunJobs([&](unsigned int j) {
      double y = 1.0 - 2.0 * (j + 0.5) / kHeight;
      double dy = -2.0 / kHeight;
      for (unsigned int i = 0; i < kWidth; ++i) {
        double x = 2.0 * (i + 0.5) / kWidth - 1.0;
        double dx = 2.0 / kWidth;

        Direction view_ray(
            model_from_clip_[0] * x + model_from_clip_[1] * y +
                model_from_clip_[2],
            model_from_clip_[3] * x + model_from_clip_[4] * y +
                model_from_clip_[5],
            model_from_clip_[6] * x + model_from_clip_[7] * y +
                model_from_clip_[8]);

        Direction view_ray_diff(
            model_from_clip_[0] * dx + model_from_clip_[1] * dy,
            model_from_clip_[3] * dx + model_from_clip_[4] * dy,
            model_from_clip_[6] * dx + model_from_clip_[7] * dy);

        RadianceSpectrum radiance = GetViewRayRadiance(view_ray, view_ray_diff);

        double r, g, b;
        if (use_luminance_) {
          Luminance x = kMaxLuminousEfficacy * Integral(radiance * cie_x_bar);
          Luminance y = kMaxLuminousEfficacy * Integral(radiance * cie_y_bar);
          Luminance z = kMaxLuminousEfficacy * Integral(radiance * cie_z_bar);
          r = (XYZ_TO_SRGB[0] * x + XYZ_TO_SRGB[1] * y + XYZ_TO_SRGB[2] * z).to(
              cd_per_square_meter);
          g = (XYZ_TO_SRGB[3] * x + XYZ_TO_SRGB[4] * y + XYZ_TO_SRGB[5] * z).to(
              cd_per_square_meter);
          b = (XYZ_TO_SRGB[6] * x + XYZ_TO_SRGB[7] * y + XYZ_TO_SRGB[8] * z).to(
              cd_per_square_meter);
        } else {
          r = radiance(kLambdaR).to(watt_per_square_meter_per_sr_per_nm);
          g = radiance(kLambdaG).to(watt_per_square_meter_per_sr_per_nm);
          b = radiance(kLambdaB).to(watt_per_square_meter_per_sr_per_nm);
        }

        r = std::pow(1.0 - std::exp(-r * exposure_()), 1.0 / 2.2);
        g = std::pow(1.0 - std::exp(-g * exposure_()), 1.0 / 2.2);
        b = std::pow(1.0 - std::exp(-b * exposure_()), 1.0 / 2.2);
        unsigned int red = static_cast<unsigned int>(r * 255.0);
        unsigned int green = static_cast<unsigned int>(g * 255.0);
        unsigned int blue = static_cast<unsigned int>(b * 255.0);
        pixels[i + j * kWidth] =
            (255 << 24) | (red << 16) | (green << 8) | blue;
        progress_bar.Increment(1);
      }
    }, kHeight);
    return pixels;
  }

/*
<h4 id="comparison">Comparison methods</h4>

<p>After some images have been rendered, we want to compare them in order to
check whether two images of the same scene, rendered with different methods, are
close enough or not. This is the goal of the following method, which uses the
<a href="https://fr.wikipedia.org/wiki/Peak_Signal_to_Noise_Ratio">Peak Signal
to Noise Ratio</a> as the image difference measure.
*/

  double ComputePSNR(const unsigned int* image1, const unsigned int* image2) {
    double square_error_sum = 0.0;
    for (unsigned int j = 0; j < kHeight; ++j) {
      for (unsigned int i = 0; i < kWidth; ++i) {
        int argb1 = image1[i + j * kWidth];
        int argb2 = image2[i + j * kWidth];
        dimensional::Vector3<Number> rgb1(
            (argb1 >> 16) & 0xFF, (argb1 >> 8) & 0xFF, argb1 & 0xFF);
        dimensional::Vector3<Number> rgb2(
            (argb2 >> 16) & 0xFF, (argb2 >> 8) & 0xFF, argb2 & 0xFF);
        square_error_sum += dot(rgb1 - rgb2, rgb1 - rgb2)();
      }
    }
    double mean_square_error = sqrt(square_error_sum / (kWidth * kHeight));
    return 10.0 * log(255 * 255 / mean_square_error) / log(10.0);
  }

/*
<p>Also, in order to visually compare the images, it is useful to have an HTML
test report, showing for each test case its two images and their PSNR score
difference. For this, the following method compares two images, writes them to
disk, creates or appends a test report entry in a test report file, and finally
returns the computed PSNR.
*/

  double Compare(Image image1, Image image2, const std::string& caption,
      bool append) {
    double psnr = ComputePSNR(image1.get(), image2.get());
    WritePngArgb(name_ + "1.png", image1.get());
    WritePngArgb(name_ + "2.png", image2.get());
    std::ofstream file(std::string(kOutputDir) + "test_report.html",
        append ? std::ios_base::app : std::ios_base::trunc);
    file << "<h2>" << name_ << " (PSNR = " << psnr << "dB)</h2>" << std::endl
         << "<p>" << caption << std::endl
         << "<p><img src=\"" << name_ << "1.png\">" << std::endl
         << "<img src=\"" << name_ << "2.png\">" << std::endl;
    file.close();
    return psnr;
  }

/*
<h3 id="cases">Test cases</h3>

<p>Thanks to the above helper methods we can finally implement some test cases
to make sure that the GPU model, despite its approximations, is almost as
accurate as the full spectral, but much slower CPU model.

<p>The first test case compares the radiance computations, done on GPU vs CPU.
It uses a separate RGB texture for the single Mie scattering, which means that
the computations done on GPU and CPU, for the 3 wavelengths of the GPU model,
are exactly the same (except for the floating point precision: single precision
on GPU - with half precision in the precomputed textures, vs double precision on
CPU). As a consequence, we expect the two images to be nearly identical:
*/

  void TestRadianceSeparateTextures() {
    const std::string kCaption = "Left: GPU model, combine_textures = false. "
        "Right: CPU model. Both images show the spectral radiance at 3 "
        "predefined wavelengths (i.e. no conversion to sRGB via CIE XYZ).";
    InitGpuModel(false /* combine_textures */,
        false /* precomputed_luminance */);
    InitCpuModel();
    SetViewParameters(65.0 * deg, 90.0 * deg, false /* use_luminance */);
    ExpectLess(
        47.0, Compare(RenderGpuImage(), RenderCpuImage(), kCaption, false));
  }

/*
<p>The following test case is almost the same as the previous one, except that
we use the the combine_textures option in the GPU model. This leads to some
approximations on the GPU side, not present in the CPU model. As a consequence,
we expect a slightly larger difference between the two images, compared to the
previous test case:
*/

  void TestRadianceCombineTextures() {
    const std::string kCaption = "Left: GPU model, combine_textures = true. "
        "Right: CPU model. Both images show the spectral radiance at 3 "
        "predefined wavelengths (i.e. no conversion to sRGB via CIE XYZ).";
    InitGpuModel(true /* combine_textures */,
        false /* precomputed_luminance */);
    InitCpuModel();
    SetViewParameters(65.0 * deg, 90.0 * deg, false /* use_luminance */);
    ExpectLess(
        46.0, Compare(RenderGpuImage(), RenderCpuImage(), kCaption, true));
  }

/*
<p>The following test case is the same as the previous one, for a sunset scene.
In this case the single Mie scattering contribution is larger than in the
previous test, so we expect a larger difference between the GPU and CPU results
(due to the GPU approximations for the single Mie scattering term):
*/

  void TestRadianceCombineTexturesSunSet() {
    const std::string kCaption = "Left: GPU model, combine_textures = true. "
        "Right: CPU model. Both images show the spectral radiance at 3 "
        "predefined wavelengths (i.e. no conversion to sRGB via CIE XYZ).";
    InitGpuModel(true /* combine_textures */,
        false /* precomputed_luminance */);
    InitCpuModel();
    SetViewParameters(88.0 * deg, 90.0 * deg, false /* use_luminance */);
    ExpectLess(
        40.0, Compare(RenderGpuImage(), RenderCpuImage(), kCaption, true));
  }

/*
<p>The following test case compares the sRGB luminance computations, done on GPU
vs CPU. The GPU computations use approximations to obtain the sRGB values from
the value of the spectral radiance at only 3 wavelengths, while the CPU
computations use a "full spectral" rendering method (using the value of the
spectral radiance at 47 wavelengths between 360 and 830 nm). To evaluate the
effect of these approximations alone, we don't use the combine_textures option
on GPU (which introduces additional approximations). Similarly, we use constant
albedos instead of wavelength dependent albedos to avoid additional
approximations:
*/

  void TestLuminanceSeparateTexturesConstantAlbedo() {
    const std::string kCaption = "Left: GPU model, combine_textures = false. "
        "Right: CPU model. Both images show the sRGB luminance (radiance "
        "converted to CIE XYZ and then to sRGB - with some approximations on "
        "GPU).";
    sphere_albedo_ = DimensionlessSpectrum(0.8);
    ground_albedo_ = DimensionlessSpectrum(0.1);
    InitGpuModel(false /* combine_textures */,
        false /* precomputed_luminance */);
    InitCpuModel();
    SetViewParameters(65.0 * deg, 90.0 * deg, true /* use_luminance */);
    ExpectLess(
        40.0, Compare(RenderGpuImage(), RenderCpuImage(), kCaption, true));
  }

/*
<p>The following test case is almost the same as the previous one, except that
we use the the combine_textures option in the GPU model. This leads to some
additional approximations on the GPU side, not present in the CPU model. As a
consequence, we expect a slightly larger difference between the two images,
compared to the previous test case:
*/

  void TestLuminanceCombineTexturesConstantAlbedo() {
    const std::string kCaption = "Left: GPU model, combine_textures = true. "
        "Right: CPU model. Both images show the sRGB luminance (radiance "
        "converted to CIE XYZ and then to sRGB - with some approximations on "
        "GPU).";
    sphere_albedo_ = DimensionlessSpectrum(0.8);
    ground_albedo_ = DimensionlessSpectrum(0.1);
    InitGpuModel(true /* combine_textures */,
        false /* precomputed_luminance */);
    InitCpuModel();
    SetViewParameters(65.0 * deg, 90.0 * deg, true /* use_luminance */);
    ExpectLess(
        40.0, Compare(RenderGpuImage(), RenderCpuImage(), kCaption, true));
  }

/*
<p>The following test case is the same as the previous one, for a sunset scene.
In this case the single Mie scattering contribution is larger than in the
previous test, so we expect a larger difference between the GPU and CPU results
(due to the GPU approximations for the single Mie scattering term):
*/

  void TestLuminanceCombineTexturesConstantAlbedoSunSet() {
    const std::string kCaption = "Left: GPU model, combine_textures = true. "
        "Right: CPU model. Both images show the sRGB luminance (radiance "
        "converted to CIE XYZ and then to sRGB - with some approximations on "
        "GPU).";
    sphere_albedo_ = DimensionlessSpectrum(0.8);
    ground_albedo_ = DimensionlessSpectrum(0.1);
    InitGpuModel(true /* combine_textures */,
        false /* precomputed_luminance */);
    InitCpuModel();
    SetViewParameters(88.0 * deg, 90.0 * deg, true /* use_luminance */);
    ExpectLess(
        35.0, Compare(RenderGpuImage(), RenderCpuImage(), kCaption, true));
  }

/*
<p>The following test case compares the sRGB luminance computations, done on GPU
vs CPU, with wavelength dependent albedo values. This leads, on the GPU side, to
new approximations compared to the CPU reference model (indeed, on GPU we
convert the sky and sun radiance to sRGB, and then multiply these sRGB values by
the albedo, sampled at 3 wavelengths - while the correct method, used on CPU, is
to perform all the computations, including the albedo multiplication, in a full
spectral way, and to convert to XYZ and then to sRGB only at the very end).
*/

  void TestLuminanceCombineTexturesSpectralAlbedo() {
    const std::string kCaption = "Left: GPU model, combine_textures = true. "
        "Right: CPU model. Both images show the sRGB luminance (radiance "
        "converted to CIE XYZ and then to sRGB - with some approximations on "
        "GPU).";
    InitGpuModel(true /* combine_textures */,
        false /* precomputed_luminance */);
    InitCpuModel();
    SetViewParameters(65.0 * deg, 90.0 * deg, true /* use_luminance */);
    ExpectLess(
        38.0, Compare(RenderGpuImage(), RenderCpuImage(), kCaption, true));
  }

/*
<p>The following test case compares the sRGB luminance computations, done on GPU
vs CPU, in a "worst case" situation: combined textures on GPU and a sunset scene
(leading to large differences in the single Mie component), and wavelength
dependent albedo values (see the previous test case):
*/

  void TestLuminanceCombineTexturesSpectralAlbedoSunSet() {
    const std::string kCaption = "Left: GPU model, combine_textures = true. "
        "Right: CPU model. Both images show the sRGB luminance (radiance "
        "converted to CIE XYZ and then to sRGB - with some approximations on "
        "GPU).";
    InitGpuModel(true /* combine_textures */,
        false /* precomputed_luminance */);
    InitCpuModel();
    SetViewParameters(88.0 * deg, 90.0 * deg, true /* use_luminance */);
    ExpectLess(
        35.0, Compare(RenderGpuImage(), RenderCpuImage(), kCaption, true));
  }

/*
<p>The following test case compares the sRGB luminance computations, done on GPU
vs CPU. The GPU computations use precomputed luminance from 15 wavelenths, while
the CPU computations use a "full spectral" rendering method (using the value of
the spectral radiance at 47 wavelengths between 360 and 830 nm). To evaluate the
effect of these approximations alone, we don't use the combine_textures option
on GPU (which introduces additional approximations). Similarly, we use constant
albedos instead of wavelength dependent albedos to avoid additional
approximations:
*/

  void TestPrecomputedLuminanceSeparateTexturesConstantAlbedo() {
    const std::string kCaption = "Left: GPU model, combine_textures = false. "
        "Right: CPU model. Both images show the sRGB luminance (radiance "
        "converted to CIE XYZ and then to sRGB - using 15 wavelengths on GPU, "
        "vs 47 on CPU).";
    sphere_albedo_ = DimensionlessSpectrum(0.8);
    ground_albedo_ = DimensionlessSpectrum(0.1);
    InitGpuModel(false /* combine_textures */,
        true /* precomputed_luminance */);
    InitCpuModel();
    SetViewParameters(65.0 * deg, 90.0 * deg, true /* use_luminance */);
    ExpectLess(
        43.0, Compare(RenderGpuImage(), RenderCpuImage(), kCaption, true));
  }

/*
<p>The following test case is almost the same as the previous one, except that
we use the the combine_textures option in the GPU model. This leads to some
additional approximations on the GPU side, not present in the CPU model. As a
consequence, we expect a slightly larger difference between the two images,
compared to the previous test case:
*/

  void TestPrecomputedLuminanceCombineTexturesConstantAlbedo() {
    const std::string kCaption = "Left: GPU model, combine_textures = true. "
        "Right: CPU model. Both images show the sRGB luminance (radiance "
        "converted to CIE XYZ and then to sRGB - using 15 wavelengths on GPU, "
        "vs 47 on CPU).";
    sphere_albedo_ = DimensionlessSpectrum(0.8);
    ground_albedo_ = DimensionlessSpectrum(0.1);
    InitGpuModel(true /* combine_textures */,
        true /* precomputed_luminance */);
    InitCpuModel();
    SetViewParameters(65.0 * deg, 90.0 * deg, true /* use_luminance */);
    ExpectLess(
        43.0, Compare(RenderGpuImage(), RenderCpuImage(), kCaption, true));
  }

/*
<p>The following test case is the same as the previous one, for a sunset scene.
In this case the single Mie scattering contribution is larger than in the
previous test, so we expect a larger difference between the GPU and CPU results
(due to the GPU approximations for the single Mie scattering term):
*/

  void TestPrecomputedLuminanceCombineTexturesConstantAlbedoSunSet() {
    const std::string kCaption = "Left: GPU model, combine_textures = true. "
        "Right: CPU model. Both images show the sRGB luminance (radiance "
        "converted to CIE XYZ and then to sRGB - using 15 wavelengths on GPU, "
        "vs 47 on CPU).";
    sphere_albedo_ = DimensionlessSpectrum(0.8);
    ground_albedo_ = DimensionlessSpectrum(0.1);
    InitGpuModel(true /* combine_textures */,
        true /* precomputed_luminance */);
    InitCpuModel();
    SetViewParameters(88.0 * deg, 90.0 * deg, true /* use_luminance */);
    ExpectLess(
        40.0, Compare(RenderGpuImage(), RenderCpuImage(), kCaption, true));
  }

/*
<p>The following test case compares the sRGB luminance computations, done on GPU
vs CPU, with wavelength dependent albedo values. This leads, on the GPU side, to
new approximations compared to the CPU reference model (indeed, on GPU we
multiply the sun and sky sRGB values by the albedo, sampled at 3 wavelengths -
while the correct method, used on CPU, is to perform all the computations,
including the albedo multiplication, in a full spectral way, and to convert to
XYZ and then to sRGB only at the very end).
*/

  void TestPrecomputedLuminanceCombineTexturesSpectralAlbedo() {
    const std::string kCaption = "Left: GPU model, combine_textures = true. "
        "Right: CPU model. Both images show the sRGB luminance (radiance "
        "converted to CIE XYZ and then to sRGB - using 15 wavelengths on GPU, "
        "vs 47 on CPU).";
    InitGpuModel(true /* combine_textures */,
        true /* precomputed_luminance */);
    InitCpuModel();
    SetViewParameters(65.0 * deg, 90.0 * deg, true /* use_luminance */);
    ExpectLess(
        39.0, Compare(RenderGpuImage(), RenderCpuImage(), kCaption, true));
  }

/*
<p>The last test case compares the sRGB luminance computations, done on GPU
vs CPU, in a "worst case" situation: combined textures on GPU and a sunset scene
(leading to large differences in the single Mie component), and wavelength
dependent albedo values (see the previous test case):
*/

  void TestPrecomputedLuminanceCombineTexturesSpectralAlbedoSunSet() {
    const std::string kCaption = "Left: GPU model, combine_textures = true. "
        "Right: CPU model. Both images show the sRGB luminance (radiance "
        "converted to CIE XYZ and then to sRGB - using 15 wavelengths on GPU, "
        "vs 47 on CPU).";
    InitGpuModel(true /* combine_textures */,
        true /* precomputed_luminance */);
    InitCpuModel();
    SetViewParameters(88.0 * deg, 90.0 * deg, true /* use_luminance */);
    ExpectLess(
        40.0, Compare(RenderGpuImage(), RenderCpuImage(), kCaption, true));
  }

/*
<p> The rest of the code simply declares the fields of our test fixture class,
and registers the test cases in the test framework:
*/

 private:
  std::string name_;
  AtmosphereParameters atmosphere_parameters_;
  DimensionlessSpectrum ground_albedo_;
  DimensionlessSpectrum sphere_albedo_;
  Position earth_center_;
  dimensional::vec2 sun_size_;

  std::unique_ptr<atmosphere::Model> model_;
  std::unique_ptr<reference::Model> reference_model_;
  GLuint program_;

  std::array<float, 9> model_from_clip_;
  Position camera_;
  Number exposure_;
  bool use_luminance_;
  Direction sun_direction_;
};

namespace {

ModelTest radiance1(
    "RadianceSeparateTextures",
    &ModelTest::TestRadianceSeparateTextures);
ModelTest radiance2(
    "RadianceCombineTextures",
    &ModelTest::TestRadianceCombineTextures);
ModelTest radiance3(
    "RadianceCombineTexturesSunSet",
    &ModelTest::TestRadianceCombineTexturesSunSet);
ModelTest luminance1(
    "LuminanceSeparateTexturesConstantAlbedo",
    &ModelTest::TestLuminanceSeparateTexturesConstantAlbedo);
ModelTest luminance2(
    "LuminanceCombineTexturesConstantAlbedo",
    &ModelTest::TestLuminanceCombineTexturesConstantAlbedo);
ModelTest luminance3(
    "LuminanceCombineTexturesConstantAlbedoSunSet",
    &ModelTest::TestLuminanceCombineTexturesConstantAlbedoSunSet);
ModelTest luminance4(
    "LuminanceCombineTexturesSpectralAlbedo",
    &ModelTest::TestLuminanceCombineTexturesSpectralAlbedo);
ModelTest luminance5(
    "LuminanceCombineTexturesSpectralAlbedoSunSet",
    &ModelTest::TestLuminanceCombineTexturesSpectralAlbedoSunSet);
ModelTest precomputed_luminance1(
    "PrecomputedLuminanceSeparateTexturesConstantAlbedo",
    &ModelTest::TestPrecomputedLuminanceSeparateTexturesConstantAlbedo);
ModelTest precomputed_luminance2(
    "PrecomputedLuminanceCombineTexturesConstantAlbedo",
    &ModelTest::TestPrecomputedLuminanceCombineTexturesConstantAlbedo);
ModelTest precomputed_luminance3(
    "PrecomputedLuminanceCombineTexturesConstantAlbedoSunSet",
    &ModelTest::TestPrecomputedLuminanceCombineTexturesConstantAlbedoSunSet);
ModelTest precomputed_luminance4(
    "PrecomputedLuminanceCombineTexturesSpectralAlbedo",
    &ModelTest::TestPrecomputedLuminanceCombineTexturesSpectralAlbedo);
ModelTest precomputed_luminance5(
    "PrecomputedLuminanceCombineTexturesSpectralAlbedoSunSet",
    &ModelTest::TestPrecomputedLuminanceCombineTexturesSpectralAlbedoSunSet);

}  // anonymous namespace

}  // namespace reference
}  // namespace atmosphere
