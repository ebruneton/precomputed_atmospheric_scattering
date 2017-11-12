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

/*<h2>atmosphere/reference/definitions.h</h2>

<p>This C++ file defines the physical types and constants which are used in the
main GLSL <a href="../functions.glsl.html">functions</a> of our atmosphere
model, in such a way that they can be compiled by the C++ compiler. The <a href=
"../definitions.glsl.html">GLSL equivalent</a> of this file provides the same
types and constants in GLSL, to allow the same functions to be compiled by the
GLSL compiler.

<p>The main purpose of this C++ compilation is to check the <a href="
https://en.wikipedia.org/wiki/Dimensional_analysis#Dimensional_homogeneity"
>dimensional homogeneity</a> of the GLSL expressions (see the
<a href="../index.html">Introduction</a>). For this we define the C++ physical
types by using generic templates parameterized by physical dimensions,
inspired from <a href=
"http://www.boost.org/doc/libs/1_61_0/doc/html/boost_units.html">Boost.Unit</a>,
and provided by the following included files:
*/

#ifndef ATMOSPHERE_REFERENCE_DEFINITIONS_H_
#define ATMOSPHERE_REFERENCE_DEFINITIONS_H_

#include "atmosphere/constants.h"
#include "math/angle.h"
#include "math/binary_function.h"
#include "math/scalar.h"
#include "math/scalar_function.h"
#include "math/ternary_function.h"
#include "math/vector.h"

namespace atmosphere {
namespace reference {

/*
<h3>Physical quantities</h3>

<p>The physical quantities we need for our atmosphere model are
<a href="https://en.wikipedia.org/wiki/Radiometry">radiometric</a> and
<a href="https://en.wikipedia.org/wiki/Photometry_(optics)">photometric</a>
quantities. We start with six base quantities: angle, length, wavelength, solid
angle, power and luminous power (wavelength is also a length, but we distinguish
the two for increased clarity).
*/

typedef dimensional::Angle Angle;
typedef dimensional::Scalar<1, 0, 0, 0, 0> Length;
typedef dimensional::Scalar<0, 1, 0, 0, 0> Wavelength;
typedef dimensional::Scalar<0, 0, 1, 0, 0> SolidAngle;
typedef dimensional::Scalar<0, 0, 0, 1, 0> Power;
typedef dimensional::Scalar<0, 0, 0, 0, 1> LuminousPower;

/*
<p>From this we derive the irradiance, radiance, spectral irradiance,
spectral radiance, luminance, etc, as well pure numbers, area, volume, etc.
*/

typedef dimensional::Scalar<0, 0, 0, 0, 0> Number;
typedef dimensional::Scalar<-1, 0, 0, 0, 0> InverseLength;
typedef dimensional::Scalar<2, 0, 0, 0, 0> Area;
typedef dimensional::Scalar<3, 0, 0, 0, 0> Volume;
typedef dimensional::Scalar<-2, 0, 0, 1, 0> Irradiance;
typedef dimensional::Scalar<-2, 0, -1, 1, 0> Radiance;
typedef dimensional::Scalar<0, -1, 0, 1, 0> SpectralPower;
typedef dimensional::Scalar<-2, -1, 0, 1, 0> SpectralIrradiance;
typedef dimensional::Scalar<-2, -1, -1, 1, 0> SpectralRadiance;
typedef dimensional::Scalar<-3, -1, -1, 1, 0> SpectralRadianceDensity;
typedef dimensional::Scalar<-1, 0, 0, 0, 0> ScatteringCoefficient;
typedef dimensional::Scalar<0, 0, -1, 0, 0> InverseSolidAngle;
typedef dimensional::Scalar<-3, 0, 0, 0, 0> NumberDensity;
typedef dimensional::Scalar<0, 0, -1, 0, 1> LuminousIntensity;
typedef dimensional::Scalar<-2, 0, -1, 0, 1> Luminance;
typedef dimensional::Scalar<-2, 0, 0, 0, 1> Illuminance;

/*
<p>We also need vectors of physical quantities, mostly to represent functions
depending on the wavelength. In this case the vector elements correspond to
values of a function at some predefined wavelengths. Here we use 47 predefined
wavelengths, uniformly distributed between 360 and 830 nanometers:
*/

template<int U1, int U2, int U3, int U4, int U5>
using WavelengthFunction = dimensional::ScalarFunction<
    0, 1, 0, 0, 0, U1, U2, U3, U4, U5, 47, 360, 830>;

// A function from Wavelength to Number.
typedef WavelengthFunction<0, 0, 0, 0, 0> DimensionlessSpectrum;
// A function from Wavelength to SpectralPower.
typedef WavelengthFunction<0, -1, 0, 1, 0> PowerSpectrum;
// A function from Wavelength to SpectralIrradiance.
typedef WavelengthFunction<-2, -1, 0, 1, 0> IrradianceSpectrum;
// A function from Wavelength to SpectralRadiance.
typedef WavelengthFunction<-2, -1, -1, 1, 0> RadianceSpectrum;
// A function from Wavelength to SpectralRadianceDensity.
typedef WavelengthFunction<-3, -1, -1, 1, 0> RadianceDensitySpectrum;
// A function from Wavelength to ScaterringCoefficient.
typedef WavelengthFunction<-1, 0, 0, 0, 0> ScatteringSpectrum;

// A position in 3D (3 length values).
typedef dimensional::Vector3<Length> Position;
// A unit direction vector in 3D (3 unitless values).
typedef dimensional::Vector3<Number> Direction;
// A vector of 3 luminance values.
typedef dimensional::Vector3<Luminance> Luminance3;
// A vector of 3 illuminance values.
typedef dimensional::Vector3<Illuminance> Illuminance3;

/*
<p>Finally, we also need precomputed textures containing physical quantities in
each texel (the texture sizes are defined in
<a href="../constants.h.html"><code>constants.h</code></a>):
*/

typedef dimensional::BinaryFunction<
    TRANSMITTANCE_TEXTURE_WIDTH,
    TRANSMITTANCE_TEXTURE_HEIGHT,
    DimensionlessSpectrum> TransmittanceTexture;

template<class T>
using AbstractScatteringTexture = dimensional::TernaryFunction<
    SCATTERING_TEXTURE_WIDTH,
    SCATTERING_TEXTURE_HEIGHT,
    SCATTERING_TEXTURE_DEPTH,
    T>;

typedef AbstractScatteringTexture<IrradianceSpectrum>
    ReducedScatteringTexture;

typedef AbstractScatteringTexture<RadianceSpectrum>
    ScatteringTexture;

typedef AbstractScatteringTexture<RadianceDensitySpectrum>
    ScatteringDensityTexture;

typedef dimensional::BinaryFunction<
    IRRADIANCE_TEXTURE_WIDTH,
    IRRADIANCE_TEXTURE_HEIGHT,
    IrradianceSpectrum> IrradianceTexture;

/*
<h3>Physical units</h3>

<p>We can then define the units for our base physical quantities:
radians (rad), meter (m), nanometer (nm), steradian (sr), watt (watt) and lumen
(lm):
*/

constexpr Angle rad = dimensional::rad;
constexpr Length m = Length::Unit();
constexpr Wavelength nm = Wavelength::Unit();
constexpr SolidAngle sr = SolidAngle::Unit();
constexpr Power watt = Power::Unit();
constexpr LuminousPower lm = LuminousPower::Unit();

/*
<p>From which we can derive the units for some derived physical quantities,
as well as some derived units (degress deg, kilometer km, kilocandela kcd):
*/

constexpr double PI = dimensional::PI;
constexpr Angle pi = dimensional::pi;
constexpr Angle deg = dimensional::deg;
constexpr Length km = 1000.0 * m;
constexpr Area m2 = m * m;
constexpr Volume m3 = m * m * m;
constexpr Irradiance watt_per_square_meter = watt / m2;
constexpr Radiance watt_per_square_meter_per_sr = watt / (m2 * sr);
constexpr SpectralIrradiance watt_per_square_meter_per_nm = watt / (m2 * nm);
constexpr SpectralRadiance watt_per_square_meter_per_sr_per_nm =
    watt / (m2 * sr * nm);
constexpr SpectralRadianceDensity watt_per_cubic_meter_per_sr_per_nm =
    watt / (m3 * sr * nm);
constexpr LuminousIntensity cd = lm / sr;
constexpr LuminousIntensity kcd = 1000.0 * cd;
constexpr Luminance cd_per_square_meter = cd / m2;
constexpr Luminance kcd_per_square_meter = kcd / m2;

/*
<h3>Atmosphere parameters</h3>

<p>Using the above types, we can now define the parameters of our atmosphere
model. We start with the definition of density profiles, which are needed for
parameters that depend on the altitude:
*/

// An atmosphere layer of width 'width', and whose density is defined as
//   'exp_term' * exp('exp_scale' * h) + 'linear_term' * h + 'constant_term',
// clamped to [0,1], and where h is the altitude.
struct DensityProfileLayer {
  DensityProfileLayer() :
      DensityProfileLayer(0.0 * m, 0.0, 0.0 / m, 0.0 / m, 0.0) {}
  DensityProfileLayer(Length width, Number exp_term, InverseLength exp_scale,
                      InverseLength linear_term, Number constant_term)
      : width(width), exp_term(exp_term), exp_scale(exp_scale),
        linear_term(linear_term), constant_term(constant_term) {
  }
  Length width;
  Number exp_term;
  InverseLength exp_scale;
  InverseLength linear_term;
  Number constant_term;
};

// An atmosphere density profile made of several layers on top of each other
// (from bottom to top). The width of the last layer is ignored, i.e. it always
// extend to the top atmosphere boundary. The profile values vary between 0
// (null density) to 1 (maximum density).
struct DensityProfile {
  DensityProfileLayer layers[2];
};

struct AtmosphereParameters {
  // The solar irradiance at the top of the atmosphere.
  IrradianceSpectrum solar_irradiance;
  // The sun's angular radius. Warning: the implementation uses approximations
  // that are valid only if this angle is smaller than 0.1 radians.
  Angle sun_angular_radius;
  // The distance between the planet center and the bottom of the atmosphere.
  Length bottom_radius;
  // The distance between the planet center and the top of the atmosphere.
  Length top_radius;
  // The density profile of air molecules, i.e. a function from altitude to
  // dimensionless values between 0 (null density) and 1 (maximum density).
  DensityProfile rayleigh_density;
  // The scattering coefficient of air molecules at the altitude where their
  // density is maximum (usually the bottom of the atmosphere), as a function of
  // wavelength. The scattering coefficient at altitude h is equal to
  // 'rayleigh_scattering' times 'rayleigh_density' at this altitude.
  ScatteringSpectrum rayleigh_scattering;
  // The density profile of aerosols, i.e. a function from altitude to
  // dimensionless values between 0 (null density) and 1 (maximum density).
  DensityProfile mie_density;
  // The scattering coefficient of aerosols at the altitude where their density
  // is maximum (usually the bottom of the atmosphere), as a function of
  // wavelength. The scattering coefficient at altitude h is equal to
  // 'mie_scattering' times 'mie_density' at this altitude.
  ScatteringSpectrum mie_scattering;
  // The extinction coefficient of aerosols at the altitude where their density
  // is maximum (usually the bottom of the atmosphere), as a function of
  // wavelength. The extinction coefficient at altitude h is equal to
  // 'mie_extinction' times 'mie_density' at this altitude.
  ScatteringSpectrum mie_extinction;
  // The asymetry parameter for the Cornette-Shanks phase function for the
  // aerosols.
  Number mie_phase_function_g;
  // The density profile of air molecules that absorb light (e.g. ozone), i.e.
  // a function from altitude to dimensionless values between 0 (null density)
  // and 1 (maximum density).
  DensityProfile absorption_density;
  // The extinction coefficient of molecules that absorb light (e.g. ozone) at
  // the altitude where their density is maximum, as a function of wavelength.
  // The extinction coefficient at altitude h is equal to
  // 'absorption_extinction' times 'absorption_density' at this altitude.
  ScatteringSpectrum absorption_extinction;
  // The average albedo of the ground.
  DimensionlessSpectrum ground_albedo;
  // The cosine of the maximum Sun zenith angle for which atmospheric scattering
  // must be precomputed (for maximum precision, use the smallest Sun zenith
  // angle yielding negligible sky light radiance values. For instance, for the
  // Earth case, 102 degrees is a good choice - yielding mu_s_min = -0.2).
  Number mu_s_min;
};

}  // namespace reference
}  // namespace atmosphere

#endif  // ATMOSPHERE_REFERENCE_DEFINITIONS_H_
