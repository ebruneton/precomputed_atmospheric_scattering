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

/*<h2>atmosphere/reference/functions_test.cc</h2>

<p>This file provides unit tests for the <a href="../functions.glsl.html">GLSL
functions</a> that implement our atmosphere model. We start with the definition
of some (arbitrary) values for the atmosphere parameters:
*/

#include "atmosphere/reference/functions.h"

#include <limits>
#include <string>

#include "atmosphere/reference/definitions.h"
#include "atmosphere/constants.h"
#include "test/test_case.h"

namespace atmosphere {
namespace reference {

constexpr double kEpsilon = 1e-3;

constexpr SpectralIrradiance kSolarIrradiance =
    123.0 * watt_per_square_meter_per_nm;
constexpr Length kBottomRadius = 1000.0 * km;
constexpr Length kTopRadius = 1500.0 * km;
constexpr Length kScaleHeight = 60.0 * km;
constexpr Length kRayleighScaleHeight = 60.0 * km;
constexpr Length kMieScaleHeight = 30.0 * km;
constexpr ScatteringCoefficient kRayleighScattering = 0.001 / km;
constexpr ScatteringCoefficient kMieScattering = 0.0015 / km;
constexpr ScatteringCoefficient kMieExtinction = 0.002 / km;
constexpr Number kGroundAlbedo = 0.1;

/*
<p>This helper function computes the cosine of the angle between the zenith and
the horizon, at some radius r. We will use it to test rays just above or below
the horizon.
*/

Number CosineOfHorizonZenithAngle(Length r) {
  assert(r >= kBottomRadius);
  return -sqrt(1.0 - (kBottomRadius / r) * (kBottomRadius / r));
}

/*
<p>Some unit tests need a precomputed texture as input, but we don't want to
precompute a whole texture for that, for efficiency reasons. Our solution is to
provide lazily computed textures instead, i.e. textures whose texels are
computed the first time we try to read them. The first type of texture we need
is a lazy transmittance texture (negative values mean "not yet computed"):
*/

class LazyTransmittanceTexture :
    public dimensional::BinaryFunction<TRANSMITTANCE_TEXTURE_WIDTH,
                                       TRANSMITTANCE_TEXTURE_HEIGHT,
                                       DimensionlessSpectrum> {
 public:
  explicit LazyTransmittanceTexture(
      const AtmosphereParameters& atmosphere_parameters)
    : BinaryFunction(DimensionlessSpectrum(-1.0)),
      atmosphere_parameters_(atmosphere_parameters) {
  }

  virtual const DimensionlessSpectrum& Get(int i, int j) const {
    int index = i + j * TRANSMITTANCE_TEXTURE_HEIGHT;
    if (value_[index][0]() < 0.0) {
      value_[index] = ComputeTransmittanceToTopAtmosphereBoundaryTexture(
          atmosphere_parameters_, vec2(i + 0.5, j + 0.5));
    }
    return value_[index];
  }

  void Clear() {
    constexpr unsigned int n =
        TRANSMITTANCE_TEXTURE_WIDTH * TRANSMITTANCE_TEXTURE_HEIGHT;
    for (unsigned int i = 0; i < n; ++i) {
      value_[i] = DimensionlessSpectrum(-1.0);
    }
  }

 private:
  const AtmosphereParameters& atmosphere_parameters_;
};

/*
<p>We also need a lazy single scattering texture:
*/

class LazySingleScatteringTexture :
    public dimensional::TernaryFunction<SCATTERING_TEXTURE_WIDTH,
                                        SCATTERING_TEXTURE_HEIGHT,
                                        SCATTERING_TEXTURE_DEPTH,
                                        IrradianceSpectrum> {
 public:
  LazySingleScatteringTexture(
      const AtmosphereParameters& atmosphere_parameters,
      const TransmittanceTexture& transmittance_texture,
      bool rayleigh)
      : TernaryFunction(IrradianceSpectrum(-watt_per_square_meter_per_nm)),
        atmosphere_parameters_(atmosphere_parameters),
        transmittance_texture_(transmittance_texture),
        rayleigh_(rayleigh) {
  }

  virtual const IrradianceSpectrum& Get(int i, int j, int k) const {
    int index =
        i + SCATTERING_TEXTURE_WIDTH * (j + SCATTERING_TEXTURE_HEIGHT * k);
    if (value_[index][0] < 0.0 * watt_per_square_meter_per_nm) {
      IrradianceSpectrum rayleigh;
      IrradianceSpectrum mie;
      ComputeSingleScatteringTexture(atmosphere_parameters_,
          transmittance_texture_, vec3(i + 0.5, j + 0.5, k + 0.5),
          rayleigh, mie);
      value_[index] = rayleigh_ ? rayleigh : mie;
    }
    return value_[index];
  }

 private:
  const AtmosphereParameters& atmosphere_parameters_;
  const TransmittanceTexture& transmittance_texture_;
  const bool rayleigh_;
};

/*
<p>a lazy multiple scattering texture, for step 1:
*/

class LazyScatteringDensityTexture :
    public dimensional::TernaryFunction<SCATTERING_TEXTURE_WIDTH,
                                        SCATTERING_TEXTURE_HEIGHT,
                                        SCATTERING_TEXTURE_DEPTH,
                                        RadianceDensitySpectrum> {
 public:
  LazyScatteringDensityTexture(
      const AtmosphereParameters& atmosphere_parameters,
      const TransmittanceTexture& transmittance_texture,
      const ReducedScatteringTexture& single_rayleigh_scattering_texture,
      const ReducedScatteringTexture& single_mie_scattering_texture,
      const ScatteringTexture& multiple_scattering_texture,
      const IrradianceTexture& irradiance_texture,
      const int order)
      : TernaryFunction(
            RadianceDensitySpectrum(-watt_per_cubic_meter_per_sr_per_nm)),
        atmosphere_parameters_(atmosphere_parameters),
        transmittance_texture_(transmittance_texture),
        single_rayleigh_scattering_texture_(single_rayleigh_scattering_texture),
        single_mie_scattering_texture_(single_mie_scattering_texture),
        multiple_scattering_texture_(multiple_scattering_texture),
        irradiance_texture_(irradiance_texture),
        order_(order) {
  }

  virtual const RadianceDensitySpectrum& Get(int i, int j, int k) const {
    int index =
        i + SCATTERING_TEXTURE_WIDTH * (j + SCATTERING_TEXTURE_HEIGHT * k);
    if (value_[index][0] < 0.0 * watt_per_cubic_meter_per_sr_per_nm) {
      value_[index] = ComputeScatteringDensityTexture(
          atmosphere_parameters_, transmittance_texture_,
          single_rayleigh_scattering_texture_, single_mie_scattering_texture_,
          multiple_scattering_texture_, irradiance_texture_,
          vec3(i + 0.5, j + 0.5, k + 0.5), order_);
    }
    return value_[index];
  }

 private:
  const AtmosphereParameters& atmosphere_parameters_;
  const TransmittanceTexture& transmittance_texture_;
  const ReducedScatteringTexture& single_rayleigh_scattering_texture_;
  const ReducedScatteringTexture& single_mie_scattering_texture_;
  const ScatteringTexture& multiple_scattering_texture_;
  const IrradianceTexture& irradiance_texture_;
  const int order_;
};


/*
<p>and step 2 of the multiple scattering computations:
*/

class LazyMultipleScatteringTexture :
    public dimensional::TernaryFunction<SCATTERING_TEXTURE_WIDTH,
                                        SCATTERING_TEXTURE_HEIGHT,
                                        SCATTERING_TEXTURE_DEPTH,
                                        RadianceSpectrum> {
 public:
  LazyMultipleScatteringTexture(
      const AtmosphereParameters& atmosphere_parameters,
      const TransmittanceTexture& transmittance_texture,
      const ScatteringDensityTexture& scattering_density_texture)
      : TernaryFunction(RadianceSpectrum(-watt_per_square_meter_per_sr_per_nm)),
        atmosphere_parameters_(atmosphere_parameters),
        transmittance_texture_(transmittance_texture),
        scattering_density_texture_(scattering_density_texture) {
  }

  virtual const RadianceSpectrum& Get(int i, int j, int k) const {
    int index =
        i + SCATTERING_TEXTURE_WIDTH * (j + SCATTERING_TEXTURE_HEIGHT * k);
    if (value_[index][0] < 0.0 * watt_per_square_meter_per_sr_per_nm) {
      Number ignored;
      value_[index] = ComputeMultipleScatteringTexture(atmosphere_parameters_,
          transmittance_texture_, scattering_density_texture_,
          vec3(i + 0.5, j + 0.5, k + 0.5), ignored);
    }
    return value_[index];
  }

 private:
  const AtmosphereParameters& atmosphere_parameters_;
  const TransmittanceTexture& transmittance_texture_;
  const ScatteringDensityTexture& scattering_density_texture_;
};

/*
<p>and, finally, a lazy ground irradiance texture:
*/

class LazyIndirectIrradianceTexture :
    public dimensional::BinaryFunction<IRRADIANCE_TEXTURE_WIDTH,
                                       IRRADIANCE_TEXTURE_HEIGHT,
                                       IrradianceSpectrum> {
 public:
  LazyIndirectIrradianceTexture(
      const AtmosphereParameters& atmosphere_parameters,
      const ReducedScatteringTexture& single_rayleigh_scattering_texture,
      const ReducedScatteringTexture& single_mie_scattering_texture,
      const ScatteringTexture& multiple_scattering_texture,
      int scattering_order)
      : BinaryFunction(IrradianceSpectrum(-watt_per_square_meter_per_nm)),
        atmosphere_parameters_(atmosphere_parameters),
        single_rayleigh_scattering_texture_(single_rayleigh_scattering_texture),
        single_mie_scattering_texture_(single_mie_scattering_texture),
        multiple_scattering_texture_(multiple_scattering_texture),
        scattering_order_(scattering_order) {
  }

  virtual const IrradianceSpectrum& Get(int i, int j) const {
    int index = i + j * IRRADIANCE_TEXTURE_HEIGHT;
    if (value_[index][0] < 0.0 * watt_per_square_meter_per_nm) {
      value_[index] = ComputeIndirectIrradianceTexture(atmosphere_parameters_,
          single_rayleigh_scattering_texture_,
          single_mie_scattering_texture_,
          multiple_scattering_texture_,
          vec2(i + 0.5, j + 0.5),
          scattering_order_);
    }
    return value_[index];
  }

 private:
  const AtmosphereParameters& atmosphere_parameters_;
  const ReducedScatteringTexture& single_rayleigh_scattering_texture_;
  const ReducedScatteringTexture& single_mie_scattering_texture_;
  const ScatteringTexture& multiple_scattering_texture_;
  int scattering_order_;
};

/*
<p>We can now define the unit tests themselves. Each test is an instance of the
following <code>TestCase</code> subclass, which has an
<code>atmosphere_parameters_</code> field initialized from the above constants.
Note that a new instance of this class is created for each unit test.
*/

class FunctionsTest : public dimensional::TestCase {
 public:
  template<typename T>
  FunctionsTest(const std::string& name, T test)
      : TestCase("FunctionsTest " + name, static_cast<Test>(test)) {
    atmosphere_parameters_.solar_irradiance[0] = kSolarIrradiance;
    atmosphere_parameters_.bottom_radius = kBottomRadius;
    atmosphere_parameters_.top_radius = kTopRadius;
    atmosphere_parameters_.rayleigh_density.layers[1] = DensityProfileLayer(
        0.0 * m, 1.0, -1.0 / kRayleighScaleHeight, 0.0 / m, 0.0);
    atmosphere_parameters_.rayleigh_scattering[0] = kRayleighScattering;
    atmosphere_parameters_.mie_density.layers[1] = DensityProfileLayer(
        0.0 * m, 1.0, -1.0 / kMieScaleHeight, 0.0 / m, 0.0);
    atmosphere_parameters_.mie_scattering[0] = kMieScattering;
    atmosphere_parameters_.mie_extinction[0] = kMieExtinction;
    atmosphere_parameters_.ground_albedo[0] = kGroundAlbedo;
    atmosphere_parameters_.mu_s_min = -1.0;
  }

/*
<p><i>Distance to the top atmosphere boundary</i>: check that this distance is
$r_\mathrm{top}-r$ for a vertical ray ($\mu=1$), and
$\sqrt{r_\mathrm{top}^2-r^2}$ for a horizontal ray ($\mu=0$).
*/

  void TestDistanceToTopAtmosphereBoundary() {
    constexpr Length r = kBottomRadius * 0.2 + kTopRadius * 0.8;
    // Vertical ray, looking top.
    ExpectNear(
        kTopRadius - r,
        DistanceToTopAtmosphereBoundary(atmosphere_parameters_, r, 1.0),
        1.0 * m);
    // Horizontal ray.
    ExpectNear(
        sqrt(kTopRadius * kTopRadius - r * r),
        DistanceToTopAtmosphereBoundary(atmosphere_parameters_, r, 0.0),
        1.0 * m);
  }

/*
<p><i>Intersections with the ground</i>: check that a vertical ray does not
intersect the ground, unless it is looking down. Likewise, check that a ray
looking slightly above the horizon ($\mu=\mu_{\mathrm{horiz}}+\epsilon$) does
not intersect the ground, but a ray looking slightly below the horizon
($\mu=\mu_{\mathrm{horiz}}-\epsilon$) does.
*/

  void TestRayIntersectsGround() {
    constexpr Length r = kBottomRadius * 0.9 + kTopRadius * 0.1;
    Number mu_horizon = CosineOfHorizonZenithAngle(r);
    ExpectFalse(RayIntersectsGround(atmosphere_parameters_, r, 1.0));
    ExpectFalse(RayIntersectsGround(
        atmosphere_parameters_, r, mu_horizon + kEpsilon));
    ExpectTrue(RayIntersectsGround(
        atmosphere_parameters_, r, mu_horizon - kEpsilon));
    ExpectTrue(RayIntersectsGround(atmosphere_parameters_, r, -1.0));
  }

/*
<p><i>Optical length to the top atmosphere boundary</i>: check that for a
vertical ray, looking up, the numerical integration in
<code>ComputeOpticalLengthToTopAtmosphereBoundary</code> gives the expected
result, which can be computed analytically:
$$
\int_r^{r_\mathrm{top}}
  \exp\left(-\frac{x-r_\mathrm{bottom}}{K}\right)\mathrm{d}x =
K\left[\exp\left(-\frac{r-r_\mathrm{bottom}}{K}\right)-
\exp\left(-\frac{r_\mathrm{top}-r_\mathrm{bottom}}{K}\right)\right]
$$
where $K$ is the scale height. Likewise, check that for a constant density
profile the optical length to the top atmosphere boundary is the distance to the
top atmosphere boundary (using a horizontal ray).
*/

  void TestComputeOpticalLengthToTopAtmosphereBoundary() {
    constexpr Length r = kBottomRadius * 0.2 + kTopRadius * 0.8;
    constexpr Length h_r = r - kBottomRadius;
    constexpr Length h_top = kTopRadius - kBottomRadius;
    // Vertical ray, looking top.
    ExpectNear(
        kRayleighScaleHeight * (exp(-h_r / kRayleighScaleHeight) -
            exp(-h_top / kRayleighScaleHeight)),
        ComputeOpticalLengthToTopAtmosphereBoundary(atmosphere_parameters_,
            atmosphere_parameters_.rayleigh_density, r, 1.0),
        1.0 * m);
    // Horizontal ray, no exponential density fall off.
    SetUniformAtmosphere();
    ExpectNear(
        sqrt(kTopRadius * kTopRadius - r * r),
        ComputeOpticalLengthToTopAtmosphereBoundary(atmosphere_parameters_,
            atmosphere_parameters_.rayleigh_density, r, 0.0),
        1.0 * m);
  }

/*
<p><i>Atmosphere density profiles</i>: check that density profiles with
exponentional, linear or constant density, and one or two layers, are correctly
computed.
*/

  void TestGetProfileDensity() {
    DensityProfile profile;
    // Only one layer, with exponentional density.
    profile.layers[1] =
        DensityProfileLayer(0.0 * m, 1.0, -1.0 / km, 0.0 / km, 0.0);
    ExpectEquals(exp(-2.0), GetProfileDensity(profile, 2.0 * km)());
    // Only one layer, with (clamped) affine density.
    profile.layers[1] =
        DensityProfileLayer(0.0 * m, 0.0, 0.0 / km, -0.5 / km, 1.0);
    ExpectEquals(1.0, GetProfileDensity(profile, 0.0 * km)());
    ExpectEquals(0.5, GetProfileDensity(profile, 1.0 * km)());
    ExpectEquals(0.0, GetProfileDensity(profile, 3.0 * km)());

    // Two layers, with (clamped) affine density.
    profile.layers[0] = DensityProfileLayer(
        25.0 * km, 0.0, 0.0 / km, 1.0 / (15.0 * km), -2.0 / 3.0);
    profile.layers[1] = DensityProfileLayer(
        0.0 * km, 0.0, 0.0 / km, -1.0 / (15.0 * km), 8.0 / 3.0);
    ExpectEquals(0.0, GetProfileDensity(profile, 0.0 * km)());
    ExpectNear(0.0, GetProfileDensity(profile, 10.0 * km)(), kEpsilon);
    ExpectNear(1.0, GetProfileDensity(profile, 25.0 * km)(), kEpsilon);
    ExpectNear(0.0, GetProfileDensity(profile, 40.0 * km)(), kEpsilon);
    ExpectEquals(0.0, GetProfileDensity(profile, 50.0 * km)());
  }

/*
<p><i>Transmittance to the top atmosphere boundary</i>: check that for a
vertical ray, looking up, the numerical integration in
<code>ComputeTransmittanceToTopAtmosphereBoundary</code> gives the expected
result, which can be computed analytically for an exponential profile (using the
above equation for the optical length) or for a triangular profile (such as the
one used for ozone). Likewise, check that for a horizontal ray, without Mie
scattering, and with a uniform density of air molecules, the optical depth is
the Rayleigh scattering coefficient times the distance to the top atmospere
boundary.
*/

  void TestComputeTransmittanceToTopAtmosphereBoundary() {
    constexpr Length r = kBottomRadius * 0.2 + kTopRadius * 0.8;
    constexpr Length h_r = r - kBottomRadius;
    constexpr Length h_top = kTopRadius - kBottomRadius;
    // Vertical ray, looking up.
    Number rayleigh_optical_depth = kRayleighScattering * kRayleighScaleHeight *
        (exp(-h_r / kRayleighScaleHeight) - exp(-h_top / kRayleighScaleHeight));
    Number mie_optical_depth = kMieExtinction * kMieScaleHeight *
        (exp(-h_r / kMieScaleHeight) - exp(-h_top / kMieScaleHeight));
    ExpectNear(
        exp(-(rayleigh_optical_depth + mie_optical_depth)),
        ComputeTransmittanceToTopAtmosphereBoundary(
            atmosphere_parameters_, r, 1.0)[0],
        Number(kEpsilon));
    // Vertical ray, looking up, no Rayleigh or Mie, only absorption, with a
    // triangular profile (whose integral is equal to 15km).
    atmosphere_parameters_.rayleigh_density.layers[1] = DensityProfileLayer();
    atmosphere_parameters_.mie_density.layers[1] = DensityProfileLayer();
    atmosphere_parameters_.absorption_density.layers[0] = DensityProfileLayer(
        25.0 * km, 0.0, 0.0 / km, 1.0 / (15.0 * km), -2.0 / 3.0);
    atmosphere_parameters_.absorption_density.layers[1] = DensityProfileLayer(
        0.0 * km, 0.0, 0.0 / km, -1.0 / (15.0 * km), 8.0 / 3.0);
    atmosphere_parameters_.absorption_extinction[0] = 0.02 / km;
    ExpectNear(
        exp(-Number(0.02 * 15.0)),
        ComputeTransmittanceToTopAtmosphereBoundary(
            atmosphere_parameters_, kBottomRadius, 1.0)[0],
        Number(kEpsilon));
    // Horizontal ray, uniform atmosphere without aerosols.
    SetUniformAtmosphere();
    RemoveAerosols();
    ExpectNear(
        exp(-kRayleighScattering * sqrt(kTopRadius * kTopRadius - r * r)),
        ComputeTransmittanceToTopAtmosphereBoundary(
            atmosphere_parameters_, r, 0.0)[0],
        Number(kEpsilon));
  }

/*
<p><i>Texture coordinates</i>: check that for a texture of size $n$, the center
of texel 0 (at texel coordinate $0.5/n$) is mapped to 0, and that the center of
texel $n-1$ (at texel coordinate $(n-0.5)/n$) is mapped to 1 (and vice-versa).
Finally, check that the mapping function and its inverse are really inverse of
each other (i.e. their composition should give the identity function).
*/

  void TestGetTextureCoordFromUnitRange() {
    ExpectNear(0.5 / 10.0, GetTextureCoordFromUnitRange(0.0, 10)(), kEpsilon);
    ExpectNear(9.5 / 10.0, GetTextureCoordFromUnitRange(1.0, 10)(), kEpsilon);

    ExpectNear(0.0, GetUnitRangeFromTextureCoord(0.5 / 10.0, 10)(), kEpsilon);
    ExpectNear(1.0, GetUnitRangeFromTextureCoord(9.5 / 10.0, 10)(), kEpsilon);

    ExpectNear(1.0 / 3.0, GetUnitRangeFromTextureCoord(
        GetTextureCoordFromUnitRange(1.0 / 3.0, 10), 10)(), kEpsilon);
  }

/*
<p><i>Mapping to transmittance texture coordinates</i>: check that the boundary
values of $r$ ($r_\mathrm{bottom}$ and $r_\mathrm{top}$) and $\mu$
($\mu_\mathrm{horiz}$ and $1$) are mapped to the centers of the boundary texels
of the transmittance texture.
*/

  void TestGetTransmittanceTextureUvFromRMu() {
    vec2 uv = GetTransmittanceTextureUvFromRMu(
        atmosphere_parameters_, kTopRadius, 1.0);
    ExpectNear(0.5 / TRANSMITTANCE_TEXTURE_WIDTH, uv.x(), kEpsilon);
    ExpectNear(1.0 - 0.5 / TRANSMITTANCE_TEXTURE_HEIGHT, uv.y(), kEpsilon);

    Number top_mu_horizon = CosineOfHorizonZenithAngle(kTopRadius);
    uv = GetTransmittanceTextureUvFromRMu(
        atmosphere_parameters_, kTopRadius, top_mu_horizon);
    ExpectNear(1.0 - 0.5 / TRANSMITTANCE_TEXTURE_WIDTH, uv.x(), kEpsilon);
    ExpectNear(1.0 - 0.5 / TRANSMITTANCE_TEXTURE_HEIGHT, uv.y(), kEpsilon);

    uv = GetTransmittanceTextureUvFromRMu(
        atmosphere_parameters_, kBottomRadius, 1.0);
    ExpectNear(0.5 / TRANSMITTANCE_TEXTURE_WIDTH, uv.x(), kEpsilon);
    ExpectNear(0.5 / TRANSMITTANCE_TEXTURE_HEIGHT, uv.y(), kEpsilon);

    uv = GetTransmittanceTextureUvFromRMu(
        atmosphere_parameters_, kBottomRadius, 0.0);
    ExpectNear(1.0 - 0.5 / TRANSMITTANCE_TEXTURE_WIDTH, uv.x(), kEpsilon);
    ExpectNear(0.5 / TRANSMITTANCE_TEXTURE_HEIGHT, uv.y(), kEpsilon);
  }

/*
<p><i>Mapping from transmittance texture coordinates</i>: check that the centers
of the boundary texels of the transmittance texture are mapped to the boundary
values of $r$ ($r_\mathrm{bottom}$ and $r_\mathrm{top}$) and $\mu$
($\mu_{\mathrm{horiz}}$ and $1$). Finally, check that the mapping function and
its inverse are really inverse of each other (i.e. their composition should give
the identity function).
*/

  void TestGetRMuFromTransmittanceTextureUv() {
    Length r;
    Number mu;
    GetRMuFromTransmittanceTextureUv(
        atmosphere_parameters_,
        vec2(0.5 / TRANSMITTANCE_TEXTURE_WIDTH,
             1.0 - 0.5 / TRANSMITTANCE_TEXTURE_HEIGHT),
        r, mu);
    ExpectNear(kTopRadius, r, 1.0 * m);
    ExpectNear(1.0, mu(), kEpsilon);

    GetRMuFromTransmittanceTextureUv(
        atmosphere_parameters_,
        vec2(1.0 - 0.5 / TRANSMITTANCE_TEXTURE_WIDTH,
             1.0 - 0.5 / TRANSMITTANCE_TEXTURE_HEIGHT),
        r, mu);
    ExpectNear(kTopRadius, r, 1.0 * m);
    ExpectNear(
        CosineOfHorizonZenithAngle(kTopRadius)(),
        mu(),
        kEpsilon);

    GetRMuFromTransmittanceTextureUv(
        atmosphere_parameters_,
        vec2(0.5 / TRANSMITTANCE_TEXTURE_WIDTH,
             0.5 / TRANSMITTANCE_TEXTURE_HEIGHT),
        r, mu);
    ExpectNear(kBottomRadius, r, 1.0 * m);
    ExpectNear(1.0, mu(), kEpsilon);

    GetRMuFromTransmittanceTextureUv(
        atmosphere_parameters_,
        vec2(1.0 - 0.5 / TRANSMITTANCE_TEXTURE_WIDTH,
             0.5 / TRANSMITTANCE_TEXTURE_HEIGHT),
        r, mu);
    ExpectNear(kBottomRadius, r, 1.0 * m);
    ExpectNear(0.0, mu(), kEpsilon);

    GetRMuFromTransmittanceTextureUv(
        atmosphere_parameters_,
        GetTransmittanceTextureUvFromRMu(atmosphere_parameters_,
            kBottomRadius * 0.2 + kTopRadius * 0.8, 0.25),
        r, mu);
    ExpectNear(kBottomRadius * 0.2 + kTopRadius * 0.8, r, 1.0 * m);
    ExpectNear(0.25, mu(), kEpsilon);
  }

/*
<p><i>Transmittance texture</i>: check that we get the same transmittance to the
top atmosphere boundary (more or less $\epsilon$) whether we compute it directly
with <code>ComputeTransmittanceToTopAtmosphereBoundary</code>, or via a
bilinearly interpolated lookup in the precomputed transmittance texture.
*/

  void TestGetTransmittanceToTopAtmosphereBoundary() {
    LazyTransmittanceTexture transmittance_texture(atmosphere_parameters_);

    const Length r = kBottomRadius * 0.2 + kTopRadius * 0.8;
    const Number mu = 0.4;
    ExpectNear(
        ComputeTransmittanceToTopAtmosphereBoundary(
            atmosphere_parameters_, r, mu)[0],
        GetTransmittanceToTopAtmosphereBoundary(
            atmosphere_parameters_, transmittance_texture, r, mu)[0],
        Number(kEpsilon));
  }

/*
<p><i>Transmittance texture</i>: check that <code>GetTransmittance</code> (which
combines two bilinearly interpolated lookups in the precomputed transmittance
texture) gives the expected result, in some cases where this result can be
computed analytically (when there are no aerosols and when the density of air
molecules is uniform, the optical depth is simply the Rayleigh scattering
coefficient times the distance travelled in the atmosphere).
*/

  void TestComputeAndGetTransmittance() {
    SetUniformAtmosphere();
    RemoveAerosols();
    LazyTransmittanceTexture transmittance_texture(atmosphere_parameters_);

    const Length r = kBottomRadius * 0.2 + kTopRadius * 0.8;
    const Length d = (kTopRadius - kBottomRadius) * 0.1;
    // Horizontal ray, from bottom atmosphere boundary.
    ExpectNear(
        exp(-kRayleighScattering * d),
        GetTransmittance(atmosphere_parameters_, transmittance_texture,
            kBottomRadius, 0.0, d, false /* ray_intersects_ground */)[0],
        Number(kEpsilon));
    // Almost vertical ray, looking up.
    ExpectNear(
        exp(-kRayleighScattering * d),
        GetTransmittance(
            atmosphere_parameters_, transmittance_texture, r, 0.7, d,
            false /* ray_intersects_ground */)[0],
        Number(kEpsilon));
    // Almost vertical ray, looking down.
    ExpectNear(
        exp(-kRayleighScattering * d),
        GetTransmittance(
            atmosphere_parameters_, transmittance_texture, r, -0.7, d,
            RayIntersectsGround(atmosphere_parameters_, r, -0.7))[0],
        Number(kEpsilon));
  }

/*
<p><i>Single scattering integrand</i>: check that the computation in
<code>ComputeSingleScatteringIntegrand</code> (which uses a precomputed
transmittance texture) gives the expected result, in 3 cases where this result
can be computed analytically:
<ul>
<li>vertical ray, looking up from the ground, sun at the zenith, scattering at
$r$. The integrand is then the product of:
  <ul>
  <li>the transmittance from the bottom to the top of the atmosphere. This
  involves the Rayleigh and Mie optical depths from the bottom to the top of the
  atmosphere, which have the form
  $$k_{ext}\int_{r_\mathrm{bottom}}^{r_\mathrm{top}}
      \exp\left(-\frac{x-r_\mathrm{bottom}}{K}\right)\mathrm{d}x =
      k_{ext}K\left[1-\exp\left(
          -\frac{r_\mathrm{top}-r_\mathrm{bottom}}{K}\right)\right]$$</li>
  <li>the (relative) number density at $r$:
  $\exp(-(r-r_\mathrm{bottom})/K))$.</li>
  </ul>
</li>
<li>vertical ray, looking down from the top atmosphere boundary, sun at the
zenith, scattering at $r$. The integrand is then the product of:
  <ul>
  <li>the transmittance from the top of the atmosphere to $r$, and back. This
  involves the Rayleigh and Mie optical depths from the top of the atmosphere to
  $r$ (times $2$), which have the form
  $$k_{ext}K\left[
          \exp\left(-\frac{r-r_\mathrm{bottom}}{K}\right)
          -\exp\left(-\frac{r_\mathrm{top}-r_\mathrm{bottom}}{K}\right)
      \right]$$</li>
  <li>the (relative) number density at $r$:
  $\exp(-(r-r_\mathrm{bottom})/K))$.</li>
  </ul>
</li>
<li>horizontal ray, from the ground, sun at the horizon, scattering at distance
$d$ from the viewer, without aerosols and with a uniform density of air
molecules. The integrand is then simply the transmittance for a horizontal ray
from the ground to the top atmosphere boundary (which has already been computed
in <code>TestComputeTransmittanceToTopAtmosphereBoundary</code>).</li>
</ul>
*/

  void TestComputeSingleScatteringIntegrand() {
    LazyTransmittanceTexture transmittance_texture(atmosphere_parameters_);

    // Vertical ray, from bottom to top atmosphere boundary, scattering at
    // middle of ray, scattering angle equal to 0.
    const Length h_top = kTopRadius - kBottomRadius;
    const Length h = h_top / 2.0;
    DimensionlessSpectrum rayleigh;
    DimensionlessSpectrum mie;
    ComputeSingleScatteringIntegrand(
        atmosphere_parameters_, transmittance_texture,
        kBottomRadius, 1.0, 1.0, 1.0, h, false, rayleigh, mie);
    Number rayleigh_optical_depth = kRayleighScattering * kRayleighScaleHeight *
        (1.0 - exp(-h_top / kRayleighScaleHeight));
    Number mie_optical_depth = kMieExtinction * kMieScaleHeight *
        (1.0 - exp(-h_top / kMieScaleHeight));
    ExpectNear(
        exp(-rayleigh_optical_depth - mie_optical_depth) *
            exp(-h / kRayleighScaleHeight),
        rayleigh[0],
        Number(kEpsilon));
    ExpectNear(
        exp(-rayleigh_optical_depth - mie_optical_depth) *
            exp(-h / kMieScaleHeight),
        mie[0],
        Number(kEpsilon));

    // Vertical ray, top to middle of atmosphere, scattering angle 180 degrees.
    ComputeSingleScatteringIntegrand(
        atmosphere_parameters_, transmittance_texture,
        kTopRadius, -1.0, 1.0, -1.0, h, true, rayleigh, mie);
    rayleigh_optical_depth = 2.0 * kRayleighScattering * kRayleighScaleHeight *
        (exp(-h / kRayleighScaleHeight) - exp(-h_top / kRayleighScaleHeight));
    mie_optical_depth = 2.0 * kMieExtinction * kMieScaleHeight *
        (exp(-h / kMieScaleHeight) - exp(-h_top / kMieScaleHeight));
    ExpectNear(
        exp(-rayleigh_optical_depth - mie_optical_depth) *
            exp(-h / kRayleighScaleHeight),
        rayleigh[0],
        Number(kEpsilon));
    ExpectNear(
        exp(-rayleigh_optical_depth - mie_optical_depth) *
            exp(-h / kMieScaleHeight),
        mie[0],
        Number(kEpsilon));

    // Horizontal ray, from bottom to top atmosphere boundary, scattering at
    // 50km, scattering angle equal to 0, uniform atmosphere, no aerosols.
    transmittance_texture.Clear();
    SetUniformAtmosphere();
    RemoveAerosols();
    ComputeSingleScatteringIntegrand(
        atmosphere_parameters_, transmittance_texture,
        kBottomRadius, 0.0, 0.0, 1.0, 50.0 * km, false, rayleigh, mie);
    rayleigh_optical_depth = kRayleighScattering * sqrt(
        kTopRadius * kTopRadius - kBottomRadius * kBottomRadius);
    ExpectNear(
        exp(-rayleigh_optical_depth),
        rayleigh[0],
        Number(kEpsilon));
  }

/*
<p><i>Distance to the nearest atmosphere boundary</i>: check that this distance
is $r_\mathrm{top}-r$ for a vertical ray looking up ($\mu=1$),
$\sqrt{r_\mathrm{top}^2-r^2}$ for a horizontal ray ($\mu=0$), and
$r-r_\mathrm{bottom}$ for a vertical ray looking down ($\mu=-1$).
*/

  void TestDistanceToNearestAtmosphereBoundary() {
    constexpr Length r = kBottomRadius * 0.2 + kTopRadius * 0.8;
    // Vertical ray, looking top.
    ExpectNear(
        kTopRadius - r,
        DistanceToNearestAtmosphereBoundary(
            atmosphere_parameters_, r, 1.0,
            RayIntersectsGround(atmosphere_parameters_, r, 1.0)),
        1.0 * m);
    // Horizontal ray.
    ExpectNear(
        sqrt(kTopRadius * kTopRadius - r * r),
        DistanceToNearestAtmosphereBoundary(
            atmosphere_parameters_, r, 0.0,
            RayIntersectsGround(atmosphere_parameters_, r, 0.0)),
        1.0 * m);
    // Vertical ray, looking down.
    ExpectNear(
        r - kBottomRadius,
        DistanceToNearestAtmosphereBoundary(
            atmosphere_parameters_, r, -1.0,
            RayIntersectsGround(atmosphere_parameters_, r, -1.0)),
        1.0 * m);
  }

/*
<p><i>Single scattering</i>: check that the numerical integration in
<code>ComputeSingleScattering</code> gives the expected result, in 2 cases where
the integral can be computed analytically:
<ul>
<li>vertical ray, looking up from the ground, sun at the zenith. The single
scattering integral has the form (where we omitted the solar irradiance and the
Mie terms to simplify the expression; see also
<code>TestComputeSingleScatteringIntegrand</code>):
$$
k_{sca}\int_{r_\mathrm{bottom}}^{r_\mathrm{top}}
    \exp\left(-k_{ext}K\left[1-
            \exp\left(-\frac{r_\mathrm{top}-r_\mathrm{bottom}}{K}\right)
        \right]\right)
    \exp\left(-\frac{r-r_\mathrm{bottom}}{K}\right)\mathrm{d}r
$$
which is equal to
$$
k_{sca}\exp\left(-k_{ext}K\left[1-
    \exp\left(-\frac{r_\mathrm{top}-r_\mathrm{bottom}}{K}\right)
        \right]\right)
K\left[
    1-\exp\left(-\frac{r_\mathrm{top}-r_\mathrm{bottom}}{K}\right)
\right]
$$
</li>
<li>vertical ray, looking down from the top atmosphere boundary, sun at the
zenith, no aerosols. The single scattering integral has the form (where we
omitted the solar irradiance; see also
<code>TestComputeSingleScatteringIntegrand</code>):
$$
k_{sca}\int_{r_\mathrm{bottom}}^{r_\mathrm{top}}
    \exp\left(-2k_{sca}K\left[
            \exp\left(-\frac{r-r_\mathrm{bottom}}{K}\right)-
            \exp\left(-\frac{r_\mathrm{top}-r_\mathrm{bottom}}{K}\right)
        \right]\right)\\
    \exp\left(-\frac{r-r_\mathrm{bottom}}{K}\right)\mathrm{d}r
$$
which, using the change of variables $u=\exp(-(r-r_\mathrm{bottom})/K)$,
gives
$$
\frac{1}{2}\left(1-\exp\left(-2k_{sca}K\left(
        1-\exp\left(-\frac{r_\mathrm{top}-r_\mathrm{bottom}}{K}\right)
    \right)\right)\right)
$$
</ul>
*/

  void TestComputeSingleScattering() {
    LazyTransmittanceTexture transmittance_texture(atmosphere_parameters_);

    // Vertical ray, from bottom atmosphere boundary, scattering angle 0.
    const Length h_top = kTopRadius - kBottomRadius;
    IrradianceSpectrum rayleigh;
    IrradianceSpectrum mie;
    ComputeSingleScattering(
        atmosphere_parameters_, transmittance_texture,
        kBottomRadius, 1.0, 1.0, 1.0, false, rayleigh, mie);
    Number rayleigh_optical_depth = kRayleighScattering * kRayleighScaleHeight *
        (1.0 - exp(-h_top / kRayleighScaleHeight));
    Number mie_optical_depth = kMieExtinction * kMieScaleHeight *
        (1.0 - exp(-h_top / kMieScaleHeight));
    // The relative error is about 1% here.
    ExpectNear(
        Number(1.0),
        rayleigh[0] / (kSolarIrradiance * rayleigh_optical_depth *
            exp(-rayleigh_optical_depth - mie_optical_depth)),
        Number(10.0 * kEpsilon));
    ExpectNear(
        Number(1.0),
        mie[0] / (kSolarIrradiance * mie_optical_depth * kMieScattering /
            kMieExtinction * exp(-rayleigh_optical_depth - mie_optical_depth)),
        Number(10.0 * kEpsilon));

    // Vertical ray, from top atmosphere boundary, scattering angle 180 degrees,
    // no aerosols.
    transmittance_texture.Clear();
    RemoveAerosols();
    ComputeSingleScattering(
        atmosphere_parameters_, transmittance_texture,
        kTopRadius, -1.0, 1.0, -1.0, true, rayleigh, mie);
    ExpectNear(
        Number(1.0),
        rayleigh[0] / (kSolarIrradiance *
            0.5 * (1.0 - exp(-2.0 * kRayleighScaleHeight * kRayleighScattering *
                (1.0 - exp(-h_top / kRayleighScaleHeight))))),
        Number(2.0 * kEpsilon));
    ExpectNear(0.0, mie[0].to(watt_per_square_meter_per_nm), kEpsilon);
  }

/*
<p><i>Rayleigh and Mie phase functions</i>: check that the integral of these
phase functions over all solid angles gives $1$.
*/

  void TestPhaseFunctions() {
    Number rayleigh_integral = 0.0;
    Number mie_integral = 0.0;
    const unsigned int N = 100;
    for (unsigned int i = 0; i < N; ++i) {
      Angle theta = (i + 0.5) * pi / N;
      SolidAngle domega = sin(theta) * (PI / N) * (2.0 * PI) * sr;
      rayleigh_integral =
          rayleigh_integral + RayleighPhaseFunction(cos(theta)) * domega;
      mie_integral =
          mie_integral + MiePhaseFunction(0.8, cos(theta)) * domega;
    }
    ExpectNear(1.0, rayleigh_integral(), 2.0 * kEpsilon);
    ExpectNear(1.0, mie_integral(), 2.0 * kEpsilon);
  }

/*
<p><i>Mapping to scattering texture coordinates</i>: check that the boundary
values of $r$ ($r_\mathrm{top}$, $r_\mathrm{bottom}$), $\mu$ ($-1$,
$\mu_{\mathrm{horiz}}$ and $1$), and $\mu_s$ and $\nu$ ($-1$ and $1$), are
mapped to the centers of the boundary texels of the scattering texture.
*/

  void TestGetScatteringTextureUvwzFromRMuMuSNu() {
    ExpectNear(
        0.5 / SCATTERING_TEXTURE_R_SIZE,
        GetScatteringTextureUvwzFromRMuMuSNu(
            atmosphere_parameters_, kBottomRadius, 0.0, 0.0, 0.0, false).w(),
        kEpsilon);
    ExpectNear(
        1.0 - 0.5 / SCATTERING_TEXTURE_R_SIZE,
        GetScatteringTextureUvwzFromRMuMuSNu(
            atmosphere_parameters_, kTopRadius, 0.0, 0.0, 0.0, false).w(),
        kEpsilon);

    const Length r = (kTopRadius + kBottomRadius) / 2.0;
    const Number mu = CosineOfHorizonZenithAngle(r);
    ExpectNear(
        0.5 / SCATTERING_TEXTURE_MU_SIZE,
        GetScatteringTextureUvwzFromRMuMuSNu(
            atmosphere_parameters_, r, mu, 0.0, 0.0, true).z(),
        kEpsilon);
    ExpectNear(
        1.0 - 0.5 / SCATTERING_TEXTURE_MU_SIZE,
        GetScatteringTextureUvwzFromRMuMuSNu(
            atmosphere_parameters_, r, mu, 0.0, 0.0, false).z(),
        kEpsilon);
    ExpectTrue(GetScatteringTextureUvwzFromRMuMuSNu(
        atmosphere_parameters_, r, -1.0, 0.0, 0.0, true).z() < 0.5);
    ExpectTrue(GetScatteringTextureUvwzFromRMuMuSNu(
        atmosphere_parameters_, r, 1.0, 0.0, 0.0, false).z() > 0.5);

    ExpectNear(
        0.5 / SCATTERING_TEXTURE_MU_S_SIZE,
        GetScatteringTextureUvwzFromRMuMuSNu(
            atmosphere_parameters_, kBottomRadius, 0.0, -1.0, 0.0, false).y(),
        kEpsilon);
    ExpectNear(
        1.0 - 0.5 / SCATTERING_TEXTURE_MU_S_SIZE,
        GetScatteringTextureUvwzFromRMuMuSNu(
            atmosphere_parameters_, kBottomRadius, 0.0, 1.0, 0.0, false).y(),
        kEpsilon);

    ExpectNear(
        0.5 / SCATTERING_TEXTURE_MU_S_SIZE,
        GetScatteringTextureUvwzFromRMuMuSNu(
            atmosphere_parameters_, kTopRadius, 0.0, -1.0, 0.0, false).y(),
        kEpsilon);
    ExpectNear(
        1.0 - 0.5 / SCATTERING_TEXTURE_MU_S_SIZE,
        GetScatteringTextureUvwzFromRMuMuSNu(
            atmosphere_parameters_, kTopRadius, 0.0, 1.0, 0.0, false).y(),
        kEpsilon);

    ExpectNear(
        0.0,
        GetScatteringTextureUvwzFromRMuMuSNu(
            atmosphere_parameters_, kBottomRadius, 0.0, 0.0, -1.0, false).x(),
        kEpsilon);
    ExpectNear(
        1.0,
        GetScatteringTextureUvwzFromRMuMuSNu(
            atmosphere_parameters_, kBottomRadius, 0.0, 0.0, 1.0, false).x(),
        kEpsilon);
  }

/*
<p><i>Mapping from scattering texture coordinates</i>: check that the centers of
the boundary texels of the scattering texture are mapped to the boundary values
of $r$ ($r_\mathrm{top}$, $r_\mathrm{bottom}$), $\mu$ ($-1$,
$\mu_{\mathrm{horiz}}$ and $1$), and $\mu_s$ and $\nu$ ($-1$ and $1$). Finally,
check that the mapping function and its inverse are really inverse of each other
(i.e. their composition should give the identity function).
*/

  void TestGetRMuMuSNuFromScatteringTextureUvwz() {
    Length r;
    Number mu;
    Number mu_s;
    Number nu;
    bool ray_r_mu_intersects_ground;
    GetRMuMuSNuFromScatteringTextureUvwz(atmosphere_parameters_,
        vec4(0.0,
             0.5 / SCATTERING_TEXTURE_MU_S_SIZE,
             0.5 / SCATTERING_TEXTURE_MU_SIZE,
             0.5 / SCATTERING_TEXTURE_R_SIZE),
        r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    ExpectNear(kBottomRadius, r, 1.0 * m);
    GetRMuMuSNuFromScatteringTextureUvwz(atmosphere_parameters_,
        vec4(0.0,
             0.5 / SCATTERING_TEXTURE_MU_S_SIZE,
             0.5 / SCATTERING_TEXTURE_MU_SIZE,
             1.0 - 0.5 / SCATTERING_TEXTURE_R_SIZE),
        r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    ExpectNear(kTopRadius, r, 1.0 * m);

    GetRMuMuSNuFromScatteringTextureUvwz(atmosphere_parameters_,
        vec4(0.0,
             0.5 / SCATTERING_TEXTURE_MU_S_SIZE,
             0.5 / SCATTERING_TEXTURE_MU_SIZE + kEpsilon,
             0.5),
        r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    const Number mu_horizon = CosineOfHorizonZenithAngle(r);
    ExpectNear(mu_horizon, mu, Number(kEpsilon));
    ExpectTrue(mu <= mu_horizon);
    ExpectTrue(ray_r_mu_intersects_ground);
    GetRMuMuSNuFromScatteringTextureUvwz(atmosphere_parameters_,
        vec4(0.0,
             0.5 / SCATTERING_TEXTURE_MU_S_SIZE,
             1.0 - 0.5 / SCATTERING_TEXTURE_MU_SIZE - kEpsilon,
             0.5),
        r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    ExpectNear(mu_horizon, mu, Number(5.0 * kEpsilon));
    ExpectTrue(mu >= mu_horizon);
    ExpectFalse(ray_r_mu_intersects_ground);

    GetRMuMuSNuFromScatteringTextureUvwz(atmosphere_parameters_,
        vec4(0.0,
             0.5 / SCATTERING_TEXTURE_MU_S_SIZE,
             0.5 / SCATTERING_TEXTURE_MU_SIZE,
             0.5 / SCATTERING_TEXTURE_R_SIZE),
        r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    ExpectNear(-1.0, mu_s(), kEpsilon);
    GetRMuMuSNuFromScatteringTextureUvwz(atmosphere_parameters_,
        vec4(0.0,
             1.0 - 0.5 / SCATTERING_TEXTURE_MU_S_SIZE,
             0.5 / SCATTERING_TEXTURE_MU_SIZE,
             0.5 / SCATTERING_TEXTURE_R_SIZE),
        r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    ExpectNear(1.0, mu_s(), kEpsilon);

    GetRMuMuSNuFromScatteringTextureUvwz(atmosphere_parameters_,
        vec4(0.0,
             0.5 / SCATTERING_TEXTURE_MU_S_SIZE,
             0.5 / SCATTERING_TEXTURE_MU_SIZE,
             0.5 / SCATTERING_TEXTURE_R_SIZE),
        r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    ExpectNear(-1.0, nu(), kEpsilon);
    GetRMuMuSNuFromScatteringTextureUvwz(atmosphere_parameters_,
        vec4(1.0,
             0.5 / SCATTERING_TEXTURE_MU_S_SIZE,
             0.5 / SCATTERING_TEXTURE_MU_SIZE,
             0.5 / SCATTERING_TEXTURE_R_SIZE),
        r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    ExpectNear(1.0, nu(), kEpsilon);

    GetRMuMuSNuFromScatteringTextureUvwz(atmosphere_parameters_,
        GetScatteringTextureUvwzFromRMuMuSNu(atmosphere_parameters_,
            kBottomRadius, -1.0, 1.0, -1.0, true),
        r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    ExpectNear(kBottomRadius, r, 1.0 * m);
    ExpectNear(-1.0, mu(), kEpsilon);
    ExpectNear(1.0, mu_s(), kEpsilon);
    ExpectNear(-1.0, nu(), kEpsilon);
    ExpectTrue(ray_r_mu_intersects_ground);

    GetRMuMuSNuFromScatteringTextureUvwz(atmosphere_parameters_,
        GetScatteringTextureUvwzFromRMuMuSNu(atmosphere_parameters_,
            kTopRadius, -1.0, 1.0, -1.0, true),
        r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    ExpectNear(kTopRadius, r, 1.0 * m);
    ExpectNear(-1.0, mu(), kEpsilon);
    ExpectNear(1.0, mu_s(), kEpsilon);
    ExpectNear(-1.0, nu(), kEpsilon);
    ExpectTrue(ray_r_mu_intersects_ground);

    GetRMuMuSNuFromScatteringTextureUvwz(atmosphere_parameters_,
        GetScatteringTextureUvwzFromRMuMuSNu(atmosphere_parameters_,
            (kBottomRadius + kTopRadius) / 2.0, 0.2, 0.3, 0.4, false),
        r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    ExpectNear((kBottomRadius + kTopRadius) / 2.0, r, 1.0 * m);
    ExpectNear(0.2, mu(), kEpsilon);
    ExpectNear(0.3, mu_s(), kEpsilon);
    ExpectNear(0.4, nu(), kEpsilon);
    ExpectFalse(ray_r_mu_intersects_ground);
  }

/*
<p><i>Single scattering texture</i>: check that we get the same single
scattering value (more or less $\epsilon$) whether we compute it directly
with <code>ComputeSingleScattering</code>, or via a
quadrilinearly interpolated lookup in the precomputed single scattering texture.
*/

  void TestComputeAndGetSingleScattering() {
    LazyTransmittanceTexture transmittance_texture(atmosphere_parameters_);
    LazySingleScatteringTexture single_rayleigh_scattering_texture(
        atmosphere_parameters_, transmittance_texture, true);
    LazySingleScatteringTexture single_mie_scattering_texture(
        atmosphere_parameters_, transmittance_texture, false);

    // Vertical ray, from bottom atmosphere boundary, scattering angle 0.
    IrradianceSpectrum rayleigh = GetScattering(
        atmosphere_parameters_, single_rayleigh_scattering_texture,
        kBottomRadius, 1.0, 1.0, 1.0, false);
    IrradianceSpectrum mie = GetScattering(
        atmosphere_parameters_, single_mie_scattering_texture,
        kBottomRadius, 1.0, 1.0, 1.0, false);
    IrradianceSpectrum expected_rayleigh;
    IrradianceSpectrum expected_mie;
    ComputeSingleScattering(
        atmosphere_parameters_, transmittance_texture,
        kBottomRadius, 1.0, 1.0, 1.0, false, expected_rayleigh, expected_mie);
    ExpectNear(1.0, (rayleigh / expected_rayleigh)[0](), kEpsilon);
    ExpectNear(1.0, (mie / expected_mie)[0](), kEpsilon);

    // Vertical ray, from top atmosphere boundary, scattering angle 180 degrees.
    rayleigh = GetScattering(
        atmosphere_parameters_, single_rayleigh_scattering_texture,
        kTopRadius, -1.0, 1.0, -1.0, true);
    mie = GetScattering(
        atmosphere_parameters_, single_mie_scattering_texture,
        kTopRadius, -1.0, 1.0, -1.0, true);
    ComputeSingleScattering(
        atmosphere_parameters_, transmittance_texture,
        kTopRadius, -1.0, 1.0, -1.0, true, expected_rayleigh, expected_mie);
    ExpectNear(1.0, (rayleigh / expected_rayleigh)[0](), kEpsilon);
    ExpectNear(1.0, (mie / expected_mie)[0](), kEpsilon);

    // Horizontal ray, from bottom of atmosphere, scattering angle 90 degrees.
    rayleigh = GetScattering(
        atmosphere_parameters_, single_rayleigh_scattering_texture,
        kBottomRadius, 0.0, 0.0, 0.0, false);
    mie = GetScattering(
        atmosphere_parameters_, single_mie_scattering_texture,
        kBottomRadius, 0.0, 0.0, 0.0, false);
    ComputeSingleScattering(
        atmosphere_parameters_, transmittance_texture,
        kBottomRadius, 0.0, 0.0, 0.0, false, expected_rayleigh, expected_mie);
    // The relative error is quite large in this case, i.e. between 6 to 8%.
    ExpectNear(1.0, (rayleigh / expected_rayleigh)[0](), 1e-1);
    ExpectNear(1.0, (mie / expected_mie)[0](), 1e-1);

    // Ray just above the horizon, sun at the zenith.
    Number mu = CosineOfHorizonZenithAngle(kTopRadius);
    rayleigh = GetScattering(
        atmosphere_parameters_, single_rayleigh_scattering_texture,
        kTopRadius, mu, 1.0, mu, false);
    mie = GetScattering(
        atmosphere_parameters_, single_mie_scattering_texture,
        kTopRadius, mu, 1.0, mu, false);
    ComputeSingleScattering(
        atmosphere_parameters_, transmittance_texture,
        kTopRadius, mu, 1.0, mu, false, expected_rayleigh, expected_mie);
    ExpectNear(1.0, (rayleigh / expected_rayleigh)[0](), kEpsilon);
    ExpectNear(1.0, (mie / expected_mie)[0](), kEpsilon);
  }

/*
<p><i>Multiple scattering, step 1</i>: check that the numerical integration in
<code>ComputeScatteringDensity</code> gives the expected result in two cases
where the integral can be computed analytically:
<ul>
<li>if the incident radiance from the $(n-1)$-th order is the same in all
directions, and if the ground albedo is $0$. In this case the scattered radiance
is the same in all directions, and equal to the scattering coefficient times the
incident radiance (the effect of the phase functions cancels out because we
integrate over all directions, and the integral of a phase function over all
directions is 1).</li>
<li>if the incident radiance from the $(n-1)$-th order is null, the ground
irradiance is the same everywhere, and the transmittance is 1. Then, at the
ground level, the radiance scattered in a horizontal direction is the scattering
coefficient times the ground irradiance, multiplied by the ground albedo,
divided by $\pi$ (the Lambertian BRDF), and divided by 2 (because the incident
radiance covers only the bottom hemisphere, and because of the symetries of the
phase functions and of the viewing conditions - especially the choice of a
horizontal direction for the scattered radiance).</li>
</ul>
*/

  void TestComputeScatteringDensity() {
    RadianceSpectrum kRadiance(13.0 * watt_per_square_meter_per_sr_per_nm);
    TransmittanceTexture full_transmittance(DimensionlessSpectrum(1.0));
    ReducedScatteringTexture no_single_scattering(
        IrradianceSpectrum(0.0 * watt_per_square_meter_per_nm));
    ScatteringTexture uniform_multiple_scattering(kRadiance);
    IrradianceTexture no_irradiance(
        IrradianceSpectrum(0.0 * watt_per_square_meter_per_nm));

    RadianceDensitySpectrum scattering_density = ComputeScatteringDensity(
        atmosphere_parameters_, full_transmittance,  no_single_scattering,
        no_single_scattering, uniform_multiple_scattering, no_irradiance,
        kBottomRadius, 0.0, 0.0, 1.0, 3);
    SpectralRadianceDensity kExpectedScatteringDensity =
        (kRayleighScattering + kMieScattering) * kRadiance[0];
    ExpectNear(
        1.0,
        (scattering_density[0] / kExpectedScatteringDensity)(),
        2.0 * kEpsilon);

    IrradianceSpectrum kIrradiance(13.0 * watt_per_square_meter_per_nm);
    IrradianceTexture uniform_irradiance(kIrradiance);
    ScatteringTexture no_multiple_scattering(
        RadianceSpectrum(0.0 * watt_per_square_meter_per_sr_per_nm));
    scattering_density = ComputeScatteringDensity(
        atmosphere_parameters_, full_transmittance, no_single_scattering,
        no_single_scattering, no_multiple_scattering, uniform_irradiance,
        kBottomRadius, 0.0, 0.0, 1.0, 3);
    kExpectedScatteringDensity = (kRayleighScattering + kMieScattering) *
        kGroundAlbedo / (2.0 * PI * sr) * kIrradiance[0];
    ExpectNear(
        1.0,
        (scattering_density[0] / kExpectedScatteringDensity)(),
        2.0 * kEpsilon);
  }

/*
<p><i>Multiple scattering, step 2</i>: check that the numerical integration in
<code>ComputeMultipleScattering</code> gives the expected result in some cases
where the integral can be computed analytically. If the radiance computed from
step 1 is the same everywhere and in all directions, and if the transmittance
is 1, then the numerical integration in <code>ComputeMultipleScattering</code>
should simply be equal to the radiance times the distance to the nearest
atmosphere boundary. We check this below for a vertical ray looking bottom, and
for a ray looking at the horizon.
*/

  void TestComputeMultipleScattering() {
    RadianceDensitySpectrum kRadianceDensity(
        0.17 * watt_per_cubic_meter_per_sr_per_nm);
    TransmittanceTexture full_transmittance(DimensionlessSpectrum(1.0));
    ScatteringDensityTexture uniform_scattering_density(kRadianceDensity);

    // Vertical ray, looking bottom.
    Length r = kBottomRadius * 0.2 + kTopRadius * 0.8;
    Length distance_to_ground = r - kBottomRadius;
    ExpectNear(
        kRadianceDensity[0] * distance_to_ground,
        ComputeMultipleScattering(atmosphere_parameters_, full_transmittance,
            uniform_scattering_density, r, -1.0, 1.0, -1.0, true)[0],
        kRadianceDensity[0] * distance_to_ground * kEpsilon);

    // Ray just below the horizon.
    Number mu = CosineOfHorizonZenithAngle(kTopRadius);
    Length distance_to_horizon =
        sqrt(kTopRadius * kTopRadius - kBottomRadius * kBottomRadius);
    ExpectNear(
        kRadianceDensity[0] * distance_to_horizon,
        ComputeMultipleScattering(atmosphere_parameters_, full_transmittance,
            uniform_scattering_density, kTopRadius, mu, 1.0, mu, true)[0],
        kRadianceDensity[0] * distance_to_horizon * kEpsilon);
  }

/*
<p><i>Multiple scattering texture, step 1</i>: check that we get the same result
for the first step of the multiple scattering computation, whether we compute
it directly, or via a quadrilinearly interpolated lookup in a precomputed
texture. For this, we use the same test cases as in
<code>TestComputeScatteringDensity</code>, where the end result can be computed
analytically.
*/

  void TestComputeAndGetScatteringDensity() {
    RadianceSpectrum kRadiance(13.0 * watt_per_square_meter_per_sr_per_nm);
    TransmittanceTexture full_transmittance(DimensionlessSpectrum(1.0));
    ReducedScatteringTexture no_single_scattering(
        IrradianceSpectrum(0.0 * watt_per_square_meter_per_nm));
    ScatteringTexture uniform_multiple_scattering(kRadiance);
    IrradianceTexture no_irradiance(
        IrradianceSpectrum(0.0 * watt_per_square_meter_per_nm));
    LazyScatteringDensityTexture multiple_scattering1(atmosphere_parameters_,
        full_transmittance, no_single_scattering, no_single_scattering,
        uniform_multiple_scattering, no_irradiance, 3);

    RadianceDensitySpectrum scattering_density = GetScattering(
        atmosphere_parameters_, multiple_scattering1,
        kBottomRadius, 0.0, 0.0, 1.0, false);
    SpectralRadianceDensity kExpectedScatteringDensity =
        (kRayleighScattering + kMieScattering) * kRadiance[0];
    ExpectNear(
        1.0,
        (scattering_density[0] / kExpectedScatteringDensity)(),
        2.0 * kEpsilon);

    IrradianceSpectrum kIrradiance(13.0 * watt_per_square_meter_per_nm);
    IrradianceTexture uniform_irradiance(kIrradiance);
    ScatteringTexture no_multiple_scattering(
        RadianceSpectrum(0.0 * watt_per_square_meter_per_sr_per_nm));

    LazyScatteringDensityTexture multiple_scattering2(atmosphere_parameters_,
        full_transmittance, no_single_scattering, no_single_scattering,
        no_multiple_scattering, uniform_irradiance, 3);
    scattering_density = GetScattering(
        atmosphere_parameters_, multiple_scattering2,
        kBottomRadius, 0.0, 0.0, 1.0, false);
    kExpectedScatteringDensity = (kRayleighScattering + kMieScattering) *
        kGroundAlbedo / (2.0 * PI * sr) * kIrradiance[0];
    ExpectNear(
        1.0,
        (scattering_density[0] / kExpectedScatteringDensity)(),
        2.0 * kEpsilon);
  }

/*
<p><i>Multiple scattering texture, step 2</i>: check that we get the same result
for the second step of the multiple scattering computation, whether we compute
it directly, or via a quadrilinearly interpolated lookup in a precomputed
texture. For this, we use the same test cases as in
<code>TestComputeMultipleScattering</code>, where the end result can be computed
analytically.
*/

  void TestComputeAndGetMultipleScattering() {
    RadianceDensitySpectrum kRadianceDensity(
        0.17 * watt_per_cubic_meter_per_sr_per_nm);
    TransmittanceTexture full_transmittance(DimensionlessSpectrum(1.0));
    ScatteringDensityTexture uniform_scattering_density(kRadianceDensity);
    LazyMultipleScatteringTexture multiple_scattering(atmosphere_parameters_,
        full_transmittance, uniform_scattering_density);

    // Vertical ray, looking bottom.
    Length r = kBottomRadius * 0.2 + kTopRadius * 0.8;
    Length distance_to_ground = r - kBottomRadius;
    ExpectNear(
        kRadianceDensity[0] * distance_to_ground,
        GetScattering(atmosphere_parameters_, multiple_scattering,
            r, -1.0, 1.0, -1.0, true)[0],
        kRadianceDensity[0] * distance_to_ground * kEpsilon);

    // Ray just below the horizon.
    Number mu = CosineOfHorizonZenithAngle(kTopRadius);
    Length distance_to_horizon =
        sqrt(kTopRadius * kTopRadius - kBottomRadius * kBottomRadius);
    ExpectNear(
        kRadianceDensity[0] * distance_to_horizon,
        GetScattering(atmosphere_parameters_, multiple_scattering,
            kTopRadius, mu, 1.0, mu, true)[0],
        kRadianceDensity[0] * distance_to_horizon * kEpsilon);
  }

/*
<p><i>Ground irradiance, indirect case</i>: check that the numerical integration
in <code>ComputeIndirectIrradiance</code> gives the expected result in a case
where this result can be computed analytically. More precisely, if the sky
radiance is the same everywhere and in all directions, then the ground
irradiance should be equal to $\pi$ times the sky radiance (because
$\int_0^{2\pi}\int_0^{\pi/2}\cos\theta\sin\theta\mathrm{d}\theta\mathrm{d}\phi=
\pi$, which is also why the Lambertian BRDF is $1/\pi$).
*/

  void TestComputeIndirectIrradiance() {
    ReducedScatteringTexture no_single_scattering;
    ScatteringTexture uniform_multiple_scattering(
        RadianceSpectrum(1.0 * watt_per_square_meter_per_sr_per_nm));
    IrradianceSpectrum irradiance = ComputeIndirectIrradiance(
        atmosphere_parameters_, no_single_scattering, no_single_scattering,
        uniform_multiple_scattering, kBottomRadius, 1.0, 2);
    // The relative error is about 1% here.
    ExpectNear(
        PI,
        irradiance[0].to(watt_per_square_meter_per_nm),
        10.0 * kEpsilon);
  }

/*
<p><i>Mapping to ground irradiance texture coordinates</i>: check that the
boundary values of $r$ ($r_\mathrm{bottom}$ and $r_\mathrm{top}$) and $\mu$
($-1$ and $1$) are mapped to the centers of the boundary texels of the ground
irradiance texture.
*/

  void TestGetIrradianceTextureUvFromRMuS() {
    ExpectNear(
        0.5 / IRRADIANCE_TEXTURE_HEIGHT,
        GetIrradianceTextureUvFromRMuS(
            atmosphere_parameters_, kBottomRadius, 0.0).y(),
        kEpsilon);
    ExpectNear(
        1.0 - 0.5 / IRRADIANCE_TEXTURE_HEIGHT,
        GetIrradianceTextureUvFromRMuS(
            atmosphere_parameters_, kTopRadius, 0.0).y(),
        kEpsilon);
    ExpectNear(
        0.5 / IRRADIANCE_TEXTURE_WIDTH,
        GetIrradianceTextureUvFromRMuS(
            atmosphere_parameters_, kBottomRadius, -1.0).x(),
        kEpsilon);
    ExpectNear(
        1.0 - 0.5 / IRRADIANCE_TEXTURE_WIDTH,
        GetIrradianceTextureUvFromRMuS(
            atmosphere_parameters_, kBottomRadius, 1.0).x(),
        kEpsilon);
  }

/*
<p><i>Mapping from ground irradiance texture coordinates</i>: check that the
centers of the boundary texels of the ground irradiance texture are mapped to
the boundary values of $r$ ($r_\mathrm{bottom}$ and $r_\mathrm{top}$) and $\mu$
($-1$ and $1$).
*/

  void TestGetRMuSFromIrradianceTextureUv() {
    Length r;
    Number mu_s;
    GetRMuSFromIrradianceTextureUv(atmosphere_parameters_,
        vec2(0.5 / IRRADIANCE_TEXTURE_WIDTH,
             0.5 / IRRADIANCE_TEXTURE_HEIGHT),
        r, mu_s);
    ExpectNear(kBottomRadius, r, 1.0 * m);
    GetRMuSFromIrradianceTextureUv(atmosphere_parameters_,
        vec2(0.5 / IRRADIANCE_TEXTURE_WIDTH,
             1.0 - 0.5 / IRRADIANCE_TEXTURE_HEIGHT),
        r, mu_s);
    ExpectNear(kTopRadius, r, 1.0 * m);
    GetRMuSFromIrradianceTextureUv(atmosphere_parameters_,
        vec2(0.5 / IRRADIANCE_TEXTURE_WIDTH,
             0.5 / IRRADIANCE_TEXTURE_HEIGHT),
        r, mu_s);
    ExpectNear(-1.0, mu_s(), kEpsilon);
    GetRMuSFromIrradianceTextureUv(atmosphere_parameters_,
        vec2(1.0 - 0.5 / IRRADIANCE_TEXTURE_WIDTH,
             0.5 / IRRADIANCE_TEXTURE_HEIGHT),
        r, mu_s);
    ExpectNear(1.0, mu_s(), kEpsilon);
  }

/*
<p><i>Ground irradiance texture</i>: check that we get the same ground
irradiance whether we compute it directly with
<code>ComputeIndirectIrradiance</code>, or via a bilinearly interpolated lookup
in the precomputed ground irradiance texture.
*/

  void TestComputeAndGetIrradiance() {
    ReducedScatteringTexture no_single_scattering(
        IrradianceSpectrum(0.0 * watt_per_square_meter_per_nm));
    ScatteringTexture fake_multiple_scattering;
    for (unsigned int x = 0; x < fake_multiple_scattering.size_x(); ++x) {
      for (unsigned int y = 0; y < fake_multiple_scattering.size_y(); ++y) {
        for (unsigned int z = 0; z < fake_multiple_scattering.size_z(); ++z) {
          double v = z + fake_multiple_scattering.size_z() *
              (y + fake_multiple_scattering.size_y() * x);
          fake_multiple_scattering.Set(x, y, z,
              RadianceSpectrum(v * watt_per_square_meter_per_sr_per_nm));
        }
      }
    }

    Length r = kBottomRadius * 0.8 + kTopRadius * 0.2;
    Number mu_s = 0.25;
    int scattering_order = 2;
    LazyIndirectIrradianceTexture irradiance_texture(atmosphere_parameters_,
        no_single_scattering, no_single_scattering, fake_multiple_scattering,
        scattering_order);
    ExpectNear(
        1.0,
        (GetIrradiance(atmosphere_parameters_, irradiance_texture, r, mu_s) /
            ComputeIndirectIrradiance(atmosphere_parameters_,
                no_single_scattering, no_single_scattering,
                fake_multiple_scattering, r, mu_s, scattering_order))[0](),
        kEpsilon);
    ExpectNotNear(
        1.0,
        (GetIrradiance(atmosphere_parameters_, irradiance_texture, r, mu_s) /
            ComputeIndirectIrradiance(atmosphere_parameters_,
                no_single_scattering, no_single_scattering,
                fake_multiple_scattering, r, 0.5, scattering_order))[0](),
        kEpsilon);
  }

/*
<p>And that's it for the unit tests! We just need to implement the two methods
that we used above to set a uniform density of air molecules and aerosols, and
to remove aerosols:
*/

 private:
  void SetUniformAtmosphere() {
    atmosphere_parameters_.rayleigh_density.layers[0] = DensityProfileLayer();
    atmosphere_parameters_.rayleigh_density.layers[1] =
        DensityProfileLayer(0.0 * m, 0.0, 0.0 / m, 0.0 / m, 1.0);
    atmosphere_parameters_.mie_density.layers[0] = DensityProfileLayer();
    atmosphere_parameters_.mie_density.layers[1] =
        DensityProfileLayer(0.0 * m, 0.0, 0.0 / m, 0.0 / m, 1.0);
    atmosphere_parameters_.absorption_density.layers[0] = DensityProfileLayer();
    atmosphere_parameters_.absorption_density.layers[1] = DensityProfileLayer();
  }

  void RemoveAerosols() {
    atmosphere_parameters_.mie_scattering[0] = 0.0 / km;
    atmosphere_parameters_.mie_extinction[0] = 0.0 / km;
  }

  AtmosphereParameters atmosphere_parameters_;
};

/*
<p>Finally, we need to create an instance of each of the above test cases (which
has the side effect of registering these instances in the <code>TestCase</code>
class, which can then run all the tests in its <code>RunAllTests</code> static
method).
*/

namespace {

FunctionsTest distance_to_top_atmosphere_boundary(
    "DistanceToTopAtmosphereBoundary",
    &FunctionsTest::TestDistanceToTopAtmosphereBoundary);
FunctionsTest ray_intersects_ground(
    "RayIntersectsGround",
    &FunctionsTest::TestRayIntersectsGround);
FunctionsTest get_profile_density(
    "GetProfileDensity",
    &FunctionsTest::TestGetProfileDensity);
FunctionsTest compute_optical_length_to_top_atmosphere_boundary(
    "ComputeOpticalLengthToTopAtmosphereBoundary",
    &FunctionsTest::TestComputeOpticalLengthToTopAtmosphereBoundary);
FunctionsTest compute_transmittance_to_top_atmosphere_boundary(
    "ComputeTransmittanceToTopAtmosphereBoundary",
    &FunctionsTest::TestComputeTransmittanceToTopAtmosphereBoundary);
FunctionsTest get_texture_coord_from_unit_range(
    "GetTextureCoordFromUnitRange",
    &FunctionsTest::TestGetTextureCoordFromUnitRange);
FunctionsTest get_transmittance_texture_uv_from_rmu(
    "GetTransmittanceTextureUvFromRMu",
    &FunctionsTest::TestGetTransmittanceTextureUvFromRMu);
FunctionsTest get_rmu_from_transmittance_texture_uv(
    "GetRMuFromTransmittanceTextureUv",
    &FunctionsTest::TestGetRMuFromTransmittanceTextureUv);
FunctionsTest get_transmittance_to_top_atmosphere_boundary(
    "GetTransmittanceToTopAtmosphereBoundary",
    &FunctionsTest::TestGetTransmittanceToTopAtmosphereBoundary);
FunctionsTest compute_and_get_transmittance(
    "ComputeAndGetTransmittance",
    &FunctionsTest::TestComputeAndGetTransmittance);

FunctionsTest compute_single_scattering_integrand(
    "ComputeSingleScatteringIntegrand",
    &FunctionsTest::TestComputeSingleScatteringIntegrand);
FunctionsTest distance_to_nearest_atmosphere_boundary(
    "DistanceToNearestAtmosphereBoundary",
    &FunctionsTest::TestDistanceToNearestAtmosphereBoundary);
FunctionsTest compute_single_scattering(
    "ComputeSingleScattering",
    &FunctionsTest::TestComputeSingleScattering);
FunctionsTest phase_functions(
    "PhaseFunctions",
    &FunctionsTest::TestPhaseFunctions);
FunctionsTest get_scattering_texture_uvwz_from_rmumusnu(
    "GetScatteringTextureUvwzFromRMuMuSNu",
    &FunctionsTest::TestGetScatteringTextureUvwzFromRMuMuSNu);
FunctionsTest get_rmumusnu_from_scattering_texture_uvwz(
    "GetRMuMuSNuFromScatteringTextureUvwz",
    &FunctionsTest::TestGetRMuMuSNuFromScatteringTextureUvwz);
FunctionsTest compute_and_get_scattering(
    "ComputeAndGetSingleScattering",
    &FunctionsTest::TestComputeAndGetSingleScattering);

FunctionsTest compute_scattering_density(
    "ComputeScatteringDensity",
    &FunctionsTest::TestComputeScatteringDensity);
FunctionsTest compute_multiple_scattering(
    "ComputeMultipleScattering",
    &FunctionsTest::TestComputeMultipleScattering);
FunctionsTest compute_and_get_scattering_density(
    "ComputeAndGetScatteringDensity",
    &FunctionsTest::TestComputeAndGetScatteringDensity);
FunctionsTest compute_and_get_multiple_scattering(
    "ComputeAndGetMultipleScattering",
    &FunctionsTest::TestComputeAndGetMultipleScattering);

FunctionsTest compute_indirect_irradiance(
    "ComputeIndirectIrradiance",
    &FunctionsTest::TestComputeIndirectIrradiance);
FunctionsTest get_irradiance_texture_uv_from_rmus(
    "GetIrradianceTextureUvFromRMuS",
    &FunctionsTest::TestGetIrradianceTextureUvFromRMuS);
FunctionsTest get_rmus_from_irradiance_texture_uv(
    "GetRMuSFromIrradianceTextureUv",
    &FunctionsTest::TestGetRMuSFromIrradianceTextureUv);
FunctionsTest get_irradiance(
    "GetComputeAndGetIrradiance",
    &FunctionsTest::TestComputeAndGetIrradiance);

}  // anonymous namespace

}  // namespace reference
}  // namespace atmosphere
