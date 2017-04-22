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

/*<h2>atmosphere/reference/functions.h</h2>

<p>This file provides a C++ header for the <a href="../functions.glsl.html">GLSL
functions</a> that implement our atmosphere model. The C++ "implementation" is
provided in <a href="functions.cc.html">functions.cc</a> (this file simply
includes the GLSL file after defining the macros it depends on). The
documentation is provided in the GLSL file.
*/

#ifndef ATMOSPHERE_REFERENCE_FUNCTIONS_H_
#define ATMOSPHERE_REFERENCE_FUNCTIONS_H_

#include "atmosphere/reference/definitions.h"

namespace atmosphere {
namespace reference {

typedef dimensional::vec2 vec2;
typedef dimensional::vec3 vec3;
typedef dimensional::vec4 vec4;

// Transmittance.

Length DistanceToTopAtmosphereBoundary(
    const AtmosphereParameters& atmosphere, Length r, Number mu);

Length DistanceToBottomAtmosphereBoundary(
    const AtmosphereParameters& atmosphere, Length r, Number mu);

bool RayIntersectsGround(
    const AtmosphereParameters& atmosphere, Length r, Number mu);

Number GetLayerDensity(const DensityProfileLayer& layer, Length altitude);

Number GetProfileDensity(const DensityProfile& profile, Length altitude);

Length ComputeOpticalLengthToTopAtmosphereBoundary(
    const AtmosphereParameters& atmosphere, const DensityProfile& profile,
    Length r, Number mu);

DimensionlessSpectrum ComputeTransmittanceToTopAtmosphereBoundary(
    const AtmosphereParameters& atmosphere, Length r, Number mu);

Number GetTextureCoordFromUnitRange(Number x, int texture_size);

Number GetUnitRangeFromTextureCoord(Number u, int texture_size);

vec2 GetTransmittanceTextureUvFromRMu(const AtmosphereParameters& atmosphere,
    Length r, Number mu);

void GetRMuFromTransmittanceTextureUv(const AtmosphereParameters& atmosphere,
    const vec2& uv, Length& r, Number& mu);

DimensionlessSpectrum ComputeTransmittanceToTopAtmosphereBoundaryTexture(
    const AtmosphereParameters& atmosphere, const vec2& gl_frag_coord);

DimensionlessSpectrum GetTransmittanceToTopAtmosphereBoundary(
    const AtmosphereParameters& atmosphere,
    const TransmittanceTexture& transmittance_texture,
    Length r, Number mu);

DimensionlessSpectrum GetTransmittance(
    const AtmosphereParameters& atmosphere,
    const TransmittanceTexture& transmittance_texture,
    Length r, Number mu, Length d, bool ray_r_mu_intersects_ground);

// Single scattering.

void ComputeSingleScatteringIntegrand(
    const AtmosphereParameters& atmosphere,
    const TransmittanceTexture& transmittance_texture,
    Length r, Number mu, Number mu_s, Number nu, Length d,
    bool ray_r_mu_intersects_ground,
    DimensionlessSpectrum& rayleigh, DimensionlessSpectrum& mie);

Length DistanceToNearestAtmosphereBoundary(
    const AtmosphereParameters& atmosphere, Length r, Number mu,
    bool ray_r_mu_intersects_ground);

void ComputeSingleScattering(
    const AtmosphereParameters& atmosphere,
    const TransmittanceTexture& transmittance_texture,
    Length r, Number mu, Number mu_s, Number nu,
    bool ray_r_mu_intersects_ground,
    IrradianceSpectrum& rayleigh, IrradianceSpectrum& mie);

InverseSolidAngle RayleighPhaseFunction(Number nu);
InverseSolidAngle MiePhaseFunction(Number g, Number nu);

vec4 GetScatteringTextureUvwzFromRMuMuSNu(
    const AtmosphereParameters& atmosphere,
    Length r, Number mu, Number mu_s, Number nu,
    bool ray_r_mu_intersects_ground);

void GetRMuMuSNuFromScatteringTextureUvwz(
    const AtmosphereParameters& atmosphere, const vec4& uvwz,
    Length& r, Number& mu, Number& mu_s, Number& nu,
    bool& ray_r_mu_intersects_ground);

void ComputeSingleScatteringTexture(const AtmosphereParameters& atmosphere,
    const TransmittanceTexture& transmittance_texture,
    const vec3& gl_frag_coord, IrradianceSpectrum& rayleigh,
    IrradianceSpectrum& mie);

template<class T>
T GetScattering(
    const AtmosphereParameters& atmosphere,
    const AbstractScatteringTexture<T>& scattering_texture,
    Length r, Number mu, Number mu_s, Number nu,
    bool ray_r_mu_intersects_ground);

RadianceSpectrum GetScattering(
    const AtmosphereParameters& atmosphere,
    const ReducedScatteringTexture& single_rayleigh_scattering_texture,
    const ReducedScatteringTexture& single_mie_scattering_texture,
    const ScatteringTexture& multiple_scattering_texture,
    Length r, Number mu, Number mu_s, Number nu,
    bool ray_r_mu_intersects_ground,
    int scattering_order);

// Multiple scattering.

RadianceDensitySpectrum ComputeScatteringDensity(
    const AtmosphereParameters& atmosphere,
    const TransmittanceTexture& transmittance_texture,
    const ReducedScatteringTexture& single_rayleigh_scattering_texture,
    const ReducedScatteringTexture& single_mie_scattering_texture,
    const ScatteringTexture& multiple_scattering_texture,
    const IrradianceTexture& irradiance_texture,
    Length r, Number mu, Number mu_s, Number nu,
    int scattering_order);

RadianceSpectrum ComputeMultipleScattering(
    const AtmosphereParameters& atmosphere,
    const TransmittanceTexture& transmittance_texture,
    const ScatteringDensityTexture& scattering_density_texture,
    Length r, Number mu, Number mu_s, Number nu,
    bool ray_r_mu_intersects_ground);

RadianceDensitySpectrum ComputeScatteringDensityTexture(
    const AtmosphereParameters& atmosphere,
    const TransmittanceTexture& transmittance_texture,
    const ReducedScatteringTexture& single_rayleigh_scattering_texture,
    const ReducedScatteringTexture& single_mie_scattering_texture,
    const ScatteringTexture& multiple_scattering_texture,
    const IrradianceTexture& irradiance_texture,
    const vec3& gl_frag_coord, int scattering_order);

RadianceSpectrum ComputeMultipleScatteringTexture(
    const AtmosphereParameters& atmosphere,
    const TransmittanceTexture& transmittance_texture,
    const ScatteringDensityTexture& scattering_density_texture,
    const vec3& gl_frag_coord, Number& nu);

// Ground irradiance.

IrradianceSpectrum ComputeDirectIrradiance(
    const AtmosphereParameters& atmosphere,
    const TransmittanceTexture& transmittance_texture,
    Length r, Number mu_s);

IrradianceSpectrum ComputeIndirectIrradiance(
    const AtmosphereParameters& atmosphere,
    const ReducedScatteringTexture& single_rayleigh_scattering_texture,
    const ReducedScatteringTexture& single_mie_scattering_texture,
    const ScatteringTexture& multiple_scattering_texture,
    Length r, Number mu_s, int scattering_order);

vec2 GetIrradianceTextureUvFromRMuS(const AtmosphereParameters& atmosphere,
    Length r, Number mu_s);

void GetRMuSFromIrradianceTextureUv(const AtmosphereParameters& atmosphere,
    const vec2& uv, Length& r, Number& mu_s);

IrradianceSpectrum ComputeDirectIrradianceTexture(
    const AtmosphereParameters& atmosphere,
    const TransmittanceTexture& transmittance_texture,
    const vec2& gl_frag_coord);

IrradianceSpectrum ComputeIndirectIrradianceTexture(
    const AtmosphereParameters& atmosphere,
    const ReducedScatteringTexture& single_rayleigh_scattering_texture,
    const ReducedScatteringTexture& single_mie_scattering_texture,
    const ScatteringTexture& multiple_scattering_texture,
    const vec2& gl_frag_coord, int scattering_order);

IrradianceSpectrum GetIrradiance(
    const AtmosphereParameters& atmosphere,
    const IrradianceTexture& irradiance_texture,
    Length r, Number mu_s);

// Rendering.

RadianceSpectrum GetSkyRadiance(
    const AtmosphereParameters& atmosphere,
    const TransmittanceTexture& transmittance_texture,
    const ReducedScatteringTexture& scattering_texture,
    const ReducedScatteringTexture& single_mie_scattering_texture,
    Position camera, const Direction& view_ray, Length shadow_length,
    const Direction& sun_direction, DimensionlessSpectrum& transmittance);

RadianceSpectrum GetSkyRadianceToPoint(
    const AtmosphereParameters& atmosphere,
    const TransmittanceTexture& transmittance_texture,
    const ReducedScatteringTexture& scattering_texture,
    const ReducedScatteringTexture& single_mie_scattering_texture,
    Position camera, const Position& point, Length shadow_length,
    const Direction& sun_direction, DimensionlessSpectrum& transmittance);

IrradianceSpectrum GetSunAndSkyIrradiance(
    const AtmosphereParameters& atmosphere,
    const TransmittanceTexture& transmittance_texture,
    const IrradianceTexture& irradiance_texture,
    const Position& point, const Direction& normal,
    const Direction& sun_direction, IrradianceSpectrum& sky_irradiance);

}  // namespace reference
}  // namespace atmosphere

#endif  // ATMOSPHERE_REFERENCE_FUNCTIONS_H_
