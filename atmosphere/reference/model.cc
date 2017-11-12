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

/*<h2>atmosphere/reference/model.cc</h2>

<p>This file implements our atmosphere model on CPU. Its main role is to
precompute the transmittance, scattering and irradiance textures. The C++
functions to precompute them are provided in
<a href="functions.h.html">functions.h</a>, but they are not sufficient.
They must be called on each texel of the precomputed textures, and these
textures must be computed in the correct order, with the correct input an output
textures, to precompute each scattering order in sequence, as described in
Algorithm 4.1 of <a href="https://hal.inria.fr/inria-00288758/en">our paper</a>.
This is the role of the following code.

<p>We start by including the files we need:
*/

#include "atmosphere/reference/model.h"

#include "atmosphere/reference/functions.h"
#include "util/progress_bar.h"

/*
<p>The constructor of the <code>Model</code> class allocates the precomputed
textures, but does not initialize them.
*/

namespace atmosphere {
namespace reference {

Model::Model(const AtmosphereParameters& atmosphere,
             const std::string& cache_directory)
    : atmosphere_(atmosphere),
      cache_directory_(cache_directory) {
  transmittance_texture_.reset(new TransmittanceTexture());
  scattering_texture_.reset(new ReducedScatteringTexture());
  single_mie_scattering_texture_.reset(new ReducedScatteringTexture());
  irradiance_texture_.reset(new IrradianceTexture());
}

/*
<p>The initialization is done in the following method, which first tries to load
the textures from disk, if they have already been precomputed.
*/

void Model::Init(unsigned int num_scattering_orders) {
  std::ifstream file;
  file.open(cache_directory_ + "transmittance.dat");
  if (file.good()) {
    file.close();
    transmittance_texture_->Load(cache_directory_ + "transmittance.dat");
    scattering_texture_->Load(cache_directory_ + "scattering.dat");
    single_mie_scattering_texture_->Load(
        cache_directory_ + "single_mie_scattering.dat");
    irradiance_texture_->Load(cache_directory_ + "irradiance.dat");
    return;
  }

/*
<p>If they have not already been precomputed, we must compute them here. This
computation requires some temporary textures, in particular to store the
contribution of one scattering order, which is needed to compute the next order
of scattering (the final precomputed textures store the sum of all the
scattering orders). We allocate these textures here (they are automatically
destroyed at the end of this method).
*/

  std::unique_ptr<IrradianceTexture>
      delta_irradiance_texture(new IrradianceTexture());
  std::unique_ptr<ReducedScatteringTexture>
      delta_rayleigh_scattering_texture(new ReducedScatteringTexture());
  ReducedScatteringTexture* delta_mie_scattering_texture =
      single_mie_scattering_texture_.get();
  std::unique_ptr<ScatteringDensityTexture>
      delta_scattering_density_texture(new ScatteringDensityTexture());
  std::unique_ptr<ScatteringTexture>
      delta_multiple_scattering_texture(new ScatteringTexture());

/*
<p>Since the computation phase takes several minutes, we show a progress bar to
provide feedback to the user. The following constants roughly represent the
relative duration of each computation phase, and are used to display a progress
value which is roughly proportional to the elapsed time.
*/

  constexpr unsigned int kTransmittanceProgress = 1;
  constexpr unsigned int kDirectIrradianceProgress = 1;
  constexpr unsigned int kSingleScatteringProgress = 10;
  constexpr unsigned int kScatteringDensityProgress = 100;
  constexpr unsigned int kIndirectIrradianceProgress = 10;
  constexpr unsigned int kMultipleScatteringProgress = 10;
  const unsigned int kTotalProgress =
      TRANSMITTANCE_TEXTURE_WIDTH * TRANSMITTANCE_TEXTURE_HEIGHT *
          kTransmittanceProgress +
      IRRADIANCE_TEXTURE_WIDTH * IRRADIANCE_TEXTURE_HEIGHT * (
          kDirectIrradianceProgress +
          kIndirectIrradianceProgress * (num_scattering_orders - 1)) +
      SCATTERING_TEXTURE_WIDTH * SCATTERING_TEXTURE_HEIGHT *
          SCATTERING_TEXTURE_DEPTH * (
              kSingleScatteringProgress +
              (kScatteringDensityProgress + kMultipleScatteringProgress) *
                  (num_scattering_orders - 1));

  ProgressBar progress_bar(kTotalProgress);

/*
<p>The remaining code of this method implements Algorithm 4.1 of our paper,
using several threads to speed up computations (by computing several texels of
a texture in parallel).
*/

  // Compute the transmittance, and store it in transmittance_texture_.
  RunJobs([&](unsigned int j) {
    for (unsigned int i = 0; i < TRANSMITTANCE_TEXTURE_WIDTH; ++i) {
      transmittance_texture_->Set(i, j,
          ComputeTransmittanceToTopAtmosphereBoundaryTexture(
              atmosphere_, vec2(i + 0.5, j + 0.5)));
      progress_bar.Increment(kTransmittanceProgress);
    }
  }, TRANSMITTANCE_TEXTURE_HEIGHT);

  // Compute the direct irradiance, store it in delta_irradiance_texture, and
  // initialize irradiance_texture_ with zeros (we don't want the direct
  // irradiance in irradiance_texture_, but only the irradiance from the sky).
  RunJobs([&](unsigned int j) {
    for (unsigned int i = 0; i < IRRADIANCE_TEXTURE_WIDTH; ++i) {
      delta_irradiance_texture->Set(i, j,
          ComputeDirectIrradianceTexture(
              atmosphere_, *transmittance_texture_, vec2(i + 0.5, j + 0.5)));
      irradiance_texture_->Set(
          i, j, IrradianceSpectrum(0.0 * watt_per_square_meter_per_nm));
      progress_bar.Increment(kDirectIrradianceProgress);
    }
  }, IRRADIANCE_TEXTURE_HEIGHT);

  // Compute the rayleigh and mie single scattering, and store them in
  // delta_rayleigh_scattering_texture and delta_mie_scattering_texture, as well
  // as in scattering_texture.
  RunJobs([&](unsigned int k) {
    for (unsigned int j = 0; j < SCATTERING_TEXTURE_HEIGHT; ++j) {
      for (unsigned int i = 0; i < SCATTERING_TEXTURE_WIDTH; ++i) {
        IrradianceSpectrum rayleigh;
        IrradianceSpectrum mie;
        ComputeSingleScatteringTexture(atmosphere_, *transmittance_texture_,
            vec3(i + 0.5, j + 0.5, k + 0.5), rayleigh, mie);
        delta_rayleigh_scattering_texture->Set(i, j, k, rayleigh);
        delta_mie_scattering_texture->Set(i, j, k, mie);
        scattering_texture_->Set(i, j, k, rayleigh);
        progress_bar.Increment(kSingleScatteringProgress);
      }
    }
  }, SCATTERING_TEXTURE_DEPTH);

  // Compute the 2nd, 3rd and 4th order of scattering, in sequence.
  for (unsigned int scattering_order = 2;
       scattering_order <= num_scattering_orders;
       ++scattering_order) {
    // Compute the scattering density, and store it in
    // delta_scattering_density_texture.
    RunJobs([&](unsigned int k) {
      for (unsigned int j = 0; j < SCATTERING_TEXTURE_HEIGHT; ++j) {
        for (unsigned int i = 0; i < SCATTERING_TEXTURE_WIDTH; ++i) {
          RadianceDensitySpectrum scattering_density;
          scattering_density = ComputeScatteringDensityTexture(atmosphere_,
              *transmittance_texture_, *delta_rayleigh_scattering_texture,
              *delta_mie_scattering_texture,
              *delta_multiple_scattering_texture, *delta_irradiance_texture,
              vec3(i + 0.5, j + 0.5, k + 0.5), scattering_order);
          delta_scattering_density_texture->Set(i, j, k, scattering_density);
          progress_bar.Increment(kScatteringDensityProgress);
        }
      }
    }, SCATTERING_TEXTURE_DEPTH);

    // Compute the indirect irradiance, store it in delta_irradiance_texture and
    // accumulate it in irradiance_texture_.
    RunJobs([&](unsigned int j) {
      for (unsigned int i = 0; i < IRRADIANCE_TEXTURE_WIDTH; ++i) {
        IrradianceSpectrum delta_irradiance;
        delta_irradiance = ComputeIndirectIrradianceTexture(
            atmosphere_, *delta_rayleigh_scattering_texture,
            *delta_mie_scattering_texture, *delta_multiple_scattering_texture,
            vec2(i + 0.5, j + 0.5), scattering_order - 1);
        delta_irradiance_texture->Set(i, j, delta_irradiance);
        progress_bar.Increment(kIndirectIrradianceProgress);
      }
    }, IRRADIANCE_TEXTURE_HEIGHT);
    (*irradiance_texture_) += *delta_irradiance_texture;

    // Compute the multiple scattering, store it in
    // delta_multiple_scattering_texture, and accumulate it in
    // scattering_texture_.
    RunJobs([&](unsigned int k) {
      for (unsigned int j = 0; j < SCATTERING_TEXTURE_HEIGHT; ++j) {
        for (unsigned int i = 0; i < SCATTERING_TEXTURE_WIDTH; ++i) {
          RadianceSpectrum delta_multiple_scattering;
          Number nu;
          delta_multiple_scattering = ComputeMultipleScatteringTexture(
              atmosphere_, *transmittance_texture_,
              *delta_scattering_density_texture,
              vec3(i + 0.5, j + 0.5, k + 0.5), nu);
          delta_multiple_scattering_texture->Set(
              i, j, k, delta_multiple_scattering);
          scattering_texture_->Set(i, j, k,
              scattering_texture_->Get(i, j, k) +
              delta_multiple_scattering * (1.0 / RayleighPhaseFunction(nu)));
          progress_bar.Increment(kMultipleScatteringProgress);
        }
      }
    }, SCATTERING_TEXTURE_DEPTH);
  }

  transmittance_texture_->Save(cache_directory_ + "transmittance.dat");
  scattering_texture_->Save(cache_directory_ + "scattering.dat");
  single_mie_scattering_texture_->Save(
      cache_directory_ + "single_mie_scattering.dat");
  irradiance_texture_->Save(cache_directory_ + "irradiance.dat");
}

/*
<p>Once the textures have been computed or loaded from the cache, they can be
used to compute the sky radiance and the sun and sky irradiance. The functions
for doing that are provided in <a href="functions.h.html">functions.h</a> and we
just need here to wrap them in their corresponding methods (except for the solar
radiance, which can be directly computed from the model parameters):
*/

RadianceSpectrum Model::GetSolarRadiance() const {
  SolidAngle sun_solid_angle = 2.0 * PI *
      (1.0 - cos(atmosphere_.sun_angular_radius)) * sr;
  return atmosphere_.solar_irradiance * (1.0 / sun_solid_angle);
}

RadianceSpectrum Model::GetSkyRadiance(Position camera, Direction view_ray,
    Length shadow_length, Direction sun_direction,
    DimensionlessSpectrum* transmittance) const {
  return reference::GetSkyRadiance(atmosphere_, *transmittance_texture_,
      *scattering_texture_, *single_mie_scattering_texture_,
      camera, view_ray, shadow_length, sun_direction, *transmittance);
}

RadianceSpectrum Model::GetSkyRadianceToPoint(Position camera, Position point,
    Length shadow_length, Direction sun_direction,
    DimensionlessSpectrum* transmittance) const {
  return reference::GetSkyRadianceToPoint(atmosphere_, *transmittance_texture_,
      *scattering_texture_, *single_mie_scattering_texture_,
      camera, point, shadow_length, sun_direction, *transmittance);
}

IrradianceSpectrum Model::GetSunAndSkyIrradiance(Position point,
    Direction normal, Direction sun_direction,
    IrradianceSpectrum* sky_irradiance) const {
  return reference::GetSunAndSkyIrradiance(atmosphere_, *transmittance_texture_,
      *irradiance_texture_, point, normal, sun_direction, *sky_irradiance);
}

}  // namespace reference
}  // namespace atmosphere
