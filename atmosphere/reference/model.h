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

/*<h2>atmosphere/reference/model.h</h2>

<p>This file defines the API to use our atmosphere model on CPU.
To use it:
<ul>
<li>create a <code>Model</code> instance with the desired atmosphere
parameters, and a directory where the precomputed textures can be cached,</li>
<li>call <code>Init</code> to precompute the atmosphere textures (or read
them from the cache directory if they have already been precomputed),</li>
<li>call <code>GetSolarRadiance</code>, <code>GetSkyRadiance</code>,
<code>GetSkyRadianceToPoint</code> and <code>GetSunAndSkyIrradiance</code> as
desired,</li>
<li>delete your <code>Model</code> when you no longer need it (the destructor
deletes the precomputed textures from memory).</li>
</ul>
*/

#ifndef ATMOSPHERE_REFERENCE_MODEL_H_
#define ATMOSPHERE_REFERENCE_MODEL_H_

#include <memory>
#include <string>
#include <vector>

#include "atmosphere/reference/definitions.h"

namespace atmosphere {
namespace reference {

class Model {
 public:
  Model(const AtmosphereParameters& atmosphere,
        const std::string& cache_directory);

  void Init(unsigned int num_scattering_orders = 4);

  RadianceSpectrum GetSolarRadiance() const;

  RadianceSpectrum GetSkyRadiance(Position camera, Direction view_ray,
      Length shadow_length, Direction sun_direction,
      DimensionlessSpectrum* transmittance) const;

  RadianceSpectrum GetSkyRadianceToPoint(Position camera, Position point,
      Length shadow_length, Direction sun_direction,
      DimensionlessSpectrum* transmittance) const;

  IrradianceSpectrum GetSunAndSkyIrradiance(Position p, Direction normal,
      Direction sun_direction, IrradianceSpectrum* sky_irradiance) const;

 private:
  const AtmosphereParameters atmosphere_;
  const std::string cache_directory_;
  std::unique_ptr<TransmittanceTexture> transmittance_texture_;
  std::unique_ptr<ReducedScatteringTexture> scattering_texture_;
  std::unique_ptr<ReducedScatteringTexture> single_mie_scattering_texture_;
  std::unique_ptr<IrradianceTexture> irradiance_texture_;
};

}  // namespace reference
}  // namespace atmosphere

#endif  // ATMOSPHERE_REFERENCE_MODEL_H_
