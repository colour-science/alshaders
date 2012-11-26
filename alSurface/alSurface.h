#pragma once

#include <ai.h>

struct ShaderData
{
   AtSampler* diffuse_sampler;
   AtSampler* glossy_sampler;
   AtSampler* glossy2_sampler;
   AtSampler* refraction_sampler;
   AtInt GI_diffuse_depth;
   AtInt GI_reflection_depth;
   AtInt GI_refraction_depth;
   AtInt GI_glossy_depth;
   AtInt GI_diffuse_samples;
   AtInt GI_glossy_samples;
   AtInt GI_refraction_samples;
   AtCritSec cs;
};
