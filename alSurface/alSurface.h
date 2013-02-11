#pragma once

#include <ai.h>
#include <map>

struct ShaderData
{
   AtSampler* diffuse_sampler;
   AtSampler* glossy_sampler;
   AtSampler* glossy2_sampler;
   AtSampler* refraction_sampler;
   int GI_diffuse_depth;
   int GI_reflection_depth;
   int GI_refraction_depth;
   int GI_glossy_depth;
   int GI_diffuse_samples;
   int GI_glossy_samples;
   int GI_refraction_samples;
   AtCritSec cs;
   std::map<AtNode*, int> lightGroups;
   bool specular1NormalConnected;
   bool specular2NormalConnected;
};
