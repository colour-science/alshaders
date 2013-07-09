#pragma once

#include <ai.h>
#include <map>
#include <string>

struct ShaderData
{
   AtSampler* diffuse_sampler;
   AtSampler* glossy_sampler;
   AtSampler* glossy2_sampler;
   AtSampler* refraction_sampler;
   AtSampler* backlight_sampler;
   int GI_diffuse_depth;
   int GI_reflection_depth;
   int GI_refraction_depth;
   int GI_glossy_depth;
   int GI_diffuse_samples;
   int diffuse_samples2;
   int GI_glossy_samples;
   int glossy_samples2;
   int GI_glossy2_samples;
   int glossy2_samples2;
   int GI_refraction_samples;
   int refraction_samples2;
   int diffuse_sample_offset;
   int glossy_sample_offset;
   int glossy2_sample_offset;
   int refraction_sample_offset;
   int backlight_sample_offset;
   int total_samples;
   AtCritSec cs;
   std::map<AtNode*, int> lightGroups;
   bool specular1NormalConnected;
   bool specular2NormalConnected;
   bool lightGroupsIndirect;
   bool standardAovs;
   int numLights;

   // AOV names
   std::string aov_diffuse_color;
   std::string aov_direct_diffuse;
   std::string aov_direct_diffuse_raw;
   std::string aov_indirect_diffuse;
   std::string aov_indirect_diffuse_raw;
   std::string aov_direct_backlight;
   std::string aov_indirect_backlight;
   std::string aov_direct_specular;
   std::string aov_indirect_specular;
   std::string aov_direct_specular_2;
   std::string aov_indirect_specular_2;
   std::string aov_single_scatter;
   std::string aov_sss;
   std::string aov_refraction;
   std::string aov_emission;
   std::string aov_uv;
   std::string aov_depth;
   std::string aov_light_group[8];
   std::string aov_id[8];
};
