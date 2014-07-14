#pragma once

#include <ai.h>
#include <map>
#include <string>
#include <vector>

#include "fresnel.h"

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
   std::map<AtNode*, float> shadowDensities;
   bool specular1NormalConnected;
   bool specular2NormalConnected;
   bool lightGroupsIndirect;
   bool standardAovs;
   bool transmitAovs;
   int numLights;
   bool rrTransmission;
   int rrTransmissionDepth;

   // data for doing RR
   int AA_samples;
   float AA_samples_inv;
   int total_depth;
   int* perm_table;
   int xres;

   float specular1IndirectClamp;
   float specular2IndirectClamp;
   float transmissionClamp;

   // AOV names
   std::vector<std::string> aovs;
   std::vector<std::string> aovs_rgba;

   // Fresnel
   Fresnel* fr1;
   Fresnel* fr2;

};
