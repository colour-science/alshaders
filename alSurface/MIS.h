#pragma once
#include <ai.h>
#include "alUtil.h"

struct BrdfData_wrap
{
   void* brdf_data;
   AtShaderGlobals* sg;
   float eta;
   AtVector V;
   AtVector N;
   mutable float kr;
};


AtRGB AiWardDuerMISBRDF_wrap( const void* brdf_data, const AtVector* indir )
{
   AtVector H;
   const BrdfData_wrap* brdfw = reinterpret_cast<const BrdfData_wrap*>(brdf_data);
   AiV3Normalize(H,(*indir)+brdfw->V);
   float kr = fresnel(std::max(0.0f,AiV3Dot(H,*indir)),brdfw->eta);
   return kr *  AiWardDuerMISBRDF(brdfw->brdf_data, indir);
}

float AiWardDuerMISPDF_wrap( const void* brdf_data, const AtVector* indir )
{
   const BrdfData_wrap* brdfw = reinterpret_cast<const BrdfData_wrap*>(brdf_data);
   return AiWardDuerMISPDF(brdfw->brdf_data, indir);
}

AtVector AiWardDuerMISSample_wrap( const void* brdf_data, float randx, float randy )
{
   const BrdfData_wrap* brdfw = reinterpret_cast<const BrdfData_wrap*>(brdf_data);
   return AiWardDuerMISSample(brdfw->brdf_data, randx, randy);
}

AtRGB AiCookTorranceMISBRDF_wrap( const void* brdf_data, const AtVector* indir )
{
   AtVector H;
   const BrdfData_wrap* brdfw = reinterpret_cast<const BrdfData_wrap*>(brdf_data);
   AiV3Normalize(H,(*indir)+brdfw->V);
   brdfw->kr = fresnel(std::max(0.0f,AiV3Dot(H,*indir)),brdfw->eta);
   return brdfw->kr *  AiCookTorranceMISBRDF(brdfw->brdf_data, indir);
}

float AiCookTorranceMISPDF_wrap( const void* brdf_data, const AtVector* indir )
{
   const BrdfData_wrap* brdfw = reinterpret_cast<const BrdfData_wrap*>(brdf_data);
   return AiCookTorranceMISPDF(brdfw->brdf_data, indir);
}

AtVector AiCookTorranceMISSample_wrap( const void* brdf_data, float randx, float randy )
{
   const BrdfData_wrap* brdfw = reinterpret_cast<const BrdfData_wrap*>(brdf_data);
   return AiCookTorranceMISSample(brdfw->brdf_data, randx, randy);
}


AtRGB AiOrenNayarMISBRDF_wrap( const void* brdf_data, const AtVector* indir )
{
   AtVector H;
   const BrdfData_wrap* brdfw = reinterpret_cast<const BrdfData_wrap*>(brdf_data);
   AiV3Normalize(H,(*indir)+brdfw->V);
   //float kr = fresnel(std::max(0.0f,AiV3Dot(H,*indir)),brdfw->eta);
   float kr = fresnel(std::max(0.0f,AiV3Dot(brdfw->N,brdfw->V)),brdfw->eta);
   return AiOrenNayarMISBRDF(brdfw->brdf_data, indir) * (1-kr);
}

float AiOrenNayarMISPDF_wrap( const void* brdf_data, const AtVector* indir )
{
   const BrdfData_wrap* brdfw = reinterpret_cast<const BrdfData_wrap*>(brdf_data);
   return AiOrenNayarMISPDF(brdfw->brdf_data, indir);
}

AtVector AiOrenNayarMISSample_wrap( const void* brdf_data, float randx, float randy )
{
   const BrdfData_wrap* brdfw = reinterpret_cast<const BrdfData_wrap*>(brdf_data);
   return AiOrenNayarMISSample(brdfw->brdf_data, randx, randy);
}
