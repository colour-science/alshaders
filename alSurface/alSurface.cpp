#include <iostream>
#include <ai_sampler.h>
#include <map>
#include <cassert>
#include <cstdlib>
#include <vector>

#include "alUtil.h"
#include "MIS.h"
#include "alSurface.h"
#include "aovs.h"

AI_SHADER_NODE_EXPORT_METHODS(alSurfaceMtd)

#define GlossyMISBRDF AiCookTorranceMISBRDF
#define GlossyMISPDF AiCookTorranceMISPDF
#define GlossyMISSample AiCookTorranceMISSample

#define GlossyMISBRDF_wrap AiCookTorranceMISBRDF_wrap
#define GlossyMISPDF_wrap AiCookTorranceMISPDF_wrap
#define GlossyMISSample_wrap AiCookTorranceMISSample_wrap

#define GlossyMISCreateData AiCookTorranceMISCreateData

#define NUM_LIGHT_GROUPS 8

#define NUM_ID_AOVS 8
static const char* id_names[NUM_ID_AOVS] =
{
    "id_1",
    "id_2",
    "id_3",
    "id_4",
    "id_5",
    "id_6",
    "id_7",
    "id_8",
};

inline void flipNormals(AtShaderGlobals* sg)
{
    sg->Nf = -sg->Nf;
    sg->Ngf = -sg->Ngf;
}

enum alSurfaceParams
{
    // diffuse
    p_diffuseStrength=0,
    p_diffuseColor,
    p_diffuseRoughness,

    // backlight
    p_backlightStrength,
    p_backlightColor,
    p_backlightIndirectStrength,

    p_emissionStrength,
    p_emissionColor,

    // sss
    p_sssMix,
    p_sssRadius,
    p_sssRadiusColor,
    p_sssDensityScale,

    p_ssStrength,
    p_ssBalance,
    p_ssTargetColor,
    p_ssSpecifyCoefficients,
    p_ssScattering,
    p_ssAbsorption,
    p_ssDensityScale,
    p_ssDirection,
    p_ssInScattering,

    p_diffuseExtraSamples,
    p_diffuseEnableCaustics,
    p_diffuseIndirectStrength,

    // specular
    p_specular1Strength,
    p_specular1Color,
    p_specular1Roughness,
    p_specular1Ior,
    p_specular1RoughnessDepthScale,
    p_specular1ExtraSamples,
    p_specular1Normal,
    p_specular1IndirectStrength,
    p_specular1IndirectClamp,

    p_specular2Strength,
    p_specular2Color,
    p_specular2Roughness,
    p_specular2Ior,
    p_specular2RoughnessDepthScale,
    p_specular2ExtraSamples,
    p_specular2Normal,
    p_specular2IndirectStrength,
    p_specular2IndirectClamp,

    // transmission
    p_transmissionStrength,
    p_transmissionColor,
    p_transmissionLinkToSpecular1,
    p_transmissionRoughness,
    p_transmissionIor,
    p_transmissionRoughnessDepthScale,
    p_transmissionEnableCaustics,
    p_transmissionExtraSamples,
    p_transmissionClamp,

    p_lightGroupsIndirect,

    p_id1,
    p_id2,
    p_id3,
    p_id4,
    p_id5,
    p_id6,
    p_id7,
    p_id8,

    p_aiEnableMatte,
    p_aiMatteColor,
    p_aiMatteColorA,

    p_aov_diffuse_color,
    p_aov_direct_diffuse,
    p_aov_direct_diffuse_raw,
    p_aov_indirect_diffuse,
    p_aov_indirect_diffuse_raw,
    p_aov_direct_backlight,
    p_aov_indirect_backlight,
    p_aov_direct_specular,
    p_aov_indirect_specular,
    p_aov_direct_specular_2,
    p_aov_indirect_specular_2,
    p_aov_single_scatter,
    p_aov_sss,
    p_aov_refraction,
    p_aov_emission,
    p_aov_uv,
    p_aov_depth,
    p_aov_light_group_1,
    p_aov_light_group_2,
    p_aov_light_group_3,
    p_aov_light_group_4,
    p_aov_light_group_5,
    p_aov_light_group_6,
    p_aov_light_group_7,
    p_aov_light_group_8,
    p_aov_id_1,
    p_aov_id_2,
    p_aov_id_3,
    p_aov_id_4,
    p_aov_id_5,
    p_aov_id_6,
    p_aov_id_7,
    p_aov_id_8,

    p_standardAovs,
    p_transmitAovs,
    p_rrTransmission,
    p_rrTransmissionDepth,

    p_opacity,

    p_bump
};

node_parameters
{
    AiParameterFLT("diffuseStrength", 1.0f );
    AiParameterRGB("diffuseColor", 0.5f, 0.5f, 0.5f );
    AiParameterFLT("diffuseRoughness", 0.0f );

    AiParameterFLT("backlightStrength", 0.0f );
    AiParameterRGB("backlightColor", 0.18f, 0.18f, 0.18f );
    AiParameterFLT("backlightIndirectStrength", 1.0f);

    AiParameterFLT("emissionStrength", 0.0f );
    AiParameterRGB("emissionColor", 1.0f, 1.0f, 1.0f);

    AiParameterFLT("sssMix", 0.0f );
    AiParameterFLT("sssRadius", 3.6f );
    AiParameterRGB("sssRadiusColor", .439f, .156f, .078f );
    AiMetaDataSetBool(mds, "sssRadiusColor", "always_linear", true);  // no inverse-gamma correction
    AiParameterFLT("sssDensityScale", 1.0f );

    AiParameterFLT("ssStrength", 0.0f );
    AiParameterFLT("ssBalance", 0.5f);
    AiParameterRGB("ssTargetColor", .439f, .156f, .078f);
    AiParameterBOOL("ssSpecifyCoefficients", false);
    AiParameterRGB("ssScattering", 1.0f, 1.0f, 1.0f);
    AiParameterRGB("ssAbsorption", 1.0f, 1.0f, 1.0f);
    AiParameterFLT("ssDensityScale", 1.0f);
    AiParameterFLT("ssDirection", 0.0f);
    AiParameterBOOL("ssInScattering", true);

    AiParameterINT("diffuseExtraSamples", 0);
    AiParameterBOOL("diffuseEnableCaustics", false);
    AiParameterFLT("diffuseIndirectStrength", 1.0f);

    AiParameterFLT("specular1Strength", 1.0f );
    AiParameterRGB("specular1Color", 1.0f, 1.0f, 1.0f );
    AiParameterFLT("specular1Roughness", 0.3f );
    AiParameterFLT("specular1Ior", 1.4f );
    AiParameterFLT("specular1RoughnessDepthScale", 1.5f);
    AiParameterINT("specular1ExtraSamples", 0);
    AiParameterVec("specular1Normal", 0, 0, 0);
    AiParameterFLT("specular1IndirectStrength", 1.0f );
    AiParameterFLT("specular1IndirectClamp", 0.0f );

    AiParameterFLT("specular2Strength", 0.0f );
    AiParameterRGB("specular2Color", 1.0f, 1.0f, 1.0f );
    AiParameterFLT("specular2Roughness", 0.3f );
    AiParameterFLT("specular2Ior", 1.4f );
    AiParameterFLT("specular2RoughnessDepthScale", 1.5f);
    AiParameterINT("specular2ExtraSamples", 0);
    AiParameterVec("specular2Normal", 0, 0, 0);
    AiParameterFLT("specular2IndirectStrength", 1.0f );
    AiParameterFLT("specular2IndirectClamp", 0.0f );

    AiParameterFLT("transmissionStrength", 0.0f );
    AiParameterRGB("transmissionColor", 1.0f, 1.0f, 1.0f );
    AiParameterBOOL("transmissionLinkToSpecular1", true);
    AiParameterFLT("transmissionRoughness", 0.f );
    AiParameterFLT("transmissionIor", 1.4f );
    AiParameterFLT("transmissionRoughnessDepthScale", 1.5f);
    AiParameterBOOL("transmissionEnableCaustics", true);
    AiParameterINT("transmissionExtraSamples", 0);
    AiParameterFLT("transmissionClamp", 0.0f );

    AiParameterBOOL("lightGroupsIndirect", true);

    AiParameterRGB("id1", 0.0f, 0.0f, 0.0f);
    AiParameterRGB("id2", 0.0f, 0.0f, 0.0f);
    AiParameterRGB("id3", 0.0f, 0.0f, 0.0f);
    AiParameterRGB("id4", 0.0f, 0.0f, 0.0f);
    AiParameterRGB("id5", 0.0f, 0.0f, 0.0f);
    AiParameterRGB("id6", 0.0f, 0.0f, 0.0f);
    AiParameterRGB("id7", 0.0f, 0.0f, 0.0f);
    AiParameterRGB("id8", 0.0f, 0.0f, 0.0f);

    AiParameterBOOL("aiEnableMatte", false);
    AiParameterRGB("aiMatteColor", 0.0f, 0.0f, 0.0f);
    AiParameterFlt("aiMatteColorA", 0.0f);

    AiParameterStr("aov_diffuse_color", "diffuse_color");
    AiParameterStr("aov_direct_diffuse", "direct_diffuse");
    AiParameterStr("aov_direct_diffuse_raw", "direct_diffuse_raw");
    AiParameterStr("aov_indirect_diffuse", "indirect_diffuse");
    AiParameterStr("aov_indirect_diffuse_raw", "indirect_diffuse_raw");
    AiParameterStr("aov_direct_backlight", "direct_backlight");
    AiParameterStr("aov_indirect_backlight", "indirect_backlight");
    AiParameterStr("aov_direct_specular", "direct_specular");
    AiParameterStr("aov_indirect_specular", "indirect_specular");
    AiParameterStr("aov_direct_specular_2", "direct_specular_2");
    AiParameterStr("aov_indirect_specular_2", "indirect_specular_2");
    AiParameterStr("aov_single_scatter", "single_scatter");
    AiParameterStr("aov_sss", "sss");
    AiParameterStr("aov_refraction", "refraction");
    AiParameterStr("aov_emission", "emission");
    AiParameterStr("aov_uv", "uv");
    AiParameterStr("aov_depth", "depth");
    AiParameterStr("aov_light_group_1", "light_group_1");
    AiParameterStr("aov_light_group_2", "light_group_2");
    AiParameterStr("aov_light_group_3", "light_group_3");
    AiParameterStr("aov_light_group_4", "light_group_4");
    AiParameterStr("aov_light_group_5", "light_group_5");
    AiParameterStr("aov_light_group_6", "light_group_6");
    AiParameterStr("aov_light_group_7", "light_group_7");
    AiParameterStr("aov_light_group_8", "light_group_8");
    AiParameterStr("aov_id_1", "id_1");
    AiParameterStr("aov_id_2", "id_2");
    AiParameterStr("aov_id_3", "id_3");
    AiParameterStr("aov_id_4", "id_4");
    AiParameterStr("aov_id_5", "id_5");
    AiParameterStr("aov_id_6", "id_6");
    AiParameterStr("aov_id_7", "id_7");
    AiParameterStr("aov_id_8", "id_8");

    AiParameterBool("standardCompatibleAOVs", false);
    AiParameterBool("transmitAovs", false);

    AiParameterBool("rrTransmission", false);
    AiParameterInt("rrTransmissionDepth", 1);

    AiParameterRGB("opacity", 1.0f, 1.0f, 1.0f);
}

#ifdef MSVC
#define _CRT_SECURE_NO_WARNINGS 1
#endif
node_loader
{
    if (i>0) return 0;
    node->methods     = alSurfaceMtd;
    node->output_type = AI_TYPE_RGB;
    node->name        = "alSurface";
    node->node_type   = AI_NODE_SHADER;
    strcpy(node->version, AI_VERSION);
    return true;
}

node_initialize
{
    ShaderData* data = new ShaderData;
    AiNodeSetLocalData(node,data);
    data->diffuse_sampler = NULL;
    data->glossy_sampler = NULL;
    data->glossy2_sampler = NULL;
    data->refraction_sampler = NULL;
    data->backlight_sampler = NULL;
};

node_finish
{
    if (AiNodeGetLocalData(node))
    {
        ShaderData* data = (ShaderData*) AiNodeGetLocalData(node);

        AiSamplerDestroy(data->diffuse_sampler);
        AiSamplerDestroy(data->glossy_sampler);
        AiSamplerDestroy(data->glossy2_sampler);
        AiSamplerDestroy(data->refraction_sampler);
        AiSamplerDestroy(data->backlight_sampler);

        AiNodeSetLocalData(node, NULL);
        delete data;
    }
}


node_update
{
    ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);

    // set up AOVs
    REGISTER_AOVS

    data->standardAovs = params[p_standardAovs].BOOL;
    data->transmitAovs = params[p_transmitAovs].BOOL;

    // store some options we'll reuse later
    AtNode *options   = AiUniverseGetOptions();
    data->GI_diffuse_depth = AiNodeGetInt(options, "GI_diffuse_depth");
    data->GI_reflection_depth = AiNodeGetInt(options, "GI_reflection_depth");
    data->GI_refraction_depth = AiNodeGetInt(options, "GI_refraction_depth");
    data->GI_glossy_depth = AiNodeGetInt(options, "GI_glossy_depth");
    data->GI_glossy_samples = AiNodeGetInt(options, "GI_glossy_samples")+params[p_specular1ExtraSamples].INT;
    data->glossy_samples2 = SQR(data->GI_glossy_samples);
    data->GI_glossy2_samples = AiNodeGetInt(options, "GI_glossy_samples")+params[p_specular2ExtraSamples].INT;
    data->glossy2_samples2 = SQR(data->GI_glossy2_samples);
    data->GI_diffuse_samples = AiNodeGetInt(options, "GI_diffuse_samples")+params[p_diffuseExtraSamples].INT;
    data->diffuse_samples2 = SQR(data->GI_diffuse_samples);
    data->GI_refraction_samples = AiNodeGetInt(options, "GI_refraction_samples")+params[p_transmissionExtraSamples].INT;
    data->refraction_samples2 = SQR(data->GI_refraction_samples);

    // setup samples
    AiSamplerDestroy(data->diffuse_sampler);
    AiSamplerDestroy(data->glossy_sampler);
    AiSamplerDestroy(data->glossy2_sampler);
    AiSamplerDestroy(data->refraction_sampler);
    AiSamplerDestroy(data->backlight_sampler);
    data->diffuse_sampler = AiSampler(data->GI_diffuse_samples, 2);
    data->glossy_sampler = AiSampler(data->GI_glossy_samples, 2);
    data->glossy2_sampler = AiSampler(data->GI_glossy_samples, 2);
    data->refraction_sampler = AiSampler(data->GI_refraction_samples, 2);
    data->backlight_sampler = AiSampler(data->GI_diffuse_samples, 2);

    // Get all the light nodes in the scene and try and find their light group parameter
    // we'll store this based on the light pointer for fast access during rendering
    AtNodeIterator* it = AiUniverseGetNodeIterator(AI_NODE_LIGHT);
    while (!AiNodeIteratorFinished(it))
    {
        AtNode* light = AiNodeIteratorGetNext(it);
        if (AiNodeLookUpUserParameter(light, "lightGroup"))
            data->lightGroups[light] = AiNodeGetInt(light, "lightGroup") - 1;
        else 
            data->lightGroups[light] = -1;
    }
    AiNodeIteratorDestroy(it);
    data->lightGroupsIndirect = params[p_lightGroupsIndirect].BOOL;

    // check whether the normal parameters are connected or not
    data->specular1NormalConnected = AiNodeIsLinked(node, "specular1Normal");
    data->specular2NormalConnected = AiNodeIsLinked(node, "specular2Normal");

    data->rrTransmission = params[p_rrTransmission].BOOL;
    data->rrTransmissionDepth = params[p_rrTransmissionDepth].INT;

    data->specular1IndirectClamp = params[p_specular1IndirectClamp].FLT;
    if (data->specular1IndirectClamp == 0.0f) data->specular1IndirectClamp = AI_INFINITE;
    data->specular2IndirectClamp = params[p_specular2IndirectClamp].FLT;
    if (data->specular2IndirectClamp == 0.0f) data->specular2IndirectClamp = AI_INFINITE;
    data->transmissionClamp = params[p_transmissionClamp].FLT;
    if (data->transmissionClamp == 0.0f) data->transmissionClamp = AI_INFINITE;
};

inline void runAbsorption(){

}

inline void prelighting(DiffuseBRDF *diffuse)
{
    // Compute albedos in order of bxdf stack
    float accumulatedAlbedo = 0;

//    spec2.computeAlbedo();
//    accumulatedAlbedo += spec2.albedo;

//    if(accumulatedAlbedo < 1.f){
//        spec1.computeAlbedo();
//        spec1.albedo *= 1.f-accumulatedAlbedo;
//        accumulatedAlbedo += spec1.albedo;
//    }

    if(accumulatedAlbedo < 1.f){
        diffuseBRDFcomputeAlbedo(diffuse);
        diffuse->albedo *= 1.f-accumulatedAlbedo;
        accumulatedAlbedo += diffuse->albedo;
    }

//    if(accumulatedAlbedo < 1.f){
//        sss.computeAlbedo();
//        sss.albedo *= 1.f-accumulatedAlbedo;
//        accumulatedAlbedo += sss.albedo;
//    }

//    if(accumulatedAlbedo < 1.f){
//        transmission.computeAlbedo();
//        transmission.albedo *= 1.f-accumulatedAlbedo;
//        accumulatedAlbedo += transmission.albedo;
//    }

    // initialize mis data
    diffuseBRDFinitMIS(diffuse);
}

inline void directlighting(AtShaderGlobals *sg,
                           DiffuseBRDF *diffuse)
{
    // Light loop
    AiLightsPrepare(sg);
    while(AiLightsGetSample(sg))
    {
        if (diffuse->active)
        {
            diffuse->result += AiEvaluateLightSample(sg,
                                                    diffuse->brdf_data,
                                                    AiOrenNayarMISSample,
                                                    AiOrenNayarMISBRDF,
                                                    AiOrenNayarMISPDF)
                             * diffuse->strength;
        }
    }
}

inline void indirectlighting(AtShaderGlobals *sg,
                             ShaderData *data,
                             DiffuseBRDF *diffuse)
{
    // Sample BRDFS
#if AI_VERSION_MAJOR_NUM > 0
    float samples[2];
#else
    double samples[2];
#endif
    AtRay wi_ray;
    AtScrSample scrs;
    AtVector wi;

    // build a local frame for sampling
    AtVector U, V;
    AiBuildLocalFramePolar(&U, &V, &sg->N);

    // View direction, omega_o
    AtVector wo = -sg->Rd;

    if (diffuse->active)
    {
        AtSamplerIterator* sampit = AiSamplerIterator(data->diffuse_sampler, sg);
        AiMakeRay(&wi_ray, AI_RAY_DIFFUSE, &sg->P, NULL, AI_BIG, sg);
        while (AiSamplerGetSample(sampit, samples))
        {
            sampleHemisphereCosine(&wi, samples[0], samples[1]);
            float cos_theta = wi.z;
            if (cos_theta <= 0.0f) continue;

            AiV3RotateToFrame(wi, U, V, sg->Nf);
            float p = cos_theta * float(AI_ONEOVERPI);

            // trace the ray
            wi_ray.dir = wi;
            if (AiTrace(&wi_ray, &scrs))
            {
                AtRGB f = diffuse->strength
                        * AiOrenNayarMISBRDF(diffuse->brdf_data, &wi)
                        / p;
                diffuse->result += scrs.color * f;
            }
        }
    }
}

shader_evaluate
{
    ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);

    // Initialize brdfs
    DiffuseBRDF diffuse;
    diffuseBRDFinit(&diffuse, sg);

    // Shadow ray branch
//    if(raytype == shadow){
//        runAbsorption();
//        return;
//    }

    // Bump mapping

    // Bxdf prelighting
    prelighting(&diffuse);

    // Lighting
    directlighting(sg, &diffuse);

    // Indirect lighting
    indirectlighting(sg, data, &diffuse);

    // Bxdf Postlighting

    // AOVs

    // Final result
    sg->out.RGB =    diffuse.result * diffuse.albedo;
//                    +backlighting.result * backlighting.albedo
//                    +sss.result * sss.albedo
//                    +spec1.result * spec1.albedo
//                    +spec2.result * spec2.albedo
//                    +transmission.result * transmission.albedo
//                    +emission;
}
