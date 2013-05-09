#include <iostream>
#include <ai_sampler.h>
#include <OpenEXR/ImathVec.h>
#include <OpenEXR/ImathMatrix.h>
#include <OpenEXR/ImathMatrixAlgo.h>
#include <map>

#include "alUtil.h"
#include "MIS.h"
#include "BeckmannMicrofacet.h"
#include "alSurface.h"

AI_SHADER_NODE_EXPORT_METHODS(alSurfaceMtd)

#define GlossyMISBRDF AiCookTorranceMISBRDF
#define GlossyMISPDF AiCookTorranceMISPDF
#define GlossyMISSample AiCookTorranceMISSample

#define GlossyMISBRDF_wrap AiCookTorranceMISBRDF_wrap
#define GlossyMISPDF_wrap AiCookTorranceMISPDF_wrap
#define GlossyMISSample_wrap AiCookTorranceMISSample_wrap

#define GlossyMISCreateData AiCookTorranceMISCreateData

#define NUM_LIGHT_GROUPS 8
static const char* lightGroupNames[NUM_LIGHT_GROUPS] =
{
    "light_group_1",
    "light_group_2",
    "light_group_3",
    "light_group_4",
    "light_group_5",
    "light_group_6",
    "light_group_7",
    "light_group_8"
};

#define NUM_ID_AOVS 8
static const char* idAovNames[NUM_ID_AOVS] = 
{
    "id_1",
    "id_2",
    "id_3",
    "id_4",
    "id_5",
    "id_6",
    "id_7",
    "id_8"
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

    // specular
    p_specular1Strength,
    p_specular1Color,
    p_specular1Roughness,
    p_specular1Ior,
    p_specular1RoughnessDepthScale,
    p_specular1ExtraSamples,
    p_specular1Normal,
    p_specular2Strength,
    p_specular2Color,
    p_specular2Roughness,
    p_specular2Ior,
    p_specular2RoughnessDepthScale,
    p_specular2ExtraSamples,
    p_specular2Normal,

    // transmission
    p_transmissionStrength,
    p_transmissionColor,
    p_transmissionLinkToSpecular1,
    p_transmissionRoughness,
    p_transmissionIor,
    p_transmissionRoughnessDepthScale,
    p_transmissionEnableCaustics,
    p_transmissionExtraSamples,

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

    p_bump
};

node_parameters
{
    AiParameterFLT("diffuseStrength", 1.0f );
    AiParameterRGB("diffuseColor", 0.18f, 0.18f, 0.18f );
    AiParameterFLT("diffuseRoughness", 0.0f );

    AiParameterFLT("backlightStrength", 0.0f );
    AiParameterRGB("backlightColor", 0.18f, 0.18f, 0.18f );

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

    AiParameterFLT("specular1Strength", 1.0f );
    AiParameterRGB("specular1Color", 1.0f, 1.0f, 1.0f );
    AiParameterFLT("specular1Roughness", 0.3f );
    AiParameterFLT("specular1Ior", 1.4f );
    AiParameterFLT("specular1RoughnessDepthScale", 1.5f);
    AiParameterINT("specular1ExtraSamples", 0);
    AiParameterVec("specular1Normal", 0, 0, 0);

    AiParameterFLT("specular2Strength", 0.0f );
    AiParameterRGB("specular2Color", 1.0f, 1.0f, 1.0f );
    AiParameterFLT("specular2Roughness", 0.3f );
    AiParameterFLT("specular2Ior", 1.4f );
    AiParameterFLT("specular2RoughnessDepthScale", 1.5f);
    AiParameterINT("specular2ExtraSamples", 0);
    AiParameterVec("specular2Normal", 0, 0, 0);

    AiParameterFLT("transmissionStrength", 0.0f );
    AiParameterRGB("transmissionColor", 1.0f, 1.0f, 1.0f );
    AiParameterBOOL("transmissionLinkToSpecular1", true);
    AiParameterFLT("transmissionRoughness", 0.1f );
    AiParameterFLT("transmissionIor", 1.4f );
    AiParameterFLT("transmissionRoughnessDepthScale", 1.5f);
    AiParameterBOOL("transmissionEnableCaustics", true);
    AiParameterINT("transmissionExtraSamples", 0);

    AiParameterBOOL("lightGroupsIndirect", false);

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

}


node_loader
{
    if (i>0) return 0;
    node->methods     = alSurfaceMtd;
    node->output_type = AI_TYPE_RGB;
    node->name        = "alSurface";
    node->node_type   = AI_NODE_SHADER;
    strcpy(node->version, AI_VERSION);
    return TRUE;
}

node_initialize
{
    ShaderData *data = new ShaderData;
    AiNodeSetLocalData(node,data);
    data->diffuse_sampler = NULL;
    data->glossy_sampler = NULL;
    data->glossy2_sampler = NULL;
    data->refraction_sampler = NULL;
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

        AiNodeSetLocalData(node, NULL);
        delete data;
    }
}


node_update
{
    // set up AOVs
    AiAOVRegister("diffuse_color", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("direct_diffuse", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("direct_backlight", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("direct_diffuse_raw", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("indirect_diffuse", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("indirect_backlight", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("indirect_diffuse_raw", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("direct_specular", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("indirect_specular", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("direct_specular_2", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("indirect_specular_2", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("single_scatter", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("sss", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("refraction", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("emission", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("uv", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("depth", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_1", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_2", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_3", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_4", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_5", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_6", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_7", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_8", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_1", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_2", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_3", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_4", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_5", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_6", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_7", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_8", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);

    // store some options we'll reuse later
    ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);
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

    // set up sample offsets for deep groups
    // order of execution is glossy->glossy2->diffuse->refraction
    data->glossy_sample_offset = 0;
    data->glossy2_sample_offset = data->glossy_sample_offset + data->glossy_samples2;
    data->diffuse_sample_offset = data->glossy2_sample_offset + data->glossy2_samples2;
    data->refraction_sample_offset = data->diffuse_sample_offset + data->diffuse_samples2;
    data->total_samples = data->refraction_sample_offset + data->refraction_samples2;

    // setup samples
    AiSamplerDestroy(data->diffuse_sampler);
    AiSamplerDestroy(data->glossy_sampler);
    AiSamplerDestroy(data->glossy2_sampler);
    AiSamplerDestroy(data->refraction_sampler);
    data->diffuse_sampler = AiSampler(data->GI_diffuse_samples, 2);
    data->glossy_sampler = AiSampler(data->GI_glossy_samples, 2);
    data->glossy2_sampler = AiSampler(data->GI_glossy_samples, 2);
    data->refraction_sampler = AiSampler(data->GI_refraction_samples, 2);

    // Get all the light nodes in the scene and try and find their light group parameter
    // we'll store this based on the light pointer for fast access during rendering
    AtNodeIterator* it = AiUniverseGetNodeIterator(AI_NODE_LIGHT);
    while (!AiNodeIteratorFinished(it))
    {
        AtNode* light = AiNodeIteratorGetNext(it);
        data->lightGroups[light] = AiNodeGetInt(light, "lightGroup") - 1;
    }
    AiNodeIteratorDestroy(it);
    data->lightGroupsIndirect = params[p_lightGroupsIndirect].BOOL;

    // check whether the normal parameters are connected or not
    data->specular1NormalConnected = AiNodeIsLinked(node, "specular1Normal");
    data->specular2NormalConnected = AiNodeIsLinked(node, "specular2Normal");
};


shader_evaluate
{
    ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);

    float ior = std::max(1.001f, AiShaderEvalParamFlt(p_specular1Ior));
    float roughness = AiShaderEvalParamFlt( p_specular1Roughness );
    roughness *= roughness;
    float transmissionRoughness;
    float transmissionIor;
    bool transmissionLinkToSpecular1 = AiShaderEvalParamBool(p_transmissionLinkToSpecular1);
    if (transmissionLinkToSpecular1)
    {
        transmissionRoughness = roughness;
        transmissionIor = ior;
    }
    else
    {
        transmissionRoughness = AiShaderEvalParamFlt(p_transmissionRoughness);
        transmissionRoughness *= transmissionRoughness;
        transmissionIor = std::max(1.001f, AiShaderEvalParamFlt(p_transmissionIor));
    }
    AtRGB transmissionColor = AiShaderEvalParamRGB(p_transmissionColor) * AiShaderEvalParamFlt(p_transmissionStrength);

    AtRGB ssScattering = AiShaderEvalParamRGB(p_ssScattering);
    AtRGB ssAbsorption = AiShaderEvalParamRGB(p_ssAbsorption);
    float ssDensityScale = AiShaderEvalParamFlt( p_ssDensityScale );
    float ssStrength = AiShaderEvalParamFlt( p_ssStrength );
    float ssDirection = AiShaderEvalParamFlt(p_ssDirection);
    float ssBalance = AiShaderEvalParamFlt(p_ssBalance);
    AtRGB ssTargetColor = AiShaderEvalParamRGB(p_ssTargetColor);
    bool ssSpecifyCoefficients = AiShaderEvalParamBool(p_ssSpecifyCoefficients);
    bool ssInScattering = AiShaderEvalParamBool(p_ssInScattering);

    // precalculate scattering coefficients as we'll need them for shadows etc.
    AtRGB sigma_t = AI_RGB_BLACK;
    AtRGB sigma_s = AI_RGB_BLACK;
    AtRGB sigma_a = AI_RGB_BLACK;
    if (ssStrength > IMPORTANCE_EPS)
    {
        if (ssSpecifyCoefficients)
        {
            sigma_s = ssScattering * ssDensityScale;
            sigma_a = ssAbsorption * ssDensityScale;
            sigma_t = sigma_s + sigma_a;
        }
        else
        {
            sigma_s = sigma_a = AI_RGB_WHITE - ssTargetColor;
            sigma_s *= ssBalance * ssDensityScale;
            sigma_a *= (1.0f - ssBalance) * ssDensityScale;
            sigma_t = sigma_s + sigma_a;
        }
    }

    // if it's a shadow ray, handle shadow colouring through absorption
    // algorithm based heavily on the example in Kettle
    if (sg->Rt & AI_RAY_SHADOW)
    {
        // if the object is transmissive and
        AtRGB outOpacity = AI_RGB_WHITE;
        if (maxh(transmissionColor))
        {
            // check transmission through the surface
            float costheta = AiV3Dot(sg->Nf, -sg->Rd);
            float kt = 1.0f - fresnel(costheta, 1.0f/transmissionIor);
            if (kt >= IMPORTANCE_EPS) // else surface is fully reflective
            {
                if (maxh(sigma_t) > 0.0f)
                {
                    AtPoint alsPreviousIntersection;
                    AtRGB als_sigma_t = sigma_t;
                    if (AiStateGetMsgPnt("alsPreviousIntersection", &alsPreviousIntersection))
                    {
                        AiStateGetMsgRGB("alsPrevious_sigma_t", &als_sigma_t);
                        bool doExtinction = false;
                        if (AiV3Dot(sg->N, sg->Rd) < 0.0f)
                        {
                            // ray is entering a closed volume
                            bool alsInside;
                            AiStateGetMsgBool("alsInside", &alsInside);
                            if (alsInside)
                            {
                                // ray is entering an embedded volume
                                doExtinction = true;
                            }
                            else
                            {
                                // shouldn't get here
                            }

                        }
                        else
                        {
                            // ray is exiting a closed volume
                            doExtinction = true;
                        }

                        if (doExtinction)
                        {
                            float z = AiV3Dist(sg->P, alsPreviousIntersection);
                            outOpacity.r = fast_exp(-z * als_sigma_t.r);
                            outOpacity.g = fast_exp(-z * als_sigma_t.g);
                            outOpacity.b = fast_exp(-z * als_sigma_t.b);
                            outOpacity = 1.0f - (outOpacity*kt);
                        }

                    }
                    else
                    {
                        // first intersection
                        // tell the next shader invocation that we're now inside the surface and what our extinction
                        // coefficient is
                        AiStateSetMsgRGB("alsPrevious_sigma_t", sigma_t);
                        AiStateSetMsgBool("alsInside", true);
                    }
                }
                else // no extinction, shadows are fresnel only.
                {
                    AiStateSetMsgRGB("alsPrevious_sigma_t", AI_RGB_BLACK);
                    outOpacity = 1.0f - kt;
                }
            }
        }

        // store intersection position
        AiStateSetMsgPnt("alsPreviousIntersection", sg->P);
        sg->out_opacity = outOpacity;
        return;
    }

    // Evaluate bump;
    AtVector N_orig;
    AtVector Nf_orig;
    AtRGB bump = AiShaderEvalParamRGB(p_bump);

    // Initialize parameter temporaries
    // TODO: reorganize this so we're not evaluating upstream when we don't need the parameters, e.g. in shadow rays
    AtRGB diffuseColor = AiShaderEvalParamRGB( p_diffuseColor ) * AiShaderEvalParamFlt( p_diffuseStrength );
    AtRGB backlightColor = AiShaderEvalParamRGB(p_backlightColor) * AiShaderEvalParamFlt(p_backlightStrength);
    float diffuseRoughness = AiShaderEvalParamFlt(p_diffuseRoughness);
    bool diffuseEnableCaustics = AiShaderEvalParamFlt(p_diffuseEnableCaustics);
    AtRGB emissionColor = AiShaderEvalParamRGB(p_emissionColor) * AiShaderEvalParamFlt(p_emissionStrength);
    float sssMix = AiShaderEvalParamFlt( p_sssMix );
    AtRGB sssRadiusColor = AiShaderEvalParamRGB( p_sssRadiusColor );
    float sssRadius = AiShaderEvalParamFlt( p_sssRadius );
    float sssDensityScale = AiShaderEvalParamFlt( p_sssDensityScale );
    AtRGB specular1Color = AiShaderEvalParamRGB( p_specular1Color ) * AiShaderEvalParamFlt( p_specular1Strength );
    AtRGB specular2Color = AiShaderEvalParamRGB( p_specular2Color ) * AiShaderEvalParamFlt( p_specular2Strength );
    AtVector specular1Normal = sg->Nf;
    if (data->specular1NormalConnected)
    {
        specular1Normal = AiShaderEvalParamVec(p_specular1Normal);
    }

    AtVector specular2Normal = sg->Nf;
    if (data->specular2NormalConnected)
    {
        specular2Normal = AiShaderEvalParamVec(p_specular2Normal);
    }

    float roughness2 = AiShaderEvalParamFlt( p_specular2Roughness );
    roughness2 *= roughness2;

    float eta = 1.0f / ior;
    float ior2 = std::max(1.001f, AiShaderEvalParamFlt(p_specular2Ior));
    float eta2 = 1.0f / ior2;


    float specular1RoughnessDepthScale = AiShaderEvalParamFlt(p_specular1RoughnessDepthScale);
    float specular2RoughnessDepthScale = AiShaderEvalParamFlt(p_specular2RoughnessDepthScale);
    float transmissionRoughnessDepthScale = AiShaderEvalParamFlt(p_transmissionRoughnessDepthScale);


    bool transmissionEnableCaustics = AiShaderEvalParamBool(p_transmissionEnableCaustics);

    // clamp roughnesses
    roughness = std::max(0.0001f, roughness);
    roughness2 = std::max(0.0001f, roughness2);
    transmissionRoughness = std::max(0.0001f, transmissionRoughness);

    // Initialize result temporaries
    AtRGB result_diffuseDirect = AI_RGB_BLACK;
    AtRGB result_diffuseDirectRaw = AI_RGB_BLACK;
    AtRGB result_backlightDirect = AI_RGB_BLACK;
    AtRGB result_glossyDirect = AI_RGB_BLACK;
    AtRGB result_glossy2Direct = AI_RGB_BLACK;
    AtRGB result_diffuseIndirect = AI_RGB_BLACK;
    AtRGB result_diffuseIndirectRaw = AI_RGB_BLACK;
    AtRGB result_backlightIndirect = AI_RGB_BLACK;
    AtRGB result_glossyIndirect = AI_RGB_BLACK;
    AtRGB result_glossy2Indirect = AI_RGB_BLACK;
    AtRGB result_sss = AI_RGB_BLACK;
    AtRGB result_ss = AI_RGB_BLACK;
    AtColor result_transmission = AI_RGB_BLACK;
    AtColor result_emission = AI_RGB_BLACK;
    // Set up flags to early out of calculations based on where we are in the ray tree
    bool do_diffuse = true;
    bool do_backlight = true;
    bool do_glossy = true;
    bool do_glossy2 = true;
    bool do_ss = true;
    bool do_sss = true;
    bool do_transmission = true;
    int glossy_samples = data->GI_glossy_samples;
    int diffuse_samples = data->GI_diffuse_samples;

    if ( sg->Rr_diff > data->GI_diffuse_depth || maxh(diffuseColor) < IMPORTANCE_EPS || sssMix == 1.0f)
    {
        do_diffuse = false;
    }

    if (maxh(backlightColor) < IMPORTANCE_EPS || !do_diffuse)
    {
        do_backlight = false;
    }

    if (sg->Rr_gloss > data->GI_glossy_depth
                || (sg->Rr_diff > 0)                                    // disable glossy->diffuse caustics
                || maxh(specular1Color) < IMPORTANCE_EPS                // skip evaluations that aren't important
                || (sg->Rr_refr > 1 && !transmissionEnableCaustics))    // disable glossy->transmitted caustics
    {
        do_glossy = false;
    }

    if (sg->Rr_gloss > data->GI_glossy_depth
            || (sg->Rr_diff > 0)                                    // disable glossy->diffuse caustics
            || maxh(specular2Color) < IMPORTANCE_EPS                // skip evaluations that aren't important
            || (sg->Rr_refr > 1 && !transmissionEnableCaustics))    // disable glossy->transmitted caustics
    {
        do_glossy2 = false;
    }

    // Grab the roughness from the previous surface and make sure we're slightly rougher than it to avoid glossy-glossy fireflies
    float alsPreviousRoughness = 0.05f;
    AiStateGetMsgFlt("alsPreviousRoughness", &alsPreviousRoughness);
    if (sg->Rr > 0)
    {
        roughness = std::max(roughness, alsPreviousRoughness*specular1RoughnessDepthScale);
        roughness2 = std::max(roughness2, alsPreviousRoughness*specular2RoughnessDepthScale);
        transmissionRoughness = std::max(transmissionRoughness, alsPreviousRoughness*transmissionRoughnessDepthScale);
    }

    if (sg->Rr_diff > 0 || sg->Rr_gloss > 1 || sssMix < 0.01f)
    {
        do_sss = false;
        sssMix = 0.0f;
    }

    // make sure diffuse and transmission can't sum > 1
    transmissionColor *= 1.0f - maxh(diffuseColor);

    if ((sg->Rr_diff > 0) || maxh(transmissionColor) < IMPORTANCE_EPS)
    {
        do_transmission = false;
    }

    // build a local frame for sampling
    AtVector U, V;
    AiBuildLocalFramePolar(&U, &V, &sg->N);

    AtVector wo = -sg->Rd;

    // prepare temporaries for light group calculation
    AtRGB lightGroupsDirect[NUM_LIGHT_GROUPS];
    memset(lightGroupsDirect, 0, sizeof(AtRGB)*NUM_LIGHT_GROUPS);

    // if this is a camera ray, prepare the temporary storage for deep groups
    AtRGB* deepGroupPtr = NULL;
    AtRGB result_directGroup[NUM_LIGHT_GROUPS];
    for (int i=0; i < NUM_LIGHT_GROUPS; ++i) result_directGroup[i] = AI_RGB_BLACK;
    bool doDeepGroups = data->lightGroupsIndirect;
    int idx = 0;

#define DEEP_DEBUG_DEPTH 5
#ifdef DEEP_DEBUG    
    AtRGB* deepDebugPtr = NULL;
#endif
    if (doDeepGroups && sg->Rt & AI_RAY_CAMERA)
    {
        // if this is a camera ray allocate the group storage
        deepGroupPtr = (AtRGB*)AiShaderGlobalsQuickAlloc(sg, sizeof(AtRGB)*NUM_LIGHT_GROUPS*data->total_samples);
        memset(deepGroupPtr, 0, sizeof(AtRGB)*NUM_LIGHT_GROUPS*data->total_samples);
        AiStateSetMsgPtr("als_deepGroupPtr", deepGroupPtr);

#ifdef DEEP_DEBUG
        // for debugging...
        deepDebugPtr = (AtRGB*)AiShaderGlobalsQuickAlloc(sg, sizeof(AtRGB)*DEEP_DEBUG_DEPTH * 5);
        memset(deepDebugPtr, 0, sizeof(AtRGB)*DEEP_DEBUG_DEPTH * 5);
        AiStateSetMsgPtr("als_deepDebugPtr", deepDebugPtr);
#endif
    }
    else if (doDeepGroups)
    {
        // secondary ray hit - get the pointer from the state
        // if the pointer hasn't been set we're being called from a BSDF that doesn't have deep group support
        // so don't try and do it or we'll be in (crashy) trouble!
        if (!AiStateGetMsgPtr("als_deepGroupPtr", (void**)&deepGroupPtr)) doDeepGroups = false;
        // Get the current sample index from the state
        // This wil be overriden if we're in a child event
        if (!AiStateGetMsgInt("als_sampleIndex", &idx)) doDeepGroups = false;
#ifdef DEEP_DEBUG
        AiStateGetMsgPtr("als_deepDebugPtr", (void**)&deepDebugPtr);
#endif
    }

    // Accumulator for transmission integrated according to the specular1 brdf. Will be used to attenuate diffuse,
    // glossy2, sss and transmission
    float kti = 1.0f;
    float kti2 = 1.0f;

     // storage for all deepgroup contributions
    AtRGB deepGroupsGlossy[NUM_LIGHT_GROUPS];
    AtRGB deepGroupsGlossy2[NUM_LIGHT_GROUPS];
    AtRGB deepGroupsDiffuse[NUM_LIGHT_GROUPS];
    AtRGB deepGroupsTransmission[NUM_LIGHT_GROUPS];
    memset(deepGroupsGlossy, 0, sizeof(AtRGB)*NUM_LIGHT_GROUPS);
    memset(deepGroupsGlossy2, 0, sizeof(AtRGB)*NUM_LIGHT_GROUPS);
    memset(deepGroupsDiffuse, 0, sizeof(AtRGB)*NUM_LIGHT_GROUPS);
    memset(deepGroupsTransmission, 0, sizeof(AtRGB)*NUM_LIGHT_GROUPS); 
    int count = 0;

    // Begin illumination calculation

    // Create the BRDF data structures for MIS
    AtVector Nold = sg->N;
    AtVector Nfold = sg->Nf;
    sg->N = sg->Nf = specular1Normal;
    void* mis;
    mis = GlossyMISCreateData(sg,&U,&V,roughness,roughness);
    BrdfData_wrap brdfw;
    brdfw.brdf_data = mis;
    brdfw.sg = sg;
    brdfw.eta = eta;
    brdfw.V = wo;
    brdfw.N = specular1Normal;
    brdfw.kr = 0.0f;


    sg->N = sg->Nf = specular2Normal;   
    void* mis2;
    mis2 = GlossyMISCreateData(sg,&U,&V,roughness2,roughness2);
    BrdfData_wrap brdfw2;
    brdfw2.brdf_data = mis2;
    brdfw2.sg = sg;
    brdfw2.eta = eta2;
    brdfw2.V = wo;
    brdfw2.N = specular2Normal;
    brdfw2.kr = 0.0f;

    sg->N = Nold;

    
    
    void* dmis = AiOrenNayarMISCreateData(sg, diffuseRoughness);

    if (do_backlight) sg->fhemi = false;

    flipNormals(sg);
    void* bmis = AiOrenNayarMISCreateData(sg, diffuseRoughness);
    flipNormals(sg);

    // Light loop
    AiLightsPrepare(sg);
    if (doDeepGroups || (sg->Rt & AI_RAY_CAMERA)) 
    {
        AtRGB LdiffuseDirect, LbacklightDirect, LspecularDirect, Lspecular2Direct;
        while(AiLightsGetSample(sg))
        {
            int lightGroup = data->lightGroups[sg->Lp];
            float specular_strength = AiLightGetSpecular(sg->Lp);
            float diffuse_strength = AiLightGetDiffuse(sg->Lp);
            if (do_glossy)
            {
                sg->Nf = specular1Normal;
                LspecularDirect =
                AiEvaluateLightSample(sg,&brdfw,GlossyMISSample_wrap,GlossyMISBRDF_wrap,GlossyMISPDF_wrap)
                    * specular_strength;
                if (lightGroup >= 0 && lightGroup < NUM_LIGHT_GROUPS)
                {
                    lightGroupsDirect[lightGroup] += LspecularDirect * specular1Color;
                }
                result_glossyDirect += LspecularDirect;
                sg->Nf = Nfold;
            }
            if (do_glossy2)
            {
                sg->Nf = specular2Normal;
                AtFloat r = (1.0f - brdfw.kr*maxh(specular1Color));
                Lspecular2Direct =
                AiEvaluateLightSample(sg,&brdfw2,GlossyMISSample_wrap,GlossyMISBRDF_wrap,GlossyMISPDF_wrap)
                                        * r * specular_strength;
                if (lightGroup >= 0 && lightGroup < NUM_LIGHT_GROUPS)
                {
                    lightGroupsDirect[lightGroup] += Lspecular2Direct * specular2Color;
                }
                result_glossy2Direct += Lspecular2Direct;
                sg->Nf = Nfold;
            }
            AtFloat r = (1.0f - brdfw.kr*maxh(specular1Color)) * (1.0f - brdfw2.kr*maxh(specular2Color));
            if (do_diffuse)
            {
                LdiffuseDirect =
                    AiEvaluateLightSample(sg,dmis,AiOrenNayarMISSample,AiOrenNayarMISBRDF, AiOrenNayarMISPDF)
                                        * r * diffuse_strength;
                if (lightGroup >= 0 && lightGroup < NUM_LIGHT_GROUPS)
                {
                    lightGroupsDirect[lightGroup] += LdiffuseDirect * diffuseColor;
                }
                result_diffuseDirect += LdiffuseDirect;
            }
            if (do_backlight)
            {
                flipNormals(sg);
                LbacklightDirect = 
                    AiEvaluateLightSample(sg,bmis,AiOrenNayarMISSample,AiOrenNayarMISBRDF, AiOrenNayarMISPDF)
                                        * r * diffuse_strength;
                if (lightGroup >= 0 && lightGroup < NUM_LIGHT_GROUPS)
                {
                    lightGroupsDirect[lightGroup] += LbacklightDirect * backlightColor;
                }
                result_backlightDirect += LbacklightDirect;
                flipNormals(sg);
            }
        }
    }
    else
    {
        while(AiLightsGetSample(sg))
        {
            float specular_strength = AiLightGetSpecular(sg->Lp);
            float diffuse_strength = AiLightGetDiffuse(sg->Lp);
            if (do_glossy)
            {
                result_glossyDirect +=
                AiEvaluateLightSample(sg,&brdfw,GlossyMISSample_wrap,GlossyMISBRDF_wrap,GlossyMISPDF_wrap)
                    * specular_strength;
            }
            if (do_glossy2)
            {
                result_glossy2Direct +=
                AiEvaluateLightSample(sg,&brdfw2,GlossyMISSample_wrap,GlossyMISBRDF_wrap,GlossyMISPDF_wrap)
                                        * (1.0f - brdfw.kr*maxh(specular1Color)) 
                                        * specular_strength;
            }
            if (do_diffuse)
            {
                result_diffuseDirect +=
                AiEvaluateLightSample(sg,dmis,AiOrenNayarMISSample,AiOrenNayarMISBRDF, AiOrenNayarMISPDF)
                                        * (1.0f - brdfw.kr*maxh(specular1Color))
                                        * (1.0f - brdfw2.kr*maxh(specular2Color))
                                        * diffuse_strength;
            }
            if (do_backlight)
            {
                flipNormals(sg);
                result_backlightDirect +=
                AiEvaluateLightSample(sg,bmis,AiOrenNayarMISSample,AiOrenNayarMISBRDF, AiOrenNayarMISPDF)
                    * (1.0f - brdfw.kr*maxh(specular1Color))
                    * (1.0f - brdfw2.kr*maxh(specular2Color))
                    * diffuse_strength;
                flipNormals(sg);
            }
        }
    }

    sg->fhemi = true;

    // Multiply by the colors
    result_diffuseDirectRaw = result_diffuseDirect;
    result_diffuseDirect *= diffuseColor;
    result_backlightDirect *= backlightColor;
    result_glossyDirect *= specular1Color;
    result_glossy2Direct *= specular2Color;       

    // Sample BRDFS
    double samples[2];
    AtRay wi_ray;
    AtVector wi;
    AtScrSample scrs;
    AtVector H;
    float kr=1, kt=1;
    AtRGB dgL[NUM_LIGHT_GROUPS];

    if ((sg->Rr_gloss < data->GI_glossy_depth) && do_glossy)
    {
        AtSamplerIterator* sampit = AiSamplerIterator(data->glossy_sampler, sg);
        AiMakeRay(&wi_ray, AI_RAY_GLOSSY, &sg->P, NULL, AI_BIG, sg);
        kti = 0.0f;
        AiStateSetMsgFlt("alsPreviousRoughness", roughness);
        sg->Nf = specular1Normal;
        while(AiSamplerGetSample(sampit, samples))
        {
            wi = GlossyMISSample(mis, samples[0], samples[1]);
            if (AiV3Dot(wi,specular1Normal) > 0.0f)
            {
                // if we're in a camera ray, pass the sample index down to the child SG
                if (doDeepGroups && (sg->Rt & AI_RAY_CAMERA))
                {
                    idx = count;
                    AiStateSetMsgInt("als_sampleIndex", idx);
                }

                // get half-angle vector for fresnel
                wi_ray.dir = wi;
                AiV3Normalize(H, wi+brdfw.V);
                kr = fresnel(std::max(0.0f,AiV3Dot(H,wi)),eta);
                kti += kr;
                if (kr > IMPORTANCE_EPS) // only trace a ray if it's going to matter
                {
                    AiTrace(&wi_ray, &scrs);
                    AtRGB f = GlossyMISBRDF(mis, &wi) / GlossyMISPDF(mis, &wi) * kr;
                    result_glossyIndirect += scrs.color * f;

                    // accumulate the lightgroup contributions calculated by the child shader
                    if (doDeepGroups)
                    {
                        for (int i=0; i < NUM_LIGHT_GROUPS; ++i)
                        {
                            deepGroupsGlossy[i] += deepGroupPtr[i*data->total_samples+idx] * f;
                        }
                    }
                }
            }
            sg->Nf = Nfold;
            count++;
        }
        result_glossyIndirect *= AiSamplerGetSampleInvCount(sampit);
        kti *= AiSamplerGetSampleInvCount(sampit);
        kti = 1.0f - kti*maxh(specular1Color);
        result_glossyIndirect *= specular1Color;

        if (doDeepGroups)
        {
            for (int i=0; i < NUM_LIGHT_GROUPS; ++i)
            {
                deepGroupsGlossy[i] *= AiSamplerGetSampleInvCount(sampit) * specular1Color;
            }
        }

    } // if (do_glossy)

    if ((sg->Rr_gloss < data->GI_glossy_depth) && do_glossy2)
    {
        AtSamplerIterator* sampit = AiSamplerIterator(data->glossy2_sampler, sg);
        AiMakeRay(&wi_ray, AI_RAY_GLOSSY, &sg->P, NULL, AI_BIG, sg);
        kti2 = 0.0f;
        AiStateSetMsgFlt("alsPreviousRoughness", roughness2);
        sg->Nf = specular2Normal;
        while(AiSamplerGetSample(sampit, samples))
        {
            wi = GlossyMISSample(mis2, samples[0], samples[1]);
            if (AiV3Dot(wi,specular2Normal) > 0.0f)
            {
                // if we're in a camera ray, pass the sample index down to the child SG
                if (doDeepGroups && (sg->Rt & AI_RAY_CAMERA))
                {
                    idx = count;
                    AiStateSetMsgInt("als_sampleIndex", idx);
                }

                wi_ray.dir = wi;
                AiV3Normalize(H, wi+brdfw2.V);
                // add the fresnel for this layer
                kr = fresnel(std::max(0.0f,AiV3Dot(H,wi)),eta2);
                if (kr > IMPORTANCE_EPS) // only trace a ray if it's going to matter
                {
                    AiTrace(&wi_ray, &scrs);
                    AtRGB f = GlossyMISBRDF(mis2, &wi) / GlossyMISPDF(mis2, &wi) * kr * kti;
                    result_glossy2Indirect += scrs.color*f;
                    kti2 += kr; 
                    
                    // accumulate the lightgroup contributions calculated by the child shader
                    if (doDeepGroups)
                    {
                        for (int i=0; i < NUM_LIGHT_GROUPS; ++i)
                        {
                            deepGroupsGlossy2[i] += deepGroupPtr[i*data->total_samples+idx] * f;
                        }
                    }
                }
            }
            sg->Nf = Nfold;
            count++;
        }
        result_glossy2Indirect*= AiSamplerGetSampleInvCount(sampit);
        kti2 *= AiSamplerGetSampleInvCount(sampit);
        kti2 = 1.0f - kti2*maxh(specular2Color);
        result_glossy2Indirect *= specular2Color;

        if (doDeepGroups)
        {
            for (int i=0; i < NUM_LIGHT_GROUPS; ++i)
            {
                deepGroupsGlossy2[i] *= AiSamplerGetSampleInvCount(sampit) * specular2Color;
            }
        }
    } // if (do_glossy2)

    if ((sg->Rr_diff < data->GI_diffuse_depth) && kti*kti2*maxh(diffuseColor) > IMPORTANCE_EPS)
    {
        float kr = kti*kti2;
        AtSamplerIterator* sampit = AiSamplerIterator(data->diffuse_sampler, sg);
        AiMakeRay(&wi_ray, AI_RAY_DIFFUSE, &sg->P, NULL, AI_BIG, sg);
        while (AiSamplerGetSample(sampit, samples))
        {
            // if we're in a camera ray, pass the sample index down to the child SG
            if (doDeepGroups && (sg->Rt & AI_RAY_CAMERA))
            {
                idx = count;
                AiStateSetMsgInt("als_sampleIndex", idx);
            }

            // cosine hemisphere sampling as O-N sampling does not work outside of a light loop
            float stheta = sqrtf(samples[0]);
            float phi = AI_PITIMES2 * samples[1];
            wi.x = stheta * cosf(phi);
            wi.y = stheta * sinf(phi);
            wi.z = sqrtf(1.0f - samples[0]);
            AiV3RotateToFrame(wi, U, V, sg->Nf);

            float cos_theta = AiV3Dot(wi, sg->Nf);
            if (cos_theta <= 0.0f) continue;

            float p = cos_theta * AI_ONEOVERPI;
            
            // trace the ray
            wi_ray.dir = wi;
            AiTrace(&wi_ray, &scrs);
            AtRGB f = kr * AiOrenNayarMISBRDF(dmis, &wi) / p;
            result_diffuseIndirectRaw += scrs.color * f;

            // accumulate the lightgroup contributions calculated by the child shader
            if (doDeepGroups)
            {
                for (int i=0; i < NUM_LIGHT_GROUPS; ++i)
                {
                    deepGroupsDiffuse[i] += deepGroupPtr[i*data->total_samples+idx] * f;
                }
            }

            count++;
        }
        result_diffuseIndirectRaw *= AiSamplerGetSampleInvCount(sampit);
        result_diffuseIndirect = result_diffuseIndirectRaw * diffuseColor;

        if (doDeepGroups)
        {
            for (int i=0; i < NUM_LIGHT_GROUPS; ++i)
            {
                deepGroupsDiffuse[i] *= AiSamplerGetSampleInvCount(sampit) * diffuseColor;
            }
        }
    } // if (do_diffuse)

    // TODO: get backlight working with deep groups
    if (do_backlight)
    {
        flipNormals(sg);
        result_backlightIndirect = AiIndirectDiffuse(&sg->Nf, sg) * backlightColor * kti * kti2;
        flipNormals(sg);
    }

    // Emission
    result_emission = emissionColor;

    // Diffusion multiple scattering
    if (do_sss)
    {
        result_sss = AiSSSPointCloudLookupCubic(sg, sssRadius*sssRadiusColor*sssDensityScale) * diffuseColor * kti * kti2;
    }

    // blend sss and diffuse
    result_diffuseDirect *= (1-sssMix);
    result_diffuseIndirect *= (1-sssMix);
    result_backlightDirect *= (1-sssMix);
    result_backlightIndirect *= (1-sssMix);
    result_sss *= sssMix;

    // Refraction
    if (do_transmission)
    {
        double samples[2];
        float n1, n2;
        float kt;
        AtRay wi_ray;
        AtScrSample sample;
        AtVector wi, R;
        bool inside;
        AtRGB sigma_t = sigma_s + sigma_a;
        AtRGB sigma_s_prime = sigma_s*(1.0f-ssDirection);
        AtRGB sigma_t_prime = (sigma_s_prime + sigma_a);
        AtRGB mfp = AI_RGB_WHITE / sigma_t_prime;

        AtSamplerIterator* sampit = AiSamplerIterator(data->refraction_sampler, sg);

        AiMakeRay(&wi_ray, AI_RAY_REFRACTED, &sg->P, NULL, AI_BIG, sg);

        while (AiSamplerGetSample(sampit, samples))
        {
            // if we're in a camera ray, pass the sample index down to the child SG
            if (doDeepGroups && (sg->Rt & AI_RAY_CAMERA))
            {
                idx = count;
                AiStateSetMsgInt("als_sampleIndex", idx);
            }

            // generate a microfacet normal, m
            // eq. 35,36
            float alpha2 = transmissionRoughness*transmissionRoughness;
            float tanThetaM = sqrtf(-alpha2 * logf(1.0f - float(samples[0])));
            float cosThetaM = 1.0f / sqrtf(1.0f + tanThetaM * tanThetaM);
            float sinThetaM = cosThetaM * tanThetaM;
            float phiM = 2.0f * float(AI_PI) * float(samples[1]);
            AtVector m = (cosf(phiM) * sinThetaM) * U +
                         (sinf(phiM) * sinThetaM) * V +
                                       cosThetaM  * sg->N;


            // get the refracted direction given m
            kt = 1.0f - fresnel(transmissionIor, m, wo, R, wi, inside);

            if (kt > IMPORTANCE_EPS) // if not TIR
            {
                // eq. 33
                float cosThetaM2 = cosThetaM * cosThetaM;
                float tanThetaM2 = tanThetaM * tanThetaM;
                float cosThetaM4 = cosThetaM2 * cosThetaM2;
                float D = fast_exp(-tanThetaM2 / alpha2) / (float(AI_PI) * alpha2 *  cosThetaM4);
                // eq. 24
                float pm = D * cosThetaM;
                // eval BRDF*cosNI
                float cosNI = AiV3Dot(sg->N, wi); // N.wi
                float cosNO = AiV3Dot(sg->N, wo);
                // eq. 26, 27: now calculate G1(i,m) and G1(o,m)
                float ao = 1 / (roughness * sqrtf((1.0f - cosNO * cosNO) / (cosNO * cosNO)));
                float ai = 1 / (roughness * sqrtf((1.0f - cosNI * cosNI) / (cosNI * cosNI)));
                float G1o = ao < 1.6f ? (3.535f * ao + 2.181f * ao * ao) / (1 + 2.276f * ao + 2.577f * ao * ao) : 1.0f;
                float G1i = ai < 1.6f ? (3.535f * ai + 2.181f * ai * ai) / (1 + 2.276f * ai + 2.577f * ai * ai) : 1.0f;
                float G = G1o * G1i;
                // eq. 21
                float cosHI = AiV3Dot(m, wi); // m.wi
                float cosHO = AiV3Dot(m, wo); // m.wo
                float Ht2 = transmissionIor * cosHI + cosHO;
                Ht2 *= Ht2;
                float brdf = (fabsf(cosHI * cosHO) * (transmissionIor * transmissionIor) * (G * D)) / fabsf(cosNO * Ht2);
                // eq. 38 and eq. 17
                float pdf = pm * (transmissionIor * transmissionIor) * fabsf(cosHI) / Ht2;

                wi_ray.dir = wi;
                AiTrace(&wi_ray, &sample);
                AtRGB transmittance = AI_RGB_WHITE;
                if (maxh(sigma_t) > 0.0f && !inside)
                {
                    transmittance.r = fast_exp(-sample.z * sigma_t.r);
                    transmittance.g = fast_exp(-sample.z * sigma_t.g);
                    transmittance.b = fast_exp(-sample.z * sigma_t.b);
                }
                AtRGB f = brdf/pdf * transmittance;
                result_transmission += sample.color * f;
                // accumulate the lightgroup contributions calculated by the child shader
                if (doDeepGroups)
                {
                    for (int i=0; i < NUM_LIGHT_GROUPS; ++i)
                    {
                        deepGroupsTransmission[i] += deepGroupPtr[i*data->total_samples+idx] * f;
                    }
                }

                // single scattering
                if (ssStrength > IMPORTANCE_EPS && maxh(sigma_s_prime) > 0.0f && !inside && ssInScattering)
                {
                    AtVector N = sg->N;
                    sg->N = m;
                    result_ss += AiSSSTraceSingleScatter(sg,bssrdfbrdf(sigma_s_prime/sigma_t_prime),mfp,ssDirection,transmissionIor) * ssStrength;
                    sg->N = N;
                }
            }
            else if (AiV3IsZero(wi)) // total internal reflection
            {
                wi_ray.dir = R;
                AiTrace(&wi_ray, &sample);
                AtRGB transmittance = AI_RGB_WHITE;
                if (maxh(sigma_t) > 0.0f && !inside)
                {
                    transmittance.r = fast_exp(-sample.z * sigma_t.r);
                    transmittance.g = fast_exp(-sample.z * sigma_t.g);
                    transmittance.b = fast_exp(-sample.z * sigma_t.b);
                }
                result_transmission += sample.color * transmittance;
                // accumulate the lightgroup contributions calculated by the child shader
                if (doDeepGroups)
                {
                    for (int i=0; i < NUM_LIGHT_GROUPS; ++i)
                    {
                        deepGroupsTransmission[i] += deepGroupPtr[i*data->total_samples+idx] * transmittance;
                    }
                }
            }

            count++;
        }

        result_transmission *= AiSamplerGetSampleInvCount(sampit) * transmissionColor * kti * kti2;
        result_ss *= AiSamplerGetSampleInvCount(sampit) * transmissionColor * kti * kti2;

        if (doDeepGroups)
        {
            for (int i=0; i < NUM_LIGHT_GROUPS; ++i)
            {
                deepGroupsTransmission[i] *= AiSamplerGetSampleInvCount(sampit) * transmissionColor * kti * kti2;
            }
        }
    }

    
#ifdef DEEP_DEBUG
    deepDebugPtr[sg->Rr * 5 + 0] = deepGroupsDiffuse[0]; 
    deepDebugPtr[sg->Rr * 5 + 1] =  deepGroupsGlossy[0]; 
    deepDebugPtr[sg->Rr * 5 + 2] = deepGroupsGlossy2[0]; 
    deepDebugPtr[sg->Rr * 5 + 3] = deepGroupsTransmission[0]; 
    deepDebugPtr[sg->Rr * 5 + 4] = lightGroupsDirect[0];
#endif
    // Now accumulate the deep group brdf results onto the relevant samples
    if (sg->Rt & AI_RAY_CAMERA)
    {
        if (doDeepGroups)
        {
            AtRGB deepGroups[NUM_LIGHT_GROUPS];
            memset(deepGroups, 0, sizeof(AtRGB)*NUM_LIGHT_GROUPS);
            for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
            {
                deepGroups[i] = deepGroupsDiffuse[i] + deepGroupsGlossy[i] + deepGroupsGlossy2[i] + deepGroupsTransmission[i] + lightGroupsDirect[i];
                AiAOVSetRGB(sg, lightGroupNames[i], deepGroups[i]);
            }
        }
        else
        {
            for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
            {
                AiAOVSetRGB(sg, lightGroupNames[i], lightGroupsDirect[i]);
            }
        }
    }
    else if (doDeepGroups)
    {
        int idx;
        AiStateGetMsgInt("als_sampleIndex", &idx);
        for (int i=0; i < NUM_LIGHT_GROUPS; ++i)
        {
            deepGroupPtr[i*data->total_samples+idx] = deepGroupsDiffuse[i] + deepGroupsGlossy[i] + deepGroupsGlossy2[i] + deepGroupsTransmission[i] 
                                                        + lightGroupsDirect[i];

            
        }
    }

    if (sg->Rt & AI_RAY_CAMERA)
    {
        // write AOVs
        AiAOVSetRGB(sg, "diffuse_color", diffuseColor);
        AiAOVSetRGB(sg, "direct_diffuse", result_diffuseDirect);
        AiAOVSetRGB(sg, "direct_diffuse_raw", result_diffuseDirectRaw);
        AiAOVSetRGB(sg, "direct_backlight", result_backlightDirect);
        AiAOVSetRGB(sg, "sss", result_sss);
        AiAOVSetRGB(sg, "direct_specular", result_glossyDirect);
        AiAOVSetRGB(sg, "direct_specular_2", result_glossy2Direct);
        AiAOVSetRGB(sg, "indirect_diffuse", result_diffuseIndirect);
        AiAOVSetRGB(sg, "indirect_diffuse_raw", result_diffuseIndirectRaw);
        AiAOVSetRGB(sg, "indirect_backlight", result_backlightIndirect);
        AiAOVSetRGB(sg, "indirect_specular", result_glossyIndirect);
        AiAOVSetRGB(sg, "indirect_specular_2", result_glossy2Indirect);
        AiAOVSetRGB(sg, "single_scatter", result_ss);
        AiAOVSetRGB(sg, "refraction", result_transmission);
        AiAOVSetRGB(sg, "emission", result_emission);

        // write IDs
        for (int i=0; i < NUM_ID_AOVS; ++i)
        {
            AtRGB tmp;
            // check if output is enabled first in case we have an expensive network upstream
            if (AiAOVEnabled(idAovNames[i], AI_TYPE_RGB))
            {
                tmp = AiShaderEvalParamRGB(p_id1 + i);
                AiAOVSetRGB(sg, idAovNames[i], tmp);
            }
        }

#ifdef DEEP_DEBUG
        char deepDebugName[32];
        for (int i=0; i < DEEP_DEBUG_DEPTH; ++i)
        {
            sprintf(deepDebugName, "deep_diffuse_%d", i);
            AiAOVSetRGB(sg, deepDebugName, deepDebugPtr[i*5+0]);
            sprintf(deepDebugName, "deep_glossy_%d", i);
            AiAOVSetRGB(sg, deepDebugName, deepDebugPtr[i*5+1]);
            sprintf(deepDebugName, "deep_glossy2_%d", i);
            AiAOVSetRGB(sg, deepDebugName, deepDebugPtr[i*5+2]);
            sprintf(deepDebugName, "deep_transmission_%d", i);
            AiAOVSetRGB(sg, deepDebugName, deepDebugPtr[i*5+3]);
            sprintf(deepDebugName, "deep_direct_%d", i);
            AiAOVSetRGB(sg, deepDebugName, deepDebugPtr[i*5+4]);
        }
#endif

        // write data AOVs
        AtRGB uv = AiColorCreate(sg->u, sg->v, 0.0f);
        AiAOVSetRGB(sg, "uv", uv);
        AtRGB depth = AiColorCreate(sg->Rl, AiV3Dot(sg->Nf, wo), 0.0f);
        AiAOVSetRGB(sg, "depth", depth);

    }

    // Sum final result from temporaries
    //
    sg->out.RGB =    result_diffuseDirect
                    +result_backlightDirect
                    +result_sss
                    +result_glossyDirect
                    +result_glossy2Direct
                    +result_diffuseIndirect
                    +result_backlightIndirect
                    +result_glossyIndirect
                    +result_glossy2Indirect
                    +result_ss
                    +result_transmission
                    +result_emission;
}
