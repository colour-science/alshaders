// Hair shader based on ISHair: Importance Sampling for Hair Scattering by Ou et. al 2012
// http://www.cs.dartmouth.edu/~ouj/site/Research/Entries/2012/6/21_ISHair__Importance_Sampling_for_Hair_Scattering.html

#include <ai.h>
#include "alUtil.h"
#include <vector>
#include <algorithm>

// hard-code IOR for now
// this is completely buggered
#define IOR 1.6f
#define FRESNEL_SAMPLES 1024
float hairFresnel(float phi, float ior)
{
    float roots[3] = {0};
    float pi3 = AI_PI*AI_PI*AI_PI;
    float rPerp = 1.0f;
    float rParal = 1.0f;
    solveCubic(0.0f, 0.0f, -2.0f, -phi, roots);
    float gamma_i = fabsf(roots[0]);
    if (gamma_i > AI_PIOVER2 )
    {
            gamma_i = AI_PI - gamma_i;
    }

    float sin_gamma = sinf(gamma_i);
    float etaPerp = sqrtf(ior*ior - sin_gamma*sin_gamma) / cosf(gamma_i);
    float inv_etaPerp = 1.0f/etaPerp;
    float etaParal = (ior*ior)*inv_etaPerp;
    float inv_etaParal = 1.0f/etaParal;

    // perp
    float a = inv_etaPerp*sinf(gamma_i);
    a *= a;
    if ( a <= 1 )
    {
            float b = etaPerp*sqrtf(1-a);
            float c = cosf(gamma_i);
            rPerp =  ( c - b ) / ( c + b );
            rPerp *= rPerp;
            rPerp = std::min(1.0f, rPerp);
    }

    // parl
    float d = inv_etaParal*sinf(gamma_i);
    d *= d;
    if ( d <= 1 )
    {
            float e = sqrtf(1-d);
            float f = etaParal*cosf(gamma_i);
            rParal = (e - f) / (e + f);
            rParal *= rParal;
            rParal = std::min(1.0f, rParal);
    }

    return 0.5f * (rPerp + rParal);
}


struct ShaderData
{
    ShaderData()
    : sampler_diffuse(NULL), sampler_R(NULL), sampler_TT(NULL), sampler_TRT(NULL), sampler_g(NULL)
    {
        invFresnelSamples = 1.0f / FRESNEL_SAMPLES;
        for (int i=0; i < FRESNEL_SAMPLES; ++i)
        {
            float phi = (AI_PITIMES2 * float(i) * invFresnelSamples) - AI_PI;
            kr[i] = hairFresnel(phi, IOR);
        }
    }
    ~ShaderData()
    {
        AiSamplerDestroy(sampler_diffuse);
        AiSamplerDestroy(sampler_R);
        AiSamplerDestroy(sampler_TT);
        AiSamplerDestroy(sampler_TRT);
        AiSamplerDestroy(sampler_g);
    }

    void update(int diffuse_samples, int glossy_samples)
    {
        AiSamplerDestroy(sampler_diffuse);
        sampler_diffuse = AiSampler(diffuse_samples, 2);
        AiSamplerDestroy(sampler_R);
        sampler_R = AiSampler(glossy_samples, 2);
        AiSamplerDestroy(sampler_TT);
        sampler_TT = AiSampler(glossy_samples, 2);
        AiSamplerDestroy(sampler_TRT);
        sampler_TRT = AiSampler(glossy_samples, 2);
        AiSamplerDestroy(sampler_g);
        sampler_g = AiSampler(glossy_samples, 2);
    }

    float fresnelLookup(float phi)
    {
        int idx = std::min(FRESNEL_SAMPLES-1, static_cast<int>(floorf((phi + AI_PI) * AI_ONEOVER2PI * FRESNEL_SAMPLES)));
        return kr[idx];
    }

    AtSampler* sampler_diffuse;
    AtSampler* sampler_R;
    AtSampler* sampler_TT;
    AtSampler* sampler_TRT;
    AtSampler* sampler_g;

    float kr[FRESNEL_SAMPLES];
    float invFresnelSamples;
};

AI_SHADER_NODE_EXPORT_METHODS(alHair)

enum alHairParams
{
    p_ior,
    p_diffuseStrength,
    p_diffuseColor,
    p_specularShift,
    p_specularWidth,
    p_extraSamples,
    p_specular1Strength,
    p_specular1Color,
    p_specular2Strength,
    p_specular2Color,
    p_glintStrength,
    p_glintRolloff,
    p_glintSeparation,
    p_transmissionStrength,
    p_transmissionColor,
    p_transmissionRolloff,
    p_opacity
};

node_parameters
{
    AiParameterFlt("ior", 1.55f);
    AiParameterFlt("diffuseStrength", 0.2f);
    AiParameterRGB("diffuseColor", 0.31f, 0.08f, 0.005f);
    AiParameterFlt("specularShift", 7.0f);
    AiParameterFlt("specularWidth", 5.0f);
    AiParameterFlt("extraSamples", 0);
    AiParameterFlt("specular1Strength", 1.0f);
    AiParameterRGB("specular1Color", 1.0f, 1.0f, 1.0f);
    AiParameterFlt("specular2Strength", 1.0f);
    AiParameterRGB("specular2Color", 0.31f, 0.08f, 0.005f);
    AiParameterFlt("glintStrength", 2.0f);
    AiParameterFlt("glintRolloff", 5.0f);
    AiParameterFlt("glintSeparation", 35.0f);
    AiParameterFlt("transmissionStrength", 1.0f);
    AiParameterRGB("transmissionColor", 0.92f, 0.7f, 0.64f);
    AiParameterFlt("transmissionRolloff", 30.0f);
    AiParameterRGB("opacity", 1.0f, 1.0f, 1.0f);
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alHair;
   node->output_type = AI_TYPE_RGB;
   node->name        = "alHair";
   node->node_type   = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return TRUE;
}

node_initialize
{
    ShaderData* data = new ShaderData;
    AiNodeSetLocalData(node, data);
}

node_finish
{
    if (AiNodeGetLocalData(node))
    {
        ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);
        delete data;
    }
}

node_update
{
    ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);
    AtNode *options   = AiUniverseGetOptions();
    int diffuse_samples = AiNodeGetInt(options, "GI_diffuse_samples");
    int glossy_samples = std::max(0, AiNodeGetInt(options, "GI_glossy_samples") + params[p_extraSamples].INT);
    data->update(diffuse_samples, glossy_samples);

    AiAOVRegister("specularDirect", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("specularIndirect", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("diffuseDirect", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("diffuseIndirect", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("specular2Direct", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("specular2Indirect", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("transmissionDirect", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("transmissionIndirect", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("glintDirect", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("glintIndirect", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
}

float g(float beta, float alpha, float theta_h)
{
    float n = theta_h-alpha;
    return expf(-(n*n)/(2.0f*beta*beta));
}

#define PIOVER4 0.7853981633974483f

void AB(float theta_r, float alpha, float beta, float& A, float& B)
{
    A = atanf((PIOVER4 + theta_r*0.5f - alpha) / beta);
    B = atanf((-PIOVER4 + theta_r*0.5f - alpha) / beta);
}

float sampleLong(double u, float theta_r, float alpha, float beta, float A, float B)
{
    float t = beta * tanf(u*(A-B) + B);
    float theta_h = t + alpha;
    return clamp( -0.4999f * AI_PI, 0.4999f * AI_PI, (2.0f*theta_h - theta_r));
}

shader_evaluate
{
    // get opacity first
    AtRGB opacity = AiShaderEvalParamRGB(p_opacity);
    float geo_opacity = 1.0f;
    if (AiUDataGetFlt("geo_opacity", &geo_opacity))
    {
        opacity *= geo_opacity;
    }
    // early out if in shadow ray or fully transparent
    if ((sg->Rt & AI_RAY_SHADOW) || AiShaderGlobalsApplyOpacity(sg, opacity))
    {
        return;
    }

    // Initialize result temporaries
    AtRGB result_diffuse_direct = AI_RGB_BLACK;
    AtRGB result_diffuse_indirect = AI_RGB_BLACK;
    AtRGB result_R_direct = AI_RGB_BLACK;
    AtRGB result_R_indirect = AI_RGB_BLACK;
    AtRGB result_TT_direct = AI_RGB_BLACK;
    AtRGB result_TT_indirect = AI_RGB_BLACK;
    AtRGB result_TRT_direct = AI_RGB_BLACK;
    AtRGB result_TRT_indirect = AI_RGB_BLACK;
    AtRGB result_TRTg_direct = AI_RGB_BLACK;
    AtRGB result_TRTg_indirect = AI_RGB_BLACK;

    // Get a local coordinate frame based on the hair fibre direction
    AtVector U = AiV3Normalize(sg->dPdv);
    AtVector V = AiV3Cross(U, sg->N);
    AtVector W = AiV3Cross(U, V);

    // Get the spherical angles of the exitant direction relative to the hair fibre
    AtVector wo = -sg->Rd;
    float theta_r = AI_PIOVER2 - sphericalTheta(wo, U);
    float phi_r = sphericalPhi(wo, V, W);

    // Get a random value per curve
    AtUInt32 curve_id = 0;
    float cn = 1.0f;
    if (AiUDataGetUInt("curve_id", &curve_id))
    {
        AtPoint2 p; p.x = float(curve_id); p.y = 0.0f;
        cn = AiCellNoise2(p);
    }

    // Get shader data
    ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);

    // Get parameters
    float eta = 1.0f / AiShaderEvalParamFlt(p_ior);

    float beta_R = AiShaderEvalParamFlt(p_specularWidth) * AI_DTOR;
    float alpha_R = -AiShaderEvalParamFlt(p_specularShift) * AI_DTOR;

    float beta_TT = beta_R * 0.5f;
    float alpha_TT = alpha_R * 0.5f;

    float beta_TRT = beta_R * 2.0f;
    float alpha_TRT = alpha_R * 1.5f;

    float gamma_TT = AiShaderEvalParamFlt(p_transmissionRolloff) * AI_DTOR;
    float gamma_g = AiShaderEvalParamFlt(p_glintRolloff) * AI_DTOR;
    float phi_g = lerp(30.0f*AI_DTOR, 45.0f*AI_DTOR, cn);

    AtRGB diffuseColor = AiShaderEvalParamRGB(p_diffuseColor) * AiShaderEvalParamFlt(p_diffuseStrength);
    AtRGB specular1Color = AiShaderEvalParamRGB(p_specular1Color) * AiShaderEvalParamFlt(p_specular1Strength);
    AtRGB specular2Color = AiShaderEvalParamRGB(p_specular2Color) * AiShaderEvalParamFlt(p_specular2Strength);
    AtRGB transmissionColor = AiShaderEvalParamRGB(p_transmissionColor) * AiShaderEvalParamFlt(p_transmissionStrength);
    float glintStrength = AiShaderEvalParamFlt(p_glintStrength) * cn;

    bool do_glossy = true;
    bool do_R = true, do_TT = true, do_TRT = true, do_g = true;
    if (sg->Rr > 0)
    {
        do_R = true;
        do_TT = true;
        do_TRT = true;
        do_g = false;
    }
    if (sg->Rr_diff > 0)
    {
        do_glossy = false;
    }

    // Direct lighting loop
    // Tell Arnold we want the full sphere for lighting.
    sg->fhemi = false;
    AiLightsPrepare(sg);
    while (AiLightsGetSample(sg))
    {
        // Diffuse
        float tl = AiV3Dot(sg->Ld, U);
        result_diffuse_direct += sqrtf(1.0f - tl*tl) * sg->Li * sg->we * AI_ONEOVER2PI;

        if (do_glossy)
        {
            // Get angle measures. See Section 3 in Ou et. al.
            float theta_i = AI_PIOVER2 - sphericalTheta(sg->Ld, U);
            float cos_theta_i = fabsf(cosf(theta_i));
            float phi_i = sphericalPhi(sg->Ld, V, W);
            float phi = phi_r - phi_i;
            if (phi < -AI_PI) phi += AI_PITIMES2;
            if (phi > AI_PI) phi -= AI_PITIMES2;
            float theta_h = (theta_r + theta_i)*0.5f;
            float theta_d = (theta_r - theta_i)*0.5f;
            float cos_theta_d = cosf(theta_d);
            float inv_cos_theta_d2 = std::max(0.001f, 1.0f/(cos_theta_d*cos_theta_d));

            // Precalculate invariants across all lobes
            AtRGB L = sg->Li * sg->we * cos_theta_i * inv_cos_theta_d2 * AI_ONEOVER2PI;

            if (1)//maxh(L) > IMPORTANCE_EPS)
            {
                //float kr = data->fresnelLookup(phi);
                //float kt = 1.0f - kr;
                // Calculate longitudinal and azimuthal functions. See Section 3.1 in Ou et. al.
                float Mr = g(beta_R, alpha_R, theta_h);
                float Nr = cosf(phi*0.5f);

                float Mtt = g(beta_TT, alpha_TT, theta_h);
                float Ntt = g(gamma_TT, 0.0f, AI_PI-phi);

                float Mtrt = g(beta_TRT, alpha_TRT, theta_h);
                float Ntrt = cosf(phi*0.5f);

                float Ng = g(gamma_g, 0.0f, fabsf(phi) - phi_g);

                // Sum result temporaries for each lobe
                result_R_direct += L * Mr * Nr;// * kr;
                result_TT_direct += L * Mtt * Ntt;// * kt*kt;
                result_TRT_direct += L * Mtrt * Ntrt;// * kt*kt;
                if (do_g) result_TRTg_direct += L * Mtrt * Ng;//* kt*kt;
                
            }
        }
    }

    // reset this.
    sg->fhemi = true;

    // Multiply by user-defined reflectance
    result_diffuse_direct *= diffuseColor;
    result_R_direct *= specular1Color;
    result_TT_direct *= transmissionColor;
    result_TRT_direct *= specular2Color;
    result_TRTg_direct *= specular2Color * glintStrength;

    // Now sample each lobe
    double samples[2];
    AtRay wi_ray;
    float theta_i, phi, phi_i;
    float pdf;
    AtVector wi;
    AtScrSample scrs;
    // Diffuse
    {
        AtSamplerIterator* sampit = AiSamplerIterator(data->sampler_diffuse, sg);
        AiMakeRay(&wi_ray, AI_RAY_DIFFUSE, &sg->P, NULL, AI_BIG, sg);
        while(AiSamplerGetSample(sampit, samples))
        {
            wi = uniformSampleSphere(samples[0], samples[1]);
            wi_ray.dir = wi;
            AiTrace(&wi_ray, &scrs);
            float tl = AiV3Dot(wi, U);
            result_diffuse_indirect += sqrtf(1.0f - tl*tl) * scrs.color * AI_ONEOVER2PI;
        }
        result_diffuse_indirect *= AiSamplerGetSampleInvCount(sampit);
        result_diffuse_indirect *= diffuseColor;
    }

    if (do_glossy)
    {
        AtSamplerIterator* sampit = AiSamplerIterator(data->sampler_R, sg);
        AiMakeRay(&wi_ray, AI_RAY_GLOSSY, &sg->P, NULL, AI_BIG, sg);
        float theta_i, phi_i;

        // precalculate AB invariants
        float A_R, B_R;
        float A_TRT, B_TRT;
        AB(theta_r, alpha_R, beta_R, A_R, B_R);
        AB(theta_r, alpha_TRT, beta_TRT, A_TRT, B_TRT);

        // sample R
        while (AiSamplerGetSample(sampit, samples))
        {
            float theta_i = sampleLong(samples[0], theta_r, alpha_R, beta_R, A_R, B_R);
            float cos_theta_i = std::max(cosf(theta_i), 0.0001f);
            float theta_h = (theta_r+theta_i)*0.5f;
            float t = theta_h-alpha_R;
            float pdf_long = (1.0f / (2.0f*cos_theta_i*(A_R-B_R))) * (beta_R / (t*t + beta_R*beta_R));
            float phi = 2.0f * asinf(2.0f*samples[1] - 1.0f);
            float phi_i = phi_r - phi;
            if (phi_i < -AI_PI) phi_i += AI_PITIMES2;
            if (phi_i > AI_PI) phi_i -= AI_PITIMES2;
            float cosphi2 = cosf(phi*0.5f);
            float pdf_phi = cosphi2*0.25f;
            float theta_d = (theta_r - theta_i)*0.5f;
            float cos_theta_d = cosf(theta_d);
            float inv_cos_theta_d2 = std::max(0.001f, 1.0f/(cos_theta_d*cos_theta_d));

            float R = cos_theta_i * inv_cos_theta_d2 * g(beta_R, alpha_R, theta_h) * cosphi2 * AI_ONEOVER2PI / (pdf_phi * pdf_long);

            if (R > IMPORTANCE_EPS)
            {
                sphericalDirection(theta_i, phi_i, V, W, U, wi_ray.dir);
                AiTrace(&wi_ray, &scrs);
                //float kr = data->fresnelLookup(phi);
                result_R_indirect += scrs.color * R;
            }
        }
        result_R_indirect *= AiSamplerGetSampleInvCount(sampit) * specular1Color;

        // sample TRT
        sampit = AiSamplerIterator(data->sampler_TRT, sg);
        while (AiSamplerGetSample(sampit, samples))
        {
            float theta_i = sampleLong(samples[0], theta_r, alpha_TRT, beta_TRT, A_TRT, B_TRT);
            float cos_theta_i = cosf(theta_i);
            float theta_h = (theta_r+theta_i)*0.5f;
            float t = theta_h-alpha_TRT;
            float pdf_long = (1.0f / (2.0f*cos_theta_i*(A_TRT-B_TRT))) * (beta_TRT / (t*t + beta_TRT*beta_TRT));
            float phi = 2.0f * asinf(2.0f*samples[1] - 1.0f);
            float phi_i = phi_r - phi;
            if (phi_i < -AI_PI) phi_i += AI_PITIMES2;
            if (phi_i > AI_PI) phi_i -= AI_PITIMES2;
            float cosphi2 = cosf(phi*0.5f);
            float pdf_phi = cosphi2*0.25f;
            float theta_d = (theta_r - theta_i)*0.5f;
            float cos_theta_d = cosf(theta_d);
            float inv_cos_theta_d2 = std::max(0.001f, 1.0f/(cos_theta_d*cos_theta_d));

            float R = cos_theta_i * inv_cos_theta_d2 * g(beta_TRT, alpha_TRT, theta_h) * cosphi2 * AI_ONEOVER2PI / (pdf_phi * pdf_long);

            if (R > IMPORTANCE_EPS)
            {
                sphericalDirection(theta_i, phi_i, V, W, U, wi_ray.dir);
                AiTrace(&wi_ray, &scrs);
                float kt = 1.0f - data->fresnelLookup(phi);
                result_TRT_indirect += scrs.color * R;
            }
        }
        result_TRT_indirect *= AiSamplerGetSampleInvCount(sampit) * specular2Color;
    }

    if (sg->Rt & AI_RAY_CAMERA)
    {
        AiAOVSetRGB(sg, "diffuseDirect", result_diffuse_direct);
        AiAOVSetRGB(sg, "diffuseIndirect", result_diffuse_indirect);
        AiAOVSetRGB(sg, "specularDirect", result_R_direct);
        AiAOVSetRGB(sg, "specularIndirect", result_R_indirect);
        AiAOVSetRGB(sg, "specular2Direct", result_TRT_direct);
        AiAOVSetRGB(sg, "specular2Indirect", result_TRT_indirect);
        AiAOVSetRGB(sg, "glintDirect", result_TRTg_direct);
        AiAOVSetRGB(sg, "glintIndirect", result_TRTg_indirect);
        AiAOVSetRGB(sg, "transmissionDirect", result_TT_direct);
        AiAOVSetRGB(sg, "transmissionIndirect", result_TT_indirect);
    }

    sg->out.RGB =   result_diffuse_direct +
                    result_diffuse_indirect +
                    result_R_direct +
                    result_R_indirect +
                    result_TT_direct +
                    result_TT_indirect +
                    result_TRT_direct +
                    result_TRT_indirect +
                    result_TRTg_direct +
                    result_TRTg_indirect;
}


