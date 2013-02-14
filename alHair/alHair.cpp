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
    return fast_exp(-(n*n)/(2.0f*beta*beta));
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

struct HairBsdf
{
    HairBsdf(AtNode* node, AtShaderGlobals* sg)
    {
        // Get a local coordinate frame based on the hair fibre direction
        U = AiV3Normalize(sg->dPdv);
        V = AiV3Cross(U, sg->N);
        W = AiV3Cross(U, V);

        // Get the spherical angles of the exitant direction relative to the hair fibre
        AtVector wo = -sg->Rd;
        theta_r = AI_PIOVER2 - sphericalTheta(wo, U);
        phi_r = sphericalPhi(wo, V, W);

        // Get a random value per curve
        AtUInt32 curve_id = 0;
        cn = 1.0f;
        if (AiUDataGetUInt("curve_id", &curve_id))
        {
            AtPoint2 p; p.x = float(curve_id); p.y = 0.0f;
            cn = AiCellNoise2(p);
        }

        beta_R = AiShaderEvalParamFlt(p_specularWidth) * AI_DTOR;
        alpha_R = -AiShaderEvalParamFlt(p_specularShift) * AI_DTOR;

        beta_TT = beta_R * 0.5f;
        alpha_TT = alpha_R * 0.5f;

        beta_TRT = beta_R * 2.0f;
        alpha_TRT = alpha_R * 1.5f;

        gamma_TT = AiShaderEvalParamFlt(p_transmissionRolloff) * AI_DTOR;
        gamma_g = AiShaderEvalParamFlt(p_glintRolloff) * AI_DTOR;
        phi_g = lerp(30.0f*AI_DTOR, 45.0f*AI_DTOR, cn);

        AB(theta_r, alpha_R, beta_R, A_R, B_R);
        AB(theta_r, alpha_TT, beta_TT, A_TT, B_TT);
        AB(theta_r, alpha_TRT, beta_TRT, A_TRT, B_TRT);

        diffuseColor = AiShaderEvalParamRGB(p_diffuseColor) * AiShaderEvalParamFlt(p_diffuseStrength);
        specular1Color = AiShaderEvalParamRGB(p_specular1Color) * AiShaderEvalParamFlt(p_specular1Strength);
        specular2Color = AiShaderEvalParamRGB(p_specular2Color) * AiShaderEvalParamFlt(p_specular2Strength);
        transmissionColor = AiShaderEvalParamRGB(p_transmissionColor) * AiShaderEvalParamFlt(p_transmissionStrength);
        glintStrength = AiShaderEvalParamFlt(p_glintStrength) * cn;

        result_diffuse_direct = AI_RGB_BLACK;
        result_diffuse_indirect = AI_RGB_BLACK;
        result_R_direct = AI_RGB_BLACK;
        result_R_indirect = AI_RGB_BLACK;
        result_TT_direct = AI_RGB_BLACK;
        result_TT_indirect = AI_RGB_BLACK;
        result_TRT_direct = AI_RGB_BLACK;
        result_TRT_indirect = AI_RGB_BLACK;
        result_TRTg_direct = AI_RGB_BLACK;
        result_TRTg_indirect = AI_RGB_BLACK;
    }

    void prepareDirectSample(AtVector wi)
    {
        // Get angle measures. See Section 3 in Ou et. al.
        theta_i = AI_PIOVER2 - sphericalTheta(wi, U);
        cos_theta_i = fabsf(cosf(theta_i));
        phi_i = sphericalPhi(wi, V, W);
        phi = phi_r - phi_i;
        if (phi < -AI_PI) phi += AI_PITIMES2;
        if (phi > AI_PI) phi -= AI_PITIMES2;
        cosphi2 = cosf(phi*0.5f);
        theta_h = (theta_r + theta_i)*0.5f;
        theta_d = (theta_r - theta_i)*0.5f;
        cos_theta_d = cosf(theta_d);
        inv_cos_theta_d2 = std::max(0.001f, 1.0f/(cos_theta_d*cos_theta_d));
        invariant = cos_theta_i * inv_cos_theta_d2 * AI_ONEOVER2PI;
    }

    inline float bsdfR()
    {
        float Mr = g(beta_R, alpha_R, theta_h);
        float Nr = cosphi2;
        return Mr * Nr;
    }

    inline float bsdfTT()
    {
        float Mtt = g(beta_TT, alpha_TT, theta_h);
        float Ntt = g(gamma_TT, 0.0f, AI_PI-phi);
        return Mtt * Ntt;
    }

    inline float bsdfTRT()
    {
        float Mtrt = g(beta_TRT, alpha_TRT, theta_h);
        float Ntrt = cosphi2;
        return Mtrt * Ntrt;
    }

    inline float bsdfg()
    {
        float Mtrt = g(beta_TRT, alpha_TRT, theta_h);
        float Ng = g(gamma_g, 0.0f, fabsf(phi) - phi_g);
        return Mtrt * Ng;
    }

    inline void FglossyDirect(AtRGB Li, float weight)
    {
        AtRGB L = Li * weight * invariant;
        if (maxh(Li*invariant) > IMPORTANCE_EPS)
        {
            result_R_direct += L * bsdfR();
            result_TT_direct += L * bsdfTT();
            result_TRT_direct += L * bsdfTRT();
            result_TRTg_direct += L * bsdfg();
        }
    }

    inline void FdiffuseDirect(const AtRGB& Li, float weight)
    {
        result_diffuse_direct += cos_theta_i * Li * weight * AI_ONEOVER2PI;
    }

    inline void scaleDirect()
    {
        result_diffuse_direct *= diffuseColor;
        result_R_direct *= specular1Color;
        result_TT_direct *= transmissionColor;
        result_TRT_direct *= specular2Color;
        result_TRTg_direct *= specular2Color * glintStrength;
    }

    inline void scaleR(float weight)
    {
        result_R_indirect *= specular1Color * weight;   
    }

    inline void scaleTT(float weight)
    {
        result_TT_indirect *= transmissionColor * weight;
    }

    inline void scaleTRT(float weight)
    {
         result_TRT_indirect *= specular2Color * weight;
    }

    inline void scaleg(float weight)
    {
        result_TRTg_indirect *= specular2Color * glintStrength * weight;
    }

    inline void scaleDiffuse(float weight)
    {
        result_diffuse_indirect *= diffuseColor * weight;
    }

    inline void sampleDiffuse(float u1, float u2, AtVector& wi)
    {
        wi = uniformSampleSphere(u1, u2);
    }

    inline void FdiffuseIndirect(const AtRGB& Li, const AtVector& wi)
    {
        float tl = AiV3Dot(wi, U);
        result_diffuse_indirect += Li * sqrtf(1.0f - tl*tl) * AI_ONEOVER2PI;
    }

    inline void prepareIndirectSample(AtVector& wi)
    {
        cos_theta_i = std::max(cosf(theta_i), 0.0001f);
        theta_h = (theta_r+theta_i)*0.5f;
        theta_d = (theta_r - theta_i)*0.5f;
        cos_theta_d = cosf(theta_d);
        inv_cos_theta_d2 = std::max(0.001f, 1.0f/(cos_theta_d*cos_theta_d));
        phi_i = phi_r - phi;
        if (phi_i < -AI_PI) phi_i += AI_PITIMES2;
        if (phi_i > AI_PI) phi_i -= AI_PITIMES2;
        sphericalDirection(theta_i, phi_i, V, W, U, wi);
        invariant = cos_theta_i * inv_cos_theta_d2 * AI_ONEOVER2PI;
    }

    inline bool importantR()
    {
        return maxh(invariant*specular1Color) > IMPORTANCE_EPS;
    }

    inline bool importantTT()
    {
        return maxh(invariant*transmissionColor) > IMPORTANCE_EPS;
    }

    inline bool importantTRT()
    {
        return maxh(invariant*specular2Color) > IMPORTANCE_EPS;
    }

    inline bool importantg()
    {
        return maxh(invariant*specular2Color*glintStrength) > IMPORTANCE_EPS;
    }

    inline void sampleR(float u1, float u2)
    {
        theta_i = sampleLong(u1, theta_r, alpha_R, beta_R, A_R, B_R);
        phi = 2.0f * asinf(2.0f*u2 - 1.0f);
    }

    inline float pdfR()
    {
        float t = theta_h-alpha_R;
        float pdf_theta = (1.0f / (2.0f*cos_theta_i*(A_R-B_R))) * (beta_R / (t*t + beta_R*beta_R));
        cosphi2 = cosf(phi*0.5f);
        float pdf_phi = cosphi2*0.25f;
        return pdf_theta * pdf_phi;
    }

    inline void FR(const AtRGB& Li)
    {
        result_R_indirect += Li * invariant * bsdfR() / pdfR();
    }

    inline void sampleTRT(float u1, float u2)
    {
        theta_i = sampleLong(u1, theta_r, alpha_TRT, beta_TRT, A_TRT, B_TRT);
        phi = 2.0f * asinf(2.0f*u2 - 1.0f);
    }

    inline float pdfTRT()
    {
        float t = theta_h-alpha_TRT;
        float pdf_theta = (1.0f / (2.0f*cos_theta_i*(A_TRT-B_TRT))) * (beta_TRT / (t*t + beta_TRT*beta_TRT));
        cosphi2 = cosf(phi*0.5f);
        float pdf_phi = cosphi2*0.25f;
        return pdf_theta * pdf_phi;
    }

    inline void FTRT(const AtRGB& Li)
    {
        result_TRT_indirect += Li * invariant * bsdfTRT() / pdfTRT();
    }

    inline void sampleTT(float u1, float u2)
    {
        theta_i = sampleLong(u1, theta_r, alpha_TT, beta_TT, A_TT, B_TT);
        C_TT = 2.0f * atanf(AI_PI/gamma_TT);
        phi = gamma_TT * tanf(C_TT * (u2-0.5f)) + AI_PI;
    }

    inline float pdfTT()
    {
        float t = theta_h-alpha_TT;
        float pdf_theta = (1.0f / (2.0f*cos_theta_i*(A_TT-B_TT))) * (beta_TT / (t*t + beta_TT*beta_TT));
        float p = phi-AI_PI;
        float pdf_phi = (1.0f / C_TT) * (gamma_TT / (p*p + gamma_TT*gamma_TT));
        return pdf_theta * pdf_phi;
    }

    inline void FTT(const AtRGB& Li)
    {
        result_TT_indirect += Li * invariant * bsdfTT() / pdfTT();
    }

    inline void sampleg(float u1, float u2)
    {
        theta_i = sampleLong(u1, theta_r, alpha_TRT, beta_TRT, A_TRT, B_TRT);
        float sign;
        if (u2 < 0.5f)
        {
            sign = 1.0f;
            u2 = 2.0f * u2;
        }
        else
        {
            sign = -1.0f;
            u2 = 2.0f * (1.0f-u2);
        }
        Cg = atanf((AI_PIOVER2 - phi_g)/gamma_g);
        Dg = atanf(-phi_g/gamma_g);
        phi = gamma_g * tanf(u2*(Cg-Dg)+Dg) + phi_g;
        phi *= sign;
    }

    inline float pdfg()
    {
        float t = theta_h-alpha_TRT;
        float pdf_theta = (1.0f / (2.0f*cos_theta_i*(A_TRT-B_TRT))) * (beta_TRT / (t*t + beta_TRT*beta_TRT));
        float p = fabs(phi) - phi_g;
        float pdf_phi = (1.0f / (2.0f * (Cg-Dg))) * (gamma_g / (p*p + gamma_g*gamma_g));
        return pdf_theta * pdf_phi;
    }

    inline void Fg(const AtRGB& Li)
    {
        result_TRTg_indirect += Li * invariant * bsdfg() / pdfg();
    }

    AtVector U, V, W; //< local coordinate frame
    float theta_r; //< exitant spherical theta
    float phi_r; //< existant spherical phi
    float cn; //< random value per curve in [0,1)

    float beta_R; //< R width
    float alpha_R; //< R shift
    float beta_TT; //< TT width
    float alpha_TT; //< TT shift
    float beta_TRT; //< TRT width
    float alpha_TRT; //< TRT shift
    float gamma_TT; //< TT rolloff
    float gamma_g; //< g rolloff
    float phi_g; //< g separation

    AtRGB diffuseColor; 
    AtRGB specular1Color;
    AtRGB specular2Color;
    AtRGB transmissionColor;
    float glintStrength;

    float theta_i;
    float cos_theta_i;
    float phi_i;
    float phi;
    float cosphi2;
    float theta_h;
    float theta_d;
    float cos_theta_d;
    float inv_cos_theta_d2;
    float invariant;

    float A_R;
    float B_R;
    float A_TT;
    float B_TT;
    float C_TT;
    float A_TRT;
    float B_TRT;
    float Cg;
    float Dg;

    AtRGB result_diffuse_direct;
    AtRGB result_diffuse_indirect;
    AtRGB result_R_direct;
    AtRGB result_R_indirect;
    AtRGB result_TT_direct;
    AtRGB result_TT_indirect;
    AtRGB result_TRT_direct;
    AtRGB result_TRT_indirect;
    AtRGB result_TRTg_direct;
    AtRGB result_TRTg_indirect;
};

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
    
    // Get shader data
    ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);

    // Create HairBsdf object 
    HairBsdf hb(node, sg);
    
    bool do_diffuse = true;
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
        if (do_diffuse)
        {
            hb.FdiffuseDirect(sg->Li, sg->we);
        }

        if (do_glossy)
        {
            hb.prepareDirectSample(sg->Ld);
            hb.FglossyDirect(sg->Li, sg->we);
        }
    }

    hb.scaleDirect();

    // reset this.
    sg->fhemi = true;

    AtRay wi_ray;
    AtScrSample scrs;
    double samples[2];
    AtSamplerIterator* sampit;

    if (do_diffuse)
    {
        AtSamplerIterator* sampit = AiSamplerIterator(data->sampler_diffuse, sg);
        AiMakeRay(&wi_ray, AI_RAY_DIFFUSE, &sg->P, NULL, AI_BIG, sg);
        while(AiSamplerGetSample(sampit, samples))
        {
            hb.sampleDiffuse(samples[0], samples[1], wi_ray.dir);
            AiTrace(&wi_ray, &scrs);
            hb.FdiffuseIndirect(scrs.color, wi_ray.dir);
        }
        hb.scaleDiffuse(AiSamplerGetSampleInvCount(sampit));
    }

    if (do_glossy)
    {
        
        AiMakeRay(&wi_ray, AI_RAY_GLOSSY, &sg->P, NULL, AI_BIG, sg);

        sampit = AiSamplerIterator(data->sampler_R, sg);
        while(AiSamplerGetSample(sampit, samples))
        {
            hb.sampleR(samples[0], samples[1]);
            hb.prepareIndirectSample(wi_ray.dir);
            if (hb.importantR())
            {
                AiTrace(&wi_ray, &scrs);
                hb.FR(scrs.color);
            }
        }
        hb.scaleR(AiSamplerGetSampleInvCount(sampit));

        sampit = AiSamplerIterator(data->sampler_TRT, sg);
        while(AiSamplerGetSample(sampit, samples))
        {
            hb.sampleTRT(samples[0], samples[1]);
            hb.prepareIndirectSample(wi_ray.dir);
            if (hb.importantTRT())
            {
                AiTrace(&wi_ray, &scrs);
                hb.FTRT(scrs.color);
            }
        }
        hb.scaleTRT(AiSamplerGetSampleInvCount(sampit));

        sampit = AiSamplerIterator(data->sampler_TT, sg);
        while(AiSamplerGetSample(sampit, samples))
        {
            hb.sampleTT(samples[0], samples[1]);
            hb.prepareIndirectSample(wi_ray.dir);
            if (hb.importantTT())
            {
                AiTrace(&wi_ray, &scrs);
                hb.FTT(scrs.color);
            }
        }
        hb.scaleTT(AiSamplerGetSampleInvCount(sampit));

        sampit = AiSamplerIterator(data->sampler_g, sg);
        while(AiSamplerGetSample(sampit, samples))
        {
            hb.sampleg(samples[0], samples[1]);
            hb.prepareIndirectSample(wi_ray.dir);
            if (hb.importantg())
            {
                AiTrace(&wi_ray, &scrs);
                hb.Fg(scrs.color);
            }
        }
        hb.scaleg(AiSamplerGetSampleInvCount(sampit));
    }

    if (sg->Rt & AI_RAY_CAMERA)
    {
        AiAOVSetRGB(sg, "diffuseDirect", hb.result_diffuse_direct);
        AiAOVSetRGB(sg, "diffuseIndirect", hb.result_diffuse_indirect);
        AiAOVSetRGB(sg, "specularDirect", hb.result_R_direct);
        AiAOVSetRGB(sg, "specularIndirect", hb.result_R_indirect);
        AiAOVSetRGB(sg, "specular2Direct", hb.result_TRT_direct);
        AiAOVSetRGB(sg, "specular2Indirect", hb.result_TRT_indirect);
        AiAOVSetRGB(sg, "glintDirect", hb.result_TRTg_direct);
        AiAOVSetRGB(sg, "glintIndirect", hb.result_TRTg_indirect);
        AiAOVSetRGB(sg, "transmissionDirect", hb.result_TT_direct);
        AiAOVSetRGB(sg, "transmissionIndirect", hb.result_TT_indirect);
    }

    sg->out.RGB =   hb.result_diffuse_direct +
                    hb.result_diffuse_indirect +
                    hb.result_R_direct +
                    hb.result_R_indirect +
                    hb.result_TT_direct +
                    hb.result_TT_indirect +
                    hb.result_TRT_direct +
                    hb.result_TRT_indirect +
                    hb.result_TRTg_direct +
                    hb.result_TRTg_indirect;
}


