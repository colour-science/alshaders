// Hair shader based on ISHair: Importance Sampling for Hair Scattering by Ou et. al 2012
// http://www.cs.dartmouth.edu/~ouj/site/Research/Entries/2012/6/21_ISHair__Importance_Sampling_for_Hair_Scattering.html

#include <ai.h>
#include "alUtil.h"
#include <vector>
#include <algorithm>

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

#define PIOVER4 0.7853981633974483f
#define ONEOVER4PI 0.07957747154594767
#define FOURPI 12.566370614359172



struct HairBsdf
{
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

    HairBsdf(AtNode* n, AtShaderGlobals* sg, ShaderData* d) :
    node(n), data(d)
    {
        depth = sg->Rr;

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

    /// Parameter evaluation. This should be called after opacity() and before anything else.
    inline void evaluateParameters(AtShaderGlobals* sg)
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

        do_diffuse = do_glossy = true;
        if (sg->Rr_diff > 0)
        {
            do_glossy = false;
        }
    }

    inline void AB(float theta_r, float alpha, float beta, float& A, float& B)
    {
        A = atanf((PIOVER4 + theta_r*0.5f - alpha) / beta);
        B = atanf((-PIOVER4 + theta_r*0.5f - alpha) / beta);
    }

    inline float sampleLong(double u, float theta_r, float alpha, float beta, float A, float B)
    {
        float t = beta * tanf(u*(A-B) + B);
        float theta_h = t + alpha;
        return clamp( -0.4999f * AI_PI, 0.4999f * AI_PI, (2.0f*theta_h - theta_r));
    }

    /// Precalculate invariants that will be used across all lobes. This function must be called before any of the bsdf functions during the direct lighting loop
    /// @param wi The incident direction 
    inline void prepareDirectSample(AtVector wi)
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

    /// Precalculate invariants that will be used across all lobes. This function must be called before any of the bsdf functions during the indirect sampling loop
    /// @param wi The incident direction 
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
        if (phi < -AI_PI) phi += AI_PITIMES2;
        if (phi > AI_PI) phi -= AI_PITIMES2;
        cosphi2 = cosf(phi*0.5f);
        sphericalDirection(theta_i, phi_i, V, W, U, wi);
        invariant = cos_theta_i * inv_cos_theta_d2 * AI_ONEOVER2PI;
    }

    /// Gaussian with offset
    inline float g(float beta, float alpha, float theta_h)
    {
        float n = theta_h-alpha;
        return fast_exp(-(n*n)/(2.0f*beta*beta));
    }

    /// Scattering of the R lobe
    inline float bsdfR()
    {
        float Mr = g(beta_R, alpha_R, theta_h);
        float Nr = cosphi2;
        return Mr * Nr;
    }

    /// Scattering of the TT lobe
    inline float bsdfTT()
    {
        float Mtt = g(beta_TT, alpha_TT, theta_h);
        float Ntt = g(gamma_TT, 0.0f, AI_PI-phi);
        return Mtt * Ntt;
    }

    /// Scatterng of the TRT lobe
    inline float bsdfTRT()
    {
        float Mtrt = g(beta_TRT, alpha_TRT, theta_h);
        float Ntrt = cosphi2;
        return Mtrt * Ntrt;
    }

    /// Scattering of the glint lobes
    inline float bsdfg()
    {
        float Mtrt = g(beta_TRT, alpha_TRT, theta_h);
        float Ng = g(gamma_g, 0.0f, fabsf(phi) - phi_g);
        return Mtrt * Ng;
    }

    /// Sample according to the diffuse bsdf (just uniform spherical sampling for now)
    inline void sampleDiffuse(float u1, float u2, AtVector& wi)
    {
        wi = uniformSampleSphere(u1, u2);
    }

    /// Sample according to the R lobe
    inline void sampleR(float u1, float u2)
    {
        theta_i = sampleLong(u1, theta_r, alpha_R, beta_R, A_R, B_R);
        phi = 2.0f * asinf(clamp(2.0f*u2 - 1.0f, -1.0f, 1.0f));
    }

    /// PDF of the R lobe
    inline float pdfR()
    {
        float t = theta_h-alpha_R;
        float pdf_theta = (1.0f / (2.0f*cos_theta_i*(A_R-B_R))) * (beta_R / (t*t + beta_R*beta_R));
        float pdf_phi = cosphi2*0.25f;
        return pdf_theta * pdf_phi;
    }

    ///  Sample according to the TRT lobe
    inline void sampleTRT(float u1, float u2)
    {
        theta_i = sampleLong(u1, theta_r, alpha_TRT, beta_TRT, A_TRT, B_TRT);
        phi = 2.0f * asinf(clamp(2.0f*u2 - 1.0f, -1.0f, 1.0f));
    }

    /// PDF of the TRT lobe
    inline float pdfTRT()
    {
        float t = theta_h-alpha_TRT;
        float pdf_theta = (1.0f / (2.0f*cos_theta_i*(A_TRT-B_TRT))) * (beta_TRT / (t*t + beta_TRT*beta_TRT));
        float pdf_phi = cosphi2*0.25f;
        return pdf_theta * pdf_phi;
    }

    /// Sample according to the TT lobe
    inline void sampleTT(float u1, float u2)
    {
        theta_i = sampleLong(u1, theta_r, alpha_TT, beta_TT, A_TT, B_TT);
        C_TT = 2.0f * atanf(AI_PI/gamma_TT);
        phi = gamma_TT * tanf(C_TT * (u2-0.5f)) + AI_PI;
    }

    /// PDF of the TT lobe
    inline float pdfTT()
    {
        float t = theta_h-alpha_TT;
        float pdf_theta = (1.0f / (2.0f*cos_theta_i*(A_TT-B_TT))) * (beta_TT / (t*t + beta_TT*beta_TT));
        float p = phi-AI_PI;
        float pdf_phi = (1.0f / C_TT) * (gamma_TT / (p*p + gamma_TT*gamma_TT));
        return pdf_theta * pdf_phi;
    }

    /// Sample according to the glints lobes
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

    /// PDF of the glints lobes
    inline float pdfg()
    {
        float t = theta_h-alpha_TRT;
        float pdf_theta = (1.0f / (2.0f*cos_theta_i*(A_TRT-B_TRT))) * (beta_TRT / (t*t + beta_TRT*beta_TRT));
        float p = fabs(phi) - phi_g;
        float pdf_phi = (1.0f / (2.0f * (Cg-Dg))) * (gamma_g / (p*p + gamma_g*gamma_g));
        return pdf_theta * pdf_phi;
    }

    /// Uniformly sample all the glossy lobes and accumulate the result
    inline void sampleGlossyUniform(double u1, double u2)
    {
        // Choose a lobe to sample based on which quadrant we are in
        if (depth < 1)
        {
            if (u1 < 0.5 && u2 < 0.5)
            {
                u1 = 2.0 * u1;
                u2 = 2.0 * u2;
                sampleR(u1, u2);
            }
            else if (u1 >= 0.5 && u2 < 0.5)
            {
                u1 = 2.0 * (1.0 - u1);
                u2 = 2.0 * u2;
                sampleTT(u1, u2);
            }
            else if (u1 < 0.5 && u2 >= 0.5)
            {
                u1 = 2.0 * u1;
                u2 = 2.0 * (1.0 - u2);
                sampleTRT(u1, u2);
            }
            else
            {
                u1 = 2.0 * (1.0 - u1);
                u2 = 2.0 * (1.0 - u2);
                sampleg(u1, u2);
            }
        }
        else
        {
            if (u1 < 0.5)
            {
                u1 = 2.0 * u1;
                sampleR(u1, u2);
            }
            else if (u1 >= 0.5)
            {
                u1 = 2.0 * (1.0 - u1);
                sampleTT(u1, u2);
            }
        }

        // precalculate some stuff
        prepareIndirectSample(wi_ray.dir);
        AtFloat p = invariant / pdfUniform();
        AtScrSample scrs;

        if (p > IMPORTANCE_EPS*0.1)
        {
            // trace our ray
            AiTrace(&wi_ray, &scrs);

            // calculate result
            if (depth < 1)
            {
                result_R_indirect += scrs.color * bsdfR() * p;
                result_TT_indirect += scrs.color * bsdfTT() * p;
                result_TRT_indirect += scrs.color * bsdfTRT() * p;
                result_TRTg_indirect += scrs.color * bsdfg() * p;
            }
            else
            {
                result_R_indirect += scrs.color * bsdfR() * p;
                result_TT_indirect += scrs.color * bsdfTT() * p;
                result_TRT_indirect += scrs.color * bsdfTRT() * p;
                result_TRTg_indirect += scrs.color * bsdfg() * p;
            }
        }
    }

    /// PDF for uniformly sampling all glossy lobes
    inline float pdfUniform()
    {
        float pdf;
        if (depth < 1) pdf = (pdfR() + pdfTT() + pdfTRT() + pdfg())* 0.25f;
        else pdf = (pdfR() + pdfTT()) * 0.5f;

        return pdf;
    }

    /// MIS diffuse sampling
    static AtVector HairDiffuseSample(const void* brdf_data, float u1, float u2)
    {
        return uniformSampleSphere(u1, u2);
    }

    /// MIS diffuse BSDF
    static AtRGB HairDiffuseBsdf(const void* brdf_data, const AtVector* wi)
    {
        const HairBsdf* hb = reinterpret_cast<const HairBsdf*>(brdf_data);
        float tl = AiV3Dot(*wi, hb->U);
        return rgb(sqrtf(1.0f - tl*tl)* AI_ONEOVER2PI);
    }

    /// MIS diffuse PDF
    static float HairDiffusePdf(const void* brdf_data, const AtVector* wi)
    {
        return ONEOVER4PI;
    }

    /// Integrate the direct illumination for all diffuse and glossy lobes
    inline void integrateDirect(AtShaderGlobals* sg)
    {
        // Tell Arnold we want the full sphere for lighting.
        sg->fhemi = false;
        AiLightsPrepare(sg);
        while (AiLightsGetSample(sg))
        {
            if (do_diffuse)
            {
                result_diffuse_direct += AiEvaluateLightSample(sg, this, HairBsdf::HairDiffuseSample, HairBsdf::HairDiffuseBsdf, HairBsdf::HairDiffusePdf);
            }

            if (do_glossy)
            {
                prepareDirectSample(sg->Ld);
                AtRGB L = sg->Li * sg->we * invariant;
                if (maxh(L) > IMPORTANCE_EPS)
                {
                    result_R_direct += L * bsdfR();
                    result_TT_direct += L * bsdfTT();
                    result_TRT_direct += L * bsdfTRT();
                    result_TRTg_direct += L * bsdfg();
                }
            }
        }
        result_diffuse_direct *= diffuseColor;
        result_R_direct *= specular1Color;
        result_TT_direct *= transmissionColor;
        result_TRT_direct *= specular2Color;
        result_TRTg_direct *= specular2Color * glintStrength;
        sg->fhemi = true;
    }

    /// Integrate the indirect illumination for all diffuse and glossy lobes
    inline void integrateIndirect(AtShaderGlobals* sg)
    {
        if (do_diffuse)
        {
            AtSamplerIterator* sampit = AiSamplerIterator(data->sampler_diffuse, sg);
            AiMakeRay(&wi_ray, AI_RAY_DIFFUSE, &sg->P, NULL, AI_BIG, sg);
            while(AiSamplerGetSample(sampit, samples))
            {
                sampleDiffuse(samples[0], samples[1], wi_ray.dir);
                AiTrace(&wi_ray, &scrs);
                float tl = AiV3Dot(wi_ray.dir, U);
                result_diffuse_indirect += scrs.color * sqrtf(1.0f - tl*tl);
            }
            // The FOURPI factor comes from moving the pdf 1/4*pi pdf out of the loop
            result_diffuse_indirect *= diffuseColor * AiSamplerGetSampleInvCount(sampit) * AI_ONEOVER2PI * FOURPI;
        }

        if (do_glossy)
        {
            AiMakeRay(&wi_ray, AI_RAY_GLOSSY, &sg->P, NULL, AI_BIG, sg);
            
            sampit = AiSamplerIterator(data->sampler_R, sg);
            while(AiSamplerGetSample(sampit, samples))
            {
                sampleGlossyUniform(samples[0], samples[1]);
            }
            float weight = AiSamplerGetSampleInvCount(sampit);
            result_R_indirect *= specular1Color * weight;
            result_TT_indirect *= transmissionColor * weight;
            result_TRT_indirect *= specular2Color * weight;
            result_TRTg_indirect *= specular2Color * glintStrength * weight;
        }
    }

    inline void writeResult(AtShaderGlobals* sg)
    {
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

    /// Calculate opacity of this shading point.
    /// @return true if we should early-out, false otherwise
    /// @param sg Shader globals
    inline bool opacity(AtShaderGlobals* sg)
    {
        AtRGB opacity = AiShaderEvalParamRGB(p_opacity);
        float geo_opacity = 1.0f;
        if (AiUDataGetFlt("geo_opacity", &geo_opacity))
        {
            opacity *= geo_opacity;
        }

        //if (!(sg->Rt & AI_RAY_SHADOW || sg->Rt & AI_RAY_CAMERA))
        if (sg->transp_index > 5)
        {
            opacity = 1.0f;
        }

        // early out if in shadow ray or fully transparent
        return ((sg->Rt & AI_RAY_SHADOW) || AiShaderGlobalsApplyOpacity(sg, opacity));
    }

    AtVector U, V, W;   //< local coordinate frame
    float theta_r;      //< exitant spherical theta
    float phi_r;        //< existant spherical phi
    float cn;           //< random value per curve in [0,1)

    float beta_R;       //< R width
    float alpha_R;      //< R shift
    float beta_TT;      //< TT width
    float alpha_TT;     //< TT shift
    float beta_TRT;     //< TRT width
    float alpha_TRT;    //< TRT shift
    float gamma_TT;     //< TT rolloff
    float gamma_g;      //< g rolloff
    float phi_g;        //< g separation

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

    bool do_diffuse;
    bool do_glossy;

    int depth;

    AtNode* node;
    ShaderData* data;
    AtRay wi_ray;
    AtScrSample scrs;
    double samples[2];
    AtSamplerIterator* sampit;
    
};


node_parameters
{
    AiParameterFlt("ior", 1.55f);
    AiParameterFlt("diffuseStrength", 0.2f);
    AiParameterRGB("diffuseColor", 0.31f, 0.08f, 0.005f);
    AiParameterFlt("specularShift", 7.0f);
    AiParameterFlt("specularWidth", 5.0f);
    AiParameterInt("extraSamples", 0);
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
    HairBsdf::ShaderData* data = new HairBsdf::ShaderData;
    AiNodeSetLocalData(node, data);
}

node_finish
{
    if (AiNodeGetLocalData(node))
    {
        HairBsdf::ShaderData* data = (HairBsdf::ShaderData*)AiNodeGetLocalData(node);
        delete data;
    }
}

node_update
{
    HairBsdf::ShaderData* data = (HairBsdf::ShaderData*)AiNodeGetLocalData(node);
    AtNode *options   = AiUniverseGetOptions();
    int diffuse_samples = AiNodeGetInt(options, "GI_diffuse_samples");
    std::cerr << "extra samples: " << params[p_extraSamples].INT << std::endl;
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





shader_evaluate
{

    // Get shader data
    HairBsdf::ShaderData* data = (HairBsdf::ShaderData*)AiNodeGetLocalData(node);

    // Create HairBsdf object 
    HairBsdf hb(node, sg, data);

    // Do opacity and early-out if possible
    if (hb.opacity(sg)) return;

    // Get parameters
    hb.evaluateParameters(sg);

    // Do direct illumination
    hb.integrateDirect(sg);

    // Do indirect illumination
    hb.integrateIndirect(sg);

    // Writeshader result
    hb.writeResult(sg);
}


