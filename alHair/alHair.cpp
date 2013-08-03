// Hair shader based on 
// [1] ISHair: Importance Sampling for Hair Scattering by Ou et al. 2012
// http://www.cs.dartmouth.edu/~ouj/site/Research/Entries/2012/6/21_ISHair__Importance_Sampling_for_Hair_Scattering.html
// [2] Dual Scattering Approximation For Fast Multiple Scattering in Hair by Zinke et al. 2008


#include <ai.h>
#include "alUtil.h"
#include "stats.h"
#include <vector>
#include <algorithm>

#define DEBUG_LUTS
#define DEBUG_FRESNEL
#ifdef DEBUG_LUTS
#include "exr.h"
#endif

AI_SHADER_NODE_EXPORT_METHODS(alHair);

enum alHairParams
{
    p_ior,
    p_hairDensity,
    p_hairColor,
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
    p_opacity,
    p_dualDepth,
    p_densityFront,
    p_densityBack
};

// hard-code IOR for now
// this is completely buggered
#define IOR 1.6f
#define FRESNEL_SAMPLES 1024

float fresnel(float n2, float n1, float theta_i)
{
    theta_i = fabsf(theta_i);
    if (theta_i > AI_PIOVER2) theta_i = AI_PI - theta_i;
    float cos_theta_i = cosf(theta_i);
    float sin_theta_i = sinf(theta_i);
    float r = (n1/n2)*sin_theta_i;
    r *= r;
    float Rs = 1.0f;
    float Rp = 1.0f;

    if (r <= 1.0f)
    {
        r = sqrtf(1.0f-r);
        float n1_cos_theta_i = n1*cos_theta_i;
        float n2_cos_theta_i = n2*cos_theta_i;

        Rs = (n1_cos_theta_i - n2*r) / (n1_cos_theta_i + n2*r);
        Rs *= Rs;

        Rp = (n1*r - n2_cos_theta_i) / (n1*r + n2_cos_theta_i);
        Rp *= Rp;
    }

    float f = std::min(1.0f, (Rs+Rp)*0.5f);
    return f;
}

void hairAttenuation(float ior, float theta_d, float phi, AtRGB absorption, AtRGB kfr[3])
{
#define ONEOVERPI3 0.032251534433199495
    // Get miller-bravais indices n' and n''
    float sin_theta_d = sinf(theta_d);
    float A = sqrtf(ior*ior - sin_theta_d*sin_theta_d);
    float n_p = A / cosf(theta_d);
    float n_pp = ior*ior / n_p;

    
    float c = asinf(1.0f/n_p);
    for (int p=0; p < 3; ++p)
    {
        // adjust phi to correct range for this component
        float phi_p = phi;
        /*
        if (p != 1)
        {
            if (phi_p > AI_PI) phi_p -= AI_PITIMES2;
            phi_p += p*AI_PI;
        }
        */
        // get roots of polynomial
        float roots[3] = {0,0,0};
        int numRoots = solveCubic(-8.0f*p*c*ONEOVERPI3, 0.0f, 6.0f*p*c*AI_ONEOVERPI - 2.0f, p*AI_PI - phi_p, roots);
        AtRGB Fr = AI_RGB_BLACK;
        if (p < 2)
        {
            float gamma_i = roots[0];
            if (p==0)
            {
                Fr = rgb(fresnel(n_p, n_pp, gamma_i));
            }
            else
            {
                float gamma_t = asinf(sinf(gamma_i)/n_p);
                float l = 2.0f * cosf(gamma_t);
                Fr = (1.0f - fresnel(n_p, n_pp, gamma_i)) * (1.0f - fresnel(1.0f/n_p, 1.0f/n_pp, gamma_t)) * exp(-absorption*l);
            }
        }
        else 
        {
            for (int i=0; i < numRoots; ++i)
            {
                float gamma_i = roots[i];
                float gamma_t = asinf(sinf(gamma_i)/n_p);
                float l = 2.0f * cosf(gamma_t);

                float Fr_TT = (1.0f - fresnel(n_p, n_pp, gamma_i)) * (1.0f - fresnel(1.0f/n_p, 1.0f/n_pp, gamma_t));
                Fr += Fr_TT * fresnel(1.0f/n_p, 1.0f/n_pp, gamma_t) * exp(-absorption*2*l);
            }
        }
        kfr[p] =  clamp(Fr, AI_RGB_BLACK, AI_RGB_WHITE);
    }
}

void hairAttenuation(float ior, float cos_theta_i, float theta_d, float phi, float phi_h, float aa, AtRGB absorption, AtRGB kfr[3])
{
    hairAttenuation(ior, theta_d, phi, absorption, kfr);
}

#define PIOVER4 0.7853981633974483f
#define ONEOVER4PI 0.07957747154594767
#define FOURPI 12.566370614359172

#define DS_NUMSTEPS 128

#define NORMALIZED_GAUSSIAN
#ifdef NORMALIZED_GAUSSIAN

/// Normalized gaussian with offset
inline float g(float beta, float alpha, float theta_h)
{
    float n = theta_h-alpha;
    return fast_exp(-(n*n)/(2.0f*beta))/sqrtf(2.0f*AI_PI*beta);
}

#else

// Normalized gaussian with offset
inline float g(float beta, float alpha, float theta_h)
{
    float n = theta_h-alpha;
    return fast_exp(-(n*n)/(2.0f*beta));
}

#endif

inline AtRGB g(AtRGB beta, AtRGB alpha, AtRGB theta_h)
{
    return AiColorCreate(
        g(beta.r, alpha.r, theta_h.r),
        g(beta.g, alpha.g, theta_h.g),
        g(beta.b, alpha.b, theta_h.b)
    );
}


/// Scattering of the R lobe
inline float bsdfR(float beta_R, float alpha_R, float theta_h, float cosphi2)
{
    float Mr = g(beta_R, alpha_R, theta_h);
    float Nr = cosphi2;
    return Mr * Nr;
}

/// Scattering of the TT lobe
inline float bsdfTT(float beta_TT, float alpha_TT, float theta_h, float gamma_TT, float phi)
{
    float Mtt = g(beta_TT, alpha_TT, theta_h);
    float Ntt = g(gamma_TT, 0.0f, AI_PI-phi);
    return Mtt * Ntt;
}

/// Scatterng of the TRT lobe
inline float bsdfTRT(float beta_TRT, float alpha_TRT, float theta_h, float cosphi2)
{
    float Mtrt = g(beta_TRT, alpha_TRT, theta_h);
    float Ntrt = cosphi2;
    return Mtrt * Ntrt;
}

/// Scattering of the glint lobes
inline float bsdfg(float beta_TRT, float alpha_TRT, float theta_h, float gamma_g, float phi, float phi_g)
{
    float Mtrt = g(beta_TRT, alpha_TRT, theta_h);
    float Ng = g(gamma_g, 0.0f, phi - phi_g);
    return Mtrt * Ng;
}

struct HairBsdf
{
    struct ShaderData
    {
        ShaderData()
        : sampler_R(NULL), sampler_TT(NULL), sampler_TRT(NULL), sampler_g(NULL)
        {}

        ~ShaderData()
        {
            AiSamplerDestroy(sampler_R);
            AiSamplerDestroy(sampler_TT);
            AiSamplerDestroy(sampler_TRT);
            AiSamplerDestroy(sampler_g);
        }

        void update(AtParamValue* params)
        {
            AtUInt32 t0 = AiMsgUtilGetElapsedTime();

            AtNode *options   = AiUniverseGetOptions();
            int glossy_samples = std::max(0, AiNodeGetInt(options, "GI_glossy_samples") + params[p_extraSamples].INT);

            ior = params[p_ior].FLT;

            AiSamplerDestroy(sampler_R);
            sampler_R = AiSampler(glossy_samples, 2);
            AiSamplerDestroy(sampler_TT);
            sampler_TT = AiSampler(glossy_samples, 2);
            AiSamplerDestroy(sampler_TRT);
            sampler_TRT = AiSampler(glossy_samples, 2);
            AiSamplerDestroy(sampler_g);
            sampler_g = AiSampler(glossy_samples, 2);

            alpha_R = -params[p_specularShift].FLT * AI_DTOR;
            beta_R = params[p_specularWidth].FLT * AI_DTOR;
            
            beta_TT = beta_R * 0.5f;
            
            alpha_TT = alpha_R * 0.5f;

            beta_TRT = beta_R * 2.0f;
            alpha_TRT = alpha_R * 1.5f;

            beta_R2 = beta_R*beta_R;
            beta_TT2 = beta_TT*beta_TT;
            beta_TRT2 = beta_TRT*beta_TRT;

            gamma_TT = params[p_transmissionRolloff].FLT * AI_DTOR;
            gamma_g = 35.0f * AI_DTOR;

            glintStrength = params[p_glintStrength].FLT;

            absorption = (AI_RGB_WHITE - min(rgb(0.9999f, 0.9999f, 0.9999f), params[p_hairColor].RGB)) * params[p_hairDensity].FLT;

            dual_depth = params[p_dualDepth].INT;

        }

        AtSampler* sampler_diffuse;
        AtSampler* sampler_R;
        AtSampler* sampler_TT;
        AtSampler* sampler_TRT;
        AtSampler* sampler_g;

        int dual_depth;

        AtRGB a_bar_f[DS_NUMSTEPS];
        AtRGB a_bar_b[DS_NUMSTEPS];
        AtRGB A_b[DS_NUMSTEPS];
        AtRGB alpha_f[DS_NUMSTEPS];
        AtRGB alpha_b[DS_NUMSTEPS];
        AtRGB beta_f[DS_NUMSTEPS];
        AtRGB beta_b[DS_NUMSTEPS];
        AtRGB sigma_b[DS_NUMSTEPS];
        AtRGB delta_b[DS_NUMSTEPS];
        AtRGB N_G_R[DS_NUMSTEPS*DS_NUMSTEPS];
        AtRGB N_G_TT[DS_NUMSTEPS*DS_NUMSTEPS];
        AtRGB N_G_TRT[DS_NUMSTEPS*DS_NUMSTEPS];
        AtRGB kf_R[DS_NUMSTEPS*DS_NUMSTEPS];
        AtRGB kf_TT[DS_NUMSTEPS*DS_NUMSTEPS];
        AtRGB kf_TRT[DS_NUMSTEPS*DS_NUMSTEPS];

        float beta_R;       //< R width
        float beta_R2;
        float alpha_R;      //< R shift
        float beta_TT;      //< TT width
        float beta_TT2;
        float alpha_TT;     //< TT shift
        float beta_TRT;     //< TRT width
        float beta_TRT2;
        float alpha_TRT;    //< TRT shift
        float gamma_TT;     //< TT rolloff
        float gamma_g;      //< g rolloff
        float phi_g;        //< g separation

        float glintStrength;

        AtRGB absorption;

        float ior;
    };

    HairBsdf(AtNode* n, AtShaderGlobals* sg, ShaderData* d) :
    node(n), data(d), numBlendHairs(2), density_front(0.7f), density_back(0.7f), _sg(sg)
    {
        depth = sg->Rr;

        // Get a local coordinate frame based on the hair fibre direction
        U = AiV3Normalize(sg->dPdv);
        V = AiV3Cross(U, sg->N);
        W = AiV3Cross(V, U);

        result_R_direct = AI_RGB_BLACK;
        result_R_indirect = AI_RGB_BLACK;
        result_TT_direct = AI_RGB_BLACK;
        result_TT_indirect = AI_RGB_BLACK;
        result_TRT_direct = AI_RGB_BLACK;
        result_TRT_indirect = AI_RGB_BLACK;
        result_TRTg_direct = AI_RGB_BLACK;
        result_TRTg_indirect = AI_RGB_BLACK;
        result_Pg_direct = AI_RGB_BLACK;
        result_Pl_direct = AI_RGB_BLACK;
        result_Pg_indirect = AI_RGB_BLACK;
        result_Pl_indirect = AI_RGB_BLACK;
        result_id1 = AI_RGB_BLACK;
        result_id2 = AI_RGB_BLACK;
        result_id3 = AI_RGB_BLACK;
        result_id4 = AI_RGB_BLACK;
        result_id5 = AI_RGB_BLACK;
        result_id6 = AI_RGB_BLACK;
        result_id7 = AI_RGB_BLACK;
        result_id8 = AI_RGB_BLACK;
    }

    /// Parameter evaluation. This should be called after opacity() and before anything else.
    inline void evaluateParameters(AtShaderGlobals* sg, ShaderData* data)
    {
        // Get the spherical angles of the exitant direction relative to the hair fibre
        wo = -sg->Rd;
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
        else
        {
            cn = 0.5f;
        }

        ior = data->ior;

        beta_R = data->beta_R;
        alpha_R = data->alpha_R;

        beta_TT = data->beta_TT;
        alpha_TT = data->alpha_TT;

        beta_TRT = data->beta_TRT;
        alpha_TRT = data->alpha_TRT;

        beta_R2 = data->beta_R2;
        beta_TT2 = data->beta_TT2;
        beta_TRT2 = data->beta_TRT2;

        gamma_TT = data->gamma_TT;
        gamma_g = data->gamma_g;
        phi_g = lerp(30.0f*AI_DTOR, 45.0f*AI_DTOR, cn);

        AB(theta_r, alpha_R, beta_R, A_R, B_R);
        AB(theta_r, alpha_TT, beta_TT, A_TT, B_TT);
        AB(theta_r, alpha_TRT, beta_TRT, A_TRT, B_TRT);

        absorption = data->absorption;
        specular1Color = AiShaderEvalParamRGB(p_specular1Color) * AiShaderEvalParamFlt(p_specular1Strength);
        specular2Color = AiShaderEvalParamRGB(p_specular2Color) * AiShaderEvalParamFlt(p_specular2Strength);
        transmissionColor = AiShaderEvalParamRGB(p_transmissionColor) * AiShaderEvalParamFlt(p_transmissionStrength);
        glintStrength = data->glintStrength;

        density_front = AiShaderEvalParamFlt(p_densityFront);
        density_back = AiShaderEvalParamFlt(p_densityBack);

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

    inline AtVector sample_R(float u1, float u2)
    {
        float theta_i = sampleLong(u1, theta_r, alpha_R, beta_R, A_R, B_R);
        float phi = 2.0f * asinf(clamp(2.0f*u2 - 1.0f, -1.0f, 1.0f)) + AI_PI;
        AtVector wi;
        sphericalDirection(theta_i, phi, V, W, U, wi);
        return wi;
    }

    struct SctGeo
    {
        SctGeo(const AtVector& w, float theta_r, float phi_r, const AtVector& U, const AtVector& V, const AtVector& W)
        {
            wi = w;
            theta_i = (AI_PIOVER2 - sphericalTheta(wi, U));
            cos_theta_i = cosf(theta_i);
            theta_h = (theta_r+theta_i)*0.5f;
            phi_i = sphericalPhi(wi, V, W);
            phi_d = phi_r - phi_i;
            if (phi_d < 0) phi_d += AI_PITIMES2;
            phi = phi_d - AI_PI;
            phi = AI_PI - fabsf(phi);
            phi_h = (phi_r+phi_i)*0.5f;
            cosphi2 = cosf(phi*0.5f);
            theta_d = (theta_r - theta_i)*0.5f;
            cos_theta_d = cosf(theta_d);
            inv_cos_theta_d2 = std::max(0.001f, 1.0f/(cos_theta_d*cos_theta_d));
        }

        AtVector wi;
        float theta_i;
        float cos_theta_i;
        float theta_h;
        float phi_i;
        float phi_d;
        float phi_h;
        float phi;
        float cosphi2;
        float theta_d;
        float cos_theta_d;
        float inv_cos_theta_d2;
    };

    inline AtRGB bsdf_R(const SctGeo& geo)
    {
        return rgb(bsdfR(beta_R2, alpha_R, geo.theta_h, geo.cosphi2)) * geo.cos_theta_i * geo.inv_cos_theta_d2 * AI_ONEOVER2PI;
    }

    inline float pdf_R(const SctGeo& geo)
    {
        float t = geo.theta_h-alpha_R;
        float pdf_theta = (1.0f / (2.0f*geo.cos_theta_i*(A_R-B_R))) * (beta_R / (t*t + beta_R2));
        float pdf_phi = geo.cosphi2*0.25f;
        return pdf_theta * pdf_phi;
    }

    inline AtVector sample_TRT(float u1, float u2)
    {
        float theta_i = sampleLong(u1, theta_r, alpha_TRT, beta_TRT, A_R, B_R);
        float phi = 2.0f * asinf(clamp(2.0f*u2 - 1.0f, -1.0f, 1.0f)) + AI_PI;
        AtVector wi;
        sphericalDirection(theta_i, phi, V, W, U, wi);
        return wi;
    }

    inline AtRGB bsdf_TRT(const SctGeo& geo)
    {
        return rgb(bsdfR(beta_TRT2, alpha_TRT, geo.theta_h, geo.cosphi2)) * geo.cos_theta_i * geo.inv_cos_theta_d2 * AI_ONEOVER2PI;;
    }

    inline float pdf_TRT(const SctGeo& geo)
    {
        float t = geo.theta_h-alpha_TRT;
        float pdf_theta = (1.0f / (2.0f*geo.cos_theta_i*(A_R-B_R))) * (beta_TRT / (t*t + beta_TRT2));
        float pdf_phi = geo.cosphi2*0.25f;
        return pdf_theta * pdf_phi;
    }

    inline AtVector sample_TRTg(float u1, float u2)
    {
        float theta_i = sampleLong(u1, theta_r, alpha_TRT, beta_TRT, A_TRT, B_TRT);
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
        float phi = gamma_g * tanf(u2*(Cg-Dg)+Dg) + phi_g;
        phi *= sign;

        AtVector wi;
        sphericalDirection(theta_i, phi, V, W, U, wi);
        return wi;
    }

    inline AtRGB bsdf_TRTg(const SctGeo& geo)
    {
        return rgb(bsdfg(beta_TRT2, alpha_TRT, geo.theta_h, gamma_g, geo.phi, phi_g)) * geo.cos_theta_i * geo.inv_cos_theta_d2 * AI_ONEOVER2PI;;
    }

    inline float pdf_TRTg(const SctGeo& geo)
    {   
        float t = geo.theta_h-alpha_TRT;
        float pdf_theta = (1.0f / (2.0f*geo.cos_theta_i*(A_TRT-B_TRT))) * (beta_TRT / (t*t + beta_TRT2));
        float p = fabs(geo.phi) - phi_g;
        float Cg = atanf((AI_PIOVER2 - phi_g)/gamma_g);
        float Dg = atanf(-phi_g/gamma_g);
        float pdf_phi = (1.0f / (2.0f * (Cg-Dg))) * (gamma_g / (p*p + gamma_g*gamma_g));
        return pdf_theta * pdf_phi;
    }

    inline AtVector sample_TT(float u1, float u2)
    {
        float theta_i = sampleLong(u1, theta_r, alpha_TT, beta_TT, A_TT, B_TT);
        C_TT = 2.0f * atanf(AI_PI/gamma_TT);
        float phi = gamma_TT * tanf(C_TT * (u2-0.5f)) + AI_PI;
        AtVector wi;
        sphericalDirection(theta_i, phi, V, W, U, wi);
        return wi;
    }

    inline AtRGB bsdf_TT(const SctGeo& geo)
    {
        return rgb(bsdfTT(beta_R2, alpha_R, geo.theta_h, gamma_TT, geo.phi_d)) * geo.cos_theta_i * geo.inv_cos_theta_d2 * AI_ONEOVER2PI;
    }

    inline float pdf_TT(const SctGeo& geo)
    {
        float t = geo.theta_h-alpha_TT;
        float pdf_theta = (1.0f / (2.0f*geo.cos_theta_i*(A_TT-B_TT))) * (beta_TT / (t*t + beta_TT2));
        float p = geo.phi-AI_PI;
        float C_TT = 2.0f * atanf(AI_PI/gamma_TT);
        float pdf_phi = (1.0f / C_TT) * (gamma_TT / (p*p + gamma_TT*gamma_TT));
        
        return pdf_theta * pdf_phi;
    }

    static AtVector HairGlossySample(const void* brdf_data, float u1, float u2)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        if (u1 < 0.5f && u2 < 0.5f) 
        {
            u1 *= 2.0f;
            u2 *= 2.0f;
            return hb->sample_R(u1, u2);
        }
        else if (u1 < 0.5f && u2 >= 0.5f)
        {
            u1 *= 2.0f;
            u2 = 2.0f * (1.0f-u2);
            return hb->sample_TT(u1, u2);
        }
        else if (u1 >= 0.5f && u2 < 0.5f)
        {
            u1 = 2.0f * (1.0f-u1);
            u2 *= 2.0f;
            return hb->sample_TRT(u1, u2);
        }
        else
        {
            u2 = 2.0f * (1.0f-u1);
            u2 = 2.0f * (1.0f-u2);
            return hb->sample_TRTg(u1, u2);
        }
    }

    static AtRGB HairGlossyBsdf(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        AtRGB result = AI_RGB_BLACK;

        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);

        AtRGB kfr[3];// = {rgb(0.1), rgb(0.9), rgb(0.1)};
        hairAttenuation(hb->ior, geo.theta_d, geo.phi_d, hb->absorption, kfr);
        result += hb->bsdf_R(geo) * kfr[0] * hb->specular1Color;// * rgb(1,0,0);
        result += hb->bsdf_TT(geo) * kfr[1] * hb->transmissionColor;// * rgb(0,1,0);
        result += hb->bsdf_TRT(geo) * kfr[2] * hb->specular2Color;// * rgb(0,0,1);
        result += hb->bsdf_TRTg(geo) * kfr[2] * hb->specular2Color;// * rgb(0,0,1);

        return result;
    }

    static float HairGlossyPdf(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;

        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);
        return (hb->pdf_R(geo) + hb->pdf_TT(geo) + hb->pdf_TRT(geo) + hb->pdf_TRTg(geo))*0.25f;
    }

   
    /// Integrate the direct illumination for all diffuse and glossy lobes
    inline void integrateDirect(AtShaderGlobals* sg)
    {
        if (!do_glossy) return;
        // Tell Arnold we want the full sphere for lighting.
        sg->fhemi = false;
        //sg->skip_shadow = true;
        AiLightsPrepare(sg);
        AtRGB kfr[4];
        while (AiLightsGetSample(sg))
        {
            result_R_direct += AiEvaluateLightSample(sg, this, HairGlossySample, HairGlossyBsdf, HairGlossyPdf);
        }
        //sg->skip_shadow = false;
        //result_R_direct *= specular1Color;// * ONEOVER4PI;
        sg->fhemi = true;        
    }


    /// Integrate the indirect illumination for all diffuse and glossy lobes
    inline void integrateIndirect(AtShaderGlobals* sg)
    {
        AiMakeRay(&wi_ray, AI_RAY_GLOSSY, &sg->P, NULL, AI_BIG, sg);
        sampit = AiSamplerIterator(data->sampler_R, sg);
        while(AiSamplerGetSample(sampit, samples))
        {
            wi_ray.dir = HairGlossySample(this, samples[0], samples[1]);

            AtScrSample scrs;

            // trace our ray
            if (AiTrace(&wi_ray, &scrs))
            {
                // calculate result
                float p = HairGlossyPdf(this, &wi_ray.dir);
                result_R_indirect += scrs.color * HairGlossyBsdf(this, &wi_ray.dir) / p;
            }
        }
        float weight = AiSamplerGetSampleInvCount(sampit);
        result_R_indirect *= weight * specular1Color;
    }



    inline void writeResult(AtShaderGlobals* sg)
    {
        if (sg->Rt & AI_RAY_CAMERA)
        {
            AiAOVSetRGB(sg, "direct_specular", result_R_direct);
            AiAOVSetRGB(sg, "indirect_specular", result_R_indirect);
            AiAOVSetRGB(sg, "direct_specular2", result_TRT_direct);
            AiAOVSetRGB(sg, "indirect_specular2", result_TRT_indirect);
            AiAOVSetRGB(sg, "direct_glint", result_TRTg_direct);
            AiAOVSetRGB(sg, "indirect_glint", result_TRTg_indirect);
            AiAOVSetRGB(sg, "direct_transmission", result_TT_direct);
            AiAOVSetRGB(sg, "indirect_transmission", result_TT_indirect);
            AiAOVSetRGB(sg, "direct_global", result_Pg_direct);
            AiAOVSetRGB(sg, "direct_local", result_Pl_direct);
            AiAOVSetRGB(sg, "indirect_global", result_Pg_indirect);
            AiAOVSetRGB(sg, "indirect_local", result_Pl_indirect);
            AiAOVSetRGB(sg, "id_1", result_id1);
            AiAOVSetRGB(sg, "id_2", result_id2);
            AiAOVSetRGB(sg, "id_3", result_id3);
            AiAOVSetRGB(sg, "id_4", result_id4);
            AiAOVSetRGB(sg, "id_5", result_id5);
            AiAOVSetRGB(sg, "id_6", result_id6);
            AiAOVSetRGB(sg, "id_7", result_id7);
            AiAOVSetRGB(sg, "id_8", result_id8);


        }

        sg->out.RGB =   result_R_direct +
                        result_R_indirect +
                        result_TT_direct +
                        result_TT_indirect +
                        result_TRT_direct +
                        result_TRT_indirect +
                        result_TRTg_direct +
                        result_TRTg_indirect +
                        result_Pg_direct +
                        result_Pl_direct +
                        result_Pg_indirect +
                        result_Pl_indirect;
        sg->out_opacity = AI_RGB_WHITE;
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

    float ior;
    float beta_R;       //< R width
    float beta_R2;
    float alpha_R;      //< R shift
    float beta_TT;      //< TT width
    float beta_TT2;
    float alpha_TT;     //< TT shift
    float beta_TRT;     //< TRT width
    float beta_TRT2;
    float alpha_TRT;    //< TRT shift
    float gamma_TT;     //< TT rolloff
    float gamma_g;      //< g rolloff
    float phi_g;        //< g separation

    AtRGB diffuseColor;
    AtRGB specular1Color;
    AtRGB specular2Color;
    AtRGB transmissionColor;
    float glintStrength;

    float A_R;
    float B_R;
    float A_TT;
    float B_TT;
    float C_TT;
    float A_TRT;
    float B_TRT;
    float Cg;
    float Dg;

    // dual-scattering parameters
    int numBlendHairs;
    float density_front;
    float density_back;

    AtRGB result_R_direct;
    AtRGB result_R_indirect;
    AtRGB result_TT_direct;
    AtRGB result_TT_indirect;
    AtRGB result_TRT_direct;
    AtRGB result_TRT_indirect;
    AtRGB result_TRTg_direct;
    AtRGB result_TRTg_indirect;

    // dual-scattering aovs
    AtRGB result_Pg_direct;
    AtRGB result_Pl_direct;
    AtRGB result_Pg_indirect;
    AtRGB result_Pl_indirect;

    // IDs
    AtRGB result_id1;
    AtRGB result_id2;
    AtRGB result_id3;
    AtRGB result_id4;
    AtRGB result_id5;
    AtRGB result_id6;
    AtRGB result_id7;
    AtRGB result_id8;

    bool do_diffuse;
    bool do_glossy;

    int depth;

    AtNode* node;
    ShaderData* data;
    AtRay wi_ray;
    AtScrSample scrs;
    double samples[2];
    AtSamplerIterator* sampit;

    AtVector wo;

    AtRGB absorption;

    AtShaderGlobals* _sg;
};


node_parameters
{
    AiParameterFlt("ior", 1.55f);
    AiParameterFlt("hairColorDensity", 1.f);
    AiParameterRGB("hairColor", 0.97f, 0.93f, 0.85f);
    AiParameterFlt("specularShift", 7.0f);
    AiParameterFlt("specularWidth", 5.0f);
    AiParameterInt("extraSamples", 0);
    AiParameterFlt("specular1Strength", 1.0f);
    AiParameterRGB("specular1Color", 1.0f, 1.0f, 1.0f);
    AiParameterFlt("specular2Strength", 1.0f);
    AiParameterRGB("specular2Color",1.0f, 1.0f, 1.0f);
    AiParameterFlt("glintStrength", 2.0f);
    AiParameterFlt("glintRolloff", 5.0f);
    AiParameterFlt("glintSeparation", 35.0f);
    AiParameterFlt("transmissionStrength", 1.0f);
    AiParameterRGB("transmissionColor", 1.0f, 1.0f, 1.0f);
    AiParameterFlt("transmissionRolloff", 30.0f);
    AiParameterRGB("opacity", 1.0f, 1.0f, 1.0f);
    AiParameterInt("dualDepth", 0);
    AiParameterFlt("densityFront", 0.7f);
    AiParameterFlt("densityBack", 0.7f);
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
    data->update(params);

    AiAOVRegister("direct_specular", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("indirect_specular", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("direct_specular2", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("indirect_specular2", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("direct_transmission", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("indirect_transmission", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("direct_glint", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("direct_global", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("indirect_global", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("direct_local", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("indirect_local", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_1", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_2", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_3", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_4", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_5", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_6", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_7", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_8", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
}

#define DUAL

shader_evaluate
{
    
    // Get shader data
    HairBsdf::ShaderData* data = (HairBsdf::ShaderData*)AiNodeGetLocalData(node);

    // Create HairBsdf object 
    HairBsdf hb(node, sg, data);
    // Get parameters
    hb.evaluateParameters(sg, data);

    AtRGB opacity = AiShaderEvalParamRGB(p_opacity);
    float geo_opacity = 1.0f;
    if (AiUDataGetFlt("geo_opacity", &geo_opacity))
    {
        opacity *= geo_opacity;
    }

    // early-out regardless if we're in a shadow ray, or if opacity is zero
    if (sg->Rt & AI_RAY_SHADOW || AiShaderGlobalsApplyOpacity(sg, opacity)) return; 

    hb.integrateDirect(sg);
    hb.integrateIndirect(sg);

    // Write shader result
    hb.writeResult(sg);
}


