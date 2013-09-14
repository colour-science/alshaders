// Hair shader based on 
// [1] ISHair: Importance Sampling for Hair Scattering by Ou et al. 2012
// http://www.cs.dartmouth.edu/~ouj/site/Research/Entries/2012/6/21_ISHair__Importance_Sampling_for_Hair_Scattering.html
// [2] Dual Scattering Approximation For Fast Multiple Scattering in Hair by Zinke et al. 2008


#include <ai.h>
#include "alUtil.h"
#include "stats.h"
#include "scattering.h"
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

#define B_WIDTH_SCALE 2.0f
#define F_WIDTH_SCALE 4.0f

struct HairBsdf
{
    struct ShaderData
    {
        ShaderData()
        : sampler_R(NULL), sampler_TT(NULL), sampler_TRT(NULL), sampler_g(NULL), ds(NULL)
        {}

        ~ShaderData()
        {
            AiSamplerDestroy(sampler_R);
            AiSamplerDestroy(sampler_TT);
            AiSamplerDestroy(sampler_TRT);
            AiSamplerDestroy(sampler_g);
            delete ds;
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

            dual_depth = params[p_dualDepth].INT;

            sampleLobesIndividually = false;

            delete ds;
            ds = new DualScattering;
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

        bool sampleLobesIndividually;

        DualScattering* ds;
        float colorNorm;
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

        sp.ior = data->ior;

        sp.beta_R = data->beta_R;
        sp.alpha_R = data->alpha_R;

        sp.beta_TT = data->beta_TT;
        sp.alpha_TT = data->alpha_TT;

        sp.beta_TRT = data->beta_TRT;
        sp.alpha_TRT = data->alpha_TRT;

        sp.beta_R2 = data->beta_R2;
        sp.beta_TT2 = data->beta_TT2;
        sp.beta_TRT2 = data->beta_TRT2;

        sp.gamma_TT = data->gamma_TT;
        sp.gamma_g = data->gamma_g;
        sp.phi_g = lerp(30.0f*AI_DTOR, 45.0f*AI_DTOR, cn);

        AB(theta_r, sp.alpha_R, sp.beta_R, A_R, B_R);
        AB(theta_r, sp.alpha_TT, sp.beta_TT, A_TT, B_TT);
        AB(theta_r, sp.alpha_TRT, sp.beta_TRT, A_TRT, B_TRT);
        AB(theta_r, sp.alpha_R, sp.beta_R*B_WIDTH_SCALE, A_b, B_b);
        AB(theta_r, sp.alpha_TT, sp.beta_TT*F_WIDTH_SCALE, A_f, B_f);

        AtRGB hairColor = min(rgb(0.99f, 0.99f, 0.99f), pow(AiShaderEvalParamRGB(p_hairColor), 0.15f));
        colorNorm = maxh(hairColor);

        sp.absorption = (AI_RGB_WHITE - hairColor);

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

    inline float sampleLong(double u, float theta_r, float alpha, float beta, float A, float B)
    {
        float t = beta * tanf(u*(A-B) + B);
        float theta_h = t + alpha;
        return clamp( -0.4999f * AI_PI, 0.4999f * AI_PI, (2.0f*theta_h - theta_r));
    }

    inline AtVector sample_R(float u1, float u2)
    {
        float theta_i = sampleLong(u1, theta_r, sp.alpha_R, sp.beta_R, A_R, B_R);
        float phi = 2.0f * asinf(clamp(2.0f*u2 - 1.0f, -1.0f, 1.0f));// + AI_PI;
        float phi_i = phi_r - phi;
        AtVector wi;
        sphericalDirection(theta_i, phi_i, V, W, U, wi);
        return wi;
    }

    inline AtRGB bsdf_R(const SctGeo& geo)
    {
        return rgb(bsdfR(sp.beta_R2, sp.alpha_R, geo.theta_h, geo.cosphi2)) * geo.cos_theta_i * geo.inv_cos_theta_d2 * AI_ONEOVER2PI;
    }

    inline float pdf_R(const SctGeo& geo)
    {
        float t = geo.theta_h-sp.alpha_R;
        float pdf_theta = (1.0f / (2.0f*geo.cos_theta_i*(A_R-B_R))) * (sp.beta_R / (t*t + sp.beta_R2));
        float pdf_phi = geo.cosphi2*0.25f;
        return pdf_theta * pdf_phi;
    }

    inline AtVector sample_TRT(float u1, float u2)
    {
        float theta_i = sampleLong(u1, theta_r, sp.alpha_TRT, sp.beta_TRT, A_R, B_R);
        float phi = 2.0f * asinf(clamp(2.0f*u2 - 1.0f, -1.0f, 1.0f));// + AI_PI;
        float phi_i = phi_r - phi;
        AtVector wi;
        sphericalDirection(theta_i, phi_i, V, W, U, wi);
        return wi;
    }

    inline AtRGB bsdf_TRT(const SctGeo& geo)
    {
        return rgb(bsdfR(sp.beta_TRT2, sp.alpha_TRT, geo.theta_h, geo.cosphi2)) * geo.cos_theta_i * geo.inv_cos_theta_d2 * AI_ONEOVER2PI;;
    }

    inline float pdf_TRT(const SctGeo& geo)
    {
        float t = geo.theta_h-sp.alpha_TRT;
        float pdf_theta = (1.0f / (2.0f*geo.cos_theta_i*(A_R-B_R))) * (sp.beta_TRT / (t*t + sp.beta_TRT2));
        float pdf_phi = geo.cosphi2*0.25f;
        return pdf_theta * pdf_phi;
    }

    inline AtVector sample_TRTg(float u1, float u2)
    {
        float theta_i = sampleLong(u1, theta_r, sp.alpha_TRT, sp.beta_TRT, A_TRT, B_TRT);
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
        Cg = atanf((AI_PIOVER2 - sp.phi_g)/sp.gamma_g);
        Dg = atanf(-sp.phi_g/sp.gamma_g);
        float phi = sp.gamma_g * tanf(u2*(Cg-Dg)+Dg) + sp.phi_g;
        phi *= sign;
        float phi_i = phi_r - phi;

        AtVector wi;
        sphericalDirection(theta_i, phi_i, V, W, U, wi);
        return wi;
    }

    inline AtRGB bsdf_TRTg(const SctGeo& geo)
    {
        return rgb(bsdfg(sp.beta_TRT2, sp.alpha_TRT, geo.theta_h, sp.gamma_g, geo.phi, sp.phi_g)) * geo.cos_theta_i * geo.inv_cos_theta_d2 * AI_ONEOVER2PI;;
    }

    inline float pdf_TRTg(const SctGeo& geo)
    {   
        float t = geo.theta_h-sp.alpha_TRT;
        float pdf_theta = (1.0f / (2.0f*geo.cos_theta_i*(A_TRT-B_TRT))) * (sp.beta_TRT / (t*t + sp.beta_TRT2));
        float p = fabs(geo.phi) - sp.phi_g;
        float Cg = atanf((AI_PIOVER2 - sp.phi_g)/sp.gamma_g);
        float Dg = atanf(-sp.phi_g/sp.gamma_g);
        float pdf_phi = (1.0f / (2.0f * (Cg-Dg))) * (sp.gamma_g / (p*p + sp.gamma_g*sp.gamma_g));
        return pdf_theta * pdf_phi;
    }

    inline AtVector sample_TT(float u1, float u2)
    {
        float theta_i = sampleLong(u1, theta_r, sp.alpha_TT, sp.beta_TT, A_TT, B_TT);
        C_TT = 2.0f * atanf(AI_PI/sp.gamma_TT);
        float phi = sp.gamma_TT * tanf(C_TT * (u2-0.5f)) + AI_PI;
        float phi_i = phi_r - phi;
        AtVector wi;
        sphericalDirection(theta_i, phi_i, V, W, U, wi);
        return wi;
    }

    inline AtRGB bsdf_TT(const SctGeo& geo)
    {
        return rgb(bsdfTT(sp.beta_R2, sp.alpha_R, geo.theta_h, sp.gamma_TT, geo.phi_d)) * geo.cos_theta_i * geo.inv_cos_theta_d2 * AI_ONEOVER2PI;
    }

    inline float pdf_TT(const SctGeo& geo)
    {
        float t = geo.theta_h-sp.alpha_TT;
        float pdf_theta = (1.0f / (2.0f*geo.cos_theta_i*(A_TT-B_TT))) * (sp.beta_TT / (t*t + sp.beta_TT2));
        float p = geo.phi-AI_PI;
        float C_TT = 2.0f * atanf(AI_PI/sp.gamma_TT);
        float pdf_phi = (1.0f / C_TT) * (sp.gamma_TT / (p*p + sp.gamma_TT*sp.gamma_TT));
        
        return pdf_theta * pdf_phi;
    }

    inline AtVector sample_Sb(float u1, float u2)
    {
        float theta_i = sampleLong(u1, theta_r, sp.alpha_R, sp.beta_R*B_WIDTH_SCALE, A_b, B_b);
        float phi = 2.0f * asinf(clamp(2.0f*u2 - 1.0f, -1.0f, 1.0f));// + AI_PI;
        float phi_i = phi_r - phi;
        AtVector wi;
        sphericalDirection(theta_i, phi_i, V, W, U, wi);
        return wi;
    }

    inline float pdf_Sb(const SctGeo& geo)
    {
        float t = geo.theta_h-sp.alpha_R*B_WIDTH_SCALE;
        float pdf_theta = (1.0f / (2.0f*geo.cos_theta_i*(A_b-B_b))) * (sp.beta_R*B_WIDTH_SCALE / (t*t + sp.beta_R2*B_WIDTH_SCALE*B_WIDTH_SCALE));
        float pdf_phi = geo.cosphi2*0.25f;
        return pdf_theta * pdf_phi;
    }

    inline AtVector sample_Sf(float u1, float u2)
    {
        float theta_i = sampleLong(u1, theta_r, sp.alpha_TT, sp.beta_TT*F_WIDTH_SCALE, A_f, B_f);
        C_TT = 2.0f * atanf(AI_PI/(sp.gamma_TT*F_WIDTH_SCALE));
        float phi = sp.gamma_TT*F_WIDTH_SCALE * tanf(C_TT * (u2-0.5f)) + AI_PI;
        float phi_i = phi_r - phi;
        AtVector wi;
        sphericalDirection(theta_i, phi_i, V, W, U, wi);
        return wi;
    }

    inline float pdf_Sf(const SctGeo& geo)
    {
        float t = geo.theta_h-sp.alpha_TT;
        float pdf_theta = (1.0f / (2.0f*geo.cos_theta_i*(A_f-B_f))) * (sp.beta_TT*F_WIDTH_SCALE / (t*t + sp.beta_TT2*F_WIDTH_SCALE*F_WIDTH_SCALE));
        float p = geo.phi-AI_PI;
        float C_TT = 2.0f * atanf(AI_PI/(sp.gamma_TT*F_WIDTH_SCALE));
        float pdf_phi = (1.0f / C_TT) * (sp.gamma_TT*F_WIDTH_SCALE / (p*p + sp.gamma_TT*sp.gamma_TT*F_WIDTH_SCALE*F_WIDTH_SCALE));
        
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

    static AtVector Hair_Sample_R(const void* brdf_data, float u1, float u2)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        return hb->sample_R(u1, u2);
    }

    static AtVector Hair_Sample_TT(const void* brdf_data, float u1, float u2)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        return hb->sample_TT(u1, u2);
    }

    static AtVector Hair_Sample_TRT(const void* brdf_data, float u1, float u2)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        return hb->sample_TRT(u1, u2);
    }

    static AtVector Hair_Sample_TRTg(const void* brdf_data, float u1, float u2)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        return hb->sample_TRTg(u1, u2);
    }

    static AtRGB HairGlossyBsdf(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        AtRGB result = AI_RGB_BLACK;

        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);

        AtRGB kfr[3];
        hairAttenuation(hb->sp.ior, geo.theta_d, geo.phi_d, hb->sp.absorption, kfr);
 
        result += hb->bsdf_R(geo) * kfr[0] * hb->specular1Color;
        result += hb->bsdf_TT(geo) * kfr[1] * hb->transmissionColor;
        result += hb->bsdf_TRT(geo) * kfr[2] * hb->specular2Color;
        result += hb->bsdf_TRTg(geo) * kfr[2] * hb->specular2Color;

        return result;
    }

    static AtRGB Hair_Bsdf_R(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);
        AtRGB kfr[3];
        hairAttenuation(hb->sp.ior, geo.theta_d, geo.phi_d, hb->sp.absorption, kfr);
        return hb->bsdf_R(geo) * kfr[0] * hb->specular1Color;
    }
    static AtRGB Hair_Bsdf_TT(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);
        AtRGB kfr[3];
        hairAttenuation(hb->sp.ior, geo.theta_d, geo.phi_d, hb->sp.absorption, kfr);
        return hb->bsdf_TT(geo) * kfr[1] * hb->transmissionColor;
    }
    static AtRGB Hair_Bsdf_TRT(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);
        AtRGB kfr[3];
        hairAttenuation(hb->sp.ior, geo.theta_d, geo.phi_d, hb->sp.absorption, kfr);
        return hb->bsdf_TRT(geo) * kfr[2] * hb->specular2Color;
    }
    static AtRGB Hair_Bsdf_TRTg(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);
        AtRGB kfr[3];
        hairAttenuation(hb->sp.ior, geo.theta_d, geo.phi_d, hb->sp.absorption, kfr);
        return hb->bsdf_TRTg(geo) * kfr[2] * hb->specular2Color;
    }

    static float HairGlossyPdf(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;

        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);
        return (hb->pdf_R(geo) + hb->pdf_TT(geo) + hb->pdf_TRT(geo) + hb->pdf_TRTg(geo))*0.25f;
    }

    static float Hair_Pdf_R(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);
        return hb->pdf_R(geo);
    }

    static float Hair_Pdf_TT(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);
        return hb->pdf_TT(geo);
    }

    static float Hair_Pdf_TRT(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);
        return hb->pdf_TRT(geo);
    }

    static float Hair_Pdf_TRTg(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);
        return hb->pdf_TRTg(geo);
    }

   
    /// Integrate the direct illumination for all diffuse and glossy lobes
    inline void integrateDirectMis(AtShaderGlobals* sg)
    {
        if (!do_glossy) return;
        // Tell Arnold we want the full sphere for lighting.
        sg->fhemi = false;
        //sg->skip_shadow = true;
        AiLightsPrepare(sg);

        if ((sg->Rt & AI_RAY_CAMERA) && data->sampleLobesIndividually)
        {
            while (AiLightsGetSample(sg))
            {
                result_R_direct += AiEvaluateLightSample(sg, this, Hair_Sample_R, Hair_Bsdf_R, Hair_Pdf_R);
                result_TT_direct += AiEvaluateLightSample(sg, this, Hair_Sample_TT, Hair_Bsdf_TT, Hair_Pdf_TT);
                result_TRT_direct += AiEvaluateLightSample(sg, this, Hair_Sample_TRT, Hair_Bsdf_TRT, Hair_Pdf_TRT);
                result_TRTg_direct += AiEvaluateLightSample(sg, this, Hair_Sample_TRTg, Hair_Bsdf_TRTg, Hair_Pdf_TRTg);
            }
        }
        else
        {
            while (AiLightsGetSample(sg))
            {
                result_R_direct += AiEvaluateLightSample(sg, this, HairGlossySample, HairGlossyBsdf, HairGlossyPdf);
            }
        }
        sg->fhemi = true;        
    }

    inline void integrateDirectDual(AtShaderGlobals* sg)
    {
        bool do_glossy = true;
        if (sg->Rt & AI_RAY_DIFFUSE) do_glossy = false;

        int als_hairNumIntersections = 0;
        AtRGB als_T_f = AI_RGB_BLACK;
        AtRGB als_sigma_bar_f = AI_RGB_BLACK;
        bool old_hemi = sg->fhemi;
        sg->fhemi = false;
        bool old_skipshadow = sg->skip_shadow;
        sg->skip_shadow = true;
        AiLightsPrepare(sg);
        AtRay ray;
        AtScrSample scrs;
        AiMakeRay(&ray, AI_RAY_SHADOW, &(sg->P), NULL, AI_BIG, sg);
        AtRGB kfr[3];
        AiStateSetMsgInt("als_raytype", ALS_RAY_DUAL);
        while (AiLightsGetSample(sg))
        {
            SctGeo geo(sg->Ld, theta_r, phi_r, U, V, W);
            
            AiStateSetMsgInt("als_hairNumIntersections", 0);
            AiStateSetMsgRGB("als_T_f", AI_RGB_WHITE);
            AiStateSetMsgRGB("als_sigma_bar_f", AI_RGB_BLACK);
            ray.dir = sg->Ld;
            AiTrace(&ray, &scrs);
            AiStateGetMsgInt("als_hairNumIntersections", &als_hairNumIntersections);
            AiStateGetMsgRGB("als_T_f", &als_T_f);
            AiStateGetMsgRGB("als_sigma_bar_f", &als_sigma_bar_f);

            AtRGB f_direct_back = data->ds->direct_back_scatter(sp, geo.theta_i, geo.theta_h, als_sigma_bar_f) * geo.inv_cos_theta_d2;
            
            float directFraction = 1.0f - std::min(als_hairNumIntersections, numBlendHairs)/float(numBlendHairs);
            AtRGB occlusion = AI_RGB_WHITE - scrs.opacity;

            if (directFraction > 0.0f)
            {
                
                result_Pl_direct += sg->Li * sg->we * occlusion * density_back * f_direct_back * geo.cos_theta_i * directFraction;

                if (!AiIsFinite(result_Pl_direct))
                {
                    std::cerr << VAR(result_Pl_direct) << std::endl;
                    std::cerr << VAR(f_direct_back) << std::endl;
                    std::cerr << VAR(als_T_f) << std::endl;
                    std::cerr << VAR(als_sigma_bar_f) << std::endl;
                }
                
                if (do_glossy)
                {
                    AtRGB kfr[3];
                    hairAttenuation(sp.ior, geo.theta_d, geo.phi_d, sp.absorption, kfr);
                    kfr[2] *= lerp(0.1f, 3.0f, SQR(cn)); //< glints
                    AtRGB L = sg->Li * sg->we * occlusion * directFraction;
                    
                    result_R_direct += L * bsdf_R(geo) * kfr[0];
                    result_TT_direct += L * bsdf_TT(geo) * kfr[1];
                    result_TRT_direct += L * bsdf_TRT(geo) * kfr[2];
                    result_TRTg_direct += L * bsdf_TRTg(geo) * kfr[2];
                }

            }
            if (directFraction < 1.0f)
            {
                AtRGB T_f = als_T_f;
                AtRGB S_f = g(als_sigma_bar_f, AI_RGB_BLACK, rgb(geo.theta_h)) / (AI_PI * geo.cos_theta_d);

                AtRGB f_s_scatter = AI_RGB_BLACK;
        
                if (geo.phi >= AI_PIOVER2 ) // forward scattering directions only
                {
                    f_s_scatter = data->ds->forward_scatter(sp, geo.theta_h, geo.theta_d, geo.phi, als_sigma_bar_f);
                    f_s_scatter *= S_f;
                }
                
                AtRGB f_scatter_back = data->ds->back_scatter(sp, geo.theta_i, geo.theta_h) * geo.inv_cos_theta_d2;

                AtRGB F_scatter = T_f * density_front * (f_s_scatter + AI_PI*density_back*f_scatter_back);

                result_Pg_direct += sg->Li * sg->we * occlusion * F_scatter * geo.cos_theta_i * (1.0f-directFraction);

                if (!AiIsFinite(result_Pg_direct))
                {
                    std::cerr << VAR(result_Pg_direct) << std::endl;
                    std::cerr << VAR(F_scatter) << std::endl;
                    std::cerr << VAR(f_scatter_back) << std::endl;
                    std::cerr << VAR(T_f) << std::endl;
                    std::cerr << VAR(S_f) << std::endl;
                }

            }
        } // END light loop
        AiStateSetMsgInt("als_raytype", ALS_RAY_UNDEFINED);

        result_R_direct *= specular1Color;
        result_TT_direct *= transmissionColor;
        result_TRT_direct *= specular2Color;
        result_TRTg_direct *= specular2Color * glintStrength;

        sg->fhemi = old_hemi;
        sg->skip_shadow = old_skipshadow;
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
                result_R_indirect += scrs.color * Hair_Bsdf_R(this, &wi_ray.dir) / p;
                result_TT_indirect += scrs.color * Hair_Bsdf_TT(this, &wi_ray.dir) / p;
                result_TRT_indirect += scrs.color * Hair_Bsdf_TRT(this, &wi_ray.dir) / p;
                result_TRTg_indirect += scrs.color * Hair_Bsdf_TRTg(this, &wi_ray.dir) / p;
            }
        }
        float weight = AiSamplerGetSampleInvCount(sampit);
        result_R_indirect *= weight * AI_PI; //< TODO: factor of pi?
        result_TT_indirect *= weight * AI_PI; //< TODO: factor of pi?
        result_TRT_indirect *= weight * AI_PI; //< TODO: factor of pi?
        result_TRTg_indirect *= weight * AI_PI; //< TODO: factor of pi?
        
    }

    inline void integrateIndirectDual(AtShaderGlobals* sg)
    {
        AiMakeRay(&wi_ray, AI_RAY_GLOSSY, &sg->P, NULL, AI_BIG, sg);
        float weight;

        bool do_glossy = true;
        if (sg->Rt & AI_RAY_DIFFUSE || sg->Rr > 0) do_glossy = false;

        if (do_glossy)
        {
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
                    result_R_indirect += scrs.color * Hair_Bsdf_R(this, &wi_ray.dir) / p;
                    result_TT_indirect += scrs.color * Hair_Bsdf_TT(this, &wi_ray.dir) / p;
                    result_TRT_indirect += scrs.color * Hair_Bsdf_TRT(this, &wi_ray.dir) / p;
                    result_TRTg_indirect += scrs.color * Hair_Bsdf_TRTg(this, &wi_ray.dir) / p;
                }
            }
            weight = AiSamplerGetSampleInvCount(sampit);
            result_R_indirect *= weight * AI_PI; //< TODO: factor of pi?
            result_TT_indirect *= weight * AI_PI; //< TODO: factor of pi?
            result_TRT_indirect *= weight * AI_PI; //< TODO: factor of pi?
            result_TRTg_indirect *= weight * AI_PI; //< TODO: factor of pi?
        }

        int als_hairNumIntersections = 0;
        AtRGB als_T_f = AI_RGB_BLACK;
        AtRGB als_sigma_bar_f = AI_RGB_BLACK;
        sampit = AiSamplerIterator(data->sampler_TT, sg);
        AiStateSetMsgInt("als_raytype", ALS_RAY_DUAL);
        while (AiSamplerGetSample(sampit, samples))
        {
            if (samples[0] < 0.5f)
            {
                samples[0] *= 2.0f;
                wi_ray.dir = sample_Sb(samples[0], samples[1]);
            }
            else
            {
                samples[0] = 1.0f - (2.0f * samples[0]);
                wi_ray.dir = sample_Sf(samples[0], samples[1]);
            }

            SctGeo geo(wi_ray.dir, theta_r, phi_r, U, V, W);

            float p = 1.0f / ((pdf_Sf(geo) + pdf_Sb(geo)) * 0.5f);

            AiStateSetMsgInt("als_hairNumIntersections", 0);
            AiStateSetMsgRGB("als_T_f", AI_RGB_WHITE);
            AiStateSetMsgRGB("als_sigma_bar_f", AI_RGB_BLACK);
            AiTrace(&wi_ray, &scrs);
            AiStateGetMsgInt("als_hairNumIntersections", &als_hairNumIntersections);
            AiStateGetMsgRGB("als_T_f", &als_T_f);
            AiStateGetMsgRGB("als_sigma_bar_f", &als_sigma_bar_f);

            AtRGB f_direct_back = data->ds->direct_back_scatter(sp, geo.theta_i, geo.theta_h, als_sigma_bar_f) * geo.inv_cos_theta_d2;
            
            float directFraction = 1.0f - std::min(als_hairNumIntersections, numBlendHairs)/float(numBlendHairs);
            AtRGB occlusion = AI_RGB_WHITE - scrs.opacity;

            if (directFraction > 0.0f)
            {
                
                result_Pl_indirect += scrs.color * density_back * f_direct_back * geo.cos_theta_i * directFraction;

            }
            if (directFraction < 1.0f)
            {
                AtRGB T_f = als_T_f;
                AtRGB S_f = g(als_sigma_bar_f, AI_RGB_BLACK, rgb(geo.theta_h)) / (AI_PI * geo.cos_theta_d);

                AtRGB f_s_scatter = AI_RGB_BLACK;
        
                if (geo.phi >= AI_PIOVER2 ) // forward scattering directions only
                {
                    f_s_scatter = data->ds->forward_scatter(sp, geo.theta_h, geo.theta_d, geo.phi, als_sigma_bar_f);
                    f_s_scatter *= S_f;
                }
                
                AtRGB f_scatter_back = data->ds->back_scatter(sp, geo.theta_i, geo.theta_h) * geo.inv_cos_theta_d2;

                AtRGB F_scatter = T_f * density_front * (f_s_scatter + AI_PI*density_back*f_scatter_back);

                result_Pg_indirect += scrs.color * F_scatter * geo.cos_theta_i * (1.0f-directFraction);

            }

        }
        AiStateSetMsgInt("als_raytype", ALS_RAY_UNDEFINED);
        weight = AiSamplerGetSampleInvCount(sampit);
        result_Pl_indirect *= weight * AI_PI;
        result_Pg_indirect *= weight * AI_PI;
    }

#if 0
    inline void integrateIndirectDual(AtShaderGlobals* sg)
    {
        if (sg->Rt & AI_RAY_DIFFUSE) return;

        AiMakeRay(&wi_ray, AI_RAY_GLOSSY, &sg->P, NULL, AI_BIG, sg);
        AtScrSample scrs;
        sampit = AiSamplerIterator(data->sampler_R, sg);
        int als_hairNumIntersections = 0;
        AtRGB als_T_f = AI_RGB_BLACK;
        AtRGB als_sigma_bar_f = AI_RGB_BLACK;
        
        sampit = AiSamplerIterator(data->sampler_R, sg);
        AiStateSetMsgInt("als_raytype", ALS_RAY_DUAL);
        while(AiSamplerGetSample(sampit, samples))
        {
            wi_ray.dir = HairGlossySample(this, samples[0], samples[1]);
            float p = 1.0f / HairGlossyPdf(this, &wi_ray.dir);

            AtScrSample scrs;

            AiStateSetMsgInt("als_hairNumIntersections", 0);
            AiStateSetMsgRGB("als_T_f", AI_RGB_WHITE);
            AiStateSetMsgRGB("als_sigma_bar_f", AI_RGB_BLACK);
            AiTrace(&wi_ray, &scrs);
            AiStateGetMsgInt("als_hairNumIntersections", &als_hairNumIntersections);
            AiStateGetMsgRGB("als_T_f", &als_T_f);
            AiStateGetMsgRGB("als_sigma_bar_f", &als_sigma_bar_f);

            SctGeo geo(wi_ray.dir, theta_r, phi_r, U, V, W);

            AtRGB f_direct_back = data->ds->direct_back_scatter(sp, geo.theta_i, geo.theta_h, als_sigma_bar_f) * geo.inv_cos_theta_d2;
            float directFraction = 1.0f - std::min(als_hairNumIntersections, numBlendHairs)/float(numBlendHairs);
            if (directFraction > 0.0f)
            {
                AtRGB kfr[3];
                hairAttenuation(sp.ior, geo.theta_d, geo.phi_d, sp.absorption, kfr);
                result_Pl_indirect += scrs.color * density_back * f_direct_back * geo.cos_theta_i * directFraction;
                AtRGB L = scrs.color * p * directFraction;
                result_R_indirect += L * Hair_Bsdf_R(this, &wi_ray.dir) * kfr[0];
                result_TT_indirect += L * Hair_Bsdf_TT(this, &wi_ray.dir) * kfr[1];
                result_TRT_indirect += L * Hair_Bsdf_TRT(this, &wi_ray.dir) * kfr[2];
                result_TRTg_indirect += L * Hair_Bsdf_TRTg(this, &wi_ray.dir) * kfr[2];

            }
            if (directFraction < 1.0f)
            {
                AtRGB T_f = als_T_f;
                AtRGB S_f = g(als_sigma_bar_f, AI_RGB_BLACK, rgb(geo.theta_h)) / (AI_PI * geo.cos_theta_d);

                AtRGB f_s_scatter = AI_RGB_BLACK;

                if (geo.phi >= AI_PIOVER2) // forward scattering directions only
                {
                    f_s_scatter = data->ds->forward_scatter(sp, geo.theta_h, geo.theta_d, geo.phi, als_sigma_bar_f);
                    f_s_scatter *= S_f;
                }

                AtRGB f_scatter_back = data->ds->back_scatter(sp, geo.theta_i, geo.theta_h) * geo.inv_cos_theta_d2;

                float directFraction = 0.0f;
                AtRGB F_scatter = T_f * density_front * (f_s_scatter + AI_PI*density_back*f_scatter_back);

                result_Pg_indirect += scrs.color * F_scatter * geo.cos_theta_i * p * (1.0f-directFraction); 
            }
        }
        AiStateSetMsgInt("als_raytype", ALS_RAY_UNDEFINED);
        float weight = AiSamplerGetSampleInvCount(sampit);
        result_R_indirect *= weight * AI_PI; //< TODO: factor of pi?
        result_TT_indirect *= weight * AI_PI; //< TODO: factor of pi?
        result_TRT_indirect *= weight * AI_PI; //< TODO: factor of pi?
        result_TRTg_indirect *= weight * AI_PI; //< TODO: factor of pi?
        result_Pl_indirect *= weight * AI_PI; //< TODO: factor of pi?
        result_Pg_indirect *= weight * AI_PI; //< TODO: factor of pi?
        
    }
#endif

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

    ScatteringParams sp; //< varying scattering parameters

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

    float A_f;
    float A_b;
    float B_f;
    float B_b;

    // dual-scattering parameters
    int numBlendHairs;
    float density_front;
    float density_back;
    float colorNorm;

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

    int als_hairNumIntersections = 0;
    AtRGB als_T_f = AI_RGB_BLACK;
    AtRGB als_sigma_bar_f = AI_RGB_BLACK;
    bool do_dual = false;
    if (sg->Rr >= data->dual_depth) do_dual = true;

    int als_raytype = ALS_RAY_UNDEFINED;
    AiStateGetMsgInt("als_raytype", &als_raytype);
    
    if (do_dual && als_raytype == ALS_RAY_DUAL)
    {
        if (AiStateGetMsgInt("als_hairNumIntersections", &als_hairNumIntersections) 
            && AiStateGetMsgRGB("als_T_f", &als_T_f)
            && AiStateGetMsgRGB("als_sigma_bar_f", &als_sigma_bar_f))
        {
            float theta_i = AI_PIOVER2 - sphericalTheta(sg->Rd, hb.U);
            int idx = int((theta_i / AI_PI + 0.5f)*(DS_NUMSTEPS-1));

            //als_T_f *= hb.data->a_bar_f[idx];
            als_T_f *=  hb.data->ds->forward_attenuation(hb.sp, theta_i);
            AiStateSetMsgRGB("als_T_f", als_T_f);

            als_sigma_bar_f += hb.sp.beta_R2 + hb.sp.beta_TRT2 + hb.sp.beta_TT2;
            AiStateSetMsgRGB("als_sigma_bar_f", als_sigma_bar_f);

            als_hairNumIntersections++;
            AiStateSetMsgInt("als_hairNumIntersections", als_hairNumIntersections);
            sg->out_opacity = AI_RGB_BLACK;
        }
        else
        {
            sg->out_opacity = AI_RGB_WHITE;
        }

        return; // early out
    }

    // early-out regardless if we're in a shadow ray, or if opacity is zero
    if (sg->Rt & AI_RAY_SHADOW || AiShaderGlobalsApplyOpacity(sg, opacity)) return; 

    if (sg->Rr_gloss > data->dual_depth) return;
    
    // calculate scattering explicitly up to the dual depth cutoff.
    // in other words, do a brute force path trace for x=dual-depth bounces, then fall back to dual scattering for the rest.
    if (do_dual)
    {
        hb.integrateDirectDual(sg);
        hb.integrateIndirectDual(sg);
    }
    else
    {
        hb.integrateDirectMis(sg);
        hb.integrateIndirect(sg);
    }

    // Write shader result
    hb.writeResult(sg);

    //sg->out_opacity = rgb(0.5f);
}


