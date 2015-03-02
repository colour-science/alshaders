#pragma once

#include <ai.h>
#include "alUtil.h"
#include <cassert>

#define SSS_MAX_SAMPLES 8

struct ScatteringParamsDiffusion
{
    ScatteringParamsDiffusion(const AtRGB& _sigma_s_prime, const AtRGB& _sigma_a, float _g, float _eta)
    : sigma_s_prime(_sigma_s_prime), sigma_a(_sigma_a), g(_g), eta(_eta)
    {
        Fdr = -1.44f/(eta*eta) + 0.71f/eta + 0.668f + 0.0636f*eta;
        float A = (1.0f + Fdr)/(1.0f - Fdr);
        sigma_t_prime = sigma_s_prime + sigma_a;
        alpha_prime = sigma_s_prime / sigma_t_prime;
        AtRGB l = AI_RGB_WHITE / sigma_t_prime;
        AtRGB D = l / 3.0f;

        sigma_tr = sqrt(3.0f*sigma_a*sigma_t_prime);
        zb = 2.0f * A * D;
        zr = AI_RGB_WHITE / sigma_t_prime;
        zv = zr + (2.0f*zb);

        AtRGB sq = sqrt(3.0f * (AI_RGB_WHITE - alpha_prime));
        albedo = alpha_prime * 0.5f * (AI_RGB_WHITE + fast_exp(-A*sq*4.0f/3.0f)) * fast_exp(-sq);
        float total = albedo.r + albedo.g + albedo.b;
        albedo_norm = albedo / total;
    }

    AtRGB sigma_s_prime;
    AtRGB sigma_t_prime;
    AtRGB sigma_a;
    float g;
    float eta;
    float Fdr;
    AtRGB sigma_tr;
    AtRGB zr;
    AtRGB zv;
    AtRGB zb;
    AtRGB alpha_prime;
    AtRGB albedo;
    AtRGB albedo_norm;
};

struct DiffusionSample
{
    AtVector P;
    AtVector N;
    AtVector Ng;
    AtRGB R;
    AtRGB L;

};

struct DiffusionMessageData
{
    int sss_depth;
    float maxdist;
    AtVector Po;
    AtVector No;
    AtVector U;
    AtVector V;
    ScatteringParamsDiffusion sp;
    AtShaderGlobals* sg;
    DiffusionSample samples[SSS_MAX_SAMPLES];
};

inline float diffusionSampleDisk(float u1, float u2, float sigma, float& dx, float& dy, float& r)
{
    r = -logf(1.0f - std::min(0.999999f,u1)) / sigma;
    float phi = u2 * AI_PITIMES2;
    dx = r * cosf(phi);
    dy = r * sinf(phi);

    // return pdf
    return sigma * fast_exp(-sigma * r);
}

inline float diffusionPdf(float r, float sigma)
{
    return sigma * fast_exp(-sigma * r);
}

inline float dipoleProfileRd(float r, float sigma_tr, float zr, float zv)
{
    float dr = sqrtf(r*r + zr*zr);
    float dv = sqrtf(r*r + zv*zv);
    dr = std::max(dr, zr);
    dv = std::max(dv, zv);

    float inv_dr3 = 1.0f / (dr*dr*dr);
    float inv_dv3 = 1.0f / (dv*dv*dv);

    float sigma_tr_dr = sigma_tr * dr;
    float sigma_tr_dv = sigma_tr * dv;

    return zr*(sigma_tr_dr+1.0f) * fast_exp(-sigma_tr_dr) * inv_dr3 + zv*(sigma_tr_dv+1.0f) * fast_exp(-sigma_tr_dv) * inv_dv3;
}

inline AtRGB dipoleProfileRd(float r, const AtRGB& sigma_tr, const AtRGB& zr, const AtRGB& zv)
{
    return rgb(
        dipoleProfileRd(r, sigma_tr.r, zr.r, zv.r),
        dipoleProfileRd(r, sigma_tr.g, zr.g, zv.g),
        dipoleProfileRd(r, sigma_tr.b, zr.b, zv.b)
    );
}

inline AtRGB dipole(float r, const AtVector& N, const AtVector& Nx, const ScatteringParamsDiffusion& sp)
{
    return dipoleProfileRd(r, sp.sigma_tr, sp.zr, sp.zv) * sp.alpha_prime;
}


void alsIrradiateSample(AtShaderGlobals* sg, DiffusionMessageData* dmd)
{
    AiStateSetMsgInt("als_raytype", ALS_RAY_UNDEFINED);

    // void* brdf_data = AiOrenNayarMISCreateData(sg, 0.0f);
    // sg->fhemi = false;
    AiLightsPrepare(sg);
    AtRGB result_diffuse = AI_RGB_BLACK;
    AtUInt32 old_fi = sg->fi;
    while (AiLightsGetSample(sg))
    {
        result_diffuse += sg->Li * MAX(AiV3Dot(sg->Ld, sg->N), 0.0f) * AI_ONEOVERPI * sg->we;
        // can't use MIS here because Arnold cocks up the shadowing ;__;
        // result_diffuse += AiEvaluateLightSample(sg, brdf_data, AiOrenNayarMISSample, AiOrenNayarMISBRDF, AiOrenNayarMISPDF);
    }
    dmd->samples[dmd->sss_depth].L = result_diffuse;

    dmd->samples[dmd->sss_depth].N = sg->N;
    dmd->samples[dmd->sss_depth].Ng = sg->Ng;
    dmd->samples[dmd->sss_depth].P = sg->P;

    dmd->sss_depth++;

    dmd->maxdist -= AiV3Length(sg->P - dmd->Po);

    if (dmd->sss_depth < SSS_MAX_SAMPLES && dmd->maxdist > 0.0f)
    {
        AiStateSetMsgInt("als_raytype", ALS_RAY_SSS);
        
        AtRay ray;
        AtScrSample scrs;
        sg->Rr--;
        AiMakeRay(&ray, AI_RAY_SUBSURFACE, &sg->P, &sg->Rd, dmd->maxdist, sg);
        AiTrace(&ray, &scrs);
    }
}

AtRGB alsDiffusion(AtShaderGlobals* sg, DiffusionMessageData* dmd, AtSampler* sss_sampler, float sssDensityScale)
{
    AtVector U, V;
    AiBuildLocalFrameShirley(&U, &V, &sg->Ng);

    // Get scattering parameters from supplied scattering colour and mfp
    // AtRGB sigma_s_prime, sigma_a;
    // alphaInversion(sssRadiusColor, sssRadius, sigma_s_prime, sigma_a);
    AtRGB sigma_s_prime = rgb(1.09f, 1.59f, 1.79f);
    AtRGB sigma_a = rgb(0.013, 0.07f, 0.145f);

    // The 10 / AI_PIOVER2 factor is there to roughly match the behaviour of the cubic. And if you're
    // going to use a magic number it should always involve pi somewhere...
    sigma_s_prime *= sssDensityScale*20.0f / AI_PIOVER2;
    sigma_a *= sssDensityScale*20.0f / AI_PIOVER2;
    float eta = 1.4f;

    ScatteringParamsDiffusion sp(sigma_s_prime, sigma_a, 0.7f, eta);

    // Find the max component of the mfp
    float l = maxh(sp.zr);

    // trick Arnold into thinking we're shooting from a different face than we actually are so he doesn't ignore intersections
    AtUInt32 old_fi = sg->fi;
    sg->fi = UINT_MAX;
    
    
    AtRGB result_sss = AI_RGB_BLACK;
    
    // Set our maximum sample distance to be some multiple of the mfp
    float R_max = l * 25.0f;
    
    AtRGB Rd_sum = AI_RGB_BLACK;
    int samplesTaken = 0;
    float samples[2];
    AtRay wi_ray;
    AtScrSample scrs;
    AtSamplerIterator* sampit = AiSamplerIterator(sss_sampler, sg);
    while (AiSamplerGetSample(sampit, samples))
    {
        // TODO: replace with a better sampling scheme
        //wi_ray.dir = uniformSampleSphere(samples[0], samples[1]);
        float dx, dy;

        AtVector Wsss, Usss, Vsss, Usss_1, Vsss_1, Usss_2, Vsss_2;
        float c_axis = 1.0f, c_axis_1, c_axis_2;;
        if (samples[0] < 0.5f)
        {
            samples[0] *= 2.0f;
            c_axis = 0.5f;
            c_axis_1 = 0.25f;
            c_axis_2 = 0.25f;
            Wsss = sg->Ng;
            
            Usss = U;
            Vsss = V;

            Usss_1 = sg->Ng;
            Vsss_1 = V;

            Usss_2 = U;
            Vsss_2 = sg->Ng;
        }
        else if (samples[0] < 0.75f)
        {
            samples[0] = (samples[0] - 0.5f) * 4.0f;
            c_axis = 0.25f;
            c_axis_1 = 0.5f;
            c_axis_2 = 0.25f;
            Wsss = U;
            
            Usss = sg->Ng;
            Vsss = V;

            Usss_1 = U;
            Vsss_1 = V;

            Usss_2 = U;
            Vsss_2 = sg->Ng;
        }
        else
        {
            samples[0] = (1.0f-samples[0])* 4.0f;
            c_axis = 0.25f;
            c_axis_1 = 0.25f;
            c_axis_2 = 0.5f;
            Wsss = V;
            
            Usss = U;
            Vsss = sg->Ng;

            Usss_1 = sg->Ng;
            Vsss_1 = V;

            Usss_2 = U;
            Vsss_2 = V;
        }

        AtVector Pd;
        float pdf_disk_a0_c0, pdf_disk_1, pdf_disk_2;
        float r_disk;
        float c_disk = 1.0f, c_disk_1, c_disk_2;
        float sigma, sigma_1, sigma_2;
        if (samples[1] < sp.albedo_norm.r)
        {
            samples[1] /= sp.albedo_norm.r;
            sigma = sp.sigma_tr.r;
            sigma_1 = sp.sigma_tr.g;
            sigma_2 = sp.sigma_tr.b;
            pdf_disk_a0_c0 = diffusionSampleDisk(samples[0], samples[1], sigma, dx, dy, r_disk);

            c_disk = sp.albedo_norm.r;
            c_disk_1 = sp.albedo_norm.g;
            c_disk_2 = sp.albedo_norm.b;
        }
        else if (samples[1] < (sp.albedo_norm.r + sp.albedo_norm.g))
        {
            samples[1] -= sp.albedo_norm.r;
            samples[1] /= sp.albedo_norm.g;
            sigma = sp.sigma_tr.g;
            sigma_1 = sp.sigma_tr.r;
            sigma_2 = sp.sigma_tr.b;
            pdf_disk_a0_c0 = diffusionSampleDisk(samples[0], samples[1], sp.sigma_tr.g, dx, dy, r_disk);
            c_disk = sp.albedo_norm.g;
            c_disk_1 = sp.albedo_norm.r;
            c_disk_2 = sp.albedo_norm.b;
        }
        else
        {
            samples[1] = 1.0f - samples[1];
            samples[1] /= sp.albedo_norm.b;
            sigma = sp.sigma_tr.b;
            sigma_1 = sp.sigma_tr.g;
            sigma_2 = sp.sigma_tr.r;
            pdf_disk_a0_c0 = diffusionSampleDisk(samples[0], samples[1], sp.sigma_tr.b, dx, dy, r_disk);   
            c_disk = sp.albedo_norm.b;
            c_disk_1 = sp.albedo_norm.g;
            c_disk_2 = sp.albedo_norm.r;
        }

        AtVector dir = -Wsss;
        float dz = R_max;
        AtPoint origin = sg->P + Wsss*(dz*0.5f) + Usss * dx + Vsss * dy;
        float maxdist = R_max * 2.0f;

        AiMakeRay(&wi_ray, AI_RAY_SUBSURFACE, &origin, &dir, maxdist, sg);

        AiStateSetMsgInt("als_raytype", ALS_RAY_SSS);
        dmd->sss_depth = 0;
        dmd->maxdist = R_max;
        dmd->Po = sg->P;
        dmd->No = sg->N;
        if (AiTrace(&wi_ray, &scrs))
        {            
            for (int i=0; i < dmd->sss_depth; ++i)
            {
                AtRGB Rd;
                
                AtVector R = dmd->samples[i].P - sg->P;
                float r = AiV3Length(R);
                
                float geom = fabsf(AiV3Dot(dmd->samples[i].Ng, Wsss));
                float geom_1 = fabsf(AiV3Dot(dmd->samples[i].Ng, Usss));
                float geom_2 = fabsf(AiV3Dot(dmd->samples[i].Ng, Vsss));

                // if (data->numExtraPoles == 0) Rd = dipole(r, sg->N, dmd->samples[i].N, sp);
                // else Rd = multipole(r, sg->N, dmd->samples[i].N, sp, .1f, data->numExtraPoles);

                Rd = dipole(r, sg->N, dmd->samples[i].N, sp);

                float r_u_1 = AiV3Dot(R, Usss_1);
                float r_v_1 = AiV3Dot(R, Vsss_1);
                float r_1 = sqrtf(SQR(r_u_1)+SQR(r_v_1));

                float r_u_2 = AiV3Dot(R, Usss_2);
                float r_v_2 = AiV3Dot(R, Vsss_2);
                float r_2 = sqrtf(SQR(r_u_2)+SQR(r_v_2));

                float pdf_disk_a0_c0 = diffusionPdf(r_disk, sigma) * geom;
                float pdf_disk_a0_c1 = diffusionPdf(r_disk, sigma_1) * geom;
                float pdf_disk_a0_c2 = diffusionPdf(r_disk, sigma_2) * geom;

                float pdf_disk_a1_c0 = diffusionPdf(r_1, sigma) * geom_1;
                float pdf_disk_a1_c1 = diffusionPdf(r_1, sigma_1) * geom_1;
                float pdf_disk_a1_c2 = diffusionPdf(r_1, sigma_2) * geom_1;

                float pdf_disk_a2_c0 = diffusionPdf(r_2, sigma) * geom_2;
                float pdf_disk_a2_c1 = diffusionPdf(r_2, sigma_1) * geom_2;
                float pdf_disk_a2_c2 = diffusionPdf(r_2, sigma_2) * geom_2;

                float pdf_sum = 
                    pdf_disk_a0_c0 * c_disk * c_axis +
                    pdf_disk_a0_c1 * c_disk_1 * c_axis +
                    pdf_disk_a0_c2 * c_disk_2 * c_axis +
                    pdf_disk_a1_c0 * c_disk * c_axis_1 +
                    pdf_disk_a1_c1 * c_disk_1 * c_axis_1 +
                    pdf_disk_a1_c2 * c_disk_2 * c_axis_1 +
                    pdf_disk_a2_c0 * c_disk * c_axis_2 +
                    pdf_disk_a2_c1 * c_disk_1 * c_axis_2 +
                    pdf_disk_a2_c2 * c_disk_2 * c_axis_2;

                result_sss += dmd->samples[i].L * Rd * r_disk / pdf_sum;
                Rd_sum += Rd * r_disk / pdf_sum;
                
            }
            
        }
    }
    float w = AiSamplerGetSampleInvCount(sampit) * 0.5f;
    result_sss *= w;
    Rd_sum *= w;

    // if (minh(Rd_sum) != 0.0f) result_sss *= (sp.albedo/Rd_sum); //< Peter Kutz's dark-edge normalization trick
    result_sss /= sp.albedo;
    sg->fi = old_fi;

    // Optimization hack: do a regular indirect diffuse and colour it like the subsurface instead of allowing sss rays to
    // continue for another diffuse bounce.
    /*
    if (!sssDoIndirect)
    {
        AtRGB result_sss_indirect = AiIndirectDiffuse(&sg->N, sg); 
        if (!data->sss_normalize) result_sss_indirect *= sp.albedo;
        result_sss += result_sss_indirect;
    }
    */

    return result_sss;
}
