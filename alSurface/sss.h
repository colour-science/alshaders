#pragma once

#include <ai.h>
#include "alUtil.h"
#include <cassert>

#define SSS_MAX_SAMPLES 8

struct ScatteringParamsDipole
{
    ScatteringParamsDipole(const AtRGB& _sigma_s_prime, const AtRGB& _sigma_a, float _g, float _eta)
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

struct ScatteringParamsDirectional
{
    ScatteringParamsDirectional(AtRGB sigma_s_prime, AtRGB sigma_a, float g) :
    g(g), eta(1.3f), 
    C_phi(0.175626f),
    C_phi_inv(1.03668f),
    C_E(0.27735f),
    A(2.05736f)
    {
        AtRGB sigma_s = sigma_s_prime/(1.0f-g);

        sigma_t_prime = sigma_s_prime + sigma_a;
        sigma_t = sigma_s + sigma_a;

        alpha_prime = sigma_s_prime / sigma_t_prime;

        float apsum = alpha_prime.r + alpha_prime.g + alpha_prime.b;
        albedo_norm = alpha_prime / apsum;

        AtRGB sq = sqrt(3.0f * (AI_RGB_WHITE - alpha_prime));
        albedo = alpha_prime * 0.5f * (AI_RGB_WHITE + fast_exp(-A*sq*4.0f/3.0f)) * fast_exp(-sq);

        D = (2*sigma_t_prime) / (3*SQR(sigma_t_prime));
        sigma_tr = sqrt(sigma_a / D);
        de = 2.131 * D / sqrt(alpha_prime);
        zr = AI_RGB_WHITE / sigma_t_prime;
    }

    float g;
    float eta;
    float C_phi;
    float C_phi_inv;
    float C_E;
    float A;
    AtRGB de;
    AtRGB sigma_t;
    AtRGB sigma_t_prime;
    AtRGB sigma_tr;
    AtRGB D;
    AtRGB zr;
    AtRGB alpha_prime;
    AtRGB albedo_norm;
    AtRGB albedo;
};

struct DiffusionSample
{
    AtVector P;     //< sampled point
    AtVector N;     //< normal at sample
    AtVector Ng;    //< geometric normal at sample
    AtRGB R;        //< diffusion result at sample
    AtRGB Rd;       //< convolved diffusion result
    AtRGB E;        //< irradiance at sample
    AtVector S;     //< vector from shading point to sample
    float r;        //< distance from shading point to sample
    float b;        //< bounce attenuation factor
};

struct DiffusionMessageData
{
    int sss_depth;
    float maxdist;
    AtVector Po;
    AtVector No;
    AtVector U;
    AtVector V;
    ScatteringParamsDipole sp;
    AtShaderGlobals* sg;
    DiffusionSample samples[SSS_MAX_SAMPLES];
};

struct DirectionalMessageData
{
    int sss_depth;
    float maxdist;
    AtVector Po;
    AtVector No;
    AtVector U;
    AtVector V;
    AtVector wo;
    ScatteringParamsDirectional sp;
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

inline AtRGB dipole(float r, const AtVector& N, const AtVector& Nx, const ScatteringParamsDipole& sp)
{
    return dipoleProfileRd(r, sp.sigma_tr, sp.zr, sp.zv) * sp.alpha_prime;
}

// Directional dipole profile evaluation
// see: http://www.ci.i.u-tokyo.ac.jp/~hachisuka/dirpole.pdf
// and: http://www.ci.i.u-tokyo.ac.jp/~hachisuka/dirpole.cpp
inline float Sp_d(const AtVector& x, const AtVector& w, const float r, const AtVector& n, const AtRGB& sigma_tr, const AtRGB& D, const float Cp_norm, 
                    const float Cp, const float Ce, const int j) 
{
    // evaluate the profile
    const float s_tr_r = sigma_tr[j] * r;
    const float s_tr_r_one = 1.0f + s_tr_r;
    const float x_dot_w = AiV3Dot(x, w);
    const float r_sqr = r * r;

    const float t0 = Cp_norm * (1.0f / (4.0f * AI_PI * AI_PI)) * expf(-s_tr_r) / (r * r_sqr);
    const float t1 = r_sqr / D[j] + 3.0f * s_tr_r_one * x_dot_w;
    const float t2 = 3.0f * D[j] * s_tr_r_one * AiV3Dot(w, n);
    const float t3 = (s_tr_r_one + 3.0f * D[j] * (3.0f * s_tr_r_one + s_tr_r * s_tr_r) / r_sqr * x_dot_w) * AiV3Dot(x, n);

    return t0 * (Cp * t1 - Ce * (t2 - t3));
}

float directionalDipoleRd(AtPoint xi, AtVector ni, AtPoint xo, AtVector no, 
                          AtVector wi, AtVector wo, float eta, AtRGB de, AtRGB D, AtRGB sigma_t, AtRGB sigma_tr, 
                          float A, float C_phi, float C_phi_inv, float C_E, unsigned int j)
{
    // distance
    AtVector xoxi = xo - xi;
    float r = AiV3Length(xoxi);

    // modified normal
    AtVector ni_s = AiV3Cross(AiV3Normalize(xoxi), AiV3Normalize(AiV3Cross(ni, xoxi)));

    // directions of ray sources
    float nnt = 1.0f / eta;
    float ddn = -AiV3Dot(wi, ni);
    AtVector wr = AiV3Normalize(wi * -nnt - ni * (ddn * nnt + sqrtf(1.0f - nnt*nnt * (1.0f - ddn*ddn))));
    AtVector wv = wr - ni_s * (2.0f * AiV3Dot(wr, ni_s));

    // distance to real sources
    const float cos_beta = -sqrtf((r * r - AiV3Dot(xoxi, wr) * AiV3Dot(xoxi, wr)) / (r * r + de[j] * de[j]));
    float dr;
    const float mu0 = -AiV3Dot(no, wr);
    if (mu0 > 0.0) 
    {
        dr = sqrtf((D[j] * mu0) * ((D[j] * mu0) - de[j] * cos_beta * 2.0) + r * r);
    } 
    else 
    {
        dr = sqrtf(1.0f / (3.0f * sigma_t[j] * 3.0f * sigma_t[j]) + r * r);
    }

    AtVector xoxv = xo - (xi + ni_s * (2.0f * A * de[j]));
    float dv = AiV3Length(xoxv);

    float result = Sp_d(xoxi, wr, dr, no, sigma_tr, D, C_phi_inv, C_phi, C_E, j) 
                 - Sp_d(xoxv, wv, dv, no, sigma_tr, D, C_phi_inv, C_phi, C_E, j);

    return std::max(0.0f, result);
}

AtRGB directionalDipole(AtPoint xi, AtVector ni, AtPoint xo, AtVector no, AtVector wi, AtVector wo, 
                        const ScatteringParamsDirectional& sp)
{
    return rgb(
        directionalDipoleRd(xi, ni, xo, no, wi, wo, sp.eta, sp.de, sp.D, sp.sigma_t, sp.sigma_tr, sp.A, sp.C_phi, sp.C_phi_inv, sp.C_E, 0),
        directionalDipoleRd(xi, ni, xo, no, wi, wo, sp.eta, sp.de, sp.D, sp.sigma_t, sp.sigma_tr, sp.A, sp.C_phi, sp.C_phi_inv, sp.C_E, 1),
        directionalDipoleRd(xi, ni, xo, no, wi, wo, sp.eta, sp.de, sp.D, sp.sigma_t, sp.sigma_tr, sp.A, sp.C_phi, sp.C_phi_inv, sp.C_E, 2)
    );
}

void alsIrradiateSample(AtShaderGlobals* sg, DirectionalMessageData* dmd)
{
    AiStateSetMsgInt("als_raytype", ALS_RAY_UNDEFINED);
    DiffusionSample& samp = dmd->samples[dmd->sss_depth];
    // void* brdf_data = AiOrenNayarMISCreateData(sg, 0.0f);
    // sg->fhemi = false;
    AiLightsPrepare(sg);
    AtRGB result_diffuse = AI_RGB_BLACK;
    AtUInt32 old_fi = sg->fi;
    samp.Rd = AI_RGB_BLACK;
    while (AiLightsGetSample(sg))
    {
        result_diffuse += sg->Li * MAX(AiV3Dot(sg->Ld, sg->N), 0.0f) * AI_ONEOVERPI * sg->we;
        // can't use MIS here because Arnold cocks up the shadowing ;__;
        // result_diffuse += AiEvaluateLightSample(sg, brdf_data, AiOrenNayarMISSample, AiOrenNayarMISBRDF, AiOrenNayarMISPDF);
        samp.Rd += directionalDipole(sg->P, sg->N, dmd->Po, dmd->No, sg->Ld, dmd->wo, dmd->sp)
                    * sg->Li * MAX(AiV3Dot(sg->Ld, sg->N), 0.0f) * sg->we;
    }

    // result_diffuse += AiIndirectDiffuse(&sg->N, sg);

    

    samp.E = result_diffuse;
    samp.S = sg->P - dmd->Po;
    samp.r = AiV3Length(samp.S);
    if (!AiColorIsZero(result_diffuse))
    {
        // samp.R = dipole(samp.r, dmd->No, sg->N, dmd->sp);
        samp.b = 1.0f ;//- (1.0f - AiV3Dot(sg->N, dmd->No) * (1.0f - AiV3Dot(dmd->No, AiV3Normalize(samp.S))))*0.25f;
    }

    samp.N = sg->N;
    samp.Ng = sg->Ng;
    samp.P = sg->P;

    dmd->sss_depth++;

    dmd->maxdist -= samp.r;

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

AtRGB alsDiffusion(AtShaderGlobals* sg, DirectionalMessageData* dmd, AtSampler* sss_sampler, 
                     AtRGB sssRadiusColor, float sssRadius, float sssDensityScale)
{
    AtVector U, V;
    AiBuildLocalFrameShirley(&U, &V, &sg->Ng);

    // Get scattering parameters from supplied scattering colour and mfp
    AtRGB sigma_s_prime, sigma_a;
#if 1 
    alphaInversion(sssRadiusColor, sssRadius, sigma_s_prime, sigma_a);
    sigma_s_prime *= sssDensityScale * AI_PIOVER2;
    sigma_a *= sssDensityScale * AI_PIOVER2;
#else
    sigma_s_prime = rgb(1.09f, 1.59f, 1.79f);
    sigma_a = rgb(0.013, 0.07f, 0.145f);

    // The 10 / AI_PIOVER2 factor is there to roughly match the behaviour of the cubic. And if you're
    // going to use a magic number it should always involve pi somewhere...
    sigma_s_prime *= sssDensityScale*20.0f / AI_PIOVER2;
    sigma_a *= sssDensityScale*20.0f / AI_PIOVER2;
#endif    
    float eta = 1.3f;

    ScatteringParamsDirectional sp(sigma_s_prime, sigma_a, 0.0f);
    dmd->sp = sp;
#if 0
    static bool first = true;
    if (first)
    {
        printf("alpha_prime: (%f, %f, %f)\n", sp.alpha_prime.r, sp.alpha_prime.g, sp.alpha_prime.b);
        printf("de: (%f, %f, %f)\n", sp.de.r, sp.de.g, sp.de.b);
        printf("D: (%f, %f, %f)\n", sp.D.r, sp.D.g, sp.D.b);
        printf("sigma_t_prime: (%f, %f, %f)\n", sp.sigma_t_prime.r, sp.sigma_t_prime.g, sp.sigma_t_prime.b);
        AtRGB st = exp( - sp.sigma_t_prime * .002f);
        printf("st: (%f, %f, %f)\n", st.r, st.g, st.b);
        first = false;
    }
#endif 

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
    dmd->wo = -sg->Rd;
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
                if (AiColorIsZero(dmd->samples[i].E)) continue;

                float geom = fabsf(AiV3Dot(dmd->samples[i].Ng, Wsss));
                float geom_1 = fabsf(AiV3Dot(dmd->samples[i].Ng, Usss));
                float geom_2 = fabsf(AiV3Dot(dmd->samples[i].Ng, Vsss));

                float r_u_1 = AiV3Dot(dmd->samples[i].S, Usss_1);
                float r_v_1 = AiV3Dot(dmd->samples[i].S, Vsss_1);
                float r_1 = sqrtf(SQR(r_u_1)+SQR(r_v_1));

                float r_u_2 = AiV3Dot(dmd->samples[i].S, Usss_2);
                float r_v_2 = AiV3Dot(dmd->samples[i].S, Vsss_2);
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

                // result_sss += dmd->samples[i].E * dmd->samples[i].R * dmd->samples[i].b * r_disk / pdf_sum;
                result_sss += dmd->samples[i].Rd * dmd->samples[i].b * r_disk / pdf_sum;
                // Rd_sum += dmd->samples[i].R * dmd->samples[i].b * r_disk / pdf_sum;
                
            }
            
        }
    }
    float w = AiSamplerGetSampleInvCount(sampit) * 0.5f;
    result_sss *= w;
    // Rd_sum *= w;

    // if (!AiColorIsZero(Rd_sum)) result_sss *= (sp.albedo/Rd_sum); //< Peter Kutz's dark-edge normalization trick
    // result_sss /= sp.albedo;
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
