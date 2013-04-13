// Hair shader based on 
// [1] ISHair: Importance Sampling for Hair Scattering by Ou et al. 2012
// http://www.cs.dartmouth.edu/~ouj/site/Research/Entries/2012/6/21_ISHair__Importance_Sampling_for_Hair_Scattering.html
// [2] Dual Scattering Approximation For Fast Multiple Scattering in Hair by Zinke et al. 2008


#include <ai.h>
#include "alUtil.h"
#include <vector>
#include <algorithm>

AI_SHADER_NODE_EXPORT_METHODS(alHair)

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
    p_densityFront,
    p_densityBack
};

// hard-code IOR for now
// this is completely buggered
#define IOR 1.6f
#define FRESNEL_SAMPLES 1024

float fresnel(float incidenceAngle, float etaPerp, float etaParal, float invert)
{
        float n1, n2;
        float rPerp = 1;
        float rParal = 1;


        float angle = abs(incidenceAngle);
        if (angle > AI_PIOVER2)
        {
            angle = AI_PI - angle;
        }

        if ( invert )
        {
            n1 = etaPerp;
            n2 = 1;
        }
        else
        {
            n1 = 1;
            n2 = etaPerp;
        }

        // Perpendicular light reflectance
        float a = (n1/n2)*sin(angle);
        a *= a;
        if ( a <= 1 )
        {

            float b = n2*sqrt(1-a);
            float c = n1*cos(angle);
            rPerp =  ( c - b ) / ( c + b );
            rPerp *= rPerp;
            rPerp = std::min(1.0f, rPerp );
        }
        if ( invert )
        {
            n1 = etaParal;
            n2 = 1;
        }
        else
        {
            n1 = 1;
            n2 = etaParal;
        }
        // Parallel light reflectance
        float d = (n1/n2)*sin(angle);
        d *= d;
        if ( d <= 1 )
        {

            float e = n1*sqrt(1-d);
            float f = n2*cos(angle);
            rParal = ( e - f ) / ( e + f );
            rParal *= rParal;
            rParal = std::min( 1.0f, rParal );
        }
        return 0.5 * (rPerp + rParal);
}



void hairAttenuation(float ior, float cos_theta_i, float theta_d, float phi, AtRGB absorption, AtRGB kfr[3])
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
        if (p != 1)
        {
            if (phi_p > AI_PI) phi_p -= AI_PITIMES2;
            phi_p += p*AI_PI;
        }
        // get roots of polynomial
        float roots[3] = {0,0,0};
        int numRoots = solveCubic(-8.0f*p*c*ONEOVERPI3, 0.0f, 6.0f*p*c*AI_ONEOVERPI - 2.0f, p*AI_PI - phi_p, roots);
        AtRGB Fr = AI_RGB_BLACK;
        if (p < 2)
        {
            float gamma_i = roots[0];
            if (p==0)
            {
                Fr = rgb(fresnel(gamma_i, n_p, n_pp, false));
            }
            else
            {
                float gamma_t = asinf(sinf(gamma_i)/n_p);
                float theta_t = acosf(std::min(1.0f, (n_p/ior)*cos_theta_i));
                float cos_theta_t = cosf(theta_t);
                float l = 2.0f * cosf(gamma_t) / std::max(0.0001f, cos_theta_t);
                Fr = (1.0f - fresnel(gamma_i, n_p, n_pp, false)) * (1.0f - fresnel(gamma_t, n_p, n_pp, true)) * exp(-absorption*l);
            }
        }
        else 
        {
            for (int i=0; i < numRoots; ++i)
            {
                float gamma_i = roots[i];
                float gamma_t = asinf(sinf(gamma_i)/n_p);
                float theta_t = acosf(std::min(1.0f, (n_p/ior)*cos_theta_i));
                float cos_theta_t = cosf(theta_t);
                float l = 2.0f * cosf(gamma_t) / std::max(0.0001f, cos_theta_t);

                float iFr = fresnel(gamma_t, n_p, n_pp, true);
                Fr += (1.0f - fresnel(gamma_i, n_p, n_pp, false)) * (1.0f - iFr) * iFr * exp(-absorption*l);
            }
        }
        kfr[p] = Fr;
    }
    //kfr[0] = kfr[1] = kfr[2] = 0.1f;
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
    struct ScatteringParameters
    {
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

        AtRGB specular1Color;
        AtRGB specular2Color;
        AtRGB transmissionColor;
        float glintStrength;

        AtRGB absorption;
    };

    struct ShaderData
    {
        ShaderData()
        : sampler_diffuse(NULL), sampler_R(NULL), sampler_TT(NULL), sampler_TRT(NULL), sampler_g(NULL)
        {}

        ~ShaderData()
        {
            AiSamplerDestroy(sampler_diffuse);
            AiSamplerDestroy(sampler_R);
            AiSamplerDestroy(sampler_TT);
            AiSamplerDestroy(sampler_TRT);
            AiSamplerDestroy(sampler_g);
        }

        void update(const ScatteringParameters& sp, int diffuse_samples, int glossy_samples)
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

            // precalculate scattering LUTs
            float theta_i_step = AI_PIOVER2 / DS_NUMSTEPS;
            float theta_r_step = AI_PI/ DS_NUMSTEPS;
            float phi_step = AI_PIOVER2 / DS_NUMSTEPS;
            int idx = 0;
            memset(a_bar_f, 0, sizeof(AtRGB)*DS_NUMSTEPS);
            memset(a_bar_b, 0, sizeof(AtRGB)*DS_NUMSTEPS);
            memset(alpha_f, 0, sizeof(AtRGB)*DS_NUMSTEPS);
            memset(alpha_b, 0, sizeof(AtRGB)*DS_NUMSTEPS);
            memset(beta_f, 0, sizeof(AtRGB)*DS_NUMSTEPS);
            memset(beta_b, 0, sizeof(AtRGB)*DS_NUMSTEPS);
            AtRGB kfr[3];
            for (float theta_i = 0; theta_i < AI_PIOVER2; theta_i+=theta_i_step)
            {
                float cos_theta_i = cosf(theta_i);
                for (float theta_r = -AI_PIOVER2; theta_r < AI_PIOVER2; theta_r += theta_r_step)
                {
                    float theta_h = (theta_i + theta_r) * 0.5f;
                    float theta_d = (theta_i - theta_r) * 0.5f;
                    // integrate over half the domain each time and double after as it's symmetrical
                    // forward scattering
                    for (float phi = AI_PIOVER2; phi < AI_PI; phi += phi_step)
                    {
                        hairAttenuation(1.55, cos_theta_i, theta_d, phi, sp.absorption, kfr);
                        float cosphi2 = cosf(phi*0.5f);
                        
                        // [2] eq (6)
                        // Compute average forward-scattering attenuation, i.e. total forward-scattered radiance
                        // TODO: Should be multiplying by cos_theta_i here too?
                        // {
                        AtRGB f_R = bsdfR(sp.beta_R2, sp.alpha_R, theta_h, cosphi2) * sp.specular1Color * kfr[0] * cos_theta_i;
                        AtRGB f_TT = bsdfTT(sp.beta_TT2, sp.alpha_TT, theta_h, sp.gamma_TT, phi) * sp.transmissionColor * kfr[1] * cos_theta_i;
                        AtRGB f_TRT = (bsdfTRT(sp.beta_TRT2, sp.alpha_TRT, theta_h, cosphi2) * sp.specular2Color
                                    + bsdfg(sp.beta_TRT2, sp.alpha_TRT, theta_h, sp.gamma_g, phi, sp.phi_g) * sp.specular2Color * sp.glintStrength)
                                    * kfr[2] * cos_theta_i;
                        a_bar_f[idx] += f_R + f_TT + f_TRT;
                        
                        // }
                        alpha_f[idx] += f_R*sp.alpha_R + f_TT*sp.alpha_TT + f_TRT*sp.alpha_TRT;
                        beta_f[idx] += f_R*sp.beta_R2 + f_TT*sp.beta_TT2 + f_TRT*sp.beta_TRT2;
                    }

                    // backward scattering
                    for (float phi = 0.0f; phi < AI_PIOVER2; phi += phi_step)
                    {
                        hairAttenuation(1.55, cosf(theta_i), theta_d, phi, sp.absorption, kfr);
                        float cosphi2 = cosf(phi*0.5f);
                        // [2] eq (11)
                        // Compute average back-scattering attenuation, i.e. total back-scattered radiance
                        // TODO: should be multiplying by cos_theta_i here too..?
                        // {
                        AtRGB f_R = bsdfR(sp.beta_R2, sp.alpha_R, theta_h, cosphi2) * sp.specular1Color * kfr[0] * cos_theta_i;
                        AtRGB f_TT = bsdfTT(sp.beta_TT2, sp.alpha_TT, theta_h, sp.gamma_TT, phi) * sp.transmissionColor * kfr[1] * cos_theta_i;
                        AtRGB f_TRT = (bsdfTRT(sp.beta_TRT2, sp.alpha_TRT, theta_h, cosphi2) * sp.specular2Color
                                    + bsdfg(sp.beta_TRT2, sp.alpha_TRT, theta_h, sp.gamma_g, phi, sp.phi_g) * sp.specular2Color * sp.glintStrength) 
                                    * kfr[2] * cos_theta_i;
                        a_bar_b[idx] += f_R + f_TT + f_TRT;
                        // }
                        alpha_b[idx] += f_R*sp.alpha_R + f_TT*sp.alpha_TT + f_TRT*sp.alpha_TRT;
                        beta_b[idx] += f_R*sp.beta_R2 + f_TT*sp.beta_TT2 + f_TRT*sp.beta_TRT2;
                    }
                }

                alpha_f[idx] /= a_bar_f[idx];
                alpha_b[idx] /= a_bar_b[idx];

                beta_f[idx] /= a_bar_f[idx];
                beta_b[idx] /= a_bar_b[idx];

                a_bar_f[idx] *= 2.0f / AI_PI * theta_r_step * phi_step;
                a_bar_b[idx] *= 2.0f / AI_PI * theta_r_step * phi_step;

                idx++;
            }

            memset(A_b, 0, sizeof(AtRGB)*DS_NUMSTEPS);
            memset(delta_b, 0, sizeof(AtRGB)*DS_NUMSTEPS);
            memset(sigma_b, 0, sizeof(AtRGB)*DS_NUMSTEPS);
            for (int i=0; i < DS_NUMSTEPS; ++i)
            {
                AtRGB af2 = a_bar_f[i]*a_bar_f[i];
                AtRGB ab2 = a_bar_b[i]*a_bar_b[i];
                AtRGB omaf2 = AI_RGB_WHITE - af2;

                // [2] eq (14)
                // Average back-scattered attenuation for up to 3 scattering events
                A_b[i] = (a_bar_b[i]*af2)/omaf2 + (ab2*a_bar_b[i]*af2)/(omaf2*omaf2);

                // [2] eq. (16)
                // Average back-scattering longitudinal shift for up to 3 scattering events
                delta_b[i] = alpha_b[i] * (1.0f - (2.0f*a_bar_b[i]*a_bar_b[i])/(omaf2*omaf2)) 
                            + alpha_f[i] * ((2.0f*omaf2*omaf2) + 4.0f * af2 * ab2) / (omaf2*omaf2*omaf2);

                // [2] eq. (17)
                // Average longitudinal variance for up to 3 scattering events
                AtRGB rtbfbb = sqrt(2.0f*beta_f[i] + beta_b[i]);
                sigma_b[i] = (1.0f + 0.7f*af2) * (a_bar_b[i]*rtbfbb + ab2*a_bar_b[i]*rtbfbb) / (a_bar_b[i] + ab2*a_bar_b[i] * (2.0f*beta_f[i] + 3.0f*beta_b[i]));
            }

            idx = 0;
            memset(N_G_R, 0, sizeof(AtRGB)*DS_NUMSTEPS);
            memset(N_G_TT, 0, sizeof(AtRGB)*DS_NUMSTEPS);
            memset(N_G_TRT, 0, sizeof(AtRGB)*DS_NUMSTEPS);
            for (float theta_d = 0; theta_d < AI_PIOVER2; theta_d+=theta_i_step)
            {
                for (float phi = AI_PIOVER2; phi < AI_PI; phi += phi_step)
                {
                    // [2] eq. (25)
                    // BCSDF due to forward scattering
                    // TODO: including a constant cos_theta_i here just isn't cool.
                    // {
                    hairAttenuation(1.55, 0.1, theta_d, phi, sp.absorption, kfr);

                    N_G_R[idx] += rgb(cosf(phi*0.5f)) * kfr[0];
                    N_G_TT[idx] += rgb(g(sp.gamma_TT, 0.0f, AI_PI-phi)) * kfr[1];
                    N_G_TRT[idx] += rgb(cosf(phi*0.5f) + g(sp.gamma_g, 0.0f, phi - sp.phi_g)) * kfr[2];
                    // }
                }

                N_G_R[idx] *= 2.0 * AI_ONEOVERPI * phi_step;
                N_G_TT[idx] *= 2.0 * AI_ONEOVERPI * phi_step;
                N_G_TRT[idx] *= 2.0 * AI_ONEOVERPI * phi_step;

                ++idx;
            }
            
        }

        float fresnelLookup(float phi)
        {
            int idx = std::min(FRESNEL_SAMPLES-1, static_cast<int>(floorf((phi * AI_ONEOVERPI * FRESNEL_SAMPLES))));
            return kr[idx];
        }

        AtSampler* sampler_diffuse;
        AtSampler* sampler_R;
        AtSampler* sampler_TT;
        AtSampler* sampler_TRT;
        AtSampler* sampler_g;

        float kr[FRESNEL_SAMPLES];
        float invFresnelSamples;

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
    };

    HairBsdf(AtNode* n, AtShaderGlobals* sg, ShaderData* d) :
    node(n), data(d), numBlendHairs(2), density_front(0.7f), density_back(0.7f)
    {
        depth = sg->Rr;

        // Get a local coordinate frame based on the hair fibre direction
        U = AiV3Normalize(sg->dPdv);
        V = AiV3Cross(U, sg->N);
        W = AiV3Cross(V, U);

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
    inline void evaluateParameters(AtShaderGlobals* sg)
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

        beta_R = AiShaderEvalParamFlt(p_specularWidth) * AI_DTOR;
        alpha_R = -AiShaderEvalParamFlt(p_specularShift) * AI_DTOR;

        beta_TT = beta_R * 0.5f;
        alpha_TT = alpha_R * 0.5f;

        beta_TRT = beta_R * 2.0f;
        alpha_TRT = alpha_R * 1.5f;

        beta_R2 = beta_R*beta_R;
        beta_TT2 = beta_TT*beta_TT;
        beta_TRT2 = beta_TRT*beta_TRT;

        gamma_TT = AiShaderEvalParamFlt(p_transmissionRolloff) * AI_DTOR;
        gamma_g = AiShaderEvalParamFlt(p_glintRolloff) * AI_DTOR;
        phi_g = lerp(30.0f*AI_DTOR, 45.0f*AI_DTOR, cn);

        AB(theta_r, alpha_R, beta_R, A_R, B_R);
        AB(theta_r, alpha_TT, beta_TT, A_TT, B_TT);
        AB(theta_r, alpha_TRT, beta_TRT, A_TRT, B_TRT);

        AtRGB hairColor = AiShaderEvalParamRGB(p_hairColor);
        float hairColorDensity = AiShaderEvalParamFlt(p_hairDensity);
        absorption = (AI_RGB_WHITE - min(rgb(0.9999f, 0.9999f, 0.9999f), hairColor)) * hairColorDensity;
        specular1Color = AiShaderEvalParamRGB(p_specular1Color) * AiShaderEvalParamFlt(p_specular1Strength);
        specular2Color = AiShaderEvalParamRGB(p_specular2Color) * AiShaderEvalParamFlt(p_specular2Strength);
        transmissionColor = AiShaderEvalParamRGB(p_transmissionColor) * AiShaderEvalParamFlt(p_transmissionStrength);
        glintStrength = AiShaderEvalParamFlt(p_glintStrength) * cn;

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

    /// Sample according to the diffuse bsdf (just uniform spherical sampling for now)
    inline void sampleDiffuse(float u1, float u2, AtVector& wi)
    {
        wi = uniformSampleSphere(u1, u2);
    }

    /// Sample according to the R lobe
    inline void sampleR(float u1, float u2)
    {
        theta_i = sampleLong(u1, theta_r, alpha_R, beta_R, A_R, B_R);
        phi = 2.0f * asinf(clamp(2.0f*u2 - 1.0f, -1.0f, 1.0f)) + AI_PI;
    }

    /// PDF of the R lobe
    inline float pdfR()
    {
        float t = theta_h-alpha_R;
        float pdf_theta = (1.0f / (2.0f*cos_theta_i*(A_R-B_R))) * (beta_R / (t*t + beta_R2));
        float pdf_phi = cosphi2*0.25f;
        return pdf_theta * pdf_phi;
    }

    ///  Sample according to the TRT lobe
    inline void sampleTRT(float u1, float u2)
    {
        theta_i = sampleLong(u1, theta_r, alpha_TRT, beta_TRT, A_TRT, B_TRT);
        phi = 2.0f * asinf(clamp(2.0f*u2 - 1.0f, -1.0f, 1.0f)) + AI_PI;
    }

    /// PDF of the TRT lobe
    inline float pdfTRT()
    {
        float t = theta_h-alpha_TRT;
        float pdf_theta = (1.0f / (2.0f*cos_theta_i*(A_TRT-B_TRT))) * (beta_TRT / (t*t + beta_TRT2));
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
        float pdf_theta = (1.0f / (2.0f*cos_theta_i*(A_TT-B_TT))) * (beta_TT / (t*t + beta_TT2));
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
        float pdf_theta = (1.0f / (2.0f*cos_theta_i*(A_TRT-B_TRT))) * (beta_TRT / (t*t + beta_TRT2));
        float p = fabs(phi) - phi_g;
        float pdf_phi = (1.0f / (2.0f * (Cg-Dg))) * (gamma_g / (p*p + gamma_g*gamma_g));
        return pdf_theta * pdf_phi;
    }

    inline void sampleUniform(double u1, double u2)
    {
        // Choose a lobe to sample based on which quadrant we are in
        if (1) // (depth < 1)
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
    }

    /// Uniformly sample all the glossy lobes and accumulate the result
    inline void integrateSingleUniform(double u1, double u2)
    {
        sampleUniform(u1, u2);

        // precalculate some stuff
        prepareIndirectSample(wi_ray.dir);
        AtFloat p = invariant / pdfUniform();
        AtScrSample scrs;

        if (p > IMPORTANCE_EPS*0.1)
        {
            // trace our ray
            AiTrace(&wi_ray, &scrs);

            AtRGB kfr[3];
            hairAttenuation(1.55, cos_theta_i, fabsf(theta_d), phi_d, absorption, kfr);

            // calculate result
            result_R_indirect += scrs.color * bsdfR(beta_R2, alpha_R, theta_h, cosphi2) * p * kfr[0];
            result_TT_indirect += scrs.color * bsdfTT(beta_TT2, alpha_TT, theta_h, gamma_TT, phi) * p * kfr[1];
            result_TRT_indirect += scrs.color * bsdfTRT(beta_TRT2, alpha_TRT, theta_h, cosphi2) * p * kfr[2];
            result_TRTg_indirect += scrs.color * bsdfg(beta_TRT2, alpha_TRT, theta_h, gamma_g, phi, phi_g) * p * kfr[2];

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
        return rgb(sqrtf(1.0f - tl*tl)* ONEOVER4PI);
    }

    /// MIS diffuse PDF
    static float HairDiffusePdf(const void* brdf_data, const AtVector* wi)
    {
        return ONEOVER4PI;
    }

    /// Precalculate invariants that will be used across all lobes. This function must be called before any of the bsdf functions during the direct lighting loop
    /// @param wi The incident direction 
    inline void prepareDirectSample(AtVector wi)
    {
        // Get angle measures. See Section 3 in Ou et. al.
        theta_i = (AI_PIOVER2 - sphericalTheta(wi, U));
        cos_theta_i = fabsf(cosf(theta_i));
        phi_i = sphericalPhi(wi, V, W);
        phi_d = phi_r - phi_i;
        if (phi_d < 0) phi_d += AI_PITIMES2;
        phi = phi_d - AI_PI;
        phi = AI_PI - fabsf(phi);

        cosphi2 = fabsf(cosf(phi*0.5f));
        theta_h = (theta_r + theta_i)*0.5f;
        theta_htt = theta_h;
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
        theta_htt = theta_h;
        theta_d = (theta_r - theta_i)*0.5f;
        cos_theta_d = cosf(theta_d);
        inv_cos_theta_d2 = std::max(0.001f, 1.0f/(cos_theta_d*cos_theta_d));
        phi_d = phi;
        if (phi_d < 0) phi_d += AI_PITIMES2;
        phi_i = phi_r - phi;
        //if (phi_i < -AI_PI) phi_i += AI_PITIMES2;
        //if (phi_i > AI_PI) phi_i -= AI_PITIMES2;
        phi = phi_d - AI_PI;
        phi = AI_PI - fabsf(phi);
        cosphi2 = cosf(phi*0.5f);
        sphericalDirection(theta_i, phi_i, V, W, U, wi);
        invariant = cos_theta_i * inv_cos_theta_d2 * AI_ONEOVER2PI;
    }

    /// Integrate the direct illumination for all diffuse and glossy lobes
    inline void integrateDirect(AtShaderGlobals* sg)
    {
        // Tell Arnold we want the full sphere for lighting.
        sg->fhemi = false;
        AiLightsPrepare(sg);
        AtRGB kfr[3];
        while (AiLightsGetSample(sg))
        {
            prepareDirectSample(sg->Ld);
            AtRGB L = sg->Li * sg->we * invariant;
            if (maxh(L) > IMPORTANCE_EPS)
            {
                hairAttenuation(1.55, cos_theta_i, fabsf(theta_d), phi_d, absorption, kfr);

                result_R_direct += L * bsdfR(beta_R2, alpha_R, theta_h, cosphi2) * kfr[0];
                result_TT_direct += L * bsdfTT(beta_TT2, alpha_TT, theta_htt, gamma_TT, phi) * kfr[1];
                result_TRT_direct += L * bsdfTRT(beta_TRT2, alpha_TRT, theta_h, cosphi2) * kfr[2];
                result_TRTg_direct += L * bsdfg(beta_TRT2, alpha_TRT, theta_h, gamma_g, phi, phi_g) * kfr[2];
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
        AiMakeRay(&wi_ray, AI_RAY_GLOSSY, &sg->P, NULL, AI_BIG, sg);
        
        sampit = AiSamplerIterator(data->sampler_R, sg);
        while(AiSamplerGetSample(sampit, samples))
        {
            integrateSingleUniform(samples[0], samples[1]);
        }
        float weight = AiSamplerGetSampleInvCount(sampit);
        result_R_indirect *= specular1Color * weight;
        result_TT_indirect *= transmissionColor * weight;
        result_TRT_indirect *= specular2Color * weight;
        result_TRTg_indirect *= specular2Color * glintStrength * weight;
        result_id2 *= weight;
    }

    /// Integrate the direct illumination with dual scattering for all diffuse and glossy lobes
    inline void integrateDirectDual(AtShaderGlobals* sg)
    {
        if (sg->Rt & AI_RAY_DIFFUSE)
        {
            // if we're in a diffuse ray then don't do the full brdf, just do a kajiya kay diffuse approximation
            // or the result will be noisy as hell (hair caustics!)
            sg->fhemi = false;
            AiLightsPrepare(sg);
            while (AiLightsGetSample(sg))
            {
                result_diffuse_direct += AiEvaluateLightSample(sg, this, HairBsdf::HairDiffuseSample, HairBsdf::HairDiffuseBsdf, HairBsdf::HairDiffusePdf);
            }
            result_diffuse_direct *= specular2Color;
            sg->fhemi = true;
        }
        else
        {
            int als_hairNumIntersections = 0;
            AtRGB als_T_f = AI_RGB_BLACK;
            AtRGB als_sigma_bar_f = AI_RGB_BLACK;
            sg->fhemi = false;
            sg->skip_shadow = true;
            AiLightsPrepare(sg);
            AtRay ray;
            AtScrSample scrs;
            AiMakeRay(&ray, AI_RAY_SHADOW, &(sg->P), NULL, AI_BIG, sg);
            AtRGB kfr[3];
            while (AiLightsGetSample(sg))
            {
                prepareDirectSample(sg->Ld);

                AiStateSetMsgInt("als_hairNumIntersections", 0);
                AiStateSetMsgRGB("als_T_f", AI_RGB_WHITE);
                AiStateSetMsgRGB("als_sigma_bar_f", AI_RGB_BLACK);
                ray.dir = sg->Ld;
                AiTrace(&ray, &scrs);
                AiStateGetMsgInt("als_hairNumIntersections", &als_hairNumIntersections);
                AiStateGetMsgRGB("als_T_f", &als_T_f);
                AiStateGetMsgRGB("als_sigma_bar_f", &als_sigma_bar_f);

                int idx = int((theta_i / AI_PIOVER2 + 0.5f)*(DS_NUMSTEPS-1));
                AtRGB theta_hr = rgb(theta_h);
                AtRGB f_direct_back = 2.0f * data->A_b[idx] * g(data->sigma_b[idx] + als_sigma_bar_f, data->delta_b[idx], theta_hr) * AI_ONEOVERPI * inv_cos_theta_d2;
                AtRGB occlusion = AI_RGB_WHITE - scrs.opacity; 
                float directFraction = 1.0f - std::min(als_hairNumIntersections, numBlendHairs)/float(numBlendHairs);

                if (directFraction > 0.0f)
                {
                    AtRGB kfr[3];
                    hairAttenuation(1.55, cos_theta_i, fabsf(theta_d), phi_d, absorption, kfr);
                    result_Pl_direct += sg->Li * sg->we * occlusion * density_back * f_direct_back * cos_theta_i * directFraction;
                    AtRGB L = sg->Li * sg->we * invariant * occlusion * directFraction;
                    
                    result_R_direct += L * bsdfR(beta_R2, alpha_R, theta_h, cosphi2) * kfr[0];
                    result_TT_direct += L * bsdfTT(beta_TT2, alpha_TT, theta_htt, gamma_TT, phi) * kfr[1];
                    result_TRT_direct += L * bsdfTRT(beta_TRT2, alpha_TRT, theta_h, cosphi2) * kfr[2];
                    result_TRTg_direct += L * bsdfg(beta_TRT2, alpha_TRT, theta_h, gamma_g, phi, phi_g) * kfr[2];
                }
                if (directFraction < 1.0f)
                {
                    AtRGB T_f = als_T_f;
                    AtRGB S_f = g(als_sigma_bar_f, AI_RGB_BLACK, rgb(theta_h)) / (AI_PI * cos_theta_d);

                    AtRGB f_s_scatter = AI_RGB_BLACK;
                    int ngidx = (fabsf(theta_d) / AI_PIOVER2) * DS_NUMSTEPS;
                    if (phi >= AI_PIOVER2 ) // forward scattering directions only
                    {
                        f_s_scatter = g(beta_R2+als_sigma_bar_f, rgb(alpha_R), theta_hr) * data->N_G_R[ngidx]
                                        + g(beta_TT2+als_sigma_bar_f, rgb(alpha_TT), theta_hr) * data->N_G_TT[ngidx]
                                        + g(beta_TRT2+als_sigma_bar_f, rgb(alpha_TRT), theta_hr) * data->N_G_TRT[ngidx];
                        // TODO: without this suppression term we get a hard line in f_s_scatter, which probably means we've made
                        // a mistake somewhere further up the chain...                
                        //f_s_scatter *= phi_c;
                        f_s_scatter *= S_f;
                    }
                    
                    AtRGB f_scatter_back = 2.0f * data->A_b[idx] * g(data->sigma_b[idx], data->delta_b[idx], theta_hr) * AI_ONEOVERPI * inv_cos_theta_d2;

                    
                    AtRGB F_scatter = T_f * density_front * (f_s_scatter + AI_PI*density_back*f_scatter_back);

                    result_Pg_direct += sg->Li * sg->we * occlusion * F_scatter * cos_theta_i * (1.0f-directFraction);

                    result_id1 += f_s_scatter;
                    result_id2 += f_scatter_back;
                    result_id3 += T_f;
                    result_id4 += data->sigma_b[idx];
                    result_id5 += data->delta_b[idx];
                }
            } // END light loop

            result_R_direct *= specular1Color;
            result_TT_direct *= transmissionColor;
            result_TRT_direct *= specular2Color;
            result_TRTg_direct *= specular2Color * glintStrength;

            sg->fhemi = true;
            sg->skip_shadow = false;
        }
    }

    /// Integrate the indirect illumination for dual scattering
    inline void integrateIndirectDual(AtShaderGlobals* sg)
    {
        if (sg->Rt & AI_RAY_DIFFUSE) return;

        AiMakeRay(&wi_ray, AI_RAY_GLOSSY, &sg->P, NULL, AI_BIG, sg);
        AtScrSample scrs;
        sampit = AiSamplerIterator(data->sampler_R, sg);
        int als_hairNumIntersections = 0;
        AtRGB als_T_f = AI_RGB_BLACK;
        AtRGB als_sigma_bar_f = AI_RGB_BLACK;
        while(AiSamplerGetSample(sampit, samples))
        {
            sampleUniform(samples[0], samples[1]);
            prepareIndirectSample(wi_ray.dir);
            AtFloat p = 1.0f / pdfUniform();

            AiStateSetMsgInt("als_hairNumIntersections", 0);
            AiStateSetMsgRGB("als_T_f", AI_RGB_WHITE);
            AiStateSetMsgRGB("als_sigma_bar_f", AI_RGB_BLACK);
            AiTrace(&wi_ray, &scrs);
            AiStateGetMsgInt("als_hairNumIntersections", &als_hairNumIntersections);
            AiStateGetMsgRGB("als_T_f", &als_T_f);
            AiStateGetMsgRGB("als_sigma_bar_f", &als_sigma_bar_f);

            int idx = int((fabs(theta_i) / AI_PIOVER2)*DS_NUMSTEPS);
            AtRGB theta_hr = rgb(theta_h);
            AtRGB f_direct_back = 2.0f * data->A_b[idx] * g(data->sigma_b[idx] + als_sigma_bar_f, data->delta_b[idx], theta_hr) * AI_ONEOVERPI 
                                    * inv_cos_theta_d2;
            float directFraction = 1.0f - std::min(als_hairNumIntersections, numBlendHairs)/float(numBlendHairs);
            if (directFraction > 0.0f)
            {
                AtRGB kfr[3];
                hairAttenuation(1.55, cos_theta_i, fabsf(theta_d), phi_d, absorption, kfr);
                result_Pl_indirect += scrs.color * density_back * f_direct_back * cos_theta_i * directFraction;
                AtRGB L = scrs.color * invariant * p * directFraction;
                result_R_indirect += L * bsdfR(beta_R2, alpha_R, theta_h, cosphi2) * kfr[0];
                result_TT_indirect += L * bsdfTT(beta_TT2, alpha_TT, theta_htt, gamma_TT, phi) * kfr[1];
                result_TRT_indirect += L * bsdfTRT(beta_TRT2, alpha_TRT, theta_h, cosphi2) * kfr[2];
                result_TRTg_indirect += L * bsdfg(beta_TRT2, alpha_TRT, theta_h, gamma_g, phi, phi_g) * kfr[2];

            }
            if (directFraction < 1.0f)
            {
                AtRGB T_f = als_T_f;
                AtRGB S_f = g(als_sigma_bar_f, AI_RGB_BLACK, rgb(theta_h)) / (AI_PI * cos_theta_d);

                AtRGB f_s_scatter = AI_RGB_BLACK;

                if (phi >= AI_PIOVER2) // forward scattering directions only
                {
                    int ngidx = (fabsf(theta_d) / AI_PIOVER2) * DS_NUMSTEPS;
                    f_s_scatter = g(beta_R2+als_sigma_bar_f, rgb(alpha_R), theta_hr) * data->N_G_R[ngidx]
                                    + g(beta_TT2+als_sigma_bar_f, rgb(alpha_TT), theta_hr) * data->N_G_TT[ngidx]
                                    + g(beta_TRT2+als_sigma_bar_f, rgb(alpha_TRT), theta_hr) * data->N_G_TRT[ngidx];

                    // TODO: without this suppression term we get a hard line in f_s_scatter, which probably means we've made
                    // a mistake somewhere further up the chain...                
                   //f_s_scatter *= phi_c;
                    f_s_scatter *= S_f;
                }
                
                AtRGB f_scatter_back = 2.0f * data->A_b[idx] * g(data->sigma_b[idx], data->delta_b[idx], theta_hr) * AI_ONEOVERPI 
                                        * inv_cos_theta_d2;

                float directFraction = 0.0f;
                AtRGB F_scatter = T_f * density_front * (f_s_scatter + AI_PI*density_back*f_scatter_back);

                result_Pg_indirect += scrs.color * F_scatter * cos_theta_i * p * (1.0f-directFraction); 
            }

        }
        float weight = AiSamplerGetSampleInvCount(sampit);
        result_Pg_indirect *= weight;
        result_Pl_indirect *= weight;
        result_R_indirect *= specular1Color * weight;
        result_TT_indirect *= transmissionColor * weight;
        result_TRT_indirect *= specular2Color * weight;
        result_TRTg_indirect *= specular2Color * glintStrength * weight;
    }  

    inline void writeResult(AtShaderGlobals* sg)
    {
        if (sg->Rt & AI_RAY_CAMERA)
        {
            AiAOVSetRGB(sg, "direct_diffuse", result_diffuse_direct);
            AiAOVSetRGB(sg, "indirect_diffuse", result_diffuse_indirect);
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

        sg->out.RGB =   result_diffuse_direct +
                        result_diffuse_indirect +
                        result_R_direct +
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

    // dual-scattering parameters
    int numBlendHairs;
    float density_front;
    float density_back;

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

    AtFloat theta_htt;
    AtRGB absorption;
    float phi_d;
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
    AtNode *options   = AiUniverseGetOptions();
    int diffuse_samples = AiNodeGetInt(options, "GI_diffuse_samples");
    int glossy_samples = std::max(0, AiNodeGetInt(options, "GI_glossy_samples") + params[p_extraSamples].INT);
    

    HairBsdf::ScatteringParameters sp;
    sp.alpha_R = params[p_specularShift].FLT * AI_DTOR;
    sp.beta_R = params[p_specularWidth].FLT * AI_DTOR;
    
    sp.beta_TT = sp.beta_R * 0.5f;
    
    sp.alpha_TT = sp.alpha_R * 0.5f;

    sp.beta_TRT = sp.beta_R * 2.0f;
    sp.alpha_TRT = sp.alpha_R * 1.5f;

    sp.beta_R2 = sp.beta_R*sp.beta_R;
    sp.beta_TT2 = sp.beta_TT*sp.beta_TT;
    sp.beta_TRT2 = sp.beta_TRT*sp.beta_TRT;

    sp.gamma_TT = params[p_transmissionRolloff].FLT * AI_DTOR;
    sp.gamma_g = 35.0f * AI_DTOR;

    sp.specular1Color = params[p_specular1Color].RGB * params[p_specular1Strength].FLT;
    sp.specular2Color = params[p_specular2Color].RGB * params[p_specular2Strength].FLT;
    sp.transmissionColor = params[p_transmissionColor].RGB * params[p_transmissionStrength].FLT;
    sp.glintStrength = params[p_glintStrength].FLT;

    sp.absorption = (AI_RGB_WHITE - min(rgb(0.9999f, 0.9999f, 0.9999f), params[p_hairColor].RGB)) * params[p_hairDensity].FLT;

    data->update(sp, diffuse_samples, glossy_samples);

    AiAOVRegister("direct_specular", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("indirect_specular", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("direct_diffuse", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("indirect_diffuse", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
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





shader_evaluate
{
    
    // Get shader data
    HairBsdf::ShaderData* data = (HairBsdf::ShaderData*)AiNodeGetLocalData(node);

    // Create HairBsdf object 
    HairBsdf hb(node, sg, data);
    // Get parameters
    hb.evaluateParameters(sg);

    int als_hairNumIntersections = 0;
    AtRGB als_T_f = AI_RGB_BLACK;
    AtRGB als_sigma_bar_f = AI_RGB_BLACK;
    AtRGB opacity = AiShaderEvalParamRGB(p_opacity);
    float geo_opacity = 1.0f;
    if (AiUDataGetFlt("geo_opacity", &geo_opacity))
    {
        opacity *= geo_opacity;
    }
    if (sg->Rt & AI_RAY_SHADOW || sg->Rt & AI_RAY_GLOSSY)
    {
        if (AiStateGetMsgInt("als_hairNumIntersections", &als_hairNumIntersections) 
            && AiStateGetMsgRGB("als_T_f", &als_T_f)
            && AiStateGetMsgRGB("als_sigma_bar_f", &als_sigma_bar_f))
        {
            float theta_i = AI_PIOVER2 - sphericalTheta(sg->Rd, hb.U);
            int idx = int((theta_i / AI_PIOVER2 + 0.5f)*(DS_NUMSTEPS-1));

            als_T_f *= hb.data->a_bar_f[idx];
            AiStateSetMsgRGB("als_T_f", als_T_f);

            als_sigma_bar_f += hb.beta_R2 + hb.beta_TRT2 + hb.beta_TT2;
            AiStateSetMsgRGB("als_sigma_bar_f", als_sigma_bar_f);

            als_hairNumIntersections++;
            AiStateSetMsgInt("als_hairNumIntersections", als_hairNumIntersections);
            sg->out_opacity = AI_RGB_BLACK;

            return; // early out
        }
        else
        {
            sg->out_opacity = AI_RGB_WHITE;
        }

        // early-out regardless if we're in a shadow ray, or if opacity is zero
        if (sg->Rt & AI_RAY_SHADOW || AiShaderGlobalsApplyOpacity(sg, opacity)) return; 
    }

    // Do direct illumination
    hb.integrateDirectDual(sg);
    //hb.integrateDirect(sg);

    // Do indirect illumination
    hb.integrateIndirectDual(sg);
    //hb.integrateIndirect(sg);

    // Write shader result
    hb.writeResult(sg);


}


