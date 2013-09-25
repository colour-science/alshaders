// scattering.h
#pragma once

#include <sstream>

#define PIOVER4 0.7853981633974483f
#define ONEOVER4PI 0.07957747154594767
#define FOURPI 12.566370614359172

//#define LUT_NAN_CHECK
#define LUT_DEBUG

#define DS_NUMSTEPS 64
#define NG_NUMSTEPS 32
#define BS_NUMSTEPS 256

void hairAttenuation(float ior, float theta_d, float phi, const AtRGB& absorption, AtRGB kfr[3])
{
    float x = clamp(1.0f - cosf(theta_d), 0.f, 1.f);
    float x3 = powf(x, 5.0f);
    float x5 = powf(x, 7.0f);
    float y = clamp(phi / AI_PI, 0.f, 1.f);
    float y5 = powf(y, 5.0f);
    float cf_front = lerp(0.02f, 0.5f, x3);
    float cf_back = lerp(0.04f, 0.9f, x5);
    kfr[0] = rgb(lerp(cf_front, cf_back, y5));
    kfr[1] = lerp(0.0f, 0.9f, 1.0f - x5) * fast_exp(-absorption * 2.0f);
    kfr[2] = lerp(0.1f, 0.2f, x5) * fast_exp(-absorption * 2.0f * lerp(1.5f, 2.0f, x3));
}

void hairAttenuation(float ior, float theta_d, float phi, float absorption, float kfr[3])
{
    float x = 1.0f - cosf(theta_d);
    float x3 = powf(x, 5.0f);
    float x5 = powf(x, 7.0f);
    float y = phi / AI_PI;
    float y5 = powf(y, 5.0f);
    float cf_front = lerp(0.02f, 0.5f, x3);
    float cf_back = lerp(0.04f, 0.9f, x5);
    kfr[0] = lerp(cf_front, cf_back, y5);
    kfr[1] = lerp(0.0f, 0.9f, 1.0f - x5) * fast_exp(-absorption * 2.0f);
    kfr[2] = lerp(0.1f, 0.2f, x5) * fast_exp(-absorption * 2.0f * lerp(1.5f, 2.0f, x3));
}

/// Normalized gaussian with offset
inline float g(float beta, float alpha, float theta_h)
{
    float n = theta_h-alpha;
    return fast_exp(-(n*n)/(2.0f*beta))/sqrtf(2.0f*AI_PI*beta);
}

inline AtRGB g(AtRGB beta, AtRGB alpha, AtRGB theta_h)
{
    return AiColorCreate(
        g(beta.r, alpha.r, theta_h.r),
        g(beta.g, alpha.g, theta_h.g),
        g(beta.b, alpha.b, theta_h.b)
    );
}

inline float g(float x, float v)
{
    return fast_exp(-SQR(x) / (2.0f*v)) / sqrtf(2.0f*AI_PI*v);
}

inline AtRGB g(AtRGB x, AtRGB v)
{
    return rgb(
        g(x.r, v.r),
        g(x.g, v.g),
        g(x.b, v.b)
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

struct ScatteringParams
{
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
    AtRGB absorption;
    AtRGB dabsorption;
    float shape;
};

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
        inv_cos_theta_d2 = 1.0f/std::max(0.001f, (cos_theta_d*cos_theta_d));
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

struct ScatteringLut
{
    ScatteringLut(float ior, float alpha_R, float alpha_TT, float alpha_TRT, float beta_R2, float beta_TT2, float beta_TRT2, 
                    float gamma_TT, float gamma_g, float phi_g, float absorption, float shape)
    {
        // precalculate scattering LUTs
        float phi_step = AI_PI / (BS_NUMSTEPS-1);
        float theta_h_step = AI_PI / (BS_NUMSTEPS-1);
        int x, y=0;
        float kfr[3] = {1,1,1};
        for (float p = 0; p < BS_NUMSTEPS; ++p)
        {
            float phi = p * phi_step;
            x = 0;
            for (float t=0; t < BS_NUMSTEPS; ++t)
            {
                float theta_h = t * theta_h_step - AI_PIOVER2;
                float cosphi2 = cosf(phi*0.5f);
                int idx = y*BS_NUMSTEPS+x;

                b_R[idx] = bsdfR(beta_R2, alpha_R, theta_h, cosphi2);
                b_TT[idx] = bsdfTT(beta_TT2, alpha_TT, theta_h, gamma_TT, phi);
                b_TRT[idx] = bsdfTRT(beta_TRT2, alpha_TRT, theta_h, cosphi2)
                                + bsdfg(beta_TRT2, alpha_TRT, theta_h, gamma_g, phi, 35.0f);

                hairAttenuation(ior, theta_h, phi, absorption, kfr);
                k_R[idx] = kfr[0];
                k_TT[idx] = kfr[1];
                k_TRT[idx] = kfr[2];
                x++;
            }
            y++;
        }

        float theta_i_step = AI_PI / (DS_NUMSTEPS-1);
        float theta_r_step = AI_PI / (DS_NUMSTEPS-1);    //< reduce the inner loop size for speed
        phi_step = AI_PIOVER2 / (DS_NUMSTEPS-1);        //< reduce the inner loop size for speed
        int idx = 0;
        memset(a_bar_f, 0, sizeof(float)*DS_NUMSTEPS);
        memset(a_bar_b, 0, sizeof(float)*DS_NUMSTEPS);
        memset(alpha_f, 0, sizeof(float)*DS_NUMSTEPS);
        memset(alpha_b, 0, sizeof(float)*DS_NUMSTEPS);
        memset(beta_f, 0, sizeof(float)*DS_NUMSTEPS);
        memset(beta_b, 0, sizeof(float)*DS_NUMSTEPS);
        
        for (float ti = 0; ti < DS_NUMSTEPS; ++ti)
        {
            float theta_i = ti * theta_i_step - AI_PIOVER2;
            float cos_theta_i = cosf(theta_i);
            for (float tr = 0; tr < DS_NUMSTEPS; ++tr)
            {
                float theta_r = tr * theta_r_step - AI_PIOVER2;
                float theta_h = (theta_i + theta_r) * 0.5f;
                float theta_d = (theta_i - theta_r) * 0.5f;
                int i_th = (int)((theta_h*AI_ONEOVERPI + 0.5f)*(BS_NUMSTEPS-1));
                int i_td = (int)((theta_d*AI_ONEOVERPI + 0.5f)*(BS_NUMSTEPS-1));

                // integrate over half the domain each time and double after as it's symmetrical
                // forward scattering
                for (float phi = AI_PIOVER2; phi <= AI_PI; phi += phi_step)
                {
                    float cosphi2 = cosf(phi*0.5f);
                    
                    // [2] eq (6)
                    // Compute average forward-scattering attenuation, i.e. total forward-scattered radiance
                    // TODO: Should be multiplying by cos_theta_i here too?
                    // {
                    int i_p = (int)(phi*AI_ONEOVERPI*(BS_NUMSTEPS-1));
                    int i_bs = i_p*BS_NUMSTEPS+i_th;
                    int i_k = i_p*BS_NUMSTEPS+i_td;
                    
                    float f_R = b_R[i_bs] * k_R[i_k] * cos_theta_i;
                    float f_TT = b_TT[i_bs] * k_TT[i_k] * cos_theta_i;
                    float f_TRT = b_TRT[i_bs] * k_TRT[i_k] * cos_theta_i;
 
                    float f = f_R + f_TT + f_TRT;
                    a_bar_f[idx] += f;
                    // }
                    alpha_f[idx] += (f_R*alpha_R + f_TT*alpha_TT + f_TRT*alpha_TRT);
                    beta_f[idx] += (f_R*beta_R2 + f_TT*beta_TT2 + f_TRT*beta_TRT2);

#ifdef LUT_NAN_CHECK
                    if (!AiIsFinite(f))
                    {
                        std::cerr << VAR(f) << std::endl;
                        std::cerr << VAR(f_R) << std::endl;;
                        std::cerr << VAR(f_TT) << std::endl;;
                        std::cerr << VAR(f_TRT) << std::endl;;
                        std::cerr << VAR(theta_i) << std::endl;;
                        std::cerr << VAR(theta_r) << std::endl;;
                        std::cerr << VAR(theta_h) << std::endl;;
                        std::cerr << VAR(theta_d) << std::endl;;
                        std::cerr << VAR(phi) << std::endl;
                        std::cerr << VAR(i_bs) << std::endl;
                        std::cerr << VAR(i_k) << std::endl;
                    }
#endif
                }

                // backward scattering
                for (float phi = 0.0f; phi <= AI_PIOVER2; phi += phi_step)
                {
                    float cosphi2 = cosf(phi*0.5f);
                    // [2] eq (11)
                    // Compute average back-scattering attenuation, i.e. total back-scattered radiance
                    // TODO: should be multiplying by cos_theta_i here too..?
                    // {
                    int i_p = (int)(phi/AI_PI*(BS_NUMSTEPS-1));
                    int i_bs = i_p*BS_NUMSTEPS+i_th;
                    int i_k = i_p*BS_NUMSTEPS+i_td;
                    
                    float f_R = b_R[i_bs] * k_R[i_k] * cos_theta_i;;
                    float f_TT = b_TT[i_bs] * k_TT[i_k] * cos_theta_i;;
                    float f_TRT = b_TRT[i_bs] * k_TRT[i_k] * cos_theta_i;;
                    
                    float b = f_R + f_TT + f_TRT;
                    a_bar_b[idx] += b;
                    // }
                    alpha_b[idx] += (f_R*alpha_R + f_TT*alpha_TT + f_TRT*alpha_TRT);
                    beta_b[idx] += (f_R*beta_R2 + f_TT*beta_TT2 + f_TRT*beta_TRT2);

#ifdef LUT_NAN_CHECK
                    if (!AiIsFinite(b))
                    {
                        std::cerr << VAR(b) << std::endl;
                        std::cerr << VAR(f_R) << std::endl;;
                        std::cerr << VAR(f_TT) << std::endl;;
                        std::cerr << VAR(f_TRT) << std::endl;;
                        std::cerr << VAR(theta_i) << std::endl;;
                        std::cerr << VAR(theta_r) << std::endl;;
                        std::cerr << VAR(theta_h) << std::endl;;
                        std::cerr << VAR(theta_d) << std::endl;;
                        std::cerr << VAR(phi) << std::endl;;
                        std::cerr << VAR(i_bs) << std::endl;
                        std::cerr << VAR(i_k) << std::endl;
                    }
#endif
                }
            }

            alpha_f[idx] /= a_bar_f[idx];
            alpha_b[idx] /= a_bar_b[idx];

            beta_f[idx] /= a_bar_f[idx];
            beta_b[idx] /= a_bar_b[idx];

            a_bar_f[idx] *= 2.0f * AI_ONEOVERPI * theta_r_step * phi_step;
            a_bar_b[idx] *= 2.0f * AI_ONEOVERPI * theta_r_step * phi_step;

            idx++;
        }

        memset(A_b, 0, sizeof(float)*DS_NUMSTEPS);
        memset(delta_b, 0, sizeof(float)*DS_NUMSTEPS);
        memset(sigma_b, 0, sizeof(float)*DS_NUMSTEPS);
        for (int i=0; i < DS_NUMSTEPS; ++i)
        {
            //a_bar_f[i] = a_bar_f[i], 0.99f);
            //a_bar_b[i] = std::min(a_bar_b[i], 0.99f);
            float af2 = a_bar_f[i]*a_bar_f[i];
            float ab2 = a_bar_b[i]*a_bar_b[i];
            float ab3 = a_bar_b[i] * ab2;
            float omaf2 = 1.0f - af2;

            // [2] eq (14)
            // Average back-scattered attenuation for up to 3 scattering events
            A_b[i] = (a_bar_b[i]*af2)/omaf2 + (ab3*af2)/(omaf2*omaf2*omaf2);

            // [2] eq. (16)
            // Average back-scattering longitudinal shift for up to 3 scattering events
            delta_b[i] = alpha_b[i] * (1.0f - (2.0f*ab2)/(omaf2*omaf2)) 
                        + alpha_f[i] * ((2.0f*omaf2*omaf2) + 4.0f * af2 * ab2) / (omaf2*omaf2*omaf2);

            // [2] eq. (17)
            // Average longitudinal variance for up to 3 scattering events
            sigma_b[i] = (1.0f + 0.7f*af2) * (a_bar_b[i]*sqrtf(2.0f*beta_f[i] + beta_b[i]) + ab3*sqrtf(2.0f*beta_f[i] + 3.0f*beta_b[i])) / (a_bar_b[i] + ab3 * (2.0f*beta_f[i] + 3.0f*beta_b[i]));
            sigma_b[i] = std::min(sigma_b[i], 1.0f);
#ifdef LUT_NAN_CHECK
            if (!AiIsFinite(A_b[i]) || !AiIsFinite(delta_b[i]) || !AiIsFinite(sigma_b[i]))
            {
                std::cerr << VAR(A_b) << std::endl;;
                std::cerr << VAR(delta_b) << std::endl;;
                std::cerr << VAR(sigma_b) << std::endl;;
                std::cerr << VAR(af2) << std::endl;
                std::cerr << VAR(ab2) << std::endl;;
                std::cerr << VAR(ab3) << std::endl;;
                std::cerr << VAR(omaf2) << std::endl;;
            }
#endif
        }

        idx = 0;
        memset(N_G_R, 0, sizeof(float)*NG_NUMSTEPS*NG_NUMSTEPS);
        memset(N_G_TT, 0, sizeof(float)*NG_NUMSTEPS*NG_NUMSTEPS);
        memset(N_G_TRT, 0, sizeof(float)*NG_NUMSTEPS*NG_NUMSTEPS);
        
        float theta_d_step = AI_PI / (NG_NUMSTEPS-1);
        phi_step = AI_PIOVER2 / (NG_NUMSTEPS-1);
        float phi_i_step = AI_PI / (NG_NUMSTEPS-1);

        // phi [pi/2, PI]
        for (int p = 0; p < NG_NUMSTEPS; ++p)        
        {
            float phi = p * phi_step + AI_PIOVER2;
            // theta_d [-pi/2, pi/2]
            for (int t=0; t < NG_NUMSTEPS; ++t)
            {
                float theta_d = t * theta_d_step - AI_PIOVER2;
                float x5 = powf(1.0f-cosf(theta_d), 5);
                int ng_idx = t*NG_NUMSTEPS+p;
                for (float phi_i = -AI_PIOVER2; phi_i <= AI_PIOVER2; phi_i += phi_i_step)
                {
                    // [2] eq. (25)
                    // BCSDF due to forward scattering
                    // {
                    
                    float phi_d = phi - phi_i;
                    if (phi_d > AI_PI) phi_d = AI_PI - (phi_d-AI_PI);
                    int i_p = (int)(phi_d/AI_PI*(BS_NUMSTEPS-1));
                    int i_td = (int)(theta_d*AI_ONEOVERPI+0.5f*(BS_NUMSTEPS-1));
                    int i_k = i_p*(BS_NUMSTEPS)+i_td;

                    float cosphi2 = cosf(phi_d*0.5f);
                    
                    N_G_R[ng_idx] += cosphi2 * k_R[i_k];
                    N_G_TT[ng_idx] += g(gamma_TT, 0.0f, AI_PI-phi_d)* k_TT[i_k];
                    N_G_TRT[ng_idx] += (cosphi2 + g(gamma_g, 0.0f, phi_d - phi_g)) * k_TRT[i_k];
                    
                    // }
                }

                // TODO: 1/pi factor here?
                N_G_R[ng_idx] *= phi_i_step;
                N_G_TT[ng_idx] *= phi_i_step;
                N_G_TRT[ng_idx] *= phi_i_step;
            }
        }

#ifdef LUT_DEBUG
    #include "exr.h"
    std::stringstream ss;

    ss << "/tmp/a_bar_f_" << absorption << ".exr";
    writeThickFloatEXR(ss.str().c_str(), a_bar_f, DS_NUMSTEPS);
    ss.str(std::string());

    ss << "/tmp/a_bar_b_" << absorption << ".exr";
    writeThickFloatEXR(ss.str().c_str(), a_bar_b, DS_NUMSTEPS);
    ss.str(std::string());

    ss << "/tmp/A_b_" << absorption << ".exr";
    writeThickFloatEXR(ss.str().c_str(), A_b, DS_NUMSTEPS);
    ss.str(std::string());

    ss << "/tmp/alpha_f_" << absorption << ".exr";
    writeThickFloatEXR(ss.str().c_str(), alpha_f, DS_NUMSTEPS);
    ss.str(std::string());

    ss << "/tmp/alpha_b_" << absorption << ".exr";
    writeThickFloatEXR(ss.str().c_str(), alpha_b, DS_NUMSTEPS);
    ss.str(std::string());

    ss << "/tmp/beta_f_" << absorption << ".exr";
    writeThickFloatEXR(ss.str().c_str(), beta_f, DS_NUMSTEPS);
    ss.str(std::string());

    ss << "/tmp/beta_b_" << absorption << ".exr";
    writeThickFloatEXR(ss.str().c_str(), beta_b, DS_NUMSTEPS);
    ss.str(std::string());

    ss << "/tmp/sigma_b_" << absorption << ".exr";
    writeThickFloatEXR(ss.str().c_str(), sigma_b, DS_NUMSTEPS);
    ss.str(std::string());

    ss << "/tmp/delta_b_" << absorption << ".exr";
    writeThickFloatEXR(ss.str().c_str(), delta_b, DS_NUMSTEPS);
    ss.str(std::string());

    ss << "/tmp/N_G_R_" << absorption << ".exr";
    writeFloatEXR(ss.str().c_str(), N_G_R, NG_NUMSTEPS, NG_NUMSTEPS);
    ss.str(std::string());

    ss << "/tmp/N_G_TT_" << absorption << ".exr";
    writeFloatEXR(ss.str().c_str(), N_G_TT, NG_NUMSTEPS, NG_NUMSTEPS);
    ss.str(std::string());

    ss << "/tmp/N_G_TRT_" << absorption << ".exr";
    writeFloatEXR(ss.str().c_str(), N_G_TRT, NG_NUMSTEPS, NG_NUMSTEPS);
    ss.str(std::string());

    ss << "/tmp/b_R_" << absorption << ".exr";
    writeFloatEXR(ss.str().c_str(), b_R, BS_NUMSTEPS, BS_NUMSTEPS);
    ss.str(std::string());

    ss << "/tmp/b_TT_" << absorption << ".exr";
    writeFloatEXR(ss.str().c_str(), b_TT, BS_NUMSTEPS, BS_NUMSTEPS);
    ss.str(std::string());

    ss << "/tmp/b_TRT_" << absorption << ".exr";
    writeFloatEXR(ss.str().c_str(), b_TRT, BS_NUMSTEPS, BS_NUMSTEPS);
    ss.str(std::string());

    ss << "/tmp/k_R_" << absorption << ".exr";
    writeFloatEXR(ss.str().c_str(), k_R, BS_NUMSTEPS, BS_NUMSTEPS);
    ss.str(std::string());

    ss << "/tmp/k_TT_" << absorption << ".exr";
    writeFloatEXR(ss.str().c_str(), k_TT, BS_NUMSTEPS, BS_NUMSTEPS);
    ss.str(std::string());

    ss << "/tmp/k_TRT_" << absorption << ".exr";
    writeFloatEXR(ss.str().c_str(), k_TRT, BS_NUMSTEPS, BS_NUMSTEPS);
    ss.str(std::string());

#endif        
    }

    inline float get_a_bar_f(float theta)
    {
        float t = (theta/AI_PI + 0.5f)*(DS_NUMSTEPS-1);
        int idx = int(t);
        int idx1 = std::min(idx+1, DS_NUMSTEPS-1);
        t -= idx;
        return a_bar_f[idx]*(1.0f-t) + a_bar_f[idx1]*t;
    }

    inline float get_A_b(float theta)
    {
        float t = (theta/AI_PI + 0.5f)*(DS_NUMSTEPS-1);
        int idx = int(t);
        int idx1 = std::min(idx+1, DS_NUMSTEPS-1);
        t -= idx;
        return A_b[idx]*(1.0f-t) + A_b[idx1]*t;
    }

    inline float get_delta_b(float theta)
    {
        float t = (theta/AI_PI + 0.5f)*(DS_NUMSTEPS-1);
        int idx = int(t);
        int idx1 = std::min(idx+1, DS_NUMSTEPS-1);
        t -= idx;
        return delta_b[idx]*(1.0f-t) + delta_b[idx1]*t;
    }

    inline float get_sigma_b(float theta)
    {
        float t = (theta/AI_PI + 0.5f)*(DS_NUMSTEPS-1);
        int idx = int(t);
        int idx1 = std::min(idx+1, DS_NUMSTEPS-1);
        t -= idx;
        return sigma_b[idx]*(1.0f-t) + sigma_b[idx1]*t;
    }

    float a_bar_f[DS_NUMSTEPS];
    float a_bar_b[DS_NUMSTEPS];
    float A_b[DS_NUMSTEPS];
    float alpha_f[DS_NUMSTEPS];
    float alpha_b[DS_NUMSTEPS];
    float beta_f[DS_NUMSTEPS];
    float beta_b[DS_NUMSTEPS];
    float sigma_b[DS_NUMSTEPS];
    float delta_b[DS_NUMSTEPS];
    float N_G_R[NG_NUMSTEPS*NG_NUMSTEPS];
    float N_G_TT[NG_NUMSTEPS*NG_NUMSTEPS];
    float N_G_TRT[NG_NUMSTEPS*NG_NUMSTEPS];
    float b_R[BS_NUMSTEPS*BS_NUMSTEPS];
    float b_TT[BS_NUMSTEPS*BS_NUMSTEPS];
    float b_TRT[BS_NUMSTEPS*BS_NUMSTEPS];
    float k_R[BS_NUMSTEPS*BS_NUMSTEPS];
    float k_TT[BS_NUMSTEPS*BS_NUMSTEPS];
    float k_TRT[BS_NUMSTEPS*BS_NUMSTEPS];
};

#define DS_MASTER_LUT_SZ 100
#define A_MAX 1.0f
struct DualScattering
{
    DualScattering()
    : _cachemisses(0.0), _cachelookups(0.0)
    {
        AiCritSecInit(&_cs);
        for (int i=0; i < DS_MASTER_LUT_SZ; ++i)
        {
            _luts[i] = NULL;
        }
    }

    void reset()
    {
        if (_luts[0])
        {
            AiCritSecEnter(&_cs);
                if (_luts[0])
                {
                    for (int i=0; i < DS_MASTER_LUT_SZ; ++i)
                    {
                        delete _luts[i];
                        _luts[i] = NULL;
                    }
                }
            AiCritSecLeave(&_cs);
        }
    }

    inline ScatteringLut* getLut(float a, const ScatteringParams& sp)
    {
        _cachelookups++;
        // See if we have a dslut struct already created for this colour, if not create it.
        int idx = (int)(floorf(a/A_MAX * DS_MASTER_LUT_SZ));
        if (!_luts[idx])
        {
            AiCritSecEnter(&_cs);
                _cachemisses++;
                if (!_luts[idx])
                {
                    float t0 = AiMsgUtilGetElapsedTime();
                    _luts[idx] = new ScatteringLut(sp.ior, sp.alpha_R, sp.alpha_TT, sp.alpha_TRT, sp.beta_R2, sp.beta_TT2,
                                                        sp.beta_TRT2, sp.gamma_TT, sp.gamma_g, sp.phi_g, float(idx)/float(DS_MASTER_LUT_SZ) * A_MAX, sp.shape); 
                    std::cerr << "generated lut at " << idx << " with absorption " << float(idx)/float(DS_MASTER_LUT_SZ) * A_MAX << std::endl;
                    _lutgen_time += AiMsgUtilGetElapsedTime() - t0;
                }
            
            AiCritSecLeave(&_cs);
        }
        return _luts[idx];
    }

    inline AtRGB forward_attenuation(const ScatteringParams& sp, float theta)
    {        
        AtRGB result;
        
        ScatteringLut* lr = getLut(sp.dabsorption.r, sp);
        ScatteringLut* lg = getLut(sp.dabsorption.g, sp);
        ScatteringLut* lb = getLut(sp.dabsorption.b, sp);
        result.r = lr->get_a_bar_f(theta);
        result.g = lg->get_a_bar_f(theta);
        result.b = lb->get_a_bar_f(theta);

        return result;
    }

    inline AtRGB f_back_direct(const ScatteringParams& sp, const SctGeo& geo)
    {
        ScatteringLut* lr = getLut(sp.dabsorption.r, sp);
        ScatteringLut* lg = getLut(sp.dabsorption.g, sp);
        ScatteringLut* lb = getLut(sp.dabsorption.b, sp);

        AtRGB result;
        float c = 2.0f * AI_ONEOVERPI * geo.inv_cos_theta_d2;
        result.r = lr->get_A_b(geo.theta_d) * g(geo.theta_h - lr->get_delta_b(geo.theta_d), SQR(lr->get_sigma_b(geo.theta_d))) * c;
        result.g = lg->get_A_b(geo.theta_d) * g(geo.theta_h - lg->get_delta_b(geo.theta_d), SQR(lg->get_sigma_b(geo.theta_d))) * c;
        result.b = lb->get_A_b(geo.theta_d) * g(geo.theta_h - lb->get_delta_b(geo.theta_d), SQR(lb->get_sigma_b(geo.theta_d))) * c;

        return result;
    }

    inline AtRGB f_back_scatter(const ScatteringParams& sp, const SctGeo& geo, const AtRGB& sigma_f2)
    {
        ScatteringLut* lr = getLut(sp.dabsorption.r, sp);
        ScatteringLut* lg = getLut(sp.dabsorption.g, sp);
        ScatteringLut* lb = getLut(sp.dabsorption.b, sp);

        AtRGB result;
        float c = 2.0f * AI_ONEOVERPI * geo.inv_cos_theta_d2;
        result.r = lr->get_A_b(geo.theta_d) * g(geo.theta_h - lr->get_delta_b(geo.theta_d), SQR(lr->get_sigma_b(geo.theta_d)) + sigma_f2.r) * c;
        result.g = lg->get_A_b(geo.theta_d) * g(geo.theta_h - lg->get_delta_b(geo.theta_d), SQR(lg->get_sigma_b(geo.theta_d)) + sigma_f2.g) * c;
        result.b = lb->get_A_b(geo.theta_d) * g(geo.theta_h - lb->get_delta_b(geo.theta_d), SQR(lb->get_sigma_b(geo.theta_d)) + sigma_f2.b) * c;

        return result;
    }

    inline AtRGB f_s_scatter(const ScatteringParams& sp, const SctGeo& geo, const AtRGB& sigma_f2)
    {
        int ng_x = int((geo.phi-AI_PIOVER2)/AI_PIOVER2) * (NG_NUMSTEPS-1);
        int ng_y = (geo.theta_d * AI_ONEOVERPI + 0.5f) * (NG_NUMSTEPS-1);
        int ng_idx = ng_y*NG_NUMSTEPS+ng_x;

        ScatteringLut* lr = getLut(sp.dabsorption.r, sp);
        ScatteringLut* lg = getLut(sp.dabsorption.g, sp);
        ScatteringLut* lb = getLut(sp.dabsorption.b, sp);

        AtRGB result = AI_RGB_BLACK;
        result.r =  g(geo.theta_h - sp.alpha_R, sp.beta_R2 + sigma_f2.r) * lr->N_G_R[ng_idx]
                    + g(geo.theta_h - sp.alpha_TT, sp.beta_TT2 + sigma_f2.r) * lr->N_G_TT[ng_idx]
                    + g(geo.theta_h - sp.alpha_TRT, sp.beta_TRT2 + sigma_f2.r) * lr->N_G_TRT[ng_idx];

        result.g =  g(geo.theta_h - sp.alpha_R, sp.beta_R2 + sigma_f2.g) * lr->N_G_R[ng_idx]
                    + g(geo.theta_h - sp.alpha_TT, sp.beta_TT2 + sigma_f2.g) * lr->N_G_TT[ng_idx]
                    + g(geo.theta_h - sp.alpha_TRT, sp.beta_TRT2 + sigma_f2.g) * lr->N_G_TRT[ng_idx];

        result.b =  g(geo.theta_h - sp.alpha_R, sp.beta_R2 + sigma_f2.b) * lr->N_G_R[ng_idx]
                    + g(geo.theta_h - sp.alpha_TT, sp.beta_TT2 + sigma_f2.b) * lr->N_G_TT[ng_idx]
                    + g(geo.theta_h - sp.alpha_TRT, sp.beta_TRT2 + sigma_f2.b) * lr->N_G_TRT[ng_idx];

        return result * geo.inv_cos_theta_d2;
    }

    ~DualScattering()
    {
        AiCritSecClose(&_cs);
        for (int i=0; i < DS_MASTER_LUT_SZ; ++i)
        {
            delete _luts[i];
        }
        if (_cachelookups > 0.0)
        {
            AiMsgInfo("[alHair] total lutgen time: %.2fs", _lutgen_time/1000.0f);
            AiMsgInfo("[alHair] cache hit rate: %.2f", (1.0-_cachemisses/_cachelookups)*100.0);
        }
    }

    ScatteringLut* _luts[DS_MASTER_LUT_SZ];
    AtCritSec _cs;
    float _lutgen_time;
    float _norm;
    double _cachemisses, _cachelookups;
};