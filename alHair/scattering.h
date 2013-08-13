// scattering.h
#pragma once

#include "exr.h"

#define PIOVER4 0.7853981633974483f
#define ONEOVER4PI 0.07957747154594767
#define FOURPI 12.566370614359172

#define DS_NUMSTEPS 128
#define NG_NUMSTEPS 32
#define BS_NUMSTEPS 256

float fresnel(float incidenceAngle, float etaPerp, float etaParal, float invert)
{

    float n1, n2;
    float rPerp = 1;
    float rParal = 1;

    etaPerp = std::max(1.0f, etaPerp);
    etaParal = std::max(1.0f, etaParal);

    float angle = fabsf(incidenceAngle);
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

        float b = n2*sqrtf(1.0f-a);
        float c = n1*cosf(angle);
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
    float d = (n1/n2)*sinf(angle);
    d *= d;
    if ( d <= 1 )
    {

        float e = n1*sqrtf(1-d);
        float f = n2*cosf(angle);
        rParal = ( e - f ) / ( e + f );
        rParal *= rParal;
        rParal = std::min( 1.0f, rParal );
    }
    return 0.5 * (rPerp + rParal);
}

#define HAIR_RADIUS 1.0f

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
        if (p != 1)
        {
            if (phi_p > AI_PI) phi_p -= AI_PITIMES2;
            phi_p += p*AI_PI;
        }
        // get roots of polynomial
        float roots[3] = {0,0,0};
        int numRoots = solveCubic(-8.0f*p*c*ONEOVERPI3, 0.0f, 6.0f*p*c*AI_ONEOVERPI - 2.0f, p*AI_PI - phi_p, roots);
        AtRGB Fr = AI_RGB_BLACK;
        kfr[2] = kfr[3] = AI_RGB_BLACK;
        if (p < 2)
        {
            float gamma_i = roots[0];
            if (p==0)
            {
                kfr[0] = rgb(fresnel(gamma_i, n_p, n_pp, false));
            }
            else
            {
                float gamma_t = asinf(sinf(gamma_i)/n_p);
                float l = 2.0f * cosf(gamma_t);
                kfr[1] = (1.0f - fresnel(gamma_i, n_p, n_pp, false)) * (1.0f - fresnel(gamma_t, n_p, n_pp, true)) * exp(-absorption*l);
            }
        }
        else 
        {
            for (int i=0; i < numRoots; ++i)
            {
                float gamma_i = roots[i];
                float gamma_t = asinf(sinf(gamma_i)/n_p);
                float l = 2.0f * cosf(gamma_t);
                float iFr = fresnel(gamma_t, n_p, n_pp, true);
                kfr[2] += (1.0f - fresnel(gamma_i, n_p, n_pp, false)) * (1.0f - iFr) * iFr * exp(-absorption*2*l);
            }
        }
        //kfr[p] = Fr;
    }
}

void hairAttenuation(float ior, float theta_d, float phi, float absorption, float kfr[3])
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
        float Fr = 0.0f;
        kfr[2] = 0.0f;
        if (p < 2)
        {
            float gamma_i = roots[0];
            if (p==0)
            {
                kfr[0] = fresnel(gamma_i, n_p, n_pp, false);
            }
            else
            {
                float gamma_t = asinf(sinf(gamma_i)/n_p);
                float l = 2.0f * cosf(gamma_t);
                kfr[1] = (1.0f - fresnel(gamma_i, n_p, n_pp, false)) * (1.0f - fresnel(gamma_t, n_p, n_pp, true)) * exp(-absorption*l);
            }
        }
        else 
        {
            for (int i=0; i < numRoots; ++i)
            {
                float gamma_i = roots[i];
                float gamma_t = asinf(sinf(gamma_i)/n_p);
                float l = 2.0f * cosf(gamma_t);
                float iFr = fresnel(gamma_t, n_p, n_pp, true);
                kfr[2] += (1.0f - fresnel(gamma_i, n_p, n_pp, false)) * (1.0f - iFr) * iFr * exp(-absorption*2*l);
            }
        }
        //kfr[p] = Fr;
    }
}

void hairAttenuation(float ior, float cos_theta_i, float theta_d, float phi, float phi_h, float aa, AtRGB absorption, AtRGB kfr[3])
{
    hairAttenuation(ior, theta_d, phi, absorption, kfr);
}

void hairAttenuation(float ior, float cos_theta_i, float theta_d, float phi, float phi_h, float aa, float absorption, float kfr[3])
{
    hairAttenuation(ior, theta_d, phi, absorption, kfr);
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
};

struct ScatteringLut
{
    ScatteringLut(float ior, float alpha_R, float alpha_TT, float alpha_TRT, float beta_R2, float beta_TT2, float beta_TRT2, 
                    float gamma_TT, float gamma_g, float phi_g, float absorption)
    {
        float glintStrength = 1.0f; //< TODO: better value for this?

        // precalculate scattering LUTs
        float phi_step = AI_PI / BS_NUMSTEPS;
        float theta_h_step = AI_PI / BS_NUMSTEPS;
        int x, y=0;
        float kfr[3] = {1,1,1};
        for (float phi=0.0f; phi < AI_PI; phi += phi_step)
        {
            x = 0;
            for (float theta_h = -AI_PIOVER2; theta_h < AI_PIOVER2; theta_h += theta_h_step)
            {
                float cosphi2 = cosf(phi*0.5f);
                int idx = y*BS_NUMSTEPS+x;

                b_R[idx] = bsdfR(beta_R2, alpha_R, theta_h, cosphi2);
                b_TT[idx] = bsdfTT(beta_TT2, alpha_TT, theta_h, gamma_TT, phi);
                b_TRT[idx] = bsdfTRT(beta_TRT2, alpha_TRT, theta_h, cosphi2)
                                + bsdfg(beta_TRT2, alpha_TRT, theta_h, gamma_g, phi, phi_g);

                hairAttenuation(ior, theta_h, phi, absorption, kfr);
                k_R[idx] = kfr[0];
                k_TT[idx] = kfr[1];
                k_TRT[idx] = kfr[2];
                x++;
            }
            y++;
        }

        float theta_i_step = AI_PI / DS_NUMSTEPS;
        float theta_r_step = AI_PI/ (DS_NUMSTEPS/2);    //< reduce the inner loop size for speed
        phi_step = AI_PIOVER2 / (DS_NUMSTEPS/2);        //< reduce the inner loop size for speed
        int idx = 0;
        memset(a_bar_f, 0, sizeof(float)*DS_NUMSTEPS);
        memset(a_bar_b, 0, sizeof(float)*DS_NUMSTEPS);
        memset(alpha_f, 0, sizeof(float)*DS_NUMSTEPS);
        memset(alpha_b, 0, sizeof(float)*DS_NUMSTEPS);
        memset(beta_f, 0, sizeof(float)*DS_NUMSTEPS);
        memset(beta_b, 0, sizeof(float)*DS_NUMSTEPS);
        

        for (float theta_i = -AI_PIOVER2; theta_i < AI_PIOVER2; theta_i+=theta_i_step)
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
                    //hairAttenuation(ior, cos_theta_i, theta_d, phi, phi, 1.0f, absorption, kfr);
                    float cosphi2 = cosf(phi*0.5f);
                    
                    // [2] eq (6)
                    // Compute average forward-scattering attenuation, i.e. total forward-scattered radiance
                    // TODO: Should be multiplying by cos_theta_i here too?
                    // {
                    int i_p = (int)(phi/AI_PI*BS_NUMSTEPS);
                    int i_th = (int)((theta_h*AI_ONEOVERPI + 0.5f)*BS_NUMSTEPS);
                    int i_bs = i_p*BS_NUMSTEPS+i_th;
                    int i_td = (int)((theta_d*AI_ONEOVERPI + 0.5f)*BS_NUMSTEPS);
                    int i_k = i_p*BS_NUMSTEPS+i_td;
                    float f_R = b_R[i_bs] * k_R[i_k] * cos_theta_i;
                    float f_TT = b_TT[i_bs] * k_TT[i_k] * cos_theta_i;
                    float f_TRT = b_TRT[i_bs] * k_TRT[i_k] * cos_theta_i;
                    float f = f_R + f_TT + f_TRT;
                    a_bar_f[idx] += f;
                    // }
                    alpha_f[idx] += f_R*alpha_R + f_TT*alpha_TT + f_TRT*alpha_TRT;
                    beta_f[idx] += f_R*beta_R2 + f_TT*beta_TT2 + f_TRT*beta_TRT2;

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
                for (float phi = 0.0f; phi < AI_PIOVER2; phi += phi_step)
                {
                    //hairAttenuation(ior, cosf(theta_i), theta_d, phi, phi, 1.0f, absorption, kfr);
                    float cosphi2 = cosf(phi*0.5f);
                    // [2] eq (11)
                    // Compute average back-scattering attenuation, i.e. total back-scattered radiance
                    // TODO: should be multiplying by cos_theta_i here too..?
                    // {
                    int i_p = (int)(phi/AI_PI*BS_NUMSTEPS);
                    int i_th = (int)((theta_h*AI_ONEOVERPI + 0.5f)*BS_NUMSTEPS);
                    int i_bs = i_p*BS_NUMSTEPS+i_th;
                    int i_td = (int)((theta_d*AI_ONEOVERPI + 0.5f)*BS_NUMSTEPS);
                    int i_k = i_p*BS_NUMSTEPS+i_td;
                    float f_R = b_R[i_bs] * k_R[i_k] * cos_theta_i;
                    float f_TT = b_TT[i_bs] * k_TT[i_k] * cos_theta_i;
                    float f_TRT = b_TRT[i_bs] * k_TRT[i_k] * cos_theta_i;
                    float b = f_R + f_TT + f_TRT;
                    a_bar_b[idx] += b;
                    // }
                    alpha_b[idx] += f_R*alpha_R + f_TT*alpha_TT + f_TRT*alpha_TRT;
                    beta_b[idx] += f_R*beta_R2 + f_TT*beta_TT2 + f_TRT*beta_TRT2;

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
            //AtRGB rtbfbb = sqrt(2.0f*beta_f[i] + beta_b[i]);
            sigma_b[i] = (1.0f + 0.7f*af2) * (a_bar_b[i]*sqrtf(2.0f*beta_f[i] + beta_b[i]) + ab3*sqrtf(2.0f*beta_f[i] + 3.0f*beta_b[i])) / (a_bar_b[i] + ab3 * (2.0f*beta_f[i] + 3.0f*beta_b[i]));
        }

        idx = 0;
        memset(N_G_R, 0, sizeof(float)*NG_NUMSTEPS*NG_NUMSTEPS);
        memset(N_G_TT, 0, sizeof(float)*NG_NUMSTEPS*NG_NUMSTEPS);
        memset(N_G_TRT, 0, sizeof(float)*NG_NUMSTEPS*NG_NUMSTEPS);
        
        float theta_d_step = AI_PI / NG_NUMSTEPS;
        phi_step = AI_PIOVER2 / NG_NUMSTEPS;
        float phi_i_step = AI_PI / NG_NUMSTEPS;

        // phi [pi/2, PI]
        for (int p = 0; p < NG_NUMSTEPS; ++p)        
        {
            float phi = p * phi_step + AI_PIOVER2;
            // theta_d [-pi/2, pi/2]
            for (int t=0; t < NG_NUMSTEPS; ++t)
            {
                float theta_d = t * theta_d_step - AI_PIOVER2;
                int ng_idx = t*NG_NUMSTEPS+p;
                for (float phi_i = -AI_PIOVER2; phi_i < AI_PIOVER2; phi_i += phi_i_step)
                {
                    // [2] eq. (25)
                    // BCSDF due to forward scattering
                    // {
                    //hairAttenuation(ior, theta_d, phi, absorption, kfr);
                    float phi_d = phi - phi_i;
                    if (phi_d > AI_PI) phi_d = AI_PI - (phi_d-AI_PI);
                    int i_p = (int)(phi_d/AI_PI*BS_NUMSTEPS);
                    int i_td = (int)(fabsf(theta_d)/AI_PIOVER2*BS_NUMSTEPS);
                    int i_k = i_p*BS_NUMSTEPS+i_td;

                    float cosphi2 = cosf(phi_d*0.5f);
                    N_G_R[ng_idx] += cosphi2 * k_R[i_k];
                    N_G_TT[ng_idx] += g(gamma_TT, 0.0f, AI_PI-phi_d)* k_TT[i_k];
                    N_G_TRT[ng_idx] += cosphi2 + g(gamma_g, 0.0f, phi_d - phi_g) * k_TRT[i_k];
                }

                N_G_R[ng_idx] *= AI_ONEOVERPI * phi_i_step;
                N_G_TT[ng_idx] *= AI_ONEOVERPI * phi_i_step;
                N_G_TRT[ng_idx] *= AI_ONEOVERPI * phi_i_step;
            }
        }
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
        // See if we have a dslut struct already created for this colour, if not create it.
        int idx = (int)(floorf(a/A_MAX * DS_MASTER_LUT_SZ));
        if (!_luts[idx])
        {
            AiCritSecEnter(&_cs);
            
                if (!_luts[idx])
                {
                    float t0 = AiMsgUtilGetElapsedTime();
                    _luts[idx] = new ScatteringLut(sp.ior, sp.alpha_R, sp.alpha_TT, sp.alpha_TRT, sp.beta_R2, sp.beta_TT2,
                                                        sp.beta_TRT2, sp.gamma_TT, sp.gamma_g, sp.phi_g, float(idx)/float(DS_MASTER_LUT_SZ) * A_MAX); 
                    std::cerr << "generated lut at " << idx << " with absorption " << float(idx)/float(DS_MASTER_LUT_SZ) * A_MAX << std::endl;
                    _lutgen_time += AiMsgUtilGetElapsedTime() - t0;
                }

            
            
            AiCritSecLeave(&_cs);
        }
        return _luts[idx];
    }

    AtRGB forward_attenuation(const ScatteringParams& sp, float theta)
    {        
        AtRGB result;
        int idx = int((theta/AI_PI + 0.5f)*(DS_NUMSTEPS-1));
        ScatteringLut* lr = getLut(sp.absorption.r, sp);
        ScatteringLut* lg = getLut(sp.absorption.g, sp);
        ScatteringLut* lb = getLut(sp.absorption.b, sp);
        result.r = lr->a_bar_f[idx];
        result.g = lg->a_bar_f[idx];
        result.b = lb->a_bar_f[idx];

        return result;
    }

    AtRGB direct_back_scatter(const ScatteringParams& sp, float theta, float theta_h, const AtRGB& als_sigma_bar_f)
    {
        AtRGB result;
        ScatteringLut* lr = getLut(sp.absorption.r, sp);
        ScatteringLut* lg = getLut(sp.absorption.g, sp);
        ScatteringLut* lb = getLut(sp.absorption.b, sp);

        int idx = int((theta/AI_PI + 0.5f)*(DS_NUMSTEPS-1));

        result.r = 2.0f * lr->A_b[idx] * g(lr->sigma_b[idx] + als_sigma_bar_f.r, lr->delta_b[idx], theta_h) * AI_ONEOVERPI;
        result.g = 2.0f * lg->A_b[idx] * g(lg->sigma_b[idx] + als_sigma_bar_f.g, lg->delta_b[idx], theta_h) * AI_ONEOVERPI;
        result.b = 2.0f * lb->A_b[idx] * g(lb->sigma_b[idx] + als_sigma_bar_f.b, lb->delta_b[idx], theta_h) * AI_ONEOVERPI;

        return result;
    }

    AtRGB back_scatter(const ScatteringParams& sp, float theta, float theta_h)
    {
    	AtRGB result;
        ScatteringLut* lr = getLut(sp.absorption.r, sp);
        ScatteringLut* lg = getLut(sp.absorption.g, sp);
        ScatteringLut* lb = getLut(sp.absorption.b, sp);

        int idx = int((theta/AI_PI + 0.5f)*(DS_NUMSTEPS-1));

        result.r = 2.0f * lr->A_b[idx] * g(lr->sigma_b[idx], lr->delta_b[idx], theta_h) * AI_ONEOVERPI;
        result.g = 2.0f * lg->A_b[idx] * g(lg->sigma_b[idx], lg->delta_b[idx], theta_h) * AI_ONEOVERPI;
        result.b = 2.0f * lb->A_b[idx] * g(lb->sigma_b[idx], lb->delta_b[idx], theta_h) * AI_ONEOVERPI;

        return result;
    }

    AtRGB forward_scatter(const ScatteringParams& sp, float theta_h, float theta_d, float phi, const AtRGB& als_sigma_bar_f)
    {
    	int ng_x = int((phi-AI_PIOVER2)/AI_PIOVER2) * (NG_NUMSTEPS-1);
        int ng_y = (theta_d * AI_ONEOVERPI + 0.5f) * (NG_NUMSTEPS-1);
        int ng_idx = ng_y*NG_NUMSTEPS+ng_x;

        ScatteringLut* lr = getLut(sp.absorption.r, sp);
        ScatteringLut* lg = getLut(sp.absorption.g, sp);
        ScatteringLut* lb = getLut(sp.absorption.b, sp);

        AtRGB result;
    	result.r = 	g(sp.beta_R2+als_sigma_bar_f.r, sp.alpha_R, theta_h) * lr->N_G_R[ng_idx]
                  + g(sp.beta_TT2+als_sigma_bar_f.r, sp.alpha_TT, theta_h) * lr->N_G_TT[ng_idx]
                  + g(sp.beta_TRT2+als_sigma_bar_f.r, sp.alpha_TRT, theta_h) * lr->N_G_TRT[ng_idx];

        result.g = 	g(sp.beta_R2+als_sigma_bar_f.g, sp.alpha_R, theta_h) * lg->N_G_R[ng_idx]
                  + g(sp.beta_TT2+als_sigma_bar_f.g, sp.alpha_TT, theta_h) * lg->N_G_TT[ng_idx]
                  + g(sp.beta_TRT2+als_sigma_bar_f.g, sp.alpha_TRT, theta_h) * lg->N_G_TRT[ng_idx];

        result.b = 	g(sp.beta_R2+als_sigma_bar_f.b, sp.alpha_R, theta_h) * lb->N_G_R[ng_idx]
                  + g(sp.beta_TT2+als_sigma_bar_f.b, sp.alpha_TT, theta_h) * lb->N_G_TT[ng_idx]
                  + g(sp.beta_TRT2+als_sigma_bar_f.b, sp.alpha_TRT, theta_h) * lb->N_G_TRT[ng_idx];

        return result;
    }



    ~DualScattering()
    {
        AiCritSecClose(&_cs);
        for (int i=0; i < DS_MASTER_LUT_SZ; ++i)
        {
            delete _luts[i];
        }
        AiMsgInfo("[alHair] total lutgen time: %.2fs", _lutgen_time/1000.0f);
    }

    ScatteringLut* _luts[DS_MASTER_LUT_SZ];
    AtCritSec _cs;
    float _lutgen_time;
};