#include "alUtil.h"

AtColor beckmannMicrofacetTransmission(AtShaderGlobals* sg, const AtVector& Z, const AtVector& X, const AtVector& Y,
                                        const AtVector& wo, AtSampler* sampler, float roughness, float eta,
                                        AtRGB sigma_s, AtRGB sigma_a, float g,
                                        float ssScale, bool inScattering, AtRGB& ss_result)
{
    double samples[2];
    float n1, n2;
    float kt;
    AtRay wi_ray;
    AtScrSample sample;
    AtVector wi, R;
    bool inside;
    AtRGB result = AI_RGB_BLACK;
    AtRGB sigma_t = sigma_s + sigma_a;
    AtRGB sigma_s_prime = sigma_s*(1.0f-g);
    AtRGB sigma_t_prime = (sigma_s_prime + sigma_a);
    AtRGB mfp = AI_RGB_WHITE / sigma_t_prime;

    AtSamplerIterator* sampit = AiSamplerIterator(sampler, sg);

    AiMakeRay(&wi_ray, AI_RAY_REFRACTED, &sg->P, NULL, AI_BIG, sg);

    while (AiSamplerGetSample(sampit, samples))
    {
        // generate a microfacet normal, m
        // eq. 35,36
        float alpha2 = roughness*roughness;
        float tanThetaM = sqrtf(-alpha2 * logf(1.0f - float(samples[0])));
        float cosThetaM = 1.0f / sqrtf(1.0f + tanThetaM * tanThetaM);
        float sinThetaM = cosThetaM * tanThetaM;
        float phiM = 2.0f * float(AI_PI) * float(samples[1]);
        AtVector m = (cosf(phiM) * sinThetaM) * X +
                     (sinf(phiM) * sinThetaM) * Y +
                                   cosThetaM  * Z;


        // get the refracted direction given m
        kt = 1.0f - fresnel(eta, m, wo, R, wi, inside);

        if (kt > IMPORTANCE_EPS) // if the final contribution is actually going to matter to the result
        {
            // eq. 33

            float cosThetaM2 = cosThetaM * cosThetaM;
            float tanThetaM2 = tanThetaM * tanThetaM;
            float cosThetaM4 = cosThetaM2 * cosThetaM2;
            float D = expf(-tanThetaM2 / alpha2) / (float(AI_PI) * alpha2 *  cosThetaM4);
            // eq. 24
            float pm = D * cosThetaM;
            // eval BRDF*cosNI
            float cosNI = AiV3Dot(Z, wi); // N.wi
            float cosNO = AiV3Dot(Z, wo);
            // eq. 26, 27: now calculate G1(i,m) and G1(o,m)
            float ao = 1 / (roughness * sqrtf((1.0f - cosNO * cosNO) / (cosNO * cosNO)));
            float ai = 1 / (roughness * sqrtf((1.0f - cosNI * cosNI) / (cosNI * cosNI)));
            float G1o = ao < 1.6f ? (3.535f * ao + 2.181f * ao * ao) / (1 + 2.276f * ao + 2.577f * ao * ao) : 1.0f;
            float G1i = ai < 1.6f ? (3.535f * ai + 2.181f * ai * ai) / (1 + 2.276f * ai + 2.577f * ai * ai) : 1.0f;
            float G = G1o * G1i;
            // eq. 21
            float cosHI = AiV3Dot(m, wi); // m.wi
            float cosHO = AiV3Dot(m, wo); // m.wo
            float Ht2 = eta * cosHI + cosHO;
            Ht2 *= Ht2;
            float brdf = (fabsf(cosHI * cosHO) * (eta * eta) * (G * D)) / fabsf(cosNO * Ht2);
            // eq. 38 and eq. 17
            float pdf = pm * (eta * eta) * fabsf(cosHI) / Ht2;

            wi_ray.dir = wi;
            AiTrace(&wi_ray, &sample);
            AtRGB transmittance = AI_RGB_WHITE;
            if (maxh(sigma_t) > 0.0f && !inside)
            {
                transmittance.r = expf(-sample.z * sigma_t.r);
                transmittance.g = expf(-sample.z * sigma_t.g);
                transmittance.b = expf(-sample.z * sigma_t.b);
            }
            result += sample.color * brdf/pdf * transmittance;

            // single scattering
            if (ssScale > IMPORTANCE_EPS && maxh(sigma_s_prime) > 0.0f && !inside && inScattering)
            {
                AtVector N = sg->N;
                sg->N = m;
                ss_result += AiSSSTraceSingleScatter(sg,bssrdfbrdf(sigma_s_prime/sigma_t_prime),mfp,g,eta) * ssScale * kt;
                sg->N = N;
            }


        }
        else if (AiV3IsZero(wi)) // total internal reflection
        {
            wi_ray.dir = R;
            AiTrace(&wi_ray, &sample);
            AtRGB transmittance = AI_RGB_WHITE;
            if (maxh(sigma_t) > 0.0f && !inside)
            {
                transmittance.r = expf(-sample.z * sigma_t.r);
                transmittance.g = expf(-sample.z * sigma_t.g);
                transmittance.b = expf(-sample.z * sigma_t.b);
            }
            result += sample.color * transmittance;
        }
    }

    result *= AiSamplerGetSampleInvCount(sampit);
    ss_result *= AiSamplerGetSampleInvCount(sampit);

    return result;
}
