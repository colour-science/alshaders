#pragma once
#include <ai.h>

struct ShaderData;

void microfacetRefraction(AtShaderGlobals *sg,
							const ShaderData *data,
							const AtFloat rf_IOR,
							AtFloat rf_roughness,
							AtColor &refr_result);


inline void make_orthonormals (const AtVector &N, AtVector &a, AtVector &b)
{
    if (N.x != N.y || N.x != N.z)
        {
                a.x = N.z-N.y;
                a.y = N.x-N.z;
                a.z = N.y-N.x;
                // (1,1,1) x N
        }
    else
        {
        a.x = N.z-N.y;
                a.y = N.x+N.z;
                a.z = -N.y-N.x;
                // (-1,1,1) x N
        }
    AiV3Normalize(a, a);
    AiV3Cross(b, N, a);
}

inline AtVector beckmann_sample (const AtFloat &roughness, const AtVector &N, const AtDouble rnd[2], AtFloat &cosThetaM, AtFloat &tanThetaM, AtFloat &alpha2)
{
        AtVector X, Y, Z = N;
        make_orthonormals(Z, X, Y);
        // generate a random microfacet normal m
        // eq. 35,36:
        // we take advantage of cos(atan(x)) == 1/sqrt(1+x^2)
        //                  and sin(atan(x)) == x/sqrt(1+x^2)

        alpha2 = roughness * roughness;
        tanThetaM = sqrtf(-alpha2 * logf(1.0f - (AtFloat)rnd[0]));
        cosThetaM = 1.0f / sqrtf(1.0f + tanThetaM * tanThetaM);
        AtFloat sinThetaM = cosThetaM * tanThetaM;
        AtFloat phiM = 2 * AtFloat(AI_PI) * (AtFloat)rnd[1];
        return  X *(cosf(phiM) * sinThetaM) +
                        Y * (sinf(phiM) * sinThetaM) +
                        Z * cosThetaM;
}

inline AtFloat beckmann_distribution (AtFloat &cosThetaM, AtFloat &tanThetaM, AtFloat &alpha2)
{
    AtFloat cosThetaM2 = cosThetaM * cosThetaM;
    AtFloat tanThetaM2 = tanThetaM * tanThetaM;
    AtFloat cosThetaM4 = cosThetaM2 * cosThetaM2;
    return expf(-tanThetaM2 / alpha2) / ((AtFloat)AI_PI * alpha2 * cosThetaM4);
}

inline AtFloat beckmann_shadowing (const AtFloat &roughness, AtFloat &cosNO, AtFloat &cosNI)
{
        AtFloat ao = 1 / (roughness * sqrtf((1.0f - cosNO * cosNO) / (cosNO * cosNO)));
        AtFloat ai = 1 / (roughness * sqrtf((1.0f - cosNI * cosNI) / (cosNI * cosNI)));
        AtFloat G1o = ao < 1.6f ? (3.535f * ao + 2.181f * ao * ao) / (1 + 2.276f * ao + 2.577f * ao * ao) : 1.0f;
        AtFloat G1i = ai < 1.6f ? (3.535f * ai + 2.181f * ai * ai) / (1 + 2.276f * ai + 2.577f * ai * ai) : 1.0f;
        return G1o * G1i;
}

inline AtFloat remap(AtFloat value, const AtFloat min, const AtFloat max, const AtFloat target_min, const AtFloat target_max)
{
        AtFloat start_rng = max - min;
        AtFloat end_rng = target_max - target_min;
        return (((value-min)*end_rng)/start_rng)+target_min;
}
