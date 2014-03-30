#pragma once
#include <ai.h>
#include "alUtil.h"

void sampleHemisphereCosine(AtVector *vec, float r1, float r2)
{
    float stheta = sqrtf(float(r1));
    float phi = float(AI_PITIMES2 * r2);
    vec->x = stheta * cosf(phi);
    vec->y = stheta * sinf(phi);
    vec->z = sqrtf(1.0f - float(r1));
}

struct DiffuseBRDF
{
    void *brdf_data;
    const AtShaderGlobals *sg;

    AtColor result;
    bool active;

    float albedo;
    AtColor strength;
    AtColor color;

};


void diffuseBRDFinit(DiffuseBRDF *brdf, const AtShaderGlobals *sg){
    brdf->sg = sg;
    brdf->active = true;
    brdf->result = AI_RGB_BLACK;
    brdf->albedo = 0.f;
    brdf->strength = AI_RGB_BLUE;
}

void diffuseBRDFcomputeAlbedo(DiffuseBRDF *brdf){
    brdf->albedo = 1.f;
}

void diffuseBRDFinitMIS(DiffuseBRDF *brdf){
    brdf->brdf_data = AiOrenNayarMISCreateData(brdf->sg, 0.f);
}

//struct SpecularBRDF
//{
//    void *brdf_data;

//    AtColor result;
//    bool active;

//    float albedo;
//    AtColor strength;
//    float metallic;
//    float roughness;
//    AtColor eta;

//};

//struct TransmisionBRDF
//{
//    void *brdf_data;

//    bool active;
//    AtColor result;

//    float albedo;
//    AtColor strength;
//    float roughness;
//    AtColor eta;

//};

