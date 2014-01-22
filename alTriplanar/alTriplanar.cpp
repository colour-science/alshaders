#include "Remap.h"
#include <ai.h>

AI_SHADER_NODE_EXPORT_METHODS(alTriplanar)

enum TriplanarSpaceEnum
{
	NS_WORLD = 0,
	NS_OBJECT,
	NS_PREF,
};

static const char* triplanarSpaceNames[] =
{
	"world",
	"object",
	"Pref",
	NULL
};

enum TriplanarTiling
{
    TM_REGULAR = 0,
    TM_CELLNOISE
};

static const char* triplanarTilingNames[] =
{
    "regular",
    "cellnoise",
    NULL
};

enum alTriplanarParams
{
    p_space,
    p_tiling,
    p_texture,
    p_blendSoftness,
    p_cellSoftness,
    p_scalex,
    p_scaley,
    p_scalez,
    p_offsetx,
    p_offsety,
    p_offsetz,
    p_rotx,
    p_roty,
    p_rotz,
    p_rotjitterx,
    p_rotjittery,
    p_rotjitterz
};

node_parameters
{
	AiParameterENUM("space", 0, triplanarSpaceNames);
    AiParameterENUM("tiling", 0, triplanarTilingNames);
    AiParameterSTR("texture", "");
    AiParameterFLT("blendSoftness", 0.1);
    AiParameterFLT("cellSoftness", 0.1);
    AiParameterFLT("scalex", 1.0f);
    AiParameterFLT("scaley", 1.0f);
    AiParameterFLT("scalez", 1.0f);
    AiParameterFLT("offsetx", 0.0f);
    AiParameterFLT("offsety", 0.0f);
    AiParameterFLT("offsetz", 0.0f);
    AiParameterFLT("rotx", 0.0f);
    AiParameterFLT("roty", 0.0f);
    AiParameterFLT("rotz", 0.0f);
    AiParameterFLT("rotjitterx", 1.0f);
    AiParameterFLT("rotjittery", 1.0f);
    AiParameterFLT("rotjitterz", 1.0f);
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alTriplanar;
   node->output_type = AI_TYPE_RGB;
   node->name        = "alTriplanar";
   node->node_type   = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return true;
}

struct ShaderData
{
    AtTextureHandle* texturehandle;
    AtTextureParams *textureparams;
};

node_initialize
{
    ShaderData *data = new ShaderData;
    const char *texname = params[p_texture].STR;
    data->texturehandle = AiTextureHandleCreate(texname);
    data->textureparams = new AtTextureParams;
    AiTextureParamsSetDefaults(data->textureparams);
    AiNodeSetLocalData(node, data);
}

node_finish
{
    ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);
    AiTextureHandleDestroy(data->texturehandle);
    delete data->textureparams;
    delete data;
}

node_update
{

}

AtPoint getProjectionGeometry(const AtNode* node, const AtShaderGlobals *sg, int space, AtPoint *P, AtVector *N){
    switch (space)
    {
    case NS_OBJECT:
        *P = sg->Po;
        *N = AiShaderGlobalsTransformNormal(sg, sg->Ng, AI_WORLD_TO_OBJECT);
        break;
    case NS_PREF:
        if (!(AiUDataGetPnt("Pref", P) && AiUDataGetVec("Nref", N))){
            // TODO: Output warning about not finding the correct data.
            *P = sg->Po;
            *N = AiShaderGlobalsTransformNormal(sg, sg->Ng, AI_WORLD_TO_OBJECT);
        }
        break;
    default:
        *P = sg->P;
        *N = sg->Ng;
        break;
    }
}

void computeBlendWeights(const AtVector N, int space, float blendSoftness, float *weights){
    weights[0] = fabsf(N.x);
    weights[1] = fabsf(N.y);
    weights[2] = fabsf(N.z);
    float weightsum = 0.f;
    for(int i=0; i<3; ++i){
        weights[i] = weights[i] - (1.f-blendSoftness)/2.f;
        weights[i] = MAX(weights[i], 0.00f);
        weightsum += weights[i];
    }
    if(weightsum){
        for(int i=0; i<3; ++i){
            weights[i] /= weightsum;
        }
    }
}

inline void rotateUVs(AtPoint &P, float degrees){
    AtVector orientVectorX;
    const double d2r = 1. / 360. * M_PI * 2;
    double phi = d2r * degrees;
    orientVectorX.x = cosf(phi);
    orientVectorX.y = sinf(phi);
    orientVectorX.z = 0.f;

    AtVector orientVectorZ;
    orientVectorZ.x = 0.f;
    orientVectorZ.y = 0.f;
    orientVectorZ.z = 1.f;

    AtVector orientVectorY = AiV3Cross(orientVectorX, orientVectorZ);

    AiV3RotateToFrame(P, orientVectorX, orientVectorY, orientVectorZ);
}

inline AtRGB tileRegular(const AtPoint &P, const AtPoint &scale, const AtPoint &offset,
                         float *weights, const AtPoint &rot, AtShaderGlobals *sg,
                         AtTextureHandle *handle, AtTextureParams *params){
    AtRGBA textureResult[3];
    bool textureAccessX = false;
    bool textureAccessY = false;
    bool textureAccessZ = false;

    // lookup X
    AtPoint ProjP;
    ProjP.x = (P.z + 123.94 + offset.x) * scale.x ;
    ProjP.y = (P.y + 87.22 + offset.x) * scale.x;
    ProjP.z = 0.;
    rotateUVs(ProjP, rot.x);

    sg->u = ProjP.x;
    sg->v = ProjP.y;
    textureResult[0] = AiTextureHandleAccess(sg, handle, params, &textureAccessX);

    // lookup Y
    ProjP.x = (P.x + 74.1 + offset.y) * scale.y;
    ProjP.y = (P.z + 9.2 + offset.y) * scale.y;
    ProjP.z = 0.;
    rotateUVs(ProjP, rot.y);

    sg->u = ProjP.x;
    sg->v = ProjP.y;
    textureResult[1] = AiTextureHandleAccess(sg, handle, params, &textureAccessY);

    // lookup Z
    ProjP.x = (P.x + 123.94 + offset.z) * scale.z;
    ProjP.y = (P.y + 87.22 + offset.z) * scale.z;
    ProjP.z = 0.;
    rotateUVs(ProjP, rot.z);

    sg->u = ProjP.x;
    sg->v = ProjP.y;
    textureResult[2] = AiTextureHandleAccess(sg, handle, params, &textureAccessZ);

    if(textureAccessX && textureAccessY && textureAccessZ){
        AtRGB result = AI_RGB_BLACK;
        result += textureResult[0].rgb() * weights[0];
        result += textureResult[1].rgb() * weights[1];
        result += textureResult[2].rgb() * weights[2];
        return result;
    }
    else {
        // Something went wrong during lookup.
        // TODO: Log the error
        return AiColorCreate(1.f, 0.0f, 0.f);
    }
}

inline bool lookupCellNoise(float u, float v, float cellSoftness,
                            float rot, float rotjitter,
                            AtShaderGlobals *sg, AtTextureHandle *handle,
                            AtTextureParams *params, AtRGBA *textureResult){
    AtPoint P;
    P.x = u;
    P.y = v;
    P.z = 0.f;

    int samples = (cellSoftness == 0) ? 1 : 3;
    // run cellnoise
    float weights[samples];
    float f[samples];
    AtVector delta[samples];
    AtUInt32 id[samples];
    AiCellular(P, samples, 1, 1.92, 1, f, delta, id);

    if(samples == 1){
        weights[0] = 1.f;
    } else {
        // find closest cell
        float closestDistance = 100000.f;
        float distances[samples];
        for(int i=0; i<samples; ++i){
            distances[i] = AiV3Length(delta[i]);
            closestDistance = MIN(distances[i], closestDistance);
        }

        float weightsum = 0.f;
        for(int i=0; i<samples; ++i){
            float diff = distances[i] - closestDistance;
            weights[i] = cellSoftness - diff;
            weights[i] = MAX(0.f, weights[i]);
            weightsum += weights[i];
        }
        for(int i=0; i<samples; ++i){
            weights[i] /= weightsum;
        }
    }

    bool success = false;
    *textureResult = AI_RGBA_BLACK;
    for(int i=0; i<samples; ++i){
        if(weights[i] > 0.f){
            // pick direction for orientation
            AtVector orientVectorX;
            double jitter = (random(id[i])-0.5) * rotjitter;
            double phi = modulo(rot/360. + jitter, 1.f) * M_PI * 2.;
            orientVectorX.x = cosf(phi);
            orientVectorX.y = sinf(phi);
            orientVectorX.z = 0.f;

            AtVector orientVectorZ;
            orientVectorZ.x = 0.f;
            orientVectorZ.y = 0.f;
            orientVectorZ.z = 1.f;

            AtVector orientVectorY = AiV3Cross(orientVectorX, orientVectorZ);

            AiV3RotateToFrame(delta[i], orientVectorX, orientVectorY, orientVectorZ);

            // find new uv coordinates, set the center of the cell to be 0.5/0.5;
            sg->u = delta[i].x * 0.75 - 0.5;
            sg->v = delta[i].y * 0.75 - 0.5;

            // texture lookup
            bool currentSuccess = false;
            *textureResult += AiTextureHandleAccess(sg, handle, params, &success) * weights[i];
            success += currentSuccess;
        }
    }

    return success;
}

inline AtRGB tileCellnoise(const AtPoint &P, const AtPoint &scale, const AtPoint &offset,
                           float *weights, float cellSoftness,
                           const AtPoint &rot, const AtPoint &rotjitter, AtShaderGlobals *sg,
                           AtTextureHandle *handle, AtTextureParams *params){
    AtRGBA textureResult[3];
    bool textureAccessX = lookupCellNoise((P.y + offset.x) * scale.x,
                                          (P.z + offset.x) * scale.x,
                                          cellSoftness,
                                          rot.x,
                                          rotjitter.x,
                                          sg,
                                          handle,
                                          params,
                                          &textureResult[0]
                                          );
    bool textureAccessY = lookupCellNoise((P.x + offset.y) * scale.y,
                                          (P.z + offset.y) * scale.y,
                                          cellSoftness,
                                          rot.y,
                                          rotjitter.y,
                                          sg,
                                          handle,
                                          params,
                                          &textureResult[1]
                                          );
    bool textureAccessZ = lookupCellNoise((P.y + offset.z) * scale.z,
                                          (P.x + offset.z) * scale.z,
                                          cellSoftness,
                                          rot.z,
                                          rotjitter.z,
                                          sg,
                                          handle,
                                          params,
                                          &textureResult[2]
                                          );


    if(textureAccessX && textureAccessY && textureAccessZ){
        AtRGB result = AI_RGB_BLACK;
        result += textureResult[0].rgb() * weights[0];
        result += textureResult[1].rgb() * weights[1];
        result += textureResult[2].rgb() * weights[2];
        return result;
    }
    else {
        // Something went wrong during lookup.
        // TODO: Log the error
        return AiColorCreate(1.f, 0.0f, 0.f);
    }
}

shader_evaluate
{
        // get shader parameters
	int space = AiShaderEvalParamInt(p_space);

    int tiling = AiShaderEvalParamInt(p_tiling);

    float blendSoftness = AiShaderEvalParamFlt(p_blendSoftness);
    blendSoftness = CLAMP(blendSoftness, 0.f, 1.f);

    float cellSoftness = AiShaderEvalParamFlt(p_cellSoftness);
    cellSoftness = CLAMP(cellSoftness, 0.f, 1.f);

    AtPoint scale;
    scale.x = 1.f/AiShaderEvalParamFlt(p_scalex);
    scale.y = 1.f/AiShaderEvalParamFlt(p_scaley);
    scale.z = 1.f/AiShaderEvalParamFlt(p_scalez);

    AtPoint offset;
    offset.x = AiShaderEvalParamFlt(p_offsetx);
    offset.y = AiShaderEvalParamFlt(p_offsety);
    offset.z = AiShaderEvalParamFlt(p_offsetz);

    AtPoint rot;
    rot.x = AiShaderEvalParamFlt(p_rotx);
    rot.y = AiShaderEvalParamFlt(p_roty);
    rot.z = AiShaderEvalParamFlt(p_rotz);

    AtPoint rotjitter;
    rotjitter.x = AiShaderEvalParamFlt(p_rotjitterx);
    rotjitter.y = AiShaderEvalParamFlt(p_rotjittery);
    rotjitter.z = AiShaderEvalParamFlt(p_rotjitterz);

        // get local data
    ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);

        // set up P and blend weights
    AtPoint P;
    AtPoint N;
    getProjectionGeometry(node, sg, space, &P, &N);
    float weights[3];
    computeBlendWeights(N, space, blendSoftness, weights);

        // compute texture values
    AtRGB result = AI_RGB_RED;
    switch(tiling){
        case TM_CELLNOISE:
            result = tileCellnoise(P, scale, offset, weights, cellSoftness, rot, rotjitter, sg, data->texturehandle, data->textureparams);
            break;
        case TM_REGULAR:
            result = tileRegular(P, scale, offset, weights, rot, sg, data->texturehandle, data->textureparams);
            break;
        default:
            // TODO: We should never end up here. Log the error to inform the shader writer.
            result = AI_RGB_BLUE;
            break;
    }
    sg->out.RGB = result;
}


