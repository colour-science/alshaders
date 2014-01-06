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
    p_scale,
	p_offset,
    p_P
};

node_parameters
{
	AiParameterENUM("space", 0, triplanarSpaceNames);
    AiParameterENUM("tiling", 0, triplanarTilingNames);
    AiParameterSTR("texture", "");
    AiParameterFLT("blendSoftness", 0.1);
    AiParameterFLT("scale", 1.0f);
    AiParameterFLT("offset", 0.0f);
	AiParameterPnt("P", 0.0f, 0.0f, 0.0f);
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

AtPoint getProjectionP(const AtNode* node, const AtShaderGlobals *sg, const AtPoint &Pin, int space, float scale, float offset){
    AtPoint P;
    if (AiNodeIsLinked(node, "P"))
    {
        P = Pin;
    }
    else
    {
        switch (space)
        {
        case NS_OBJECT:
            P = sg->Po;
            break;
        case NS_PREF:
            if (!AiUDataGetPnt("Pref", &P))
                P = sg->Po;
            break;
        default:
            P = sg->P;
            break;
        }
    }
    return P / scale + offset;
}

void computeBlendWeights(const AtShaderGlobals *sg, int space, float blendSoftness, float *weights){
    AtPoint N =
    weights[0] = fabsf(sg->N.x);
    weights[1] = fabsf(sg->N.y);
    weights[2] = fabsf(sg->N.z);
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

inline AtRGB tileRegular(const AtPoint &P, float *weights, AtShaderGlobals *sg,
AtTextureHandle *handle, AtTextureParams *params){
    AtRGBA textureResult[3];
    bool textureAccessX = false;
    bool textureAccessY = false;
    bool textureAccessZ = false;

    // lookup X
    sg->u = P.z + 123.94;
    sg->v = P.y + 87.22;
    textureResult[0] = AiTextureHandleAccess(sg, handle, params, &textureAccessX);

    // lookup Y
    sg->u = P.z + 9.2;
    sg->v = P.x + 74.1;
    textureResult[1] = AiTextureHandleAccess(sg, handle, params, &textureAccessY);

    // lookup Z
    sg->u = P.x - 12.55;
    sg->v = P.y + 6;
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

inline bool lookupCellNoise(float u, float v,  AtRGBA *textureResult, AtShaderGlobals *sg,
                            AtTextureHandle *handle, AtTextureParams *params){
    AtPoint P;
    P.x = u;
    P.y = v;
    P.z = 0.f;

    // run cellnoise
    float f;
    AtVector delta;
    AtUInt32 id;
    AiCellular(P, 1, 1, 1.92, 1, &f, &delta, &id);

    // pick random direction for orientation
    AtVector orientVectorX;
    double phi = random(id) * M_PI * 2;
    orientVectorX.x = cosf(phi);
    orientVectorX.y = sinf(phi);
    orientVectorX.z = 0.f;

    AtVector orientVectorZ;
    orientVectorZ.x = 0.f;
    orientVectorZ.y = 0.f;
    orientVectorZ.z = 1.f;

    AtVector orientVectorY = AiV3Cross(orientVectorX, orientVectorZ);

    AiV3RotateToFrame(delta, orientVectorX, orientVectorY, orientVectorZ);

    // find new uv coordinates
    sg->u = delta.x;
    sg->v = delta.y;

    bool success = false;
    *textureResult = AiTextureHandleAccess(sg, handle, params, &success);

    return success;
}

inline AtRGB tileCellnoise(const AtPoint &P, float *weights, AtShaderGlobals *sg,
                           AtTextureHandle *handle, AtTextureParams *params){
    AtRGBA textureResult[3];
    bool textureAccessX = lookupCellNoise(P.y,
                                          P.z,
                                          &textureResult[0],
                                          sg,
                                          handle,
                                          params);
    bool textureAccessY = lookupCellNoise(P.x,
                                          P.z,
                                          &textureResult[1],
                                          sg,
                                          handle,
                                          params);
    bool textureAccessZ = lookupCellNoise(P.y,
                                          P.x,
                                          &textureResult[2],
                                          sg,
                                          handle,
                                          params);


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
    float scale = AiShaderEvalParamFlt(p_scale);
	float offset = AiShaderEvalParamFlt(p_offset);
    AtPoint Pin = AiShaderEvalParamPnt(p_P);

        // get local data
    ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);

        // set up P and blend weights
    AtPoint P = getProjectionP(node, sg, Pin, space, scale, offset);
    float weights[3];
    computeBlendWeights(sg, space, blendSoftness, weights);

        // compute texture values
    AtRGB result = AI_RGB_RED;
    switch(tiling){
        case TM_CELLNOISE:
            result = tileCellnoise(P, weights, sg, data->texturehandle, data->textureparams);
            break;
        case TM_REGULAR:
            result = tileRegular(P, weights, sg, data->texturehandle, data->textureparams);
            break;
        default:
            // TODO: We should never end up here. Log the error to inform the shader writer.
            result = AI_RGB_BLUE;
            break;
    }
    sg->out.RGB = result;
}


