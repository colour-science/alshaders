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
    AiParameterFLT("cellSoftness", 0.1);
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
    weights[0] = fabsf(sg->Ng.x);
    weights[1] = fabsf(sg->Ng.y);
    weights[2] = fabsf(sg->Ng.z);
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

inline bool lookupCellNoise(float u, float v, float cellSoftness, AtShaderGlobals *sg,
                            AtTextureHandle *handle, AtTextureParams *params, AtRGBA *textureResult){
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
            // pick random direction for orientation
            AtVector orientVectorX;
            double phi = random(id[i]) * M_PI * 2;
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

inline AtRGB tileCellnoise(const AtPoint &P, float *weights, float cellSoftness, AtShaderGlobals *sg,
                           AtTextureHandle *handle, AtTextureParams *params){
    AtRGBA textureResult[3];
    bool textureAccessX = lookupCellNoise(P.y,
                                          P.z,
                                          cellSoftness,
                                          sg,
                                          handle,
                                          params,
                                          &textureResult[0]
                                          );
    bool textureAccessY = lookupCellNoise(P.x,
                                          P.z,
                                          cellSoftness,
                                          sg,
                                          handle,
                                          params,
                                          &textureResult[1]
                                          );
    bool textureAccessZ = lookupCellNoise(P.y,
                                          P.x,
                                          cellSoftness,
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
            result = tileCellnoise(P, weights, cellSoftness, sg, data->texturehandle, data->textureparams);
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


