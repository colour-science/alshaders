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

inline bool tileCellNoise(float u, float v,  AtRGBA *textureResult, AtShaderGlobals *sg,
                          AtTextureHandle *handle, AtTextureParams *params){
        // lookup X
    AtPoint PX;
    PX.x = u;
    PX.y = v;
    PX.z = 0.f;

    // run cellnoise
    float f;
    AtVector delta;
    AtUInt32 id;
    AiCellular(PX, 1, 1, 1.92, 1, &f, &delta, &id);

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

shader_evaluate
{
	int space = AiShaderEvalParamInt(p_space);
    int tiling = AiShaderEvalParamInt(p_tiling);
    float blendSoftness = AiShaderEvalParamFlt(p_blendSoftness);
    blendSoftness = CLAMP(blendSoftness, 0.f, 1.f);
    float scale = AiShaderEvalParamFlt(p_scale);
	float offset = AiShaderEvalParamFlt(p_offset);
    AtPoint Pin = AiShaderEvalParamPnt(p_P);

    ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);

    // choose what space we want to calculate in
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

        // compute blend weights
    AtRGB textureResult = AI_RGB_BLACK;
    float weights[3] = {fabsf(sg->N.x), fabsf(sg->N.y), fabsf(sg->N.z)};
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


    if(tiling == TM_CELLNOISE){

        AtRGBA textureResult[3];
        bool textureAccessX = tileCellNoise((P.y + offset)/scale,
                                            (P.z + offset)/scale,
                                            &textureResult[0],
                                            sg,
                                            data->texturehandle,
                                            data->textureparams);
        bool textureAccessY = tileCellNoise((P.x + offset)/scale,
                                            (P.z + offset)/scale,
                                            &textureResult[1],
                                            sg,
                                            data->texturehandle,
                                            data->textureparams);
        bool textureAccessZ = tileCellNoise((P.y + offset)/scale,
                                            (P.x + offset)/scale,
                                            &textureResult[2],
                                            sg,
                                            data->texturehandle,
                                            data->textureparams);


        if(textureAccessX && textureAccessY && textureAccessZ){
            AtRGB result = AI_RGB_BLACK;
            result += textureResult[0].rgb() * weights[0];
            result += textureResult[1].rgb() * weights[1];
            result += textureResult[2].rgb() * weights[2];
            sg->out.RGB = result;
        }
        else {
            sg->out.RGB = AiColorCreate(1.f, 0.0f, 0.f);
        }

    } else {
        // regular tiling
        AtRGBA textureResultX;
        AtRGBA textureResultY;
        AtRGBA textureResultZ;
        bool textureAccessX = false;
        bool textureAccessY = false;
        bool textureAccessZ = false;

        // lookup X
        sg->u = (P.z + 123.94) / scale + offset;
        sg->v = (P.y + 87.22) / scale + offset;
        textureResultX = AiTextureHandleAccess(sg, data->texturehandle, data->textureparams, &textureAccessX);

        // lookup Y
        sg->u = (P.z + 9.2) / scale + offset;
        sg->v = (P.x + 74.1) / scale + offset;
        textureResultY = AiTextureHandleAccess(sg, data->texturehandle, data->textureparams, &textureAccessY);

        // lookup Z
        sg->u = (P.x - 12.55) / scale + offset;
        sg->v = (P.y + 6) / scale + offset;
        textureResultZ = AiTextureHandleAccess(sg, data->texturehandle, data->textureparams, &textureAccessZ);

        if(textureAccessX && textureAccessY && textureAccessZ){
            textureResult += textureResultX.rgb() * weights[0];
            textureResult += textureResultY.rgb() * weights[1];
            textureResult += textureResultZ.rgb() * weights[2];
            sg->out.RGB = textureResult;
        }
        else {
            sg->out.RGB = AiColorCreate(1.f, 0.0f, 0.f);
        }
    }

}


