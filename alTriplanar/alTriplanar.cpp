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

enum alTriplanarParams
{
	p_space,
    p_texture,
    p_blendSoftness,
    p_scale,
	p_offset,
    p_P
};

node_parameters
{
	AiParameterENUM("space", 0, triplanarSpaceNames);
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

shader_evaluate
{
	int space = AiShaderEvalParamInt(p_space);
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
        // blend result
        AtRGB textureResult = AI_RGB_BLACK;
        float weights[3] = {fabsf(sg->N.x), fabsf(sg->N.y), fabsf(sg->N.z)};
        float sum = 0.f;
        int highest = 0;
        for(int i=0; i<3; ++i){
            weights[i] = weights[i] - (1.f-blendSoftness)/2.f;
            weights[i] = MAX(weights[i], 0.00f);
            sum += weights[i];
        }
        if(sum){
            for(int i=0; i<3; ++i){
                weights[i] /= sum;
            }
        }
        textureResult += textureResultX.rgb() * weights[0];
        textureResult += textureResultY.rgb() * weights[1];
        textureResult += textureResultZ.rgb() * weights[2];
        sg->out.RGB = textureResult;
    }
    else {
        sg->out.RGB = AiColorCreate(1.f, 0.0f, 0.f);
    }

}


