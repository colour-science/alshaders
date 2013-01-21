#include "alUtil.h"
#include <ai.h>

AI_SHADER_NODE_EXPORT_METHODS(alLayer)

enum alLayerParams
{
	p_layer1,
	p_layer2,
	p_mix,
	p_debug
};

enum DebugModes
{
	kOff = 0,
	kLayer1,
	kLayer2,
	kMixer
};

static const char* debugModeNames[] =
{
	"off",
	"layer1",
	"layer2",
	"mixer",
	NULL
};

node_parameters
{
	AiParameterRGB("layer1", 0.0f, 0.0f, 0.0f);
	AiParameterRGB("layer2", 0.0f, 0.0f, 0.0f);
	AiParameterFlt("mix", 0.0f);
	AiParameterEnum("debug", kOff, debugModeNames);
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alLayer;
   node->output_type = AI_TYPE_RGB;
   node->name        = "alLayer";
   node->node_type   = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return TRUE;
}

node_initialize
{

}

node_finish
{

}

#define NUM_AOVs 20
static const char* AOVs[NUM_AOVs] = {
	"diffuseDirect",
	"diffuseIndirect",
	"specularDirect",
	"specularIndirect",
	"specular2Direct",
	"specular2Indirect",
	"singleScatter",
	"multiScatter",
	"transmission",
	"emission",
	"uv",
	"depth",
	"lightGroup1",
	"lightGroup2",
	"lightGroup3",
	"lightGroup4",
	"lightGroup5",
	"lightGroup6",
	"lightGroup7",
	"lightGroup8"
};


node_update
{
	AiAOVRegister("diffuseDirect", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("diffuseIndirect", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("specularDirect", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("specularIndirect", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("specular2Direct", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("specular2Indirect", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("singleScatter", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("multiScatter", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("transmission", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("emission", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("uv", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("depth", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("lightGroup1", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("lightGroup2", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("lightGroup3", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("lightGroup4", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("lightGroup5", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("lightGroup6", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("lightGroup7", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("lightGroup8", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
}

shader_evaluate
{
	AtRGB result = AI_RGB_BLACK;
	
	AtFloat mix = AiShaderEvalParamFlt(p_mix);
	int debug = AiShaderEvalParamEnum(p_debug);
	if (debug == kMixer)
	{
		result = AiColorCreate(mix, mix, mix);
	}
	else
	{
		if (debug == kLayer1) mix = 0.0f;
		else if (debug == kLayer2) mix = 1.0f;


		if (mix >= (1.0f-IMPORTANCE_EPS))
		{
			result = AiShaderEvalParamRGB(p_layer2);
		}
		else if (mix <= IMPORTANCE_EPS)
		{
			result = AiShaderEvalParamRGB(p_layer1);
		}
		else
		{
			AtRGB tmp[NUM_AOVs];
			AtRGB layer1 = AiShaderEvalParamRGB(p_layer1);
			for (int i=0; i < NUM_AOVs; ++i)
			{
				if (!AiAOVGetRGB(sg, AOVs[i], tmp[i]))
				{
					tmp[i] = AI_RGB_BLACK;
				}
			}
			AtRGB layer2 = AiShaderEvalParamRGB(p_layer2);
			result = lerp(layer1, layer2, mix);
			for (int i=0; i < NUM_AOVs; ++i)
			{
				AtRGB tmp2;
				if (!AiAOVGetRGB(sg, AOVs[i], tmp2))
				{
					tmp2 = AI_RGB_BLACK;
				}
				AiAOVSetRGB(sg, AOVs[i], lerp(tmp[i], tmp2, mix));
			}
		}
	}

	sg->out.RGB = result;
}


