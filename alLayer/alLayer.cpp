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

#define NUM_AOVs 31
static const char* AOVs[NUM_AOVs] = {
	"diffuse_color",
	"direct_diffuse",
	"direct_diffuse_raw",
	"indirect_diffuse",
	"indirect_diffuse_raw",
	"direct_specular",
	"indirect_specular",
	"direct_specular_2",
	"indirect_specular_2",
	"single_scatter",
	"sss",
	"refraction",
	"emission",
	"uv",
	"depth",
	"light_group_1",
	"light_group_2",
	"light_group_3",
	"light_group_4",
	"light_group_5",
	"light_group_6",
	"light_group_7",
	"light_group_8",
	"id_1",
	"id_2",
	"id_3",
	"id_4",
	"id_5",
	"id_6",
	"id_7",
	"id_8"
};


node_update
{
	AiAOVRegister("diffuse_color", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("direct_diffuse", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("indirect_diffuse", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("direct_diffuse_raw", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("indirect_diffuse_raw", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("direct_specular", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("indirect_specular", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("direct_specular_2", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("indirect_specular_2", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("single_scatter", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("sss", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("refraction", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("emission", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("uv", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("depth", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("light_group_1", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("light_group_2", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("light_group_3", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("light_group_4", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("light_group_5", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("light_group_6", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("light_group_7", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("light_group_8", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
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
			if (sg->Rt & AI_RAY_CAMERA) // handle aovs
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
			else // just layer the results
			{
				AtRGB layer1 = AiShaderEvalParamRGB(p_layer1);
				AtRGB layer2 = AiShaderEvalParamRGB(p_layer2);
				result = lerp(layer1, layer2, mix);
			}
		}
	}

	sg->out.RGB = result;
}


