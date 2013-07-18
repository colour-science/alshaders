#include "alUtil.h"
#include <ai.h>
#include <vector>
#include <string>
#include <cassert>

AI_SHADER_NODE_EXPORT_METHODS(alLayer)

struct ShaderData
{
	// AOV names
   std::vector<std::string> aovs;

   bool standardAovs;
};

enum alLayerParams
{
	p_layer1,
	p_layer2,
	p_mix,
	p_debug,


	p_aov_diffuse_color,
    p_aov_direct_diffuse,
    p_aov_direct_diffuse_raw,
    p_aov_indirect_diffuse,
    p_aov_indirect_diffuse_raw,
    p_aov_direct_backlight,
    p_aov_indirect_backlight,
    p_aov_direct_specular,
    p_aov_indirect_specular,
    p_aov_direct_specular_2,
    p_aov_indirect_specular_2,
    p_aov_single_scatter,
    p_aov_sss,
    p_aov_refraction,
    p_aov_emission,
    p_aov_uv,
    p_aov_depth,
    p_aov_light_group_1,
    p_aov_light_group_2,
    p_aov_light_group_3,
    p_aov_light_group_4,
    p_aov_light_group_5,
    p_aov_light_group_6,
    p_aov_light_group_7,
    p_aov_light_group_8,
    p_aov_id_1,
    p_aov_id_2,
    p_aov_id_3,
    p_aov_id_4,
    p_aov_id_5,
    p_aov_id_6,
    p_aov_id_7,
    p_aov_id_8,

    p_standardAovs
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

	AiParameterStr("aov_diffuse_color", "diffuse_color");
    AiParameterStr("aov_direct_diffuse", "direct_diffuse");
    AiParameterStr("aov_direct_diffuse_raw", "direct_diffuse_raw");
    AiParameterStr("aov_indirect_diffuse", "indirect_diffuse");
    AiParameterStr("aov_indirect_diffuse_raw", "indirect_diffuse_raw");
    AiParameterStr("aov_direct_backlight", "direct_backlight");
    AiParameterStr("aov_indirect_backlight", "indirect_backlight");
    AiParameterStr("aov_direct_specular", "direct_specular");
    AiParameterStr("aov_indirect_specular", "indirect_specular");
    AiParameterStr("aov_direct_specular_2", "direct_specular_2");
    AiParameterStr("aov_indirect_specular_2", "indirect_specular_2");
    AiParameterStr("aov_single_scatter", "single_scatter");
    AiParameterStr("aov_sss", "sss");
    AiParameterStr("aov_refraction", "refraction");
    AiParameterStr("aov_emission", "emission");
    AiParameterStr("aov_uv", "uv");
    AiParameterStr("aov_depth", "depth");
    AiParameterStr("aov_light_group_1", "light_group_1");
    AiParameterStr("aov_light_group_2", "light_group_2");
    AiParameterStr("aov_light_group_3", "light_group_3");
    AiParameterStr("aov_light_group_4", "light_group_4");
    AiParameterStr("aov_light_group_5", "light_group_5");
    AiParameterStr("aov_light_group_6", "light_group_6");
    AiParameterStr("aov_light_group_7", "light_group_7");
    AiParameterStr("aov_light_group_8", "light_group_8");
    AiParameterStr("aov_id_1", "id_1");
    AiParameterStr("aov_id_2", "id_2");
    AiParameterStr("aov_id_3", "id_3");
    AiParameterStr("aov_id_4", "id_4");
    AiParameterStr("aov_id_5", "id_5");
    AiParameterStr("aov_id_6", "id_6");
    AiParameterStr("aov_id_7", "id_7");
    AiParameterStr("aov_id_8", "id_8");

    AiParameterBool("standardCompatibleAOVs", false);
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
	ShaderData* data = new ShaderData;
    AiNodeSetLocalData(node,data);
}

node_finish
{
	if (AiNodeGetLocalData(node))
    {
        ShaderData* data = (ShaderData*) AiNodeGetLocalData(node);
        AiNodeSetLocalData(node, NULL);
        delete data;
    }
}

#define NUM_AOVs 33

node_update
{
	ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);

    // set up AOVs
    data->aovs.clear();
    data->aovs.push_back(params[p_aov_diffuse_color].STR);
    data->aovs.push_back(params[p_aov_direct_diffuse].STR);
    data->aovs.push_back(params[p_aov_direct_diffuse_raw].STR);
    data->aovs.push_back(params[p_aov_indirect_diffuse].STR);
    data->aovs.push_back(params[p_aov_indirect_diffuse_raw].STR);
    data->aovs.push_back(params[p_aov_direct_specular].STR);
    data->aovs.push_back(params[p_aov_indirect_specular].STR);
    data->aovs.push_back(params[p_aov_direct_specular_2].STR);
    data->aovs.push_back(params[p_aov_indirect_specular_2].STR);
    data->aovs.push_back(params[p_aov_direct_backlight].STR);
    data->aovs.push_back(params[p_aov_indirect_backlight].STR);
    data->aovs.push_back(params[p_aov_single_scatter].STR);
    data->aovs.push_back(params[p_aov_sss].STR);
    data->aovs.push_back(params[p_aov_refraction].STR);
    data->aovs.push_back(params[p_aov_emission].STR);
    data->aovs.push_back(params[p_aov_uv].STR);
    data->aovs.push_back(params[p_aov_depth].STR);
    data->aovs.push_back(params[p_aov_light_group_1].STR);
    data->aovs.push_back(params[p_aov_light_group_2].STR);
    data->aovs.push_back(params[p_aov_light_group_3].STR);
    data->aovs.push_back(params[p_aov_light_group_4].STR);
    data->aovs.push_back(params[p_aov_light_group_5].STR);
    data->aovs.push_back(params[p_aov_light_group_6].STR);
    data->aovs.push_back(params[p_aov_light_group_7].STR);
    data->aovs.push_back(params[p_aov_light_group_8].STR);
    data->aovs.push_back(params[p_aov_id_1].STR);
    data->aovs.push_back(params[p_aov_id_2].STR);
    data->aovs.push_back(params[p_aov_id_3].STR);
    data->aovs.push_back(params[p_aov_id_4].STR);
    data->aovs.push_back(params[p_aov_id_5].STR);
    data->aovs.push_back(params[p_aov_id_6].STR);
    data->aovs.push_back(params[p_aov_id_7].STR);
    data->aovs.push_back(params[p_aov_id_8].STR);

    assert(NUM_AOVs == data->aovs.size() && "[alLayer] NUM_AOVs does not match size of aovs array!");

    for (size_t i=0; i < data->aovs.size(); ++i)
    	AiAOVRegister(data->aovs[i].c_str(), AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);

    data->standardAovs = params[p_standardAovs].BOOL;
}

shader_evaluate
{
	ShaderData* data = (ShaderData*) AiNodeGetLocalData(node);

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
				for (size_t i=0; i < data->aovs.size(); ++i)
				{
					if (!AiAOVGetRGB(sg, data->aovs[i].c_str(), tmp[i]))
					{
						tmp[i] = AI_RGB_BLACK;
					}
				}
				AtRGB layer2 = AiShaderEvalParamRGB(p_layer2);
				result = lerp(layer1, layer2, mix);
				for (size_t i=0; i < data->aovs.size(); ++i)
				{
					AtRGB tmp2;
					if (!AiAOVGetRGB(sg, data->aovs[i].c_str(), tmp2))
					{
						tmp2 = AI_RGB_BLACK;
					}
					AiAOVSetRGB(sg, data->aovs[i].c_str(), lerp(tmp[i], tmp2, mix));
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


