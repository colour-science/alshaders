#include "alUtil.h"
#include <ai.h>

AI_SHADER_NODE_EXPORT_METHODS(alLayer)

enum alLayerParams
{
	p_layer1,
	p_layer2,
	p_mix
};

node_parameters
{
	AiParameterRGB("layer1", 0.0f, 0.0f, 0.0f);
	AiParameterRGB("layer2", 0.0f, 0.0f, 0.0f);
	AiParameterFlt("mix", 0.0f);
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

node_update
{

}

shader_evaluate
{
	AtRGB result = AI_RGB_BLACK;
	
	AtFloat mix = AiShaderEvalParamFlt(p_mix);
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
		AtRGB layer1 = AiShaderEvalParamRGB(p_layer1);
		AtRGB layer2 = AiShaderEvalParamRGB(p_layer2);
		result = lerp(layer1, layer2, mix);
	}

	sg->out.RGB = result;
}


