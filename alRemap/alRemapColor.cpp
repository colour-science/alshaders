#include "alUtil.h"
#include <ai.h>

AI_SHADER_NODE_EXPORT_METHODS(alRemapColorMtd);

enum alRemapParams
{
	p_input,
	p_gamma,
	p_exposure,
};

node_parameters
{
	AiParameterRGB("input", 0.18f, 0.18f, 0.18f);
	AiParameterFLT("gamma", 1.0f);
	AiParameterFLT("exposure", 0.f);
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alRemapColorMtd;
   node->output_type = AI_TYPE_RGB;
   node->name        = "alRemapColor";
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
	AtRGB input = AiShaderEvalParamRGB(p_input);
	AtFloat gamma = AiShaderEvalParamFlt(p_gamma);
	AtFloat exposure = AiShaderEvalParamFlt(p_exposure);

	AtRGB result = pow(input, 1.0f/gamma);
	result = result * powf(2.0f, exposure);
	sg->out.RGB = result;
}


