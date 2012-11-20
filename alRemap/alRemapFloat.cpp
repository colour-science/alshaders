#include "alUtil.h"
#include <ai.h>

AI_SHADER_NODE_EXPORT_METHODS(alRemapFloatMtd)

enum alRemapParams
{
	p_input,
	p_inputMin,
	p_inputMax,
	p_outputMin,
	p_outputMax
};

node_parameters
{
	AiParameterFLT("input", 0.0f);
	AiParameterFLT("inputMin", 0.0f);
	AiParameterFLT("inputMax", 1.0f);
	AiParameterFLT("outputMin", 0.0f);
	AiParameterFLT("outputMax", 1.0f);
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alRemapFloatMtd;
   node->output_type = AI_TYPE_FLOAT;
   node->name        = "alRemapFloat";
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
	AtFloat input = AiShaderEvalParamFlt(p_input);
	AtFloat inputMin = AiShaderEvalParamFlt(p_inputMin);
	AtFloat inputMax = AiShaderEvalParamFlt(p_inputMax);
	AtFloat outputMin = AiShaderEvalParamFlt(p_outputMin);
	AtFloat outputMax = AiShaderEvalParamFlt(p_outputMax);

	AtFloat f = (input-inputMin)/(inputMax-inputMin);
	f = lerp(outputMin, outputMax, f);
	sg->out.FLT = f;
}


