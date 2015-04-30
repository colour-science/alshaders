#include "alUtil.h"
#include "Remap.h"
#include <ai.h>

AI_SHADER_NODE_EXPORT_METHODS(alRemapFloatMtd)

enum alRemapParams
{
	p_input,
	REMAP_FLOAT_PARAM_ENUM
};

node_parameters
{
	AiParameterFLT("input", 0.0f);
	REMAP_FLOAT_PARAM_DECLARE;
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alRemapFloatMtd;
   node->output_type = AI_TYPE_FLOAT;
   node->name        = "alRemapFloat";
   node->node_type   = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return true;
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
	float input = AiShaderEvalParamFlt(p_input);

	RemapFloat r = REMAP_FLOAT_CREATE;

	sg->out.FLT = r.remap(input);
}


