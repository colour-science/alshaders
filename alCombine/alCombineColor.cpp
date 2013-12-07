#include "alUtil.h"
#include <ai.h>

AI_SHADER_NODE_EXPORT_METHODS(alCombineColorMtd)

enum alRemapParams
{
	p_input1,
	p_input2,
	p_input3,
	p_combineOp
};

enum CombineOpEnum
{
	CO_MULTIPLY=0,
	CO_ADD,
	CO_DIVIDE,
	CO_SUBTRACT,
	CO_LERP
};

static const char* combineOpNames[] =
{
	"multiply 1*2",
	"add 1+2",
	"divide 1/2",
	"subtract 1-2",
	"lerp(1, 2, 3)",
	NULL
};

node_parameters
{
	AiParameterRGB("input1", 1.0f, 1.0f, 1.0f);
	AiParameterRGB("input2", 1.0f, 1.0f, 1.0f);
	AiParameterFLT("input3", 0.0f);
	AiParameterENUM("combineOp", 0, combineOpNames);
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alCombineColorMtd;
   node->output_type = AI_TYPE_RGB;
   node->name        = "alCombineColor";
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
	AtRGB input1 = AiShaderEvalParamRGB(p_input1);
	AtRGB input2 = AiShaderEvalParamRGB(p_input2);
	float input3 = AiShaderEvalParamFlt(p_input3);
	int combineOp = AiShaderEvalParamInt(p_combineOp);

	AtRGB f = input1;
	switch(combineOp)
	{
	case CO_MULTIPLY:
		f = input1*input2;
		break;
	case CO_ADD:
		f = input1+input2;
		break;
	case CO_DIVIDE:
		f = input1/input2;
		break;
	case CO_SUBTRACT:
		f = input1-input2;
		break;
	case CO_LERP:
		f = lerp(input1, input2, input3);
		break;
	default:
		f = input1;
		break;
	}

	sg->out.RGB = f;
}


