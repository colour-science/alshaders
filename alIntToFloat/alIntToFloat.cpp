#include "alUtil.h"
#include <ai.h>

AI_SHADER_NODE_EXPORT_METHODS(alIntToFloatMtd)

enum alSwitchParams
{
    p_inputA
};

node_parameters
{
    AiParameterINT("input", 0);
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alSwitchColorMtd;
   node->output_type = AI_TYPE_FLOAT;
   node->name        = "alIntToFloat";
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
    sg->out.FLT = AiShaderEvalParamINT(p_input);
}


