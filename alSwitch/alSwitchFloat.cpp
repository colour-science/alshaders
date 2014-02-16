#include "alUtil.h"
#include <ai.h>

AI_SHADER_NODE_EXPORT_METHODS(alSwitchFloatMtd)

enum alSwitchParams
{
    p_inputA,
    p_inputB,
    p_mix,
    p_threshold
};

node_parameters
{
    AiParameterFLT("inputA", 0.f);
    AiParameterFLT("inputB", 1.f);
    AiParameterFLT("mix", 1.0f);
    AiParameterFLT("threshold", 0.5f);
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alSwitchFloatMtd;
   node->output_type = AI_TYPE_FLOAT;
   node->name        = "alSwitchFloat";
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
    float mix = AiShaderEvalParamFlt(p_mix);
    float threshold = AiShaderEvalParamFlt(p_threshold);

    if(mix <= threshold){
        sg->out.FLT = AiShaderEvalParamFlt(p_inputA);
    } else {
        sg->out.FLT = AiShaderEvalParamFlt(p_inputB);
    }
}


