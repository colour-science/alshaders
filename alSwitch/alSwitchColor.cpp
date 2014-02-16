#include "alUtil.h"
#include <ai.h>

AI_SHADER_NODE_EXPORT_METHODS(alSwitchColorMtd)

enum alSwitchParams
{
    p_inputA,
    p_inputB,
    p_mix,
    p_threshold
};

node_parameters
{
    AiParameterRGB("inputA", 0.f, 0.f, 0.f);
    AiParameterRGB("inputB", 1.f, 1.f, 1.f);
    AiParameterFLT("mix", 1.0f);
    AiParameterFLT("threshold", 0.5f);
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alSwitchColorMtd;
   node->output_type = AI_TYPE_RGB;
   node->name        = "alSwitchColor";
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
        sg->out.RGB = AiShaderEvalParamRGB(p_inputA);
    } else {
        sg->out.RGB = AiShaderEvalParamRGB(p_inputB);
    }
}


