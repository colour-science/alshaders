#include "alUtil.h"
#include <ai.h>

AI_SHADER_NODE_EXPORT_METHODS(alSwitchFloatMtd)

enum alSwitchParams
{
    p_inputA,
    p_inputB,
    p_inputC,
    p_inputD,
    p_inputE,
    p_inputF,
    p_inputG,
    p_inputH,
    p_mix,
    p_threshold
};

node_parameters
{
    AiParameterFLT("inputA", 0.0f);
    AiParameterFLT("inputB", 1.0f);
    AiParameterFLT("inputC", .15f);
    AiParameterFLT("inputD", .30f);
    AiParameterFLT("inputE", .45f);
    AiParameterFLT("inputF", .60f);
    AiParameterFLT("inputG", .75f);
    AiParameterFLT("inputH", .90f);
    AiParameterFLT("mixer", 1.0f);
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

    int input = floorf(mix);
    if (mix-input >= threshold) input++;
    input = clamp(input, 0, 7);

    float result = 0.0f;

    switch(input)
    {
    case 0:
        result = AiShaderEvalParamFlt(p_inputA);
        break;
    case 1:
        result = AiShaderEvalParamFlt(p_inputB);
        break;
    case 2:
        result = AiShaderEvalParamFlt(p_inputC);
        break;
    case 3:
        result = AiShaderEvalParamFlt(p_inputD);
        break;
    case 4:
        result = AiShaderEvalParamFlt(p_inputE);
        break;
    case 5:
        result = AiShaderEvalParamFlt(p_inputF);
        break;
    case 6:
        result = AiShaderEvalParamFlt(p_inputG);
        break;
    case 7:
        result = AiShaderEvalParamFlt(p_inputH);
        break;
    default:
        // should never get here
        result = 0.0f;
        break;
    }

    sg->out.FLT = result;
}


