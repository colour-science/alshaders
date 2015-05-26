#include <ai.h>
#include <cstring>
#include "tea.h"
#include "alUtil.h"

AI_SHADER_NODE_EXPORT_METHODS(alFlake)

enum alFlakeParams
{
   p_amount = 0,
   p_size,
   p_divergence,
   p_P
};

node_parameters
{
   AiParameterFlt("amount", 0.7f);
   AiParameterFlt("size", 0.01f);
   AiParameterFlt("divergence", 0.5f);
   AiParameterPnt("P", 0.0f, 0.0, 0.0);
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alFlake;
   node->output_type = AI_TYPE_VECTOR;
   node->name        = "alFlake";
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
   AtVector result = AiVector(0.0f, 0.0f, 1.0f);

   float amount = clamp(AiShaderEvalParamFlt(p_amount), 0.0f, 1.0f);
   float size = AiShaderEvalParamFlt(p_size);
   float divergence = clamp(AiShaderEvalParamFlt(p_divergence), 0.0f, 1.0f);

   AtPoint P = sg->P;
   if (AiNodeIsLinked(node, "P"))
   {
      P = AiShaderEvalParamPnt(p_P);
   }

   P /= size;

   // get a cellular id
   AtUInt32 id;
   AiCellular(P, 1, 1.0f, 2.0f, 1.0f, NULL, NULL, &id);

   // get two random numbers based on the id
   float u1 = sampleTEAFloat(id, 0);
   if (u1 < amount)
   {
      u1 /= amount;
      float u2 = sampleTEAFloat(id, 1);

      // get a random vector in the hemisphere
      AtVector d = uniformSampleHemisphere(u1, u2);
      // swap y and z to make a normal map
      std::swap(d.y, d.z);

      // blend it in 
      result = lerp(result, d, divergence);
   }
   sg->out.VEC = result * .5f + AiVector(0.5f, 0.5f, 0.5f);
}


