#include <ai.h>

AI_SHADER_NODE_EXPORT_METHODS(alLightGroups)

enum alLightGroupParams
{
    p_color
};

struct ShaderData
{
    AtSampler* sampler;
};

node_parameters
{
    AiParameterRGB("color", 0.18f, 0.18f, 0.18f);
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alLightGroups;
   node->output_type = AI_TYPE_RGB;
   node->name        = "alLightGroups";
   node->node_type   = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return TRUE;
}

node_initialize
{
    ShaderData *data = new ShaderData;
    AiNodeSetLocalData(node,data);
    data->sampler = NULL;
}

node_finish
{
    if (AiNodeGetLocalData(node))
    {
        ShaderData* data = (ShaderData*) AiNodeGetLocalData(node);

        AiSamplerDestroy(data->sampler);
        AiNodeSetLocalData(node, NULL);
        delete data;
    }
}

node_update
{
    ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);
    AtNode *options   = AiUniverseGetOptions();
    AiSamplerDestroy(data->sampler);
    data->sampler = AiSampler(AiNodeGetInt(options, "GI_diffuse_samples"), 2);
}

shader_evaluate
{
    
    AtRGB result_direct = AI_RGB_BLACK;
    AtRGB result_indirect = AI_RGB_BLACK;
    
    AtRGB color = AiShaderEvalParamRGB(p_color);

    ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);

    AtRGB* deepGroupPtr = NULL;
    if (sg->Rt & AI_RAY_CAMERA)
    {
        // if this is a camera ray allocate the group storage
        deepGroupPtr = (AtRGB*)AiShaderGlobalsQuickAlloc(sg, sizeof(AtRGB));
        memset(deepGroupPtr, 0, sizeof(AtRGB));
        AiStateSetMsgPtr("als_deepGroupPtr", deepGroupPtr);
    }
    else
    {
        // secondary ray hit - get the pointer from the state
        AiStateGetMsgPtr("als_deepGroupPtr", (void**)&deepGroupPtr);
    }
    
    void* mis = AiOrenNayarMISCreateData(sg, 0);
    AtRGB L = AI_RGB_BLACK;
    AiLightsPrepare(sg);
    while (AiLightsGetSample(sg))
    {
        L = AiEvaluateLightSample(sg, mis,AiOrenNayarMISSample,AiOrenNayarMISBRDF, AiOrenNayarMISPDF);

        result_direct += L * color;
    }
    
    double u[2];
    AtRay wir;
    AtScrSample scrs;
    AtSamplerIterator* sit = AiSamplerIterator(data->sampler, sg);
    AiMakeRay(&wir, AI_RAY_DIFFUSE, &sg->P, NULL, AI_BIG, sg);
    while (AiSamplerGetSample(sit, u))
    {
        wir.dir = AiOrenNayarMISSample(mis, u[0], u[1]);
        float p = AiOrenNayarMISPDF(mis, &wir.dir);

        if (p > 0)
        {
            AiTrace(&wir, &scrs);
            AtRGB f = AiOrenNayarMISBRDF(mis, &wir.dir) / p;
            L = scrs.color * f * color;
            result_indirect += L;
            *deepGroupPtr = *deepGroupPtr * f * color + result_direct;
        }
    }
    result_indirect *= AiSamplerGetSampleInvCount(sit);
    

    //sg->out.RGB = result_direct + result_indirect;
    sg->out.RGB = *deepGroupPtr;
}