#include <ai.h>
#include <map>

AI_SHADER_NODE_EXPORT_METHODS(alLightGroups)

#define NUM_LIGHT_GROUPS 8
static const char* lightGroupNames[NUM_LIGHT_GROUPS] =
{
    "lightGroup1",
    "lightGroup2",
    "lightGroup3",
    "lightGroup4",
    "lightGroup5",
    "lightGroup6",
    "lightGroup7",
    "lightGroup8"
};

enum alLightGroupParams
{
    p_color
};

struct ShaderData
{
    AtSampler* sampler;
    std::map<AtNode*, int> lightGroups;
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
    // Get all the light nodes in the scene and try and find their light group parameter
    // we'll store this based on the light pointer for fast access during rendering
    AtNodeIterator* it = AiUniverseGetNodeIterator(AI_NODE_LIGHT);
    while (!AiNodeIteratorFinished(it))
    {
        AtNode* light = AiNodeIteratorGetNext(it);
        data->lightGroups[light] = AiNodeGetInt(light, "lightGroup") - 1;
    }
    AiNodeIteratorDestroy(it);
}

shader_evaluate
{
    
    AtRGB result_direct = AI_RGB_BLACK;
    AtRGB result_indirect = AI_RGB_BLACK;
    
    AtRGB color = AiShaderEvalParamRGB(p_color);

    ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);

    AtRGB* deepGroupPtr = NULL;
    AtRGB result_directGroup[NUM_LIGHT_GROUPS];
    for (int i=0; i < NUM_LIGHT_GROUPS; ++i) result_directGroup[i] = AI_RGB_BLACK;
    if (sg->Rt & AI_RAY_CAMERA)
    {
        // if this is a camera ray allocate the group storage
        deepGroupPtr = (AtRGB*)AiShaderGlobalsQuickAlloc(sg, sizeof(AtRGB)*NUM_LIGHT_GROUPS);
        memset(deepGroupPtr, 0, sizeof(AtRGB)*NUM_LIGHT_GROUPS);
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
        int lightGroup = data->lightGroups[sg->Lp];
        if (lightGroup >= 0 && lightGroup < NUM_LIGHT_GROUPS)
            result_directGroup[lightGroup] += L * color;
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
            for (int i=0; i < NUM_LIGHT_GROUPS; ++i)
            {
                deepGroupPtr[i] = deepGroupPtr[i] * f * color + result_directGroup[i];    
            }
            
        }
    }
    result_indirect *= AiSamplerGetSampleInvCount(sit);
    
    for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
    {
        AiAOVSetRGB(sg, lightGroupNames[i], deepGroupPtr[i]);
    }

    sg->out.RGB = result_direct + result_indirect;
    //sg->out.RGB = *deepGroupPtr;
}