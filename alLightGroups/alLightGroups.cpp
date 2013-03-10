#include <ai.h>
#include <map>

AI_SHADER_NODE_EXPORT_METHODS(alLightGroups)

#define NUM_LIGHT_GROUPS 8
static const char* lightGroupNames[] =
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

//#define DO_LIGHT_GROUPS

enum alLightGroupParams
{
    p_color
};

struct ShaderData
{
    AtSampler* sampler;
    int samples;
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
    data->samples = AiNodeGetInt(options, "GI_diffuse_samples");
    data->sampler = AiSampler(data->samples, 2);
    data->samples = SQR(data->samples);
#ifdef DO_LIGHT_GROUPS
    // Get all the light nodes in the scene and try and find their light group parameter
    // we'll store this based on the light pointer for fast access during rendering
    AtNodeIterator* it = AiUniverseGetNodeIterator(AI_NODE_LIGHT);
    while (!AiNodeIteratorFinished(it))
    {
        AtNode* light = AiNodeIteratorGetNext(it);
        data->lightGroups[light] = AiNodeGetInt(light, "lightGroup") - 1;
    }
    AiNodeIteratorDestroy(it);
#endif
}

shader_evaluate
{
    
    AtRGB result_direct = AI_RGB_BLACK;
    AtRGB result_indirect = AI_RGB_BLACK;
    
    AtRGB color = AiShaderEvalParamRGB(p_color);

    ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);

#ifdef DO_LIGHT_GROUPS
    AtRGB* deepGroupPtr = NULL;
    AtRGB result_directGroup[NUM_LIGHT_GROUPS];
    for (int i=0; i < NUM_LIGHT_GROUPS; ++i) result_directGroup[i] = AI_RGB_BLACK;
    if (sg->Rt & AI_RAY_CAMERA)
    {
        // if this is a camera ray allocate the group storage
        deepGroupPtr = (AtRGB*)AiShaderGlobalsQuickAlloc(sg, sizeof(AtRGB)*NUM_LIGHT_GROUPS*data->samples);
        memset(deepGroupPtr, 0, sizeof(AtRGB)*NUM_LIGHT_GROUPS*data->samples);
        AiStateSetMsgPtr("als_deepGroupPtr", deepGroupPtr);
    }
    else
    {
        // secondary ray hit - get the pointer from the state
        AiStateGetMsgPtr("als_deepGroupPtr", (void**)&deepGroupPtr);
    }
#endif

    void* mis = AiOrenNayarMISCreateData(sg, 0);
    AtRGB L = AI_RGB_BLACK;
    AiLightsPrepare(sg);
    while (AiLightsGetSample(sg))
    {
        L = AiEvaluateLightSample(sg, mis,AiOrenNayarMISSample,AiOrenNayarMISBRDF, AiOrenNayarMISPDF);

        result_direct += L * color;

#ifdef DO_LIGHT_GROUPS        
        int lightGroup = data->lightGroups[sg->Lp];
        if (lightGroup >= 0 && lightGroup < NUM_LIGHT_GROUPS)
            result_directGroup[lightGroup] += L * color;
#endif
    }

    
    double u[2];
    AtRay wir;
    AtScrSample scrs;
    AtSamplerIterator* sit = AiSamplerIterator(data->sampler, sg);
    AiMakeRay(&wir, AI_RAY_DIFFUSE, &sg->P, NULL, AI_BIG, sg);
    int count = 0;
    int idx;
    while (AiSamplerGetSample(sit, u))
    {
        wir.dir = AiOrenNayarMISSample(mis, u[0], u[1]);
        float p = AiOrenNayarMISPDF(mis, &wir.dir);

        if (p > 0)
        {
#ifdef DO_LIGHT_GROUPS
            // if we're in a camera ray, pass the sample index down to the child SG
            if (sg->Rt & AI_RAY_CAMERA)
            {
                idx = count;
                AiStateSetMsgInt("als_sampleIndex", idx);
            }
            else
            {
                AiStateGetMsgInt("als_sampleIndex", &idx);
            }
            AiTrace(&wir, &scrs);
            AtRGB f = AiOrenNayarMISBRDF(mis, &wir.dir) / p;
            L = scrs.color * f * color;
            result_indirect += L;

            for (int i=0; i < NUM_LIGHT_GROUPS; ++i)
            {
                deepGroupPtr[i*data->samples+idx] = deepGroupPtr[i*data->samples+idx] * f * color + result_directGroup[i];    
            }
#else
           AiTrace(&wir, &scrs);
            AtRGB f = AiOrenNayarMISBRDF(mis, &wir.dir) / p;
            L = scrs.color * f * color;
            result_indirect += L; 
 #endif           
        }

        count++;
    }
    result_indirect *= AiSamplerGetSampleInvCount(sit);

#ifdef DO_LIGHT_GROUPS    
    if (sg->Rt & AI_RAY_CAMERA)
    {
        AtRGB deepGroups[NUM_LIGHT_GROUPS];
        memset(deepGroups, 0, sizeof(AtRGB)*NUM_LIGHT_GROUPS);
        for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
        {
            for (int s=0; s < data->samples; ++s)
            {
                deepGroups[i] += deepGroupPtr[i*data->samples+s];
            }
            deepGroups[i] /= float(data->samples);
            AiAOVSetRGB(sg, lightGroupNames[i], deepGroups[i]);
        }
    }
#endif

    sg->out.RGB = result_direct + result_indirect;
}