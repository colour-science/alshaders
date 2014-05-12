#include "Remap.h"
#include <ai.h>

AI_SHADER_NODE_EXPORT_METHODS(alCurvatureMtd)

struct ShaderData
{
	AtSampler* sampler;
};

enum alCurvatureParams
{
	p_samples,
	p_sampleOffset,
	p_sampleRadius,
	p_color1,
	p_color2,
	REMAP_FLOAT_PARAM_ENUM
};

node_parameters
{
	AiParameterINT("samples", 2);
	AiParameterFLT("sampleOffset", 1.0f);
	AiParameterFLT("sampleRadius", 1.0f);
	AiParameterRGB("color1", 0.0f, 0.0f, 0.0f);
	AiParameterRGB("color2", 1.0f, 1.0f, 1.0f);
	REMAP_FLOAT_PARAM_DECLARE;
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alCurvatureMtd;
   node->output_type = AI_TYPE_RGB;
   node->name        = "alCurvature";
   node->node_type   = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return true;
}

node_initialize
{
	ShaderData *data = (ShaderData*) AiMalloc(sizeof(ShaderData));
	AiNodeSetLocalData(node,data);
	data->sampler = NULL;
}

node_finish
{
	if (AiNodeGetLocalData(node))
	{
		ShaderData* data = (ShaderData*) AiNodeGetLocalData(node);
		AiSamplerDestroy(data->sampler);

		AiFree((void*) data);
		AiNodeSetLocalData(node, NULL);
	}
}

node_update
{
	ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);
	AiSamplerDestroy(data->sampler);
	data->sampler = AiSampler(params[p_samples].INT, 2);
	
}

shader_evaluate
{
	float result = 0.0f;
	float sampleOffset = AiShaderEvalParamFlt(p_sampleOffset);
	float sampleRadius = AiShaderEvalParamFlt(p_sampleRadius);
	AtRGB color1 = AiShaderEvalParamRGB(p_color1);
	AtRGB color2 = AiShaderEvalParamRGB(p_color2);
	

	
	// build a local frame for sampling
	AtVector U, V;
	AiBuildLocalFramePolar(&U, &V, &sg->N);
	
	ShaderData* data = (ShaderData*) AiNodeGetLocalData(node);
	AtSamplerIterator* sampit = AiSamplerIterator(data->sampler, sg);
	AtRay ray;
	AtShaderGlobals* hitpoint = AiShaderGlobals();
#if AI_VERSION_MAJOR_NUM > 0
    float samples[2];
#else
    double samples[2];
#endif
	float du, dv;
	AtVector dir = -sg->N;
	float count = 0.0f;
	while (AiSamplerGetSample(sampit, samples))
	{
		// sample a disk above the normal to get the src points
		//concentricSampleDisk(samples[0], samples[1], du, dv);
		//AtPoint srcpoint = sg->P + du*U*sampleRadius + dv*V*sampleRadius + sampleOffset*sg->N;
		AtPoint srcpoint = sg->P + sg->N*sampleOffset + uniformSampleSphere(samples[0], samples[1]);

		// trace straight back down
		AiMakeRay(&ray, AI_RAY_GENERIC, &srcpoint, &dir, sampleRadius, sg);
		AiTraceProbe(&ray, hitpoint);
		if (hitpoint)
		{
			AtVector L = AiV3Normalize(hitpoint->P - sg->P);
			result += AiV3Dot(hitpoint->N, L) * fabs(AiV3Dot(hitpoint->N, sg->N));
			//result += AiV3Dot(hitpoint->N, sg->N);
			count++;
		}		
	}
	if (count)
		result /= count;
		
	RemapFloat r = REMAP_FLOAT_CREATE;
	result = r.remap(result);

	sg->out.RGB = AiColorCreate(std::max(result, 0.0f), std::max(-result, 0.0f), 0.0f);
	AiShaderGlobalsDestroy(hitpoint);
}


