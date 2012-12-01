#include "Irawan.h"
#include <ai.h>
#include <cstring>

AI_SHADER_NODE_EXPORT_METHODS(alClothMtd)

struct ShaderData
{
	ShaderData(): irawan(NULL){}
	~ShaderData() {delete irawan;}
	IrawanBRDF* irawan;
};

enum alClothParams
{
	p_preset = 0,
	p_specularScale,
	p_repeatU,
	p_repeatV,
};

const char* presetNames[] =
{
	"denim",
	"silk charmeuse",
	"cotton twill",
	"wool garbadine",
	"polyester lining cloth",
	"silk shantung",
	NULL
};

node_parameters
{
	AiParameterENUM("preset", 0, presetNames);
	AiParameterFLT("specularScale", 1.0f );
	AiParameterFLT("repeatU", 10.0f );
	AiParameterFLT("repeatV", 10.0f );
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alClothMtd;
   node->output_type = AI_TYPE_RGB;
   node->name        = "alCloth";
   node->node_type   = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return TRUE;
}

node_initialize
{
	ShaderData* data = new ShaderData;
	AiNodeSetLocalData(node, data);
}

node_finish
{
	if (AiNodeGetLocalData(node))
	{
		ShaderData* data = (ShaderData*) AiNodeGetLocalData(node);
		delete data;
	}
}

node_update
{
	ShaderData* data = (ShaderData*) AiNodeGetLocalData(node);

	if (!data->irawan)
	{
		WeavePattern* wv = createWeavePreset(params[p_preset].INT);
		data->irawan = new IrawanBRDF(*wv);
	}
}

shader_evaluate
{
	ShaderData* data = (ShaderData*) AiNodeGetLocalData(node);

	AtFloat specularScale = AiShaderEvalParamFlt(p_specularScale);
	AtFloat repeatU = AiShaderEvalParamFlt(p_repeatU);
	AtFloat repeatV = AiShaderEvalParamFlt(p_repeatV);

	IrawanBRDF& iw = *(data->irawan);

	// build a local frame for sampling
	AtVector U, V;
	AiBuildLocalFramePolar(&U, &V, &sg->N);

	iw.setFrame(sg->N, U, V);
	iw.setRepeat(repeatU, repeatV);

	AtVector wo = iw.worldToLocal(-sg->Rd);

	// Light loop
	AiLightsPrepare(sg);
	AtRGB result = AI_RGB_BLACK;
	while(AiLightsGetSample(sg))
	{
		AtVector wi = iw.worldToLocal(sg->Ld);
		result += sg->Li * iw.eval(wi, wo, sg->u, sg->v) * sg->we * AiV3Dot(sg->N, sg->Ld)
						*specularScale;
	}

	sg->out.RGB = result;
}


