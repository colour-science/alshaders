#include "Irawan.h"
#include <ai.h>
#include <cstring>

AI_SHADER_NODE_EXPORT_METHODS(alClothMtd)

struct ShaderData
{
	ShaderData(): irawan(NULL), diffuse_sampler(NULL), glossy_sampler(NULL){}
	~ShaderData() {delete irawan;}
	IrawanBRDF* irawan;
	AtSampler* diffuse_sampler;
	AtSampler* glossy_sampler;

};

enum alClothParams
{
	p_preset = 0,
	p_diffuseScale,
	p_specularScale,
	p_repeatU,
	p_repeatV,
	p_warpDiffuseColor,
	p_weftDiffuseColor,
	p_warpSpecularColor,
	p_weftSpecularColor
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
	AiParameterFLT("diffuseScale", 1.0f );
	AiParameterFLT("specularScale", 1.0f );
	AiParameterFLT("repeatU", 10.0f );
	AiParameterFLT("repeatV", 10.0f );
	AiParameterRGB("warpDiffuseColor", 0.6f, 0.6f, 0.6f);
	AiParameterRGB("weftDiffuseColor", 0.6f, 0.6f, 0.6f);
	AiParameterRGB("warpSpecularColor", 0.6f, 0.6f, 0.6f);
	AiParameterRGB("weftSpecularColor", 0.6f, 0.6f, 0.6f);
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
		AiSamplerDestroy(data->diffuse_sampler);
		AiSamplerDestroy(data->glossy_sampler);
		delete data;
		AiNodeSetLocalData(node, NULL);
	}
}

node_update
{
	ShaderData* data = (ShaderData*) AiNodeGetLocalData(node);
	AtNode *options   = AiUniverseGetOptions();

	// create cloth BRDF based on preset
	// TODO: add weave file support
	if (!data->irawan)
	{
		WeavePattern* wv = createWeavePreset(params[p_preset].INT);
		data->irawan = new IrawanBRDF(*wv);
	}

	// set up samplers
	AtInt diffuseSamples = AiNodeGetInt(options, "GI_diffuse_samples");
	AtInt glossySamples = AiNodeGetInt(options, "GI_glossy_samples");
	AiSamplerDestroy(data->diffuse_sampler);
	AiSamplerDestroy(data->glossy_sampler);
	data->diffuse_sampler = AiSampler(diffuseSamples, 2);
	data->glossy_sampler = AiSampler(glossySamples, 2);
}

shader_evaluate
{
	ShaderData* data = (ShaderData*) AiNodeGetLocalData(node);

	AtFloat diffuseScale = AiShaderEvalParamFlt(p_diffuseScale);
	AtFloat specularScale = AiShaderEvalParamFlt(p_specularScale);
	AtFloat repeatU = AiShaderEvalParamFlt(p_repeatU);
	AtFloat repeatV = AiShaderEvalParamFlt(p_repeatV);
	AtRGB warpDiffuseColor = AiShaderEvalParamRGB(p_warpDiffuseColor);
	AtRGB weftDiffuseColor = AiShaderEvalParamRGB(p_weftDiffuseColor);
	AtRGB warpSpecularColor = AiShaderEvalParamRGB(p_warpSpecularColor);
	AtRGB weftSpecularColor = AiShaderEvalParamRGB(p_weftSpecularColor);

	IrawanBRDF& iw = *(data->irawan);

	// Figure out what parts of calculate
	bool do_diffuse=true, do_glossy=true;
	if ( sg->Rr_diff > 0)
	{
		do_glossy = false;
	}

	// build a local frame for sampling
	AtVector U, V;
	AiBuildLocalFramePolar(&U, &V, &sg->Nf);

	iw.setFrame(sg->Nf, U, V);
	iw.setRepeat(repeatU, repeatV);
	iw.setColors(warpDiffuseColor, weftDiffuseColor, warpSpecularColor, weftSpecularColor);

	AtVector wo = iw.worldToLocal(-sg->Rd);

	// initialize accumulators
	AtRGB result_diffuseDirect = AI_RGB_BLACK;
	AtRGB result_diffuseIndirect = AI_RGB_BLACK;
	AtRGB result_specularDirect = AI_RGB_BLACK;
	AtRGB result_specularIndirect = AI_RGB_BLACK;

	// Light loop
	AiLightsPrepare(sg);
	AtRGB tempDiffuse;
	while(AiLightsGetSample(sg))
	{
		AtVector wi = iw.worldToLocal(sg->Ld);

		if (do_glossy)
		{
			result_specularDirect += sg->Li * iw.eval(wi, wo, sg->u, sg->v, tempDiffuse) * sg->we * AiV3Dot(sg->Nf, sg->Ld)
						*specularScale;
			result_diffuseDirect += sg->Li * tempDiffuse * sg->we * AiV3Dot(sg->Nf, sg->Ld) * diffuseScale;
		}
	}

	// brdf sampling
	double samples[2];
	AtRay wi_ray;
	AtVector wi, wi_l;
	AtScrSample scrs;
	if (do_glossy)
	{
		AtSamplerIterator* sampit = AiSamplerIterator(data->glossy_sampler, sg);
		AiMakeRay(&wi_ray, AI_RAY_GLOSSY, &sg->P, NULL, AI_BIG, sg);
		AtFloat count=0.0f;
		while(AiSamplerGetSample(sampit, samples))
		{
			wi_l = iw.sample(samples[0], samples[1]);
			wi = iw.localToWorld(wi_l);
			wi_ray.dir = wi;
			AiTrace(&wi_ray, &scrs);
			result_specularIndirect += scrs.color * iw.eval(wi_l, wo, sg->u, sg->v, tempDiffuse);
			result_diffuseIndirect += scrs.color * tempDiffuse;
														// only doing cosine hemisphere sampling
														// so we can leave this out for now
														// * fabs(wi_l.z) / iw.pdf(wi_l, wo);
			count++;
		}
		result_specularIndirect /= count;
		result_diffuseIndirect /= count;
	}

	sg->out.RGB = result_diffuseDirect + result_diffuseIndirect +
					result_specularDirect + result_specularIndirect;
}


