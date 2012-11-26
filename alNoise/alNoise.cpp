#include "Remap.h"
#include <ai.h>

AI_SHADER_NODE_EXPORT_METHODS(alNoise)

enum NoiseSpaceEnum
{
	NS_WORLD = 0,
	NS_OBJECT,
	NS_PREF,
	NS_UV
};

static const char* noiseSpaceNames[] = {
		"world",
		"object",
		"Pref",
		"UV",
		NULL
};

enum alNoiseParams
{
	p_space,
	p_time,
	p_octaves,
	p_distortion,
	p_lacunarity,
	REMAP_FLOAT_PARAM_ENUM,
	p_color1,
	p_color2,
	p_P
};

node_parameters
{
	AiParameterENUM("space", 0, noiseSpaceNames);
	AiParameterFLT("time", 0.0f);
	AiParameterINT("octaves", 8);
	AiParameterFLT("distortion", 0.5f);
	AiParameterFLT("lacunarity", 2.0f);
	REMAP_FLOAT_PARAM_DECLARE;
	AiParameterRGB("color1", 0.0f, 0.0f, 0.0f);
	AiParameterRGB("color2", 1.0f, 1.0f, 1.0f);
	AiParameterPnt("P", 0.0f, 0.0f, 0.0f);
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alNoise;
   node->output_type = AI_TYPE_RGB;
   node->name        = "alNoise";
   node->node_type   = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return TRUE;
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
	int space = AiShaderEvalParamInt(p_space);
	int time = AiShaderEvalParamFlt(p_time);
	int octaves = AiShaderEvalParamInt(p_octaves);
	AtFloat distortion = AiShaderEvalParamFlt(p_distortion);
	AtFloat lacunarity = AiShaderEvalParamFlt(p_lacunarity);
	AtRGB color1 = AiShaderEvalParamRGB(p_color1);
	AtRGB color2 = AiShaderEvalParamRGB(p_color2);
	AtPoint Pin = AiShaderEvalParamPnt(p_P);

	// choose what space we want to calculate in
	AtPoint P;
	if (AiNodeIsLinked(node, "P"))
	{
		P = Pin;
	}
	else
	{
		switch (space)
		{
		case NS_OBJECT:
			P = sg->Po;
			break;
		case NS_UV:
			P.x = sg->u;
			P.y = sg->v;
			P.z = 0.0f;
			break;
		case NS_PREF:
			if (!AiUDataGetPnt("Pref", &P))
				P = sg->Po;
			break;
		default:
			P = sg->P;
			break;
		}
	}

	float n = AiNoise4(P, time, octaves, distortion, lacunarity);

	RemapFloat r = REMAP_FLOAT_CREATE;
	n = r.remap(n);

	sg->out.RGB = AiColorLerp(n, color1, color2);
}


