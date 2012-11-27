#include "Remap.h"
#include <ai.h>

AI_SHADER_NODE_EXPORT_METHODS(alFractal)

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
	p_frequency,
	p_time,
	p_octaves,
	p_distortion,
	p_lacunarity,
	p_gain,
	p_turbulent,
	p_ridged,
	p_ridgeOffset,
	REMAP_FLOAT_PARAM_ENUM,
	p_color1,
	p_color2,
	p_P
};

node_parameters
{
	AiParameterENUM("space", 0, noiseSpaceNames);
	AiParameterFLT("frequency", 1.0f);
	AiParameterFLT("time", 0.0f);
	AiParameterINT("octaves", 8);
	AiParameterFLT("distortion", 0.0f);
	AiParameterFLT("lacunarity", 2.121f);
	AiParameterFLT("gain", 0.5f);
	AiParameterBOOL("turbulent", false);
	AiParameterBOOL("ridged", false);
	AiParameterFLT("ridgeOffset", 0.0f);
	REMAP_FLOAT_PARAM_DECLARE;
	AiParameterRGB("color1", 0.0f, 0.0f, 0.0f);
	AiParameterRGB("color2", 1.0f, 1.0f, 1.0f);
	AiParameterPnt("P", 0.0f, 0.0f, 0.0f);
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alFractal;
   node->output_type = AI_TYPE_RGB;
   node->name        = "alFractal";
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
	AtFloat frequency = AiShaderEvalParamFlt(p_frequency);
	AtFloat time = AiShaderEvalParamFlt(p_time);
	int octaves = AiShaderEvalParamInt(p_octaves);
	AtFloat distortion = AiShaderEvalParamFlt(p_distortion);
	AtFloat lacunarity = AiShaderEvalParamFlt(p_lacunarity);
	AtFloat gain = AiShaderEvalParamFlt(p_gain);
	bool turbulent = AiShaderEvalParamBool(p_turbulent);
	bool ridged = AiShaderEvalParamBool(p_ridged);
	AtFloat ridgeOffset = AiShaderEvalParamFlt(p_ridgeOffset);
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

	P *= frequency;

	AtFloat n = 0.0f;
	AtFloat amp = 1.0f;
	AtFloat weight = 1;
	float v;
	for (int i=0; i < octaves; ++i)
	{
		AtPoint PP = P;
		if (distortion != 0.0f)
			PP += distortion * AiVNoise3(P, 1, 0, 0);
		v = AiPerlin3(PP);
		if (turbulent) v = fabsf(v);
		if (ridged)
		{
			v = ridgeOffset - v;
			v *= v;
			v *= weight;
			weight = v * 2;
			weight = clamp(weight, 0.0f, 1.0f);
		}
		n += v * amp;
		amp *= gain;

		P *= lacunarity;
	}

	RemapFloat r = REMAP_FLOAT_CREATE;
	n = r.remap(n);

	sg->out.RGB = AiColorLerp(n, color1, color2);
}


