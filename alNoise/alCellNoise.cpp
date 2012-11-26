#include "Remap.h"
#include <ai.h>

AI_SHADER_NODE_EXPORT_METHODS(alCellNoise)

enum CellNoiseSpaceEnum
{
	NS_WORLD = 0,
	NS_OBJECT,
	NS_PREF,
	NS_UV
};

static const char* cellnoiseSpaceNames[] = {
	"world",
	"object",
	"Pref",
	"UV",
	NULL
};

enum alCellNoiseParams
{
	p_space,
	p_frequency,
	p_octaves,
	p_randomness,
	p_lacunarity,
	p_f1w,
	p_f2w,
	p_f3w,
	p_f4w,
	p_mynkowskiShape,
	p_color1,
	p_color2,
	p_P,
	REMAP_FLOAT_PARAM_ENUM
};

node_parameters
{
	AiParameterENUM("space", 0, cellnoiseSpaceNames);
	AiParameterFLT("frequency", 1.0f);
	AiParameterINT("octaves", 1);
	AiParameterFLT("randomness", 1.0f);
	AiParameterFLT("lacunarity", 1.92f);
	AiParameterFlt("f1w", -1.0f);
	AiParameterFlt("f2w", 1.0f);
	AiParameterFlt("f3w", 0.0f);
	AiParameterFlt("f4w", 0.0f);
	AiParameterFLT("mynkowskiShape", 2.0f);
	AiParameterRGB("color1", 0.0f, 0.0f, 0.0f);
	AiParameterRGB("color2", 1.0f, 1.0f, 1.0f);
	AiParameterPnt("P", 0.0f, 0.0f, 0.0f);
	REMAP_FLOAT_PARAM_DECLARE;
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alCellNoise;
   node->output_type = AI_TYPE_RGB;
   node->name        = "alCellNoise";
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
	int octaves = AiShaderEvalParamInt(p_octaves);
	AtFloat randomness = AiShaderEvalParamFlt(p_randomness);
	AtFloat lacunarity = AiShaderEvalParamFlt(p_lacunarity);
	AtFloat f1w = AiShaderEvalParamFlt(p_f1w);
	AtFloat f2w = AiShaderEvalParamFlt(p_f2w);
	AtFloat f3w = AiShaderEvalParamFlt(p_f3w);
	AtFloat f4w = AiShaderEvalParamFlt(p_f4w);
	AtFloat ms = AiShaderEvalParamFlt(p_mynkowskiShape);
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
	// scale the space
	P *= frequency;

	// get cellular result
	AtFloat F[4];
	AtVector delta[4];
	AtUInt32 ID[4];
	AiCellular(P, 4, octaves, lacunarity, randomness, F, delta, ID);

	// check what distance metric we're using
	if (ms != 2.0f)
	{
		// Use Mynkowski distance metric instead of the distances computed by Arnold
		AtFloat ims = 1.0f / ms;
		for (int i=0; i < 4; ++i)
		{
			F[i] = powf(fabsf(delta[i].x), ms) + powf(fabsf(delta[i].y), ms) + powf(fabsf(delta[i].z), ms);
			F[i] = powf(F[i], ims);
		}
	}

	// weight the feature distances
	float n = F[0]*f1w + F[1]*f2w + F[2]*f3w + F[3]*f4w;

	// normalize for the number of octaves
	n /= float(octaves);

	RemapFloat r = REMAP_FLOAT_CREATE;
	n = r.remap(n);

	sg->out.RGB = AiColorLerp(n, color1, color2);
}


