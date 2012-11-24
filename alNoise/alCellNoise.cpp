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

enum CellNoiseModeEnum
{
	CN_F1 = 0,
	CN_F2,
	CN_F3,
	CN_F4,
	CN_F2_MINUS_F1,
	CN_F1_PLUS_F2,
	CN_F1_TIMES_F2
};

static const char* cellnoiseModeNames[] = {
	"F1",
	"F2",
	"F3",
	"F4",
	"F2 - F1",
	"F1 + F2",
	"F1 * F2",
	NULL
};

enum alNoiseParams
{
	p_space,
	p_frequency,
	p_octaves,
	p_randomness,
	p_lacunarity,
	p_mode,
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
	AiParameterENUM("mode", 0, cellnoiseModeNames);
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
	int mode = AiShaderEvalParamInt(p_mode);
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
		default:
			P = sg->P;
			break;
		}
	}
	// scale the space
	P *= frequency;

	// get cellular result
	AtFloat F[5];
	AtVector delta[5];
	AtUInt32 ID[5];
	AiCellular(P, 5, octaves, lacunarity, randomness, F, delta, ID);

	// check what distance metric we're using
	if (ms != 2.0f)
	{
		// Use Mynkowski distance metric instead of the distances computed by Arnold
		AtFloat ims = 1.0f / ms;
		for (int i=0; i < 5; ++i)
		{
			F[i] = powf(fabsf(delta[i].x), ms) + powf(fabsf(delta[i].y), ms) + powf(fabsf(delta[i].z), ms);
			F[i] = powf(F[i], ims);
		}
	}

	// choose what to do with the feature distances
	float n;
	switch (mode)
	{
	case CN_F2:
		n = F[1] * 0.5f;
		break;
	case CN_F3:
		n = F[2] * 0.5f;
			break;
	case CN_F4:
		n = F[3] * 0.5f;
			break;
	case CN_F2_MINUS_F1:
		n = (F[1] - F[0])*0.5f;
		break;
	case CN_F1_PLUS_F2:
		n = (F[0] + F[1]) * 0.25f;
		break;
	case CN_F1_TIMES_F2:
		n = (F[0] * F[1]) * 0.25f;
		break;
	default:
		n = F[0] * 0.5f;
		break;
	}

	// normalize for the number of octaves
	n /= float(octaves);

	RemapFloat r = REMAP_FLOAT_CREATE;
	n = r.remap(n);

	sg->out.RGB = AiColorLerp(n, color1, color2);
}


