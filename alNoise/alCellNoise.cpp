#include "Remap.h"
#include <ai.h>
#include <cmath>
#include <cstring>

AI_SHADER_NODE_EXPORT_METHODS(alCellNoise)

enum CellNoiseSpaceEnum
{
	NS_WORLD = 0,
	NS_OBJECT,
	NS_PREF,
	NS_UV
};

static const char* cellnoiseSpaceNames[] =
{
	"world",
	"object",
	"Pref",
	"UV",
	NULL
};

enum CellNoiseModeEnum
{
	CN_FEATURES = 0,
	CN_CHIPS
};

static const char* cellNoiseModeNames[] =
{
	"features",
	"chips",
	NULL
};

enum alCellNoiseParams
{
	p_space,
	p_frequency,
	p_mode,
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
	p_chipColor1,
	p_chipProb1,
	p_chipColor2,
	p_chipProb2,
	p_chipColor3,
	p_chipProb3,
	p_chipColor4,
	p_chipProb4,
	p_chipColor5,
	p_chipProb5,
	p_P,
	REMAP_FLOAT_PARAM_ENUM
};

node_parameters
{
	AiParameterENUM("space", 0, cellnoiseSpaceNames);
	AiParameterFLT("frequency", 1.0f);
	AiParameterENUM("mode", 0, cellNoiseModeNames);
	AiParameterINT("octaves", 1);
	AiParameterFLT("randomness", 1.0f);
	AiParameterFLT("lacunarity", 2.121f);
	AiParameterFlt("f1w", -1.0f);
	AiParameterFlt("f2w", 1.0f);
	AiParameterFlt("f3w", 0.0f);
	AiParameterFlt("f4w", 0.0f);
	AiParameterFLT("mynkowskiShape", 2.0f);
	AiParameterRGB("color1", 0.0f, 0.0f, 0.0f);
	AiParameterRGB("color2", 1.0f, 1.0f, 1.0f);
	AiParameterRGB("chipColor1", 0.0f, 0.0f, 0.0f);
	AiParameterFLT("chipProb1", 0.2f);
	AiParameterRGB("chipColor2", 0.0f, 0.0f, 0.0f);
	AiParameterFLT("chipProb2", 0.4f);
	AiParameterRGB("chipColor3", 0.0f, 0.0f, 0.0f);
	AiParameterFLT("chipProb3", 0.6f);
	AiParameterRGB("chipColor4", 0.0f, 0.0f, 0.0f);
	AiParameterFLT("chipProb4", 0.8f);
	AiParameterRGB("chipColor5", 0.0f, 0.0f, 0.0f);
	AiParameterFLT("chipProb5", 1.0f);
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
   ::strcpy(node->version, AI_VERSION);
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
	int mode = AiShaderEvalParamInt(p_mode);
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

	AtRGB chipColor1 = AiShaderEvalParamRGB(p_chipColor1);
	AtFloat chipProb1 = AiShaderEvalParamFlt(p_chipProb1);
	AtRGB chipColor2 = AiShaderEvalParamRGB(p_chipColor2);
	AtFloat chipProb2 = AiShaderEvalParamFlt(p_chipProb2);
	AtRGB chipColor3 = AiShaderEvalParamRGB(p_chipColor3);
	AtFloat chipProb3 = AiShaderEvalParamFlt(p_chipProb3);
	AtRGB chipColor4 = AiShaderEvalParamRGB(p_chipColor4);
	AtFloat chipProb4 = AiShaderEvalParamFlt(p_chipProb4);
	AtRGB chipColor5 = AiShaderEvalParamRGB(p_chipColor5);
	AtFloat chipProb5 = AiShaderEvalParamFlt(p_chipProb5);


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

	if (mode == CN_FEATURES)
	{
		// check what distance metric we're using
		if (ms != 2.0f)
		{
			// Use Mynkowski distance metric instead of the distances computed by Arnold
			AtFloat ims = 1.0f / ms;
			for (int i=0; i < 4; ++i)
			{
				F[i] = pow(fabs(delta[i].x), ms) + pow(fabs(delta[i].y), ms) + pow(fabs(delta[i].z), ms);
				F[i] = pow(F[i], ims);
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
	else if (mode == CN_CHIPS)
	{
		double rr = random(ID[0]);

		AtRGB chipColor;
		if (rr < chipProb1)
		{
			chipColor = chipColor1;
		}
		else if (rr < chipProb2)
		{
			chipColor = chipColor2;
		}
		else if (rr < chipProb3)
		{
			chipColor = chipColor3;
		}
		else if (rr < chipProb4)
		{
			chipColor = chipColor4;
		}
		else
		{
			chipColor = chipColor5;
		}

		sg->out.RGB = chipColor;
	}
}


