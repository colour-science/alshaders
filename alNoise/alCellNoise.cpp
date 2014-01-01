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
	p_smoothChips,
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
	AiParameterBOOL("smoothChips", false);
	AiParameterRGB("chipColor1", .383f, .318f, .252f);
	AiParameterFLT("chipProb1", 0.2f);
	AiParameterRGB("chipColor2", .383f, .191f, 0.01f);
	AiParameterFLT("chipProb2", 0.4f);
	AiParameterRGB("chipColor3", .635f, .612f, .563f);
	AiParameterFLT("chipProb3", 0.6f);
	AiParameterRGB("chipColor4", .509f, .361f, .213f);
	AiParameterFLT("chipProb4", 0.8f);
	AiParameterRGB("chipColor5", .593f, .472f, .248f);
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
   return true;
}


struct ShaderData
{
	int space;
	float frequency;
	int mode;
	int octaves;
	float randomness;
	float lacunarity;
	float f1w;
	float f2w;
	float f3w;
	float f4w;
	float mynkowskiShape;
	bool smoothChips;
	float chipProb1;
	float chipProb2;
	float chipProb3;
	float chipProb4;
	float chipProb5;
};

node_initialize
{
	ShaderData* data = new ShaderData;
	AiNodeSetLocalData(node, data);
}

node_finish
{
	ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);
	delete data;
}

node_update
{
	ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);
	data->space = params[p_space].INT;
	data->frequency = params[p_frequency].FLT;
	data->mode = params[p_mode].INT;
	data->octaves = params[p_octaves].INT;
	data->randomness= params[p_randomness].FLT;
	data->lacunarity = params[p_lacunarity].FLT;
	data->f1w = params[p_f1w].FLT;
	data->f2w = params[p_f2w].FLT;
	data->f3w = params[p_f3w].FLT;
	data->f4w = params[p_f4w].FLT;
	data->mynkowskiShape = params[p_mynkowskiShape].FLT;
	data->smoothChips = params[p_smoothChips].BOOL;
	data->chipProb1 = params[p_chipProb1].FLT;
	data->chipProb2 = params[p_chipProb2].FLT;
	data->chipProb3 = params[p_chipProb3].FLT;
	data->chipProb4 = params[p_chipProb4].FLT;
	data->chipProb5 = params[p_chipProb5].FLT;
}

shader_evaluate
{
	ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);
	AtRGB color1 = AiShaderEvalParamRGB(p_color1);
	AtRGB color2 = AiShaderEvalParamRGB(p_color2);

	AtRGB chipColor1 = AiShaderEvalParamRGB(p_chipColor1);
	AtRGB chipColor2 = AiShaderEvalParamRGB(p_chipColor2);
	AtRGB chipColor3 = AiShaderEvalParamRGB(p_chipColor3);
	AtRGB chipColor4 = AiShaderEvalParamRGB(p_chipColor4);
	AtRGB chipColor5 = AiShaderEvalParamRGB(p_chipColor5);

	AtPoint Pin = AiShaderEvalParamPnt(p_P);

	// choose what space we want to calculate in
	AtPoint P;
	if (AiNodeIsLinked(node, "P"))
	{
		P = Pin;
	}
	else
	{
		switch (data->space)
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
	P *= data->frequency;

	// get cellular result
	float F[4];
	AtVector delta[4];
	AtUInt32 ID[4];
	AiCellular(P, 4, data->octaves, data->lacunarity, data->randomness, F, delta, ID);

	if (data->mode == CN_FEATURES)
	{
		// check what distance metric we're using
		if (data->mynkowskiShape != 2.0f)
		{
			// Use Mynkowski distance metric instead of the distances computed by Arnold
			float ims = 1.0f / data->mynkowskiShape;
			for (int i=0; i < 4; ++i)
			{
				F[i] = pow(fabs(delta[i].x), data->mynkowskiShape) + pow(fabs(delta[i].y), data->mynkowskiShape) 
						+ pow(fabs(delta[i].z), data->mynkowskiShape);
				F[i] = pow(F[i], ims);
			}
		}

		// weight the feature distances
		float n = F[0]*data->f1w + F[1]*data->f2w + F[2]*data->f3w + F[3]*data->f4w;

		// normalize for the number of octaves
		n /= float(data->octaves);

		RemapFloat r = REMAP_FLOAT_CREATE;
		n = r.remap(n);

		sg->out.RGB = AiColorLerp(n, color1, color2);
	}
	else if (data->mode == CN_CHIPS)
	{
		double rr = random(ID[0]);

		AtRGB chipColor;
		if (rr < data->chipProb1)
		{
			chipColor = chipColor1;
		}
		else if (rr < data->chipProb2)
		{
			if (data->smoothChips) chipColor = lerp(chipColor1, chipColor2, (rr-data->chipProb1)/(data->chipProb2-data->chipProb1));
			else chipColor = chipColor2;
		}
		else if (rr < data->chipProb3)
		{
			if (data->smoothChips) chipColor = lerp(chipColor2, chipColor3, (rr-data->chipProb2)/(data->chipProb3-data->chipProb2));
			else chipColor = chipColor3;
		}
		else if (rr < data->chipProb4)
		{
			if (data->smoothChips) chipColor = lerp(chipColor3, chipColor4, (rr-data->chipProb3)/(data->chipProb4-data->chipProb3));
			else chipColor = chipColor4;
		}
		else
		{
			if (data->smoothChips) chipColor = lerp(chipColor4, chipColor5, (rr-data->chipProb4)/(data->chipProb5-data->chipProb4));
			else chipColor = chipColor5;
		}

		sg->out.RGB = chipColor;
	}
}


