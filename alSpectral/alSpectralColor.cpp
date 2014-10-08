#include <ai.h>
#include "SPD.h"
#include "alUtil.h"

AI_SHADER_NODE_EXPORT_METHODS(alSpectralColorMtd)

enum alSpectralColorParams
{
	p_mode,
	p_space
};

enum ColorEnum
{
	dark_skin,
	light_skin,
	blue_sky,
	foliage,
	blue_flower,
	bluish_green,
	orange,
	purplish_blue,
	moderate_red,
	purple,
	yellow_green,
	orange_yellow,
	blue,
	green,
	red,
	yellow,
	magenta,
	cyan,
	white_95,
	neutral_8,
	neutral_65,
	neutral_5,
	neutral_35,
	black_2,
	
	checker
};

static const char* ColorNames[] = 
{
	"dark_skin",
	"light_skin",
	"blue_sky",
	"foliage",
	"blue_flower",
	"bluish_green",
	"orange",
	"purplish_blue",
	"moderate_red",
	"purple",
	"yellow_green",
	"orange_yellow",
	"blue",
	"green",
	"red",
	"yellow",
	"magenta",
	"cyan",
	"white_95",
	"neutral_8",
	"neutral_65",
	"neutral_5",
	"neutral_35",
	"black_2",
	"checker",

	NULL
};

enum SpaceEnum
{
	kRec709 = 0,
	kRec2020,
	kACES,
	kACES_D65,
	kP3,
	kP3_D65
};

static const char* SpaceNames[] = 
{
	"rec709",
	"rec2020",
	"aces",
	"aces_D65",
	"p3",
	"p3_D65",
	NULL
};

node_parameters
{
	AiParameterEnum("mode", white_95, ColorNames);
	AiParameterEnum("space", kRec709, SpaceNames);
}

node_initialize
{
	AiNodeSetLocalData(node, NULL);
}

node_finish
{
	SPD* spd = (SPD*)AiNodeGetLocalData(node);
	delete spd;
}

node_update
{
	SPD* spd = (SPD*)AiNodeGetLocalData(node);
	delete spd;


}

shader_evaluate
{
	
	sg->out.RGB = AI_RGB_BLACK;

	SPD* spd = (SPD*)AiShaderGlobalsQuickAlloc(sg, sizeof(SPD));

	int mode = AiShaderEvalParamInt(p_mode);
	int space = AiShaderEvalParamInt(p_space);

	if (mode < checker)
	{
		int c = SPD::kDARK_SKIN + mode;
		spd->set(c);
	}
	else
	{
		// make colorchecker pattern
		float u = fmodf(sg->u, 1.0f) * 6;
		float v = fmodf(sg->v, 1.0f) * 4;

		int x = int(u);
		int y = int(v);

		int cell = y*6+x;
		spd->set(SPD::kDARK_SKIN + cell);

		if (fmodf(u + 0.05f, 1.0f) < 0.1f) spd->set(SPD::kBLACK);
		if (fmodf(v + 0.05f, 1.0f) < 0.1f) spd->set(SPD::kBLACK);
	}

	AiStateSetMsgPtr("als_spd_color", spd);
	if (space == kRec709)
	{
		sg->out.RGB = clamp(CsRec709.xyzToRgb(spd->xyz(SPD::kD65)) * 100.0f, AI_RGB_BLACK, AI_RGB_WHITE);
	}
	else if (space == kRec2020)
	{
		sg->out.RGB = clamp(CsRec2020.xyzToRgb(spd->xyz(SPD::kD65)) * 100.0f, AI_RGB_BLACK, AI_RGB_WHITE);
	}
	else if (space == kACES)
	{
		sg->out.RGB = clamp(CsACES.xyzToRgb(spd->xyz(SPD::kD60)) * 100.0f, AI_RGB_BLACK, AI_RGB_WHITE);
	}
	else if (space == kACES_D65)
	{
		sg->out.RGB = clamp(CsACES_D65.xyzToRgb(spd->xyz(SPD::kD65)) * 100.0f, AI_RGB_BLACK, AI_RGB_WHITE);
	}
	else if (space == kP3)
	{
		sg->out.RGB = clamp(CsP3.xyzToRgb(spd->xyz(SPD::kD65)) * 100.0f, AI_RGB_BLACK, AI_RGB_WHITE);
	}
	else
	{
		sg->out.RGB = clamp(CsP3_D65.xyzToRgb(spd->xyz(SPD::kD65)) * 100.0f, AI_RGB_BLACK, AI_RGB_WHITE);	
	}
}

node_loader
{
	if (i>0) return 0;
    node->methods     = alSpectralColorMtd;
    node->output_type = AI_TYPE_RGB;
    node->name        = "alSpectralColor";
    node->node_type   = AI_NODE_SHADER;
    strcpy(node->version, AI_VERSION);
    return true;
}