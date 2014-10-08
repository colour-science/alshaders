#include <ai.h>
#include <cstring>
#include "SPD.h"
#include "alUtil.h"

AI_SHADER_NODE_EXPORT_METHODS(alSpectralMtd)

enum alSpectralParams
{
	p_mode,
	p_color,
	p_space,
	p_scale,
	p_white_point
};

enum ModeEnum
{
	kCONSTANT = 0,
	kDIFFUSE
};

static const char* ModeNames[] = 
{
	"constant",
	"diffuse",
	NULL
};

enum SpaceEnum
{
	kRec709 = 0,
	kRec2020,
	kACES,
	kP3,
	kXYZ,
	kRGB,
};

static const char* SpaceNames[] = 
{
	"rec709",
	"rec2020",
	"aces",
	"p3",
	"xyz",
	"rgb",
	NULL
};

node_parameters
{
	AiParameterEnum("mode", kDIFFUSE, ModeNames);
	AiParameterRGB("color", 0.0f, 0.0f, 0.0f);
	AiParameterEnum("space", kRec709, SpaceNames);
	AiParameterFlt("scale", 1.0f);
}

struct ShaderData
{
	ShaderData() : diffuse_sampler(NULL) {}
	AtSampler* diffuse_sampler;
};

node_initialize
{
	ShaderData* data = new ShaderData();
	AiNodeSetLocalData(node, data);
}

node_finish
{
	ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);
	AiSamplerDestroy(data->diffuse_sampler);
	delete data;
}

node_update
{
	AtNode *options = AiUniverseGetOptions();
	ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);
	AiSamplerDestroy(data->diffuse_sampler);
	data->diffuse_sampler = AiSampler(AiNodeGetInt(options, "GI_diffuse_samples"), 2);
}

shader_evaluate
{
	ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);
	int mode = AiShaderEvalParamInt(p_mode);
	
	int space = AiShaderEvalParamInt(p_space);
	float scale = AiShaderEvalParamFlt(p_scale);

	AtRGB rgb_color = AiShaderEvalParamRGB(p_color);
	SPD* spd_color = NULL;
	AiStateGetMsgPtr("als_spd_color", (void**)&spd_color);

	SPD spd_result(SPD::kBLACK);
	SPD spd_indirect(SPD::kBLACK);
	SPD* spd_accum = NULL;
	if (sg->Rt & AI_RAY_CAMERA)
	{
		spd_accum = (SPD*)AiShaderGlobalsQuickAlloc(sg, sizeof(SPD));
		spd_accum->set(SPD::kBLACK);
		AiStateSetMsgPtr("als_spd_accum", spd_accum);
	}
	else
	{
		AiStateGetMsgPtr("als_spd_accum",(void**)&spd_accum);
	}

	if (!spd_accum)
	{
		// something went badly wrong. 
		// set cyan for error and exit
		sg->out.RGB = rgb(0, 1, 1);
		return;
	}

	// either we had nothing connected, or something else has gone wrong upstream
	// set magenta to signify error
	if (!spd_color)
	{
		sg->out.RGB = rgb(1.0f, 0.0f, 1.0f);
		return;
	}

	AtRGB rgb_result = AI_RGB_BLACK;
	
	if (mode == kCONSTANT)
	{
		// constant shading, just assign the color
		spd_result = *spd_color;
	}
	else
	{
		void* mis = AiOrenNayarMISCreateData(sg, 0.0f);

		// diffuse shading
		SPD spd_d65(SPD::kD65);
		AiLightsPrepare(sg);
		while(AiLightsGetSample(sg))
        {
        	rgb_result += sg->Li * sg->we * MAX(AiV3Dot(sg->Nf, sg->Ld), 0.0f) * rgb_color;
        	float intensity = sg->Li.g;
        	float weight = intensity * sg->we * MAX(AiV3Dot(sg->Nf, sg->Ld), 0.0f);
        	spd_result = spd_result + (*spd_color * spd_d65 * weight) * scale;
        }

        // indirect diffuse
        // build a local frame for sampling
	    AtVector U, V;
	    if (!AiV3isZero(sg->dPdu))
	    {
	        // we have valid a valid dPdu derivative, construct V 
	        AtVector Utmp = AiV3Normalize(sg->dPdu);
	        V = AiV3Normalize(AiV3Cross(sg->Nf, Utmp));
	        U = AiV3Cross(V, sg->Nf);
	    }
	    else
	    {
	        AiBuildLocalFramePolar(&U, &V, &sg->Nf);
	    }
        AtSamplerIterator* sampit = AiSamplerIterator(data->diffuse_sampler, sg);
        AtRay wi_ray;
        AtScrSample scrs;
        AiMakeRay(&wi_ray, AI_RAY_DIFFUSE, &sg->P, NULL, AI_BIG, sg);
        float samples[2];
        AtVector wi;
        AtRGB rgb_indirect = AI_RGB_BLACK;
        while (AiSamplerGetSample(sampit, samples))
        {
            // cosine hemisphere sampling as O-N sampling does not work outside of a light loop
            float stheta = sqrtf(float(samples[0]));
            float phi = float(AI_PITIMES2 * samples[1]);
            wi.x = stheta * cosf(phi);
            wi.y = stheta * sinf(phi);
            wi.z = sqrtf(1.0f - float(samples[0]));
            AiV3RotateToFrame(wi, U, V, sg->Nf);

            float cos_theta = AiV3Dot(wi, sg->Nf);
            if (cos_theta <= 0.0f) continue;

            float p = cos_theta * float(AI_ONEOVERPI);
            
            // trace the ray
            wi_ray.dir = wi;
            AtRGB f = AiOrenNayarMISBRDF(mis, &wi) / p;
            if (AiTrace(&wi_ray, &scrs))
            {
            	spd_indirect = spd_indirect + (*spd_accum) * f.g;
            	rgb_indirect += f * scrs.color;
            }
        }
        spd_indirect = spd_indirect * (*spd_color) * AiSamplerGetSampleInvCount(sampit);
        spd_result = spd_result + spd_indirect;
        rgb_indirect *= rgb_color * AiSamplerGetSampleInvCount(sampit);
        rgb_result += rgb_indirect;
	}

	if (sg->Rt & AI_RAY_CAMERA)
	{
		if (space == kRec709)
		{
			sg->out.RGB = CsRec709.xyzToRgb(spd_result.xyz(SPD::kWHITE));
		}
		else if (space == kRec2020)
		{
			sg->out.RGB = CsRec2020.xyzToRgb(spd_result.xyz(SPD::kWHITE));
		}
		else if (space == kACES)
		{
			sg->out.RGB = CsACES.xyzToRgb(spd_result.xyz(SPD::kWHITE));
		}
		else if (space == kP3)
		{
			sg->out.RGB = CsP3_D65.xyzToRgb(spd_result.xyz(SPD::kWHITE));
		}
		else if (space == kXYZ)
		{
			sg->out.RGB = spd_result.xyz(SPD::kWHITE);
		}
		else
		{
			sg->out.RGB = rgb_result;
		}

	}
	else
	{
		*spd_accum = spd_result;
		sg->out.RGB = rgb_result;
	}
}

node_loader
{
	if (i>0) return 0;
    node->methods     = alSpectralMtd;
    node->output_type = AI_TYPE_RGB;
    node->name        = "alSpectral";
    node->node_type   = AI_NODE_SHADER;
    strcpy(node->version, AI_VERSION);
    return true;
}