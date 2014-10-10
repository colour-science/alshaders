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
	p_roughness,
	p_inputspace
};

enum ModeEnum
{
	kCONSTANT = 0,
	kDIFFUSE,
	kGOLD,
	kSILVER,
	kCOPPER,
	kTUNGSTEN
};

static const char* ModeNames[] = 
{
	"constant",
	"diffuse",
	"gold",
	"silver",
	"copper",
	"tungsten",
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
	AiParameterFlt("roughness", 0.0f);
	AiParameterEnum("inputspace", kRec709, SpaceNames);
}

struct ShaderData
{
	ShaderData() : diffuse_sampler(NULL), glossy_sampler(NULL)
	 {}
	AtSampler* diffuse_sampler;
	AtSampler* glossy_sampler;
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
	AiSamplerDestroy(data->glossy_sampler);
	delete data;
}

node_update
{
	AtNode *options = AiUniverseGetOptions();
	ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);
	AiSamplerDestroy(data->diffuse_sampler);
	AiSamplerDestroy(data->glossy_sampler);
	data->diffuse_sampler = AiSampler(AiNodeGetInt(options, "GI_diffuse_samples"), 2);
	data->glossy_sampler = AiSampler(AiNodeGetInt(options, "GI_glossy_samples"), 2);
}

shader_evaluate
{
	ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);
	int mode = AiShaderEvalParamInt(p_mode);
	
	int space = AiShaderEvalParamInt(p_space);
	int inputspace = AiShaderEvalParamInt(p_inputspace);
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

	bool do_glossy = !(sg->Rt & AI_RAY_DIFFUSE);

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

	if (mode == kCONSTANT)
	{
		// constant shading, just assign the color
		spd_result = *spd_color;
	}
	else if (mode == kDIFFUSE)
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
	else if (do_glossy)
	{
		// glossy reflection
		float roughness = AiShaderEvalParamFlt(p_roughness);
		void* mis = AiCookTorranceMISCreateData(sg, &U, &V, roughness, roughness);

		ColorSpace* cs;
		if (inputspace == kRec2020)
		{
			cs = &CsRec2020;
		}
		else if (inputspace == kACES)
		{
			cs = &CsACES_D65;
		}
		else if (inputspace == kP3)
		{
			cs = &CsP3_D65;
		}
		else
		{
			cs = &CsRec709;
		}

		SPD spd_d65(SPD::kD65);
		SPD kr;
		AtVector H, wi, wo;
		wo = -sg->Rd;
		AiLightsPrepare(sg);
		while(AiLightsGetSample(sg))
        {
        	wi = sg->Ld;
        	AtRGB f = AiCookTorranceMISBRDF(mis, &wi);
        	
        	H = AiV3Normalize(wi+wo);
        	float intensity = sg->Li.g;
        	
        	float cos_theta = MAX(AiV3Dot(H, wi), 0.0f);
        	SPD::fresnel(cos_theta, mode - kGOLD, kr);
        	AtRGB kr_rgb = cs->xyzToRgb(kr.xyz(SPD::kD65)) * 100.0f;
        	float weight = intensity * sg->we;


        	spd_result = spd_result + (spd_d65 * weight) * kr * f[1] * scale;

        	rgb_result += sg->Li * sg->we * f * kr_rgb;
        }

        // glossy indirect
        AtSamplerIterator* sampit = AiSamplerIterator(data->glossy_sampler, sg);
        AtRay wi_ray;
        AtScrSample scrs;
        AiMakeRay(&wi_ray, AI_RAY_GLOSSY, &sg->P, NULL, AI_BIG, sg);
        float samples[2];
        AtRGB rgb_indirect = AI_RGB_BLACK;
        while (AiSamplerGetSample(sampit, samples))
        {
            wi = AiCookTorranceMISSample(mis, samples[0], samples[1]);
            AtRGB f = AiCookTorranceMISBRDF(mis, &wi);
            float pdf = AiCookTorranceMISPDF(mis, &wi);

            H = AiV3Normalize(wi+wo);
            float cos_theta = MAX(AiV3Dot(H, wi), 0.0f);
        	SPD::fresnel(cos_theta, mode - kGOLD, kr);

        	AtRGB kr_rgb = cs->xyzToRgb(kr.xyz(SPD::kD65)) * 100.0f;

            // trace the ray
            wi_ray.dir = wi;
            if (AiTrace(&wi_ray, &scrs))
            {
            	spd_indirect = spd_indirect + (*spd_accum) * (f.g/pdf) * kr;
            	rgb_indirect += f/pdf * scrs.color * kr_rgb;
            }
        }
        spd_indirect = spd_indirect * AiSamplerGetSampleInvCount(sampit);
        spd_result = spd_result + spd_indirect;
        rgb_indirect *= AiSamplerGetSampleInvCount(sampit);
        rgb_result += rgb_indirect;
	}

	if (sg->Rt & AI_RAY_CAMERA)
	{
		SPD::Constant c;
		if (mode == kCONSTANT)
		{
			c = SPD::kD65;
		}
		else
		{
			c = SPD::kWHITE;
		}

		if (space == kRec709)
		{
			sg->out.RGB = CsRec709.xyzToRgb(spd_result.xyz(c));
		}
		else if (space == kRec2020)
		{
			sg->out.RGB = CsRec2020.xyzToRgb(spd_result.xyz(c));
		}
		else if (space == kACES)
		{
			sg->out.RGB = CsACES_D65.xyzToRgb(spd_result.xyz(c));
		}
		else if (space == kP3)
		{
			sg->out.RGB = CsP3_D65.xyzToRgb(spd_result.xyz(c));
		}
		else if (space == kXYZ)
		{
			sg->out.RGB = spd_result.xyz(c);
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