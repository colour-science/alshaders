#include <ai.h>
#include <cstring>
#include <iostream>
#include <ai_sampler.h>
#include <OpenEXR/ImathVec.h>
#include <OpenEXR/ImathMatrix.h>
#include <OpenEXR/ImathMatrixAlgo.h>

#include "alUtil.h"
#include "MIS.h"
#include "BeckmannMicrofacet.h"
#include "alSurface.h"
#include "Shadows.h"

AI_SHADER_NODE_EXPORT_METHODS(alSurfaceMtd)

#define GlossyMISBRDF AiCookTorranceMISBRDF
#define GlossyMISPDF AiCookTorranceMISPDF
#define GlossyMISSample AiCookTorranceMISSample

#define GlossyMISBRDF_wrap AiCookTorranceMISBRDF_wrap
#define GlossyMISPDF_wrap AiCookTorranceMISPDF_wrap
#define GlossyMISSample_wrap AiCookTorranceMISSample_wrap

#define GlossyMISCreateData AiCookTorranceMISCreateData

enum alSurfaceParams
{
	// diffuse
	p_diffuseScale=0,
	p_diffuseColor,

	// sss
	p_sssMix,
	p_sssRadius,
	p_sssRadiusColor,
	p_sssScale,

	p_ssScale,
	p_ssRadius,
	p_ssRadiusColor,

	// specular
	p_specular1Scale,
	p_specular1Color,
	p_specular1Roughness,
	p_specular1Ior,
	p_specular2Scale,
	p_specular2Color,
	p_specular2Roughness,
	p_specular2Ior,

	// transmission
	p_transmissionScale,
	p_transmissionColor,
	p_transmissionLinkToSpecular1,
	p_transmissionRoughness,
	p_transmissionIor,
	p_absorptionEnable,
	p_absorptionDensity,
	p_absorptionColor,

	p_bump
};

node_parameters
{
	AiParameterFLT( "diffuseScale", 1.0f );
	AiParameterRGB( "diffuseColor", 0.18f, 0.18f, 0.18f );

	AiParameterFLT( "sssMix", 0.0f );
	AiParameterFLT( "sssRadius", 3.6f );
	AiParameterRGB( "sssRadiusColor", .439f, .156f, .078f );
	AiMetaDataSetBool(mds, "sssRadiusColor", "always_linear", true);  // no inverse-gamma correction
	AiParameterFLT( "sssScale", 10.0f );

	AiParameterFLT( "ssScale", 0.0f );
	AiParameterFLT( "ssRadius", 3.6f );
	AiParameterRGB( "ssRadiusColor", .439f, .156f, .078f );
	AiMetaDataSetBool(mds, "ssRadiusColor", "always_linear", true);  // no inverse-gamma correction

	AiParameterFLT( "specular1Scale", 1.0f );
	AiParameterRGB( "specular1Color", 1.0f, 1.0f, 1.0f );
	AiParameterFLT( "specular1Roughness", 0.3f );
	AiParameterFLT( "specular1Ior", 1.4f );

	AiParameterFLT( "specular2Scale", 1.0f );
	AiParameterRGB( "specular2Color", 1.0f, 1.0f, 1.0f );
	AiParameterFLT( "specular2Roughness", 0.3f );
	AiParameterFLT( "specular2Ior", 1.4f );

	AiParameterFLT("transmissionScale", 0.0f );
	AiParameterRGB("transmissionColor", 1.0f, 1.0f, 1.0f );
	AiParameterBOOL("transmissionLinkToSpecular1", true);
	AiParameterFLT("transmissionRoughness", 0.1f );
	AiParameterFLT("transmissionIor", 1.4f );
	AiParameterBOOL("absorptionEnable", false);
	AiParameterFLT("absorptionDensity", 1.0f);
	AiParameterRGB("absorptionColor", 1.0f, 1.0f, 1.0f);
}


node_loader
{
	if (i>0) return 0;
	node->methods     = alSurfaceMtd;
	node->output_type = AI_TYPE_RGB;
	node->name        = "alSurface";
	node->node_type   = AI_NODE_SHADER;
	strcpy(node->version, AI_VERSION);
	return TRUE;
}

node_initialize
{
	ShaderData *data = (ShaderData*) AiMalloc(sizeof(ShaderData));
	AiNodeSetLocalData(node,data);
	data->diffuse_sampler = NULL;
	data->glossy_sampler = NULL;
	data->refraction_sampler = NULL;
};

node_finish
{
	if (AiNodeGetLocalData(node))
	{
		ShaderData* data = (ShaderData*) AiNodeGetLocalData(node);

		AiSamplerDestroy(data->diffuse_sampler);
		AiSamplerDestroy(data->glossy_sampler);
		AiSamplerDestroy(data->refraction_sampler);

		AiFree((void*) data);
		AiNodeSetLocalData(node, NULL);
		}
}


node_update
{
	ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);
	AtNode *options   = AiUniverseGetOptions();
	data->GI_diffuse_depth = AiNodeGetInt(options, "GI_diffuse_depth");
	data->GI_reflection_depth = AiNodeGetInt(options, "GI_reflection_depth");
	data->GI_refraction_depth = AiNodeGetInt(options, "GI_refraction_depth");
	data->GI_glossy_depth = AiNodeGetInt(options, "GI_glossy_depth");
	data->GI_glossy_samples = AiNodeGetInt(options, "GI_glossy_samples");
	data->GI_diffuse_samples = AiNodeGetInt(options, "GI_diffuse_samples");

	// setup samples
	AiSamplerDestroy(data->diffuse_sampler);
	AiSamplerDestroy(data->glossy_sampler);
	AiSamplerDestroy(data->refraction_sampler);
	data->diffuse_sampler = AiSampler(data->GI_diffuse_samples, 2);
	data->glossy_sampler = AiSampler(data->GI_glossy_samples, 2);
	data->refraction_sampler = AiSampler(AiNodeGetInt(options, "GI_refraction_samples"), 2);
};


shader_evaluate
{
	ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);


	// Evaluate bump;
	AtVector N_orig;
	AtVector Nf_orig;
	AtRGB bump = AiShaderEvalParamRGB( p_bump );

	// Initialize parameter temporaries
	AtRGB diffuseColor = AiShaderEvalParamRGB( p_diffuseColor ) * AiShaderEvalParamFlt( p_diffuseScale );
	AtFloat sssMix = AiShaderEvalParamFlt( p_sssMix );
	AtRGB sssRadiusColor = AiShaderEvalParamRGB( p_sssRadiusColor );
	float sssRadius = AiShaderEvalParamFlt( p_sssRadius );
	float sssScale = AiShaderEvalParamFlt( p_sssScale );
	AtRGB specular1Color = AiShaderEvalParamRGB( p_specular1Color ) * AiShaderEvalParamFlt( p_specular1Scale );
	AtFloat roughness = AiShaderEvalParamFlt( p_specular1Roughness );
	roughness *= roughness;
	AtFloat ior = AiShaderEvalParamFlt( p_specular1Ior );
	AtFloat eta = 1.0f / ior;

	AtRGB ssRadiusColor = AiShaderEvalParamRGB( p_ssRadiusColor );
	float ssRadius = AiShaderEvalParamFlt( p_ssRadius );
	float ssScale = AiShaderEvalParamFlt( p_ssScale );

	AtRGB transmissionColor = AiShaderEvalParamRGB(p_transmissionColor) * AiShaderEvalParamFlt(p_transmissionScale);

	AtFloat transmissionRoughness;
	AtFloat transmissionIor;
	bool transmissionLinkToSpecular1 = AiShaderEvalParamBool(p_transmissionLinkToSpecular1);
	if (transmissionLinkToSpecular1)
	{
		transmissionRoughness = roughness;
		transmissionIor = ior;
	}
	else
	{
		transmissionRoughness = AiShaderEvalParamFlt(p_transmissionRoughness);
		transmissionRoughness *= transmissionRoughness;
		transmissionIor = AiShaderEvalParamFlt(p_transmissionIor);
	}

	AtRGB absorption = AI_RGB_BLACK;
	if (AiShaderEvalParamBool(p_absorptionEnable))
	{
		absorption = (AI_RGB_WHITE - AiShaderEvalParamRGB(p_absorptionColor));
		absorption = max(rgb(AI_EPSILON), absorption) * AiShaderEvalParamFlt(p_absorptionDensity);
	}

	if (sg->Rt & AI_RAY_SHADOW)
	{
		//Kettle_shadows(fresnel(costheta, 1.0f/1.5f), 1.0f, false, AI_RGB_WHITE, sg, node);
		float costheta = AiV3Dot(sg->Nf, -sg->Rd);
		sg->out_opacity = fresnel(costheta, 1.0f/transmissionIor);
		return;
	}

	// Initialize result temporaries
	AtRGB result_diffuseDirect = AI_RGB_BLACK;
	AtRGB result_glossyDirect = AI_RGB_BLACK;
	AtRGB result_diffuseIndirect = AI_RGB_BLACK;
	AtRGB result_glossyIndirect = AI_RGB_BLACK;
	AtRGB result_sss = AI_RGB_BLACK;
	AtRGB result_ss = AI_RGB_BLACK;
	AtColor	result_transmission = AI_RGB_BLACK;
	// Set up flags to early out of calculations based on where we are in the ray tree
	bool do_diffuse = true;
	bool do_glossy = true;
	bool do_ss = true;
	bool do_sss = true;
	bool do_transmission = true;
	AtInt glossy_samples = data->GI_glossy_samples;
	AtInt diffuse_samples = data->GI_diffuse_samples;

	if ( sg->Rr_diff > data->GI_diffuse_depth || maxh(diffuseColor) < 0.01 )
	{
		do_diffuse = false;
	}


	if (sg->Rr_gloss > data->GI_glossy_depth || sg->Rr_diff > 0 || maxh(specular1Color) < 0.01)
	{
		do_glossy = false;
	}


	if ( sg->Rr_diff > 0 || sg->Rr_gloss > 0 || ssScale < 0.01f )
	{
		do_ss = false;
	}

	if ( sg->Rr_diff > 0 || sg->Rr_gloss > 1 || sssMix < 0.01f || !do_diffuse )
	{
		do_sss = false;
		sssMix = 0.0f;
	}

	if (sg->Rr_diff >  0 || maxh(transmissionColor) < 0.01)
	{
		do_transmission = false;
	}

	// build a local frame for sampling
	AtVector U, V;
	AiBuildLocalFramePolar(&U, &V, &sg->N);

	AtVector wo = -sg->Rd;

	// Begin illumination calculation
	if ( do_diffuse || do_glossy )
	{
		// Create the BRDF data structures for MIS
		void* mis;
		mis = GlossyMISCreateData(sg,&U,&V,roughness,roughness);
		BrdfData_wrap brdfw;
		brdfw.brdf_data = mis;
		brdfw.sg = sg;
		brdfw.eta = eta;
		brdfw.V = wo;
		brdfw.N = sg->N;

		void* dmis;
		dmis = AiOrenNayarMISCreateData(sg, 0.0f);
		BrdfData_wrap brdfd;
		brdfd.brdf_data = dmis;
		brdfd.sg = sg;
		brdfd.eta = eta;
		brdfd.V = wo;
		brdfd.N = sg->N;


		// Light loop
		AiLightsPrepare(sg);
		while(AiLightsGetSample(sg))
		{
			if (do_diffuse)
			{
				AtRGB Li = AiEvaluateLightSample(sg,&brdfd,AiOrenNayarMISSample_wrap,AiOrenNayarMISBRDF_wrap, AiOrenNayarMISPDF_wrap);
				result_diffuseDirect += Li*diffuseColor;
			}
			if (do_glossy)
			{
				result_glossyDirect +=
				AiEvaluateLightSample(sg,&brdfw,GlossyMISSample_wrap,GlossyMISBRDF_wrap,GlossyMISPDF_wrap)
				*specular1Color;
			}
		}

		// Sample BRDFS
		double samples[2];
		AtRay wi_ray;
		AtVector wi;
		AtScrSample scrs;
		AtVector H;
		float kr=1, kt=1;

		if (do_glossy)
		{
			glossy_samples *= glossy_samples;
			AtSamplerIterator* sampit = AiSamplerIterator(data->glossy_sampler, sg);
			AiMakeRay(&wi_ray, AI_RAY_GLOSSY, &sg->P, NULL, AI_BIG, sg);
			AtInt count=0;
			while(AiSamplerGetSample(sampit, samples))
			{
				wi = GlossyMISSample(mis, samples[0], samples[1]);
				if (AiV3Dot(wi,sg->Nf) > 0.0f)
				{
					wi_ray.dir = wi;
					AiV3Normalize(H, wi+brdfw.V);
					kr = fresnel(std::max(0.0f,AiV3Dot(H,wi)),eta);
					if (kr > IMPORTANCE_EPS) // only trace a ray if it's going to matter
					{
						AiTrace(&wi_ray, &scrs);
						result_glossyIndirect +=
						scrs.color*GlossyMISBRDF(mis, &wi) / GlossyMISPDF(mis, &wi) * kr * specular1Color;
					}
				}
				count++;
			}
			if (count) result_glossyIndirect /= float(count);
		} // if (do_glossy)

		if ( do_diffuse )
		{
			diffuse_samples *= diffuse_samples;
			AtSamplerIterator* sampit = AiSamplerIterator(data->diffuse_sampler, sg);
			AiMakeRay( &wi_ray, AI_RAY_DIFFUSE, &sg->P, NULL, AI_BIG, sg );
			AtInt count=0;
			while(AiSamplerGetSample(sampit, samples))
			{
				Imath::V3f wi_im = cosineSampleHemisphere(samples[0], samples[1]);
				Imath::V3f N_im( sg->N.x, sg->N.y, sg->N.z );
				Imath::M44f rm = rotationMatrix(Imath::V3f(0,1,0),N_im);
				rm.multDirMatrix(wi_im,wi_im);
				wi.x = wi_im.x;wi.y=wi_im.y;wi.z=wi_im.z;
				wi_ray.dir = wi;
				AiV3Normalize(H, wi+brdfw.V);
				kt = 1.0f - fresnel(std::max(0.0f,AiV3Dot(H,wi)),eta);
				if (kt > IMPORTANCE_EPS)
				{
					AiTrace(&wi_ray, &scrs);
					result_diffuseIndirect += scrs.color * kt * AI_ONEOVERPI;
				}
				count++;
			}
			if (count) result_diffuseIndirect /= float(count);
			result_diffuseIndirect *= diffuseColor;
		} // if (do_diffuse)

	} // if (do_diffuse || do_glossy)

	// Diffusion multiple scattering
	if ( do_sss )
	{
		result_sss = AiSSSPointCloudLookupCubic(sg, sssRadius*sssRadiusColor) * diffuseColor;
	}

	// Single-scattering
	if (do_ss)
	{
		AtRGB sigma_s_prime, sigma_a;
		alphaInversion( ssRadiusColor*ssRadius, ssRadius, sigma_s_prime, sigma_a );
		AtRGB sigma_t_prime = (sigma_s_prime+sigma_a);
		AtRGB mfp = AI_RGB_WHITE / sigma_t_prime;
		AtRGB alpha_prime = sigma_s_prime / sigma_t_prime;
		result_ss = AiSSSTraceSingleScatter(sg,AI_RGB_WHITE,mfp,0.0f,1.3) / brdf(alpha_prime) * ssScale;
	}

	// blend sss and direct diffuse
	result_diffuseDirect *= (1-sssMix);
	result_diffuseIndirect *= (1-sssMix);
	result_sss *= sssMix;

	// Refraction
	if (do_transmission)
	{
		//microfacetRefraction(sg, data, transmissionIor, transmissionRoughness, absorption, result_transmission);
		result_transmission = beckmannMicrofacetTransmission(sg, sg->N, U, V, wo, data->refraction_sampler,
																transmissionRoughness, transmissionIor, absorption);
	}

	if (sg->Rt & AI_RAY_CAMERA)
	{
		// write AOVs
		AiAOVSetRGB(sg, "diffuseDirect", result_diffuseDirect);
		AiAOVSetRGB(sg, "multiScatter", result_sss);
		AiAOVSetRGB(sg, "specularDirect", result_glossyDirect);
		AiAOVSetRGB(sg, "diffuseIndirect", result_diffuseIndirect);
		AiAOVSetRGB(sg, "specularIndirect", result_glossyIndirect);
		AiAOVSetRGB(sg, "singleScatter", result_ss);
		AiAOVSetRGB(sg, "transmission", result_transmission);
	}

	// Sum final result from temporaries
	//
	sg->out.RGB =  	 result_diffuseDirect
					+result_sss
					+result_glossyDirect
					+result_diffuseIndirect
					+result_glossyIndirect
					+result_ss
					+result_transmission;
}
