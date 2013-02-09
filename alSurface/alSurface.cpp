#include <iostream>
#include <ai_sampler.h>
#include <OpenEXR/ImathVec.h>
#include <OpenEXR/ImathMatrix.h>
#include <OpenEXR/ImathMatrixAlgo.h>
#include <map>

#include "alUtil.h"
#include "MIS.h"
#include "BeckmannMicrofacet.h"
#include "alSurface.h"

AI_SHADER_NODE_EXPORT_METHODS(alSurfaceMtd)

#define GlossyMISBRDF AiCookTorranceMISBRDF
#define GlossyMISPDF AiCookTorranceMISPDF
#define GlossyMISSample AiCookTorranceMISSample

#define GlossyMISBRDF_wrap AiCookTorranceMISBRDF_wrap
#define GlossyMISPDF_wrap AiCookTorranceMISPDF_wrap
#define GlossyMISSample_wrap AiCookTorranceMISSample_wrap

#define GlossyMISCreateData AiCookTorranceMISCreateData

#define NUM_LIGHT_GROUPS 8
static const char* lightGroupNames[NUM_LIGHT_GROUPS] =
{
	"lightGroup1",
	"lightGroup2",
	"lightGroup3",
	"lightGroup4",
	"lightGroup5",
	"lightGroup6",
	"lightGroup7",
	"lightGroup8"
};

enum alSurfaceParams
{
	// diffuse
	p_diffuseStrength=0,
	p_diffuseColor,
	p_diffuseRoughness,
	p_emissionStrength,
	p_emissionColor,

	// sss
	p_sssMix,
	p_sssRadius,
	p_sssRadiusColor,
	p_sssDensityScale,

	p_ssStrength,
	p_ssBalance,
	p_ssTargetColor,
	p_ssSpecifyCoefficients,
	p_ssScattering,
	p_ssAbsorption,
	p_ssDensityScale,
	p_ssDirection,
	p_ssInScattering,

	p_diffuseExtraSamples,
	p_diffuseEnableCaustics,

	// specular
	p_specular1Strength,
	p_specular1Color,
	p_specular1Roughness,
	p_specular1Ior,
	p_specular1RoughnessDepthScale,
	p_specular1ExtraSamples,
	p_specular1Normal,
	p_specular2Strength,
	p_specular2Color,
	p_specular2Roughness,
	p_specular2Ior,
	p_specular2RoughnessDepthScale,
	p_specular2ExtraSamples,
	p_specular2Normal,

	// transmission
	p_transmissionStrength,
	p_transmissionColor,
	p_transmissionLinkToSpecular1,
	p_transmissionRoughness,
	p_transmissionIor,
	p_transmissionRoughnessDepthScale,
	p_transmissionEnableCaustics,
	p_transmissionExtraSamples,


	p_bump
};

node_parameters
{
	AiParameterFLT("diffuseStrength", 1.0f );
	AiParameterRGB("diffuseColor", 0.18f, 0.18f, 0.18f );
	AiParameterFLT("diffuseRoughness", 0.0f );
	AiParameterFLT("emissionStrength", 0.0f );
	AiParameterRGB("emissionColor", 1.0f, 1.0f, 1.0f);

	AiParameterFLT("sssMix", 0.0f );
	AiParameterFLT("sssRadius", 3.6f );
	AiParameterRGB("sssRadiusColor", .439f, .156f, .078f );
	AiMetaDataSetBool(mds, "sssRadiusColor", "always_linear", true);  // no inverse-gamma correction
	AiParameterFLT("sssDensityScale", 1.0f );

	AiParameterFLT("ssStrength", 0.0f );
	AiParameterFLT("ssBalance", 0.5f);
	AiParameterRGB("ssTargetColor", .439f, .156f, .078f);
	AiParameterBOOL("ssSpecifyCoefficients", false);
	AiParameterRGB("ssScattering", 1.0f, 1.0f, 1.0f);
	AiParameterRGB("ssAbsorption", 1.0f, 1.0f, 1.0f);
	AiParameterFLT("ssDensityScale", 1.0f);
	AiParameterFLT("ssDirection", 0.0f);
	AiParameterBOOL("ssInScattering", true);

	AiParameterINT("diffuseExtraSamples", 0);
	AiParameterBOOL("diffuseEnableCaustics", false);

	AiParameterFLT("specular1Strength", 1.0f );
	AiParameterRGB("specular1Color", 1.0f, 1.0f, 1.0f );
	AiParameterFLT("specular1Roughness", 0.3f );
	AiParameterFLT("specular1Ior", 1.4f );
	AiParameterFLT("specular1RoughnessDepthScale", 1.0f);
	AiParameterINT("specular1ExtraSamples", 0);
	AiParameterVec("specular1Normal", 0, 0, 0);

	AiParameterFLT("specular2Strength", 0.0f );
	AiParameterRGB("specular2Color", 1.0f, 1.0f, 1.0f );
	AiParameterFLT("specular2Roughness", 0.3f );
	AiParameterFLT("specular2Ior", 1.4f );
	AiParameterFLT("specular2RoughnessDepthScale", 1.0f);
	AiParameterINT("specular2ExtraSamples", 0);
	AiParameterVec("specular2Normal", 0, 0, 0);

	AiParameterFLT("transmissionStrength", 0.0f );
	AiParameterRGB("transmissionColor", 1.0f, 1.0f, 1.0f );
	AiParameterBOOL("transmissionLinkToSpecular1", true);
	AiParameterFLT("transmissionRoughness", 0.1f );
	AiParameterFLT("transmissionIor", 1.4f );
	AiParameterFLT("transmissionRoughnessDepthScale", 1.0f);
	AiParameterBOOL("transmissionEnableCaustics", true);
	AiParameterINT("transmissionExtraSamples", 0);

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
	ShaderData *data = new ShaderData;
	AiNodeSetLocalData(node,data);
	data->diffuse_sampler = NULL;
	data->glossy_sampler = NULL;
	data->glossy2_sampler = NULL;
	data->refraction_sampler = NULL;
};

node_finish
{
	if (AiNodeGetLocalData(node))
	{
		ShaderData* data = (ShaderData*) AiNodeGetLocalData(node);

		AiSamplerDestroy(data->diffuse_sampler);
		AiSamplerDestroy(data->glossy_sampler);
		AiSamplerDestroy(data->glossy2_sampler);
		AiSamplerDestroy(data->refraction_sampler);

		AiNodeSetLocalData(node, NULL);
		delete data;
	}
}


node_update
{
	// set up AOVs
	AiAOVRegister("diffuseDirect", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("diffuseIndirect", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("specularDirect", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("specularIndirect", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("specular2Direct", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("specular2Indirect", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("singleScatter", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("multiScatter", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("transmissionIndirect", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("emission", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("uv", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("depth", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("lightGroup1", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("lightGroup2", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("lightGroup3", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("lightGroup4", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("lightGroup5", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("lightGroup6", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("lightGroup7", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("lightGroup8", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);

	// store some options we'll reuse later
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
	AiSamplerDestroy(data->glossy2_sampler);
	AiSamplerDestroy(data->refraction_sampler);
	data->diffuse_sampler = AiSampler(data->GI_diffuse_samples+params[p_diffuseExtraSamples].INT, 2);
	data->glossy_sampler = AiSampler(data->GI_glossy_samples+params[p_specular1ExtraSamples].INT, 2);
	data->glossy2_sampler = AiSampler(data->GI_glossy_samples+params[p_specular2ExtraSamples].INT, 2);
	data->refraction_sampler = AiSampler(AiNodeGetInt(options, "GI_refraction_samples")+params[p_transmissionExtraSamples].INT, 2);

	// Get all the light nodes in the scene and try and find their light group parameter
	// we'll store this based on the light pointer for fast access during rendering
	AtNodeIterator* it = AiUniverseGetNodeIterator(AI_NODE_LIGHT);
	while (!AiNodeIteratorFinished(it))
	{
		AtNode* light = AiNodeIteratorGetNext(it);
		data->lightGroups[light] = AiNodeGetInt(light, "lightGroup") - 1;
	}
	AiNodeIteratorDestroy(it);
};


shader_evaluate
{
	ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);

	AtFloat ior = std::max(1.001f, AiShaderEvalParamFlt(p_specular1Ior));
	AtFloat roughness = AiShaderEvalParamFlt( p_specular1Roughness );
	roughness *= roughness;
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
		transmissionIor = std::max(1.001f, AiShaderEvalParamFlt(p_transmissionIor));
	}
	AtRGB transmissionColor = AiShaderEvalParamRGB(p_transmissionColor) * AiShaderEvalParamFlt(p_transmissionStrength);

	AtRGB ssScattering = AiShaderEvalParamRGB(p_ssScattering);
	AtRGB ssAbsorption = AiShaderEvalParamRGB(p_ssAbsorption);
	AtFloat ssDensityScale = AiShaderEvalParamFlt( p_ssDensityScale );
	AtFloat ssStrength = AiShaderEvalParamFlt( p_ssStrength );
	AtFloat ssDirection = AiShaderEvalParamFlt(p_ssDirection);
	AtFloat ssBalance = AiShaderEvalParamFlt(p_ssBalance);
	AtRGB ssTargetColor = AiShaderEvalParamRGB(p_ssTargetColor);
	bool ssSpecifyCoefficients = AiShaderEvalParamBool(p_ssSpecifyCoefficients);
	bool ssInScattering = AiShaderEvalParamBool(p_ssInScattering);

	// precalculate scattering coefficients as we'll need them for shadows etc.
	AtRGB sigma_t = AI_RGB_BLACK;
	AtRGB sigma_s = AI_RGB_BLACK;
	AtRGB sigma_a = AI_RGB_BLACK;
	if (ssStrength > IMPORTANCE_EPS)
	{
		if (ssSpecifyCoefficients)
		{
			sigma_s = ssScattering * ssDensityScale;
			sigma_a = ssAbsorption * ssDensityScale;
			sigma_t = sigma_s + sigma_a;
		}
		else
		{
			sigma_s = sigma_a = AI_RGB_WHITE - ssTargetColor;
			sigma_s *= ssBalance * ssDensityScale;
			sigma_a *= (1.0f - ssBalance) * ssDensityScale;
			sigma_t = sigma_s + sigma_a;
		}
	}

	// if it's a shadow ray, handle shadow colouring through absorption
	// algorithm based heavily on the example in Kettle
	if (sg->Rt & AI_RAY_SHADOW)
	{
		// if the object is transmissive and
		AtRGB outOpacity = AI_RGB_WHITE;
		if (maxh(transmissionColor))
		{
			// check transmission through the surface
			AtFloat costheta = AiV3Dot(sg->Nf, -sg->Rd);
			AtFloat kt = 1.0f - fresnel(costheta, 1.0f/transmissionIor);
			if (kt >= IMPORTANCE_EPS) // else surface is fully reflective
			{
				if (maxh(sigma_t) > 0.0f)
				{
					AtPoint alsPreviousIntersection;
					AtRGB als_sigma_t = sigma_t;
					if (AiStateGetMsgPnt("alsPreviousIntersection", &alsPreviousIntersection))
					{
						AiStateGetMsgRGB("alsPrevious_sigma_t", &als_sigma_t);
						bool doExtinction = false;
						if (AiV3Dot(sg->N, sg->Rd) < 0.0f)
						{
							// ray is entering a closed volume
							bool alsInside;
							AiStateGetMsgBool("alsInside", &alsInside);
							if (alsInside)
							{
								// ray is entering an embedded volume
								doExtinction = true;
							}
							else
							{
								// shouldn't get here
							}

						}
						else
						{
							// ray is exiting a closed volume
							doExtinction = true;
						}

						if (doExtinction)
						{
							AtFloat z = AiV3Dist(sg->P, alsPreviousIntersection);
							outOpacity.r = expf(-z * als_sigma_t.r);
							outOpacity.g = expf(-z * als_sigma_t.g);
							outOpacity.b = expf(-z * als_sigma_t.b);
							outOpacity = 1.0f - (outOpacity*kt);
						}

					}
					else
					{
						// first intersection
						// tell the next shader invocation that we're now inside the surface and what our extinction
						// coefficient is
						AiStateSetMsgRGB("alsPrevious_sigma_t", sigma_t);
						AiStateSetMsgBool("alsInside", true);
					}
				}
				else // no extinction, shadows are fresnel only.
				{
					AiStateSetMsgRGB("alsPrevious_sigma_t", AI_RGB_BLACK);
					outOpacity = 1.0f - kt;
				}
			}
		}

		// store intersection position
		AiStateSetMsgPnt("alsPreviousIntersection", sg->P);
		sg->out_opacity = outOpacity;
		return;
	}

	// Evaluate bump;
	AtVector N_orig;
	AtVector Nf_orig;
	AtRGB bump = AiShaderEvalParamRGB( p_bump );

	// Initialize parameter temporaries
	// TODO: reorganize this so we're not evaluating upstream when we don't need the parameters, e.g. in shadow rays
	AtRGB diffuseColor = AiShaderEvalParamRGB( p_diffuseColor ) * AiShaderEvalParamFlt( p_diffuseStrength );
	AtFloat diffuseRoughness = AiShaderEvalParamFlt(p_diffuseRoughness);
	bool diffuseEnableCaustics = AiShaderEvalParamFlt(p_diffuseEnableCaustics);
	AtRGB emissionColor = AiShaderEvalParamRGB(p_emissionColor) * AiShaderEvalParamFlt(p_emissionStrength);
	AtFloat sssMix = AiShaderEvalParamFlt( p_sssMix );
	AtRGB sssRadiusColor = AiShaderEvalParamRGB( p_sssRadiusColor );
	AtFloat sssRadius = AiShaderEvalParamFlt( p_sssRadius );
	AtFloat sssDensityScale = AiShaderEvalParamFlt( p_sssDensityScale );
	AtRGB specular1Color = AiShaderEvalParamRGB( p_specular1Color ) * AiShaderEvalParamFlt( p_specular1Strength );
	AtRGB specular2Color = AiShaderEvalParamRGB( p_specular2Color ) * AiShaderEvalParamFlt( p_specular2Strength );
	AtVector specular1Normal = sg->N;
	if (AiNodeIsLinked(node, "specular1Normal"))
	{
		specular1Normal = AiShaderEvalParamVec(p_specular1Normal);
	}

	AtVector specular2Normal = sg->N;
	if (AiNodeIsLinked(node, "specular2Normal"))
	{
		specular2Normal = AiShaderEvalParamVec(p_specular2Normal);
	}

	AtFloat roughness2 = AiShaderEvalParamFlt( p_specular2Roughness );
	roughness2 *= roughness2;

	AtFloat eta = 1.0f / ior;
	AtFloat ior2 = std::max(1.001f, AiShaderEvalParamFlt(p_specular2Ior));
	AtFloat eta2 = 1.0f / ior2;


	AtFloat specular1RoughnessDepthScale = AiShaderEvalParamFlt(p_specular1RoughnessDepthScale);
	AtFloat specular2RoughnessDepthScale = AiShaderEvalParamFlt(p_specular2RoughnessDepthScale);
	AtFloat transmissionRoughnessDepthScale = AiShaderEvalParamFlt(p_transmissionRoughnessDepthScale);


	bool transmissionEnableCaustics = AiShaderEvalParamBool(p_transmissionEnableCaustics);

	// clamp roughnesses
	// TODO: fall back to single-ray solution when roughness is 0
	roughness = std::max(0.0001f, roughness);
	roughness2 = std::max(0.0001f, roughness2);
	transmissionRoughness = std::max(0.0001f, transmissionRoughness);

	// Initialize result temporaries
	AtRGB result_diffuseDirect = AI_RGB_BLACK;
	AtRGB result_glossyDirect = AI_RGB_BLACK;
	AtRGB result_glossy2Direct = AI_RGB_BLACK;
	AtRGB result_diffuseIndirect = AI_RGB_BLACK;
	AtRGB result_glossyIndirect = AI_RGB_BLACK;
	AtRGB result_glossy2Indirect = AI_RGB_BLACK;
	AtRGB result_sss = AI_RGB_BLACK;
	AtRGB result_ss = AI_RGB_BLACK;
	AtColor	result_transmission = AI_RGB_BLACK;
	AtColor result_emission = AI_RGB_BLACK;
	// Set up flags to early out of calculations based on where we are in the ray tree
	bool do_diffuse = true;
	bool do_glossy = true;
	bool do_glossy2 = true;
	bool do_ss = true;
	bool do_sss = true;
	bool do_transmission = true;
	AtInt glossy_samples = data->GI_glossy_samples;
	AtInt diffuse_samples = data->GI_diffuse_samples;

	if ( sg->Rr_diff > data->GI_diffuse_depth || maxh(diffuseColor) < IMPORTANCE_EPS || sssMix == 1.0f)
	{
		do_diffuse = false;
	}


	if (sg->Rr_gloss > data->GI_glossy_depth
				|| (sg->Rr_diff > 0)			// disable glossy->diffuse caustics
				|| maxh(specular1Color) < IMPORTANCE_EPS				// skip evaluations that aren't important
				|| (sg->Rr_refr > 1 && !transmissionEnableCaustics))	// disable glossy->transmitted caustics
	{
		do_glossy = false;
	}

	if (sg->Rr_gloss > data->GI_glossy_depth
			|| (sg->Rr_diff > 0)			// disable glossy->diffuse caustics
			|| maxh(specular2Color) < IMPORTANCE_EPS				// skip evaluations that aren't important
			|| (sg->Rr_refr > 1 && !transmissionEnableCaustics))	// disable glossy->transmitted caustics
	{
		do_glossy2 = false;
	}

	if (sg->Rr_gloss > 0 && do_glossy)
	{
		roughness *= powf(specular1RoughnessDepthScale, sg->Rr_gloss);
		roughness2 *= powf(specular2RoughnessDepthScale, sg->Rr_gloss);
	}
	if (sg->Rr_refr > 0 && do_transmission)
	{
		transmissionRoughness *= powf(transmissionRoughnessDepthScale, sg->Rr_refr);
	}

	if (sg->Rr_diff > 0 || sg->Rr_gloss > 1 || sssMix < 0.01f)
	{
		do_sss = false;
		sssMix = 0.0f;
	}

	// make sure diffuse and transmission can't sum > 1
	transmissionColor *= 1.0f - maxh(diffuseColor);

	if ((sg->Rr_diff > 0) || maxh(transmissionColor) < IMPORTANCE_EPS)
	{
		do_transmission = false;
	}

	// build a local frame for sampling
	AtVector U, V;
	AiBuildLocalFramePolar(&U, &V, &sg->N);

	AtVector wo = -sg->Rd;

	// prepare temporaries for light group calculation
	AtRGB lightGroupDiffuse[NUM_LIGHT_GROUPS];
	memset(lightGroupDiffuse, 0, sizeof(AtRGB)*NUM_LIGHT_GROUPS);
	AtRGB lightGroupSpecular[NUM_LIGHT_GROUPS];
	memset(lightGroupSpecular, 0, sizeof(AtRGB)*NUM_LIGHT_GROUPS);
	AtRGB lightGroupSpecular2[NUM_LIGHT_GROUPS];
	memset(lightGroupSpecular2, 0, sizeof(AtRGB)*NUM_LIGHT_GROUPS);

	// Accumulator for transmission integrated according to the specular1 brdf. Will be used to attenuate diffuse,
	// glossy2, sss and transmission
	AtFloat kti = 1.0f;
	AtFloat kti2 = 1.0f;

	// Begin illumination calculation
	if (do_diffuse || do_glossy || do_glossy2)
	{
		// Create the BRDF data structures for MIS
		AtVector Nold = sg->N;
		sg->N = sg->Nf = specular1Normal;
		void* mis;
		mis = GlossyMISCreateData(sg,&U,&V,roughness,roughness);
		BrdfData_wrap brdfw;
		brdfw.brdf_data = mis;
		brdfw.sg = sg;
		brdfw.eta = eta;
		brdfw.V = wo;
		brdfw.N = specular1Normal;
		brdfw.kr = 0.0f;


		sg->N = sg->Nf = specular2Normal;	
		void* mis2;
		mis2 = GlossyMISCreateData(sg,&U,&V,roughness2,roughness2);
		BrdfData_wrap brdfw2;
		brdfw2.brdf_data = mis2;
		brdfw2.sg = sg;
		brdfw2.eta = eta2;
		brdfw2.V = wo;
		brdfw2.N = specular2Normal;

		sg->N = Nold;

		void* dmis;
		dmis = AiOrenNayarMISCreateData(sg, diffuseRoughness);
		BrdfData_wrap brdfd;
		brdfd.brdf_data = dmis;
		brdfd.sg = sg;
		brdfd.eta = eta;
		brdfd.V = wo;
		brdfd.N = sg->N;

		// Light loop
		AiLightsPrepare(sg);
		if (sg->Rt & AI_RAY_CAMERA)
		{
			AtRGB LdiffuseDirect, LspecularDirect, Lspecular2Direct;
			while(AiLightsGetSample(sg))
			{
				int lightGroup = data->lightGroups[sg->Lp];
				if (do_glossy)
				{
					sg->N = sg->Nf = specular1Normal;
					LspecularDirect =
					AiEvaluateLightSample(sg,&brdfw,GlossyMISSample_wrap,GlossyMISBRDF_wrap,GlossyMISPDF_wrap);
					if (lightGroup >= 0 && lightGroup < NUM_LIGHT_GROUPS)
					{
						lightGroupSpecular[lightGroup] += LspecularDirect;
					}
					result_glossyDirect += LspecularDirect;
					sg->N = Nold;
				}
				if (do_glossy2)
				{
					sg->N = sg->Nf = specular2Normal;
					Lspecular2Direct =
					AiEvaluateLightSample(sg,&brdfw2,GlossyMISSample_wrap,GlossyMISBRDF_wrap,GlossyMISPDF_wrap)
											* (1.0f - brdfw.kr*maxh(specular1Color));
					if (lightGroup >= 0 && lightGroup < NUM_LIGHT_GROUPS)
					{
						lightGroupSpecular2[lightGroup] += Lspecular2Direct;
					}
					result_glossy2Direct += Lspecular2Direct;
					sg->N = Nold;
				}
				if (do_diffuse)
				{
					LdiffuseDirect =
					AiEvaluateLightSample(sg,dmis,AiOrenNayarMISSample,AiOrenNayarMISBRDF, AiOrenNayarMISPDF)
											* (1.0f - brdfw.kr*maxh(specular1Color))
											* (1.0f - brdfw2.kr*maxh(specular2Color));
					if (lightGroup >= 0 && lightGroup < NUM_LIGHT_GROUPS)
					{
						lightGroupDiffuse[lightGroup] += LdiffuseDirect;
					}
					result_diffuseDirect += LdiffuseDirect;
				}
			}
			for (int i=0; i < NUM_LIGHT_GROUPS; ++i)
			{
				lightGroupDiffuse[i] *= diffuseColor;
				lightGroupSpecular[i] *= specular1Color;
				lightGroupSpecular2[i] *= specular2Color;
			}
		}
		else
		{
			while(AiLightsGetSample(sg))
			{
				if (do_glossy)
				{
					result_glossyDirect +=
					AiEvaluateLightSample(sg,&brdfw,GlossyMISSample_wrap,GlossyMISBRDF_wrap,GlossyMISPDF_wrap);
				}
				if (do_glossy2)
				{
					result_glossy2Direct +=
					AiEvaluateLightSample(sg,&brdfw2,GlossyMISSample_wrap,GlossyMISBRDF_wrap,GlossyMISPDF_wrap)
											* (1.0f - brdfw.kr*maxh(specular1Color));
				}
				if (do_diffuse)
				{
					result_diffuseDirect +=
					AiEvaluateLightSample(sg,&brdfd,AiOrenNayarMISSample_wrap,AiOrenNayarMISBRDF_wrap, AiOrenNayarMISPDF_wrap)
											* (1.0f - brdfw.kr*maxh(specular1Color))
											* (1.0f - brdfw2.kr*maxh(specular2Color));
				}
			}
		}

		// Multiply by the colors
		result_diffuseDirect *= diffuseColor;
		result_glossyDirect *= specular1Color;
		result_glossy2Direct *= specular2Color;


		// Sample BRDFS
		double samples[2];
		AtRay wi_ray;
		AtVector wi;
		AtScrSample scrs;
		AtVector H;
		float kr=1, kt=1;

		if (do_glossy)
		{
			AtSamplerIterator* sampit = AiSamplerIterator(data->glossy_sampler, sg);
			AiMakeRay(&wi_ray, AI_RAY_GLOSSY, &sg->P, NULL, AI_BIG, sg);
			kti = 0.0f;
			while(AiSamplerGetSample(sampit, samples))
			{
				wi = GlossyMISSample(mis, samples[0], samples[1]);
				if (AiV3Dot(wi,specular1Normal) > 0.0f)
				{
					// get half-angle vector for fresnel
					wi_ray.dir = wi;
					AiV3Normalize(H, wi+brdfw.V);
					kr = fresnel(std::max(0.0f,AiV3Dot(H,wi)),eta);
					kti += kr;
					if (kr > IMPORTANCE_EPS) // only trace a ray if it's going to matter
					{
						AiTrace(&wi_ray, &scrs);
						result_glossyIndirect +=
						scrs.color*GlossyMISBRDF(mis, &wi) / GlossyMISPDF(mis, &wi) * kr;
					}
				}
			}
			result_glossyIndirect *= AiSamplerGetSampleInvCount(sampit);
			kti *= AiSamplerGetSampleInvCount(sampit);
			kti = 1.0f - kti*maxh(specular1Color);
			result_glossyIndirect *= specular1Color;
		} // if (do_glossy)

		if (do_glossy2)
		{
			AtSamplerIterator* sampit = AiSamplerIterator(data->glossy2_sampler, sg);
			AiMakeRay(&wi_ray, AI_RAY_GLOSSY, &sg->P, NULL, AI_BIG, sg);
			kti2 = 0.0f;
			while(AiSamplerGetSample(sampit, samples))
			{
				wi = GlossyMISSample(mis2, samples[0], samples[1]);
				if (AiV3Dot(wi,specular2Normal) > 0.0f)
				{
					wi_ray.dir = wi;
					AiV3Normalize(H, wi+brdfw2.V);
					// add the fresnel for this layer
					kr = fresnel(std::max(0.0f,AiV3Dot(H,wi)),eta2);
					if (kr > IMPORTANCE_EPS) // only trace a ray if it's going to matter
					{
						AiTrace(&wi_ray, &scrs);
						result_glossy2Indirect +=
						scrs.color*GlossyMISBRDF(mis2, &wi) / GlossyMISPDF(mis2, &wi) * kr * kti;
						kti2 += kr; 
					}
				}
			}
			result_glossy2Indirect*= AiSamplerGetSampleInvCount(sampit);
			kti2 *= AiSamplerGetSampleInvCount(sampit);
			kti2 = 1.0f - kti2*maxh(specular2Color);
			result_glossy2Indirect *= specular2Color;
		} // if (do_glossy2)

		if ( do_diffuse )
		{
			result_diffuseIndirect = AiOrenNayarIntegrate(&sg->Nf, sg, diffuseRoughness) * diffuseColor * kti * kti2;
		} // if (do_diffuse)

	} // if (do_diffuse || do_glossy)

	// Emission
	result_emission = emissionColor;

	// Diffusion multiple scattering
	if (do_sss)
	{
		result_sss = AiSSSPointCloudLookupCubic(sg, sssRadius*sssRadiusColor*sssDensityScale) * diffuseColor * kti * kti2;
	}


	// blend sss and direct diffuse
	result_diffuseDirect *= (1-sssMix);
	result_diffuseIndirect *= (1-sssMix);
	result_sss *= sssMix;

	// Refraction
	if (do_transmission)
	{
		result_transmission = beckmannMicrofacetTransmission(sg, sg->N, U, V, wo, data->refraction_sampler,
																transmissionRoughness, transmissionIor,
																sigma_s, sigma_a,
																ssDirection, ssStrength, ssInScattering, result_ss) * kti * kti2 * transmissionColor;
	}

	if (sg->Rt & AI_RAY_CAMERA)
	{
		// write AOVs
		AiAOVSetRGB(sg, "diffuseDirect", result_diffuseDirect);
		AiAOVSetRGB(sg, "multiScatter", result_sss);
		AiAOVSetRGB(sg, "specularDirect", result_glossyDirect);
		AiAOVSetRGB(sg, "specular2Direct", result_glossy2Direct);
		AiAOVSetRGB(sg, "diffuseIndirect", result_diffuseIndirect);
		AiAOVSetRGB(sg, "specularIndirect", result_glossyIndirect);
		AiAOVSetRGB(sg, "specular2Indirect", result_glossy2Indirect);
		AiAOVSetRGB(sg, "singleScatter", result_ss);
		AiAOVSetRGB(sg, "transmissionIndirect", result_transmission);
		AiAOVSetRGB(sg, "emission", result_emission);

		// write light groups
		AtRGB lightGroup[NUM_LIGHT_GROUPS];
		for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
		{
			lightGroup[i] = lightGroupDiffuse[i] + lightGroupSpecular[i] + lightGroupSpecular2[i];
			AiAOVSetRGB(sg, lightGroupNames[i], lightGroup[i]);
		}

		// write data AOVs
		AtRGB uv = AiColorCreate(sg->u, sg->v, 0.0f);
		AiAOVSetRGB(sg, "uv", uv);
		AtRGB depth = AiColorCreate(sg->Rl, AiV3Dot(sg->Nf, wo), 0.0f);
		AiAOVSetRGB(sg, "depth", depth);

	}

	// Sum final result from temporaries
	//
	sg->out.RGB =  	 result_diffuseDirect
					+result_sss
					+result_glossyDirect
					+result_glossy2Direct
					+result_diffuseIndirect
					+result_glossyIndirect
					+result_glossy2Indirect
					+result_ss
					+result_transmission
					+result_emission;
}
