// Hair shader based on ISHair: Importance Sampling for Hair Scattering by Ou et. al 2012
// http://www.cs.dartmouth.edu/~ouj/site/Research/Entries/2012/6/21_ISHair__Importance_Sampling_for_Hair_Scattering.html

#include <ai.h>
#include "alUtil.h"

struct ShaderData
{
	ShaderData()
	: sampler_R(NULL), sampler_TT(NULL), sampler_TRT(NULL), sampler_g(NULL)
	{}
	~ShaderData()
	{
		AiSamplerDestroy(sampler_R);
		AiSamplerDestroy(sampler_TT);
		AiSamplerDestroy(sampler_TRT);
		AiSamplerDestroy(sampler_g);
	}

	void update(int glossy_samples)
	{
		AiSamplerDestroy(sampler_R);
		sampler_R = AiSampler(glossy_samples, 2);
		AiSamplerDestroy(sampler_TT);
		sampler_TT = AiSampler(glossy_samples, 2);
		AiSamplerDestroy(sampler_TRT);
		sampler_TRT = AiSampler(glossy_samples, 2);
		AiSamplerDestroy(sampler_g);
		sampler_g = AiSampler(glossy_samples, 2);
	}

	AtSampler* sampler_R;
	AtSampler* sampler_TT;
	AtSampler* sampler_TRT;
	AtSampler* sampler_g;
};

AI_SHADER_NODE_EXPORT_METHODS(alHair)

enum alHairParams
{
	p_ior,
	p_specularShift,
	p_specularWidth,
	p_specular1Strength,
	p_specular1Color,
	p_specular2Strength,
	p_specular2Color,
	p_glintStrength,
	p_glintRolloff,
	p_glintSeparation,
	p_transmissionStrength,
	p_transmissionColor,
	p_transmissionRolloff
};

node_parameters
{
	AiParameterFlt("ior", 1.55f);
	AiParameterFlt("specularShift", 7.0f);
	AiParameterFlt("specularWidth", 7.0f);
	AiParameterFlt("specular1Strength", 1.0f);
	AiParameterRGB("specular1Color", 1.0f, 1.0f, 1.0f);
	AiParameterFlt("specular2Strength", 1.0f);
	AiParameterRGB("specular2Color", 0.31f, 0.08f, 0.005f);
	AiParameterFlt("glintStrength", 2.0f);
	AiParameterFlt("glintRolloff", 15.0f);
	AiParameterFlt("glintSeparation", 35.0f);
	AiParameterFlt("transmissionStrength", 1.0f);
	AiParameterRGB("transmissionColor", 0.92f, 0.7f, 0.64f);
	AiParameterFlt("transmissionRolloff", 30.0f);
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alHair;
   node->output_type = AI_TYPE_RGB;
   node->name        = "alHair";
   node->node_type   = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return TRUE;
}

node_initialize
{
	ShaderData* data = new ShaderData;
	AiNodeSetLocalData(node, data);
}

node_finish
{
	if (AiNodeGetLocalData(node))
	{
		ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);
		delete data;
	}
}

node_update
{
	ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);
	AtNode *options   = AiUniverseGetOptions();
	int glossy_samples = AiNodeGetInt(options, "GI_glossy_samples");
	data->update(glossy_samples);

	AiAOVRegister("specularDirect", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("specularIndirect", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
}

AtFloat g(AtFloat beta, AtFloat alpha, AtFloat theta_h)
{
	AtFloat n = theta_h-alpha;
	return expf(-(n*n)/(2.0f*beta*beta));
}

// returns theta
AtFloat longitudinalSample(AtFloat u, AtFloat theta_r, AtFloat alpha, AtFloat beta)
{
	AtFloat A = atanf((AI_PI*0.25f + theta_r*0.5f - alpha) / beta);
	AtFloat B = atanf((-AI_PI*0.25f + theta_r*0.5f - alpha) / beta);
	return 2.0f*beta*tanf(u*(A-B)+B) + 2.0f*alpha - theta_r;
}

// returns p(theta)
AtFloat longitudinalPdf(AtFloat theta_i, AtFloat theta_r, AtFloat alpha, AtFloat beta)
{
	AtFloat A = atanf((AI_PI*0.25f + theta_r*0.5f - alpha) / beta);
	AtFloat B = atanf((-AI_PI*0.25f + theta_r*0.5f - alpha) / beta);
	AtFloat theta_h = (theta_r + theta_i)*0.5f;
	AtFloat x = 1.0f / (2.0f * cosf(theta_i) * (A-B));
	AtFloat y = beta / ((theta_h-alpha)*(theta_h-alpha)+beta*beta);
	return x*y;
}

// returns phi
// get phi_i = phi_r - phi
AtFloat NrSample(AtFloat u)
{
	return 2.0f * asinf(2.0f*u - 1);
}

AtFloat NrPdf(AtFloat phi)
{
	return cosf(phi*0.5f)*0.25f;
}

// returns phi
// get phi_i = phi_r - phi
AtFloat NttSample(AtFloat u, AtFloat gamma_TT)
{
	AtFloat C_TT = 2.0f*atanf(AI_PI/gamma_TT);
	return gamma_TT*tanf(C_TT*(u - 0.5f)) + AI_PI;
}

AtFloat NttPdf(AtFloat phi, AtFloat gamma_TT)
{
	AtFloat C_TT = 2.0f*atanf(AI_PI/gamma_TT);
	return (gamma_TT/((phi - AI_PI)*(phi - AI_PI) + gamma_TT*gamma_TT)) / C_TT;
}

AtFloat NgSample(AtFloat u, AtFloat gamma_g, AtFloat phi_g)
{
	AtFloat sign = 1.0f;
	if (u < 0.5f)
	{
		u = 2.0f*u;
	}
	else
	{
		u = 2.0f*(1 - u);
		sign = -1.0f;
	}

	AtFloat Cg = atanf((AI_PIOVER2 - phi_g)/gamma_g);
	AtFloat Dg = atanf(-phi_g/gamma_g);

	return sign * (gamma_g*tanf(u*(Cg - Dg) + Dg) + phi_g);
}

AtFloat NgPdf(AtFloat phi, AtFloat gamma_g, AtFloat phi_g)
{
	AtFloat Cg = atanf((AI_PIOVER2 - phi_g)/gamma_g);
	AtFloat Dg = atanf(-phi_g/gamma_g);
	AtFloat x = 0.5f*(Cg - Dg);
	AtFloat y = gamma_g / ((fabsf(phi) - phi_g)*(fabsf(phi) - phi_g) + gamma_g*gamma_g);
	return x*y;
}

shader_evaluate
{
	// Initialize result temporaries
	AtRGB result_diffuse_direct = AI_RGB_BLACK;
	AtRGB result_R_direct = AI_RGB_BLACK;
	AtRGB result_R_indirect = AI_RGB_BLACK;
	AtRGB result_TT_direct = AI_RGB_BLACK;
	AtRGB result_TT_indirect = AI_RGB_BLACK;
	AtRGB result_TRT_direct = AI_RGB_BLACK;
	AtRGB result_TRT_indirect = AI_RGB_BLACK;
	AtRGB result_TRTg_direct = AI_RGB_BLACK;
	AtRGB result_TRTg_indirect = AI_RGB_BLACK;

	// Get a local coordinate frame based on the hair fibre direction
	AtVector U = AiV3Normalize(sg->dPdv);
	AtVector up; AiV3Create(up, 0.0f, 1.0f, 0.0f);
	AtVector V = AiV3Cross(U, up);
	AtVector W = AiV3Cross(U, V);

	// Get the spherical angles of the exitant direction relative to the hair fibre
	AtVector wo = -sg->Rd;
	AtFloat theta_r = AI_PIOVER2 - sphericalTheta(wo, U);
	AtFloat phi_r = sphericalPhi(wo, V, W);

	// Get a random value per curve
	AtUInt32 curve_id = 0;
	AtFloat phi_offset = 0.0f;
	AiUDataGetUInt("curve_id", &curve_id);
	AtPoint2 p; p.x = AtFloat(curve_id); p.y = 0.0f;
	AtFloat cn = AiCellNoise2(p);
	// Map the random value to an offset for phi. This will be used by the glints to simulate the hair fibre twisting.
	phi_offset = cn * AI_PITIMES2;

	// Get shader data
	ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);

	// Get parameters
	AtFloat eta = 1.0f / AiShaderEvalParamFlt(p_ior);

	AtFloat beta_R = AiShaderEvalParamFlt(p_specularWidth) * AI_DTOR;
	AtFloat alpha_R = -AiShaderEvalParamFlt(p_specularShift) * AI_DTOR;

	AtFloat beta_TT = beta_R * 0.5f;
	AtFloat alpha_TT = alpha_R * 0.5f;

	AtFloat beta_TRT = beta_R * 2.0f;
	AtFloat alpha_TRT = alpha_R * 1.5f;

	AtFloat gamma_TT = AiShaderEvalParamFlt(p_transmissionRolloff) * AI_DTOR;
	AtFloat gamma_g = AiShaderEvalParamFlt(p_glintRolloff) * AI_DTOR;
	AtFloat phi_g = AiShaderEvalParamFlt(p_glintSeparation) * AI_DTOR;

	AtRGB specular1Color = AiShaderEvalParamRGB(p_specular1Color) * AiShaderEvalParamFlt(p_specular1Strength);
	AtRGB specular2Color = AiShaderEvalParamRGB(p_specular2Color) * AiShaderEvalParamFlt(p_specular2Strength);
	AtRGB transmissionColor = AiShaderEvalParamRGB(p_transmissionColor) * AiShaderEvalParamFlt(p_transmissionStrength);
	AtFloat glintStrength = AiShaderEvalParamFlt(p_glintStrength);

	// Direct lighting loop
	AiLightsPrepare(sg);
	while (AiLightsGetSample(sg))
	{
		// Get angle measures. See Section 3 in Ou et. al.
		AtFloat theta_i = AI_PIOVER2 - sphericalTheta(sg->Ld, U);
		AtFloat cos_theta_i = fabsf(cosf(theta_i));
		AtFloat phi_i = sphericalPhi(sg->Ld, V, W);
		AtFloat phi = phi_r - phi_i;
		if (phi < -AI_PI) phi += AI_PITIMES2;
		if (phi > AI_PI) phi -= AI_PITIMES2;
		AtFloat theta_h = (theta_r + theta_i)*0.5f;
		AtFloat theta_d = (theta_r - theta_i)*0.5f;
		AtFloat cos_theta_d = cosf(theta_d);
		AtFloat inv_cos_theta_d2 = 1.0f/(cos_theta_d*cos_theta_d);

		// fresnel
		AtFloat kr = fresnel(cos_theta_i, eta);

		// Calculate longitudinal and azimuthal functions. See Section 3.1 in Ou et. al.
		AtFloat Mr = g(beta_R, alpha_R, theta_h);
		AtFloat Nr = cosf(phi*0.5f);

		AtFloat Mtt = g(beta_TT, alpha_TT, theta_h);
		AtFloat Ntt = g(gamma_TT, 0.0f, AI_PI-phi);

		AtFloat Mtrt = g(beta_TRT, alpha_TRT, theta_h);
		AtFloat Ntrt = cosf(phi*0.5f);

		AtFloat phi_o = phi + phi_offset;
		if (phi_o < -AI_PI) phi_o += AI_PITIMES2;
		if (phi_o > AI_PI) phi_o -= AI_PITIMES2;
		AtFloat Ng = g(gamma_g, 0.0f, fabsf(phi_o) - phi_g);

		// Precalculate invariants across all lobes
		AtRGB L = sg->Li * sg->we * cos_theta_i * inv_cos_theta_d2 * AI_ONEOVER2PI;

		// Sum result temporaries for each lobe
		result_R_direct += L * Mr * Nr * kr;
		result_TT_direct += L * Mtt * Ntt * (1.0f-kr);
		result_TRT_direct += L * Mtrt * Ntrt * kr;
		result_TRTg_direct += L * Mtrt * Ng * kr;
	}

	// Multiply by user-defined reflectance
	result_R_direct *= specular1Color;
	result_TT_direct *= transmissionColor;
	result_TRT_direct *= specular2Color;
	result_TRTg_direct *= specular2Color * glintStrength;

	// Now sample each lobe
	double samples[2];
	AtRay wi_ray;
	AtFloat theta_i, phi, phi_i;
	AtFloat pdf;
	AtVector wi;
	AtScrSample scrs;
	// R
	{
		AtSamplerIterator* sampit = AiSamplerIterator(data->sampler_R, sg);
		AiMakeRay(&wi_ray, AI_RAY_GLOSSY, &sg->P, NULL, AI_BIG, sg);
		AtInt count=0;
		while(AiSamplerGetSample(sampit, samples))
		{
			theta_i = longitudinalSample(samples[0], theta_r, alpha_R, beta_R);
			phi = NrSample(samples[1]);
			phi_i = phi_r - phi;
			pdf = longitudinalPdf(theta_i, theta_r, alpha_R, beta_R) * NrPdf(phi);

			AtFloat theta_h = (theta_r + theta_i)*0.5f;

			sphericalDirection(theta_i, phi_i, V, W, U, wi);
			wi_ray.dir = wi;

			AiTrace(&wi_ray, &scrs);
			AtFloat Mr = g(beta_R, alpha_R, theta_h);
			AtFloat Nr = cosf(phi*0.5f);
			AtFloat theta_d = (theta_r - theta_i)*0.5f;
			AtFloat cos_theta_d = cosf(theta_d);
			AtFloat cos_theta_i = cosf(theta_i);
			AtFloat inv_cos_theta_d2 = 1.0f/(cos_theta_d*cos_theta_d);
			AtFloat kr = fresnel(cos_theta_i, eta);

			result_R_indirect += scrs.color*Mr*Nr*kr*cos_theta_i*inv_cos_theta_d2*AI_ONEOVER2PI / pdf;
		}
		if (count) result_R_indirect /= AtFloat(count);
		result_R_indirect *= specular1Color;
	}
	// TT
	{
		AtSamplerIterator* sampit = AiSamplerIterator(data->sampler_TT, sg);
		AiMakeRay(&wi_ray, AI_RAY_GLOSSY, &sg->P, NULL, AI_BIG, sg);
		AtInt count=0;
		while(AiSamplerGetSample(sampit, samples))
		{
			theta_i = longitudinalSample(samples[0], theta_r, alpha_TT, beta_TT);
			phi = NttSample(samples[1], gamma_TT);
			phi_i = phi_r - phi;
			pdf = longitudinalPdf(theta_i, theta_r, alpha_TT, beta_TT) * NttPdf(phi, gamma_TT);

			AtFloat theta_h = (theta_r + theta_i)*0.5f;

			sphericalDirection(theta_i, phi_i, V, W, U, wi);
			wi_ray.dir = wi;

			AiTrace(&wi_ray, &scrs);
			AtFloat Mtt = g(beta_TT, alpha_TT, theta_h);
			AtFloat Ntt = g(gamma_TT, 0.0f, AI_PI-phi);
			AtFloat theta_d = (theta_r - theta_i)*0.5f;
			AtFloat cos_theta_d = cosf(theta_d);
			AtFloat cos_theta_i = cosf(theta_i);
			AtFloat inv_cos_theta_d2 = 1.0f/(cos_theta_d*cos_theta_d);
			AtFloat kr = fresnel(fabsf(cos_theta_i), eta);

			result_TT_indirect += scrs.color*Mtt*Ntt*(1.0f-kr)*cos_theta_i*inv_cos_theta_d2*AI_ONEOVER2PI / pdf;
		}
		if (count) result_TT_indirect /= AtFloat(count);
		result_TT_indirect *= transmissionColor;
	}
	// TRT
	{
		AtSamplerIterator* sampit = AiSamplerIterator(data->sampler_TRT, sg);
		AiMakeRay(&wi_ray, AI_RAY_GLOSSY, &sg->P, NULL, AI_BIG, sg);
		AtInt count=0;
		while(AiSamplerGetSample(sampit, samples))
		{
			theta_i = longitudinalSample(samples[0], theta_r, alpha_TRT, beta_TRT);
			phi = NrSample(samples[1]);
			phi_i = phi_r - phi;
			pdf = longitudinalPdf(theta_i, theta_r, alpha_TRT, beta_TRT) * NrPdf(phi);

			AtFloat theta_h = (theta_r + theta_i)*0.5f;

			sphericalDirection(theta_i, phi_i, V, W, U, wi);
			wi_ray.dir = wi;

			AiTrace(&wi_ray, &scrs);
			AtFloat Mtrt = g(beta_TRT, alpha_TRT, theta_h);
			AtFloat Ntrt = cosf(phi*0.5f);
			AtFloat theta_d = (theta_r - theta_i)*0.5f;
			AtFloat cos_theta_d = cosf(theta_d);
			AtFloat cos_theta_i = cosf(theta_i);
			AtFloat inv_cos_theta_d2 = 1.0f/(cos_theta_d*cos_theta_d);
			AtFloat kr = fresnel(fabsf(cos_theta_i), eta);

			result_TRT_indirect += scrs.color*Mtrt*Ntrt*kr*cos_theta_i*inv_cos_theta_d2*AI_ONEOVER2PI / pdf;
		}
		if (count) result_TRT_indirect /= AtFloat(count);
		result_TRT_indirect *= specular2Color;
	}
	// TRTg
	{
		AtSamplerIterator* sampit = AiSamplerIterator(data->sampler_g, sg);
		AiMakeRay(&wi_ray, AI_RAY_GLOSSY, &sg->P, NULL, AI_BIG, sg);
		AtInt count=0;
		while(AiSamplerGetSample(sampit, samples))
		{
			theta_i = longitudinalSample(samples[0], theta_r, alpha_TRT, beta_TRT);
			phi = NgSample(samples[1], gamma_g, phi_g);
			phi_i = phi_r - phi;
			pdf = longitudinalPdf(theta_i, theta_r, alpha_TRT, beta_TRT) * NgPdf(phi, gamma_g, phi_g);

			AtFloat theta_h = (theta_r + theta_i)*0.5f;

			sphericalDirection(theta_i, phi_i, V, W, U, wi);
			wi_ray.dir = wi;

			AiTrace(&wi_ray, &scrs);
			AtFloat Mtrt = g(beta_TRT, alpha_TRT, theta_h);
			AtFloat phi_o = phi + phi_offset;
			if (phi_o < -AI_PI) phi_o += AI_PITIMES2;
			if (phi_o > AI_PI) phi_o -= AI_PITIMES2;
			AtFloat Ng = g(gamma_g, 0.0f, fabsf(phi_o) - phi_g);
			AtFloat theta_d = (theta_r - theta_i)*0.5f;
			AtFloat cos_theta_d = cosf(theta_d);
			AtFloat cos_theta_i = cosf(theta_i);
			AtFloat inv_cos_theta_d2 = 1.0f/(cos_theta_d*cos_theta_d);
			AtFloat kr = fresnel(fabsf(cos_theta_i), eta);

			result_TRTg_indirect += scrs.color*Mtrt*Ng*kr*cos_theta_i*inv_cos_theta_d2*AI_ONEOVER2PI / pdf;
		}
		if (count) result_TRTg_indirect /= AtFloat(count);
		result_TRTg_indirect *= specular2Color * glintStrength;
	}

	AiAOVSetRGB(sg, "specularDirect", result_R_direct);
	AiAOVSetRGB(sg, "specularIndirect", result_R_indirect);
	AiAOVSetRGB(sg, "specular2Direct", result_TRT_direct + result_TRTg_direct);
	AiAOVSetRGB(sg, "specular2Indirect", result_TRT_indirect + result_TRTg_indirect);
	AiAOVSetRGB(sg, "transmissionDirect", result_TT_direct);
	AiAOVSetRGB(sg, "transmissionIndirect", result_TT_indirect);

	sg->out.RGB = 	result_diffuse_direct +
					result_R_direct +
					result_R_indirect +
					result_TT_direct +
					result_TT_indirect +
					result_TRT_direct +
					result_TRT_indirect +
					result_TRTg_direct +
					result_TRTg_indirect;
}


