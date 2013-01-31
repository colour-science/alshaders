// Hair shader based on ISHair: Importance Sampling for Hair Scattering by Ou et. al 2012
// http://www.cs.dartmouth.edu/~ouj/site/Research/Entries/2012/6/21_ISHair__Importance_Sampling_for_Hair_Scattering.html

#include <ai.h>
#include "alUtil.h"
#include <vector>
#include <algorithm>

#define SAMPLE_ALL 0

struct ShaderData
{
	ShaderData()
	: sampler_diffuse(NULL), sampler_R(NULL), sampler_TT(NULL), sampler_TRT(NULL), sampler_g(NULL)
	{}
	~ShaderData()
	{
		AiSamplerDestroy(sampler_diffuse);
		AiSamplerDestroy(sampler_R);
		AiSamplerDestroy(sampler_TT);
		AiSamplerDestroy(sampler_TRT);
		AiSamplerDestroy(sampler_g);
	}

	void update(int diffuse_samples, int glossy_samples)
	{
		AiSamplerDestroy(sampler_diffuse);
		sampler_diffuse = AiSampler(diffuse_samples, 2);
		AiSamplerDestroy(sampler_R);
		sampler_R = AiSampler(glossy_samples, 2);
		AiSamplerDestroy(sampler_TT);
		sampler_TT = AiSampler(glossy_samples, 2);
		AiSamplerDestroy(sampler_TRT);
		sampler_TRT = AiSampler(glossy_samples, 2);
		AiSamplerDestroy(sampler_g);
		sampler_g = AiSampler(glossy_samples, 2);
	}

	AtSampler* sampler_diffuse;
	AtSampler* sampler_R;
	AtSampler* sampler_TT;
	AtSampler* sampler_TRT;
	AtSampler* sampler_g;
};

AI_SHADER_NODE_EXPORT_METHODS(alHair)

enum alHairParams
{
	p_ior,
	p_diffuseStrength,
	p_diffuseColor,
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
	AiParameterFlt("diffuseStrength", 1.0f);
	AiParameterRGB("diffuseColor", 0.31f, 0.08f, 0.005f);
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
	int diffuse_samples = AiNodeGetInt(options, "GI_diffuse_samples");
	int glossy_samples = AiNodeGetInt(options, "GI_glossy_samples");
	data->update(diffuse_samples, glossy_samples);

	AiAOVRegister("specularDirect", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
	AiAOVRegister("specularIndirect", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
}

AtFloat g(AtFloat beta, AtFloat alpha, AtFloat theta_h)
{
	AtFloat n = theta_h-alpha;
	return expf(-(n*n)/(2.0f*beta*beta));
}

void AB(AtFloat theta_r, AtFloat alpha, AtFloat beta, AtFloat& A, AtFloat& B)
{
	A = atanf(((AI_PIOVER2+theta_r)*0.5f - alpha) / beta);
	B = atanf(((theta_r-AI_PIOVER2)*0.5f - alpha) / beta);
}

AtFloat sampleLong(double u, AtFloat theta_r, AtFloat alpha, AtFloat beta)
{
	AtFloat A, B;
	AB(theta_r, alpha, beta, A, B);

	AtFloat t = beta * tanf(u*(A-B) + B);
	AtFloat theta_h = t + alpha;
	return clamp( -0.4999f * AI_PI, 0.4999f * AI_PI, (2.0f*theta_h - theta_r));
}


void sampleSphere(double u[2], const AtVector& U, const AtVector& V, const AtVector& W, AtFloat& theta_i, AtFloat& phi_i)
{
	AtVector w = uniformSampleSphere(u[0], u[1]);
	theta_i = AI_PIOVER2 - sphericalTheta(w, U);
	phi_i = sphericalPhi(w, V, W);
}

AtFloat spherePdf()
{
	return AI_ONEOVERPI * 0.25f;
}

void sample_R(double u[2], AtFloat theta_r, AtFloat phi_r, AtFloat alpha, AtFloat beta,
				AtFloat& theta_i, AtFloat& phi_i, AtFloat& brdf, const AtVector& U, const AtVector& V, const AtVector& W)
{

	theta_i = sampleLong(u[0], theta_r, alpha, beta);
	AtFloat theta_h = (theta_i+theta_r)*0.5f;

	AtFloat phi = 2.0f*asinf(1.0f - 2.0f*u[1]);
	phi_i = phi_r - phi;

	/*
	sampleSphere(u, U, V, W, theta_i, phi_i);
	AtFloat theta_h = (theta_i+theta_r)*0.5f;
	AtFloat phi = phi_r - phi_i;
	*/
	if (phi < -AI_PI) phi += AI_PITIMES2;
	if (phi > AI_PI) phi -= AI_PITIMES2;

	AtFloat Mr = g(beta, alpha, theta_h);
	AtFloat Nr = cosf(phi*0.5f);
	AtFloat theta_d = (theta_r - theta_i)*0.5f;
	AtFloat cos_theta_d = cosf(theta_d);
	AtFloat cos_theta_i = cosf(theta_i);
	AtFloat inv_cos_theta_d2 = std::max(0.001f, 1.0f/(cos_theta_d*cos_theta_d));
	AtFloat kr = 1;//fresnel(cos_theta_i, eta);

	brdf = Mr*Nr*kr*cos_theta_i*inv_cos_theta_d2*AI_ONEOVER2PI;
}

void sample_TT(double u[2], AtFloat theta_r, AtFloat phi_r, AtFloat alpha, AtFloat beta, AtFloat gamma,
				AtFloat& theta_i, AtFloat& phi_i, AtFloat& brdf, const AtVector& U, const AtVector& V, const AtVector& W)
{

	theta_i = sampleLong(u[0], theta_r, alpha, beta);
	AtFloat c = 2.0f * atanf(AI_PIOVER2 / gamma);
	AtFloat p = gamma * tanf((u[1] - 0.5f)*c);
	AtFloat phi = p + AI_PI;
	phi_i = phi_r - phi;

	/*
	sampleSphere(u, U, V, W, theta_i, phi_i);
	AtFloat phi = phi_r - phi_i;
	*/
	if (phi < -AI_PI) phi += AI_PITIMES2;
	if (phi > AI_PI) phi -= AI_PITIMES2;

	AtFloat theta_h = (theta_i+theta_r)*0.5f;
	AtFloat Mtt = g(beta, alpha, theta_h);
	AtFloat Ntt = g(gamma, 0.0f, AI_PI-phi);
	AtFloat theta_d = (theta_r - theta_i)*0.5f;
	AtFloat cos_theta_d = cosf(theta_d);
	AtFloat cos_theta_i = fabs(cosf(theta_i));
	AtFloat inv_cos_theta_d2 = std::max(0.001f, 1.0f/(cos_theta_d*cos_theta_d));
	AtFloat kt = 1;//1.0f - fresnel(fabsf(cos_theta_i), eta);

	brdf = Mtt*Ntt*kt*cos_theta_i*inv_cos_theta_d2*AI_ONEOVER2PI;
}


void sample_g(double u[2], AtFloat theta_r, AtFloat phi_r, AtFloat alpha, AtFloat beta, AtFloat gamma, AtFloat phi_g,
				AtFloat phi_offset, AtFloat& theta_i, AtFloat& phi_i, AtFloat& brdf, const AtVector& U, const AtVector& V, const AtVector& W)
{
	theta_i = sampleLong(u[0], theta_r, alpha, beta);
	AtFloat c = atanf((AI_PIOVER2-phi_g)/gamma);
	AtFloat d = atanf(-phi_g/gamma);

/*
	AtFloat sign;
	if (u[1] < 0.5)
	{
		u[1] = 2.0f*u[1];
		sign = 1.0f;
	}
	else
	{
		u[1] = 2.0f*(1.0f-u[1]);
		sign = -1.0f;
	}

	AtFloat p = gamma * tanf(u[1]*(c-d) + d);
	AtFloat phi = sign * (p + phi_g);
	phi_i = phi_r - phi;
	if (phi_i < -AI_PI) phi_i += AI_PITIMES2;
	if (phi_i > AI_PI) phi_i -= AI_PITIMES2;


*/
	// FIXME: g sampling broken. use trt sampling instead
	AtFloat phi = 2.0f*asinf(1.0f - 2.0f*u[1]);
	phi_i = phi_r - phi;


	/*
	sampleSphere(u, U, V, W, theta_i, phi_i);
	phi = phi_r - phi_i;
	*/
	if (phi < -AI_PI) phi += AI_PITIMES2;
	if (phi > AI_PI) phi -= AI_PITIMES2;

	AtFloat theta_h = (theta_i+theta_r)*0.5f;
	AtFloat Mtrt = g(beta, alpha, theta_h);
	AtFloat Ng = g(gamma, 0.0f, fabsf(phi) - phi_g);
	AtFloat theta_d = (theta_r - theta_i)*0.5f;
	AtFloat cos_theta_d = cosf(theta_d);
	AtFloat cos_theta_i = cosf(theta_i);
	AtFloat inv_cos_theta_d2 = std::max(0.001f, 1.0f/(cos_theta_d*cos_theta_d));
	AtFloat kr = 1;//fresnel(fabsf(cos_theta_i), eta);

	brdf = Mtrt*Ng*kr*cos_theta_i*inv_cos_theta_d2*AI_ONEOVER2PI;
	//brdf = kr*cos_theta_i*inv_cos_theta_d2*AI_ONEOVER2PI;;
	//brdf = 1;
}

AtFloat pdfLong(AtFloat theta_r, AtFloat theta_i, AtFloat alpha, AtFloat beta)
{
	AtFloat A, B;
	AB(theta_r, alpha, beta, A, B);
	AtFloat theta_h = (theta_r+theta_i)*0.5f;
	AtFloat t = theta_h - alpha;
	return beta / ((t*t + beta*beta) * std::max(0.001f, (2.0f * (A-B) * cosf(theta_i))));
	//return 1.0f;
	//return beta / (t*t);
	//return beta / fabsf(t*t);
	//return t*t;
	//return beta;
	//return fabsf(beta / (t*t));
}

AtFloat pdf_R(AtFloat theta_r, AtFloat phi_r, AtFloat theta_i, AtFloat phi_i, AtFloat alpha, AtFloat beta)
{
	AtFloat phi_pdf = cosf((phi_r - phi_i) * 0.5f) * 0.25f;
	return phi_pdf * pdfLong(theta_r, theta_i, alpha, beta);
	//return 1;
}

AtFloat pdf_TT(AtFloat theta_r, AtFloat phi_r, AtFloat theta_i, AtFloat phi_i, AtFloat alpha, AtFloat beta, AtFloat gamma)
{
	AtFloat c = 2.0f * atanf(AI_PIOVER2 / gamma);
	AtFloat phi = fabsf(phi_r - phi_i);
	if (phi < -AI_PI) phi += AI_PITIMES2;
	if (phi > AI_PI) phi -= AI_PITIMES2;
	AtFloat p = phi + AI_PI;
	AtFloat phi_pdf = (gamma / (p*p + gamma*gamma)) / c;
	return phi_pdf * pdfLong(theta_r, theta_i, alpha, beta);
	//return 1;
}

AtFloat pdf_g(AtFloat theta_r, AtFloat phi_r, AtFloat theta_i, AtFloat phi_i, AtFloat alpha, AtFloat beta,
				AtFloat gamma, AtFloat phi_g)
{
	AtFloat c = atanf((AI_PIOVER2-phi_g)/gamma);
	AtFloat d = atanf(-phi_g/gamma);
	AtFloat phi = fabsf(phi_r - phi_i);
	AtFloat p = phi - phi_g;
	AtFloat phi_pdf = gamma / (p*p + gamma*gamma) / (2.0f*fabsf((c-d)));
	return phi_pdf * pdfLong(theta_r, theta_i, alpha, beta);
	//return 1.0f;
}

shader_evaluate
{
	// Initialize result temporaries
	AtRGB result_diffuse_direct = AI_RGB_BLACK;
	AtRGB result_diffuse_indirect = AI_RGB_BLACK;
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
	//AtFloat phi_g = AiShaderEvalParamFlt(p_glintSeparation) * AI_DTOR;
	AtFloat phi_g = lerp(30.0f*AI_DTOR, 45.0f*AI_DTOR, cn);

	AtRGB diffuseColor = AiShaderEvalParamRGB(p_diffuseColor) * AiShaderEvalParamFlt(p_diffuseStrength);
	AtRGB specular1Color = AiShaderEvalParamRGB(p_specular1Color) * AiShaderEvalParamFlt(p_specular1Strength);
	AtRGB specular2Color = AiShaderEvalParamRGB(p_specular2Color) * AiShaderEvalParamFlt(p_specular2Strength);
	AtRGB transmissionColor = AiShaderEvalParamRGB(p_transmissionColor) * AiShaderEvalParamFlt(p_transmissionStrength);
	AtFloat glintStrength = AiShaderEvalParamFlt(p_glintStrength);

	bool do_glossy = true;
	bool do_R = true, do_TT = true, do_TRT = true, do_g = true;
	if (sg->Rr > 0)
	{
		do_R = true;
		do_TT = true;
		do_TRT = true;
		do_g = false;
	}
	else if (sg->Rr > 0)
	{
		do_glossy = false;
	}

	// Direct lighting loop
	// Tell Arnold we want the full sphere for lighting.
	sg->fhemi = false;
	AiLightsPrepare(sg);
	while (AiLightsGetSample(sg))
	{
		// Diffuse
		AtFloat tl = AiV3Dot(sg->Ld, U);
		result_diffuse_direct += sqrtf(1.0f - tl*tl) * sg->Li * sg->we * AI_ONEOVER2PI;

		if (do_glossy)
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
			AtFloat inv_cos_theta_d2 = std::max(0.001f, 1.0f/(cos_theta_d*cos_theta_d));

			// Precalculate invariants across all lobes
			AtRGB L = sg->Li * sg->we * cos_theta_i * inv_cos_theta_d2 * AI_ONEOVER2PI;

			if (maxh(L) > IMPORTANCE_EPS)
			{
				// fresnel
				AtFloat kr = 1.0f;//fresnel(cos_theta_i, eta);
				AtFloat kt = 1.0f;

				// Calculate longitudinal and azimuthal functions. See Section 3.1 in Ou et. al.
				AtFloat Mr = g(beta_R, alpha_R, theta_h);
				AtFloat Nr = cosf(phi*0.5f);

				AtFloat Mtt = g(beta_TT, alpha_TT, theta_h);
				AtFloat Ntt = g(gamma_TT, 0.0f, AI_PI-phi);

				AtFloat Mtrt = g(beta_TRT, alpha_TRT, theta_h);
				AtFloat Ntrt = cosf(phi*0.5f);

				AtFloat Ng = g(gamma_g, 0.0f, fabsf(phi) - phi_g);

				// Sum result temporaries for each lobe
				if (do_R) result_R_direct += L * Mr * Nr * kr;
				if (do_TT) result_TT_direct += L * Mtt * Ntt * kt;
				if (do_TRT) result_TRT_direct += L * Mtrt * Ntrt * kr;
				if (do_g) result_TRTg_direct += L * Mtrt * Ng * kr;
			}
		}
	}

	// reset this.
	sg->fhemi = true;

	// Multiply by user-defined reflectance
	result_diffuse_direct *= diffuseColor;
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
	// Diffuse
	{
		AtSamplerIterator* sampit = AiSamplerIterator(data->sampler_diffuse, sg);
		AiMakeRay(&wi_ray, AI_RAY_DIFFUSE, &sg->P, NULL, AI_BIG, sg);
		AtInt count=0;
		while(AiSamplerGetSample(sampit, samples))
		{
			wi = uniformSampleSphere(samples[0], samples[1]);
			wi_ray.dir = wi;
			AiTrace(&wi_ray, &scrs);
			AtFloat tl = AiV3Dot(wi, U);
			AtFloat kr = fresnel(fabsf(tl), eta);
			result_diffuse_indirect += sqrtf(1.0f - tl*tl) * scrs.color * AI_ONEOVER2PI;
		}
		result_diffuse_indirect *= AiSamplerGetSampleInvCount(sampit);
		result_diffuse_indirect *= diffuseColor;
	}

	if (do_glossy)
	{
		AtSamplerIterator* sampit = AiSamplerIterator(data->sampler_R, sg);
		AiMakeRay(&wi_ray, AI_RAY_GLOSSY, &sg->P, NULL, AI_BIG, sg);
		AtFloat theta_i, phi_i;
		AtFloat brdf;
		AtRGB* sum;
		AtRGB I;
		AtRGB tmp;
		while(AiSamplerGetSample(sampit, samples))
		{
			if (samples[0] < 0.5 && samples[1] < 0.5) // R
			{
				samples[0] = 2.0 * samples[0];
				samples[1] = 2.0 * samples[1];
				sample_R(samples, theta_r, phi_r, alpha_R, beta_R, theta_i, phi_i, brdf, U, V, W);
				if (do_R ) sum = &result_R_indirect; else sum = &tmp;
				I = specular1Color;
			}
			else if (samples[0] >= 0.5 && samples[1] < 0.5) // TT
			{
				samples[0] = 2.0 * (1.0 - samples[0]);
				samples[1] = 2.0 * samples[1];
				sample_TT(samples, theta_r, phi_r, alpha_TT, beta_TT, gamma_TT, theta_i, phi_i, brdf, U, V, W);
				sum = &result_TT_indirect;
				if (do_TT) sum = &result_TT_indirect; else sum = &tmp;
				I = transmissionColor;
			}
			else if (samples[0] < 0.5 && samples[1] >= 0.5) // TRT
			{
				samples[0] = 2.0 * samples[0];
				samples[1] = 2.0 * (1.0 - samples[1]);
				sample_R(samples, theta_r, phi_r, alpha_TRT, beta_TRT, theta_i, phi_i, brdf, U, V, W);
				sum = &result_TRT_indirect;
				if (do_TRT) sum = &result_TRT_indirect; else sum = &tmp;
				I = specular2Color;
			}
			else // G
			{
				samples[0] = 2.0 * (1.0 - samples[0]);
				samples[1] = 2.0 * (1.0 - samples[1]);
				sample_g(samples, theta_r, phi_r, alpha_TRT, beta_TRT, gamma_g, phi_g, phi_offset, theta_i, phi_i, brdf, U, V, W);
				sum = &result_TRTg_indirect;
				if (do_g) sum = &result_TRTg_indirect; else sum = &tmp;
				I = specular2Color * glintStrength;
			}

			sphericalDirection(theta_i, phi_i, V, W, U, wi);
			wi_ray.dir = wi;
			AtFloat pdf = 	pdf_R(theta_r, phi_r, theta_i, phi_i, alpha_R, beta_R)
						+ 	pdf_TT(theta_r, phi_r, theta_i, phi_i, alpha_TT, beta_TT, gamma_TT)
						+	pdf_R(theta_r, phi_r, theta_i, phi_i, alpha_TRT, beta_TRT)
						//+ 	pdf_g(theta_r, phi_r, theta_i, phi_i, alpha_TRT, beta_TRT, phi_g, gamma_g);
						+	pdf_R(theta_r, phi_r, theta_i, phi_i, alpha_TRT, beta_TRT);
			pdf *= 0.25f;
			AiTrace(&wi_ray, &scrs);
			(*sum) += (scrs.color * brdf * I * AI_PI) / (pdf);
		}
		AtFloat in = AiSamplerGetSampleInvCount(sampit);
		result_R_indirect *= in;
		result_TT_indirect *= in;
		result_TRT_indirect *= in;
		result_TRTg_indirect *= in;
	}

	if (sg->Rt & AI_RAY_CAMERA)
	{
		AiAOVSetRGB(sg, "diffuseDirect", result_diffuse_direct);
		AiAOVSetRGB(sg, "diffuseIndirect", result_diffuse_indirect);
		AiAOVSetRGB(sg, "specularDirect", result_R_direct);
		AiAOVSetRGB(sg, "specularIndirect", result_R_indirect);
		AiAOVSetRGB(sg, "specular2Direct", result_TRT_direct);
		AiAOVSetRGB(sg, "specular2Indirect", result_TRT_indirect);
		AiAOVSetRGB(sg, "glintDirect", result_TRTg_direct);
		AiAOVSetRGB(sg, "glintIndirect", result_TRTg_indirect);
		AiAOVSetRGB(sg, "transmissionDirect", result_TT_direct);
		AiAOVSetRGB(sg, "transmissionIndirect", result_TT_indirect);
	}

	sg->out.RGB = 	result_diffuse_direct +
					result_diffuse_indirect +
					result_R_direct +
					result_R_indirect +
					result_TT_direct +
					result_TT_indirect +
					result_TRT_direct +
					result_TRT_indirect +
					result_TRTg_direct +
					result_TRTg_indirect;
}


