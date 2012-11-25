#include "alUtil.h"

AtColor beckmannMicrofacetTransmission(AtShaderGlobals* sg, const AtVector& Z, const AtVector& X, const AtVector& Y,
										const AtVector& wo, AtSampler* sampler, AtFloat roughness, AtFloat eta,
										AtColor absorption)
{
	AtFloat count = 0.0f;
	double samples[2];
	AtFloat n1, n2;
	AtFloat kt;
	AtRay wi_ray;
	AtScrSample sample;
	AtVector wi, R;
	bool inside;
	AtColor result = AI_RGB_BLACK;
	float z = 0.0f;

	AtSamplerIterator* sampit = AiSamplerIterator(sampler, sg);

	AiMakeRay(&wi_ray, AI_RAY_REFRACTED, &sg->P, NULL, AI_BIG, sg);

	while (AiSamplerGetSample(sampit, samples))
	{
		// generate a microfacet normal, m
		// eq. 35,36
		AtFloat alpha2 = roughness*roughness;
		AtFloat tanThetaM = sqrtf(-alpha2 * logf(1.0f - AtFloat(samples[0])));
		AtFloat cosThetaM = 1.0f / sqrtf(1.0f + tanThetaM * tanThetaM);
		AtFloat sinThetaM = cosThetaM * tanThetaM;
		AtFloat phiM = 2.0f * AtFloat(AI_PI) * AtFloat(samples[1]);
		AtVector m = (cosf(phiM) * sinThetaM) * X +
				 	 (sinf(phiM) * sinThetaM) * Y +
							   	   cosThetaM  * Z;


		// get the refracted direction given m
		kt = 1.0f - fresnel(eta, m, wo, R, wi, inside);

		if (kt > IMPORTANCE_EPS) // if the final contribution is actually going to matter to the result
		{
			// eq. 33

			AtFloat cosThetaM2 = cosThetaM * cosThetaM;
			AtFloat tanThetaM2 = tanThetaM * tanThetaM;
			AtFloat cosThetaM4 = cosThetaM2 * cosThetaM2;
			AtFloat D = expf(-tanThetaM2 / alpha2) / (AtFloat(AI_PI) * alpha2 *  cosThetaM4);
			// eq. 24
			AtFloat pm = D * cosThetaM;
			// eval BRDF*cosNI
			AtFloat cosNI = AiV3Dot(Z, wi); // N.wi
			AtFloat cosNO = AiV3Dot(Z, wo);
			// eq. 26, 27: now calculate G1(i,m) and G1(o,m)
			AtFloat ao = 1 / (roughness * sqrtf((1.0f - cosNO * cosNO) / (cosNO * cosNO)));
			AtFloat ai = 1 / (roughness * sqrtf((1.0f - cosNI * cosNI) / (cosNI * cosNI)));
			AtFloat G1o = ao < 1.6f ? (3.535f * ao + 2.181f * ao * ao) / (1 + 2.276f * ao + 2.577f * ao * ao) : 1.0f;
			AtFloat G1i = ai < 1.6f ? (3.535f * ai + 2.181f * ai * ai) / (1 + 2.276f * ai + 2.577f * ai * ai) : 1.0f;
			AtFloat G = G1o * G1i;
			// eq. 21
			AtFloat cosHI = AiV3Dot(m, wi); // m.wi
			AtFloat cosHO = AiV3Dot(m, wo); // m.wo
			AtFloat Ht2 = eta * cosHI + cosHO;
			Ht2 *= Ht2;
			AtFloat brdf = (fabsf(cosHI * cosHO) * (eta * eta) * (G * D)) / fabsf(cosNO * Ht2);
			// eq. 38 and eq. 17
			AtFloat pdf = pm * (eta * eta) * fabsf(cosHI) / Ht2;

			wi_ray.dir = wi;
			AiTrace(&wi_ray, &sample);
			result += sample.color * kt * brdf/pdf;
			z += AtFloat(sample.z);
		}
		else if (AiV3IsZero(wi)) // total internal reflection
		{
			wi_ray.dir = R;
			AiTrace(&wi_ray, &sample);
			result += sample.color;
			z += AtFloat(sample.z);
		}

		count++;
	}

	result /= count;
	z /= count;

	if (maxh(absorption)>0.0f && !inside)
	{
			AtColor transmittance;
			transmittance.r = expf(absorption.r * -z);
			transmittance.g = expf(absorption.g * -z);
			transmittance.b = expf(absorption.b * -z);

			result *= transmittance;
	}

	return result;
}
