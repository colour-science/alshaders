#include "alUtil.h"
#include <iostream>

AtVector squareToCosineHemisphere(AtFloat u1, AtFloat u2)
{
   AtVector ret;
   concentricSampleDisk(u1, u2, ret.x, ret.y);
   ret.z = sqrtf(std::max(0.0f, 1.0f - ret.x*ret.x - ret.y*ret.y));
   return ret;
}

AtVector sample (const AtVector& N, const AtVector& wo, AtFloat roughness, AtFloat r1, AtFloat r2)
{
	AtFloat roughness_squared = roughness*roughness;
	AtFloat tanThetaM = sqrtf(-roughness_squared * logf(1 - r2));
	AtFloat cosThetaM = 1 / sqrtf(1 + tanThetaM * tanThetaM);
	AtFloat sinThetaM = cosThetaM * tanThetaM;
	AtFloat phiM = 2 * (float)AI_PI * r1;

	AtVector X; X.x = 1; X.y = 0; X.z = 0;
	AtVector Y; X.x = 0; Y.y = 1; Y.z = 0;
	AtVector Z; Z.x = 0; Z.y = 0; Z.z = 1;

	AtVector m =	X *(cosf(phiM) * sinThetaM) +
			Y * (sinf(phiM) * sinThetaM) +
			Z * cosThetaM;

	AtVector omega = m * 2 * AiV3Dot(m, wo) - wo;
	//AtVector eye = -wo;
	//AiReflect(&eye, &m, &omega);

	return omega;
}

 AtFloat pdf (const AtVector& N, const AtVector& wi, const AtVector& wo, AtFloat roughness)
{
	AtVector H = wi + wo;

	AiV3Normalize(H, H);

	AtFloat cosThetaM = AiV3Dot(N, H);
	if (cosThetaM <= 0)
	{
	  return 0.0;
	}
	else
	{
		AtFloat r2 = roughness*roughness;
		AtFloat cosThetaM2 = cosThetaM * cosThetaM;
		AtFloat tanThetaM2 = (1 - cosThetaM2) / cosThetaM2;
		AtFloat cosThetaM4 = cosThetaM2 * cosThetaM2;
		AtFloat D = expf(-tanThetaM2 /r2) / (AtFloat(AI_PI) * r2 *  cosThetaM4);
	return (D * cosThetaM * 0.25f) / AiV3Dot(H, wo);
	}
}

AtFloat brdf (const AtVector& N, const AtVector& wi, const AtVector& wo, AtFloat roughness, AtFloat eta)
{
    AtFloat cosNO = AiV3Dot(N, wo);
    AtFloat cosNI = AiV3Dot(N, wi);
    if (cosNO > 0 && cosNI > 0) {
		// get half vector
		AtVector H = wi + wo;
		AtVector Hr;
		AiV3Normalize(Hr, H);

		AtFloat cosThetaM = AiV3Dot(N, Hr);

		if (cosThetaM <=0)
		{
			return 0.0f;
		}
		else
		{
			AtFloat r2 = roughness*roughness;
			AtFloat cosThetaM2 = cosThetaM * cosThetaM;
			AtFloat tanThetaM2 = (1 - cosThetaM2) / cosThetaM2;
			AtFloat cosThetaM4 = cosThetaM2 * cosThetaM2;
			AtFloat D = expf(-tanThetaM2 / r2) / ((float)AI_PI * r2 *  cosThetaM4);

			AtFloat ao = 1 / (roughness * sqrtf((1 - cosNO * cosNO) / (cosNO * cosNO)));
			AtFloat ai = 1 / (roughness * sqrtf((1 - cosNI * cosNI) / (cosNI * cosNI)));
			AtFloat G1o = ao < 1.6f ? (3.535f * ao + 2.181f * ao * ao) / (1 + 2.276f * ao + 2.577f * ao * ao) : 1.0f;
			AtFloat G1i = ai < 1.6f ? (3.535f * ai + 2.181f * ai * ai) / (1 + 2.276f * ai + 2.577f * ai * ai) : 1.0f;
			AtFloat G = G1o * G1i;

			//AtFloat coshi = AiV3Dot(Hr, wi);
			//AtFloat kr = fresnel(coshi, eta);

			//return (G * D) * 0.25f * kr / cosNO;
			return (G * D) * 0.25f / cosNO;
		}
    }
	return 0.0f;
};

int main(int argc, char** argv)
{

	AtRGB tmp;
	AtVector N; N.x = 0; N.y = 0; N.z = 1;
	AtUInt32 ns = 1000000;

	AtUInt32 nos = 128;
	AtFloat* result = new AtFloat[nos];
	memset(result, 0, nos*sizeof(AtFloat));
	AtFloat step = 1.0f / nos;
	AtFloat y = 0;
	AtVector wo;
	wo.x = wo.y = wo.z = 0.0f;
	AtFloat roughness = 0.5f;
	AtFloat eta = 1.0f / 1.5f;
	for (AtUInt32 j = 0; j < nos; ++j)
	{
		wo.z = 1.0f - wo.y*wo.y;

		for (AtUInt32 i = 0; i < ns; ++i)
		{
			AtVector wi = sample(N, wo, roughness, drand48(), drand48());
			float p = pdf(N, wi, wo, roughness);
			if (p > 0.0f)
				result[j] += brdf(N, wi, wo, roughness, eta) / p;

		}
		result[j] /= AtFloat(ns);
		std::cout << wo.z << " " << result[j] << std::endl;
		wo.y += step;
	}
}
