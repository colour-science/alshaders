#include "alUtil.h"

int main(int argc, char** argv)
{

	AtUInt32 ns = 100000;
	AtRGB result = AI_RGB_BLACK;
	AtRGB tmp;
	for (AtUInt32 i=0; i < ns; ++i)
	{
		AtVector wi = squareToCosineHemisphere(drand48(), drand48());
		AtVector wo = squareToCosineHemisphere(drand48(), drand48());
		result += eval(wi, wo, drand48(), drand48(), tmp) / cosTheta(wo);
	}
}
