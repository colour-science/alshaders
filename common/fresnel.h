#pragma once
#include <ai.h>
#include "alUtil.h"

class Fresnel
{
public:
	virtual ~Fresnel(){};
	virtual AtRGB kr(float cos_theta)=0;

	float _eta;
};

class FresnelDielectric : public Fresnel
{
public:
	FresnelDielectric(float eta) { _eta = eta; }
	virtual ~FresnelDielectric(){}

	virtual AtRGB kr(float cos_theta)
	{
		return rgb(fresnel(cos_theta, _eta));
	}
};

#define FRCOND_STEPS 90
class FresnelConductor : public Fresnel
{
public:
	FresnelConductor();
	virtual ~FresnelConductor(){}
	void setMaterial(int material);
	virtual AtRGB kr(float cos_theta);
private:
	float* _data;
};

enum ConductorMaterial
{
	kAluminium=0,
	kChrome,
	kCopper,
	kGold,
	kSilver,
	kPlatinum,
	kTitanium,
	kTungsten	
};