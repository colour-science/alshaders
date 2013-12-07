#include "Color.h"
#include "alUtil.h"
#include <ai.h>

AI_SHADER_NODE_EXPORT_METHODS(alBlackbody)

struct BlackbodySpectrum
{
	BlackbodySpectrum(const float temperature)
	: temp(temperature)
	{}

	float operator()(float wavelength) const
	{
		 double lambda = wavelength * 1e-9;
		 return (3.74183e-16 * pow(lambda, -5.0)) / (exp(1.4388e-2 / (lambda * temp)) - 1.0);
	}

	double temp;
};

struct ShaderData
{
	BlackbodySpectrum bs;
	AtRGB result;
};

enum alBlackbodyParams
{
	p_temperature,
	p_strength,
	p_physicalIntensity,
	p_physicalExposure
};

node_parameters
{
	AiParameterFLT("temperature", 1000.0f);
	AiParameterFLT("strength", 1.0f);
	AiParameterFLT("physicalIntensity", 1.0f);
	AiParameterFLT("physicalExposure", -16.0f);
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alBlackbody;
   node->output_type = AI_TYPE_RGB;
   node->name        = "alBlackbody";
   node->node_type   = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return true;
}

node_initialize
{
	ShaderData *data = (ShaderData*) AiMalloc(sizeof(ShaderData));
	AiNodeSetLocalData(node,data);
}

node_finish
{
	if (AiNodeGetLocalData(node))
	{
		ShaderData* data = (ShaderData*) AiNodeGetLocalData(node);
		AiFree((void*) data);
		AiNodeSetLocalData(node, NULL);
	}
}

node_update
{
	ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);
	if (!AiNodeIsLinked(node, "temperature"))
	{
		// temperature parameter is not connected, precalculate spectrum
		data->bs.temp = params[p_temperature].FLT;
		AtColor xyz = spectrumToXyz(data->bs);
		data->result = xyzToRgb(CsRec709, xyz);
	}
}



shader_evaluate
{
	float temperature = AiShaderEvalParamFlt(p_temperature);
	float physicalIntensity = AiShaderEvalParamFlt(p_physicalIntensity);
	float strength = AiShaderEvalParamFlt(p_strength);
	float exposure = AiShaderEvalParamFlt(p_physicalExposure);

	AtRGB result = AI_RGB_BLACK;

	ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);
	if (AiNodeIsLinked(node, "temperature"))
	{
		BlackbodySpectrum bs(temperature);
		AtColor xyz = spectrumToXyz(bs);
		result = xyzToRgb(CsRec709, xyz);
	}
	else
	{
		result = data->result;
	}

	result = max(AI_RGB_BLACK, result);

	if (physicalIntensity > 0.0f)
	{
		result *= lerp(1.0f, pow(temperature, 4) * 5.67e-8 * powf(2.0f, exposure), physicalIntensity);
	}

	sg->out.RGB = result * strength ;
}


