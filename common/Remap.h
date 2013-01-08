#pragma once

#include <alUtil.h>
#include <ai.h>

#define REMAP_FLOAT_PARAM_ENUM 	\
	p_inputMin,					\
	p_inputMax,					\
	p_contrast,					\
	p_contrastPivot,			\
	p_contrastSoftClip,			\
	p_biasVal,						\
	p_gainVal,						\
	p_outputMin,				\
	p_outputMax,				\
	p_clampEnable,				\
	p_threshold,				\
	p_clampMin,					\
	p_clampMax					\

#define REMAP_FLOAT_PARAM_DECLARE 				\
	AiParameterFLT("inputMin", 0.0f);			\
	AiParameterFLT("inputMax", 1.0f);			\
	AiParameterFLT("contrast", 1.0f);			\
	AiParameterFLT("contrastPivot", 0.5f);		\
	AiParameterFLT("contrastSoftClip", 0.0f);	\
	AiParameterFLT("bias", 0.5f);				\
	AiParameterFLT("gain", 0.5f);				\
	AiParameterFLT("outputMin", 0.0f);			\
	AiParameterFLT("outputMax", 1.0f);			\
	AiParameterBOOL("clampEnable", false);		\
	AiParameterBOOL("threshold", false);		\
	AiParameterFLT("clampMin", 0.0f);			\
	AiParameterFLT("clampMax", 1.0f);			\

struct RemapFloat
{
public:
	RemapFloat(AtFloat imn, AtFloat imx, AtFloat ct, AtFloat ctp, AtFloat ctsc, AtFloat bs, AtFloat gn, AtFloat omn, AtFloat omx,
				bool ce, bool t, AtFloat cmn, AtFloat cmx) :
		inputMin(imn),
		inputMax(imx),
		contrastVal(ct),
		contrastPivot(ctp),
		contrastSoftClip(ctsc),
		bias(bs),
		gain(gn),
		outputMin(omn),
		outputMax(omx),
		clampEnable(ce),
		threshold(t),
		clampMin(cmn),
		clampMax(cmx)
	{}

	AtFloat remap(AtFloat input)
	{
		AtFloat f = (input-inputMin)/(inputMax-inputMin);
		f = contrast(f, contrastVal, contrastPivot, contrastSoftClip);
		f = biasandgain(f, bias, gain);
		f = lerp(outputMin, outputMax, f);
		if (clampEnable)
		{
			f = std::min(clampMax, f);
			f = std::max(clampMin, f);
			if (threshold)
			{
				f = (f-clampMin)/(clampMax-clampMin);
			}
		}
		return f;
	}

	AtFloat inputMin;
	AtFloat inputMax;
	AtFloat contrastVal;
	AtFloat contrastPivot;
	AtFloat contrastSoftClip;
	AtFloat bias;
	AtFloat gain;
	AtFloat outputMin;
	AtFloat outputMax;
	bool clampEnable;
	bool threshold;
	AtFloat clampMin;
	AtFloat clampMax;
};

#define REMAP_FLOAT_CREATE							\
	RemapFloat( 									\
		AiShaderEvalParamFlt(p_inputMin), 			\
		AiShaderEvalParamFlt(p_inputMax),			\
		AiShaderEvalParamFlt(p_contrast),			\
		AiShaderEvalParamFlt(p_contrastPivot),		\
		AiShaderEvalParamFlt(p_contrastSoftClip),	\
		AiShaderEvalParamFlt(p_biasVal),				\
		AiShaderEvalParamFlt(p_gainVal),				\
		AiShaderEvalParamFlt(p_outputMin),			\
		AiShaderEvalParamFlt(p_outputMax),			\
		AiShaderEvalParamBool(p_clampEnable),		\
		AiShaderEvalParamBool(p_threshold),		\
		AiShaderEvalParamFlt(p_clampMin),			\
		AiShaderEvalParamFlt(p_clampMax)			\
	)
