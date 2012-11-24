#pragma once

#include <alUtil.h>
#include <ai.h>

#define REMAP_FLOAT_PARAM_ENUM 	\
	p_inputMin,					\
	p_inputMax,					\
	p_contrast,					\
	p_contrastPivot,			\
	p_contrastSoftClip,			\
	p_outputMin,				\
	p_outputMax					\

#define REMAP_FLOAT_PARAM_DECLARE 				\
	AiParameterFLT("inputMin", 0.0f);			\
	AiParameterFLT("inputMax", 1.0f);			\
	AiParameterFLT("contrast", 1.0f);			\
	AiParameterFLT("contrastPivot", 0.5f);		\
	AiParameterFLT("contrastSoftClip", 0.0f);	\
	AiParameterFLT("outputMin", 0.0f);			\
	AiParameterFLT("outputMax", 1.0f);			\

struct RemapFloat
{
public:
	RemapFloat(AtFloat imn, AtFloat imx, AtFloat ct, AtFloat ctp, AtFloat ctsc, AtFloat omn, AtFloat omx) :
		inputMin(imn),
		inputMax(imx),
		contrastVal(ct),
		contrastPivot(ctp),
		contrastSoftClip(ctsc),
		outputMin(omn),
		outputMax(omx)
	{}

	AtFloat remap(AtFloat input)
	{
		AtFloat f = (input-inputMin)/(inputMax-inputMin);
		f = contrast(f, contrastVal, contrastPivot, contrastSoftClip);
		f = lerp(outputMin, outputMax, f);
		return f;
	}

	AtFloat inputMin;
	AtFloat inputMax;
	AtFloat contrastVal;
	AtFloat contrastPivot;
	AtFloat contrastSoftClip;
	AtFloat outputMin;
	AtFloat outputMax;
};

#define REMAP_FLOAT_CREATE							\
	RemapFloat( 									\
		AiShaderEvalParamFlt(p_inputMin), 			\
		AiShaderEvalParamFlt(p_inputMax),			\
		AiShaderEvalParamFlt(p_contrast),			\
		AiShaderEvalParamFlt(p_contrastPivot),		\
		AiShaderEvalParamFlt(p_contrastSoftClip),	\
		AiShaderEvalParamFlt(p_outputMin),			\
		AiShaderEvalParamFlt(p_outputMax)			\
	)
