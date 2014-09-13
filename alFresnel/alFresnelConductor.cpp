#include <ai.h>
#include <cstring>
#include "fresnel.h"

AI_SHADER_NODE_EXPORT_METHODS(alFresnelConductor)

node_loader
{
	if (i>0) return 0;
	node->methods     = alFresnelConductor;
	node->output_type = AI_TYPE_RGB;
	node->name        = "alFresnelConductor";
	node->node_type   = AI_NODE_SHADER;
	strcpy(node->version, AI_VERSION);
	return true;
}

enum alFresnelConductorParams
{
	p_material=0,
	p_normalize,
	p_n,
	p_k
};

static const char* alFresnelConductorMaterialNames[] = 
{
	"aluminium",
	"chrome",
	"copper",
	"gold",
	"silver",
	"platinum",
	"titanium",
	"tungsten",
	"custom",
	NULL
};

node_parameters
{
	AiParameterEnum("material", 0, alFresnelConductorMaterialNames);
	AiParameterBool("normalize", false);
	AiParameterFlt("n", 1.19781);
	AiParameterFlt("k", 7.0488);

}

node_initialize
{

}

node_update
{

}

node_finish
{

}

shader_evaluate
{
	int material = AiShaderEvalParamInt(p_material);
	float n = AiShaderEvalParamFlt(n);
	float k = AiShaderEvalParamFlt(k);
	FresnelConductor fr;
	fr.setMaterial(material, n, k);
	sg->out.RGB = fr.kr(AiV3Dot(-sg->Rd, sg->Nf), 1.f);
}