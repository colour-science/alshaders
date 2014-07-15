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
	p_material=0
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
	NULL
};

node_parameters
{
	AiParameterEnum("material", 0, alFresnelConductorMaterialNames);
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
	FresnelConductor fr;
	fr.setMaterial(material);
	sg->out.RGB = fr.kr(AiV3Dot(-sg->Rd, sg->Nf));
}