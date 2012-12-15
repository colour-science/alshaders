#include "IES.h"
#include "alUtil.h"
#include <ai.h>

AI_SHADER_NODE_EXPORT_METHODS(alIES)

struct ShaderData
{
	IESData* ies;
};

enum alIESParams
{
	p_filename
};

node_parameters
{
	AiParameterSTR("filename", "");
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alIES;
   node->output_type = AI_TYPE_RGB;
   node->name        = "alIES";
   node->node_type   = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return TRUE;
}

node_initialize
{
	ShaderData *data = (ShaderData*) AiMalloc(sizeof(ShaderData));
	data->ies = NULL;
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
	delete data->ies;
	data->ies = new IESData(params[p_filename].STR);
}

shader_evaluate
{
	// Get the light's xform matrix and its inverse
	AtMatrix lightToWorld, worldToLight;
	AiNodeGetMatrix(sg->Lp, "matrix", lightToWorld);
	AiM4Invert(lightToWorld, worldToLight);

	// get the major axis of the light in local space
	AtVector L;
	AiM4VectorByMatrixMult(&L, worldToLight, &sg->Ld);
	// just in case the user scaled the light
	L = AiV3Normalize(L);

	float phi = sphericalPhi(L);
	float theta = sphericalTheta(L);

	ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);

	if (data->ies->isValid())
	{
		float i = data->ies->lookup(phi, theta);
		sg->Liu *= AiColorCreate(i,i,i);

		//AtRGB red = AiColorCreate(1,0,0);
		//AtRGB green = AiColorCreate(0,1,0);
		//sg->Liu = AiColorLerp(i, red, green);
		//sg->Liu = AiColorCreate(i,i,i);

	}

}


