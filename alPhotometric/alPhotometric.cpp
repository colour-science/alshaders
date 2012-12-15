#include "Photometric.h"
#include "alUtil.h"
#include <ai.h>

AI_SHADER_NODE_EXPORT_METHODS(alPhotometric)

struct ShaderData
{
	PhotometricData* photometric;
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
   node->methods     = alPhotometric;
   node->output_type = AI_TYPE_RGB;
   node->name        = "alPhotometric";
   node->node_type   = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return TRUE;
}

node_initialize
{
	ShaderData *data = (ShaderData*) AiMalloc(sizeof(ShaderData));
	data->photometric = NULL;
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
	delete data->photometric;
	data->photometric = new PhotometricData(params[p_filename].STR);
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

	if (data->photometric->isValid())
	{
		float i = data->photometric->lookup(phi, theta);
		sg->Liu *= AiColorCreate(i,i,i);

	}

}


