#include <ai.h>
#include <alUtil.h>

#define NUM_LAYERS 8

AI_SHADER_NODE_EXPORT_METHODS(alLayerColor)

const char* param_names[] = 
{
	"layer1",
	"layer1a",
	"layer2",
	"layer2a",
	"layer3",
	"layer3a",
	"layer4",
	"layer4a",
	"layer5",
	"layer5a",
	"layer6",
	"layer6a",
	"layer7",
	"layer7a",
	"layer8",
	"layer8a"
};

node_parameters
{
	AiParameterRGB("layer1", 0.0f, 0.0f, 0.0f);
	AiParameterFlt("layer1a", 0.0f);
	AiParameterRGB("layer2", 0.0f, 0.0f, 0.0f);
	AiParameterFlt("layer2a", 0.0f);
	AiParameterRGB("layer3", 0.0f, 0.0f, 0.0f);
	AiParameterFlt("layer3a", 0.0f);
	AiParameterRGB("layer4", 0.0f, 0.0f, 0.0f);
	AiParameterFlt("layer4a", 0.0f);
	AiParameterRGB("layer5", 0.0f, 0.0f, 0.0f);
	AiParameterFlt("layer5a", 0.0f);
	AiParameterRGB("layer6", 0.0f, 0.0f, 0.0f);
	AiParameterFlt("layer6a", 0.0f);
	AiParameterRGB("layer7", 0.0f, 0.0f, 0.0f);
	AiParameterFlt("layer7a", 0.0f);
	AiParameterRGB("layer8", 0.0f, 0.0f, 0.0f);
	AiParameterFlt("layer8a", 0.0f);
}

node_loader
{
	if (i>0) return false;
	node->methods     = alLayerColor;
	node->output_type = AI_TYPE_RGB;
	node->name        = "alLayerColor";
	node->node_type   = AI_NODE_SHADER;
	strcpy(node->version, AI_VERSION);
	return true;
}

struct ShaderData
{
	int max_layer;
};

node_initialize
{
	ShaderData* data = (ShaderData*) AiMalloc(sizeof(ShaderData));
	AiNodeSetLocalData(node, data);
}

node_finish
{
	ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);
	AiFree(data);
	AiNodeSetLocalData(node, NULL);
}

node_update
{
	ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);
	
	data->max_layer = -1;
	// check to see what the top layer should be
	for (int i=0; i < NUM_LAYERS; ++i)
	{
		// if the alpha is either linked or non-zero, layer is active
		if (AiNodeIsLinked(node, param_names[i*2+1]) ||
			params[i*2+1].FLT != 0.0f)
		{
			data->max_layer = i;
		}
	}

}

shader_evaluate
{
	ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);

	AtRGB result = AI_RGB_BLACK;

	for (int i=0; i <= data->max_layer; ++i)
	{
		AtRGB layerVal = AiShaderEvalParamRGB(i*2);
		float layerAlpha = AiShaderEvalParamFlt(i*2+1);
		result = lerp(result, layerVal, layerAlpha);
	}

	sg->out.RGB = result;
}