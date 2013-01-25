#include <ai.h>

AI_SHADER_NODE_EXPORT_METHODS(alInputVector)

enum alInputVectorParams
{
	p_input,
	p_type,
	p_userName,
	p_vector,
	p_matrix
};

enum Inputs
{
	IN_P = 0,
	IN_PO,
	IN_USER,
	IN_CUSTOM
};

static const char* InputNames[] = {
	"P",
	"Po",
	"User",
	"Custom",
	NULL
};

enum Types
{
	T_POINT = 0,
	T_VECTOR = 1
};

static const char* TypeNames[] =
{
	"Point",
	"Vector",
	NULL
};

node_parameters
{
	AiParameterEnum("input", IN_P, InputNames);
	AiParameterEnum("type", T_POINT, TypeNames);
	AiParameterStr("userName", "");
	AiParameterVec("vector", 0.0f, 0.0f, 0.0f);
	AtFloat mtx[4][4];
	AiM4Identity(mtx);
	AiParameterMtx("matrix", mtx);
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alInputVector;
   node->output_type = AI_TYPE_VECTOR;
   node->name        = "alInputVector";
   node->node_type   = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return TRUE;
}

node_initialize
{

}

node_finish
{

}

node_update
{

}

shader_evaluate
{
	int input = AiShaderEvalParamInt(p_input);
	int type = AiShaderEvalParamInt(p_type);
	const char* userName = AiShaderEvalParamStr(p_userName);
	AtVector vector = AiShaderEvalParamVec(p_vector);
	AtMatrix& mtx = *(AiShaderEvalParamMtx(p_matrix));
	AtVector result;

	// first select the input vector to use
	switch(input)
	{
	case IN_P:
		vector = sg->P;
		break;
	case IN_PO:
		vector = sg->Po;
		break;
	case IN_USER:
		AiUDataGetPnt(userName, &vector);
		break;
	default:
		break;
	}

	// now perform the appropriate xform based on the type
	if (type == T_POINT)
	{
		AiM4PointByMatrixMult(&result, mtx, &vector);
	}
	else
	{
		AiM4VectorByMatrixMult(&result, mtx, &vector);
	}

	sg->out.VEC = result;
}


