#include <ai.h>
#include <string>
#include <map>
#include <algorithm>
#include <vector>
#include "filters.h"


///////////////////////////////////////////////
//
//		Filter node things
//
///////////////////////////////////////////////

AI_FILTER_NODE_EXPORT_METHODS(cryptomatte_filterMtd)

enum cryptomatte_filterParams {
	p_width,
	p_rank,
	p_filter,
	p_mode,
};


enum modeEnum {
	p_mode_double_rgba,
};

static const char* modeEnumNames[] = {
	"double_rgba",
	NULL
};

node_parameters {
	AiMetaDataSetStr(mds, NULL, "maya.attr_prefix", "filter_");
	AiMetaDataSetStr(mds, NULL, "maya.translator", "cryptomatteFilter");
	AiMetaDataSetInt(mds, NULL, "maya.id", 0x00116420);

	AiParameterFlt("width", 2.0);
	AiParameterInt("rank", 0);
	AiParameterEnum("filter", p_filter_gaussian, filterEnumNames);
	AiParameterEnum("mode", p_mode_double_rgba, modeEnumNames);
}

node_loader {
	if (i>0) return 0;
	node->methods     = cryptomatte_filterMtd;
	node->output_type = AI_TYPE_RGBA;
	node->name        = "cryptomatte_filter";
	node->node_type   = AI_NODE_FILTER;
	strcpy(node->version, AI_VERSION);
	return true;
}

node_initialize {
	static const char* necessary_aovs[] = {
		"FLOAT Z", 
		NULL 
	};
	AiFilterInitialize(node, true, necessary_aovs, NULL);
}

node_finish {
	AiFilterDestroy(node);
}

node_update {
	AtShaderGlobals shader_globals;
	AtShaderGlobals *sg = &shader_globals;
	if (AiShaderEvalParamEnum(p_filter) == p_filter_box) {
		AiFilterUpdate(node, 1.0f);
	} else {		
		float width = AiShaderEvalParamFlt(p_width);
		AiFilterUpdate(node, width);
	}
}

filter_output_type {
	return AI_TYPE_RGBA;
}


///////////////////////////////////////////////
//
//		Sample-Weight map Class and type definitions
//
///////////////////////////////////////////////

class compareTail {
public:
	bool operator()(const std::pair<float, float> x, const std::pair<float, float> y) {
		return x.second > y.second;
	}
};

typedef std::map<float,float> 			sw_map_type ;
typedef std::map<float,float>::iterator sw_map_iterator_type;

void write_to_samples_map(sw_map_type * vals, float hash, float sample_weight) {
   (*vals)[hash] += sample_weight;
}

///////////////////////////////////////////////
//
//		Filter proper
//
///////////////////////////////////////////////

filter_pixel {	
	AtShaderGlobals shader_globals;
	AtShaderGlobals *sg = &shader_globals;
	AtRGBA *out_value = (AtRGBA *)data_out;
	*out_value = AI_RGBA_BLACK;

	int rank =        AiShaderEvalParamInt (p_rank);
	float width =     AiShaderEvalParamFlt (p_width);
	int main_filter = AiShaderEvalParamEnum(p_filter);

	AtNode * renderOptions = AiUniverseGetOptions();
	int auto_transparency_depth = AiNodeGetInt(renderOptions, "auto_transparency_depth");

	AtNode * camera = AiUniverseGetCamera();
	float camera_far_clip = AiNodeGetFlt(camera, "far_clip");
	
	AtRGBA val = AI_RGBA_BLACK;

	float (*filter)(AtPoint2, float);


	///////////////////////////////////////////////
	//
	//		early out for black pixels
	//
	///////////////////////////////////////////////
	
	bool early_out = true;

	while (AiAOVSampleIteratorGetNext(iterator)) {

		if (AiAOVSampleIteratorHasValue(iterator))	{
			early_out = false;	
			break;
		}

	}

	if (early_out) {
		*out_value = AI_RGBA_BLACK;
		if (rank == 0) {
			out_value->g = 1.0f;					
		}
		return;
	}
	AiAOVSampleIteratorReset(iterator);


	///////////////////////////////////////////////
	//
	//		Select Filter
	//
	///////////////////////////////////////////////

	switch (main_filter) {
		case p_filter_triangle:
			filter = &triangle;
			break;
		case p_filter_blackman_harris:
			filter = &blackman_harris;
			break;
		case p_filter_box:
			filter = &box;
			break;
		case p_filter_disk:
			filter = &disk;
			break;
		case p_filter_cone:
			filter = &cone;
			break;
		case p_filter_gaussian:
		default:
			filter = &gaussian;
			break;
	}


	///////////////////////////////////////////////
	//
	//		Set up sample-weight maps and friends
	//
	///////////////////////////////////////////////

	sw_map_type vals;


	int total_samples = 0;
	float total_weight = 0.0f;
	float sample_weight;
	AtPoint2 offset;
	AtRGBA color;

	///////////////////////////////////////////////
	//
	//		Iterate samples
	//
	///////////////////////////////////////////////

	while (AiAOVSampleIteratorGetNext(iterator)) {
		offset = AiAOVSampleIteratorGetOffset(iterator);
		sample_weight = filter(offset, width);

		if (sample_weight == 0.0f) {
			continue;
		}
		total_samples++;

		///////////////////////////////////////////////
		//
		//		Samples with no value.
		//
		///////////////////////////////////////////////

		if (!AiAOVSampleIteratorHasValue(iterator))	{			
			total_weight += sample_weight;
			write_to_samples_map(&vals, 0.0f, sample_weight);
			continue;
		}


		///////////////////////////////////////////////
		//
		//		Samples with values, depth or no depth
		//
		///////////////////////////////////////////////

		float quota = sample_weight;
		float iterative_transparency_weight = 1.0f;
		int depth_counter = 0;
		float z = 0.0f;
		bool traced_bg = false;

		color = AiAOVSampleIteratorGetRGBA(iterator);
		while (AiAOVSampleIteratorGetNextDepth(iterator)) {
			if (depth_counter > 0)
				color = AiAOVSampleIteratorGetRGBA(iterator);
			depth_counter++;
			z = AiAOVSampleIteratorGetAOVFlt(iterator, "Z");

			float sub_sample_opacity = color.b;
			color.b = 0.0f;

			// so if our sample is 80% opaque, and iterative transparency indicates that 10% opacity is left over, we'll get 8% of that sample. 
			float sub_sample_weight = sub_sample_opacity * iterative_transparency_weight * sample_weight;  
			iterative_transparency_weight *= (1.0f - sub_sample_opacity); // so if the current sub sample is 80% opaque, it means 20% of the weight will remain for the next subsample

			quota -= sub_sample_weight;
			total_weight += sub_sample_weight;
			write_to_samples_map(&vals, color.r, sub_sample_weight);

			traced_bg = traced_bg || (z == camera_far_clip);

			if (quota < 0.0f) {
				// of questionable use
				quota = 0.0f;
				break;
			}
		}
		
		AtRGBA exhaustion_color = AI_RGBA_BLACK;

		if (traced_bg) {
			// died with weight left for a good reason (hit the BG)
			exhaustion_color = AI_RGBA_BLACK;
		} else if (depth_counter >= auto_transparency_depth) {
			// died with weight left for a good reason (ran out of transparency depth)
			exhaustion_color = color;
		} else if (depth_counter == 0) {
			// this happens with volume shaders, for some reason. 
			exhaustion_color = color;
		} else {
			// died with weight left for no good reason
			exhaustion_color = AI_RGBA_BLACK;
		}

		total_weight += quota;

		write_to_samples_map(&vals, exhaustion_color.r, quota);
    }


	///////////////////////////////////////////////
	//
	//		Rank samples and make pixels
	//
	///////////////////////////////////////////////

	if (total_samples == 0) {
		return;
	}

	sw_map_iterator_type vals_iter;
	
	std::vector<std::pair<float, float> > all_vals;
	std::vector<std::pair<float, float> >::iterator all_vals_iter;

	for (vals_iter = vals.begin(); vals_iter != vals.end(); ++vals_iter) {
		all_vals.push_back(*vals_iter);
	}
	
	std::sort(all_vals.begin(), all_vals.end(), compareTail());


	int iter = 0;
	for (all_vals_iter = all_vals.begin(); all_vals_iter != all_vals.end(); ++all_vals_iter) {
		if (iter == rank) {
			out_value->r = all_vals_iter->first;
			out_value->g = (all_vals_iter->second / total_weight);
		}
		else if (iter == rank + 1) {
			out_value->b = all_vals_iter->first;
			out_value->a = (all_vals_iter->second / total_weight);
			return;
		}

		iter++;
	}
}


