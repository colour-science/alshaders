#pragma once
#include <ai.h>

void Kettle_shadows(const float fresnel,
					 const AtFloat rf_scale,
					 const bool rf_col_shadows,
 					 const AtColor rf_color,
					 AtShaderGlobals* sg,
					 const AtNode* node)
{
	//bool rf_useTrans = AiShaderEvalParamBool(p_rf_useTrans);
	// check is object is refractive
	if (rf_scale > 0.0f)
	{
		// check whether we are using absorption
		if (0) // (rf_useTrans && rf_col_shadows)
		{
			/*
			We are doing absorption shadows. shadow rays have multiple intersections
			on a single ray, so we check the facing of the normal on each intersection
			and store it, the color and orientation of the face. Using this info
			we can build a set of rules governing the absorption.

			We presume the geometric rules that objects are either fully inside
			one another or fully outside, overlapping objects will fail.

			We check whether a previous hitpoint exists, if not, then we exit out, and set
			the messages required for the next intersection.

			Otherwise, if we are entering and the previous intersection was also entering,
			we do absorption. If we are exiting, we do absorption.
			*/

			AtPoint s_hit_previous;
			AtColor absorbance = AI_RGB_WHITE;
			AtColor rf_transCol = AI_RGB_WHITE;//AiShaderEvalParamRGB(p_rf_transCol);

			if (AiStateGetMsgPnt("s_hit_previous", &s_hit_previous))
			{
				AiStateGetMsgRGB("s_trans_previous", &rf_transCol);
				AtBoolean do_absorb = false;

				if (AiV3Dot(sg->N, sg->Rd) < 0.0)
				{
					//entering

					AtBoolean old_entering;
					AiStateGetMsgBool("s_hit_entering", &old_entering);
					// if entering, and the previous state was entering,
					// then we need to process absorbtion as we're inside a mesh
					if (old_entering)
					{
						do_absorb = true;
					}
					else
					{
						// ladies and gentlemen, we are floating in space
						rf_transCol	= AI_RGB_WHITE;
					}
					// set entering state for previous
					AiStateSetMsgBool("s_hit_entering", true);
				}
				else
				{
					//exiting, always do absorbtion
					do_absorb = true;
					// set entering state for previous
					AiStateSetMsgBool("s_hit_entering", false);
					//rf_transCol	= rf_transCol;
				}

				if (do_absorb)
				{
					AtFloat distance = AiV3Dist(sg->P, s_hit_previous);
					// do absorbtion

					// invert absorbtion colour (nicety for artists)
					rf_transCol.r = 1.0f - rf_transCol.r;
					rf_transCol.g = 1.0f - rf_transCol.g;
					rf_transCol.b = 1.0f - rf_transCol.b;

					//AiColorScale(absorbance, rf_transCol, -distance * AiShaderEvalParamFlt(p_rf_transExp));

					absorbance.r = exp( absorbance.r );
					absorbance.g = exp( absorbance.g );
					absorbance.b = exp( absorbance.b );

					// invert resulting absorption colour for opacity composite
					absorbance = 1.0f - absorbance;

					if (!AiColorCorrupted(absorbance) && !AiColorIsSmall(absorbance))
						sg->out_opacity = absorbance;
					else
						sg->out_opacity = 0.0f;
				}
				else
				{
					sg->out_opacity = 0.0f;
				}
			}
			else
			{
				sg->out_opacity = 0.0f;
			}

			// store absorbtion current shader colour and mark entering as false (first intersection)
			AiStateSetMsgRGB("s_trans_previous", rf_transCol);
			AiStateSetMsgBool("s_hit_entering", false);
		}
		else
		{
			// this shader does not use absorbtion, but we still need to
			// store values for the hit as the next encountered shader might.
			// We store absorbtion colour as white ("do nothing" color)
			AiStateSetMsgRGB("s_trans_previous", AI_RGB_WHITE);
			sg->out_opacity = 0.0f;
		}

		const AtColor inverse_rf_col = (1.0f - rf_color);
		sg->out_opacity += inverse_rf_col;

		// lerp between not refractive and refractive based on refraction scale. We also blend fresnel.
		// todo: RGB fresnel in shadows? Split into three...?
		AiColorLerp(sg->out_opacity, LERP((1.0f-fresnel), 0.0f, rf_scale), AI_RGB_WHITE, sg->out_opacity);
		//store hitpoint position
		AiStateSetMsgPnt("s_hit_previous", sg->P);
		return;
	}
	// early out
	return;
};
