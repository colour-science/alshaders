/*
 * Simple hair shader, roughly based on Marschner's shading model
 */

#define _STANDALONE_MARSCHNER_   // NOTE: define this to build this shader as standalone

#ifdef _STANDALONE_MARSCHNER_
#include <ai.h>
inline float fast_expf(float v) { return fast_exp(v); }
#else
#include "../util/fast_math.h"
#endif

#include <ai.h>


enum Hair_shader_params
{
   p_opacity,
   p_root_color_multiplier,
   p_tip_color_multiplier,
   p_max_v,
   p_uparam,
   p_vparam,

   p_R_color,
   p_R_intensity,
   p_R_width,
   p_R_shift,

   p_TT_color,
   p_TT_intensity,
   p_TT_width,
   p_TT_shift,
   p_TT_perp_width,

   p_TRT_color,
   p_TRT_intensity,
   p_TRT_width,
   p_TRT_shift,

   p_glints_rel_intensity,
   p_glints_perp_width,
   p_glints_half_angle,
   p_glints_random,

   p_back_scatter_color,
   p_back_scatter_intensity,
   p_back_scatter_width_factor,
   p_back_scatter_shift_offset,

   p_forward_scatter_color,
   p_forward_scatter_intensity,

   p_override_MS_params,
   p_R_color_MS,
   p_R_intensity_MS,
   p_R_width_MS,
   p_R_shift_MS,
   p_TT_color_MS,
   p_TT_intensity_MS,
   p_TT_width_MS,
   p_TT_shift_MS,
   p_TT_perp_width_MS,
   p_TRT_color_MS,
   p_TRT_intensity_MS,
   p_TRT_width_MS,
   p_TRT_shift_MS,
   p_glints_rel_intensity_MS,
   p_glints_perp_width_MS,
   p_glints_half_angle_MS,

   p_num_hairs_fully_opaque,

   p_indirect_color,
   p_indirect_kd,
   p_indirect_n_samples,
   p_indirect_n_mint,
   p_indirect_n_maxt,
   p_indirect_n_spread,
   p_indirect_n_falloff,

   p_aov_hair_f_s_direct,
   p_aov_hair_f_s_direct_R,
   p_aov_hair_f_s_direct_TRT,
   p_aov_hair_f_s_direct_TT,
   p_aov_hair_f_forward_scatter,
   p_aov_hair_f_back_scatter,
   p_aov_hair_indirect_diffuse

};

AI_SHADER_NODE_EXPORT_METHODS(MarschnerMtd);



node_parameters
{
   AiParameterRGB ( "opacity"                    ,  .35f, .35f, .35f );
   AiParameterRGB ( "root_color_multiplier"      ,  1.0f, 1.0f, 1.0f );
   AiParameterRGB ( "tip_color_multiplier"       ,  1.0f, 1.0f, 1.0f );
   AiParameterFlt ( "max_v"                      ,  1.0f  );
   AiParameterStr ( "uparam"                     ,  NULL  );
   AiParameterStr ( "vparam"                     ,  NULL  );

   // R, Primary Reflection
   AiParameterRGB ( "primary_spec_color"         ,  1.0f, 1.0f, 0.7f );  // Color
   AiParameterFlt ( "primary_spec_intensity"     ,  1.0f  );  // Intensity
   AiParameterFlt ( "primary_spec_width"         ,  8.0f  );  // Longitudinal Width of Reflection, 5 to 10 degrees,
   AiParameterFlt ( "primary_spec_shift"         , -8.0f  );  // Longitudinal Offset of Reflection, -10 to -5 degrees

   // TT
   AiParameterRGB ( "transmit_color"             ,  1.0f, 0.6f, 0.28f );  // Color
   AiParameterFlt ( "transmit_intensity"         ,  4.0f  );  // Intensity
   AiParameterFlt ( "transmit_width"             ,  4.0f  );  // Width, should be around R_width / 2.0
   AiParameterFlt ( "transmit_shift"             ,  4.0f  );  // Longitudinal Shift, degrees should be - R_shift / 2
   AiParameterFlt ( "transmit_perp_width"        ,  10.f  );  // Perpendicular Width (across hairs), degrees

   // TRT, Secondary Reflection
   AiParameterRGB ( "secondary_spec_color"       ,  1.0f, 0.62f, 0.2f );  // Color
   AiParameterFlt ( "secondary_spec_intensity"   ,  0.5f  );  // Intensity
   AiParameterFlt ( "secondary_spec_width"       ,  15.f  );  // Longitudinal Width of Secondary Reflection, should be around R_width * 2.0
   AiParameterFlt ( "secondary_spec_shift"       ,  5.0f  );  // Longitudinal Offset of Secondary Reflection, should be - 3. / 2. * R_shift

   // TRT, Glints
   AiParameterFlt ( "glints_rel_intensity"       ,  2.5f  );  // Relative Intensity to secondary reflection, .5 to 5.0
   AiParameterFlt ( "glints_perp_width"          ,  10.f  );  // Perpendicular Width (across hairs), 10 to 25 degrees
   AiParameterFlt ( "glints_half_angle"          ,  30.f  );  // Glint Separation (across hairs), degrees (30 to 45)
   AiParameterFlt ( "glints_random"              ,  1.0f  );  // Glint Randomization value, 0 to 1

   // Backwards scattering, F_back_scatter
   AiParameterRGB ( "back_scatter_color"         ,  1.0f, .58f, .3f );  // Color
   AiParameterFlt ( "back_scatter_intensity"     ,  2.0f  );  // Intensity
   AiParameterFlt ( "back_scatter_width_factor"  ,  1.0f  );  // Width Adjust
   AiParameterFlt ( "back_scatter_shift_offset"  ,  0.0f  );  // Shift Adjust

   // Forwards Scattering, F_s_scatter
   AiParameterRGB ( "forward_scatter_color"      ,  1.0f, 1.0f, 1.0f );  // Color
   AiParameterFlt ( "forward_scatter_intensity"  ,  6.0f  );  // Intensity

   // Override multi-scattering parameters with _MS params if needed
   AiParameterBool( "override_MS_params"         ,  false );
   AiParameterRGB ( "primary_spec_color_MS"      ,  1.0f, 1.0f, 1.0f );
   AiParameterFlt ( "primary_spec_intensity_MS"  ,  1.0f  );
   AiParameterFlt ( "primary_spec_width_MS"      ,  8.0f  );
   AiParameterFlt ( "primary_spec_shift_MS"      , -8.0f  );
   AiParameterRGB ( "transmit_color_MS"          ,  1.0f, 1.0f, 1.0f );
   AiParameterFlt ( "transmit_intensity_MS"      ,  4.0f  );
   AiParameterFlt ( "transmit_width_MS"          ,  4.0f  );
   AiParameterFlt ( "transmit_shift_MS"          ,  4.0f  );
   AiParameterFlt ( "transmit_perp_width_MS"     ,  10.f  );
   AiParameterRGB ( "secondary_spec_color_MS"    ,  1.0f, 1.0f, 1.0f );
   AiParameterFlt ( "secondary_spec_intensity_MS",  0.5f  );
   AiParameterFlt ( "secondary_spec_width_MS"    ,  15.f  );
   AiParameterFlt ( "secondary_spec_shift_MS"    ,  5.0f  );
   AiParameterFlt ( "glints_rel_intensity_MS"    ,  2.5f  );
   AiParameterFlt ( "glints_perp_width_MS"       ,  10.f  );
   AiParameterFlt ( "glints_half_angle_MS"       ,  30.f  );

   // TODO: How many hairs cast an opaque shadow?
   // This would depend on hair transparency
   AiParameterFlt ( "num_hairs_fully_opaque"     ,  100.f );

   // Indirect Diffuse
   AiParameterRGB ( "indirect_color"             ,  1.0f, 1.0f, 1.0f );
   AiParameterFlt ( "indirect_kd"                ,  0.f    );
   AiParameterInt ( "indirect_n_samples"         ,  0      );
   AiParameterFlt ( "indirect_n_mint"            ,  0.f    );
   AiParameterFlt ( "indirect_n_maxt"            ,  1000.f );
   AiParameterFlt ( "indirect_n_spread"          ,  1.f    );
   AiParameterFlt ( "indirect_n_falloff"         ,  1.f    );

   AiParameterStr ( "aov_hair_f_s_direct"        , "hair_f_s_direct"        );
   AiParameterStr ( "aov_hair_f_s_direct_R"      , "hair_f_s_direct_R"      );
   AiParameterStr ( "aov_hair_f_s_direct_TRT"    , "hair_f_s_direct_TRT"    );
   AiParameterStr ( "aov_hair_f_s_direct_TT"     , "hair_f_s_direct_TT"     );
   AiParameterStr ( "aov_hair_f_forward_scatter" , "hair_f_forward_scatter" );
   AiParameterStr ( "aov_hair_f_back_scatter"    , "hair_f_back_scatter"    );
   AiParameterStr ( "aov_hair_indirect_diffuse"  , "hair_indirect_diffuse"  );

   // Synonyms
   AiMetaDataSetStr(mds, "primary_spec_color"         , "synonym", "R_color"          );
   AiMetaDataSetStr(mds, "primary_spec_intensity"     , "synonym", "R_intensity"      );
   AiMetaDataSetStr(mds, "primary_spec_width"         , "synonym", "R_width"          );
   AiMetaDataSetStr(mds, "primary_spec_shift"         , "synonym", "R_shift"          );
   AiMetaDataSetStr(mds, "transmit_color"             , "synonym", "TT_color"         );
   AiMetaDataSetStr(mds, "transmit_intensity"         , "synonym", "TT_intensity"     );
   AiMetaDataSetStr(mds, "transmit_width"             , "synonym", "TT_width"         );
   AiMetaDataSetStr(mds, "transmit_shift"             , "synonym", "TT_shift"         );
   AiMetaDataSetStr(mds, "transmit_perp_width"        , "synonym", "TT_perp_width"    );
   AiMetaDataSetStr(mds, "secondary_spec_color"       , "synonym", "TRT_color"        );
   AiMetaDataSetStr(mds, "secondary_spec_intensity"   , "synonym", "TRT_intensity"    );
   AiMetaDataSetStr(mds, "secondary_spec_with"        , "synonym", "TRT_width"        );
   AiMetaDataSetStr(mds, "secondary_spec_shift"       , "synonym", "TRT_shift"        );
   //
   AiMetaDataSetStr(mds, "primary_spec_color_MS"      , "synonym", "R_color_MS"       );
   AiMetaDataSetStr(mds, "primary_spec_intensity_MS"  , "synonym", "R_intensity_MS"   );
   AiMetaDataSetStr(mds, "primary_spec_with_MS"       , "synonym", "R_width_MS"       );
   AiMetaDataSetStr(mds, "primary_spec_shift_MS"      , "synonym", "R_shift_MS"       );
   AiMetaDataSetStr(mds, "transmit_color_MS"          , "synonym", "TT_color_MS"      );
   AiMetaDataSetStr(mds, "transmit_intensity_MS"      , "synonym", "TT_intensity_MS"  );
   AiMetaDataSetStr(mds, "transmit_width_MS"          , "synonym", "TT_width_MS"      );
   AiMetaDataSetStr(mds, "transmit_shift_MS"          , "synonym", "TT_shift_MS"      );
   AiMetaDataSetStr(mds, "transmit_perp_width_MS"     , "synonym", "TT_perp_width_MS" );
   AiMetaDataSetStr(mds, "secondary_spec_color_MS"    , "synonym", "TRT_color_MS"     );
   AiMetaDataSetStr(mds, "secondary_spec_intensity_MS", "synonym", "TRT_intensity_MS" );
   AiMetaDataSetStr(mds, "secondary_spec_with_MS"     , "synonym", "TRT_width_MS"     );
   AiMetaDataSetStr(mds, "secondary_spec_shift_MS"    , "synonym", "TRT_shift_MS"     );
}

struct ShaderData
{
   int max_diffuse_depth;
   int max_reflect_depth;
   int max_refract_depth;
   int max_glossy_depth;

   // average Properties, ideally these would be exact
   AtColor R_color;
   float   R_intensity;
   float   R_width;
   float   R_width2;
   float   R_shift;
   AtColor TT_color;
   float   TT_intensity;
   float   TT_width;
   float   TT_width2;
   float   TT_shift;
   float   TT_perp_width;
   float   TT_perp_width2;
   AtColor TRT_color;
   float   TRT_intensity;
   float   TRT_width;
   float   TRT_width2;
   float   TRT_shift;
   float   glints_rel_intensity;
   float   glints_perp_width;
   float   glints_perp_width2;
   float   glints_half_angle;

   static const int num_steps = 100;

   // average forward attenuation
   AtColor a_f[num_steps];

   // forwards and backwards scattering factors
   float d_f;
   float d_b;

   // average width for backwards scatter
   AtColor sigma_b2[num_steps];

   // average width for forward scatter
   AtColor beta_f2[num_steps];
   // average width for backwards scatter
   AtColor beta_b2[num_steps];

   // average attenuation for backwards scatter
   AtColor A_b[num_steps];

   // average longitudinal shift
   AtColor delta_b[num_steps];

   // forward scatter, azimuthal component integrals
   // Note it is just a number, because we have defined azimuthal component N*
   // as solely dependent on phi (Marschner's paper takes theta into account too
   // to compute Fresnel attenuation)
   AtColor N_R_G;
   AtColor N_TT_G;
   AtColor N_TRT_G;

   AtSampler *sampler;
};

static inline float safe_acosf(float v)
{
   return acosf((v < -1.f) ? -1.f : ((v > 1.f) ? 1.f : v));
}

static float gaussian(float x, float variance)
{
   return fast_expf(-.5f * x * x / variance);
}

static float normalized_gaussian(float x, float variance)
{
   return fast_expf(-.5f * x * x / variance) / sqrtf(2.f * AI_PI * variance);
}

// Primary highlight, Reflection (R)
static AtColor f_R(float theta, float phi, ShaderData &sd)
{
   return sd.R_color * sd.R_intensity *
      // longitudinal pseudo scattering function
      gaussian(theta - sd.R_shift, sd.R_width2) *
      // azimuthal scattering function for the primary highlight
      cosf(phi / 2.f);
}

// Transmission, Forward scattering (TT)
static AtColor f_TT(float theta, float phi, ShaderData &sd)
{
   return sd.TT_color * sd.TT_intensity *
      // longitudinal pseudo scattering function
      gaussian(theta - sd.TT_shift, sd.TT_width2) *
      // azimuthal pseudo scattering function
      gaussian(AI_PI - phi, sd.TT_perp_width2);
}

// Secondary highlight, Internal reflection (TRT)
static AtColor f_TRT(float theta, float phi, ShaderData &sd)
{
   return sd.TRT_color * sd.TRT_intensity *
      // longitudinal pseudo scattering function
      gaussian(theta - sd.TRT_shift, sd.TRT_width2) *
      // azimuthal scattering function (secondary reflection + glints)
      (cosf(phi / 2.f) + sd.glints_rel_intensity * gaussian(sd.glints_half_angle - phi,
                                                               sd.glints_perp_width2));
}

// Create precomputed tables
static void precompute_tables(AtNode *node, ShaderData *data)
{
   bool override_MS_params = AiNodeGetBool(node, "override_MS_params");
   if (override_MS_params)
   {
      // Use parameters specifically set for multiple scattering
      data->R_color     = AiNodeGetRGB(node, "primary_spec_color_MS");
      data->R_intensity = AiNodeGetFlt(node, "primary_spec_intensity_MS");
      data->R_width     = AiNodeGetFlt(node, "primary_spec_width_MS") * AI_DTOR;
      data->R_shift     = AiNodeGetFlt(node, "primary_spec_shift_MS") * AI_DTOR;
      data->R_width2    = SQR(data->R_width);

      data->TT_color       = AiNodeGetRGB(node, "transmit_color_MS");
      data->TT_intensity   = AiNodeGetFlt(node, "transmit_intensity_MS");
      data->TT_width       = AiNodeGetFlt(node, "transmit_width_MS") * AI_DTOR;
      data->TT_shift       = AiNodeGetFlt(node, "transmit_shift_MS") * AI_DTOR;
      data->TT_perp_width  = AiNodeGetFlt(node, "transmit_perp_width_MS") * AI_DTOR;
      data->TT_width2      = SQR(data->TT_width);
      data->TT_perp_width2 = SQR(data->TT_perp_width);

      data->TRT_color            = AiNodeGetRGB(node, "secondary_spec_color_MS");
      data->TRT_intensity        = AiNodeGetFlt(node, "secondary_spec_intensity_MS");
      data->TRT_width            = AiNodeGetFlt(node, "secondary_spec_width_MS") * AI_DTOR;
      data->TRT_shift            = AiNodeGetFlt(node, "secondary_spec_shift_MS") * AI_DTOR;
      data->glints_rel_intensity = AiNodeGetFlt(node, "glints_rel_intensity_MS");
      data->glints_perp_width    = AiNodeGetFlt(node, "glints_perp_width_MS") * AI_DTOR;
      data->glints_half_angle    = AiNodeGetFlt(node, "glints_half_angle_MS") * AI_DTOR;
      data->TRT_width2           = SQR(data->TRT_width);
      data->glints_perp_width2   = SQR(data->glints_perp_width);
   }
   else
   {
      // Use parameters shared with single scatter
      data->R_color     = AiNodeGetRGB(node, "primary_spec_color");
      data->R_intensity = AiNodeGetFlt(node, "primary_spec_intensity");
      data->R_width     = AiNodeGetFlt(node, "primary_spec_width") * AI_DTOR;
      data->R_shift     = AiNodeGetFlt(node, "primary_spec_shift") * AI_DTOR;
      data->R_width2    = SQR(data->R_width);

      data->TT_color       = AiNodeGetRGB(node, "transmit_color");
      data->TT_intensity   = AiNodeGetFlt(node, "transmit_intensity");
      data->TT_width       = AiNodeGetFlt(node, "transmit_width") * AI_DTOR;
      data->TT_shift       = AiNodeGetFlt(node, "transmit_shift") * AI_DTOR;
      data->TT_perp_width  = AiNodeGetFlt(node, "transmit_perp_width") * AI_DTOR;
      data->TT_width2      = SQR(data->TT_width);
      data->TT_perp_width2 = SQR(data->TT_perp_width);

      data->TRT_color            = AiNodeGetRGB(node, "secondary_spec_color");
      data->TRT_intensity        = AiNodeGetFlt(node, "secondary_spec_intensity");
      data->TRT_width            = AiNodeGetFlt(node, "secondary_spec_width") * AI_DTOR;
      data->TRT_shift            = AiNodeGetFlt(node, "secondary_spec_shift") * AI_DTOR;
      data->glints_rel_intensity = AiNodeGetFlt(node, "glints_rel_intensity");
      data->glints_perp_width    = AiNodeGetFlt(node, "glints_perp_width") * AI_DTOR;
      data->glints_half_angle    = AiNodeGetFlt(node, "glints_half_angle") * AI_DTOR;
      data->TRT_width2           = SQR(data->TRT_width);
      data->glints_perp_width2   = SQR(data->glints_perp_width);
   }

   // Hair Density Factors for back and front scattering
   data->d_b = 0.7f;
   data->d_f = 0.7f;

   int step;
   float theta, phi, theta_d, theta_h;
   float num_steps = 100;
   float phi_step = 2.f * AI_PI / data->num_steps;
   float theta_step = 2.f * AI_PI / data->num_steps;

   // average backward attenuation
   AtColor a_b[data->num_steps];
   // average forward scattering shift
   AtColor alpha_f[data->num_steps];
   // average backward scattering shift
   AtColor alpha_b[data->num_steps];

   // This is the normalization factor to make f_s be < 1 , and prevent it from blowing up
   // It is only used for the multiple scattering computations
   AtColor f_s_norm_color = AI_RGB_BLACK;
   for (theta = -AI_PIOVER2; theta < AI_PIOVER2; theta += theta_step)
   {
      for (phi = 0.f; phi < AI_PI; phi += phi_step)
      {
         f_s_norm_color += f_R(theta, phi, *data) +
                           f_TT(theta, phi, *data) +
                           f_TRT(theta, phi, *data);
      }
   }

   f_s_norm_color *= 2.f * theta_step * phi_step;
   float f_s_norm_factor = 1.f / AiColorMaxRGB(f_s_norm_color);

   // a_f is average forward scattering attenuation
   // a_b is average backward scattering attenuation
   // alpha_f is the weighted average for shifts over the front hemisphere
   // alpha_b is the weighted average for shifts over the back hemisphere
   // beta_f2 is the weighted average for squared widths (variance)
   //                 over the front hemisphere
   // beta_b2 is the weighted average for squared widths (variance)
   //                 over the back hemisphere

   // All formulas from Zinke's and Sadeghi's papers
   for (step = 0; step < num_steps; step += 1)
   {
      theta_d =  AI_PI * step / (num_steps - 1) - AI_PIOVER2;

      data->a_f[step] = AI_RGB_BLACK;
      a_b[step] = AI_RGB_BLACK;

      alpha_f[step] = AI_RGB_BLACK;
      alpha_b[step] = AI_RGB_BLACK;

      data->beta_f2[step] = AI_RGB_BLACK;
      data->beta_b2[step] = AI_RGB_BLACK;

      for (theta = -AI_PIOVER2; theta < AI_PIOVER2; theta += theta_step)
      {
         theta_h = (theta_d + theta) * .5f;

         // forwards
         for (phi = AI_PIOVER2; phi < AI_PI; phi += phi_step)
         {
            data->a_f[step] += f_R(theta_h, phi, *data) +
                               f_TT(theta_h, phi, *data) +
                               f_TRT(theta_h, phi, *data);

            alpha_f[step] += data->R_shift * f_R(theta_h, phi, *data) +
                             data->TT_shift * f_TT(theta_h, phi, *data) +
                             data->TRT_shift * f_TRT(theta_h, phi, *data);

            data->beta_f2[step] += data->R_width2 * f_R(theta_h, phi, *data) +
                                   data->TT_width2 * f_TT(theta_h, phi, *data) +
                                   data->TRT_width2 * f_TRT(theta_h, phi, *data);
         }

         // backwards
         for (phi = 0; phi < AI_PIOVER2; phi += phi_step)
         {
            a_b[step] += f_R(theta_h, phi, *data) +
                         f_TT(theta_h, phi, *data) +
                         f_TRT(theta_h, phi, *data);

            alpha_b[step] += data->R_shift * f_R(theta_h, phi, *data) +
                             data->TT_shift * f_TT(theta_h, phi, *data) +
                             data->TRT_shift * f_TRT(theta_h, phi, *data);

            data->beta_b2[step] += data->R_width2 * f_R(theta_h, phi, *data) +
                                   data->TT_width2 * f_TT(theta_h, phi, *data) +
                                   data->TRT_width2 * f_TRT(theta_h, phi, *data);
         }
      }

      alpha_f[step] /= data->a_f[step];
      alpha_b[step] /= a_b[step];

      data->beta_f2[step] /= data->a_f[step];
      data->beta_b2[step] /= a_b[step];

      data->a_f[step] *= 2.f / AI_PI * theta_step * phi_step;
      data->a_f[step] *= f_s_norm_factor;

      a_b[step] *= 2.f / AI_PI * theta_step * phi_step;
      a_b[step] *= f_s_norm_factor;
   }

   // All this could be optimized

   // See Zinke
   AtColor term1, term2, term3, term4;
   for (step = 0; step < num_steps; step += 1)
   {
      term1 = AI_RGB_WHITE - SQR(data->a_f[step]);
      term2 = term1 * term1 * term1;
      term3 = SQR(data->a_f[step]) * SQR(a_b[step]);

      data->delta_b[step] = alpha_b[step] * (AI_RGB_WHITE - 2 * SQR(a_b[step]))
                            / SQR(term1);

      data->delta_b[step] += alpha_f[step] * (2 * SQR(term1) + 4 * term3)
                             / term2;

      data->A_b[step] = a_b[step] * SQR(data->a_f[step]) / term1;
      data->A_b[step] += a_b[step] * term3 / term2;
   }

   // See Zinke
   for (step = 0; step < num_steps; step += 1)
   {
      for (int i = 0; i < 3; i++)
      {
         term1[i] = sqrtf(2.f * data->beta_f2[step][i] + data->beta_b2[step][i]);
         term2[i] = sqrtf(2.f * data->beta_f2[step][i] + 3.f * data->beta_b2[step][i]);
         term3[i] = sqrtf(data->beta_f2[step][i]);
         term4[i] = sqrtf(data->beta_b2[step][i]);
      }

      data->sigma_b2[step] = AI_RGB_WHITE + 0.7 * SQR(data->a_f[step]);
      data->sigma_b2[step] *= (a_b[step] * term1) + (SQR(a_b[step]) * a_b[step] * term2);
      data->sigma_b2[step] /= a_b[step] + SQR(a_b[step]) * a_b[step] * (2 * term3 + 3 * term4);
      data->sigma_b2[step] *= data->sigma_b2[step];
   }

   data->N_R_G   = AI_RGB_BLACK;
   data->N_TT_G  = AI_RGB_BLACK;
   data->N_TRT_G = AI_RGB_BLACK;
   // Since we are ignoring the Fresnel effects for N, this integral is the same for all thetas
   for (phi = AI_PIOVER2; phi < AI_PI; phi += phi_step)
   {
      data->N_R_G   += cosf(phi * 0.5);
      data->N_TT_G  += gaussian(AI_PI - phi, data->TT_perp_width2);
      data->N_TRT_G += cosf(phi * 0.5) +
      data->glints_rel_intensity * gaussian(data->glints_half_angle - phi,
                                            data->glints_perp_width2);
   }

   data->N_R_G   *= 2.0 * AI_ONEOVERPI * phi_step;
   data->N_TT_G  *= 2.0 * AI_ONEOVERPI * phi_step;
   data->N_TRT_G *= 2.0 * AI_ONEOVERPI * phi_step;
}



node_initialize
{
   ShaderData *data = (ShaderData*) AiMalloc(sizeof(ShaderData));
   AiNodeSetLocalData(node, data);
   data->sampler = NULL;
   precompute_tables(node, data);
}

node_update
{
   ShaderData *data = (ShaderData*) AiNodeGetLocalData(node);
   AtNode *options = AiUniverseGetOptions();

   if (data->sampler)
      AiSamplerDestroy(data->sampler);
   int samples = AiNodeGetInt(node, "indirect_n_samples");
   data->sampler = (samples != 0) ?  AiSampler(samples, 2) : NULL;

   data->max_diffuse_depth = AiNodeGetInt(options, "GI_diffuse_depth");
   data->max_reflect_depth = AiNodeGetInt(options, "GI_reflection_depth");
   data->max_refract_depth = AiNodeGetInt(options, "GI_refraction_depth");
   data->max_glossy_depth = AiNodeGetInt(options, "GI_glossy_depth");

   AiAOVRegister(AiNodeGetStr(node, "aov_hair_f_s_direct"), AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
   AiAOVRegister(AiNodeGetStr(node, "aov_hair_f_s_direct_R"), AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
   AiAOVRegister(AiNodeGetStr(node, "aov_hair_f_s_direct_TRT"), AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
   AiAOVRegister(AiNodeGetStr(node, "aov_hair_f_s_direct_TT"), AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
   AiAOVRegister(AiNodeGetStr(node, "aov_hair_f_forward_scatter"), AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
   AiAOVRegister(AiNodeGetStr(node, "aov_hair_f_back_scatter"), AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
   AiAOVRegister(AiNodeGetStr(node, "aov_hair_indirect_diffuse"), AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);

   precompute_tables(node, data);
}

node_finish
{
   ShaderData *data = (ShaderData*) AiNodeGetLocalData(node);
   AiSamplerDestroy(data->sampler);
   AiFree(data);
}

shader_evaluate
{
   // For environment checking
   ShaderData *data = (ShaderData*) AiNodeGetLocalData(node);

   AtRGB opacity = AiShaderEvalParamRGB(p_opacity);

   // This piece of user-data is automatically set by the curves node when
   // using auto-enlargement (min_pixel_width > 0)
   float geo_opacity;
   if (AiUDataGetFlt("geo_opacity", &geo_opacity))
      opacity *= geo_opacity;

   if (AiShaderGlobalsApplyOpacity(sg, opacity))
      return;

   // early out for shadow rays and totally transparent objects
   if ((sg->Rt & AI_RAY_SHADOW) || AiColorIsZero(sg->out_opacity))
      return;

   AtVector V = -sg->Rd;
   AtVector T = AiV3Normalize(sg->dPdv);

   AtColor Cdiff = AI_RGB_BLACK;
   AtColor Cspec = AI_RGB_BLACK;

   // to support texture mapping across different hairs, change the current (u,v)
   // position according to the specified userdata channels
   float oldU = sg->u;
   float oldV = sg->v;
   // FIXME: we use two separate userdata channels here instead of a single POINT2
   //        channel because of a specific dataset from a specific customer, which
   //        we like to use for performance benchmarking. once we cleanup that scene,
   //        we can switch this shader around
   // NOTE:  sg->u,v are untouched if the userdata channels don't exist (in which
   //        case AiUDataGet* returns false, which we don't bother to check here)
   float rootU = sg->u;
   float rootV = sg->v;
   AiUDataGetFlt(AiShaderEvalParamStr(p_uparam), &rootU);
   AiUDataGetFlt(AiShaderEvalParamStr(p_vparam), &rootV);
   sg->u = rootU;
   sg->v = rootV;

   AtColor root_color  = AiShaderEvalParamRGB(p_root_color_multiplier);
   AtColor tip_color   = AiShaderEvalParamRGB(p_tip_color_multiplier);
   float max_v         = AiShaderEvalParamFlt(p_max_v);

   AtColor R_color     = AiShaderEvalParamRGB(p_R_color);
   float R_intensity   = AiShaderEvalParamFlt(p_R_intensity);
   float R_width       = AI_DTOR * AiShaderEvalParamFlt(p_R_width);
   float R_shift       = AI_DTOR * AiShaderEvalParamFlt(p_R_shift);

   AtColor TT_color    = AiShaderEvalParamRGB(p_TT_color);
   float TT_intensity  = AiShaderEvalParamFlt(p_TT_intensity);
   float TT_width      = AI_DTOR * AiShaderEvalParamFlt(p_TT_width);
   float TT_shift      = AI_DTOR * AiShaderEvalParamFlt(p_TT_shift);
   float TT_perp_width = AI_DTOR *AiShaderEvalParamFlt(p_TT_perp_width);

   AtColor TRT_color   = AiShaderEvalParamRGB(p_TRT_color);
   float TRT_intensity = AiShaderEvalParamFlt(p_TRT_intensity);
   float TRT_width     = AI_DTOR * AiShaderEvalParamFlt(p_TRT_width);
   float TRT_shift     = AI_DTOR * AiShaderEvalParamFlt(p_TRT_shift);

   float glints_rel_intensity = AiShaderEvalParamFlt(p_glints_rel_intensity);
   float glints_perp_width    = AI_DTOR * AiShaderEvalParamFlt(p_glints_perp_width);
   float glints_half_angle    = AI_DTOR * AiShaderEvalParamFlt(p_glints_half_angle);
   float glints_random        = AiShaderEvalParamFlt(p_glints_random);

   AtColor back_scatter_color      = AiShaderEvalParamRGB(p_back_scatter_color);
   float back_scatter_intensity    = AiShaderEvalParamFlt(p_back_scatter_intensity);
   float back_scatter_width_factor = AiShaderEvalParamFlt(p_back_scatter_width_factor);
   float back_scatter_shift_offset = AI_DTOR * AiShaderEvalParamFlt(p_back_scatter_shift_offset);

   AtColor forward_scatter_color   = AiShaderEvalParamRGB(p_forward_scatter_color);
   float forward_scatter_intensity = AiShaderEvalParamFlt(p_forward_scatter_intensity);

   float num_hairs_fully_opaque    = AiShaderEvalParamFlt(p_num_hairs_fully_opaque);

   AtColor indirect_color   = AiShaderEvalParamRGB(p_indirect_color);
   float indirect_kd        = AiShaderEvalParamFlt(p_indirect_kd);
   int   indirect_n_samples = AiShaderEvalParamInt(p_indirect_n_samples);
   float indirect_n_mint    = AiShaderEvalParamFlt(p_indirect_n_mint);
   float indirect_n_maxt    = AiShaderEvalParamFlt(p_indirect_n_maxt);
   float indirect_n_spread  = AiShaderEvalParamFlt(p_indirect_n_spread);
   float indirect_n_falloff = AiShaderEvalParamFlt(p_indirect_n_falloff);

   // restore original (u,v)
   sg->u = oldU;
   sg->v = oldV;

   // mix root and tip colors
   AtColor color_multiplier = AiColorLerp(sg->v/max_v, root_color, tip_color);
   TT_color *= color_multiplier;
   TRT_color *= color_multiplier;
   back_scatter_color *= color_multiplier;
   forward_scatter_color *= color_multiplier;
   indirect_color  *= color_multiplier;

   // precomputed squares
   float TT_width2          = SQR(TT_width);
   float R_width2           = SQR(R_width);
   float TT_perp_width2     = SQR(TT_perp_width);
   float TRT_width2         = SQR(TRT_width);
   float glints_perp_width2 = SQR(glints_perp_width);

   // For testing glint randomization
   AtPoint2 noisePos;
   noisePos.x = rootU;
   noisePos.y = rootV;
   glints_random = LERP(glints_random, 1.0f, AiCellNoise2(noisePos));

   // We are using curves to represent tubes, the normal is meaningless for the
   // shading calculations, we need to consider all directions
   sg->fhemi = false;

   float V_dot_T = AiV3Dot(V, T);

   // View direction projected on normal plane, projected outgoing light direction
   AtVector V_proj = AiV3Normalize(V - V_dot_T * T);

   // Light direction projected on normal plane, projected incoming light direction
   AtVector Ld_proj;

   // Marschner, longitudinal outgoing light direction
   float theta_r = SGN(V_dot_T) * safe_acosf(AiV3Dot(V, V_proj));
   // Marschner, longitudinal incoming light direction
   float theta_i;

   // Marschner, longitudinal half angle between incoming and outgoing light directions: (theta_i + theta_r) / 2
   float theta_h;
   // Marschner, difference angle: (theta_i - theta_r) / 2
   float theta_d;
   float cos_theta_d, cos_theta_d2;
   float cos_theta_i;

   // Marschner, relative azimuth angle: phi_i - phi_r
   float phi;

   // Number of hairs between light and shading point
   float num_hairs_in_front;

   float T_dot_L;

   AtColor hair_f_s_direct = AI_RGB_BLACK;
   AtColor hair_f_s_direct_R = AI_RGB_BLACK;
   AtColor hair_f_s_direct_TT = AI_RGB_BLACK;
   AtColor hair_f_s_direct_TRT = AI_RGB_BLACK;
   AtColor hair_f_back_direct = AI_RGB_BLACK;
   AtColor hair_f_s_scatter = AI_RGB_BLACK;
   AtColor hair_f_back_scatter = AI_RGB_BLACK;

   // direct lighting
   AiLightsPrepare(sg);
   while (AiLightsGetSample(sg))
   {
      if (AiLightGetAffectSpecular(sg->Lp) && sg->Li != AI_RGB_BLACK && sg->we != 0.0f)
      {
         //float light_specular = AiLightGetDiffuse(sg->Lp);
         if (1)
         {
            T_dot_L = AiV3Dot(T, sg->Ld);
            Ld_proj = AiV3Normalize(sg->Ld - T_dot_L * T);

            cos_theta_i = AiV3Dot(sg->Ld, Ld_proj);
            theta_i = SGN(T_dot_L) * safe_acosf(cos_theta_i);
            theta_h = (theta_r + theta_i) * 0.5;
            theta_d = (theta_r - theta_i) * 0.5;
            cos_theta_d = cosf(theta_d);
            cos_theta_d2 = SQR(cos_theta_d);
            phi = safe_acosf(AiV3Dot(V_proj, Ld_proj));

            int index_i = ROUND((theta_i / AI_PI + .5f) * (data->num_steps - 1));

            // We use max for occlusion, because we need to determine the maximum number of hairs
            // Are between the shading point and the lihgt. It would not make sense to do it
            // per channel. This is a very rough approximation
            num_hairs_in_front = num_hairs_fully_opaque * AiColorMaxRGB(sg->Lo);

            // "acceptable approx based on production needs", see Sadeghi
            AtColor sigma_f2 = num_hairs_in_front * data->beta_f2[index_i];

            // Note also that we are ignoring Fresnel,
            // so the azimuthal scattering does not depend on theta_h at all

            AtColor f_s_direct_R  = AI_RGB_BLACK;
            AtColor f_s_direct_TT  = AI_RGB_BLACK;
            AtColor f_s_direct_TRT  = AI_RGB_BLACK;

            // Primary Reflection
            f_s_direct_R = R_intensity * R_color *
               // longitudinal pseudo scattering function
               gaussian(theta_h - R_shift, R_width2) *
               // azimuthal scattering function for the primary highlight
               cosf(phi * 0.5);

            // Transmitted
            f_s_direct_TT = TT_intensity * TT_color * color_multiplier *
               // longitudinal pseudo scattering function
               gaussian(theta_h - TT_shift, TT_width2) *
               // azimuthal pseudo scattering function
               gaussian(AI_PI - phi, TT_perp_width2);

            // Secondary Reflection and glints
            f_s_direct_TRT = TRT_intensity * TRT_color * color_multiplier *
               // longitudinal pseudo scattering function
               gaussian(theta_h - TRT_shift, TRT_width2) *
               // azimuthal scattering function (secondary reflection + glints)
               // glint angle is different for each hair strand and has a randomized value between 30o and 45o
               (cosf(phi * 0.5) + glints_rel_intensity * gaussian(glints_half_angle + AI_PIOVER2 * (glints_random - .5f) - phi,
                                                                  glints_perp_width2));

            // Note that f_s_direct is not normalized, as per the paper (but it could)
            AtColor f_s_direct  = f_s_direct_R + f_s_direct_TT + f_s_direct_TRT;

            AtColor f_s_scatter = AI_RGB_BLACK;

            // Only for _forward_ scattering directions !!
            if (phi >= AI_PIOVER2 && forward_scatter_intensity > 0.f && forward_scatter_color != AI_RGB_BLACK)
            {
               for (int i=0; i<3; i++)
               {
                  // R
                  f_s_scatter[i] += data->R_intensity * data->R_color[i] *
                     normalized_gaussian(theta_h - data->R_shift, data->R_width2 + sigma_f2[i]) *
                     data->N_R_G[i];

                  // TT
                  f_s_scatter[i] += data->TT_intensity * data->TT_color[i] *
                     normalized_gaussian(theta_h - data->TT_shift, data->TT_width2 + sigma_f2[i]) *
                     data->N_TT_G[i];

                  // TRT
                  f_s_scatter[i] += data->TRT_intensity * data->TRT_color[i] *
                     normalized_gaussian(theta_h - data->TRT_shift, data->TRT_width2 + sigma_f2[i]) *
                     data->N_TRT_G[i];

                  // Adjust
                  f_s_scatter[i] *= forward_scatter_intensity * forward_scatter_color[i];
               }
            }

            AtColor f_back_direct = AI_RGB_BLACK;
            AtColor f_back_scatter = AI_RGB_BLACK;

            if (back_scatter_intensity > 0.f && back_scatter_color != AI_RGB_BLACK)
            {
               for (int i = 0; i < 3; i++)
               {
                  f_back_direct[i] = back_scatter_intensity * back_scatter_color[i] *
                     2.f * data->A_b[index_i][i] *
                     normalized_gaussian(theta_h - data->delta_b[index_i][i] - back_scatter_shift_offset,
                                         data->sigma_b2[index_i][i] * back_scatter_width_factor)
                     * AI_ONEOVERPI / cos_theta_d2;

                  f_back_scatter[i] = back_scatter_intensity * back_scatter_color[i] *
                     2.f * data->A_b[index_i][i] *
                     normalized_gaussian(theta_h - data->delta_b[index_i][i] - back_scatter_shift_offset ,
                                         (data->sigma_b2[index_i][i] + sigma_f2[i]) * back_scatter_width_factor)
                     * AI_ONEOVERPI / cos_theta_d2;
               }
            }

            // We cannot follow what Zinke does in the paper for raytracing case
            // as it is too noisy. That is because the shading completely changes when
            // sg->Lo != 0.0f
            // To smooth that harsh transition we use powf over the non occluded fraction of
            // Light hitting the hair, that is : illuminated =  (WHITE - sg->Lo)
            // powf(illum, infinite) corresponds to what the paper says (too noisy)
            // powf(illum, >1) for more accurate but noisy results, with some smoothing
            // powf(illum, 1) is used here,
            // powf(illum, <1) could be used to let direct light penetrate deeper into the hair,
            //                 for thinner than normal hair, or to soften even more
            // Also note, we consider that T_f is already factored into sg->Li
            AtColor direct_fraction = AI_RGB_WHITE - sg->Lo;
            AtColor F_direct = direct_fraction * sg->Li * sg->we * (f_s_direct + data->d_b * f_back_direct);
            AtColor F_scatter = (AI_RGB_WHITE - direct_fraction) * sg->Li * sg->we *
               data->d_f * (f_s_scatter + AI_PI * data->d_b * f_back_scatter);

            Cspec += (F_direct + F_scatter) * cos_theta_i;

            // Write to AOVs
            hair_f_s_direct_R += direct_fraction * sg->Li * sg->we *
                        f_s_direct_R * cos_theta_i;
            hair_f_s_direct_TT += direct_fraction * sg->Li * sg->we *
                        f_s_direct_TT * cos_theta_i;
            hair_f_s_direct_TRT += direct_fraction * sg->Li * sg->we *
                        f_s_direct_TRT * cos_theta_i;
            hair_f_s_direct += direct_fraction * sg->Li * sg->we *
                        f_s_direct * cos_theta_i;
            hair_f_back_direct += direct_fraction * sg->Li * sg->we *
                        data->d_b * f_back_direct * cos_theta_i;
            hair_f_s_scatter += (AI_RGB_WHITE - direct_fraction) * sg->Li * sg->we *
                        data->d_f * f_s_scatter * cos_theta_i;
            hair_f_back_scatter += (AI_RGB_WHITE - direct_fraction) * sg->Li * sg->we *
                        data->d_f * AI_PI * data->d_b * f_back_scatter * cos_theta_i;
         }
      }
   }

   // Indirect diffuse hack
   if (sg->Rr_diff < data->max_diffuse_depth)
   {
      if (indirect_kd > 0)
      {
         if (indirect_n_samples != 0)
         {
            AtVector Nbent;
            AiOcclusion(&V, &V, sg, indirect_n_mint, indirect_n_maxt, indirect_n_spread, indirect_n_falloff, data->sampler, &Nbent);
            Cdiff += indirect_kd * indirect_color * color_multiplier * AiIndirectDiffuse(&Nbent,sg);
         }
         else
         {
            Cdiff += indirect_kd * indirect_color * color_multiplier * AiIndirectDiffuse(&V,sg);
         }
      }
   }

   AiAOVSetRGB(sg, AiShaderEvalParamStr(p_aov_hair_f_s_direct), hair_f_s_direct);
   AiAOVSetRGB(sg, AiShaderEvalParamStr(p_aov_hair_f_s_direct_R), hair_f_s_direct_R);
   AiAOVSetRGB(sg, AiShaderEvalParamStr(p_aov_hair_f_s_direct_TRT), hair_f_s_direct_TRT);
   AiAOVSetRGB(sg, AiShaderEvalParamStr(p_aov_hair_f_s_direct_TT), hair_f_s_direct_TT);
   AiAOVSetRGB(sg, AiShaderEvalParamStr(p_aov_hair_f_forward_scatter), hair_f_s_scatter);
   AiAOVSetRGB(sg, AiShaderEvalParamStr(p_aov_hair_f_back_scatter), hair_f_back_scatter + hair_f_back_direct);
   AiAOVSetRGB(sg, AiShaderEvalParamStr(p_aov_hair_indirect_diffuse), Cdiff);

   sg->out.RGB = Cspec + Cdiff;

}

node_loader
{
    if (i > 0) return FALSE;

    node->methods = MarschnerMtd;
    node->output_type = AI_TYPE_RGB;
    node->name = "aiMarschner";
    node->node_type = AI_NODE_SHADER;
    strcpy(node->version, AI_VERSION);
    return TRUE;
}
