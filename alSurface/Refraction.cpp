#include "Refraction.h"
#include "alSurface.h"
#include "alUtil.h"

static void micro_refract_integrate (AtColor &result, AtBoolean &entering, AtFloat &z_rf, AtSampler *sampler, AtShaderGlobals *sg, const AtFloat roughness, const AtFloat rf_IOR, const AtFloat rf_IOR_out, const AtBoolean do_env)
{
        /* microfacet beckman distribution cook torrance model.
        Massivly oiked from http://code.google.com/p/openshadinglanguage/source/browse/trunk/src/liboslexec/bsdf_microfacet.cpp
        */

        AtColor         accum = AI_RGB_BLACK;
        AtVector        omega_out = -sg->Rd;
        AtFloat         z = 0.0f;
        AtFloat         cosNO   = AiV3Dot(omega_out, sg->N);

        AtInt counter = 0;

        AtSamplerIterator *iter;
        iter = AiSamplerIterator(sampler, sg);

        AtDouble        rnd[2];
        AtScrSample     sample;
        AtRay           ray;

        AiMakeRay(&ray, AI_RAY_REFRACTED, &sg->P, NULL, AI_BIG, sg);

        AtFloat n1, n2;
        AtVector omega_in, mf, m;

        AtFloat         cosThetaM, tanThetaM, alpha2;

        while (AiSamplerGetSample(iter, rnd))
        {
                counter ++;
                m = beckmann_sample(roughness, sg->N, rnd, cosThetaM, tanThetaM, alpha2);
                // eq. 33
                AtFloat D = beckmann_distribution(cosThetaM, tanThetaM, alpha2);
                // eq. 24
                AtFloat pm = D * cosThetaM;

                if (AiV3Dot(m, sg->Rd) < 0.0)
                {
                        //entering
                        //AiMsgWarning("entering");
                        n1 = rf_IOR_out;
                        n2 = rf_IOR;
                        entering = true;
                        mf = m;
                }
                else
                {
                        //exiting
                        //AiMsgWarning("exiting");
                        n1 = rf_IOR;
                        n2 = rf_IOR_out;
                        entering = false;
                        AiV3Neg(mf, m);
                }
        // if (mod_refract(rf_IOR, rf_IOR_out, m, omega_out, omega_in, is_inside)) // && !is_inside)
                if (AiRefract(&sg->Rd, &mf, &omega_in, n1, n2))
                {
                        // eval BRDF*cosNI
                        AtFloat cosNI = AiV3Dot(sg->N, omega_in);
                        // eq. 26, 27: now calculate G1(i,m) and G1(o,m)
                        AtFloat G = beckmann_shadowing(roughness, cosNO, cosNI);
                        // eq. 21
                        AtFloat cosHI = AiV3Dot(mf, omega_in);
                        AtFloat cosHO = AiV3Dot(mf, omega_out);
                        AtFloat Ht2 = rf_IOR * cosHI + cosHO;
                        Ht2 *= Ht2;
                        AtVector mfn = mf;
                        AtVector oin = omega_in;
                        float eta = rf_IOR;
                        if (AiV3Dot(mfn, oin) < 0.0f)
                        {
                        	AiV3Neg(oin, oin);
                        	eta = 1.0f/eta;
                        }
                        float chi = AiV3Dot(mfn, oin);
                        float kt =  1.0f - fresnel(chi,eta);
                        //float kt = 1.0f;
                        AtFloat brdf = (fabsf(cosHI * cosHO) * (rf_IOR * rf_IOR) * (G * D)) / (fabsf(cosNO * Ht2)+(AtFloat)AI_EPSILON); // epsilon to avoid divide by zero errors at air like IOR values
                        // eq. 38 and eq. 17
                        AtFloat pdf = pm * (rf_IOR * rf_IOR) * fabsf(cosHI) / (Ht2+(AtFloat)AI_EPSILON); // epsilon to avoid divide by zero errors at air like IOR values

                        // trace ray
                        ray.dir = omega_in;
                        if (do_env)
                                AiTraceBackground(&ray, &sample);
                        else
                                AiTrace(&ray, &sample);
                        accum += sample.color * kt * (brdf/pdf);
                        z += (AtFloat)sample.z;
        }
                else
                {
                        // TIR ray
                        AtFloat cosMO = AiV3Dot(m, omega_out);
                        AiReflect(&sg->Rd, &mf, &omega_in);
                        // remainder or reflective BRDF
                        AtFloat cosNI = AiV3Dot(sg->N, omega_in);
                        // eq. 26, 27: now calculate G1(i,m) and G1(o,m)
                        AtFloat G = beckmann_shadowing(roughness, cosNO, cosNI);
                        // eq. 20: (F*G*D)/(4*in*on)
                        AtFloat brdf = (G * D) * 0.25f / cosNO;
                        AtFloat pdf = pm * 0.25f / cosMO;
                        ray.dir = omega_in;
                        // trace ray
                        ray.dir = omega_in;
                        if (do_env)
                                AiTraceBackground(&ray, &sample);
                        else
                                AiTrace(&ray, &sample);
                        accum += sample.color *(brdf/pdf);
                        z += (AtFloat)sample.z;

                }
        }
        if (counter>0){
                result += accum/(AtFloat)counter;
                z_rf += z/(AtFloat)counter;
        }
}

// ==============================
//
// kettle refraction function
//
// ==============================

void microfacetRefraction(AtShaderGlobals *sg,
                                  const ShaderData *data,
                                  const AtFloat rf_IOR,
                                  AtFloat rf_roughness,
                                  AtColor &refr_result)
{
        AtInt           max_refract                     = data->GI_refraction_depth;
        AtInt           max_gloss                       = data->GI_glossy_depth;
        AtInt           rf_brdf                         = 0; //= (AtInt)AiShaderEvalParamFlt(p_rf_brdf);
        //AtFloat         rf_roughness           			= = AiShaderEvalParamFlt(p_rf_roughness);
        // we exponent roughness here, it's too sensitive.
        //rf_roughness = powf(rf_roughness, 2.0f);
        AtBoolean       do_rf_env                       = false;
        AtColor         rf_exitCol;
        AtColor         rf_exitCol_swatch       		= AI_RGB_BLACK; //= AiShaderEvalParamRGB(p_rf_exitCol);
        AtBoolean       rf_useEnv                       = true; //= AiShaderEvalParamBool(p_rf_useEnv);
        AtFloat         rf_IOR_out                      = 1.0f; // = AiShaderEvalParamFlt(p_rf_IOR_out);
        AtFloat         refr_alpha                      = 1.0f;

        AtVector        R;
        AtRay       ray;
        AtScrSample sample;

        /*
        if (sss_pc_creation)
        {
                rf_roughness = 0.0f;
                if (max_refract > 1)
                        max_refract = 1;
        }
        */

        /*
        Limit the glossy depth if it's greater than the max reflect parameter
        and said parameter is active. No way to do this through node init,
        since glossy rays are reflection OR refraction, therefore we have
        check reflection and refraction
        */

        // maximum refraction depth reached
        if (sg->Rr_refr == max_refract)         // is mirror
                do_rf_env = true;

        // set up exit color
        if (rf_useEnv)
                rf_exitCol = AI_RGB_BLACK;
        else
                rf_exitCol = rf_exitCol_swatch;

        AtFloat z_rf = 0.0f;
        AtBoolean entering;

        // don't process if we're on an exit ray and not sampling enviroment
        if ((do_rf_env==(AtBoolean)true) && (rf_useEnv == (AtBoolean)false))
        {
                // do nothing
        }
        else
        {
                // check whether refraction is glossy, we have the BRDF option in here if we add further BRDFs in the future.
                if (rf_roughness >= AI_EPSILON)
                {
                        // two different sampler loops, to minimise logic inside the sampler
                        if (!do_rf_env)
                        {
                                // full trace
                                switch (rf_brdf)
                                {
                                case 0: // spi
                                        micro_refract_integrate(refr_result, entering, z_rf, data->refraction_sampler, sg, rf_roughness, rf_IOR, rf_IOR_out, false);
                                        break;
                                default:
                                        break;
                                }
                        }
                        // check if this needs an enviroment lookup (if it has a falloff, or if it's an env_lookup)
                        else if (rf_useEnv)
                        {
                                // env only
                                switch (rf_brdf)
                                {
                                case 0: // spi
                                        micro_refract_integrate(rf_exitCol, entering, z_rf, data->refraction_sampler, sg, rf_roughness, rf_IOR, rf_IOR_out, true);
                                        break;
                                default:
                                        break;
                                }
                        }
                }
                else
                {
                        // mirror refraction (no samples)

                        AtVector        refract_dir;
                        AtFloat         n1, n2;

                        if (AiV3Dot(sg->N, sg->Rd) < 0.0)
                        {
                                        //entering
                                        n1 = rf_IOR_out;
                                        n2 = rf_IOR;
                                        entering = true;
                        }
                        else
                        {
                                        //exiting
                                        n1 = rf_IOR;
                                        n2 = rf_IOR_out;
                                        entering = false;
                        }

                        if (AiRefract(&sg->Rd, &sg->Nf, &refract_dir, n1, n2))
                        {
                                R = refract_dir;
                                AiMakeRay(&ray, AI_RAY_REFRACTED, &sg->P, &R, AI_BIG, sg);
                                if (!do_rf_env)
                                {
                                        AiTrace(&ray, &sample);
                                        refr_result = sample.color;
                                        z_rf += (AtFloat)sample.z;
                                }
                                else if (rf_useEnv)
                                {
                                        AiTraceBackground(&ray, &sample);
                                        rf_exitCol = sample.color;
                                }
                        }
                        else    // Refraction ray failed, TIR
                        {
                                AtVector TIR;
                                AiReflectSafe(&sg->Rd, &sg->Nf, &sg->Ng, &TIR);
                                if (rf_roughness > AI_EPSILON)
                                        AiMakeRay(&ray, AI_RAY_GLOSSY, &sg->P, &TIR, AI_BIG, sg);
                                else
                                        AiMakeRay(&ray, AI_RAY_REFRACTED, &sg->P, &TIR, AI_BIG, sg);

                                if (!do_rf_env)
                                {
                                        AiTrace(&ray, &sample);
                                        refr_result = sample.color;
                                        //refr_result = AI_RGB_RED;
                                        z_rf += (AtFloat)sample.z;
                                }
                                else if (rf_useEnv)
                                {
                                        AiTraceBackground(&ray, &sample);
                                        rf_exitCol = sample.color;
                                }
                        }
                }
                // process the result

                // multiply the enviroment color by the exit color
                if (rf_useEnv)
                        rf_exitCol *= rf_exitCol_swatch;

                // process transmission, only do when inside a volume (entering)
                /*
                if (AiShaderEvalParamBool(p_rf_useTrans) && entering && (do_rf_env == false))
                {
                        AtColor rf_transCol     = AiShaderEvalParamRGB(p_rf_transCol);
                        // invert absorbtion colour (nicety for artists)
                        rf_transCol.r = 1.0f - rf_transCol.r;
                        rf_transCol.g = 1.0f - rf_transCol.g;
                        rf_transCol.b = 1.0f - rf_transCol.b;

                        AtColor absorbance;
                        AiColorScale(absorbance, rf_transCol, -z_rf * AiShaderEvalParamFlt(p_rf_transExp));

                        absorbance.r = expf( absorbance.r );
                        absorbance.g = expf( absorbance.g );
                        absorbance.b = expf( absorbance.b );

                        refr_result *= absorbance;
                }
                */
        }

        // set env colour if relavent
        if (do_rf_env)
        {
                refr_alpha = 0.0f;
                refr_result = rf_exitCol;
        }
}

