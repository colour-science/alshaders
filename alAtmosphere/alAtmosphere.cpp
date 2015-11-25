#include "spectrum.h"
#include "Color.h"
#include "alUtil.h"
#include <ai.h>

using sly::HeroSpectrum;
using sly::HeroWavelengths;

AI_SHADER_NODE_EXPORT_METHODS(alAtmosphere);

double rayleigh_scattering_coefficient(double eta, double N,
                                       double lambda_meters)
{
    return (8.0 * pow(AI_PI, 3.0) * pow((pow(eta, 2.0) - 1.0), 2.0)) /
           (3.0 * N * pow(lambda_meters, 4));
}

double planck_blackbody(double temp, double lambda_meters)
{
    return (3.74183e-16 * pow(lambda_meters, -5.0)) /
           (exp(1.4388e-2 / (lambda_meters * temp)) - 1.0);
}

AtRGB sly_to_ai(const sly::RGB& c) { return AiColor(c.r, c.g, c.b); }

struct ShaderData
{
    int AA_samples2;
};

enum alAtmosphereParams
{
    p_units,
    p_mode,
    p_camera_height,
    p_earth_center,
    p_earth_radius,
    p_atmos_radius,
    p_sun_direction,
    p_aerosol_density,
    p_exposure
};

const char* mode_names[] = {"sky", "planet", NULL};

node_parameters
{
    AiParameterFLT("units", 0.01f);
    AiParameterENUM("mode", 0, mode_names);  // 0 = sky, 1 = planet
    AiParameterFLT("camera_height", 2);
    AiParameterVEC("earth_center", 0, -6360000, 0);
    AiParameterFLT("earth_radius", 6360000);
    AiParameterFLT("atmos_radius", 6420000);
    AiParameterVEC("sun_direction", 0, -1, 0);
    AiParameterFLT("aerosol_density", 1);
    AiParameterFLT("exposure", -41.5);
}

node_loader
{
    if (i > 0) return 0;
    node->methods = alAtmosphere;
    node->output_type = AI_TYPE_RGB;
    node->name = "alAtmosphere";
    node->node_type = AI_NODE_SHADER;
    strcpy(node->version, AI_VERSION);
    return true;
}

node_initialize
{
    ShaderData* data = (ShaderData*)new ShaderData();
    AiNodeSetLocalData(node, data);
    sly::SampledSpectrum::init();
}

node_finish
{
    ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);
    delete data;
    AiNodeSetLocalData(node, NULL);
}

node_update
{
    ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);
    AtNode* options = AiUniverseGetOptions();
    data->AA_samples2 = SQR(AiNodeGetInt(options, "AA_samples"));
}

bool ray_sphere(AtVector o, AtVector d, AtVector c, float r, float& t0,
                float& t1)
{
    float r2 = SQR(r);
    AtVector l = c - o;
    float tca = AiV3Dot(l, d);
    float d2 = AiV3Dot(l, l) - SQR(tca);

    if (d2 > r2) return false;

    float thc = sqrtf(r2 - d2);
    t0 = tca - thc;
    t1 = tca + thc;

    if (t0 < 0 && t1 < 0) return false;

    if (t0 > t1) std::swap(t0, t1);

    return true;
}

shader_evaluate
{
    ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);

    // Ci is the color of the surface that we've hit (only available in volume
    // context). For shadows this (I think) represents the transmission so
    // by setting sg->out.RGB we can specify the attenuation along the ray
    AtRGB Ci = sg->Ci;
    sg->out.RGB = Ci;

    sly::ColorSpaceRGB color_space = sly::Primaries::Rec709;
    const sly::SpectralPowerDistribution& spd_white =
        HeroSpectrum::getIlluminant(color_space.wp);

    // Vo is the output volume radiance, in other words in-scattering
    sg->Vo = AI_RGB_BLACK;

    if (sg->Rt & AI_RAY_CAMERA)
    {
        float units = AiShaderEvalParamFlt(p_units);

        int mode = AiShaderEvalParamInt(p_mode);

        float earth_radius = AiShaderEvalParamFlt(p_earth_radius) * units;
        AtVector earth_center;
        if (mode == 0)
        {
            // calculate earth center from camera's height about ground.
            float camera_height = AiShaderEvalParamFlt(p_camera_height) * units;
            earth_center = sg->Ro - AiVector(0, camera_height, 0) -
                           AiVector(0, earth_radius, 0);
        }
        else
        {
            earth_center = AiShaderEvalParamVec(p_earth_center) * units;
        }

        float atmos_radius = AiShaderEvalParamFlt(p_atmos_radius) * units;
        AtVector sun_direction =
            -AiV3Normalize(AiShaderEvalParamVec(p_sun_direction));
        float aerosol_density = AiShaderEvalParamFlt(p_aerosol_density);
        float exposure = AiShaderEvalParamFlt(p_exposure);
        const float H_r = 7994 * units;
        const float H_m = 1200 * units;

#define ATMOS_SPECTRAL
#ifdef ATMOS_SPECTRAL
        HeroSpectrum beta_R;
        HeroSpectrum sun_radiance;
        float u = (float(sg->si) + drand48()) / float(data->AA_samples2);
        HeroWavelengths hw = sly::sample_hero_wavelength(u);
        for (int i = 0; i < 4; ++i)
        {
            beta_R[i] = rayleigh_scattering_coefficient(1.000276, 2.504e25,
                                                        hw[i] * 1e-9);
            sun_radiance[i] = sly::Illuminant::D65.value(hw[i]);
        }

        const float sun_intensity = 1;
#else
        HeroSpectrum beta_R;
        HeroSpectrum sun_radiance(1);
        beta_R[0] = rayleigh_scattering_coefficient(1.000276, 2.504e25, 680e-9);
        beta_R[1] = rayleigh_scattering_coefficient(1.000276, 2.504e25, 550e-9);
        beta_R[2] = rayleigh_scattering_coefficient(1.000276, 2.504e25, 440e-9);
        const float sun_intensity = 20;
#endif
        beta_R /= units;
        const float beta_M = 21e-6 * aerosol_density / units;

        float t0, t1;
        if (ray_sphere(sg->Ro, sg->Rd, earth_center, atmos_radius, t0, t1))
        {
            float mu = AiV3Dot(sun_direction, sg->Rd);
            float phase_R = 3 / (16 * AI_PI) * (1 + mu * mu);
            float g = 0.76;
            float phase_M = 3 / (8 * AI_PI) * ((1 - g * g) * (1 + mu * mu)) /
                            ((2 + g * g) * pow(1 + g * g - 2 * g * mu, 1.5));

            // clip the atmosphere to the ray segment
            t0 = std::max(0.0f, t0);
            t1 = std::min(t1, (float)sg->Rl);

            // clip the atmosphere to the ground
            float tg0, tg1;
            if (ray_sphere(sg->Ro, sg->Rd, earth_center, earth_radius, tg0,
                           tg1))
            {
                if (tg0 < 0)
                {
                    // camera is inside the earth
                    return;
                }

                if (tg1 < t1)
                {
                    // ray is pointing at the ground
                    t1 = tg1;
                }
            }

            const int num_samples = 32;
            const int num_light_samples = 8;

            float ds = (t1 - t0) / float(num_samples);

            float t = t0;
            float optical_depth_r = 0;
            float optical_depth_m = 0;
            HeroSpectrum sum_in_r;
            HeroSpectrum sum_in_m;
            for (int i = 0; i < num_samples; ++i)
            {
                AtVector position = sg->Ro + sg->Rd * (t + drand48() * ds);
                AtVector center_offset = position - earth_center;
                float height = AiV3Length(center_offset) - earth_radius;

                float h_r = exp(-height / H_r) * ds;
                float h_m = exp(-height / H_m) * ds;

                optical_depth_r += h_r;
                optical_depth_m += h_m;

                // first check if we're shadowed
                AtRay shadow_ray;
                AtScrSample scr_sample;
                AiMakeRay(&shadow_ray, AI_RAY_SHADOW, &position, &sun_direction,
                          AI_BIG, sg);
                if (AiTrace(&shadow_ray, &scr_sample))
                {
                    t += ds;
                    continue;
                }

                // check the distance from the current position to the edge of
                // the
                // armosphere
                ray_sphere(position, sun_direction, earth_center, atmos_radius,
                           t0, t1);
                float ds_l = t1 / float(num_light_samples);
                float t_l = 0;
                float optical_depth_light_r = 0;
                float optical_depth_light_m = 0;
                int j = 0;

                for (; j < num_light_samples; ++j)
                {
                    AtVector position_light =
                        position + sun_direction * (t_l + drand48() * ds_l);
                    float height_light =
                        AiV3Length(position_light - earth_center) -
                        earth_radius;
                    if (height_light < 0) break;

                    optical_depth_light_r += exp(-height_light / H_r) * ds_l;
                    optical_depth_light_m += exp(-height_light / H_m) * ds_l;

                    t_l += ds_l;
                }
                if (j == num_light_samples)
                {
                    HeroSpectrum tau_r =
                        beta_R * (optical_depth_r + optical_depth_light_r);
                    float tau_m = beta_M * 1.1f *
                                  (optical_depth_m + optical_depth_light_m);
                    HeroSpectrum tau = tau_r + tau_m;
                    HeroSpectrum attenuation = exp(-tau);
                    sum_in_r += h_r * attenuation;
                    sum_in_m += h_m * attenuation;
                }

                t += ds;
            }

            // compute total in-scattering
            HeroSpectrum inscattering_rayleigh =
                sum_in_r * phase_R * beta_R * sun_intensity * sun_radiance;
            HeroSpectrum inscattering_mie =
                sum_in_m * phase_M * beta_M * sun_intensity * sun_radiance;
            HeroSpectrum inscattering =
                inscattering_rayleigh + inscattering_mie;

            // compute attenuation from full optical depth along camera ray
            HeroSpectrum tau =
                beta_R * optical_depth_r + beta_M * 1.1f * optical_depth_m;
            HeroSpectrum attenuation = exp(-tau);
#ifdef ATMOS_SPECTRAL
            // convert inscattering and attenuation to RGB to apply
            // std::cerr << "conversion: " << hw << std::endl;
            sly::RGB inscattering_rgb = inscattering.toRGB(
                sly::Primaries::Rec709, hw, sly::WhitePoint::ONE);
            sly::RGB attenuation_rgb = inscattering.toRGB(
                sly::Primaries::Rec709, hw, sly::WhitePoint::ONE);

            AtRGB inscattering_atrgb =
                sly_to_ai(inscattering_rgb) * pow(2, exposure);
            AtRGB attenuation_atrgb =
                sly_to_ai(attenuation_rgb) * pow(2, exposure);
            ;
#else
            AtRGB inscattering_atrgb =
                AiColor(inscattering[0], inscattering[1], inscattering[2]);
            AtRGB attenuation_atrgb =
                AiColor(attenuation[0], attenuation[1], attenuation[2]);
#endif
            sg->Vo = inscattering_atrgb;
            sg->out.RGB = Ci * attenuation_atrgb;
        }
    }
}
