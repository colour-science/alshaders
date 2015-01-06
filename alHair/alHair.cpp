
// Hair shader based on 
// [1] ISHair: Importance Sampling for Hair Scattering by Ou et al. 2012
// http://www.cs.dartmouth.edu/~ouj/site/Research/Entries/2012/6/21_ISHair__Importance_Sampling_for_Hair_Scattering.html
// [2] Dual Scattering Approximation For Fast Multiple Scattering in Hair by Zinke et al. 2008


#include <ai.h>
#include "alUtil.h"
#include "stats.h"
#include "scattering.h"
#include <vector>
#include <algorithm>
#include <map>

//#define DEBUG_LUTS
#ifdef DEBUG_LUTS
#include "exr.h"
#endif

#define NUM_LIGHT_GROUPS 8

AI_SHADER_NODE_EXPORT_METHODS(alHair);

enum alHairParams
{
    p_twist,
    p_hairColor,
    p_specularShift,
    p_specularWidth,
    p_extraSamplesDiffuse,
    p_extraSamplesGlossy,
    p_diffuseStrength,
    p_diffuseColor,
    p_diffuseScatteringMix,
    p_specular1Strength,
    p_specular1Color,
    p_specular2Strength,
    p_specular2Color,
    p_glintStrength,
    p_glintRolloff,
    p_transmissionStrength,
    p_transmissionColor,
    p_transmissionRolloff,
    p_opacity,
    p_dualDepth,
    p_densityFront,
    p_densityBack,
    p_singleSaturation,
    p_multipleSaturation,
    p_specular1WidthScale,
    p_specular2WidthScale,
    p_transmissionWidthScale,
    p_specular1Shift,
    p_specular2Shift,
    p_transmissionShift,
    p_glintTexture,
    p_doMis,
    p_diffuseIndirectStrength,
    p_glossyIndirectStrength,

    p_randomTangent,
    p_randomHue,
    p_randomSaturation,

    p_uparam,
    p_vparam,

    p_aiEnableMatte,
    p_aiMatteColor,
    p_aiMatteColorA,

    p_id1,
    p_id2,
    p_id3,
    p_id4,
    p_id5,
    p_id6,
    p_id7,
    p_id8,

    p_aov_diffuse_color,
    p_aov_direct_diffuse,
    p_aov_indirect_diffuse,
    p_aov_direct_local,
    p_aov_indirect_local,
    p_aov_direct_global,
    p_aov_indirect_global,
    p_aov_direct_specular,
    p_aov_indirect_specular,
    p_aov_direct_specular_2,
    p_aov_indirect_specular_2,
    p_aov_direct_glint,
    p_aov_indirect_glint,
    p_aov_direct_transmission,
    p_aov_indirect_transmission,
    p_aov_depth,
    p_aov_light_group_1,
    p_aov_light_group_2,
    p_aov_light_group_3,
    p_aov_light_group_4,
    p_aov_light_group_5,
    p_aov_light_group_6,
    p_aov_light_group_7,
    p_aov_light_group_8,
    p_aov_id_1,
    p_aov_id_2,
    p_aov_id_3,
    p_aov_id_4,
    p_aov_id_5,
    p_aov_id_6,
    p_aov_id_7,
    p_aov_id_8
};

enum AovIndices
{
    k_diffuse_color=0,
    k_direct_diffuse,
    k_indirect_diffuse,
    k_direct_local,
    k_indirect_local,
    k_direct_global,
    k_indirect_global,
    k_direct_specular,
    k_indirect_specular,
    k_direct_specular_2,
    k_indirect_specular_2,
    k_direct_glint,
    k_indirect_glint,
    k_direct_transmission,
    k_indirect_transmission,
    k_depth,
    k_light_group_1,
    k_light_group_2,
    k_light_group_3,
    k_light_group_4,
    k_light_group_5,
    k_light_group_6,
    k_light_group_7,
    k_light_group_8,
    k_id_1,
    k_id_2,
    k_id_3,
    k_id_4,
    k_id_5,
    k_id_6,
    k_id_7,
    k_id_8
};

#define NUM_ID_AOVS 8
static const char* id_names[NUM_ID_AOVS] =
{
    "id_1",
    "id_2",
    "id_3",
    "id_4",
    "id_5",
    "id_6",
    "id_7",
    "id_8",
};

node_parameters
{
    AiParameterFlt("twist", 20.0f);
    AiParameterRGB("hairColor", 1.0f, 0.829f, 0.488f);
    AiParameterFlt("specularShift", -5.0f);
    AiParameterFlt("specularWidth", 5.0f);
    AiParameterInt("extraSamplesDiffuse", 0);
    AiParameterInt("extraSamplesGlossy", 0);
    AiParameterFlt("diffuseStrength", 1.0f);
    AiParameterRGB("diffuseColor", 1.0f, 1.0f, 1.0f);
    AiParameterFlt("diffuseScatteringMix", 1.0f);
    AiParameterFlt("specular1Strength", 1.0f);
    AiParameterRGB("specular1Color", 1.0f, 1.0f, 1.0f);
    AiParameterFlt("specular2Strength", 1.0f);
    AiParameterRGB("specular2Color",1.0f, 1.0f, 1.0f);
    AiParameterFlt("glintStrength", 2.0f);
    AiParameterFlt("glintRolloff", 5.0f);
    AiParameterFlt("transmissionStrength", 1.0f);
    AiParameterRGB("transmissionColor", 1.0f, 1.0f, 1.0f);
    AiParameterFlt("transmissionRolloff", 10.0f);
    AiParameterRGB("opacity", 1.0f, 1.0f, 1.0f);
    AiParameterInt("dualDepth", 0);
    AiParameterFlt("diffuseForward", 0.7f);
    AiParameterFlt("diffuseBack", 0.7f);
    AiParameterFlt("singleSaturation", 0.1f);
    AiParameterFlt("multipleSaturation", 0.15f);
    AiParameterFlt("specular1WidthScale", 1.0f);
    AiParameterFlt("specular2WidthScale", 1.0f);
    AiParameterFlt("transmissionWidthScale", 1.0f);
    AiParameterFlt("specular1Shift", 0.0f);
    AiParameterFlt("specular2Shift", 0.0f);
    AiParameterFlt("transmissionShift", 0.0f);
    AiParameterFlt("glintTexture", 1.0f);
    AiParameterBool("MIS", true);
    AiParameterFlt("diffuseIndirectStrength", 1.0f);
    AiParameterFlt("glossyIndirectStrength", 1.0f);

    AiParameterFlt("randomTangent", 0.0f);
    AiParameterFlt("randomHue", 0.0f);
    AiParameterFlt("randomSaturation", 0.0f);

    AiParameterStr("uparam", "uparamcoord");
    AiParameterStr("vparam", "vparamcoord");

    AiParameterBOOL("aiEnableMatte", false);
    AiParameterRGB("aiMatteColor", 0.0f, 0.0f, 0.0f);
    AiParameterFlt("aiMatteColorA", 0.0f);

    AiParameterRGB("id1", 0.0f, 0.0f, 0.0f);
    AiParameterRGB("id2", 0.0f, 0.0f, 0.0f);
    AiParameterRGB("id3", 0.0f, 0.0f, 0.0f);
    AiParameterRGB("id4", 0.0f, 0.0f, 0.0f);
    AiParameterRGB("id5", 0.0f, 0.0f, 0.0f);
    AiParameterRGB("id6", 0.0f, 0.0f, 0.0f);
    AiParameterRGB("id7", 0.0f, 0.0f, 0.0f);
    AiParameterRGB("id8", 0.0f, 0.0f, 0.0f);

    AiParameterStr("aov_diffuse_color", "diffuse_color");
    AiParameterStr("aov_direct_diffuse", "direct_diffuse");
    AiParameterStr("aov_indirect_diffuse", "indirect_diffuse");
    AiParameterStr("aov_direct_local", "direct_local");
    AiParameterStr("aov_indirect_local", "indirect_local");
    AiParameterStr("aov_direct_global", "direct_global");
    AiParameterStr("aov_indirect_global", "indirect_global");
    AiParameterStr("aov_direct_specular", "direct_specular");
    AiParameterStr("aov_indirect_specular", "indirect_specular");
    AiParameterStr("aov_direct_specular_2", "direct_specular_2");
    AiParameterStr("aov_indirect_specular_2", "indirect_specular_2");
    AiParameterStr("aov_direct_glint", "direct_glint");
    AiParameterStr("aov_indirect_glint", "indirect_glint");
    AiParameterStr("aov_direct_transmission", "direct_transmission");
    AiParameterStr("aov_indirect_transmission", "indirect_transmission");
    AiParameterStr("aov_depth", "depth");
    AiParameterStr("aov_light_group_1", "light_group_1");
    AiParameterStr("aov_light_group_2", "light_group_2");
    AiParameterStr("aov_light_group_3", "light_group_3");
    AiParameterStr("aov_light_group_4", "light_group_4");
    AiParameterStr("aov_light_group_5", "light_group_5");
    AiParameterStr("aov_light_group_6", "light_group_6");
    AiParameterStr("aov_light_group_7", "light_group_7");
    AiParameterStr("aov_light_group_8", "light_group_8");
    AiParameterStr("aov_id_1", "id_1");
    AiParameterStr("aov_id_2", "id_2");
    AiParameterStr("aov_id_3", "id_3");
    AiParameterStr("aov_id_4", "id_4");
    AiParameterStr("aov_id_5", "id_5");
    AiParameterStr("aov_id_6", "id_6");
    AiParameterStr("aov_id_7", "id_7");
    AiParameterStr("aov_id_8", "id_8");
}

#define B_WIDTH_SCALE 2.0f
#define F_WIDTH_SCALE 4.0f

struct HairBsdf
{
    struct ShaderData
    {
        ShaderData()
        : sampler_glossy(NULL), sampler_diffuse(NULL), ds(NULL)
        {}

        ~ShaderData()
        {
            AiSamplerDestroy(sampler_glossy);
            AiSamplerDestroy(sampler_diffuse);
            delete ds;
        }

        void update(AtParamValue* params)
        {
            AtUInt32 t0 = AiMsgUtilGetElapsedTime();

            AtNode *options   = AiUniverseGetOptions();
            int glossy_samples = std::max(0, AiNodeGetInt(options, "GI_glossy_samples") + params[p_extraSamplesGlossy].INT);
            int diffuse_samples = std::max(0, AiNodeGetInt(options, "GI_diffuse_samples") + params[p_extraSamplesDiffuse].INT);

            AiSamplerDestroy(sampler_glossy);
            sampler_glossy = AiSampler(glossy_samples, 2);
            AiSamplerDestroy(sampler_diffuse);
            sampler_diffuse = AiSampler(diffuse_samples, 2);

            alpha_R = -params[p_specularShift].FLT * AI_DTOR;
            beta_R = params[p_specularWidth].FLT * AI_DTOR;
            
            beta_TT = beta_R * 0.5f;
            
            alpha_TT = -alpha_R * 0.5f;

            beta_TRT = beta_R * 2.0f;
            alpha_TRT = -alpha_R * 1.5f;

            beta_R2 = beta_R*beta_R;
            beta_TT2 = beta_TT*beta_TT;
            beta_TRT2 = beta_TRT*beta_TRT;

            gamma_TT = params[p_transmissionRolloff].FLT * AI_DTOR;
            gamma_g = params[p_glintRolloff].FLT * AI_DTOR;

            glintStrength = params[p_glintStrength].FLT;

            dual_depth = params[p_dualDepth].INT;

            sampleLobesIndividually = false;

            doMis = params[p_doMis].BOOL;

            uparam = params[p_uparam].STR;
            vparam = params[p_vparam].STR;

            aovs.clear(); 
            aovs.push_back(params[p_aov_diffuse_color].STR); 
            aovs.push_back(params[p_aov_direct_diffuse].STR); 
            aovs.push_back(params[p_aov_indirect_diffuse].STR); 
            aovs.push_back(params[p_aov_direct_local].STR); 
            aovs.push_back(params[p_aov_indirect_local].STR); 
            aovs.push_back(params[p_aov_direct_global].STR); 
            aovs.push_back(params[p_aov_indirect_global].STR); 
            aovs.push_back(params[p_aov_direct_specular].STR); 
            aovs.push_back(params[p_aov_indirect_specular].STR); 
            aovs.push_back(params[p_aov_direct_specular_2].STR); 
            aovs.push_back(params[p_aov_indirect_specular_2].STR); 
            aovs.push_back(params[p_aov_direct_glint].STR); 
            aovs.push_back(params[p_aov_indirect_glint].STR); 
            aovs.push_back(params[p_aov_direct_transmission].STR); 
            aovs.push_back(params[p_aov_indirect_transmission].STR); 
            aovs.push_back(params[p_aov_depth].STR); 
            aovs.push_back(params[p_aov_light_group_1].STR); 
            aovs.push_back(params[p_aov_light_group_2].STR); 
            aovs.push_back(params[p_aov_light_group_3].STR); 
            aovs.push_back(params[p_aov_light_group_4].STR); 
            aovs.push_back(params[p_aov_light_group_5].STR); 
            aovs.push_back(params[p_aov_light_group_6].STR); 
            aovs.push_back(params[p_aov_light_group_7].STR); 
            aovs.push_back(params[p_aov_light_group_8].STR); 
            aovs.push_back(params[p_aov_id_1].STR); 
            aovs.push_back(params[p_aov_id_2].STR); 
            aovs.push_back(params[p_aov_id_3].STR); 
            aovs.push_back(params[p_aov_id_4].STR); 
            aovs.push_back(params[p_aov_id_5].STR); 
            aovs.push_back(params[p_aov_id_6].STR); 
            aovs.push_back(params[p_aov_id_7].STR); 
            aovs.push_back(params[p_aov_id_8].STR); 
            for (size_t i=0; i < aovs.size(); ++i)
                AiAOVRegister(aovs[i].c_str(), AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);

            delete ds;
            ds = new DualScattering;
        }

        AtSampler* sampler_glossy;
        AtSampler* sampler_diffuse;

        int dual_depth;

        float beta_R;       //< R width
        float beta_R2;
        float alpha_R;      //< R shift
        float beta_TT;      //< TT width
        float beta_TT2;
        float alpha_TT;     //< TT shift
        float beta_TRT;     //< TRT width
        float beta_TRT2;
        float alpha_TRT;    //< TRT shift
        float gamma_TT;     //< TT rolloff
        float gamma_g;      //< g rolloff
        float phi_g;        //< g separation

        float glintStrength;

        bool sampleLobesIndividually;

        DualScattering* ds;

        bool doMis;

        std::vector<std::string> aovs;
        std::map<AtNode*, int> lightGroups;

        std::string uparam, vparam;
    };

    HairBsdf(AtNode* n, AtShaderGlobals* sg, ShaderData* d) :
    node(n), data(d), numBlendHairs(2), density_front(0.7f), density_back(0.7f), _sg(sg)
    {
        depth = sg->Rr;

        result_R_direct = AI_RGB_BLACK;
        result_R_indirect = AI_RGB_BLACK;
        result_TT_direct = AI_RGB_BLACK;
        result_TT_indirect = AI_RGB_BLACK;
        result_TRT_direct = AI_RGB_BLACK;
        result_TRT_indirect = AI_RGB_BLACK;
        result_TRTg_direct = AI_RGB_BLACK;
        result_TRTg_indirect = AI_RGB_BLACK;
        result_Pg_direct = AI_RGB_BLACK;
        result_Pl_direct = AI_RGB_BLACK;
        result_Pg_indirect = AI_RGB_BLACK;
        result_Pl_indirect = AI_RGB_BLACK;
        result_id1 = AI_RGB_BLACK;
        result_id2 = AI_RGB_BLACK;
        result_id3 = AI_RGB_BLACK;
        result_id4 = AI_RGB_BLACK;
        result_id5 = AI_RGB_BLACK;
        result_id6 = AI_RGB_BLACK;
        result_id7 = AI_RGB_BLACK;
        result_id8 = AI_RGB_BLACK;

        
        memset(lightGroupsDirect, 0, sizeof(AtRGB)*NUM_LIGHT_GROUPS);
    }

    /// Parameter evaluation. This should be called after opacity() and before anything else.
    inline void evaluateParameters(AtShaderGlobals* sg, ShaderData* data)
    {
        // save v coord 
        float v = sg->v;

        // set the hair u, v coords so other shaders (e.g. ramp) can use them
        AtPoint2 hair_uv = AiPoint2(sg->u, sg->v);
        AiStateSetMsgPnt2("maya_ramp_uv_override", hair_uv);

        // replace uvs
        float s, t;
        if (AiUDataGetFlt(data->uparam.c_str(), &s) && AiUDataGetFlt(data->vparam.c_str(), &t))
        {
            sg->u = s;
            sg->v = t;
        }

        // Get a random value per curve
        AtUInt32 curve_idi = 0;
        float curve_id = 0.0f;
        cn = 1.0f;
        AtVector cv = aivec(0.0f);
        float randomTangent = AiShaderEvalParamFlt(p_randomTangent) * 0.3f;
        if (AiUDataGetUInt("curve_id", &curve_idi))
        {
            AtPoint2 p; p.x = float(curve_idi); p.y = 0.0f;
            cn = AiCellNoise2(p);

            AtPoint p2 = aivec(float(curve_idi)+17.0f, 0.0f, 0.0f);
            cv = (AiVCellNoise3(p2)*2.0f - aivec(1.0f));
        }
        else if (AiUDataGetFlt("curve_id", &curve_id))
        {
            AtPoint2 p; p.x = curve_id; p.y = 0.0f;
            cn = AiCellNoise2(p);

            AtPoint p2 = aivec(curve_id+17.0f, 0.0f, 0.0f);
            cv = (AiVCellNoise3(p2)*2.0f - aivec(1.0f));

        }
        else
        {
            cn = 0.5f;
        }

        // Get a local coordinate frame based on the hair fibre direction
        U = AiV3Normalize(sg->dPdv + cv*randomTangent);
        V = AiV3Cross(U, sg->N);
        W = AiV3Cross(V, U);

        // Get the spherical angles of the exitant direction relative to the hair fibre
        wo = -sg->Rd;
        theta_r = AI_PIOVER2 - sphericalTheta(wo, U);
        phi_r = sphericalPhi(wo, V, W);

        sp.beta_R = data->beta_R;
        sp.alpha_R = data->alpha_R;
        sp.alpha_R_offset = -AiShaderEvalParamFlt(p_specular1Shift);

        sp.beta_TT = data->beta_TT;
        sp.alpha_TT = data->alpha_TT;
        sp.alpha_TT_offset = -AiShaderEvalParamFlt(p_transmissionShift);

        sp.beta_TRT = data->beta_TRT;
        sp.alpha_TRT = data->alpha_TRT;
        sp.alpha_TRT_offset = -AiShaderEvalParamFlt(p_specular2Shift);

        sp.beta_R2 = data->beta_R2;
        sp.beta_TT2 = data->beta_TT2;
        sp.beta_TRT2 = data->beta_TRT2;

        sp.gamma_TT = data->gamma_TT;
        sp.gamma_g = data->gamma_g;
        float twist = AiShaderEvalParamFlt(p_twist);
        float glintTexture = AiShaderEvalParamFlt(p_glintTexture);
        sp.phi_g = lerp((35.0f-30*(1.0f-glintTexture))*AI_DTOR, (35.0f+50.0f*glintTexture)*AI_DTOR, fabsf(sinf(AI_PI*(cn + (cn+v)*(twist*cn)))));

        sp.ior = 1.55;

        AB(theta_r, sp.alpha_R, sp.beta_R, A_R, B_R);
        AB(theta_r, sp.alpha_TT, sp.beta_TT, A_TT, B_TT);
        AB(theta_r, sp.alpha_TRT, sp.beta_TRT, A_TRT, B_TRT);
        AB(theta_r, sp.alpha_R, sp.beta_R*B_WIDTH_SCALE, A_b, B_b);
        AB(theta_r, sp.alpha_TT, sp.beta_TT*F_WIDTH_SCALE, A_f, B_f);

        float singleSaturation = AiShaderEvalParamFlt(p_singleSaturation);
        float multipleSaturation = AiShaderEvalParamFlt(p_multipleSaturation);
        hairColor = clamp(AiShaderEvalParamRGB(p_hairColor), AI_RGB_BLACK, AI_RGB_WHITE);

        float randomHue = AiShaderEvalParamFlt(p_randomHue) * 0.1f;
        float randomSaturation = AiShaderEvalParamFlt(p_randomSaturation);

        if (randomHue != 0.0f || randomSaturation != 0.0f)
        {
            hairColor = rgb2hsv(hairColor);
            hairColor.r += cv.x * randomHue;
            hairColor.r = clamp(-1.0f, 1.0f, hairColor.r);
            hairColor.g += cv.y * randomSaturation;
            hairColor.g = clamp(0.0f, 1.0f, hairColor.g);
            hairColor = hsv2rgb(hairColor);
        }

        float hcmax = maxh(hairColor);
        AtRGB scol=hairColor, mcol=hairColor;
        scol = pow(scol, singleSaturation);
        mcol = pow(mcol, multipleSaturation);

        sp.absorption = clamp(AI_RGB_WHITE - scol, rgb(0.01f), rgb(0.99f));
        sp.dabsorption = clamp(AI_RGB_WHITE - mcol, rgb(0.01f), rgb(0.99f));

        diffuseColor = AiShaderEvalParamRGB(p_diffuseColor) * AiShaderEvalParamFlt(p_diffuseStrength);
        diffuseScatteringMix = AiShaderEvalParamFlt(p_diffuseScatteringMix);
        specular1Color = AiShaderEvalParamRGB(p_specular1Color) * AiShaderEvalParamFlt(p_specular1Strength);
        specular2Color = AiShaderEvalParamRGB(p_specular2Color) * AiShaderEvalParamFlt(p_specular2Strength);
        transmissionColor = AiShaderEvalParamRGB(p_transmissionColor) * AiShaderEvalParamFlt(p_transmissionStrength);
        glintStrength = AiShaderEvalParamFlt(p_glintStrength);
        

        specular1WidthScale = AiShaderEvalParamFlt(p_specular1WidthScale);
        specular2WidthScale = AiShaderEvalParamFlt(p_specular2WidthScale);
        transmissionWidthScale = AiShaderEvalParamFlt(p_transmissionWidthScale);

        diffuseIndirectStrength = AiShaderEvalParamFlt(p_diffuseIndirectStrength);
        glossyIndirectStrength = AiShaderEvalParamFlt(p_glossyIndirectStrength);

        density_front = AiShaderEvalParamFlt(p_densityFront);
        density_back = AiShaderEvalParamFlt(p_densityBack);

        if (density_back*maxh(diffuseColor) > IMPORTANCE_EPS && density_front*maxh(diffuseColor) > IMPORTANCE_EPS)
            dualImportant = true;
        else
            dualImportant = false;

        do_diffuse = do_glossy = true;
        if (sg->Rr_diff > 0)
        {
            do_glossy = false;
        }

        hairColor *= diffuseColor;

        doLightGroups = false;
        for (int i=0; i < NUM_LIGHT_GROUPS; ++i)
        {
            if (AiAOVEnabled(data->aovs[k_light_group_1+i].c_str(), AI_TYPE_RGB))
            {
                doLightGroups = true;
                break;
            }
        }
    }

    inline void AB(float theta_r, float alpha, float beta, float& A, float& B)
    {
        A = atanf((PIOVER4 + theta_r*0.5f - alpha) / beta);
        B = atanf((-PIOVER4 + theta_r*0.5f - alpha) / beta);
    }

    inline float sampleLong(double u, float theta_r, float alpha, float beta, float A, float B)
    {
        float t = 2.0f*beta * tanf(u*(A-B) + B) + 2.0f*alpha - theta_r;
        return t;//clamp( -0.4999f * AI_PI, 0.4999f * AI_PI, t); //< TODO: not handling the t == +/- AI_PIOVER2 case doesn't seem to hurt?
    }

    inline AtVector sample_R(float u1, float u2)
    {
        float theta_i = sampleLong(u1, theta_r, sp.alpha_R - sp.alpha_R_offset, sp.beta_R*specular1WidthScale, A_R, B_R);
        float phi = 2.0f * asinf(clamp(2.0f*u2 - 1.0f, -1.0f, 1.0f));
        float phi_i = phi_r - phi;
        AtVector wi;
        sphericalDirection(theta_i, phi_i, V, W, U, wi);
        return wi;
    }

    inline AtRGB bsdf_R(const SctGeo& geo)
    {
        return rgb(bsdfR(sp.beta_R2*SQR(specular1WidthScale), sp.alpha_R - sp.alpha_R_offset, geo.theta_h, geo.cosphi2)) * geo.cos_theta_i * geo.inv_cos_theta_d2 * AI_ONEOVERPI;
    }

    inline float pdf_R(const SctGeo& geo)
    {
        float t = geo.theta_h-sp.alpha_R - sp.alpha_R_offset;
        float pdf_theta = (1.0f / (2.0f*geo.cos_theta_i*(A_R-B_R))) * (sp.beta_R*specular1WidthScale / (t*t + sp.beta_R2*SQR(specular1WidthScale)));
        float pdf_phi = geo.cosphi2*0.25f;
        return pdf_theta * pdf_phi;
    }

    inline AtVector sample_TRT(float u1, float u2)
    {
        float theta_i = sampleLong(u1, theta_r, sp.alpha_TRT - sp.alpha_TRT_offset, sp.beta_TRT*specular2WidthScale, A_R, B_R);
        float phi = 2.0f * asinf(clamp(2.0f*u2 - 1.0f, -1.0f, 1.0f));
        float phi_i = phi_r - phi;
        AtVector wi;
        sphericalDirection(theta_i, phi_i, V, W, U, wi);
        return wi;
    }

    inline AtRGB bsdf_TRT(const SctGeo& geo)
    {
        return rgb(bsdfR(sp.beta_TRT2*SQR(specular2WidthScale), sp.alpha_TRT - sp.alpha_TRT_offset, geo.theta_h, geo.cosphi2)) * geo.cos_theta_i * geo.inv_cos_theta_d2 * AI_ONEOVERPI;
    }

    inline float pdf_TRT(const SctGeo& geo)
    {
        float t = geo.theta_h-sp.alpha_TRT - sp.alpha_TRT_offset;
        float pdf_theta = (1.0f / (2.0f*geo.cos_theta_i*(A_R-B_R))) * (sp.beta_TRT*specular2WidthScale / (t*t + sp.beta_TRT2*SQR(specular2WidthScale)));
        float pdf_phi = geo.cosphi2*0.25f;
        return pdf_theta * pdf_phi;
    }

    inline AtVector sample_TRTg(float u1, float u2)
    {
        float theta_i = sampleLong(u1, theta_r, sp.alpha_TRT - sp.alpha_TRT_offset, sp.beta_TRT*specular2WidthScale, A_TRT, B_TRT);
        float sign;
        if (u2 < 0.5f)
        {
            sign = 1.0f;
            u2 = 2.0f * u2;
        }
        else
        {
            sign = -1.0f;
            u2 = 2.0f * (1.0f-u2);
        }
        Cg = atanf((AI_PIOVER2 - sp.phi_g)/sp.gamma_g);
        Dg = atanf(-sp.phi_g/sp.gamma_g);
        float phi = sp.gamma_g * tanf(u2*(Cg-Dg)+Dg) + sp.phi_g;
        phi *= sign;
        float phi_i = phi_r + phi;

        AtVector wi;
        sphericalDirection(theta_i, phi_i, V, W, U, wi);
        return wi;
    }

    inline AtRGB bsdf_TRTg(const SctGeo& geo)
    {
        return rgb(bsdfg(sp.beta_TRT2*SQR(specular2WidthScale), sp.alpha_TRT - sp.alpha_TRT_offset, geo.theta_h, sp.gamma_g, geo.phi, sp.phi_g)) * geo.cos_theta_i * geo.inv_cos_theta_d2 * AI_ONEOVERPI;
    }

    inline float pdf_TRTg(const SctGeo& geo)
    {   
        float t = geo.theta_h-sp.alpha_TRT - sp.alpha_TRT_offset;
        float pdf_theta = (1.0f / (2.0f*geo.cos_theta_i*(A_TRT-B_TRT))) * (sp.beta_TRT*specular2WidthScale / (t*t + sp.beta_TRT2*SQR(specular2WidthScale)));
        float p = fabsf(geo.phi) - sp.phi_g;
        float Cg = atanf((AI_PIOVER2 - sp.phi_g)/sp.gamma_g);
        float Dg = atanf(-sp.phi_g/sp.gamma_g);
        float pdf_phi = (1.0f / (2.0f * (Cg-Dg))) * (sp.gamma_g / (p*p + sp.gamma_g*sp.gamma_g));
        return pdf_theta * pdf_phi;
    }

    inline AtVector sample_TT(float u1, float u2)
    {
        float theta_i = sampleLong(u1, theta_r, sp.alpha_TT - sp.alpha_TT_offset, sp.beta_TT*transmissionWidthScale, A_TT, B_TT);
        C_TT = 2.0f * atanf(AI_PI/sp.gamma_TT);
        float phi = sp.gamma_TT * tanf(C_TT * (u2-0.5f)) + AI_PI;
        float phi_i = phi_r - phi;
        AtVector wi;
        sphericalDirection(theta_i, phi_i, V, W, U, wi);
        return wi;
    }

    inline AtRGB bsdf_TT(const SctGeo& geo)
    {
        return rgb(bsdfTT(sp.beta_TT2*SQR(transmissionWidthScale), sp.alpha_TT - sp.alpha_TT_offset, geo.theta_h, sp.gamma_TT, geo.phi_d)) * geo.cos_theta_i * geo.inv_cos_theta_d2 * AI_ONEOVERPI;
    }

    inline float pdf_TT(const SctGeo& geo)
    {
        float t = geo.theta_h-sp.alpha_TT - sp.alpha_TT_offset;
        float pdf_theta = (1.0f / (2.0f*geo.cos_theta_i*(A_TT-B_TT))) * (sp.beta_TT*transmissionWidthScale / (t*t + sp.beta_TT2*SQR(transmissionWidthScale)));
        float p = geo.phi-AI_PI;
        float C_TT = 2.0f * atanf(AI_PI/sp.gamma_TT);
        float pdf_phi = (1.0f / C_TT) * (sp.gamma_TT / (p*p + sp.gamma_TT*sp.gamma_TT));
        
        return pdf_theta * pdf_phi;
    }

    inline AtVector sample_Sb(float u1, float u2)
    {
        float theta_i = sampleLong(u1, theta_r, sp.alpha_R, sp.beta_R*B_WIDTH_SCALE, A_b, B_b);
        float phi = 2.0f * asinf(clamp(2.0f*u2 - 1.0f, -1.0f, 1.0f));
        float phi_i = phi_r - phi;
        AtVector wi;
        sphericalDirection(theta_i, phi_i, V, W, U, wi);
        return wi;
    }

    inline float pdf_Sb(const SctGeo& geo)
    {
        float t = geo.theta_h-sp.alpha_R*B_WIDTH_SCALE;
        float pdf_theta = (1.0f / (2.0f*geo.cos_theta_i*(A_b-B_b))) * (sp.beta_R*B_WIDTH_SCALE / (t*t + sp.beta_R2*B_WIDTH_SCALE*B_WIDTH_SCALE));
        float pdf_phi = geo.cosphi2*0.25f;
        return pdf_theta * pdf_phi;
    }

    inline AtVector sample_Sf(float u1, float u2)
    {
        float theta_i = sampleLong(u1, theta_r, sp.alpha_TT, sp.beta_TT*F_WIDTH_SCALE, A_f, B_f);
        C_TT = 2.0f * atanf(AI_PI/(sp.gamma_TT*F_WIDTH_SCALE));
        float phi = sp.gamma_TT*F_WIDTH_SCALE * tanf(C_TT * (u2-0.5f)) + AI_PI;
        float phi_i = phi_r - phi;
        AtVector wi;
        sphericalDirection(theta_i, phi_i, V, W, U, wi);
        return wi;
    }

    inline float pdf_Sf(const SctGeo& geo)
    {
        float t = geo.theta_h-sp.alpha_TT;
        float pdf_theta = (1.0f / (2.0f*geo.cos_theta_i*(A_f-B_f))) * (sp.beta_TT*F_WIDTH_SCALE / (t*t + sp.beta_TT2*F_WIDTH_SCALE*F_WIDTH_SCALE));
        float p = geo.phi-AI_PI;
        float C_TT = 2.0f * atanf(AI_PI/(sp.gamma_TT*F_WIDTH_SCALE));
        float pdf_phi = (1.0f / C_TT) * (sp.gamma_TT*F_WIDTH_SCALE / (p*p + sp.gamma_TT*sp.gamma_TT*F_WIDTH_SCALE*F_WIDTH_SCALE));
        
        return pdf_theta * pdf_phi;
    }

    static AtVector Hair_Sample_R(const void* brdf_data, float u1, float u2)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        return hb->sample_R(u1, u2);
    }

    static AtVector Hair_Sample_TT(const void* brdf_data, float u1, float u2)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        return hb->sample_TT(u1, u2);
    }

    static AtVector Hair_Sample_TRT(const void* brdf_data, float u1, float u2)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        return hb->sample_TRT(u1, u2);
    }

    static AtVector Hair_Sample_TRTg(const void* brdf_data, float u1, float u2)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        return hb->sample_TRTg(u1, u2);
    }

    static AtRGB Hair_Bsdf_R(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);
        AtRGB kfr[3];
        hairAttenuation(hb->sp.ior, geo.theta_d, geo.phi_d, hb->sp.absorption, kfr);
        return hb->bsdf_R(geo) * kfr[0] * hb->specular1Color;
    }
    static AtRGB Hair_Bsdf_TT(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);
        AtRGB kfr[3];
        hairAttenuation(hb->sp.ior, geo.theta_d, geo.phi_d, hb->sp.absorption, kfr);
        return hb->bsdf_TT(geo) * kfr[1] * hb->transmissionColor;
    }
    static AtRGB Hair_Bsdf_TRT(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);
        AtRGB kfr[3];
        hairAttenuation(hb->sp.ior, geo.theta_d, geo.phi_d, hb->sp.absorption, kfr);
        return hb->bsdf_TRT(geo) * kfr[2] * hb->specular2Color;
    }
    static AtRGB Hair_Bsdf_TRTg(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);
        AtRGB kfr[3];
        hairAttenuation(hb->sp.ior, geo.theta_d, geo.phi_d, hb->sp.absorption, kfr);
        return hb->bsdf_TRTg(geo) * kfr[2] * hb->specular2Color * hb->glintStrength;
    }

    static AtRGB HairGlossyBsdf(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        AtRGB result = AI_RGB_BLACK;

        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);

        AtRGB kfr[3];
        hairAttenuation(hb->sp.ior, geo.theta_d, geo.phi_d, hb->sp.absorption, kfr);
        
        AtRGB f_R = hb->do_glossy ? hb->bsdf_R(geo) * kfr[0] * hb->specular1Color : AI_RGB_BLACK;
        AtRGB f_TT = hb->do_glossy ? hb->bsdf_TT(geo) * kfr[1] * hb->transmissionColor : AI_RGB_BLACK;
        AtRGB f_TRT = hb->do_glossy ? hb->bsdf_TRT(geo) * kfr[2] * hb->specular2Color : AI_RGB_BLACK;
        AtRGB f_TRTg = hb->do_glossy ? hb->bsdf_TRTg(geo) * kfr[2] * hb->specular2Color * hb->glintStrength : AI_RGB_BLACK;

        // Store terms needed for MIS
        if (hb->is_bsdf_sample) 
        {
            hb->L_b = hb->_sg->Li;
            hb->d_b = hb->_sg->Ldist;
            hb->f_R_b = f_R;
            hb->f_TT_b = f_TT;
            hb->f_TRT_b = f_TRT;
            hb->f_TRTg_b = f_TRTg;
            hb->pdf_bl = 1.0f / (hb->_sg->we*hb->_sg->n);
            hb->w_b = *wi;
        }
        else 
        {
            hb->L_l = hb->_sg->Li;
            hb->d_l = hb->_sg->Ldist;
            hb->f_R_l = f_R;
            hb->f_TT_l = f_TT;
            hb->f_TRT_l = f_TRT;
            hb->f_TRTg_l = f_TRTg;
            hb->pdf_ll = 1.0f / (hb->_sg->we*hb->_sg->n);
            hb->w_l = *wi;
        }

        return f_R + f_TT + f_TRT + f_TRTg;
    }

    static AtVector HairGlossySample(const void* brdf_data, float u1, float u2)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        hb->is_bsdf_sample = true;
        if (u1 < 0.5f && u2 < 0.5f) 
        {
            u1 *= 2.0f;
            u2 *= 2.0f;
            return hb->sample_R(u1, u2);
        }
        else if (u1 < 0.5f && u2 >= 0.5f)
        {
            u1 *= 2.0f;
            u2 = 2.0f * (1.0f-u2);
            return hb->sample_TT(u1, u2);
        }
        else if (u1 >= 0.5f && u2 < 0.5f)
        {
            u1 = 2.0f * (1.0f-u1);
            u2 *= 2.0f;
            return hb->sample_TRT(u1, u2);
        }
        else
        {
            u2 = 2.0f * (1.0f-u1);
            u2 = 2.0f * (1.0f-u2);
            return hb->sample_TRTg(u1, u2);
        }


    }

    static float HairGlossyPdf(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;

        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);

        float p = (hb->pdf_R(geo) + hb->pdf_TT(geo) + hb->pdf_TRT(geo) + hb->pdf_TRTg(geo))*0.25f;

        if (hb->is_bsdf_sample) 
        {
            hb->pdf_bb = p;
        }
        else 
        {
            hb->pdf_lb = p;
        }

        return p;
    }

    static float Hair_Pdf_R(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);
        return hb->pdf_R(geo);
    }

    static float Hair_Pdf_TT(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);
        return hb->pdf_TT(geo);
    }

    static float Hair_Pdf_TRT(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);
        return hb->pdf_TRT(geo);
    }

    static float Hair_Pdf_TRTg(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);
        return hb->pdf_TRTg(geo);
    }
   
    /// Integrate the direct illumination for all diffuse and glossy lobes
    inline void integrateDirectMis(AtShaderGlobals* sg)
    {
        if (sg->Rt & AI_RAY_DIFFUSE) return;
        // Tell Arnold we want the full sphere for lighting.
        sg->fhemi = false;
        //sg->skip_shadow = true;
        AiLightsPrepare(sg);

        if ((sg->Rt & AI_RAY_CAMERA) && doLightGroups)
        {
            while (AiLightsGetSample(sg))
            {
                // per-light specular and diffuse strength multipliers
                float specular_strength = AiLightGetSpecular(sg->Lp);
                if (specular_strength < IMPORTANCE_EPS) continue;

                // get the group assigned to this light from the hash table using the light's pointer
                int lightGroup = data->lightGroups[sg->Lp];

                is_bsdf_sample = false;
                pdf_bb = 0.0f;
                AiEvaluateLightSample(sg, this, HairGlossySample, HairGlossyBsdf, HairGlossyPdf);

                float w = powerHeuristic(pdf_ll, pdf_lb);
                L_l = L_l * w / (pdf_ll*sg->n);
                result_R_direct += L_l * f_R_l;
                result_TT_direct += L_l * f_TT_l;
                result_TRT_direct += L_l * f_TRT_l;
                result_TRTg_direct += L_l * f_TRTg_l;

                if (pdf_bb != 0.0f)
                {
                    w = powerHeuristic(pdf_bb, pdf_bl);
                    L_b = L_b * w / (pdf_bb*sg->n);
                    result_R_direct += L_b * f_R_b;
                    result_TT_direct += L_b * f_TT_b;
                    result_TRT_direct += L_b * f_TRT_b;
                    result_TRTg_direct += L_b * f_TRTg_b;
                }

                // if the light is assigned a valid group number, add this sample's contribution to that light group
                if (lightGroup >= 0 && lightGroup < NUM_LIGHT_GROUPS)
                {
                    lightGroupsDirect[lightGroup] += L_b * (f_R_b + f_TT_b + f_TRT_b + f_TRTg_b);
                    lightGroupsDirect[lightGroup] += L_l * (f_R_l + f_TT_l + f_TRT_l + f_TRTg_l);
                }

            }
        }
        else
        {
            while (AiLightsGetSample(sg))
            {
                // per-light specular and diffuse strength multipliers
                float specular_strength = AiLightGetSpecular(sg->Lp);
                if (specular_strength < IMPORTANCE_EPS) continue;

                is_bsdf_sample = false;
                pdf_bb = 0.0f;
                AiEvaluateLightSample(sg, this, HairGlossySample, HairGlossyBsdf, HairGlossyPdf);

                float w = powerHeuristic(pdf_ll, pdf_lb);
                L_l = L_l * w / (pdf_ll*sg->n);
                result_R_direct += L_l * f_R_l;
                result_TT_direct += L_l * f_TT_l;
                result_TRT_direct += L_l * f_TRT_l;
                result_TRTg_direct += L_l * f_TRTg_l;

                if (pdf_bb != 0.0f)
                {
                    w = powerHeuristic(pdf_bb, pdf_bl);
                    L_b = L_b * w / (pdf_bb*sg->n);
                    result_R_direct += L_b * f_R_b;
                    result_TT_direct += L_b * f_TT_b;
                    result_TRT_direct += L_b * f_TRT_b;
                    result_TRTg_direct += L_b * f_TRTg_b;
                }

            }
        }
    
        sg->fhemi = true;        
    }

    inline void integrateDirect(AtShaderGlobals* sg)
    {
        sg->fhemi = false;
        
        if (do_glossy)
        {
            AiLightsPrepare(sg);
            AtRGB kfr[3];
            if ((sg->Rt & AI_RAY_CAMERA) && doLightGroups)
            {
                while (AiLightsGetSample(sg))
                {
                    // per-light specular and diffuse strength multipliers
                    float specular_strength = AiLightGetSpecular(sg->Lp);
                    if (specular_strength < IMPORTANCE_EPS) continue;

                    // get the group assigned to this light from the hash table using the light's pointer
                    int lightGroup = data->lightGroups[sg->Lp];

                    SctGeo geo(sg->Ld, theta_r, phi_r, U, V, W);

                    AtRGB kfr[3];
                    hairAttenuation(sp.ior, geo.theta_d, geo.phi_d, sp.absorption, kfr);
                    AtRGB L = sg->Li * sg->we;
                    AtRGB f_R = L * bsdf_R(geo) * kfr[0] * specular1Color;
                    AtRGB f_TT = L * bsdf_TT(geo) * kfr[1] * transmissionColor;
                    AtRGB f_TRT = L * bsdf_TRT(geo) * kfr[2] * specular2Color;
                    AtRGB f_TRTg = L * bsdf_TRTg(geo) * kfr[2] * specular2Color * glintStrength;

                    result_R_direct += f_R;
                    result_TT_direct += f_TT;
                    result_TRT_direct += f_TRT;
                    result_TRTg_direct += f_TRTg;

                    // if the light is assigned a valid group number, add this sample's contribution to that light group
                    if (lightGroup >= 0 && lightGroup < NUM_LIGHT_GROUPS)
                    {
                        lightGroupsDirect[lightGroup] += f_R + f_TT + f_TRT + f_TRTg;
                    }

                } // END light loop
            }
            else
            {
                while (AiLightsGetSample(sg))
                {
                    // per-light specular and diffuse strength multipliers
                    float specular_strength = AiLightGetSpecular(sg->Lp);
                    if (specular_strength < IMPORTANCE_EPS) continue;

                    SctGeo geo(sg->Ld, theta_r, phi_r, U, V, W);

                    AtRGB kfr[3];
                    hairAttenuation(sp.ior, geo.theta_d, geo.phi_d, sp.absorption, kfr);
                    AtRGB L = sg->Li * sg->we;
                    
                    result_R_direct += L * bsdf_R(geo) * kfr[0];
                    result_TT_direct += L * bsdf_TT(geo) * kfr[1];
                    result_TRT_direct += L * bsdf_TRT(geo) * kfr[2];
                    result_TRTg_direct += L * bsdf_TRTg(geo) * kfr[2];

                } // END light loop
                result_R_direct *= specular1Color;
                result_TT_direct *= transmissionColor;
                result_TRT_direct *= specular2Color;
                result_TRTg_direct *= specular2Color * glintStrength;
            }
        }

        if (do_diffuse && diffuseScatteringMix < 1.0f)
        {
            AiLightsPrepare(sg);
            while (AiLightsGetSample(sg))
            {
                // per-light specular and diffuse strength multipliers
                float diffuse_strength = AiLightGetDiffuse(sg->Lp);
                if (diffuse_strength < IMPORTANCE_EPS) continue;

                float cos_theta_i = AiV3Dot(U, sg->Ld);
                float sin_theta_i = sqrtf(MAX(1.0f - SQR(cos_theta_i), 0.0f));

                AtRGB L = sg->Li * sg->we * diffuse_strength * sin_theta_i;
                result_Pl_direct += L * hairColor * (1.0f / (4*AI_PI)) * (1.0f - diffuseScatteringMix);
            }
        }
        sg->fhemi = true;

        
    }

    inline void integrateDirectDual(AtShaderGlobals* sg)
    {
        float als_hairNumIntersections = 0;
        AtRGB T_f = AI_RGB_BLACK;
        AtRGB sigma_f = AI_RGB_BLACK;
        bool old_hemi = sg->fhemi;
        sg->fhemi = false;
        bool old_skipshadow = sg->skip_shadow;
        sg->skip_shadow = true;
        AiLightsPrepare(sg);
        AtRay ray;
        AtScrSample scrs;
        AtRGB kfr[3];
        AtRGB occlusion;
        float directFraction;

        if ((sg->Rt & AI_RAY_CAMERA) && doLightGroups)
        {
            while (AiLightsGetSample(sg))
            {
                SctGeo geo(sg->Ld, theta_r, phi_r, U, V, W);

                // per-light specular and diffuse strength multipliers
                float specular_strength = AiLightGetSpecular(sg->Lp);
                float diffuse_strength = AiLightGetDiffuse(sg->Lp);
                
                // get the group assigned to this light from the hash table using the light's pointer
                int lightGroup = data->lightGroups[sg->Lp];

                directFraction = 1.0f;
                occlusion = AI_RGB_WHITE;

                if (diffuse_strength > IMPORTANCE_EPS && diffuseScatteringMix > 0.0f)
                {
                    AiStateSetMsgInt("als_raytype", ALS_RAY_DUAL);
                    AiStateSetMsgFlt("als_hairNumIntersections", 0);
                    AiStateSetMsgRGB("als_T_f", AI_RGB_WHITE);
                    AiStateSetMsgRGB("als_sigma_bar_f", AI_RGB_BLACK);
                    AiMakeRay(&ray, AI_RAY_SHADOW, &(sg->P), &(sg->Ld), sg->Ldist, sg);
                    AiTrace(&ray, &scrs);
                    AiStateGetMsgFlt("als_hairNumIntersections", &als_hairNumIntersections);
                    AiStateGetMsgRGB("als_T_f", &T_f);
                    AiStateGetMsgRGB("als_sigma_bar_f", &sigma_f);
                    
                    directFraction = 1.0f - std::min(als_hairNumIntersections, float(numBlendHairs))/float(numBlendHairs);
                    occlusion = AI_RGB_WHITE - scrs.opacity;

                    AtRGB F_direct = directFraction * density_back*data->ds->f_back_direct(sp, geo); //< do f_s_direct separately for AOVs
                    AtRGB F_scatter = (1.0f-directFraction) * T_f * density_front * 
                                        (data->ds->f_s_scatter(sp, geo, sigma_f) + AI_PI * density_back * data->ds->f_back_scatter(sp, geo, sigma_f));

                    AtRGB f_Pl = sg->Li * sg->we * occlusion * F_direct * geo.cos_theta_i * AI_ONEOVERPI * diffuseColor;
                    AtRGB f_Pg = sg->Li * sg->we * occlusion * F_scatter * geo.cos_theta_i * AI_ONEOVERPI * diffuseColor;
                    result_Pl_direct += f_Pl * diffuseScatteringMix;
                    result_Pg_direct += f_Pg * diffuseScatteringMix;

                    // if the light is assigned a valid group number, add this sample's contribution to that light group
                    if (lightGroup >= 0 && lightGroup < NUM_LIGHT_GROUPS)
                    {
                        lightGroupsDirect[lightGroup] += f_Pl;
                        lightGroupsDirect[lightGroup] += f_Pg;
                    }
                }
                else
                {
                    AiStateSetMsgInt("als_raytype", ALS_RAY_UNDEFINED);
                    AiMakeRay(&ray, AI_RAY_SHADOW, &(sg->P), &(sg->Ld), sg->Ldist, sg);
                    AiTrace(&ray, &scrs);
                    occlusion = AI_RGB_WHITE - scrs.opacity;
                }

                if (do_glossy && specular_strength > IMPORTANCE_EPS)
                {
                    AtRGB kfr[3];
                    hairAttenuation(sp.ior, geo.theta_d, geo.phi_d, sp.absorption, kfr);
                    AtRGB L = sg->Li * sg->we * occlusion * directFraction;
                    AtRGB f_R = L * bsdf_R(geo) * kfr[0] * specular1Color;
                    AtRGB f_TT = L * bsdf_TT(geo) * kfr[1] * transmissionColor;
                    AtRGB f_TRT = L * bsdf_TRT(geo) * kfr[2] * specular2Color;
                    AtRGB f_TRTg = L * bsdf_TRTg(geo) * kfr[2] * specular2Color * glintStrength;

                    result_R_direct += f_R;
                    result_TT_direct += f_TT;
                    result_TRT_direct += f_TRT;
                    result_TRTg_direct += f_TRTg;

                    // if the light is assigned a valid group number, add this sample's contribution to that light group
                    if (lightGroup >= 0 && lightGroup < NUM_LIGHT_GROUPS)
                    {
                        lightGroupsDirect[lightGroup] += f_R + f_TT + f_TRT + f_TRTg;
                    }
                }

            } // END light loop
        }
        else
        {
            while (AiLightsGetSample(sg))
            {
                SctGeo geo(sg->Ld, theta_r, phi_r, U, V, W);

                // per-light specular and diffuse strength multipliers
                float specular_strength = AiLightGetSpecular(sg->Lp);
                float diffuse_strength = AiLightGetDiffuse(sg->Lp);

                directFraction = 1.0f;
                occlusion = AI_RGB_WHITE;

                if (diffuse_strength > IMPORTANCE_EPS && diffuseScatteringMix > 0.0f)
                {
                    AiStateSetMsgInt("als_raytype", ALS_RAY_DUAL);
                    AiStateSetMsgFlt("als_hairNumIntersections", 0);
                    AiStateSetMsgRGB("als_T_f", AI_RGB_WHITE);
                    AiStateSetMsgRGB("als_sigma_bar_f", AI_RGB_BLACK);
                    AiMakeRay(&ray, AI_RAY_SHADOW, &(sg->P), &(sg->Ld), sg->Ldist, sg);
                    AiTrace(&ray, &scrs);
                    AiStateGetMsgFlt("als_hairNumIntersections", &als_hairNumIntersections);
                    AiStateGetMsgRGB("als_T_f", &T_f);
                    AiStateGetMsgRGB("als_sigma_bar_f", &sigma_f);
                    
                    directFraction = 1.0f - std::min(als_hairNumIntersections, float(numBlendHairs))/float(numBlendHairs);
                    occlusion = AI_RGB_WHITE - scrs.opacity;

                    AtRGB F_direct = directFraction * density_back*data->ds->f_back_direct(sp, geo); //< do f_s_direct separately for AOVs
                    AtRGB F_scatter = (1.0f-directFraction) * T_f * density_front * 
                                        (data->ds->f_s_scatter(sp, geo, sigma_f) + AI_PI * density_back * data->ds->f_back_scatter(sp, geo, sigma_f));

                    result_Pl_direct += sg->Li * sg->we * occlusion * F_direct * geo.cos_theta_i * AI_ONEOVERPI;
                    result_Pg_direct += sg->Li * sg->we * occlusion * F_scatter * geo.cos_theta_i * AI_ONEOVERPI;
                }
                else
                {
                    AiStateSetMsgInt("als_raytype", ALS_RAY_UNDEFINED);
                    AiMakeRay(&ray, AI_RAY_SHADOW, &(sg->P), &(sg->Ld), sg->Ldist, sg);
                    AiTrace(&ray, &scrs);
                    occlusion = AI_RGB_WHITE - scrs.opacity;
                }

                if (do_glossy)
                {
                    AtRGB kfr[3];
                    hairAttenuation(sp.ior, geo.theta_d, geo.phi_d, sp.absorption, kfr);
                    AtRGB L = sg->Li * sg->we * occlusion * directFraction;
                    
                    result_R_direct += L * bsdf_R(geo) * kfr[0];
                    result_TT_direct += L * bsdf_TT(geo) * kfr[1];
                    result_TRT_direct += L * bsdf_TRT(geo) * kfr[2];
                    result_TRTg_direct += L * bsdf_TRTg(geo) * kfr[2];
                }

            } // END light loop
            result_R_direct *= specular1Color;
            result_TT_direct *= transmissionColor;
            result_TRT_direct *= specular2Color;
            result_TRTg_direct *= specular2Color * glintStrength;

            result_Pg_direct *= diffuseColor * diffuseScatteringMix;
            result_Pl_direct *= diffuseColor * diffuseScatteringMix;
        }
        AiStateSetMsgInt("als_raytype", ALS_RAY_UNDEFINED);

         sg->skip_shadow = old_skipshadow;

        if (do_diffuse && diffuseScatteringMix < 1.0f)
        {
            AiLightsPrepare(sg);
            while (AiLightsGetSample(sg))
            {
                // per-light specular and diffuse strength multipliers
                float diffuse_strength = AiLightGetDiffuse(sg->Lp);
                if (diffuse_strength < IMPORTANCE_EPS) continue;

                float cos_theta_i = AiV3Dot(U, sg->Ld);
                float sin_theta_i = sqrtf(MAX(1.0f - SQR(cos_theta_i), 0.0f));

                AtRGB L = sg->Li * sg->we * diffuse_strength * sin_theta_i;
                result_Pl_direct += L * hairColor * (1.0f / (4*AI_PI))  * (1.0f - diffuseScatteringMix);
            }
        }

        sg->fhemi = old_hemi;
    }

    inline void integrateDirectDualMis(AtShaderGlobals* sg)
    {
        float als_hairNumIntersections = 0;
        AtRGB T_f = AI_RGB_BLACK;
        AtRGB sigma_f = AI_RGB_BLACK;
        bool old_hemi = sg->fhemi;
        sg->fhemi = false;
        bool old_skipshadow = sg->skip_shadow;
        sg->skip_shadow = true;
        AiLightsPrepare(sg);
        AtRay ray;
        AtScrSample scrs;
        float directFraction;
        AtRGB occlusion;
        AtRGB F_direct;
        AtRGB F_scatter;
        
        AtRGB kfr[3];

        if ((sg->Rt & AI_RAY_CAMERA) && doLightGroups)
        {
            while (AiLightsGetSample(sg))
            {
                // per-light specular and diffuse strength multipliers
                float specular_strength = AiLightGetSpecular(sg->Lp);
                float diffuse_strength = AiLightGetDiffuse(sg->Lp);

                // get the group assigned to this light from the hash table using the light's pointer
                int lightGroup = data->lightGroups[sg->Lp];

                directFraction = 1.0f;
                occlusion = AI_RGB_WHITE;
                F_direct = F_scatter = AI_RGB_BLACK;

                is_bsdf_sample = false;
                pdf_bb = 0.0f;
                AiEvaluateLightSample(sg, this, HairGlossySample, HairGlossyBsdf, HairGlossyPdf);
                SctGeo geo_l(w_l, theta_r, phi_r, U, V, W);
                float w = powerHeuristic(pdf_ll, pdf_lb);
                
                if (diffuse_strength > IMPORTANCE_EPS)
                {
                    AiStateSetMsgInt("als_raytype", ALS_RAY_DUAL);
                    AiStateSetMsgFlt("als_hairNumIntersections", 0);
                    AiStateSetMsgRGB("als_T_f", AI_RGB_WHITE);
                    AiStateSetMsgRGB("als_sigma_bar_f", AI_RGB_BLACK);
                    AiMakeRay(&ray, AI_RAY_SHADOW, &(sg->P), &w_l, d_l, sg);
                    AiTrace(&ray, &scrs);
                    AiStateGetMsgFlt("als_hairNumIntersections", &als_hairNumIntersections);
                    AiStateGetMsgRGB("als_T_f", &T_f);
                    AiStateGetMsgRGB("als_sigma_bar_f", &sigma_f);
                    
                    directFraction = 1.0f - std::min(als_hairNumIntersections, float(numBlendHairs))/float(numBlendHairs);
                    occlusion = AI_RGB_WHITE - scrs.opacity;

                    F_direct = directFraction * density_back*data->ds->f_back_direct(sp, geo_l); //< do f_s_direct separately for AOVs
                    F_scatter = (1.0f-directFraction) * T_f * density_front * 
                                        (data->ds->f_s_scatter(sp, geo_l, sigma_f) + AI_PI * density_back * data->ds->f_back_scatter(sp, geo_l, sigma_f));  
                }
                else
                {
                    AiStateSetMsgInt("als_raytype", ALS_RAY_UNDEFINED);
                    AiMakeRay(&ray, AI_RAY_SHADOW, &(sg->P), &w_l, d_l, sg);
                    AiTrace(&ray, &scrs);
                    occlusion = AI_RGB_WHITE - scrs.opacity;
                }
                AtRGB L = L_l * w * occlusion / (pdf_ll*sg->n);

                AtRGB f_Pl = L * F_direct * geo_l.cos_theta_i * AI_ONEOVERPI * diffuseColor * diffuseScatteringMix;
                AtRGB f_Pg = L * F_scatter * geo_l.cos_theta_i * AI_ONEOVERPI * diffuseColor * diffuseScatteringMix;
                AtRGB f_R = L * directFraction * f_R_l * specular1Color;
                AtRGB f_TT = L * directFraction * f_TT_l * transmissionColor;
                AtRGB f_TRT = L * directFraction * f_TRT_l * specular2Color;
                AtRGB f_TRTg = L * directFraction * f_TRTg_l * specular2Color * glintStrength;

                result_Pl_direct += f_Pl;
                result_Pg_direct += f_Pg;
                result_R_direct += f_R;
                result_TT_direct += f_TT;
                result_TRT_direct += f_TRT;
                result_TRTg_direct += f_TRTg;

                // if the light is assigned a valid group number, add this sample's contribution to that light group
                if (lightGroup >= 0 && lightGroup < NUM_LIGHT_GROUPS)
                {
                    lightGroupsDirect[lightGroup] += f_Pl + f_Pg + f_R + f_TT + f_TRT + f_TRTg;
                }

                if (pdf_bb != 0.0f)
                {
                    SctGeo geo_b(w_b, theta_r, phi_r, U, V, W);
                    w = powerHeuristic(pdf_bb, pdf_bl);
                    
                    if (diffuse_strength > IMPORTANCE_EPS)
                    {
                        AiStateSetMsgInt("als_raytype", ALS_RAY_DUAL);
                        AiStateSetMsgFlt("als_hairNumIntersections", 0);
                        AiStateSetMsgRGB("als_T_f", AI_RGB_WHITE);
                        AiStateSetMsgRGB("als_sigma_bar_f", AI_RGB_BLACK);
                        AiMakeRay(&ray, AI_RAY_SHADOW, &(sg->P), &w_b, d_b, sg);
                        AiTrace(&ray, &scrs);
                        AiStateGetMsgFlt("als_hairNumIntersections", &als_hairNumIntersections);
                        AiStateGetMsgRGB("als_T_f", &T_f);
                        AiStateGetMsgRGB("als_sigma_bar_f", &sigma_f);
                    
                        directFraction = 1.0f - std::min(als_hairNumIntersections, float(numBlendHairs))/float(numBlendHairs);
                        occlusion = AI_RGB_WHITE - scrs.opacity;

                        F_direct = directFraction * density_back*data->ds->f_back_direct(sp, geo_b); //< do f_s_direct separately for AOVs
                        F_scatter = (1.0f-directFraction) * T_f * density_front * 
                                            (data->ds->f_s_scatter(sp, geo_b, sigma_f) + AI_PI * density_back * data->ds->f_back_scatter(sp, geo_b, sigma_f));
                    }
                    else
                    {
                        AiStateSetMsgInt("als_raytype", ALS_RAY_UNDEFINED);
                        AiMakeRay(&ray, AI_RAY_SHADOW, &(sg->P), &w_b, d_b, sg);
                        AiTrace(&ray, &scrs);
                        occlusion = AI_RGB_WHITE - scrs.opacity;
                    }
                    
                    L = L_b * w * occlusion / (pdf_bb*sg->n);
                    AtRGB f_Pl = L * F_direct * geo_b.cos_theta_i * AI_ONEOVERPI * diffuseColor;
                    AtRGB f_Pg = L * F_scatter * geo_b.cos_theta_i * AI_ONEOVERPI * diffuseColor;
                    AtRGB f_R = L * directFraction * f_R_b * specular1Color;
                    AtRGB f_TT = L * directFraction * f_TT_b * transmissionColor;
                    AtRGB f_TRT = L * directFraction * f_TRT_b * specular2Color;
                    AtRGB f_TRTg = L * directFraction * f_TRTg_b * specular2Color * glintStrength;

                    result_Pl_direct += f_Pl * diffuseScatteringMix;
                    result_Pg_direct += f_Pg * diffuseScatteringMix;
                    result_R_direct += f_R;
                    result_TT_direct += f_TT;
                    result_TRT_direct += f_TRT;
                    result_TRTg_direct += f_TRTg;

                    // if the light is assigned a valid group number, add this sample's contribution to that light group
                    if (lightGroup >= 0 && lightGroup < NUM_LIGHT_GROUPS)
                    {
                        lightGroupsDirect[lightGroup] += f_Pl + f_Pg + f_R + f_TT + f_TRT + f_TRTg;
                    }
                }

            } // END light loop
        }
        else
        {
            while (AiLightsGetSample(sg))
            {
                // per-light specular and diffuse strength multipliers
                float specular_strength = AiLightGetSpecular(sg->Lp);
                float diffuse_strength = AiLightGetDiffuse(sg->Lp);

                directFraction = 1.0f;
                occlusion = AI_RGB_WHITE;
                F_direct = F_scatter = AI_RGB_BLACK;

                is_bsdf_sample = false;
                pdf_bb = 0.0f;
                AiEvaluateLightSample(sg, this, HairGlossySample, HairGlossyBsdf, HairGlossyPdf);
                SctGeo geo_l(w_l, theta_r, phi_r, U, V, W);
                float w = powerHeuristic(pdf_ll, pdf_lb);
                
                if (diffuse_strength > IMPORTANCE_EPS)
                {
                    AiStateSetMsgInt("als_raytype", ALS_RAY_DUAL);
                    AiStateSetMsgFlt("als_hairNumIntersections", 0);
                    AiStateSetMsgRGB("als_T_f", AI_RGB_WHITE);
                    AiStateSetMsgRGB("als_sigma_bar_f", AI_RGB_BLACK);
                    AiMakeRay(&ray, AI_RAY_SHADOW, &(sg->P), &w_l, d_l, sg);
                    AiTrace(&ray, &scrs);
                    AiStateGetMsgFlt("als_hairNumIntersections", &als_hairNumIntersections);
                    AiStateGetMsgRGB("als_T_f", &T_f);
                    AiStateGetMsgRGB("als_sigma_bar_f", &sigma_f);
                    
                    directFraction = 1.0f - std::min(als_hairNumIntersections, float(numBlendHairs))/float(numBlendHairs);
                    occlusion = AI_RGB_WHITE - scrs.opacity;

                    F_direct = directFraction * density_back*data->ds->f_back_direct(sp, geo_l); //< do f_s_direct separately for AOVs
                    F_scatter = (1.0f-directFraction) * T_f * density_front * 
                                        (data->ds->f_s_scatter(sp, geo_l, sigma_f) + AI_PI * density_back * data->ds->f_back_scatter(sp, geo_l, sigma_f));  
                }
                else
                {
                    AiStateSetMsgInt("als_raytype", ALS_RAY_UNDEFINED);
                    AiMakeRay(&ray, AI_RAY_SHADOW, &(sg->P), &w_l, d_l, sg);
                    AiTrace(&ray, &scrs);
                    occlusion = AI_RGB_WHITE - scrs.opacity;
                }
                AtRGB L = L_l * w * occlusion / (pdf_ll*sg->n);
                result_Pl_direct += L * F_direct * geo_l.cos_theta_i * AI_ONEOVERPI;
                result_Pg_direct += L * F_scatter * geo_l.cos_theta_i * AI_ONEOVERPI;
                result_R_direct += L * directFraction * f_R_l;
                result_TT_direct += L * directFraction * f_TT_l;
                result_TRT_direct += L * directFraction * f_TRT_l;
                result_TRTg_direct += L * directFraction * f_TRTg_l;

                if (pdf_bb != 0.0f)
                {
                    SctGeo geo_b(w_b, theta_r, phi_r, U, V, W);
                    w = powerHeuristic(pdf_bb, pdf_bl);
                    
                    if (diffuse_strength > IMPORTANCE_EPS)
                    {
                        AiStateSetMsgInt("als_raytype", ALS_RAY_DUAL);
                        AiStateSetMsgFlt("als_hairNumIntersections", 0);
                        AiStateSetMsgRGB("als_T_f", AI_RGB_WHITE);
                        AiStateSetMsgRGB("als_sigma_bar_f", AI_RGB_BLACK);
                        AiMakeRay(&ray, AI_RAY_SHADOW, &(sg->P), &w_b, d_b, sg);
                        AiTrace(&ray, &scrs);
                        AiStateGetMsgFlt("als_hairNumIntersections", &als_hairNumIntersections);
                        AiStateGetMsgRGB("als_T_f", &T_f);
                        AiStateGetMsgRGB("als_sigma_bar_f", &sigma_f);
                    
                        directFraction = 1.0f - std::min(als_hairNumIntersections, float(numBlendHairs))/float(numBlendHairs);
                        occlusion = AI_RGB_WHITE - scrs.opacity;

                        F_direct = directFraction * density_back*data->ds->f_back_direct(sp, geo_b); //< do f_s_direct separately for AOVs
                        F_scatter = (1.0f-directFraction) * T_f * density_front * 
                                            (data->ds->f_s_scatter(sp, geo_b, sigma_f) + AI_PI * density_back * data->ds->f_back_scatter(sp, geo_b, sigma_f));
                    }
                    else
                    {
                        AiStateSetMsgInt("als_raytype", ALS_RAY_UNDEFINED);
                        AiMakeRay(&ray, AI_RAY_SHADOW, &(sg->P), &w_b, d_b, sg);
                        AiTrace(&ray, &scrs);
                        occlusion = AI_RGB_WHITE - scrs.opacity;
                    }
                    
                    L = L_b * w * occlusion / (pdf_bb*sg->n);
                    result_Pl_direct += L * F_direct * geo_b.cos_theta_i * AI_ONEOVERPI;
                    result_Pg_direct += L * F_scatter * geo_b.cos_theta_i * AI_ONEOVERPI;
                    result_R_direct += L * directFraction * f_R_b;
                    result_TT_direct += L * directFraction * f_TT_b;
                    result_TRT_direct += L * directFraction * f_TRT_b;
                    result_TRTg_direct += L * directFraction * f_TRTg_b;
                }

            } // END light loop
            result_R_direct *= specular1Color;
            result_TT_direct *= transmissionColor;
            result_TRT_direct *= specular2Color;
            result_TRTg_direct *= specular2Color * glintStrength;

            result_Pg_direct *= diffuseColor * diffuseScatteringMix;
            result_Pl_direct *= diffuseColor * diffuseScatteringMix;
        }
        AiStateSetMsgInt("als_raytype", ALS_RAY_UNDEFINED);

        sg->skip_shadow = old_skipshadow;

        if (do_diffuse && diffuseScatteringMix < 1.0f)
        {
            AiLightsPrepare(sg);
            while (AiLightsGetSample(sg))
            {
                // per-light specular and diffuse strength multipliers
                float diffuse_strength = AiLightGetDiffuse(sg->Lp);
                if (diffuse_strength < IMPORTANCE_EPS) continue;

                float cos_theta_i = AiV3Dot(U, sg->Ld);
                float sin_theta_i = sqrtf(MAX(1.0f - SQR(cos_theta_i), 0.0f));

                AtRGB L = sg->Li * sg->we * diffuse_strength * sin_theta_i;
                result_Pl_direct += L * hairColor * (1.0f / (4*AI_PI)) * (1.0f - diffuseScatteringMix);
            }
        }

        sg->fhemi = old_hemi;
        
    }


    /// Integrate the indirect illumination for all diffuse and glossy lobes
    inline void integrateIndirect(AtShaderGlobals* sg)
    {
        if (do_glossy && glossyIndirectStrength > 0.0f)
        {
            AiMakeRay(&wi_ray, AI_RAY_GLOSSY, &sg->P, NULL, AI_BIG, sg);

            sampit = AiSamplerIterator(data->sampler_glossy, sg);
            AiStateSetMsgInt("als_raytype", ALS_RAY_HAIR);
            while(AiSamplerGetSample(sampit, samples))
            {
                wi_ray.dir = HairGlossySample(this, samples[0], samples[1]);

                AtScrSample scrs;

                // trace our ray
                if (AiTrace(&wi_ray, &scrs))
                {
                    // calculate result
                    float p = HairGlossyPdf(this, &wi_ray.dir);
                    result_R_indirect += scrs.color * Hair_Bsdf_R(this, &wi_ray.dir) / p;
                    result_TT_indirect += scrs.color * Hair_Bsdf_TT(this, &wi_ray.dir) / p;
                    result_TRT_indirect += scrs.color * Hair_Bsdf_TRT(this, &wi_ray.dir) / p;
                    result_TRTg_indirect += scrs.color * Hair_Bsdf_TRTg(this, &wi_ray.dir) / p;
                }
            }
            AiStateSetMsgInt("als_raytype", ALS_RAY_UNDEFINED);
            float weight = AiSamplerGetSampleInvCount(sampit) * glossyIndirectStrength;
            result_R_indirect *= weight; //< TODO: factor of pi?
            result_TT_indirect *= weight; //< TODO: factor of pi?
            result_TRT_indirect *= weight; //< TODO: factor of pi?
            result_TRTg_indirect *= weight; //< TODO: factor of pi?
        }

        if (diffuseScatteringMix < 1.0f && diffuseIndirectStrength > 0.0f)
        {
            result_Pl_indirect = AiIndirectDiffuse(&wo, sg) * hairColor * (1.0f-diffuseScatteringMix) * diffuseIndirectStrength;
        }
    }

    inline void integrateIndirectDual(AtShaderGlobals* sg)
    {
        AiMakeRay(&wi_ray, AI_RAY_GLOSSY, &sg->P, NULL, AI_BIG, sg);
        float weight;

        bool do_glossy = true;
        if (sg->Rt & AI_RAY_DIFFUSE || sg->Rr > 0) do_glossy = false;

        if (do_glossy && glossyIndirectStrength > 0.0f)
        {
            AiStateSetMsgInt("als_raytype", ALS_RAY_HAIR);
            sampit = AiSamplerIterator(data->sampler_glossy, sg);
            while(AiSamplerGetSample(sampit, samples))
            {
                wi_ray.dir = HairGlossySample(this, samples[0], samples[1]);

                AtScrSample scrs;

                // trace our ray
                AiStateSetMsgFlt("alsPreviousRoughness", 1.0f);
                if (AiTrace(&wi_ray, &scrs))
                {
                    // calculate result
                    float p = HairGlossyPdf(this, &wi_ray.dir);
                    result_R_indirect += scrs.color * Hair_Bsdf_R(this, &wi_ray.dir) / p;
                    result_TT_indirect += scrs.color * Hair_Bsdf_TT(this, &wi_ray.dir) / p;
                    result_TRT_indirect += scrs.color * Hair_Bsdf_TRT(this, &wi_ray.dir) / p;
                    result_TRTg_indirect += scrs.color * Hair_Bsdf_TRTg(this, &wi_ray.dir) / p;
                }
            }
            AiStateSetMsgInt("als_raytype", ALS_RAY_UNDEFINED);
            weight = AiSamplerGetSampleInvCount(sampit) * glossyIndirectStrength;
            result_R_indirect *= weight; //< TODO: factor of pi?
            result_TT_indirect *= weight; //< TODO: factor of pi?
            result_TRT_indirect *= weight; //< TODO: factor of pi?
            result_TRTg_indirect *= weight; //< TODO: factor of pi?
        }

        if (diffuseIndirectStrength > 0.0f && diffuseScatteringMix > 0)
        {
            float als_hairNumIntersections = 0;
            AtRGB T_f = AI_RGB_BLACK;
            AtRGB sigma_f = AI_RGB_BLACK;
            sampit = AiSamplerIterator(data->sampler_diffuse, sg);
            AiStateSetMsgInt("als_raytype", ALS_RAY_DUAL);
            while (AiSamplerGetSample(sampit, samples))
            {
                if (samples[0] < 0.5f)
                {
                    samples[0] *= 2.0f;
                    wi_ray.dir = sample_Sb(samples[0], samples[1]);
                }
                else
                {
                    samples[0] = 1.0f - (2.0f * samples[0]);
                    wi_ray.dir = sample_Sf(samples[0], samples[1]);
                }

                SctGeo geo(wi_ray.dir, theta_r, phi_r, U, V, W);

                float p = 1.0f / ((pdf_Sf(geo) + pdf_Sb(geo)) * 0.5f);

                AiStateSetMsgFlt("als_hairNumIntersections", 0);
                AiStateSetMsgRGB("als_T_f", AI_RGB_WHITE);
                AiStateSetMsgRGB("als_sigma_bar_f", AI_RGB_BLACK);
                AiStateSetMsgFlt("alsPreviousRoughness", 1.0f);
                AiTrace(&wi_ray, &scrs);
                AiStateGetMsgFlt("als_hairNumIntersections", &als_hairNumIntersections);
                AiStateGetMsgRGB("als_T_f", &T_f);
                AiStateGetMsgRGB("als_sigma_bar_f", &sigma_f);

                float directFraction = 1.0f - std::min(als_hairNumIntersections, float(numBlendHairs))/float(numBlendHairs);
                AtRGB occlusion = AI_RGB_WHITE - scrs.opacity;

                AtRGB F_direct = directFraction * density_back*data->ds->f_back_direct(sp, geo); //< do f_s_direct separately for AOVs
                AtRGB F_scatter = (1.0f-directFraction) * T_f * density_front * 
                                    (data->ds->f_s_scatter(sp, geo, sigma_f) + AI_PI * density_back * data->ds->f_back_scatter(sp, geo, sigma_f));

                result_Pl_indirect += scrs.color * F_direct * geo.cos_theta_i * AI_ONEOVERPI;
                result_Pg_indirect += scrs.color * F_scatter * geo.cos_theta_i * AI_ONEOVERPI;

            }
            AiStateSetMsgInt("als_raytype", ALS_RAY_UNDEFINED);
            weight = AiSamplerGetSampleInvCount(sampit) * diffuseIndirectStrength;
            result_Pl_indirect *= weight * diffuseColor * diffuseScatteringMix;
            result_Pg_indirect *= weight * diffuseColor * diffuseScatteringMix;
        }

        if (diffuseScatteringMix < 1.0f && diffuseIndirectStrength > 0.0f)
        {
            result_Pl_indirect = AiIndirectDiffuse(&wo, sg) * hairColor * (1.0f-diffuseScatteringMix) * diffuseIndirectStrength;
        }
    }

    inline void writeResult(AtShaderGlobals* sg)
    {
        if (sg->Rt & AI_RAY_CAMERA)
        {
            AiAOVSetRGB(sg, data->aovs[k_diffuse_color].c_str(), hairColor);
            if ((result_Pg_direct + result_Pl_direct) != AI_RGB_BLACK) AiAOVSetRGB(sg, data->aovs[k_direct_diffuse].c_str(), result_Pg_direct + result_Pl_direct);
            if ((result_Pg_indirect + result_Pl_indirect) != AI_RGB_BLACK) AiAOVSetRGB(sg, data->aovs[k_indirect_diffuse].c_str(), result_Pg_indirect + result_Pl_indirect);
            if (result_R_direct != AI_RGB_BLACK) AiAOVSetRGB(sg, data->aovs[k_direct_specular].c_str(), result_R_direct);
            if (result_R_indirect != AI_RGB_BLACK) AiAOVSetRGB(sg, data->aovs[k_indirect_specular].c_str(), result_R_indirect);
            if (result_TRT_direct != AI_RGB_BLACK) AiAOVSetRGB(sg, data->aovs[k_direct_specular_2].c_str(), result_TRT_direct);
            if (result_TRT_indirect != AI_RGB_BLACK) AiAOVSetRGB(sg, data->aovs[k_indirect_specular_2].c_str(), result_TRT_indirect);
            if (result_TRTg_direct != AI_RGB_BLACK) AiAOVSetRGB(sg, data->aovs[k_direct_glint].c_str(), result_TRTg_direct);
            if (result_TRTg_indirect != AI_RGB_BLACK) AiAOVSetRGB(sg, data->aovs[k_indirect_glint].c_str(), result_TRTg_indirect);
            if (result_TT_direct != AI_RGB_BLACK) AiAOVSetRGB(sg, data->aovs[k_direct_transmission].c_str(), result_TT_direct);
            if (result_TT_indirect != AI_RGB_BLACK) AiAOVSetRGB(sg, data->aovs[k_indirect_transmission].c_str(), result_TT_indirect);
            if (result_Pg_direct != AI_RGB_BLACK) AiAOVSetRGB(sg, data->aovs[k_direct_global].c_str(), result_Pg_direct);
            if (result_Pl_direct != AI_RGB_BLACK) AiAOVSetRGB(sg, data->aovs[k_direct_local].c_str(), result_Pl_direct);
            if (result_Pg_indirect != AI_RGB_BLACK) AiAOVSetRGB(sg, data->aovs[k_indirect_global].c_str(), result_Pg_indirect);
            if (result_Pl_indirect != AI_RGB_BLACK) AiAOVSetRGB(sg, data->aovs[k_indirect_local].c_str(), result_Pl_indirect);
            for (int i=0; i < NUM_ID_AOVS; ++i)
            {
                AtRGB tmp;
                // check if output is enabled first in case we have an expensive network upstream
                if (AiAOVEnabled(data->aovs[k_id_1+i].c_str(), AI_TYPE_RGB))
                {
                    tmp = AiShaderEvalParamRGB(p_id1 + i);

                    // check if we're overriding it with a per-object id
                    if (AiNodeLookUpUserParameter(sg->Op, id_names[i]))
                    {
                        tmp = AiNodeGetRGB(sg->Op, id_names[i]);
                    }

                    if (tmp != AI_RGB_BLACK)
                        AiAOVSetRGB(sg, data->aovs[k_id_1+i].c_str(), tmp);
                }
            }

            if (doLightGroups)
            {
                for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
                {
                    if (lightGroupsDirect[i] != AI_RGB_BLACK)
                        AiAOVSetRGB(sg, data->aovs[k_light_group_1+i].c_str(), lightGroupsDirect[i]);
                }
            }

        }

        sg->out.RGB =   result_R_direct +
                        result_R_indirect +
                        result_TT_direct +
                        result_TT_indirect +
                        result_TRT_direct +
                        result_TRT_indirect +
                        result_TRTg_direct +
                        result_TRTg_indirect +
                        result_Pg_direct +
                        result_Pl_direct +
                        result_Pg_indirect +
                        result_Pl_indirect;
        sg->out_opacity = AI_RGB_WHITE;
    }

    AtVector U, V, W;   //< local coordinate frame
    float theta_r;      //< exitant spherical theta
    float phi_r;        //< existant spherical phi
    float cn;           //< random value per curve in [0,1)

    ScatteringParams sp; //< varying scattering parameters

    float diffuseStrength;
    AtRGB diffuseColor;
    float diffuseScatteringMix;
    AtRGB hairColor;
    AtRGB specular1Color;
    AtRGB specular2Color;
    AtRGB transmissionColor;
    float glintStrength;
    float diffuseIndirectStrength;
    float glossyIndirectStrength;

    float A_R;
    float B_R;
    float A_TT;
    float B_TT;
    float C_TT;
    float A_TRT;
    float B_TRT;
    float Cg;
    float Dg;

    float A_f;
    float A_b;
    float B_f;
    float B_b;

    // dual-scattering parameters
    int numBlendHairs;
    float density_front;
    float density_back;
    float colorNorm;

    AtRGB result_R_direct;
    AtRGB result_R_indirect;
    AtRGB result_TT_direct;
    AtRGB result_TT_indirect;
    AtRGB result_TRT_direct;
    AtRGB result_TRT_indirect;
    AtRGB result_TRTg_direct;
    AtRGB result_TRTg_indirect;

    // dual-scattering aovs
    AtRGB result_Pg_direct;
    AtRGB result_Pl_direct;
    AtRGB result_Pg_indirect;
    AtRGB result_Pl_indirect;

    // IDs
    AtRGB result_id1;
    AtRGB result_id2;
    AtRGB result_id3;
    AtRGB result_id4;
    AtRGB result_id5;
    AtRGB result_id6;
    AtRGB result_id7;
    AtRGB result_id8;

    bool do_diffuse;
    bool do_glossy;

    int depth;

    AtNode* node;
    ShaderData* data;
    AtRay wi_ray;
    AtScrSample scrs;
#if AI_VERSION_MAJOR_NUM > 0
    float samples[2];
#else
    double samples[2];
#endif
    AtSamplerIterator* sampit;

    AtVector wo;

    AtShaderGlobals* _sg;

    float specular1WidthScale;
    float specular2WidthScale;
    float transmissionWidthScale;

    bool dualImportant;

    // MIS temporaries
    // {
    bool is_bsdf_sample;
    AtRGB L_l;
    AtRGB f_R_l;
    AtRGB f_TT_l;
    AtRGB f_TRT_l;
    AtRGB f_TRTg_l;
    float pdf_ll;
    float pdf_lb;
    float d_l;
    AtVector w_l;
    AtRGB L_b;
    AtRGB f_R_b;
    AtRGB f_TT_b;
    AtRGB f_TRT_b;
    AtRGB f_TRTg_b;
    float pdf_bb;
    float pdf_bl;
    float d_b;
    AtVector w_b;
    // }

    bool doLightGroups;
    AtRGB lightGroupsDirect[NUM_LIGHT_GROUPS];
};

node_loader
{
   if (i>0) return 0;
   node->methods     = alHair;
   node->output_type = AI_TYPE_RGB;
   node->name        = "alHair";
   node->node_type   = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return true;
}

node_initialize
{
    HairBsdf::ShaderData* data = new HairBsdf::ShaderData;
    AiNodeSetLocalData(node, data);
}

node_finish
{
    if (AiNodeGetLocalData(node))
    {
        HairBsdf::ShaderData* data = (HairBsdf::ShaderData*)AiNodeGetLocalData(node);
        delete data;
    }
}

node_update
{
    HairBsdf::ShaderData* data = (HairBsdf::ShaderData*)AiNodeGetLocalData(node);
    data->update(params);

    AiAOVRegister("diffuse_color", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("direct_specular", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("indirect_specular", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("direct_specular_2", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("indirect_specular_2", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("direct_transmission", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("indirect_transmission", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("direct_glint", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("direct_global", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("indirect_global", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("direct_local", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("indirect_local", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_1", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_2", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_3", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_4", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_5", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_6", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_7", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_8", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_1", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_2", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_3", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_4", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_5", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_6", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_7", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_8", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);

    // Get all the light nodes in the scene and try and find their light group parameter
    // we'll store this based on the light pointer for fast access during rendering
    AtNodeIterator* it = AiUniverseGetNodeIterator(AI_NODE_LIGHT);
    while (!AiNodeIteratorFinished(it))
    {
        AtNode* light = AiNodeIteratorGetNext(it);
        if (AiNodeLookUpUserParameter(light, "lightGroup"))
            data->lightGroups[light] = AiNodeGetInt(light, "lightGroup") - 1;
        else 
            data->lightGroups[light] = -1;
    }
}

shader_evaluate
{
    
    // Get shader data
    HairBsdf::ShaderData* data = (HairBsdf::ShaderData*)AiNodeGetLocalData(node);

    // Create HairBsdf object 
    HairBsdf hb(node, sg, data);
    // Get parameters
    hb.evaluateParameters(sg, data);

    AtRGB opacity = AiShaderEvalParamRGB(p_opacity);
    float geo_opacity = 1.0f;
    if (AiUDataGetFlt("geo_opacity", &geo_opacity))
    {
        opacity *= geo_opacity;
    }

    float als_hairNumIntersections = 0;
    AtRGB als_T_f = AI_RGB_BLACK;
    AtRGB als_sigma_bar_f = AI_RGB_BLACK;
    bool do_dual = false;
    if (sg->Rr >= data->dual_depth) do_dual = true;

    int als_raytype = ALS_RAY_UNDEFINED;
    AiStateGetMsgInt("als_raytype", &als_raytype);
    
    if (do_dual && als_raytype == ALS_RAY_DUAL)
    {
        if (AiStateGetMsgFlt("als_hairNumIntersections", &als_hairNumIntersections) 
            && AiStateGetMsgRGB("als_T_f", &als_T_f)
            && AiStateGetMsgRGB("als_sigma_bar_f", &als_sigma_bar_f))
        {
            float theta_i = AI_PIOVER2 - sphericalTheta(sg->Rd, hb.U);

            AtRGB T_f = hb.data->ds->forward_attenuation(hb.sp, theta_i);
            //T_f = AI_RGB_WHITE - ((AI_RGB_WHITE - T_f)*opacity); //< modify transmission to account for opacity
            als_T_f *= T_f;
            AiStateSetMsgRGB("als_T_f", als_T_f);

            als_sigma_bar_f += hb.sp.beta_R2 + hb.sp.beta_TRT2 + hb.sp.beta_TT2;
            AiStateSetMsgRGB("als_sigma_bar_f", als_sigma_bar_f);

            als_hairNumIntersections+=minh(opacity);
            AiStateSetMsgFlt("als_hairNumIntersections", als_hairNumIntersections);

            if (maxh(als_T_f) > IMPORTANCE_EPS)
                sg->out_opacity = AI_RGB_BLACK;
            else
                sg->out_opacity = AI_RGB_WHITE;
        }
        else
        {
            sg->out_opacity = AI_RGB_WHITE;
        }

        return; // early out
    }

    // early-out regardless if we're in a shadow ray, or if opacity is zero
    if (sg->Rt & AI_RAY_SHADOW || AiShaderGlobalsApplyOpacity(sg, opacity)) return; 

    // early out if we're a hair-hair glossy ray and the ray depth says we should be calculating dual scattering only
    if (sg->Rr_gloss > data->dual_depth && als_raytype == ALS_RAY_HAIR) 
    {
        sg->out.RGB = AI_RGB_BLACK;
        sg->out_opacity = AI_RGB_WHITE;
        return;
    }
    
    // calculate scattering explicitly up to the dual depth cutoff.
    // in other words, do a brute force path trace for x=dual-depth bounces, then fall back to dual scattering for the rest.
    if (do_dual && hb.dualImportant && hb.diffuseScatteringMix > 0.0f)
    {
        if (data->doMis) hb.integrateDirectDualMis(sg);
        else hb.integrateDirectDual(sg);
        hb.integrateIndirectDual(sg);
    }
    else
    {
        if (data->doMis) hb.integrateDirectMis(sg);
        else hb.integrateDirect(sg);
        hb.integrateIndirect(sg);
    }

    // Write shader result
    hb.writeResult(sg);
    sg->out_opacity = opacity;
}


