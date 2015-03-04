#include "sss.h"
#include "norm_table.h"

float ScatteringParamsDirectional::_albedo_lut[SSS_ALBEDO_LUT_SZ] = {0.004660, 0.005024, 0.005285, 0.005503, 0.005697, 0.005875, 0.006041, 0.006197, 0.006345, 0.006488, 0.006625, 0.006758, 0.006887, 0.007013, 0.007136, 0.007256, 0.007374, 0.007490, 0.007604, 0.007716, 0.007826, 0.007935, 0.008043, 0.008149, 0.008255, 0.008358, 0.008461, 0.008564, 0.008665, 0.008765, 0.008865, 0.008964, 0.009063, 0.009160, 0.009257, 0.009354, 0.009450, 0.009546, 0.009641, 0.009735, 0.009829, 0.009923, 0.010016, 0.010110, 0.010203, 0.010295, 0.010387, 0.010479, 0.010571, 0.010662, 0.010754, 0.010845, 0.010935, 0.011026, 0.011116, 0.011207, 0.011297, 0.011387, 0.011477, 0.011566, 0.011656, 0.011746, 0.011835, 0.011925, 0.012014, 0.012103, 0.012192, 0.012281, 0.012370, 0.012459, 0.012548, 0.012637, 0.012726, 0.012815, 0.012904, 0.012993, 0.013082, 0.013171, 0.013260, 0.013349, 0.013438, 0.013527, 0.013616, 0.013705, 0.013794, 0.013883, 0.013972, 0.014061, 0.014150, 0.014240, 0.014329, 0.014418, 0.014508, 0.014597, 0.014687, 0.014777, 0.014866, 0.014956, 0.015046, 0.015136, 0.015226, 0.015317, 0.015407, 0.015497, 0.015588, 0.015678, 0.015769, 0.015860, 0.015951, 0.016042, 0.016133, 0.016224, 0.016315, 0.016407, 0.016498, 0.016590, 0.016682, 0.016774, 0.016866, 0.016958, 0.017051, 0.017143, 0.017236, 0.017329, 0.017421, 0.017514, 0.017608, 0.017701, 0.017794, 0.017888, 0.017982, 0.018076, 0.018170, 0.018264, 0.018359, 0.018453, 0.018548, 0.018643, 0.018738, 0.018833, 0.018928, 0.019024, 0.019120, 0.019215, 0.019311, 0.019408, 0.019504, 0.019601, 0.019698, 0.019795, 0.019892, 0.019990, 0.020087, 0.020185, 0.020283, 0.020381, 0.020479, 0.020577, 0.020676, 0.020774, 0.020873, 0.020973, 0.021072, 0.021172, 0.021272, 0.021372, 0.021472, 0.021573, 0.021673, 0.021774, 0.021875, 0.021976, 0.022078, 0.022180, 0.022282, 0.022383, 0.022486, 0.022588, 0.022690, 0.022794, 0.022897, 0.023001, 0.023104, 0.023208, 0.023311, 0.023417, 0.023521, 0.023626, 0.023730, 0.023837, 0.023942, 0.024047, 0.024153, 0.024260, 0.024366, 0.024472, 0.024579, 0.024687, 0.024794, 0.024901, 0.025010, 0.025118, 0.025226, 0.025335, 0.025444, 0.025554, 0.025663, 0.025773, 0.025883, 0.025992, 0.026104, 0.026214, 0.026326, 0.026436, 0.026548, 0.026661, 0.026772, 0.026886, 0.026997, 0.027111, 0.027225, 0.027338, 0.027452, 0.027567, 0.027682, 0.027795, 0.027911, 0.028027, 0.028143, 0.028259, 0.028376, 0.028491, 0.028609, 0.028726, 0.028844, 0.028962, 0.029082, 0.029201, 0.029320, 0.029439, 0.029559, 0.029678, 0.029800, 0.029921, 0.030041, 0.030164, 0.030286, 0.030409, 0.030531, 0.030655, 0.030778, 0.030903, 0.031028, 0.031151, 0.031277, 0.031403};
// float ScatteringParamsDirectional::_albedo_lut_d[SSS_ALBEDO_LUT_SZ] = {0.014215, 0.014978, 0.015518, 0.015965, 0.016359, 0.016716, 0.017047, 0.017357, 0.017650, 0.017930, 0.018198, 0.018457, 0.018707, 0.018949, 0.019185, 0.019415, 0.019639, 0.019858, 0.020073, 0.020284, 0.020491, 0.020694, 0.020894, 0.021092, 0.021286, 0.021478, 0.021667, 0.021854, 0.022039, 0.022221, 0.022402, 0.022581, 0.022758, 0.022934, 0.023107, 0.023280, 0.023451, 0.023620, 0.023789, 0.023956, 0.024121, 0.024286, 0.024450, 0.024612, 0.024774, 0.024934, 0.025094, 0.025252, 0.025410, 0.025567, 0.025724, 0.025879, 0.026034, 0.026188, 0.026341, 0.026494, 0.026646, 0.026797, 0.026948, 0.027098, 0.027248, 0.027397, 0.027545, 0.027694, 0.027841, 0.027988, 0.028135, 0.028281, 0.028427, 0.028572, 0.028717, 0.028862, 0.029006, 0.029150, 0.029294, 0.029437, 0.029580, 0.029722, 0.029864, 0.030006, 0.030148, 0.030289, 0.030430, 0.030571, 0.030712, 0.030852, 0.030992, 0.031132, 0.031272, 0.031411, 0.031550, 0.031689, 0.031828, 0.031967, 0.032105, 0.032243, 0.032381, 0.032519, 0.032657, 0.032795, 0.032932, 0.033069, 0.033207, 0.033344, 0.033481, 0.033617, 0.033754, 0.033891, 0.034027, 0.034163, 0.034300, 0.034436, 0.034572, 0.034708, 0.034844, 0.034980, 0.035115, 0.035251, 0.035386, 0.035522, 0.035657, 0.035793, 0.035928, 0.036063, 0.036199, 0.036334, 0.036469, 0.036604, 0.036739, 0.036874, 0.037009, 0.037144, 0.037279, 0.037414, 0.037549, 0.037684, 0.037819, 0.037954, 0.038089, 0.038224, 0.038358, 0.038493, 0.038628, 0.038763, 0.038898, 0.039033, 0.039167, 0.039302, 0.039437, 0.039572, 0.039707, 0.039842, 0.039977, 0.040112, 0.040247, 0.040382, 0.040517, 0.040652, 0.040787, 0.040922, 0.041058, 0.041193, 0.041328, 0.041463, 0.041599, 0.041734, 0.041869, 0.042005, 0.042140, 0.042276, 0.042412, 0.042547, 0.042683, 0.042819, 0.042955, 0.043091, 0.043227, 0.043363, 0.043499, 0.043635, 0.043771, 0.043907, 0.044044, 0.044180, 0.044317, 0.044485, 0.044622, 0.044759, 0.044896, 0.045033, 0.045170, 0.045307, 0.045445, 0.045582, 0.045719, 0.045857, 0.045995, 0.046132, 0.046270, 0.046408, 0.046546, 0.046684, 0.046822, 0.046960, 0.047099, 0.047237, 0.047376, 0.047514, 0.047653, 0.047792, 0.047931, 0.048070, 0.048209, 0.048348, 0.048523, 0.048663, 0.048803, 0.048943, 0.049083, 0.049223, 0.049363, 0.049503, 0.049643, 0.049784, 0.049924, 0.050065, 0.050206, 0.050347, 0.050526, 0.050667, 0.050809, 0.050950, 0.051092, 0.051234, 0.051376, 0.051518, 0.051660, 0.051802, 0.051984, 0.052127, 0.052270, 0.052413, 0.052556, 0.052699, 0.052842, 0.052986, 0.053170, 0.053314, 0.053458, 0.053602, 0.053747, 0.053891, 0.054077, 0.054222, 0.054367, 0.054513};
ScatteringParamsDirectional::ScatteringParamsDirectional(AtRGB sigma_s_prime, AtRGB sigma_a, float g) :
g(g), eta(1.3f), 
C_phi(0.175626f),
C_phi_inv(1.03668f),
C_E(0.27735f),
_3C2(0.0611156f),
A(2.05736f)
{
    AtRGB sigma_s = sigma_s_prime/(1.0f-g);

    sigma_t_prime = sigma_s_prime + sigma_a;
    sigma_t = sigma_s + sigma_a;

    alpha_prime = sigma_s_prime / sigma_t_prime;

    float apsum = alpha_prime.r + alpha_prime.g + alpha_prime.b;
    albedo_norm = alpha_prime / apsum;

    // AtRGB sq = sqrt(3.0f * (AI_RGB_WHITE - alpha_prime));
    // albedo = alpha_prime * 0.5f * (AI_RGB_WHITE + fast_exp(-A*sq*4.0f/3.0f)) * fast_exp(-sq);



    D = (2*sigma_t_prime) / (3*SQR(sigma_t_prime));
    sigma_tr = sqrt(sigma_a / D);
    de = 2.131 * D / sqrt(alpha_prime);
    zr = AI_RGB_WHITE / sigma_t_prime;
}

ScatteringParamsDirectional::ScatteringParamsDirectional(AtRGB c, float scale, float g, bool norm_color, bool dir) :
g(g), eta(1.3f), 
C_phi(0.175626f),
C_phi_inv(1.03668f),
C_E(0.27735f),
_3C2(0.0611156f),
A(2.05736f)
{
    AtRGB c_n = c;
    if (norm_color) c_n /= maxh(c);
    AtRGB ap = rgb(
        computeAlphaPrime(c_n.r * 0.439f),
        computeAlphaPrime(c_n.g * 0.439f),
        computeAlphaPrime(c_n.b * 0.439f)
    );

    AtRGB str = AI_RGB_WHITE / c_n;

    AtRGB stp = rgb(str[0] / ( sqrt( 3 * ( 1 - ap[0] ) ) ),
                    str[1] / ( sqrt( 3 * ( 1 - ap[1] ) ) ),
                    str[2] / ( sqrt( 3 * ( 1 - ap[2] ) ) ) );

    AtRGB sigma_s_prime = stp * ap ;
    AtRGB sigma_a = stp - sigma_s_prime;

    // factor of 2.5 is eyeballed to roughly match the look of the cubic
    sigma_s_prime *= scale * AI_PI * 2.5;
    sigma_a *= scale * AI_PI * 2.5;

    AtRGB sigma_s = sigma_s_prime/(1.0f-g);

    sigma_t_prime = sigma_s_prime + sigma_a;
    sigma_t = sigma_s + sigma_a;

    alpha_prime = sigma_s_prime / sigma_t_prime;

    float apsum = alpha_prime.r + alpha_prime.g + alpha_prime.b;
    albedo_norm = alpha_prime / apsum;

    D = (2*sigma_t_prime) / (3*SQR(sigma_t_prime));
    sigma_tr = sqrt(sigma_a / D);
    de = 2.131 * D / sqrt(alpha_prime);
    zr = AI_RGB_WHITE / sigma_t_prime;

    float ld = maxh(zr);
    float maxdist = ld * SSS_MAX_RADIUS;

    
    int r_idx = int(c_n.r * (SSS_ALBEDO_LUT_SZ-1));
    int g_idx = int(c_n.g * (SSS_ALBEDO_LUT_SZ-1));
    int b_idx = int(c_n.b * (SSS_ALBEDO_LUT_SZ-1));
    // if (dir)
    // {
    //     albedo.r = _albedo_lut_d[r_idx];
    //     albedo.g = _albedo_lut_d[g_idx];
    //     albedo.b = _albedo_lut_d[b_idx];
    // }
    // else
    // {
        albedo.r = _albedo_lut[r_idx];
        albedo.g = _albedo_lut[g_idx];
        albedo.b = _albedo_lut[b_idx];
    // }
}

void alsIrradiateSample(AtShaderGlobals* sg, DirectionalMessageData* dmd, bool directional)
{
    AiStateSetMsgInt("als_raytype", ALS_RAY_UNDEFINED);
    DiffusionSample& samp = dmd->samples[dmd->sss_depth];
    // void* brdf_data = AiOrenNayarMISCreateData(sg, 0.0f);
    // sg->fhemi = false;
    AiLightsPrepare(sg);
    AtRGB result_diffuse = AI_RGB_BLACK;
    AtUInt32 old_fi = sg->fi;
    samp.Rd = AI_RGB_BLACK;
    while (AiLightsGetSample(sg))
    {
        // can't use MIS here because Arnold cocks up the shadowing ;__;
        // result_diffuse += AiEvaluateLightSample(sg, brdf_data, AiOrenNayarMISSample, AiOrenNayarMISBRDF, AiOrenNayarMISPDF);
        if (directional)
        {
            samp.Rd += directionalDipole(sg->P, sg->N, dmd->Po, dmd->No, sg->Ld, dmd->wo, dmd->sp)
                    * sg->Li * MAX(AiV3Dot(sg->Ld, sg->N), 0.0f) * AI_ONEOVERPI * sg->we;
        }
        else
        {
            result_diffuse += sg->Li * MAX(AiV3Dot(sg->Ld, sg->N), 0.0f) * AI_ONEOVERPI * sg->we;
        }

    }

    if (!directional)
        samp.Rd += result_diffuse * directionalDipole(sg->P, sg->N, dmd->Po, dmd->No, sg->N, dmd->No, dmd->sp);
    
    samp.Rd += AiIndirectDiffuse(&sg->N, sg) * directionalDipole(sg->P, sg->N, dmd->Po, dmd->No, sg->N, dmd->No, dmd->sp);

    // printf("samp.Rd: (%f, %f, %f)\n", samp.Rd.r, samp.Rd.g, samp.Rd.b);
    

    samp.S = sg->P - dmd->Po;
    samp.r = AiV3Length(samp.S);
    if (!AiColorIsZero(result_diffuse))
    {
        // samp.R = dipole(samp.r, dmd->No, sg->N, dmd->sp);
        samp.b = 1.0f ;//- (1.0f - AiV3Dot(sg->N, dmd->No) * (1.0f - AiV3Dot(dmd->No, AiV3Normalize(samp.S))))*0.25f;
    }

    samp.N = sg->N;
    samp.Ng = sg->Ng;
    samp.P = sg->P;

    dmd->sss_depth++;

    dmd->maxdist -= samp.r;

    if (dmd->sss_depth < SSS_MAX_SAMPLES && dmd->maxdist > 0.0f)
    {
        AiStateSetMsgInt("als_raytype", ALS_RAY_SSS);
        
        AtRay ray;
        AtScrSample scrs;
        sg->Rr--;
        AiMakeRay(&ray, AI_RAY_SUBSURFACE, &sg->P, &sg->Rd, dmd->maxdist, sg);
        AiTrace(&ray, &scrs);
    }
}

AtRGB alsDiffusion(AtShaderGlobals* sg, DirectionalMessageData* dmd, AtSampler* sss_sampler, 
                     AtRGB sssRadiusColor, float sssRadius, float sssDensityScale, bool directional)
{
    AtVector U, V;
    AiBuildLocalFrameShirley(&U, &V, &sg->Ng);

    // Get scattering parameters from supplied scattering colour and mfp
    // AtRGB sigma_s_prime, sigma_a;
#if 1 
    // alphaInversion(sssRadiusColor, sssRadius, sigma_s_prime, sigma_a);
    // sigma_s_prime *= sssDensityScale / sssRadius * AI_PI;
    // sigma_a *= sssDensityScale / sssRadius * AI_PI;
#else
    sigma_s_prime = rgb(1.09f, 1.59f, 1.79f);
    sigma_a = rgb(0.013, 0.07f, 0.145f);

    // The 10 / AI_PIOVER2 factor is there to roughly match the behaviour of the cubic. And if you're
    // going to use a magic number it should always involve pi somewhere...
    sigma_s_prime *= sssDensityScale*20.0f / AI_PIOVER2;
    sigma_a *= sssDensityScale*20.0f / AI_PIOVER2;
#endif    
    // float eta = 1.3f;

    // ScatteringParamsDirectional sp(sigma_s_prime, sigma_a, 0.0f);
    ScatteringParamsDirectional sp(sssRadiusColor, sssDensityScale/sssRadius, 0.0f, true, directional);
    dmd->sp = sp;
#if 1
    // static bool first = true;
    // if (0)
    // {
    //     printf("alpha_prime: (%f, %f, %f)\n", sp.alpha_prime.r, sp.alpha_prime.g, sp.alpha_prime.b);
    //     printf("de: (%f, %f, %f)\n", sp.de.r, sp.de.g, sp.de.b);
    //     printf("D: (%f, %f, %f)\n", sp.D.r, sp.D.g, sp.D.b);
    //     printf("sigma_s_prime: (%f, %f, %f)\n", sigma_s_prime.r, sigma_s_prime.g, sigma_s_prime.b);
    //     printf("sigma_a: (%f, %f, %f)\n", sigma_a.r, sigma_a.g, sigma_a.b);
    //     printf("sigma_t_prime: (%f, %f, %f)\n", sp.sigma_t_prime.r, sp.sigma_t_prime.g, sp.sigma_t_prime.b);
    //     AtRGB st = exp( - sp.sigma_t_prime * .002f);
    //     printf("st: (%f, %f, %f)\n", st.r, st.g, st.b);
    //     first = false;
    // }
#endif 

    // Find the max component of the mfp
    float l = maxh(sp.zr);

    // trick Arnold into thinking we're shooting from a different face than we actually are so he doesn't ignore intersections
    AtUInt32 old_fi = sg->fi;
    sg->fi = UINT_MAX;

    AtRGB result_sss = AI_RGB_BLACK;
    
    // Set our maximum sample distance to be some multiple of the mfp
    float R_max = l * SSS_MAX_RADIUS;
    // sp.albedo = integrateDirectional(sp, R_max, 50);
    // sp.albedo = directionalNorm(sssRadiusColor);
    // sp.albedo = AI_RGB_WHITE;
    
    AtRGB Rd_sum = AI_RGB_BLACK;
    int samplesTaken = 0;
    float samples[2];
    AtRay wi_ray;
    AtScrSample scrs;
    AtSamplerIterator* sampit = AiSamplerIterator(sss_sampler, sg);
    dmd->wo = -sg->Rd;
    while (AiSamplerGetSample(sampit, samples))
    {
        // TODO: replace with a better sampling scheme
        //wi_ray.dir = uniformSampleSphere(samples[0], samples[1]);
        float dx, dy;

        AtVector Wsss, Usss, Vsss, Usss_1, Vsss_1, Usss_2, Vsss_2;
        float c_axis = 1.0f, c_axis_1, c_axis_2;;
        if (samples[0] < 0.5f)
        {
            samples[0] *= 2.0f;
            c_axis = 0.5f;
            c_axis_1 = 0.25f;
            c_axis_2 = 0.25f;
            Wsss = sg->Ng;
            
            Usss = U;
            Vsss = V;

            Usss_1 = sg->Ng;
            Vsss_1 = V;

            Usss_2 = U;
            Vsss_2 = sg->Ng;
        }
        else if (samples[0] < 0.75f)
        {
            samples[0] = (samples[0] - 0.5f) * 4.0f;
            c_axis = 0.25f;
            c_axis_1 = 0.5f;
            c_axis_2 = 0.25f;
            Wsss = U;
            
            Usss = sg->Ng;
            Vsss = V;

            Usss_1 = U;
            Vsss_1 = V;

            Usss_2 = U;
            Vsss_2 = sg->Ng;
        }
        else
        {
            samples[0] = (1.0f-samples[0])* 4.0f;
            c_axis = 0.25f;
            c_axis_1 = 0.25f;
            c_axis_2 = 0.5f;
            Wsss = V;
            
            Usss = U;
            Vsss = sg->Ng;

            Usss_1 = sg->Ng;
            Vsss_1 = V;

            Usss_2 = U;
            Vsss_2 = V;
        }

        AtVector Pd;
        float pdf_disk_a0_c0, pdf_disk_1, pdf_disk_2;
        float r_disk;
        float c_disk = 1.0f, c_disk_1, c_disk_2;
        float sigma, sigma_1, sigma_2;
        if (samples[1] < sp.albedo_norm.r)
        {
            samples[1] /= sp.albedo_norm.r;
            sigma = sp.sigma_tr.r;
            sigma_1 = sp.sigma_tr.g;
            sigma_2 = sp.sigma_tr.b;
            pdf_disk_a0_c0 = diffusionSampleDisk(samples[0], samples[1], sigma, dx, dy, r_disk);

            c_disk = sp.albedo_norm.r;
            c_disk_1 = sp.albedo_norm.g;
            c_disk_2 = sp.albedo_norm.b;
        }
        else if (samples[1] < (sp.albedo_norm.r + sp.albedo_norm.g))
        {
            samples[1] -= sp.albedo_norm.r;
            samples[1] /= sp.albedo_norm.g;
            sigma = sp.sigma_tr.g;
            sigma_1 = sp.sigma_tr.r;
            sigma_2 = sp.sigma_tr.b;
            pdf_disk_a0_c0 = diffusionSampleDisk(samples[0], samples[1], sp.sigma_tr.g, dx, dy, r_disk);
            c_disk = sp.albedo_norm.g;
            c_disk_1 = sp.albedo_norm.r;
            c_disk_2 = sp.albedo_norm.b;
        }
        else
        {
            samples[1] = 1.0f - samples[1];
            samples[1] /= sp.albedo_norm.b;
            sigma = sp.sigma_tr.b;
            sigma_1 = sp.sigma_tr.g;
            sigma_2 = sp.sigma_tr.r;
            pdf_disk_a0_c0 = diffusionSampleDisk(samples[0], samples[1], sp.sigma_tr.b, dx, dy, r_disk);   
            c_disk = sp.albedo_norm.b;
            c_disk_1 = sp.albedo_norm.g;
            c_disk_2 = sp.albedo_norm.r;
        }

        AtVector dir = -Wsss;
        float dz = R_max;
        AtPoint origin = sg->P + Wsss*(dz*0.5f) + Usss * dx + Vsss * dy;
        float maxdist = R_max * 2.0f;

        AiMakeRay(&wi_ray, AI_RAY_SUBSURFACE, &origin, &dir, maxdist, sg);

        AiStateSetMsgInt("als_raytype", ALS_RAY_SSS);
        dmd->sss_depth = 0;
        dmd->maxdist = R_max;
        dmd->Po = sg->P;
        dmd->No = sg->N;
        if (AiTrace(&wi_ray, &scrs))
        {            
            for (int i=0; i < dmd->sss_depth; ++i)
            {                
                if (AiColorIsZero(dmd->samples[i].Rd)) continue;

                float geom = fabsf(AiV3Dot(dmd->samples[i].Ng, Wsss));
                float geom_1 = fabsf(AiV3Dot(dmd->samples[i].Ng, Usss));
                float geom_2 = fabsf(AiV3Dot(dmd->samples[i].Ng, Vsss));

                float r_u_1 = AiV3Dot(dmd->samples[i].S, Usss_1);
                float r_v_1 = AiV3Dot(dmd->samples[i].S, Vsss_1);
                float r_1 = sqrtf(SQR(r_u_1)+SQR(r_v_1));

                float r_u_2 = AiV3Dot(dmd->samples[i].S, Usss_2);
                float r_v_2 = AiV3Dot(dmd->samples[i].S, Vsss_2);
                float r_2 = sqrtf(SQR(r_u_2)+SQR(r_v_2));

                float pdf_disk_a0_c0 = diffusionPdf(r_disk, sigma) * geom;
                float pdf_disk_a0_c1 = diffusionPdf(r_disk, sigma_1) * geom;
                float pdf_disk_a0_c2 = diffusionPdf(r_disk, sigma_2) * geom;

                float pdf_disk_a1_c0 = diffusionPdf(r_1, sigma) * geom_1;
                float pdf_disk_a1_c1 = diffusionPdf(r_1, sigma_1) * geom_1;
                float pdf_disk_a1_c2 = diffusionPdf(r_1, sigma_2) * geom_1;

                float pdf_disk_a2_c0 = diffusionPdf(r_2, sigma) * geom_2;
                float pdf_disk_a2_c1 = diffusionPdf(r_2, sigma_1) * geom_2;
                float pdf_disk_a2_c2 = diffusionPdf(r_2, sigma_2) * geom_2;

                float pdf_sum = 
                    pdf_disk_a0_c0 * c_disk * c_axis +
                    pdf_disk_a0_c1 * c_disk_1 * c_axis +
                    pdf_disk_a0_c2 * c_disk_2 * c_axis +
                    pdf_disk_a1_c0 * c_disk * c_axis_1 +
                    pdf_disk_a1_c1 * c_disk_1 * c_axis_1 +
                    pdf_disk_a1_c2 * c_disk_2 * c_axis_1 +
                    pdf_disk_a2_c0 * c_disk * c_axis_2 +
                    pdf_disk_a2_c1 * c_disk_1 * c_axis_2 +
                    pdf_disk_a2_c2 * c_disk_2 * c_axis_2;

                // result_sss += dmd->samples[i].E * dmd->samples[i].R * dmd->samples[i].b * r_disk / pdf_sum;
                result_sss += dmd->samples[i].Rd * r_disk / pdf_sum;
                // Rd_sum += dmd->samples[i].R * dmd->samples[i].b * r_disk / pdf_sum;
                
            }
            
        }
    }
    float w = AiSamplerGetSampleInvCount(sampit);
    result_sss *= w;
    // Rd_sum *= w;
    // if (!AiColorIsZero(Rd_sum)) result_sss *= (sp.albedo/Rd_sum); //< Peter Kutz's dark-edge normalization trick
    result_sss /= sp.albedo;
    sg->fi = old_fi;

    // Optimization hack: do a regular indirect diffuse and colour it like the subsurface instead of allowing sss rays to
    // continue for another diffuse bounce.
    /*
    if (!sssDoIndirect)
    {
        AtRGB result_sss_indirect = AiIndirectDiffuse(&sg->N, sg); 
        if (!data->sss_normalize) result_sss_indirect *= sp.albedo;
        result_sss += result_sss_indirect;
    }
    */

    return result_sss;
}