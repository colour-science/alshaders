#include "sss.h"
#include "norm_table.h"

float ScatteringParamsDirectional::_albedo_lut[SSS_ALBEDO_LUT_SZ] = {0.004660, 0.005024, 0.005285, 0.005503, 0.005697, 0.005875, 0.006041, 0.006197, 0.006345, 0.006488, 0.006625, 0.006758, 0.006887, 0.007013, 0.007136, 0.007256, 0.007374, 0.007490, 0.007604, 0.007716, 0.007826, 0.007935, 0.008043, 0.008149, 0.008255, 0.008358, 0.008461, 0.008564, 0.008665, 0.008765, 0.008865, 0.008964, 0.009063, 0.009160, 0.009257, 0.009354, 0.009450, 0.009546, 0.009641, 0.009735, 0.009829, 0.009923, 0.010016, 0.010110, 0.010203, 0.010295, 0.010387, 0.010479, 0.010571, 0.010662, 0.010754, 0.010845, 0.010935, 0.011026, 0.011116, 0.011207, 0.011297, 0.011387, 0.011477, 0.011566, 0.011656, 0.011746, 0.011835, 0.011925, 0.012014, 0.012103, 0.012192, 0.012281, 0.012370, 0.012459, 0.012548, 0.012637, 0.012726, 0.012815, 0.012904, 0.012993, 0.013082, 0.013171, 0.013260, 0.013349, 0.013438, 0.013527, 0.013616, 0.013705, 0.013794, 0.013883, 0.013972, 0.014061, 0.014150, 0.014240, 0.014329, 0.014418, 0.014508, 0.014597, 0.014687, 0.014777, 0.014866, 0.014956, 0.015046, 0.015136, 0.015226, 0.015317, 0.015407, 0.015497, 0.015588, 0.015678, 0.015769, 0.015860, 0.015951, 0.016042, 0.016133, 0.016224, 0.016315, 0.016407, 0.016498, 0.016590, 0.016682, 0.016774, 0.016866, 0.016958, 0.017051, 0.017143, 0.017236, 0.017329, 0.017421, 0.017514, 0.017608, 0.017701, 0.017794, 0.017888, 0.017982, 0.018076, 0.018170, 0.018264, 0.018359, 0.018453, 0.018548, 0.018643, 0.018738, 0.018833, 0.018928, 0.019024, 0.019120, 0.019215, 0.019311, 0.019408, 0.019504, 0.019601, 0.019698, 0.019795, 0.019892, 0.019990, 0.020087, 0.020185, 0.020283, 0.020381, 0.020479, 0.020577, 0.020676, 0.020774, 0.020873, 0.020973, 0.021072, 0.021172, 0.021272, 0.021372, 0.021472, 0.021573, 0.021673, 0.021774, 0.021875, 0.021976, 0.022078, 0.022180, 0.022282, 0.022383, 0.022486, 0.022588, 0.022690, 0.022794, 0.022897, 0.023001, 0.023104, 0.023208, 0.023311, 0.023417, 0.023521, 0.023626, 0.023730, 0.023837, 0.023942, 0.024047, 0.024153, 0.024260, 0.024366, 0.024472, 0.024579, 0.024687, 0.024794, 0.024901, 0.025010, 0.025118, 0.025226, 0.025335, 0.025444, 0.025554, 0.025663, 0.025773, 0.025883, 0.025992, 0.026104, 0.026214, 0.026326, 0.026436, 0.026548, 0.026661, 0.026772, 0.026886, 0.026997, 0.027111, 0.027225, 0.027338, 0.027452, 0.027567, 0.027682, 0.027795, 0.027911, 0.028027, 0.028143, 0.028259, 0.028376, 0.028491, 0.028609, 0.028726, 0.028844, 0.028962, 0.029082, 0.029201, 0.029320, 0.029439, 0.029559, 0.029678, 0.029800, 0.029921, 0.030041, 0.030164, 0.030286, 0.030409, 0.030531, 0.030655, 0.030778, 0.030903, 0.031028, 0.031151, 0.031277, 0.031403};
float ScatteringParamsDirectional::_albedo_lut_d[SSS_ALBEDO_LUT_SZ] = {0.006032, 0.006421, 0.006700, 0.006933, 0.007138, 0.007327, 0.007501, 0.007664, 0.007821, 0.007969, 0.008113, 0.008252, 0.008387, 0.008518, 0.008645, 0.008770, 0.008892, 0.009012, 0.009129, 0.009246, 0.009359, 0.009472, 0.009583, 0.009692, 0.009800, 0.009906, 0.010012, 0.010120, 0.010223, 0.010323, 0.010425, 0.010528, 0.010627, 0.010726, 0.010828, 0.010926, 0.011022, 0.011122, 0.011217, 0.011311, 0.011407, 0.011505, 0.011600, 0.011694, 0.011788, 0.011882, 0.011974, 0.012067, 0.012160, 0.012252, 0.012345, 0.012437, 0.012530, 0.012622, 0.012712, 0.012800, 0.012894, 0.012984, 0.013075, 0.013164, 0.013257, 0.013346, 0.013434, 0.013523, 0.013614, 0.013702, 0.013790, 0.013882, 0.013972, 0.014061, 0.014150, 0.014238, 0.014327, 0.014417, 0.014505, 0.014593, 0.014683, 0.014773, 0.014861, 0.014949, 0.015038, 0.015124, 0.015214, 0.015301, 0.015390, 0.015480, 0.015570, 0.015658, 0.015746, 0.015836, 0.015924, 0.016011, 0.016100, 0.016190, 0.016280, 0.016368, 0.016457, 0.016544, 0.016634, 0.016722, 0.016810, 0.016899, 0.016991, 0.017081, 0.017171, 0.017261, 0.017350, 0.017439, 0.017528, 0.017617, 0.017708, 0.017798, 0.017888, 0.017978, 0.018068, 0.018157, 0.018247, 0.018337, 0.018428, 0.018521, 0.018613, 0.018703, 0.018795, 0.018886, 0.018975, 0.019066, 0.019157, 0.019250, 0.019343, 0.019434, 0.019525, 0.019618, 0.019711, 0.019803, 0.019893, 0.019987, 0.020080, 0.020173, 0.020265, 0.020358, 0.020450, 0.020544, 0.020638, 0.020733, 0.020827, 0.020922, 0.021017, 0.021110, 0.021202, 0.021298, 0.021392, 0.021487, 0.021580, 0.021676, 0.021772, 0.021868, 0.021962, 0.022059, 0.022157, 0.022253, 0.022349, 0.022445, 0.022541, 0.022638, 0.022735, 0.022832, 0.022930, 0.023027, 0.023124, 0.023223, 0.023321, 0.023418, 0.023516, 0.023615, 0.023714, 0.023813, 0.023912, 0.024012, 0.024111, 0.024210, 0.024309, 0.024408, 0.024509, 0.024609, 0.024709, 0.024826, 0.024927, 0.025028, 0.025130, 0.025231, 0.025331, 0.025435, 0.025537, 0.025640, 0.025743, 0.025846, 0.025948, 0.026052, 0.026155, 0.026258, 0.026362, 0.026466, 0.026571, 0.026675, 0.026780, 0.026884, 0.026989, 0.027095, 0.027200, 0.027304, 0.027411, 0.027516, 0.027622, 0.027729, 0.027853, 0.027960, 0.028067, 0.028174, 0.028282, 0.028390, 0.028498, 0.028608, 0.028716, 0.028824, 0.028934, 0.029043, 0.029153, 0.029264, 0.029392, 0.029502, 0.029613, 0.029724, 0.029835, 0.029945, 0.030056, 0.030168, 0.030278, 0.030390, 0.030523, 0.030636, 0.030748, 0.030861, 0.030974, 0.031087, 0.031202, 0.031317, 0.031454, 0.031569, 0.031684, 0.031800, 0.031916, 0.032031, 0.032169, 0.032285, 0.032401, 0.032518};
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

    // AtRGB sigma_s = sigma_s_prime/(1.0f-g);
    AtRGB sigma_s = sigma_s_prime;
    sigma_s_prime /= (1.0f - g);

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
    if (0)//(dir)
    {
        albedo.r = _albedo_lut_d[r_idx];
        albedo.g = _albedo_lut_d[g_idx];
        albedo.b = _albedo_lut_d[b_idx];
    }
    else
    {
        albedo.r = _albedo_lut[r_idx];
        albedo.g = _albedo_lut[g_idx];
        albedo.b = _albedo_lut[b_idx];
    }
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