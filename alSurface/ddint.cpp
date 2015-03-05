#include "alUtil.h"
#include "sss.h"

int main()
{
#if 0
   AtRGB sigma_s_prime = rgb(9.820210, 23.136164, 48.246452);
   AtRGB sigma_a = rgb(0.416882, 8.468369, 45.620243);

   ScatteringParamsDirectional sp(sigma_s_prime, sigma_a, 0.0f);

   float l = maxh(sp.zr);
   float rmax = l * 25.0f;

   AtRGB result = integrateDirectional(sp, rmax, 1000);

   printf("result: (%f, %f, %f)\n", result.r, result.g, result.b);
#endif
   /*
   int numsteps = 100;
   float step = 1.0f / float(numsteps); 
   int count = 0;
   printf("#include \"norm_table.h\"\nfloat norm_table[NORM_SZ*NORM_SZ*NORM_SZ*3] = {");
   for (float b = step*0.5f; b < 1.0f; b += step)
   {
      for (float g = step*0.5f; g < 1.0f; g += step)
      {
         for (float r = step*0.5f; r < 1.0f; r += step)
         {
            AtRGB sc = rgb(r, g, b);
            ScatteringParamsDirectional sp(sc, 10.0f, 0.0f);
            float l = maxh(sp.zr);
            float rmax = l * SSS_MAX_RADIUS;
            AtRGB result = integrateDirectional(sp, rmax, 1000);
            if (count == 0)
            {

            }
            else
            {
               printf(", ");
            }
            printf("%f, %f, %f", result.r, result.g, result.b);
            

            count++;
         }
         // printf("\n");
      }
      // printf("\n");
   }
   printf("};\n");
   */

   float step = 1.0f / SSS_ALBEDO_LUT_SZ;
   int count = 0;
   printf("_albedo_lut[SSS_ALBEDO_LUT_SZ] = {");

   for (float f=step/2; f < 1.0f; f += step)
   {
      AtRGB sc = rgb(f, f, f);
      ScatteringParamsDirectional sp(sc, 1.0f, 0.0f, false, false);
      float l = maxh(sp.zr);
      float rmax = l * SSS_MAX_RADIUS;
      AtRGB result = integrateDirectional(sp, rmax, 10000);
      if (count == 0)
      {

      }
      else
      {
         printf(", ");
      }
      printf("%f", result.r);
      count++;
   }

   printf("};\n");

   printf("_albedo_lut_d[SSS_ALBEDO_LUT_SZ] = {");
   count = 0;
   for (float f=step/2; f < 1.0f; f += step)
   {
      AtRGB sc = rgb(f, f, f);
      ScatteringParamsDirectional sp(sc, 1.0f, 0.0f, false, true);
      float l = maxh(sp.zr);
      float rmax = l * SSS_MAX_RADIUS;
      AtRGB result = integrateDirectionalHemi(sp, rmax, 1000);
      if (count == 0)
      {

      }
      else
      {
         printf(", ");
      }
      printf("%f", result.r);
      count++;
   }

   printf("};\n");

   return 0;

}