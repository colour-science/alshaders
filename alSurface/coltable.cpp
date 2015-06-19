#include "alUtil.h"
#include <cstdio>

int main()
{
   for (int i=1; i <= 100; ++i)
   {
      float ld = float(i) / 100.0f;
      float sigma_s_prime;
      float sigma_a;

      alphaInversion(ld, sigma_s_prime, sigma_a);
      float sigma_t_prime = sigma_s_prime + sigma_a;
      float alpha_prime = sigma_s_prime / sigma_t_prime;
      printf("%f,%f\n", ld, alpha_prime);
   }
}