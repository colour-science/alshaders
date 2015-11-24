import math

def beta(n, N, l):
   return (8 * math.pi**3 * (n**2 - 1)**2) / (3 * N * l**4) 
   
l_r = 680e-9
l_g = 550e-9
l_b = 440e-9
n = 1.000276
N = 2.504e25

print beta(n, N, l_r)
print beta(n, N, l_g)
print beta(n, N, l_b)

