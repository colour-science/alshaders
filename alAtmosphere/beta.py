import math

def beta(n, N, l):
   return (8 * math.pi**3 * (n**2 - 1)**2) / (3 * N * l**4)

l_r = 680e-9
l_g = 550e-9
l_b = 440e-9
n = 1.000276
N = 2.504e25

# print beta(n, N, l_r)
# print beta(n, N, l_g)
# print beta(n, N, l_b)

def planck(temp, l):
    p = (3.74183e-16 * math.pow(l, -5.0)) / (math.exp(1.4388e-2 / (l * temp)) - 1.0)
    # return p * math.pow(temp, 4) * 5.67e-8
    return p;

for i in range(300, 780, 10):
    print planck(6500, i*1e-9)
