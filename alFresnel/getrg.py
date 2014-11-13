import sys
from math import *

if len(sys.argv) != 7:
   print 'usage:\n\t %s n_r k_r n_g k_g n_b k_b\n\n\tInputs are eta and absorption at 650nm, 550nm and 450nm from refractiveindex.info' % sys.argv[0]
   sys.exit(-1)

def n_min(r):
   return (1-r)/(1+r);


def n_max(r):
   return (1.0+sqrt(r))/(1.0-sqrt(r)); 


def get_r(n,  k):
   return ((n-1.0)*(n-1.0)+k*k)/((n+1.0)*(n+1.0)+k*k)


def get_g(n,  k):
   r = get_r(n,k)
   return (n_max(r)-n)/(n_max(r)-n_min(r))

r_r = get_r(float(sys.argv[1]), float(sys.argv[2]))
r_g = get_r(float(sys.argv[3]), float(sys.argv[4]))
r_b = get_r(float(sys.argv[5]), float(sys.argv[6]))

g_r = get_g(float(sys.argv[1]), float(sys.argv[2]))
g_g = get_g(float(sys.argv[3]), float(sys.argv[4]))
g_b = get_g(float(sys.argv[5]), float(sys.argv[6]))

print 'r: %.3f %.3f %.3f' % (r_r, r_g, r_b)
print 'g: %.3f %.3f %.3f' % (g_r, g_g, g_b)