#pragma once

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <OpenEXR/ImathVec.h>

#include <ai.h>

#define IMPORTANCE_EPS 0.01

// concentricSampleDisk and cosineSampleHemisphere lifted from PBRT
/*
Copyright (c) 1998-2012, Matt Pharr and Greg Humphreys.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
* Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
inline void concentricSampleDisk(float u1, float u2, float& dx, float& dy)
{
    float r, theta;
    // Map uniform random numbers to $[-1,1]^2$
    float sx = 2 * u1 - 1;
    float sy = 2 * u2 - 1;

    // Map square to $(r,\theta)$

    // Handle degeneracy at the origin
    if (sx == 0.0 && sy == 0.0) {
        dx = 0.0;
        dy = 0.0;
        return;
    }
    if (sx >= -sy) {
        if (sx > sy) {
            // Handle first region of disk
            r = sx;
            if (sy > 0.0) theta = sy/r;
            else          theta = 8.0f + sy/r;
        }
        else {
            // Handle second region of disk
            r = sy;
            theta = 2.0f - sx/r;
        }
    }
    else {
        if (sx <= sy) {
            // Handle third region of disk
            r = -sx;
            theta = 4.0f - sy/r;
        }
        else {
            // Handle fourth region of disk
            r = -sy;
            theta = 6.0f + sx/r;
        }
    }
    theta *= M_PI / 4.f;
    dx = r * cosf(theta);
    dy = r * sinf(theta);
}

inline Imath::V3f cosineSampleHemisphere(float u1, float u2)
{
   Imath::V3f ret;
   concentricSampleDisk(u1, u2, ret.x, ret.z);
   ret.y = sqrtf(std::max(0.0f, 1.0f - ret.x*ret.x - ret.z*ret.z));
   return ret;
}

inline float maxh(const AtRGB& c)
{
   return std::max(std::max(c.r, c.g), c.b);
}

inline float minh(const AtRGB& c)
{
   return std::min(std::min(c.r, c.g ), c.b);
}

inline AtRGB max(const AtRGB& c1, const AtRGB& c2)
{
	AtRGB c;
	c.r = std::max(c1.r, c2.r);
	c.g = std::max(c1.g, c2.g);
	c.b = std::max(c1.b, c2.b);
	return c;
}

inline AtRGB min(const AtRGB& c1, const AtRGB& c2)
{
	AtRGB c;
	c.r = std::min(c1.r, c2.r);
	c.g = std::min(c1.g, c2.g);
	c.b = std::min(c1.b, c2.b);
	return c;
}

inline float clamp(float a, float mn, float mx)
{
	return std::min(std::max(a, mn), mx);
}

inline AtRGB clamp(const AtRGB& a, const AtRGB& mn, const AtRGB& mx)
{
	return min(max(a, mn), mx);
}

inline float lerp(const float a, const float b, const float t)
{
   return (1-t)*a + t*b;
}

inline AtRGB lerp(const AtRGB& a, const AtRGB& b, const float t)
{
   AtRGB r;
   r.r = lerp( a.r, b.r, t );
   r.g = lerp( a.g, b.g, t );
   r.b = lerp( a.b, b.b, t );
   return r;
}

inline AtFloat fresnel(AtFloat cosi, AtFloat etai)
{
	if (cosi >= 1.0f) return 0.0f;
	AtFloat sint = etai * sqrtf(1.0f-cosi*cosi);
	if ( sint >= 1.0f ) return 1.0f;

	AtFloat cost = sqrtf(1.0f-sint*sint);
	AtFloat pl =	(cosi - (etai * cost))
					/ (cosi + (etai * cost));
	AtFloat pp =	((etai * cosi) - cost)
					/ ((etai * cosi) + cost);
	return (pl*pl+pp*pp)*0.5f;
}

// Stolen wholesale from OSL:
// http://code.google.com/p/openshadinglanguage
/*
Copyright (c) 2009-2010 Sony Pictures Imageworks Inc., et al.
All Rights Reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
* Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.
* Neither the name of Sony Pictures Imageworks nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

inline AtFloat fresnel(AtFloat eta, const AtVector& N, const AtVector& I, AtVector& R, AtVector& T, bool& is_inside)
{
	AtFloat c = AiV3Dot(N,I);
	AtFloat neta;
    AtVector Nn;
    // compute reflection
    R = (2 * c) * N - I;
    // check which side of the surface we are on
    if (c > 0) {
        // we are on the outside of the surface, going in
        neta = 1 / eta;
        Nn   = N;
        is_inside = false;
    } else {
        // we are inside the surface,
        c  = -c;
        neta = eta;
        Nn   = -N;
        is_inside = true;
    }
    R = (2 * c) * Nn - I;
    float arg = 1 - (neta * neta * (1 - (c * c)));
    if (arg < 0) {
        T.x = T.y = T.z;
        return 1; // total internal reflection
    } else {
        float dnp = sqrtf(arg);
        float nK = (neta * c) - dnp;
        T = -(neta * I) + (nK * Nn);
        // compute Fresnel terms
        float cosTheta1 = c; // N.R
        float cosTheta2 = -AiV3Dot(Nn,T);
        float pPara = (cosTheta1 - eta * cosTheta2) / (cosTheta1 + eta * cosTheta2);
        float pPerp = (eta * cosTheta1 - cosTheta2) / (eta * cosTheta1 + cosTheta2);
        return 0.5f * (pPara * pPara + pPerp * pPerp);
    }
}

inline AtRGB sqrt(AtRGB c)
{
	c.r = sqrtf(c.r);
	c.g = sqrtf(c.g);
	c.b = sqrtf(c.b);
	return c;
}

inline AtRGB exp(AtRGB c)
{
	c.r = expf(c.r);
	c.g = expf(c.g);
	c.b = expf(c.b);
	return c;
}

inline AtRGB pow(AtRGB c, float e)
{
	c.r = powf(c.r, e);
	c.g = powf(c.g, e);
	c.b = powf(c.b, e);
	return c;
}

inline AtRGB rgb(float f)
{
	AtRGB c;
	c.r = c.g = c.b = f;
	return c;
}

inline AtRGB rgb(float r, float g, float b)
{
	AtRGB c;
	c.r = r; c.g = g; c.b = b;
	return c;
}

inline float luminance(const AtRGB& c)
{
	return c.r*0.212671 + c.g*0.715160 + c.b*0.072169;
}

inline float luminance(float f)
{
	return f;
}

template <typename T>
inline T contrast(T input, float contrast, float pivot, float softClip)
{
	if (contrast == 1.0f) return input;

	float xmin = std::max(0.0f, (pivot * (contrast-1.0f)/contrast));
	float xmax = std::min(1.0f, ((1.0f-pivot * (1.0f-contrast))/contrast));

	float xl = (pivot - xmin) * softClip + xmin;
	float al = contrast * xl + pivot * (1.0f-contrast);
	float gl = contrast * xl / al;
	float xh = xmax - (xmax-pivot) * softClip;
	float ah = 1.0f - (contrast*xh + pivot * (1.0f-contrast));
	float gh = contrast * (1.0f-xh) / ah;

	float lum = luminance(input);
	if (softClip > 0.0f)
	{
		lum = clamp(lum, 0.0f, 1.0f);
	}

	float result;
	if (softClip != 0 && lum < xl)
	{
		result = (al * pow(lum/xl, gl)) / lum;
	}
	else if (softClip != 0.0f && lum > xh)
	{
		result = (1.0f - ah * pow((1.0f-lum)/(1-xh), gh)) / lum;
	}
	else
	{
		result = (contrast*lum + pivot * (1.0f-contrast)) / lum;
	}
	result = std::max(0.0f, result);

	return result * input;
}

inline AtRGB rgb2hsv(AtRGB rgb)
{
	float mn = minh(rgb);
	float mx = maxh(rgb);
	float chroma = mx - mn;
	float h, s, v;
	v = mx;
	if (v > 0.0f)
	{
		s = chroma / mx;
	}
	else
	{
		s = 0.0f;
	}

	if (s <= 0.0f)
	{
		h = 0.0f;
	}
	else
	{
		if (rgb.r >= mx)
		{
			h = (rgb.g - rgb.b) / chroma;
		}
		else if (rgb.g >= mx)
		{
			h = 2.0f + (rgb.b - rgb.r) / chroma;
		}
		else
		{
			h = 4.0f + (rgb.r - rgb.g) / chroma;
		}
		h /= 6.0f;
		if (h < 0.0f) h += 1.0f;
	}
	return AiColorCreate(h,s,v);
}

inline AtRGB hsv2rgb(AtRGB hsv)
{
	if (hsv.g == 0.0f)
	{
		return AiColorCreate(hsv.b, hsv.b, hsv.b);
	}
	else
	{
		float h = 6.0f * (hsv.r - floorf(hsv.r));
		const float& s = hsv.g;
		const float& v = hsv.b;
		int hi = int(hi);
		float f = h - hi;
		float p = v * (1.0f - s);
		float q = v * (1.0f - s*f);
		float t = v * (1.0f - s*(1.0f-f));

		switch(hi)
		{
		case 0: return AiColorCreate(v, t, p);
		case 1: return AiColorCreate(q, v, p);
		case 2: return AiColorCreate(p, v, t);
		case 3: return AiColorCreate(p, q, v);
		case 4: return AiColorCreate(t, p, v);
		default: return AiColorCreate(v, p, q);
		}
	}

}

// For the sake of simplicity we limit eta to 1.3 so cache A here
inline float A(float eta)
{
	float Fdr = -1.440/(eta*eta) + 0.710/eta + 0.668 + 0.0636*eta;
	return (1 + Fdr)/(1 - Fdr);
}

inline AtRGB brdf( const AtRGB& _alpha_prime )
{
	AtRGB sq = sqrt( 3.0f * (AI_RGB_WHITE - _alpha_prime) );
	return _alpha_prime * 0.5f * (AI_RGB_WHITE + exp( -(A(1.3f)*rgb(4.0f/3.0f)*sq ) )) * exp( -sq );
}

void alphaInversion( const AtRGB& scatterColour, float scatterDist, AtRGB& sigma_s_prime_, AtRGB& sigma_a_ );
float alpha1_3( float x );
inline AtRGB alpha1_3(const AtRGB& c)
{
	return rgb(alpha1_3(c.r), alpha1_3(c.g), alpha1_3(c.b));
}

inline std::ostream& operator<<( std::ostream& os, AtVector v )
{
   os << "(" << v.x << "," << v.y << "," << v.z << ")";
   return os;
}

