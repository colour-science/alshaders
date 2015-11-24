#ifndef SLY_MATHUTIL_H
#define SLY_MATHUTIL_H

#include "constants.h"
#include <OpenEXR/ImathVec.h>
#include <OpenEXR/ImathMatrix.h>

namespace sly
{

using Imath::V2f;
using Imath::V3f;
using Imath::M33f;
using Imath::M44f;
using std::min;
using std::max;

inline bool solve_linear_system_22(const float A[2][2], const float B[2],
                                   float& x0, float& x1)
{
    float det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
    if (fabsf(det) < 1e-10f) return false;

    x0 = (A[1][1] * B[0] - A[0][1] * B[1]) / det;
    x1 = (A[0][0] * B[1] - A[1][0] * B[0]) / det;

    if (!std::isfinite(x0) || !std::isfinite(x1)) return false;

    return true;
}

inline float radians(float d) { return d * c_pi_over_180; }

inline M44f translate(float x, float y, float z)
{
    M44f xf;
    xf.translate(V3f(x, y, z));
    return xf;
}

inline M44f scale(float x, float y, float z)
{
    M44f xf;
    xf.scale(V3f(x, y, z));
    return xf;
}

inline float lerp(float a, float b, float t) { return (1.0f - t) * a + t * b; }

inline float clamp(const float f, const float mn = 0.0f,
                   const float mx = c_infinity)
{
    return max(min(f, mx), mn);
}

inline bool ispfinite(float f) { return std::isfinite(f) && (f >= 0); }

inline float hmax(V3f v) { return std::max(std::max(v.x, v.y), v.z); }

inline V3f normalize(V3f v) { return v.normalize(); }

inline float dot(V3f a, V3f b) { return a.dot(b); }

inline float absdot(V3f a, V3f b) { return fabsf(a.dot(b)); }

inline float clampdot(V3f a, V3f b) { return clamp(a.dot(b)); }

inline float safe_acosf(float cos_theta)
{
    return acosf(clamp(cos_theta, -1.0f, 1.0f));
}

template <typename T>
inline T sqr(T x)
{
    return x * x;
}

inline float powerHeuristic(float fp, float gp)
{
    return sqr(fp) / (sqr(fp) + sqr(gp));
}

inline V3f cross(V3f a, V3f b) { return a.cross(b); }

inline V3f multmtxpnt(M44f mtx, V3f p)
{
    V3f r;
    mtx.multVecMatrix(p, r);
    return r;
}

inline V3f multmtxvec(M44f mtx, V3f p)
{
    V3f r;
    mtx.multDirMatrix(p, r);
    return r;
}

inline V3f reflect(const V3f wo, const V3f n)
{
    return -wo + 2.0f * dot(wo, n) * n;
}

inline bool refract(const V3f wi, const V3f n, float eta, V3f& wt)
{
    // Compute $\cos \theta_\roman{t}$ using Snell's law
    float cos_theta_i = dot(n, wi);
    float sin2_theta_i = std::max(0.f, 1.f - cos_theta_i * cos_theta_i);
    float sin2_theta_t = eta * eta * sin2_theta_i;

    // Handle total internal reflection for transmission
    if (sin2_theta_t >= 1.f) return false;

    float cos_theta_t = sqrtf(1.f - sin2_theta_t);
    wt = eta * -wi + (eta * cos_theta_i - cos_theta_t) * n;
    return true;
}

// Note: Could replace this with memcpy, which gcc optimizes to the same
// assembly as the code below. I'm not sure how other compiler treat it though,
// since it's really part of the C runtime. The union seems to be portable
// enough.
static inline float uintBitsToFloat(uint32_t i)
{
    union
    {
        float f;
        uint32_t i;
    } unionHack;
    unionHack.i = i;
    return unionHack.f;
}

static inline uint32_t floatBitsToUint(float f)
{
    union
    {
        float f;
        uint32_t i;
    } unionHack;
    unionHack.f = f;
    return unionHack.i;
}

// 2x-5x faster than i/float(UINT_MAX)
static inline float normalizedUint(uint32_t i)
{
    return uintBitsToFloat((i >> 9u) | 0x3F800000u) - 1.0f;
}

////////////////////////////////////////////////////////////////////////////
// Everything below this is ripped wholesale from OIIO

/// Multiply and add: (a * b) + c
template <typename T>
inline T madd(const T& a, const T& b, const T& c)
{
    // NOTE:  in the future we may want to explicitly ask for a fused
    // multiply-add in a specialized version for float.
    // NOTE2: GCC/ICC will turn this (for float) into a FMA unless
    // explicitly asked not to, clang seems to leave the code alone.
    return a * b + c;
}

template <typename IN_TYPE, typename OUT_TYPE>
inline OUT_TYPE bit_cast(const IN_TYPE in)
{
    // NOTE: this is the only standards compliant way of doing this type of
    // casting,
    // luckily the compilers we care about know how to optimize away this idiom.
    OUT_TYPE out;
    memcpy(&out, &in, sizeof(IN_TYPE));
    return out;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//
// SAFE MATH
//
// The functions named "safe_*" are versions with various internal clamps
// or other deviations from IEEE standards with the specific intent of
// never producing NaN or Inf values or throwing exceptions. But within the
// valid range, they should be full precision and match IEEE standards.
//

/// Safe (clamping) sqrt: safe_sqrt(x<0) returns 0, not NaN.
template <typename T>
inline T safe_sqrt(T x)
{
    return x >= T(0) ? std::sqrt(x) : T(0);
}

/// Safe (clamping) inverse sqrt: safe_inversesqrt(x<=0) returns 0.
template <typename T>
inline T safe_inversesqrt(T x)
{
    return x > T(0) ? T(1) / std::sqrt(x) : T(0);
}

/// Safe (clamping) arcsine: clamp to the valid domain.
template <typename T>
inline T safe_asin(T x)
{
    if (x <= T(-1)) return T(-M_PI_2);
    if (x >= T(+1)) return T(+M_PI_2);
    return std::asin(x);
}

/// Safe (clamping) arccosine: clamp to the valid domain.
template <typename T>
inline T safe_acos(T x)
{
    if (x <= T(-1)) return T(M_PI);
    if (x >= T(+1)) return T(0);
    return std::acos(x);
}

/// Safe log2: clamp to valid domain.
template <typename T>
inline T safe_log2(T x)
{
    // match clamping from fast version
    if (x < std::numeric_limits<T>::min()) x = std::numeric_limits<T>::min();
    if (x > std::numeric_limits<T>::max()) x = std::numeric_limits<T>::max();
#if defined(OIIO_CPLUSPLUS11)
    return std::log2(x);
#else
    return log2f(x);  // punt: just use the float one
#endif
}

/// Safe log: clamp to valid domain.
template <typename T>
inline T safe_log(T x)
{
    // slightly different than fast version since clamping happens before
    // scaling
    if (x < std::numeric_limits<T>::min()) x = std::numeric_limits<T>::min();
    if (x > std::numeric_limits<T>::max()) x = std::numeric_limits<T>::max();
    return std::log(x);
}

/// Safe log10: clamp to valid domain.
template <typename T>
inline T safe_log10(T x)
{
    // slightly different than fast version since clamping happens before
    // scaling
    if (x < std::numeric_limits<T>::min()) x = std::numeric_limits<T>::min();
    if (x > std::numeric_limits<T>::max()) x = std::numeric_limits<T>::max();
    return log10f(x);
}

/// Safe logb: clamp to valid domain.
template <typename T>
inline T safe_logb(T x)
{
#if defined(OIIO_CPLUSPLUS11)
    return (x != T(0)) ? std::logb(x) : -std::numeric_limits<T>::max();
#else
    return (x != T(0)) ? logbf(x) : -std::numeric_limits<T>::max();
#endif
}

/// Safe pow: clamp the domain so it never returns Inf or NaN or has divide
/// by zero error.
template <typename T>
inline T safe_pow(T x, T y)
{
    if (y == T(0)) return T(1);
    if (x == T(0)) return T(0);
    // if x is negative, only deal with integer powers
    if ((x < T(0)) && (y != floor(y))) return T(0);
    // FIXME: this does not match "fast" variant because clamping limits are
    // different
    T r = std::pow(x, y);
    // Clamp to avoid returning Inf.
    const T big = std::numeric_limits<T>::max();
    return clamp(r, -big, big);
}

// (end of safe_* functions)
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//
// FAST & APPROXIMATE MATH
//
// The functions named "fast_*" provide a set of replacements to libm that
// are much faster at the expense of some accuracy and robust handling of
// extreme values. One design goal for these approximation was to avoid
// branches as much as possible and operate on single precision values only
// so that SIMD versions should be straightforward ports We also try to
// implement "safe" semantics (ie: clamp to valid range where possible)
// natively since wrapping these inline calls in another layer would be
// wasteful.
//
// Some functions are fast_safe_*, which is both a faster approximation as
// well as clamped input domain to ensure no NaN, Inf, or divide by zero.
//

/// Round to nearest integer, returning as an int.
inline int fast_rint(float x)
{
// used by sin/cos/tan range reduction
#if OIIO_SIMD_SSE >= 4
    // single roundps instruction on SSE4.1+ (for gcc/clang at least)
    return static_cast<int>(rintf(x));
#else
    // emulate rounding by adding/substracting 0.5
    return static_cast<int>(x + copysignf(0.5f, x));
#endif
}

inline float fast_sin(float x)
{
    // very accurate argument reduction from SLEEF
    // starts failing around x=262000
    // Results on: [-2pi,2pi]
    // Examined 2173837240 values of sin: 0.00662760244 avg ulp diff, 2 max ulp,
    // 1.19209e-07 max error
    int q = fast_rint(x * float(M_1_PI));
    float qf = q;
    x = madd(qf, -0.78515625f * 4, x);
    x = madd(qf, -0.00024187564849853515625f * 4, x);
    x = madd(qf, -3.7747668102383613586e-08f * 4, x);
    x = madd(qf, -1.2816720341285448015e-12f * 4, x);
    x = float(M_PI_2) - (float(M_PI_2) - x);  // crush denormals
    float s = x * x;
    if ((q & 1) != 0) x = -x;
    // this polynomial approximation has very low error on [-pi/2,+pi/2]
    // 1.19209e-07 max error in total over [-2pi,+2pi]
    float u = 2.6083159809786593541503e-06f;
    u = madd(u, s, -0.0001981069071916863322258f);
    u = madd(u, s, +0.00833307858556509017944336f);
    u = madd(u, s, -0.166666597127914428710938f);
    u = madd(s, u * x, x);
    // For large x, the argument reduction can fail and the polynomial can be
    // evaluated with arguments outside the valid internal. Just clamp the bad
    // values away (setting to 0.0f means no branches need to be generated).
    if (fabsf(u) > 1.0f) u = 0.0f;
    return u;
}

inline float fast_cos(float x)
{
    // same argument reduction as fast_sin
    int q = fast_rint(x * float(M_1_PI));
    float qf = q;
    x = madd(qf, -0.78515625f * 4, x);
    x = madd(qf, -0.00024187564849853515625f * 4, x);
    x = madd(qf, -3.7747668102383613586e-08f * 4, x);
    x = madd(qf, -1.2816720341285448015e-12f * 4, x);
    x = float(M_PI_2) - (float(M_PI_2) - x);  // crush denormals
    float s = x * x;
    // polynomial from SLEEF's sincosf, max error is
    // 4.33127e-07 over [-2pi,2pi] (98% of values are "exact")
    float u = -2.71811842367242206819355e-07f;
    u = madd(u, s, +2.47990446951007470488548e-05f);
    u = madd(u, s, -0.00138888787478208541870117f);
    u = madd(u, s, +0.0416666641831398010253906f);
    u = madd(u, s, -0.5f);
    u = madd(u, s, +1.0f);
    if ((q & 1) != 0) u = -u;
    if (fabsf(u) > 1.0f) u = 0.0f;
    return u;
}

inline void fast_sincos(float x, float* sine, float* cosine)
{
    // same argument reduction as fast_sin
    int q = fast_rint(x * float(M_1_PI));
    float qf = q;
    x = madd(qf, -0.78515625f * 4, x);
    x = madd(qf, -0.00024187564849853515625f * 4, x);
    x = madd(qf, -3.7747668102383613586e-08f * 4, x);
    x = madd(qf, -1.2816720341285448015e-12f * 4, x);
    x = float(M_PI_2) - (float(M_PI_2) - x);  // crush denormals
    float s = x * x;
    // NOTE: same exact polynomials as fast_sin and fast_cos above
    if ((q & 1) != 0) x = -x;
    float su = 2.6083159809786593541503e-06f;
    su = madd(su, s, -0.0001981069071916863322258f);
    su = madd(su, s, +0.00833307858556509017944336f);
    su = madd(su, s, -0.166666597127914428710938f);
    su = madd(s, su * x, x);
    float cu = -2.71811842367242206819355e-07f;
    cu = madd(cu, s, +2.47990446951007470488548e-05f);
    cu = madd(cu, s, -0.00138888787478208541870117f);
    cu = madd(cu, s, +0.0416666641831398010253906f);
    cu = madd(cu, s, -0.5f);
    cu = madd(cu, s, +1.0f);
    if ((q & 1) != 0) cu = -cu;
    if (fabsf(su) > 1.0f) su = 0.0f;
    if (fabsf(cu) > 1.0f) cu = 0.0f;
    *sine = su;
    *cosine = cu;
}

// NOTE: this approximation is only valid on [-8192.0,+8192.0], it starts
// becoming
// really poor outside of this range because the reciprocal amplifies errors
inline float fast_tan(float x)
{
    // derived from SLEEF implementation
    // note that we cannot apply the "denormal crush" trick everywhere because
    // we sometimes need to take the reciprocal of the polynomial
    int q = fast_rint(x * float(2 * M_1_PI));
    float qf = q;
    x = madd(qf, -0.78515625f * 2, x);
    x = madd(qf, -0.00024187564849853515625f * 2, x);
    x = madd(qf, -3.7747668102383613586e-08f * 2, x);
    x = madd(qf, -1.2816720341285448015e-12f * 2, x);
    if ((q & 1) == 0)
        x = float(M_PI_4) - (float(M_PI_4) - x);  // crush denormals (only if we
                                                  // aren't inverting the result
                                                  // later)
    float s = x * x;
    float u = 0.00927245803177356719970703f;
    u = madd(u, s, 0.00331984995864331722259521f);
    u = madd(u, s, 0.0242998078465461730957031f);
    u = madd(u, s, 0.0534495301544666290283203f);
    u = madd(u, s, 0.133383005857467651367188f);
    u = madd(u, s, 0.333331853151321411132812f);
    u = madd(s, u * x, x);
    if ((q & 1) != 0) u = -1.0f / u;
    return u;
}

/// Fast, approximate sin(x*M_PI) with maximum absolute error of 0.000918954611.
/// Adapted from
/// http://devmaster.net/posts/9648/fast-and-accurate-sine-cosine#comment-76773
inline float fast_sinpi(float x)
{
    // Fast trick to strip the integral part off, so our domain is [-1,1]
    const float z = x - ((x + 25165824.0f) - 25165824.0f);
    const float y = z - z * fabsf(z);
    const float Q = 3.10396624f;
    const float P = 3.584135056f;  // P = 16-4*Q
    return y * (Q + P * fabsf(y));
    /* The original article used used inferior constants for Q and P and
     * so had max error 1.091e-3.
     *
     * The optimal value for Q was determined by exhaustive search, minimizing
     * the absolute numerical error relative to
     *float(std::sin(double(phi*M_PI)))
     * over the interval [0,2] (which is where most of the invocations happen).
     *
     * The basic idea of this approximation starts with the coarse
     *approximation:
     *      sin(pi*x) ~= f(x) =  4 * (x - x * abs(x))
     *
     * This approximation always _over_ estimates the target. On the otherhand,
     *the
     * curve:
     *      sin(pi*x) ~= f(x) * abs(f(x)) / 4
     *
     * always lies _under_ the target. Thus we can simply numerically search for
     *the
     * optimal constant to LERP these curves into a more precise approximation.
     * After folding the constants together and simplifying the resulting math,
     *we
     * end up with the compact implementation below.
     *
     * NOTE: this function actually computes sin(x * pi) which avoids one or two
     * mults in many cases and guarantees exact values at integer periods.
     */
}

/// Fast approximate cos(x*M_PI) with ~0.1% absolute error.
inline float fast_cospi(float x) { return fast_sinpi(x + 0.5f); }

inline float fast_acos(float x)
{
    const float f = fabsf(x);
    const float m =
        (f < 1.0f) ? 1.0f - (1.0f - f) : 1.0f;  // clamp and crush denormals
    // based on http://www.pouet.net/topic.php?which=9132&page=2
    // 85% accurate (ulp 0)
    // Examined 2130706434 values of acos: 15.2000597 avg ulp diff, 4492 max
    // ulp,
    // 4.51803e-05 max error // without "denormal crush"
    // Examined 2130706434 values of acos: 15.2007108 avg ulp diff, 4492 max
    // ulp,
    // 4.51803e-05 max error // with "denormal crush"
    const float a =
        sqrtf(1.0f - m) *
        (1.5707963267f +
         m * (-0.213300989f + m * (0.077980478f + m * -0.02164095f)));
    return x < 0 ? float(M_PI) - a : a;
}

inline float fast_asin(float x)
{
    // based on acosf approximation above
    // max error is 4.51133e-05 (ulps are higher because we are consistently off
    // by a little amount)
    const float f = fabsf(x);
    const float m =
        (f < 1.0f) ? 1.0f - (1.0f - f) : 1.0f;  // clamp and crush denormals
    const float a =
        float(M_PI_2) -
        sqrtf(1.0f - m) *
            (1.5707963267f +
             m * (-0.213300989f + m * (0.077980478f + m * -0.02164095f)));
    return copysignf(a, x);
}

inline float fast_atan(float x)
{
    const float a = fabsf(x);
    const float k = a > 1.0f ? 1 / a : a;
    const float s = 1.0f - (1.0f - k);  // crush denormals
    const float t = s * s;
    // http://mathforum.org/library/drmath/view/62672.html
    // Examined 4278190080 values of atan: 2.36864877 avg ulp diff, 302 max ulp,
    // 6.55651e-06 max error      // (with  denormals)
    // Examined 4278190080 values of atan: 171160502 avg ulp diff, 855638016 max
    // ulp, 6.55651e-06 max error // (crush denormals)
    float r = s * madd(0.43157974f, t, 1.0f) /
              madd(madd(0.05831938f, t, 0.76443945f), t, 1.0f);
    if (a > 1.0f) r = 1.570796326794896557998982f - r;
    return copysignf(r, x);
}

inline float fast_atan2(float y, float x)
{
    // based on atan approximation above
    // the special cases around 0 and infinity were tested explicitly
    // the only case not handled correctly is x=NaN,y=0 which returns 0 instead
    // of nan
    const float a = fabsf(x);
    const float b = fabsf(y);

    const float k =
        (b == 0) ? 0.0f : ((a == b) ? 1.0f : (b > a ? a / b : b / a));
    const float s = 1.0f - (1.0f - k);  // crush denormals
    const float t = s * s;

    float r = s * madd(0.43157974f, t, 1.0f) /
              madd(madd(0.05831938f, t, 0.76443945f), t, 1.0f);

    if (b > a)
        r = 1.570796326794896557998982f - r;  // account for arg reduction
    if (bit_cast<float, unsigned>(x) & 0x80000000u)  // test sign bit of x
        r = float(M_PI) - r;
    return copysignf(r, y);
}

static inline float fast_log2(float x)
{
    // NOTE: clamp to avoid special cases and make result "safe" from large
    // negative values/nans
    if (x < std::numeric_limits<float>::min())
        x = std::numeric_limits<float>::min();
    if (x > std::numeric_limits<float>::max())
        x = std::numeric_limits<float>::max();
    // based on https://github.com/LiraNuna/glsl-sse2/blob/master/source/vec4.h
    unsigned bits = bit_cast<float, unsigned>(x);
    int exponent = int(bits >> 23) - 127;
    float f =
        bit_cast<unsigned, float>((bits & 0x007FFFFF) | 0x3f800000) - 1.0f;
    // Examined 2130706432 values of log2 on [1.17549435e-38,3.40282347e+38]:
    // 0.0797524457 avg ulp diff, 3713596 max ulp, 7.62939e-06 max error
    // ulp histogram:
    //  0  = 97.46%
    //  1  =  2.29%
    //  2  =  0.11%
    float f2 = f * f;
    float f4 = f2 * f2;
    float hi = madd(f, -0.00931049621349f, 0.05206469089414f);
    float lo = madd(f, 0.47868480909345f, -0.72116591947498f);
    hi = madd(f, hi, -0.13753123777116f);
    hi = madd(f, hi, 0.24187369696082f);
    hi = madd(f, hi, -0.34730547155299f);
    lo = madd(f, lo, 1.442689881667200f);
    return ((f4 * hi) + (f * lo)) + exponent;
}

inline float fast_log(float x)
{
    // Examined 2130706432 values of logf on [1.17549435e-38,3.40282347e+38]:
    // 0.313865375 avg ulp diff, 5148137 max ulp, 7.62939e-06 max error
    return fast_log2(x) * float(M_LN2);
}

inline float fast_log10(float x)
{
    // Examined 2130706432 values of log10f on [1.17549435e-38,3.40282347e+38]:
    // 0.631237033 avg ulp diff, 4471615 max ulp, 3.8147e-06 max error
    return fast_log2(x) * float(M_LN2 / M_LN10);
}

inline float fast_logb(float x)
{
    // don't bother with denormals
    x = fabsf(x);
    if (x < std::numeric_limits<float>::min())
        x = std::numeric_limits<float>::min();
    if (x > std::numeric_limits<float>::max())
        x = std::numeric_limits<float>::max();
    unsigned bits = bit_cast<float, unsigned>(x);
    return int(bits >> 23) - 127;
}

inline float fast_exp2(float x)
{
    // clamp to safe range for final addition
    if (x < -126.0f) x = -126.0f;
    if (x > 126.0f) x = 126.0f;
    // range reduction
    int m = int(x);
    x -= m;
    x = 1.0f - (1.0f - x);  // crush denormals (does not affect max ulps!)
    // 5th degree polynomial generated with sollya
    // Examined 2247622658 values of exp2 on [-126,126]: 2.75764912 avg ulp
    // diff,
    // 232 max ulp
    // ulp histogram:
    //  0  = 87.81%
    //  1  =  4.18%
    float r = 1.33336498402e-3f;
    r = madd(x, r, 9.810352697968e-3f);
    r = madd(x, r, 5.551834031939e-2f);
    r = madd(x, r, 0.2401793301105f);
    r = madd(x, r, 0.693144857883f);
    r = madd(x, r, 1.0f);
    // multiply by 2 ^ m by adding in the exponent
    // NOTE: left-shift of negative number is undefined behavior
    return bit_cast<unsigned, float>(bit_cast<float, unsigned>(r) +
                                     (unsigned(m) << 23));
}

inline float fast_exp(float x)
{
    // Examined 2237485550 values of exp on [-87.3300018,87.3300018]: 2.6666452
    // avg ulp diff, 230 max ulp
    return fast_exp2(x * float(1 / M_LN2));
}

/// Faster float exp than is in libm, but still 100% accurate
inline float fast_correct_exp(float x)
{
#if defined(__x86_64__) && defined(__GNU_LIBRARY__) && defined(__GLIBC__) && \
    defined(__GLIBC_MINOR__) && __GLIBC__ <= 2 && __GLIBC_MINOR__ < 16
    // On x86_64, versions of glibc < 2.16 have an issue where expf is
    // much slower than the double version.  This was fixed in glibc 2.16.
    return static_cast<float>(std::exp(static_cast<double>(x)));
#else
    return std::exp(x);
#endif
}

inline float fast_exp10(float x)
{
    // Examined 2217701018 values of exp10 on [-37.9290009,37.9290009]:
    // 2.71732409 avg ulp diff, 232 max ulp
    return fast_exp2(x * float(M_LN10 / M_LN2));
}

inline float fast_expm1(float x)
{
    if (fabsf(x) < 1e-5f)
    {
        x = 1.0f - (1.0f - x);  // crush denormals
        return madd(0.5f, x * x, x);
    }
    else
        return fast_exp(x) - 1.0f;
}

inline float fast_sinh(float x)
{
    float a = fabsf(x);
    if (a > 1.0f)
    {
        // Examined 53389559 values of sinh on [1,87.3300018]: 33.6886442 avg
        // ulp
        // diff, 178 max ulp
        float e = fast_exp(a);
        return copysignf(0.5f * e - 0.5f / e, x);
    }
    else
    {
        a = 1.0f - (1.0f - a);  // crush denorms
        float a2 = a * a;
        // degree 7 polynomial generated with sollya
        // Examined 2130706434 values of sinh on [-1,1]: 1.19209e-07 max error
        float r = 2.03945513931e-4f;
        r = madd(r, a2, 8.32990277558e-3f);
        r = madd(r, a2, 0.1666673421859f);
        r = madd(r * a, a2, a);
        return copysignf(r, x);
    }
}

inline float fast_cosh(float x)
{
    // Examined 2237485550 values of cosh on [-87.3300018,87.3300018]:
    // 1.78256726
    // avg ulp diff, 178 max ulp
    float e = fast_exp(fabsf(x));
    return 0.5f * e + 0.5f / e;
}

inline float fast_tanh(float x)
{
    // Examined 4278190080 values of tanh on [-3.40282347e+38,3.40282347e+38]:
    // 3.12924e-06 max error
    // NOTE: ulp error is high because of sub-optimal handling around the origin
    float e = fast_exp(2.0f * fabsf(x));
    return copysignf(1 - 2 / (1 + e), x);
}

inline float fast_safe_pow(float x, float y)
{
    if (y == 0) return 1.0f;  // x^0=1
    if (x == 0) return 0.0f;  // 0^y=0
    // be cheap & exact for special case of squaring and identity
    if (y == 1.0f) return x;
    if (y == 2.0f) return std::min(x * x, std::numeric_limits<float>::max());
    float sign = 1.0f;
    if (x < 0)
    {
        // if x is negative, only deal with integer powers
        // powf returns NaN for non-integers, we will return 0 instead
        int ybits = bit_cast<float, int>(y) & 0x7fffffff;
        if (ybits >= 0x4b800000)
        {
            // always even int, keep positive
        }
        else if (ybits >= 0x3f800000)
        {
            // bigger than 1, check
            int k = (ybits >> 23) - 127;  // get exponent
            int j = ybits >> (23 - k);    // shift out possible fractional bits
            if ((j << (23 - k)) ==
                ybits)  // rebuild number and check for a match
                sign = bit_cast<int, float>(
                    0x3f800000 | (j << 31));  // +1 for even, -1 for odd
            else
                return 0.0f;  // not integer
        }
        else
        {
            return 0.0f;  // not integer
        }
    }
    return sign * fast_exp2(y * fast_log2(fabsf(x)));
}

inline float fast_erf(float x)
{
    // Examined 1082130433 values of erff on [0,4]: 1.93715e-06 max error
    // Abramowitz and Stegun, 7.1.28
    const float a1 = 0.0705230784f;
    const float a2 = 0.0422820123f;
    const float a3 = 0.0092705272f;
    const float a4 = 0.0001520143f;
    const float a5 = 0.0002765672f;
    const float a6 = 0.0000430638f;
    const float a = fabsf(x);
    const float b = 1.0f - (1.0f - a);  // crush denormals
    const float r = madd(
        madd(madd(madd(madd(madd(a6, b, a5), b, a4), b, a3), b, a2), b, a1), b,
        1.0f);
    const float s = r * r;  // ^2
    const float t = s * s;  // ^4
    const float u = t * t;  // ^8
    const float v = u * u;  // ^16
    return copysignf(1.0f - 1.0f / v, x);
}

inline float fast_erfc(float x)
{
    // Examined 2164260866 values of erfcf on [-4,4]: 1.90735e-06 max error
    // ulp histogram:
    //   0  = 80.30%
    return 1.0f - fast_erf(x);
}

inline float fast_ierf(float x)
{
    // from: Approximating the erfinv function by Mike Giles
    // to avoid trouble at the limit, clamp input to 1-eps
    float a = fabsf(x);
    if (a > 0.99999994f) a = 0.99999994f;
    float w = -fast_log((1.0f - a) * (1.0f + a)), p;
    if (w < 5.0f)
    {
        w = w - 2.5f;
        p = 2.81022636e-08f;
        p = madd(p, w, 3.43273939e-07f);
        p = madd(p, w, -3.5233877e-06f);
        p = madd(p, w, -4.39150654e-06f);
        p = madd(p, w, 0.00021858087f);
        p = madd(p, w, -0.00125372503f);
        p = madd(p, w, -0.00417768164f);
        p = madd(p, w, 0.246640727f);
        p = madd(p, w, 1.50140941f);
    }
    else
    {
        w = sqrtf(w) - 3.0f;
        p = -0.000200214257f;
        p = madd(p, w, 0.000100950558f);
        p = madd(p, w, 0.00134934322f);
        p = madd(p, w, -0.00367342844f);
        p = madd(p, w, 0.00573950773f);
        p = madd(p, w, -0.0076224613f);
        p = madd(p, w, 0.00943887047f);
        p = madd(p, w, 1.00167406f);
        p = madd(p, w, 2.83297682f);
    }
    return p * x;
}

// (end of fast* functions)
////////////////////////////////////////////////////////////////////////////
}

#endif
