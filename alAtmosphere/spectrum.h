#ifndef SLY_SPECTRUM_H
#define SLY_SPECTRUM_H

#include <vector>
#include "colorimetry.h"
#include "spectrum_grid.h"
#include "constants.h"

#define HANIKE_CONVERSION

#define slyassert(x, m) assert((x) && m)
namespace sly
{
static const int SpectrumLambdaStart = 390;
static const int SpectrumLambdaEnd = 710;
static const int SpectrumSamples = 32;
static const int SpectrumBinSize =
    (SpectrumLambdaEnd - SpectrumLambdaStart) / SpectrumSamples;

#define SPECTRUM_CIE_X CMF::CIE_1931_2_x
#define SPECTRUM_CIE_Y CMF::CIE_1931_2_y
#define SPECTRUM_CIE_Z CMF::CIE_1931_2_z

enum SpectrumType
{
    kReflectance,
    kIlluminant
};

template <int NumSamples>
class CoefficientSpectrum
{
 public:
    CoefficientSpectrum(float v = 0.0f)
    {
        for (int i = 0; i < NumSamples; ++i)
        {
            _c[i] = v;
        }
    }

    CoefficientSpectrum operator+(const CoefficientSpectrum& s) const
    {
        CoefficientSpectrum r;
        for (int i = 0; i < NumSamples; ++i)
        {
            r._c[i] = _c[i] + s._c[i];
        }

        return r;
    }

    CoefficientSpectrum operator+=(const CoefficientSpectrum& s)
    {
        for (int i = 0; i < NumSamples; ++i)
        {
            _c[i] += s._c[i];
        }

        return *this;
    }

    CoefficientSpectrum operator-(const CoefficientSpectrum& s) const
    {
        CoefficientSpectrum r;
        for (int i = 0; i < NumSamples; ++i)
        {
            r._c[i] = _c[i] - s._c[i];
        }

        return r;
    }

    CoefficientSpectrum operator-=(const CoefficientSpectrum& s)
    {
        for (int i = 0; i < NumSamples; ++i)
        {
            _c[i] -= s._c[i];
        }

        return *this;
    }

    CoefficientSpectrum operator*(const CoefficientSpectrum& s) const
    {
        CoefficientSpectrum r;
        for (int i = 0; i < NumSamples; ++i)
        {
            r._c[i] = _c[i] * s._c[i];
        }

        return r;
    }

    CoefficientSpectrum operator*=(const CoefficientSpectrum& s)
    {
        for (int i = 0; i < NumSamples; ++i)
        {
            _c[i] *= s._c[i];
        }

        return *this;
    }

    CoefficientSpectrum operator/(const CoefficientSpectrum& s) const
    {
        CoefficientSpectrum r;
        for (int i = 0; i < NumSamples; ++i)
        {
            r._c[i] = _c[i] / s._c[i];
        }

        return r;
    }

    CoefficientSpectrum operator/=(const CoefficientSpectrum& s)
    {
        for (int i = 0; i < NumSamples; ++i)
        {
            _c[i] /= s._c[i];
        }

        return *this;
    }

    CoefficientSpectrum operator-() const
    {
        CoefficientSpectrum r;
        for (int i = 0; i < NumSamples; ++i)
        {
            r._c[i] = -_c[i];
        }

        return r;
    }

    bool operator==(const CoefficientSpectrum& s) const
    {
        for (int i = 0; i < NumSamples; ++i)
        {
            if (_c[i] != s._c[i]) return false;
        }
        return true;
    }

    bool operator!=(const CoefficientSpectrum& s) const
    {
        for (int i = 0; i < NumSamples; ++i)
        {
            if (_c[i] != s._c[i]) return true;
        }
        return false;
    }

    float& operator[](size_t i) { return _c[i]; }

    const float& operator[](size_t i) const { return _c[i]; }

    bool isBlack() const
    {
        for (int i = 0; i < NumSamples; ++i)
        {
            if (_c[i] > c_epsilon) return false;
        }

        return true;
    }

    void write(std::ostream& os) const
    {
        float l = SpectrumLambdaStart;
        float l_step =
            (SpectrumLambdaEnd - SpectrumLambdaStart) / SpectrumSamples;
        for (size_t i = 0; i < NumSamples; ++i, l += l_step)
        {
            os << l << " " << _c[i] << std::endl;
        }
    }

    float _c[NumSamples];
};

template <int NumSamples>
CoefficientSpectrum<NumSamples> sqrt(const CoefficientSpectrum<NumSamples>& s)
{
    CoefficientSpectrum<NumSamples> r;
    for (int i = 0; i < NumSamples; ++i)
    {
        r._c[i] = sqrtf(s._c[i]);
    }
    return r;
}

template <int NumSamples>
CoefficientSpectrum<NumSamples> pow(const CoefficientSpectrum<NumSamples>& s,
                                    float f)
{
    CoefficientSpectrum<NumSamples> r;
    for (int i = 0; i < NumSamples; ++i)
    {
        r._c[i] = powf(s._c[i], f);
    }
    return r;
}

template <int NumSamples>
CoefficientSpectrum<NumSamples> exp(const CoefficientSpectrum<NumSamples>& s)
{
    CoefficientSpectrum<NumSamples> r;
    for (int i = 0; i < NumSamples; ++i)
    {
        r._c[i] = expf(s._c[i]);
    }
    return r;
}

template <int NumSamples>
CoefficientSpectrum<NumSamples> log(const CoefficientSpectrum<NumSamples>& s)
{
    CoefficientSpectrum<NumSamples> r;
    for (int i = 0; i < NumSamples; ++i)
    {
        r._c[i] = logf(s._c[i]);
    }
    return r;
}

template <int NumSamples>
CoefficientSpectrum<NumSamples> lerp(const CoefficientSpectrum<NumSamples>& s1,
                                     const CoefficientSpectrum<NumSamples>& s2,
                                     float t)
{
    CoefficientSpectrum<NumSamples> r;
    for (int i = 0; i < NumSamples; ++i)
    {
        r._c[i] = lerp(s1._c[i], s2._c[i], t);
    }
    return r;
}

template <int NumSamples>
CoefficientSpectrum<NumSamples> clamp(const CoefficientSpectrum<NumSamples>& s,
                                      float mn = 0.0f, float mx = INFINITY)
{
    CoefficientSpectrum<NumSamples> r;
    for (int i = 0; i < NumSamples; ++i)
    {
        r._c[i] = clamp(s._c[i], mn, mx);
    }
    return r;
}

template <int NumSamples>
CoefficientSpectrum<NumSamples> min(const CoefficientSpectrum<NumSamples>& s,
                                    float f)
{
    CoefficientSpectrum<NumSamples> r;
    for (int i = 0; i < NumSamples; ++i)
    {
        r._c[i] = min(s._c[i], f);
    }
    return r;
}

template <int NumSamples>
CoefficientSpectrum<NumSamples> max(const CoefficientSpectrum<NumSamples>& s,
                                    float f)
{
    CoefficientSpectrum<NumSamples> r;
    for (int i = 0; i < NumSamples; ++i)
    {
        r._c[i] = max(s._c[i], f);
    }
    return r;
}

template <int NumSamples>
CoefficientSpectrum<NumSamples> operator*(
    float f, const CoefficientSpectrum<NumSamples>& s)
{
    return s * f;
}

template <int NumSamples>
bool isnan(const CoefficientSpectrum<NumSamples>& s)
{
    for (int i = 0; i < NumSamples; ++i)
    {
        if (std::isnan(s._c[i])) return true;
    }
    return false;
}

template <int NumSamples>
bool isfinite(const CoefficientSpectrum<NumSamples>& s)
{
    for (int i = 0; i < NumSamples; ++i)
    {
        if (!std::isfinite(s._c[i])) return false;
    }
    return true;
}

template <int NumSamples>
bool ispositive(const CoefficientSpectrum<NumSamples>& s)
{
    for (int i = 0; i < NumSamples; ++i)
    {
        if (s._c[i] < 0.0f) return false;
    }
    return true;
}

template <int NumSamples>
bool ispfinite(const CoefficientSpectrum<NumSamples>& s)
{
    return isfinite(s) && ispositive(s);
}

template <int NumSamples>
float hmax(const CoefficientSpectrum<NumSamples>& s)
{
    float mx = 0.0f;
    for (int i = 0; i < NumSamples; ++i)
    {
        mx = std::max(mx, s._c[i]);
    }
    return mx;
}

inline bool spectrumSamplesSorted(const float* lambda, int n)
{
    for (int i = 0; i < n - 1; ++i)
    {
        if (lambda[i] > lambda[i + 1]) return false;
    }
    return true;
}

inline void sortSpectrumSamples(float* lambda, float* v, int n)
{
    std::vector<std::pair<float, float> > sortVec;
    sortVec.reserve(n);
    for (int i = 0; i < n; ++i)
    {
        sortVec.push_back(std::make_pair(lambda[i], v[i]));
    }
    std::sort(sortVec.begin(), sortVec.end());
    for (int i = 0; i < n; ++i)
    {
        lambda[i] = sortVec[i].first;
        v[i] = sortVec[i].second;
    }
}

inline float interpolateSpectrumSamples(const float* lambda, const float* v,
                                        int n, float l)
{
    if (l <= lambda[0]) return v[0];
    if (l >= lambda[n - 1]) return v[n - 1];
    for (int i = 0; i < n - 1; ++i)
    {
        if (l >= lambda[i] && l <= lambda[i + 1])
        {
            float t = (l - lambda[i]) / (lambda[i + 1] - lambda[i]);
            return lerp(v[i], v[i + 1], t);
        }
    }

    slyassert(false, "should never get here");
    return (0.0f);
}

struct HeroWavelengths
{
    float lambda[4];

    float& operator[](size_t i) { return lambda[i]; }

    const float& operator[](size_t i) const { return lambda[i]; }
};

class SampledSpectrum : public CoefficientSpectrum<SpectrumSamples>
{
 public:
    SampledSpectrum(float v = 0.0f)
    {
        for (int i = 0; i < SpectrumSamples; ++i)
        {
            _c[i] = v;
        }
    }

    SampledSpectrum(float r, float g, float b, const ColorSpaceRGB& cs,
                    SpectrumType type = kReflectance)
    {
        RGB rgb(r, g, b);
        HeroWavelengths hw;
        *this = fromRGB(cs, rgb, hw, type);
    }

    SampledSpectrum(const std::vector<float>& vec)
    {
        slyassert(vec.size() == SpectrumSamples,
                  "supplied coefficient vector is wrong size");
        for (int i = 0; i < vec.size(); ++i)
        {
            _c[i] = vec[i];
        }
    }

    SampledSpectrum(const SpectralPowerDistribution& spd)
    {
        spd.interpolateOnto(_c, SpectrumLambdaStart, SpectrumLambdaEnd,
                            SpectrumSamples);
    }

    SampledSpectrum(const CoefficientSpectrum<SpectrumSamples>& s)
        : CoefficientSpectrum<SpectrumSamples>(s)
    {
    }

    static float averageSpectrumSamples(const float* lambda, const float* vals,
                                        int n, float lambdaStart,
                                        float lambdaEnd)
    {
        for (int i = 0; i < n - 1; ++i)
            slyassert(lambda[i + 1] > lambda[i], "wavelengths are not ordered");
        slyassert(lambdaStart < lambdaEnd,
                  "start wavelength is not less than end wavelength");
        // Handle cases with out-of-bounds range or single sample only
        if (lambdaEnd <= lambda[0]) return vals[0];
        if (lambdaStart >= lambda[n - 1]) return vals[n - 1];
        if (n == 1) return vals[0];
        float sum = 0.f;
        // Add contributions of constant segments before/after samples
        if (lambdaStart < lambda[0]) sum += vals[0] * (lambda[0] - lambdaStart);
        if (lambdaEnd > lambda[n - 1])
            sum += vals[n - 1] * (lambdaEnd - lambda[n - 1]);

        // Advance to first relevant wavelength segment
        int i = 0;
        while (lambdaStart > lambda[i + 1])
            ++i;
        slyassert(i + 1 < n, "advanced off the end of the SPD");

// Loop over wavelength sample segments and add contributions
#define INTERP(w, i)                                               \
    lerp(((w)-lambda[i]) / (lambda[(i) + 1] - lambda[i]), vals[i], \
         vals[(i) + 1])
#define SEG_AVG(wl0, wl1, i) (0.5f * (INTERP(wl0, i) + INTERP(wl1, i)))
        for (; i + 1 < n && lambdaEnd >= lambda[i]; ++i)
        {
            float segStart = max(lambdaStart, lambda[i]);
            float segEnd = min(lambdaEnd, lambda[i + 1]);
            sum += SEG_AVG(segStart, segEnd, i) * (segEnd - segStart);
        }
#undef INTERP
#undef SEG_AVG
        return sum / (lambdaEnd - lambdaStart);
    }

    static SampledSpectrum fromSampled(const float* lambda, const float* v,
                                       int n)
    {
        // sort samples if unordered
        if (!spectrumSamplesSorted(lambda, n))
        {
            std::vector<float> slambda(&lambda[0], &lambda[n]);
            std::vector<float> sv(&v[0], &v[n]);
            sortSpectrumSamples(&slambda[0], &sv[0], n);
            return fromSampled(&slambda[0], &sv[0], n);
        }

        SampledSpectrum r;
        for (int i = 0; i < SpectrumSamples; ++i)
        {
            // compute average value of SPD over given sample's range
            float l0 = lerp(SpectrumLambdaStart, SpectrumLambdaEnd,
                            float(i) / float(SpectrumSamples));
            float l1 = lerp(SpectrumLambdaStart, SpectrumLambdaEnd,
                            float(i + 1) / float(SpectrumSamples));
            r._c[i] = averageSpectrumSamples(lambda, v, n, l0, l1);
        }
        return r;
    }

    static SampledSpectrum fromSampled(const SampledSpectrum& s,
                                       const ColorSpaceRGB& cs, WhitePoint wp,
                                       const HeroWavelengths& hw)
    {
        return s;
    }

    float interpolate(float l) const
    {
        float l_t = (l - SpectrumLambdaStart) /
                    (SpectrumLambdaEnd - SpectrumLambdaStart);
        int l0 = static_cast<int>(floorf(l_t));
        int l1 = static_cast<int>(ceilf(l_t));
        float t = l_t - l0;

        l0 *= SpectrumSamples;
        l1 *= SpectrumSamples;

        l0 = max(0, l0);
        l1 = min(SpectrumSamples - 1, l1);

        return _c[l0] * (1.0f - t) + _c[l1] * t;
    }

    static const SampledSpectrum& getIlluminant(WhitePoint wp)
    {
        switch (wp)
        {
        case WhitePoint::D65:
            return D65;
        case WhitePoint::D60:
            return D60;
        case WhitePoint::D55:
            return D55;
        case WhitePoint::D50:
            return D50;
        default:
            return ones;
        }
    }

    XYZ toXYZ(const ColorSpaceRGB& cs, const HeroWavelengths& hw,
              WhitePoint wp) const
    {
        XYZ xyz(0.0f);
        const SampledSpectrum& illum = getIlluminant(wp);
        float N = 0.0f;
        for (int i = 0; i < SpectrumSamples; ++i)
        {
            float Me = _c[i] * illum._c[i];
            xyz.x += X._c[i] * Me;
            xyz.y += Y._c[i] * Me;
            xyz.z += Z._c[i] * Me;
            N += Y._c[i] * illum._c[i];
        }
        slyassert(isfinite(xyz.x), "result is not finite");
        slyassert(N > 0.0f, "N is 0");
        return xyz / N;
    }

    RGB toRGB(const ColorSpaceRGB& cs, const HeroWavelengths& hw,
              WhitePoint wp) const
    {
        return cs.XYZtoRGB(toXYZ(cs, hw, wp));
    }

    // Meng et al conversion method from
    // "Physically Meaningful Rendering using Tristimulus Colors"
    static SampledSpectrum fromRGB(const ColorSpaceRGB& cs_, RGB rgb,
                                   const HeroWavelengths&,
                                   SpectrumType type = kReflectance)
    {
        // the conversion method assumes using a white point of E for all
        // colours
        ColorSpaceRGB cs = cs_;
        cs.set_white_point(V2f(0.333f, 0.333f));

        if (type == kReflectance)
        {
            XYZ xyz = cs.RGBtoXYZ(rgb);
            int i = 0;
            SampledSpectrum s;
            for (float lambda = SpectrumLambdaStart; lambda < SpectrumLambdaEnd;
                 lambda += SpectrumBinSize)
            {
                s._c[i] = clamp(spectrum_xyz_to_p(lambda, &xyz.x) /
                                    equal_energy_reflectance,
                                0.0f, 1.0f);
                i++;
            }
            return s;
        }
        else
        {
            if (rgb == RGB(1.0f))
            {
                // if pure white, just return the renderer's white point
                return getIlluminant(cs_.wp);
            }
            else
            {
                XYZ xyz = cs.RGBtoXYZ(rgb);
                int i = 0;
                SampledSpectrum s;
                for (float lambda = SpectrumLambdaStart;
                     lambda < SpectrumLambdaEnd; lambda += SpectrumBinSize)
                {
                    s._c[i] = clamp(spectrum_xyz_to_p(lambda, &xyz.x) /
                                    equal_energy_reflectance);
                    i++;
                }
                //  return getIlluminant(Options::get().white_point) * s;
                return getIlluminant(cs_.wp) * s;
            }
        }
    }

    static SampledSpectrum fromXYZ(XYZ xyz, const ColorSpaceRGB& cs,
                                   const HeroWavelengths& hw,
                                   SpectrumType type = kReflectance)
    {
        RGB rgb = cs.XYZtoRGB(xyz);
        return fromRGB(cs, rgb, hw, type);
    }

    float y() const
    {
        float yy = 0.0f;
        for (int i = 0; i < SpectrumSamples; ++i)
        {
            yy += Y._c[i] * _c[i];
        }
        return yy / yint;
    }

    // this must be called once before SampledSpectrum can be used
    static void init()
    {
        // fill out our CMFs
        SPECTRUM_CIE_X.interpolateOnto(X._c, SpectrumLambdaStart,
                                       SpectrumLambdaEnd, SpectrumSamples);
        SPECTRUM_CIE_Y.interpolateOnto(Y._c, SpectrumLambdaStart,
                                       SpectrumLambdaEnd, SpectrumSamples);
        SPECTRUM_CIE_Z.interpolateOnto(Z._c, SpectrumLambdaStart,
                                       SpectrumLambdaEnd, SpectrumSamples);

        CMF::SmitsWhite.interpolateOnto(rgbRefl2SpectWhite._c,
                                        SpectrumLambdaStart, SpectrumLambdaEnd,
                                        SpectrumSamples);
        CMF::SmitsRed.interpolateOnto(rgbRefl2SpectRed._c, SpectrumLambdaStart,
                                      SpectrumLambdaEnd, SpectrumSamples);
        CMF::SmitsGreen.interpolateOnto(rgbRefl2SpectGreen._c,
                                        SpectrumLambdaStart, SpectrumLambdaEnd,
                                        SpectrumSamples);
        CMF::SmitsBlue.interpolateOnto(rgbRefl2SpectBlue._c,
                                       SpectrumLambdaStart, SpectrumLambdaEnd,
                                       SpectrumSamples);
        CMF::SmitsCyan.interpolateOnto(rgbRefl2SpectCyan._c,
                                       SpectrumLambdaStart, SpectrumLambdaEnd,
                                       SpectrumSamples);
        CMF::SmitsMagenta.interpolateOnto(rgbRefl2SpectMagenta._c,
                                          SpectrumLambdaStart,
                                          SpectrumLambdaEnd, SpectrumSamples);
        CMF::SmitsYellow.interpolateOnto(rgbRefl2SpectYellow._c,
                                         SpectrumLambdaStart, SpectrumLambdaEnd,
                                         SpectrumSamples);

        Illuminant::D50.interpolateOnto(D50._c, SpectrumLambdaStart,
                                        SpectrumLambdaEnd, SpectrumSamples);
        Illuminant::D55.interpolateOnto(D55._c, SpectrumLambdaStart,
                                        SpectrumLambdaEnd, SpectrumSamples);
        Illuminant::D60.interpolateOnto(D60._c, SpectrumLambdaStart,
                                        SpectrumLambdaEnd, SpectrumSamples);
        Illuminant::D65.interpolateOnto(D65._c, SpectrumLambdaStart,
                                        SpectrumLambdaEnd, SpectrumSamples);

        for (int i = 0; i < SpectrumSamples; ++i)
        {
            yint += Y._c[i];
        }
    }

    friend std::ostream& operator<<(std::ostream& os, const SampledSpectrum& s);

    static SampledSpectrum D50;
    static SampledSpectrum D55;
    static SampledSpectrum D60;
    static SampledSpectrum D65;
    static SampledSpectrum ones;
    static SampledSpectrum X, Y, Z;

 private:
    static SampledSpectrum rgbRefl2SpectWhite;
    static SampledSpectrum rgbRefl2SpectCyan;
    static SampledSpectrum rgbRefl2SpectMagenta;
    static SampledSpectrum rgbRefl2SpectYellow;
    static SampledSpectrum rgbRefl2SpectRed;
    static SampledSpectrum rgbRefl2SpectGreen;
    static SampledSpectrum rgbRefl2SpectBlue;
    static SampledSpectrum rgbIllum2SpectWhite;
    static SampledSpectrum rgbIllum2SpectCyan;
    static SampledSpectrum rgbIllum2SpectMagenta;
    static SampledSpectrum rgbIllum2SpectYellow;
    static SampledSpectrum rgbIllum2SpectRed;
    static SampledSpectrum rgbIllum2SpectGreen;
    static SampledSpectrum rgbIllum2SpectBlue;
    static float yint;
};

inline std::ostream& operator<<(std::ostream& os, const SampledSpectrum& s)
{
    os << "{";
    for (size_t i = 0; i < SpectrumSamples; ++i)
    {
        os << s._c[i];
        if (i != SpectrumSamples - 1)
        {
            os << ", ";
        }
    }
    os << "}";
    return os;
}

class RGBSpectrum : public CoefficientSpectrum<3>
{
 public:
    RGBSpectrum(float v = 0.0f) : CoefficientSpectrum<3>(v) {}
    RGBSpectrum(const CoefficientSpectrum<3>& s) : CoefficientSpectrum<3>(s) {}
    RGBSpectrum(const ColorSpaceRGB& cs, float r, float g, float b,
                SpectrumType type = kReflectance)
    {
        RGB rgb(r, g, b);
        HeroWavelengths hw;
        *this = fromRGB(cs, rgb, hw, type);
    }

    static RGBSpectrum fromRGB(const ColorSpaceRGB& cs, RGB rgb,
                               const HeroWavelengths&,
                               SpectrumType type = kReflectance)
    {
        RGBSpectrum s;
        s._c[0] = rgb.r;
        s._c[1] = rgb.g;
        s._c[2] = rgb.b;
        return s;
    }

    static RGBSpectrum fromXYZ(XYZ xyz, const ColorSpaceRGB& cs,
                               const HeroWavelengths& hw)
    {
        return fromRGB(cs, cs.XYZtoRGB(xyz), hw);
    }

    static RGBSpectrum fromSampled(const SampledSpectrum& ss,
                                   const ColorSpaceRGB& cs, WhitePoint wp,
                                   const HeroWavelengths hw)
    {
        XYZ xyz(0.0f);
        const SampledSpectrum& illum = SampledSpectrum::getIlluminant(wp);
        const SampledSpectrum& X = SampledSpectrum::X;
        const SampledSpectrum& Y = SampledSpectrum::Y;
        const SampledSpectrum& Z = SampledSpectrum::Z;
        float N = 0.0f;
        for (int i = 0; i < SpectrumSamples; ++i)
        {
            float Me = ss[i] * illum[i];
            xyz.x += X[i] * Me;
            xyz.y += Y[i] * Me;
            xyz.z += Z[i] * Me;
            N += Y[i] * illum[i];
        }
        slyassert(isfinite(xyz.x), "result is not finite");
        slyassert(N > 0.0f, "N is 0");

        return fromXYZ(xyz / N, cs, hw);
    }

    RGB toRGB(const ColorSpaceRGB& cs, const HeroWavelengths& hw, WhitePoint wp,
              const SampledSpectrum* illum = nullptr) const
    {
        RGB rgb;
        rgb.r = _c[0];
        rgb.g = _c[1];
        rgb.b = _c[2];
        return rgb;
    }

    XYZ toXYZ(const ColorSpaceRGB& cs, const HeroWavelengths& hw,
              WhitePoint wp) const
    {
        return cs.RGBtoXYZ(toRGB(cs, hw, wp));
    }

    friend std::ostream& operator<<(std::ostream& os, const RGBSpectrum& s);

 private:
};

inline std::ostream& operator<<(std::ostream& os, const RGBSpectrum& s)
{
    os << "(" << s._c[0] << ", " << s._c[1] << ", " << s._c[2] << ")";
    return os;
}

// fill out the lambda array to spread the remaining samples over
// the full spectral range. the hero wavelength goes in lambda[0]
inline HeroWavelengths rotate_wavelengths(float lambda_h)
{
    HeroWavelengths hw;
    hw.lambda[0] = lambda_h;
    float lambda_min = SpectrumLambdaStart;
    float lambda_max = SpectrumLambdaEnd;
    float lambda_bar = lambda_max - lambda_min;
    for (int j = 1; j < 4; ++j)
    {
        hw.lambda[j] =
            fmodf(lambda_h - lambda_min + float(j) * 0.25f * lambda_bar,
                  lambda_bar) +
            lambda_min;
    }

    return hw;
}

inline HeroWavelengths sample_hero_wavelength(float u)
{
    float lambda_min = SpectrumLambdaStart;
    float lambda_max = SpectrumLambdaEnd;
    float lambda_bar = lambda_max - lambda_min;
    float lambda_h = lambda_min + u * lambda_bar;

    HeroWavelengths hw;
    hw[0] = lambda_h;

    for (int j = 1; j < 4; ++j)
    {
        hw[j] = fmodf(lambda_h - lambda_min + float(j) * 0.25f * lambda_bar,
                      lambda_bar) +
                lambda_min;
    }

    return hw;
}

inline const SpectralPowerDistribution& spd_white()
{
    /*
        switch (Options::get().white_point)
        {
        case WhitePoint::D50:
            return Illuminant::D50;
            break;
        case WhitePoint::D55:
            return Illuminant::D55;
            break;
        case WhitePoint::D60:
            return Illuminant::D60;
            break;
        case WhitePoint::D65:
            return Illuminant::D65;
            break;
        default:
            slyassert("unknown white point SPD requested");
            return Illuminant::ones;
            break;
        }
    */
    return Illuminant::D65;
}

class HeroSpectrum : public CoefficientSpectrum<4>
{
 public:
    HeroSpectrum(float v = 0.0f) : CoefficientSpectrum<4>(v) {}
    HeroSpectrum(const CoefficientSpectrum<4>& s) : CoefficientSpectrum<4>(s) {}
    HeroSpectrum(const ColorSpaceRGB& cs, float r, float g, float b,
                 HeroWavelengths hw, SpectrumType type = kReflectance)
    {
        RGB rgb(r, g, b);
        *this = fromRGB(cs, rgb, hw, type);
    }

    // create a new HeroSpectrum by inrerpolating the given SPD at the given
    // wavelengths
    static HeroSpectrum from_spd(const SpectralPowerDistribution& spd,
                                 const HeroWavelengths& hw)
    {
        HeroSpectrum hs;
        for (int i = 0; i < 4; ++i)
        {
            hs._c[i] = spd.value(hw.lambda[i]);
        }
        return hs;
    }

    static HeroSpectrum fromRGB(const ColorSpaceRGB& cs_, RGB rgb,
                                HeroWavelengths hw,
                                SpectrumType type = kReflectance)
    {
        // create a copy of the colour space for doing the conversion and
        // set its white point to E
        ColorSpaceRGB cs = cs_;
        cs.set_white_point(V2f(0.333f, 0.333f));

        if (type == kReflectance)
        {
            XYZ xyz = cs.RGBtoXYZ(rgb);
            HeroSpectrum s;
            for (int i = 0; i < 4; ++i)
            {
                s._c[i] = clamp(spectrum_xyz_to_p(hw.lambda[i], &xyz.x) /
                                    equal_energy_reflectance,
                                0, 1);
            }
            return s;
        }
        else
        {
            const SpectralPowerDistribution& spd = spd_white();
            if (rgb == RGB(1.0f))
            {
                // if pure white, just return the renderer's white point
                return from_spd(spd, hw);
            }
            else
            {
                XYZ xyz = cs.RGBtoXYZ(rgb);
                HeroSpectrum s;
                for (int i = 0; i < 4; ++i)
                {
                    s._c[i] = clamp(spectrum_xyz_to_p(hw.lambda[i], &xyz.x) /
                                        equal_energy_reflectance,
                                    0, 1) *
                              spd.value(hw.lambda[i]);
                }
                return s;
            }
        }
    }

    RGB toRGB(const ColorSpaceRGB& cs, const HeroWavelengths& hw,
              WhitePoint wp) const
    {
        return cs.XYZtoRGB(toXYZ(cs, hw, wp));
    }

    static HeroSpectrum fromXYZ(XYZ xyz, const ColorSpaceRGB& cs,
                                const HeroWavelengths& hw,
                                SpectrumType type = kReflectance)
    {
        if (type == kReflectance)
        {
            HeroSpectrum s;
            for (int i = 0; i < 4; ++i)
            {
                s._c[i] = clamp(spectrum_xyz_to_p(hw.lambda[i], &xyz.x) /
                                    equal_energy_reflectance,
                                0, 1);
            }
            return s;
        }
        else
        {
            const SpectralPowerDistribution& spd = spd_white();
            if (xyz == XYZ(1.0f))
            {
                // if pure white, just return the renderer's white point
                return from_spd(spd, hw);
            }
            else
            {
                HeroSpectrum s;
                for (int i = 0; i < 4; ++i)
                {
                    s._c[i] = clamp(spectrum_xyz_to_p(hw.lambda[i], &xyz.x) /
                                        equal_energy_reflectance,
                                    0, 1) *
                              spd.value(hw.lambda[i]);
                }
                return s;
            }
        }
    }

    static const SpectralPowerDistribution& getIlluminant(WhitePoint wp)
    {
        switch (wp)
        {
        case WhitePoint::D50:
            return Illuminant::D50;
            break;
        case WhitePoint::D55:
            return Illuminant::D55;
            break;
        case WhitePoint::D60:
            return Illuminant::D60;
            break;
        case WhitePoint::D65:
            return Illuminant::D65;
            break;
        default:
            slyassert("unknown white point SPD requested", "");
            return Illuminant::ones;
            break;
        }
    }

    XYZ toXYZ(const ColorSpaceRGB& cs, const HeroWavelengths& hw,
              WhitePoint wp) const
    {
        XYZ xyz(0.0f);
        const SpectralPowerDistribution& illum = getIlluminant(wp);
        const SpectralPowerDistribution& X = CMF::CIE_1931_2_x;
        const SpectralPowerDistribution& Y = CMF::CIE_1931_2_y;
        const SpectralPowerDistribution& Z = CMF::CIE_1931_2_z;
        float N = 0.0f;
        for (int i = 0; i < 4; ++i)
        {
            float e_i = illum.value(hw.lambda[i]);
            float x_i = X.value(hw.lambda[i]);
            float y_i = Y.value(hw.lambda[i]);
            float z_i = Z.value(hw.lambda[i]);
            float Me = _c[i] * e_i;
            xyz.x += x_i * Me;
            xyz.y += y_i * Me;
            xyz.z += z_i * Me;
            N += y_i * e_i;
        }
        slyassert(isfinite(xyz.x), "result is not finite");
        slyassert(N > 0.0f, "N is 0");
        return xyz / N;
    }

    static HeroSpectrum fromSampled(const SampledSpectrum& ss,
                                    const ColorSpaceRGB&, WhitePoint wp,
                                    const HeroWavelengths hw)
    {
        HeroSpectrum hs;
        for (int i = 0; i < 4; ++i)
        {
            hs[i] = ss.interpolate(hw[i]);
        }
        return hs;
    }

    friend std::ostream& operator<<(std::ostream& os, const HeroSpectrum& s);
};

inline std::ostream& operator<<(std::ostream& os, const HeroSpectrum& s)
{
    os << "{";
    for (size_t i = 0; i < 4; ++i)
    {
        os << s._c[i];
        if (i != 4 - 1)
        {
            os << ", ";
        }
    }
    os << "}";
    return os;
}

#if defined(SLY_SPECTRUM_HERO)
typedef HeroSpectrum Spectrum;
#else
typedef RGBSpectrum Spectrum;
#endif
}

#endif
