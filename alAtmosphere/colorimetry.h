#ifndef SLY_COLOR_H
#define SLY_COLOR_H

#include <cmath>
#include <cassert>
#include <vector>
#include <OpenEXR/ImathVec.h>
#include <OpenEXR/ImathMatrix.h>

#define slyassert(x, m) assert((x) && m)
namespace sly
{
using Imath::V2f;
using Imath::V3f;
using Imath::M33f;

/**
 * @brief Explicit RGB tristimulus struct
 */
struct RGB
{
    RGB(float v = 0.0f) : r(v), g(v), b(v) {}
    RGB(float r, float g, float b) : r(r), g(g), b(b) {}
    float r;
    float g;
    float b;

    RGB operator+(RGB o) const
    {
        RGB c(r + o.r, g + o.g, b + o.b);
        return c;
    }

    RGB operator+=(RGB o)
    {
        r += o.r;
        g += o.g;
        b += o.b;
        return *this;
    }

    RGB operator-(RGB o) const
    {
        RGB c(r - o.r, g - o.g, b - o.b);
        return c;
    }

    RGB operator-=(RGB o)
    {
        r -= o.r;
        g -= o.g;
        b -= o.b;
        return *this;
    }

    RGB operator*(RGB o) const
    {
        RGB c(r * o.r, g * o.g, b * o.b);
        return c;
    }

    RGB operator*=(RGB o)
    {
        r *= o.r;
        g *= o.g;
        b *= o.b;
        return *this;
    }

    RGB operator/(RGB o) const
    {
        RGB c(r / o.r, g / o.g, b / o.b);
        return c;
    }

    RGB operator/=(RGB o)
    {
        r /= o.r;
        g /= o.g;
        b /= o.b;
        return *this;
    }

    bool operator==(RGB o) { return r == o.r && g == o.g && b == o.b; }
};

inline RGB operator*(float f, RGB c) { return c * f; }

inline std::ostream& operator<<(std::ostream& os, RGB c)
{
    os << "RGB(" << c.r << ", " << c.g << ", " << c.b << ")";
    return os;
}

inline RGB pow(RGB c, float e)
{
    return RGB(powf(c.r, e), powf(c.g, e), powf(c.b, e));
}

inline float lerp(float a, float b, float t) { return (1 - t) * a + t * b; }
inline RGB lerp(RGB a, RGB b, float t) { return (1 - t) * a + t * b; }

inline float hmin(RGB c) { return std::min(c.r, std::min(c.g, c.b)); }

inline float hmax(RGB c) { return std::max(c.r, std::max(c.g, c.b)); }

inline RGB max(RGB c, float f)
{
    c.r = std::max(c.r, f);
    c.g = std::max(c.g, f);
    c.b = std::max(c.b, f);
    return c;
}

inline bool ispositive(RGB c) { return hmax(c) >= 0.0f; }

inline bool isnegative(RGB c) { return !ispositive(c); }

inline bool isfinite(RGB c)
{
    return std::isfinite(c.r) && std::isfinite(c.g) && std::isfinite(c.b);
}

inline bool ispfinite(RGB c) { return isfinite(c) && ispositive(c); }

inline float sRGB_curve(float c)
{
    if (c <= 0.0031308f)
    {
        return 12.92f * c;
    }
    else
    {
        return (1.0f + 0.055f) * powf(c, 1.0f / 2.4f) - 0.055f;
    }
}

inline RGB linear_to_sRGB(RGB c)
{
    c.r = sRGB_curve(c.r);
    c.g = sRGB_curve(c.g);
    c.b = sRGB_curve(c.b);
    return c;
}

/**
 * @brief XYZ tristimulus type
 * @details NOTE: by convention the domain of XYZ is [0,100] but we choose to
 * use [0,1] instead to save an extra factor every time we use the objects
 */
typedef V3f XYZ;

enum class WhitePoint
{
    D50,
    D55,
    D60,
    D65,
    ONE
};

struct ColorSpaceRGB
{
    ColorSpaceRGB(const std::string& name_, const V2f& r, const V2f& g,
                  const V2f& b, const V2f& w)
        : name(name_), red(r), green(g), blue(b), white(w)
    {
        initialize_matrix();
    }

    ColorSpaceRGB(const std::string& name_, const V2f& r, const V2f& g,
                  const V2f& b, const V2f& w, M33f m)
        : name(name_), red(r), green(g), blue(b), white(w), m_to_rgb(m)
    {
        m_to_xyz = m_to_rgb.inverse();
    }

    // Convert an XYZ to an RGB in this color space
    RGB XYZtoRGB(const XYZ& xyz) const;

    XYZ RGBtoXYZ(const RGB& rgb) const;

    void set_white_point(const V2f& w)
    {
        white = w;
        initialize_matrix();
    }

    std::string name;  /// name of the system, e.g. rec709, ACES etc
    V2f red;           /// chromaticity of the red primary
    V2f green;         /// chromaticity of the green primary
    V2f blue;          /// chromticity of the blue primary
    V2f white;         /// chromaticity of the white point
    M33f m_to_rgb;     /// XYZ to RGB matrix
    M33f m_to_xyz;     /// RGB to XYZ matrix
    WhitePoint wp;     /// enum for the white point

 private:
    void initialize_matrix()
    {
        float xr, yr, zr, xg, yg, zg, xb, yb, zb;
        float xw, yw, zw;
        float rx, ry, rz, gx, gy, gz, bx, by, bz;
        float rw, gw, bw;

        xr = red.x;
        yr = red.y;
        zr = 1 - (xr + yr);
        xg = green.x;
        yg = green.y;
        zg = 1 - (xg + yg);
        xb = blue.x;
        yb = blue.y;
        zb = 1 - (xb + yb);

        xw = white.x;
        yw = white.y;
        zw = 1 - (xw + yw);

        // xyz -> rgb matrix, before scaling to white
        rx = (yg * zb) - (yb * zg);
        ry = (xb * zg) - (xg * zb);
        rz = (xg * yb) - (xb * yg);
        gx = (yb * zr) - (yr * zb);
        gy = (xr * zb) - (xb * zr);
        gz = (xb * yr) - (xr * yb);
        bx = (yr * zg) - (yg * zr);
        by = (xg * zr) - (xr * zg);
        bz = (xr * yg) - (xg * yr);

        // White scaling factors.
        // Dividing by yw scales the white luminance to unity, as conventional
        rw = ((rx * xw) + (ry * yw) + (rz * zw)) / yw;
        gw = ((gx * xw) + (gy * yw) + (gz * zw)) / yw;
        bw = ((bx * xw) + (by * yw) + (bz * zw)) / yw;

        // xyz -> rgb matrix, correctly scaled to white
        m_to_rgb[0][0] = rx / rw;
        m_to_rgb[0][1] = ry / rw;
        m_to_rgb[0][2] = rz / rw;
        m_to_rgb[1][0] = gx / gw;
        m_to_rgb[1][1] = gy / gw;
        m_to_rgb[1][2] = gz / gw;
        m_to_rgb[2][0] = bx / bw;
        m_to_rgb[2][1] = by / bw;
        m_to_rgb[2][2] = bz / bw;

        m_to_xyz = m_to_rgb.inverse();
    }
};

namespace Chromaticies
{
extern V2f D60;
extern V2f D65;
extern V2f E;
extern V2f ACES;
extern V2f P3;
}

namespace Primaries
{
extern ColorSpaceRGB Rec709;
extern ColorSpaceRGB Rec709E;
extern ColorSpaceRGB Rec2020;
extern ColorSpaceRGB ACES;
extern ColorSpaceRGB ACEScg;
extern ColorSpaceRGB P3;
}

/**
 * @brief Generalized spectral power distribution container
 */
class SpectralPowerDistribution
{
 public:
    explicit SpectralPowerDistribution(float start, float end, float step,
                                       float v = 0.0f)
        : _step(step)
    {
        slyassert(end >= start + step, "end of range is less than step");

        for (float l = start; l < end; l += step)
        {
            _wavelengths.push_back(l);
            _values.push_back(v);
        }
    }

    SpectralPowerDistribution(std::vector<float>&& w, std::vector<float>&& v)
    {
        slyassert(w.size() == v.size(),
                  "wavelengh and value array sizes do not match");
        slyassert(w.size() > 1, "number of wavelengths must be greater than 1");

        _wavelengths = std::move(w);
        _values = std::move(v);

        // scan through the wavelengths and determine if it's a uniform
        // distribution or not
        _step = _wavelengths[1] - _wavelengths[0];
        for (size_t i = 2; i < _wavelengths.size(); ++i)
        {
            if (_wavelengths[i] - _wavelengths[i - 1] != _step)
            {
                _step = 0.0f;
                break;
            }
        }
    }

    SpectralPowerDistribution(float start, float end, float step,
                              std::vector<float>&& v)
    {
        _step = step;
        float end_inc = end - step;
        slyassert(int((start - end_inc) / step) == v.size(),
                  "range does not match size of values array");

        _values = std::move(v);

        for (float l = start; l < end; l += step)
        {
            _wavelengths.push_back(l);
        }

        slyassert(_wavelengths.size() == _values.size(),
                  "number of wavelengths does not match number of values");
    }

    bool operator==(const SpectralPowerDistribution& rhs)
    {
        if (start() != rhs.start() || end() != rhs.end() ||
            step() != rhs.step())
            return false;

        for (size_t i = 0; i < numSamples(); ++i)
        {
            if (_wavelengths[i] != rhs._wavelengths[i]) return false;
            if (_values[i] != rhs._values[i]) return false;
        }

        return true;
    }

    bool operator!=(const SpectralPowerDistribution& rhs)
    {
        return !(*this == rhs);
    }

    float start() const { return _wavelengths[0]; }
    float end() const { return _wavelengths[_wavelengths.size() - 1] + _step; }
    float step() const { return _step; }
    size_t numSamples() const { return _wavelengths.size(); }
    bool isUniform() const { return !(_step == 0.0f); }

    /**
     * @brief Check if rhs has the same range and step size as this
     */
    bool isEqualScale(const SpectralPowerDistribution& rhs) const
    {
        if (start() != rhs.start() || end() != rhs.end() ||
            step() != rhs.step())
            return false;

        for (size_t i = 0; i < _wavelengths.size(); ++i)
        {
            if (_wavelengths[i] != rhs._wavelengths[i]) return false;
        }

        return true;
    }

    friend std::ostream& operator<<(std::ostream& os,
                                    const SpectralPowerDistribution& spd);

    void write(std::ostream& os) const
    {
        for (size_t i = 0; i < _wavelengths.size(); ++i)
        {
            os << _wavelengths[i] << " " << _values[i] << std::endl;
        }
    }

    // interpolate this SPD from another
    // currently rhs must have a wider range than this spd (no extrapolation)
    void interpolateFrom(const SpectralPowerDistribution& rhs)
    {
        size_t sz = _wavelengths.size();

        // if the other SPD is the same distribution as us, just copy it
        if (isEqualScale(rhs))
        {
            *this = rhs;
            return;
        }

        // first check that the range of the SPD we're interpolating is
        // at least as big as the range we're representing
        // i.e. we don't do extrapolation yet
        slyassert(_wavelengths[0] >= _wavelengths[0], "");
        slyassert(_wavelengths[sz - 1] <=
                      rhs._wavelengths[rhs._wavelengths.size() - 1],
                  "");
        slyassert(rhs._wavelengths.size() > 1, "");

        // foreach wavelength we want to interpolate onto
        for (size_t i = 0; i < sz; ++i)
        {
            size_t j1 = 0;
            size_t j0 = 0;
            // find the first wavelength in rhs that's past our desired
            // wavelength
            while (rhs._wavelengths[j1++] < _wavelengths[i])
                ;

            // j is now the index of the first wavelength greater than our
            // desired
            // get the preceding index
            if (j1 > 0)
            {
                j0 = j1 - 1;

                // get the t of where our wavelength sits between the rhs ones
                float t = (_wavelengths[i] - rhs._wavelengths[j0]) /
                          (rhs._wavelengths[j1] - rhs._wavelengths[j0]);

                // lerp to get our new value
                _values[i] = lerp(rhs._values[j0], rhs._values[j1], t);
            }
            else
            {
                _values[i] = rhs._values[j1];
            }
        }
    }

    /**
     * @brief Interpolate this SPD onto a simple float array of values
     * @details Used by the CoefficientSpectrum class to fill out its
     * precomputed
     * arrays
     */
    void interpolateOnto(float* v, float lambda_start, float lambda_end,
                         int sz) const
    {
        slyassert(lambda_start >= start(), "");
        slyassert(lambda_end <= end(), "");

        // if sampling is the same, just copy
        if (sz == _wavelengths.size() && lambda_start == start() &&
            lambda_end == end())
            memcpy(v, &(_wavelengths[0]), sizeof(float) * sz);

        float lambda = lambda_start;
        float lambda_step = (lambda_end - lambda_start) / float(sz);
        for (int i = 0; i < sz; ++i, lambda += lambda_step)
        {
            int j1 = 0;
            int j0 = 0;
            while (_wavelengths[j1++] < lambda)
                ;

            if (j1 > 0)
            {
                j0 = j1 - 1;
                float t = (lambda - _wavelengths[j0]) /
                          (_wavelengths[j1] - _wavelengths[j0]);

                v[i] = lerp(_values[j0], _values[j1], t);
            }
            else
            {
                v[i] = _values[0];
            }
        }
    }

    //
    float value(float lambda) const
    {
        slyassert(lambda >= start(), "");
        slyassert(lambda <= end(), "");
        slyassert(_step != 0, "");

        float s = start();
        float e = end();
        float sz = e - s;
        int num = _wavelengths.size();
        float t = (lambda - s) / (sz);
        int idx = int(t * num);
        return _values[idx];
    }

    /**
     * @brief Convert this SPD to XYZ coordinates using the supplied CMFs
     * @details NOTE: this function assumes the supplied CMFs have the same
     * scale
     * as this SPD
     */
    XYZ toXYZ(const SpectralPowerDistribution& cie_x,
              const SpectralPowerDistribution& cie_y,
              const SpectralPowerDistribution& cie_z);

    /**
     * @brief Convert this SPD to XYZ coordinates using the supplied CMFs and
     * illuminant SPD
     */
    XYZ toXYZ(const SpectralPowerDistribution& cie_x,
              const SpectralPowerDistribution& cie_y,
              const SpectralPowerDistribution& cie_z,
              const SpectralPowerDistribution& illum);

 private:
    // step size between wavelength samples. will be zero if non-uniform
    float _step;
    std::vector<float> _wavelengths;
    std::vector<float> _values;
};

inline std::ostream& operator<<(std::ostream& os,
                                const SpectralPowerDistribution& spd)
{
    //  os << "{";
    //  for (size_t i = 0; i < spd._wavelengths.size(); ++i)
    //  {
    //      os << boost::format("%.2f: %.2f\n") % spd._wavelengths[i] %
    //                spd._values[i];
    //  }
    //  os << "}";
    //  return os;
}

namespace Illuminant
{
extern const SpectralPowerDistribution D50;
extern const SpectralPowerDistribution D55;
extern const SpectralPowerDistribution D60;
extern const SpectralPowerDistribution D65;
extern const SpectralPowerDistribution ones;
}

namespace CMF
{
extern const SpectralPowerDistribution CIE_1931_2_x;
extern const SpectralPowerDistribution CIE_1931_2_y;
extern const SpectralPowerDistribution CIE_1931_2_z;
extern const SpectralPowerDistribution CIE_1964_10_x;
extern const SpectralPowerDistribution CIE_1964_10_y;
extern const SpectralPowerDistribution CIE_1964_10_z;
extern const SpectralPowerDistribution CIE_2012_2_x;
extern const SpectralPowerDistribution CIE_2012_2_y;
extern const SpectralPowerDistribution CIE_2012_2_z;
extern const SpectralPowerDistribution CIE_2012_10_x;
extern const SpectralPowerDistribution CIE_2012_10_y;
extern const SpectralPowerDistribution CIE_2012_10_z;

extern const SpectralPowerDistribution SmitsWhite;
extern const SpectralPowerDistribution SmitsRed;
extern const SpectralPowerDistribution SmitsGreen;
extern const SpectralPowerDistribution SmitsBlue;
extern const SpectralPowerDistribution SmitsCyan;
extern const SpectralPowerDistribution SmitsMagenta;
extern const SpectralPowerDistribution SmitsYellow;
}

namespace ColorChecker
{
namespace BabelAverage
{
extern const SpectralPowerDistribution dark_skin;
extern const SpectralPowerDistribution light_skin;
extern const SpectralPowerDistribution blue_sky;
extern const SpectralPowerDistribution foliage;
extern const SpectralPowerDistribution blue_flower;
extern const SpectralPowerDistribution bluish_green;
extern const SpectralPowerDistribution orange;
extern const SpectralPowerDistribution purplish_blue;
extern const SpectralPowerDistribution moderate_red;
extern const SpectralPowerDistribution purple;
extern const SpectralPowerDistribution yellow_green;
extern const SpectralPowerDistribution orange_yellow;
extern const SpectralPowerDistribution blue;
extern const SpectralPowerDistribution green;
extern const SpectralPowerDistribution red;
extern const SpectralPowerDistribution yellow;
extern const SpectralPowerDistribution magenta;
extern const SpectralPowerDistribution cyan;
extern const SpectralPowerDistribution white_95;
extern const SpectralPowerDistribution neutral_80;
extern const SpectralPowerDistribution neutral_65;
extern const SpectralPowerDistribution neutral_50;
extern const SpectralPowerDistribution neutral_35;
extern const SpectralPowerDistribution black_20;
}

namespace Ohta1997
{
extern const SpectralPowerDistribution dark_skin;
extern const SpectralPowerDistribution light_skin;
extern const SpectralPowerDistribution blue_sky;
extern const SpectralPowerDistribution foliage;
extern const SpectralPowerDistribution blue_flower;
extern const SpectralPowerDistribution bluish_green;
extern const SpectralPowerDistribution orange;
extern const SpectralPowerDistribution purplish_blue;
extern const SpectralPowerDistribution moderate_red;
extern const SpectralPowerDistribution purple;
extern const SpectralPowerDistribution yellow_green;
extern const SpectralPowerDistribution orange_yellow;
extern const SpectralPowerDistribution blue;
extern const SpectralPowerDistribution green;
extern const SpectralPowerDistribution red;
extern const SpectralPowerDistribution yellow;
extern const SpectralPowerDistribution magenta;
extern const SpectralPowerDistribution cyan;
extern const SpectralPowerDistribution white_95;
extern const SpectralPowerDistribution neutral_80;
extern const SpectralPowerDistribution neutral_65;
extern const SpectralPowerDistribution neutral_50;
extern const SpectralPowerDistribution neutral_35;
extern const SpectralPowerDistribution black_20;
}
}

inline XYZ SpectralPowerDistribution::toXYZ(
    const SpectralPowerDistribution& cie_x,
    const SpectralPowerDistribution& cie_y,
    const SpectralPowerDistribution& cie_z)
{
    XYZ xyz(0.0f);
    slyassert(isEqualScale(cie_x), "cie_x SPD is not equal scale.");
    slyassert(isEqualScale(cie_y), "cie_y SPD is not equal scale.");
    slyassert(isEqualScale(cie_z), "cie_z SPD is not equal scale.");

    for (size_t i = 0; i < _wavelengths.size(); ++i)
    {
        xyz.x += _values[i] * cie_x._values[i];
        xyz.y += _values[i] * cie_y._values[i];
        xyz.z += _values[i] * cie_z._values[i];
    }

    return xyz;
}

inline XYZ SpectralPowerDistribution::toXYZ(
    const SpectralPowerDistribution& cie_x,
    const SpectralPowerDistribution& cie_y,
    const SpectralPowerDistribution& cie_z,
    const SpectralPowerDistribution& illum)
{
    XYZ xyz(0.0f);
    slyassert(isEqualScale(cie_x), "cie_x SPD is not equal scale.");
    slyassert(isEqualScale(cie_y), "cie_y SPD is not equal scale.");
    slyassert(isEqualScale(cie_z), "cie_z SPD is not equal scale.");
    slyassert(isEqualScale(illum), "illum SPD is not equal scale.");

    float N = 0.0f;
    for (size_t i = 0; i < _values.size(); ++i)
    {
        float Me = _values[i] * illum._values[i];
        xyz.x += Me * cie_x._values[i];
        xyz.y += Me * cie_y._values[i];
        xyz.z += Me * cie_z._values[i];
        N += cie_y._values[i] * illum._values[i];
    }

    return xyz / N;
}
}

#endif
