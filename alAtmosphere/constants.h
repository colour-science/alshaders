#ifndef SLY_CONSTANTS_H
#define SLY_CONSTANTS_H

#include <climits>
#include <limits>
#include <cmath>
#include <cassert>

namespace sly
{
static const float c_epsilon = 1e-5f;
static const float c_one_minus_epsilon = 0.9999994f;
static const float c_pi = M_PI;
static const float c_two_pi = M_PI * 2.0f;
static const float c_four_pi = M_PI * 4.0f;
static const float c_pi_over_two = M_PI * 0.5f;
static const float c_pi_over_180 = M_PI / 180.0f;
static const float c_one_over_pi = 1.0f / M_PI;
static const float c_one_over_sqrt_pi = 1.0f / sqrt(M_PI);
static const float c_one_over_two_pi = 1.0f / (2.0f * M_PI);
static const float c_infinity = std::numeric_limits<float>::infinity();
static const float c_huge = std::numeric_limits<float>::max();
}

#endif
