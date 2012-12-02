#include "Irawan.h"
#include "alUtil.h"

IrawanBRDF::IrawanBRDF(const WeavePattern& pattern)
: _pattern(pattern), _initialization(false), _specularNormalization(1.0f), _repeatU(10.0f), _repeatV(10.0f)
{
	AtUInt32 ns = 100000;
	AtRGB result = AI_RGB_BLACK;
	_initialization = true;
	AtRGB tmp;
	for (AtUInt32 i=0; i < ns; ++i)
	{
		AtVector wi = squareToCosineHemisphere(drand48(), drand48());
		AtVector wo = squareToCosineHemisphere(drand48(), drand48());
		result += eval(wi, wo, drand48(), drand48(), tmp) / cosTheta(wo);
	}
	_initialization = false;

	if (maxh(result) == 0.0f)
		_specularNormalization = 0.0f;
	else
		_specularNormalization = AtFloat(ns) / (maxh(result) * AI_PI);
}

AtRGB IrawanBRDF::eval(const AtVector& wi, const AtVector& wo, AtFloat uu, AtFloat vv, AtRGB& resultDiffuse) const
{
	resultDiffuse = AI_RGB_BLACK;
	// Make sure we're in the right hemisphere for the cloth
	// return black otherwise
	if (cosTheta(wi) <= 0 || cosTheta(wo) <= 0)
		return AI_RGB_BLACK;

	// scale the UVs by the size of the pattern
	AtPoint2 uv;
	uv.x = uu * _repeatU;
	uv.y = 	(1 - vv) * _repeatV;
	AtPoint2 xy;
	xy.x = uv.x * _pattern.tileWidth;
	xy.y = uv.y * _pattern.tileHeight;

	// modulo to find our coordinates relative to a single tile of the weave pattern
	/*
	Point2i lookup(
		modulo((int) xy.x, m_pattern.tileWidth),
		modulo((int) xy.y, m_pattern.tileHeight));
	*/
	int lookupx = modulo((int)xy.x, _pattern.tileWidth);
	int lookupy = modulo((int)xy.y, _pattern.tileHeight);

	// figure out which yarn in the pattern our lookup point is on
	int yarnID = _pattern.pattern[lookupx + lookupy * _pattern.tileWidth] - 1;

	// get the yarn
	const Yarn &yarn = _pattern.yarns.at(yarnID);
	// store center of the yarn segment
	AtPoint2 center;
	center.x = ((int) xy.x / _pattern.tileWidth) * _pattern.tileWidth + yarn.centerU * _pattern.tileWidth;
	center.y = 	((int) xy.y / _pattern.tileHeight) * _pattern.tileHeight
			+ (1 - yarn.centerV) * _pattern.tileHeight;

	// transform x and y to new coordinate system with (0,0) at the
	// center of the yarn segment
	xy.x =	  xy.x - center.x;
	xy.y = - (xy.y - center.y);

	int type = yarn.type;
	AtFloat w = yarn.width;
	AtFloat l = yarn.length;

	// Get incident and exitant directions.
	// These directions should be in the local coordinate system with Z normal to the surface
	AtVector om_i = wi;
	if (om_i.z < 0.0f) om_i = -om_i;
	AtVector om_r = wo;
	if (om_r.z < 0.0f) om_r = -om_r;

	AtFloat psi = yarn.psi;
	AtFloat umax = yarn.umax;
	AtFloat kappa = yarn.kappa;

	AtFloat dUmaxOverDWarp, dUmaxOverDWeft;
	AtRGB ks, kd;
	if (type == Yarn::kWarp)
	{
		dUmaxOverDWarp = _pattern.dWarpUmaxOverDWarp;
		dUmaxOverDWeft = _pattern.dWarpUmaxOverDWeft;
		ks = _warpSpecularColor;
		kd = _warpDiffuseColor;
	}
	else
	{ // type == EWeft
		dUmaxOverDWarp = _pattern.dWeftUmaxOverDWarp;
		dUmaxOverDWeft = _pattern.dWeftUmaxOverDWeft;
		// Rotate xy, incident, and exitant directions pi/2 radian about z-axis
		AtFloat tmp = xy.x;
		xy.x = -xy.y;
		xy.y = tmp;
		tmp = om_i.x;
		om_i.x = -om_i.y;
		om_i.y = tmp;
		tmp = om_r.x;
		om_r.x = -om_r.y;
		om_r.y = tmp;
		ks = _weftSpecularColor;
		kd = _weftDiffuseColor;
	}

	// Correlated (Perlin) noise.
	AtFloat random1 = 1.0f;
	AtFloat random2 = 1.0f;

	/* Number of TEA iterations (the more, the better the
	   quality of the pseudorandom floats) */
	const int teaIterations = 8;


/*
	if (_pattern.period > 0.0f)
	{
		// generate 1 seed per yarn segment
		AtUInt32 posx = center.x, posy = center.y;

		// we'll take a chance here that Arnold's implementation matches Mitsuba's closely enough
		AtPoint2 pp1, pp2;
		pp1.x = (center.x * (_pattern.tileHeight * _repeatV
				+ sampleTEAFloat(posx, 2*posy, teaIterations)) + center.y) / _pattern.period;
		pp1.y = 0.0f;
		pp2.x = (center.y * (_pattern.tileWidth * _repeatU
				+ sampleTEAFloat(posx, 2*posy+1, teaIterations)) + center.x) / _pattern.period;
		pp2.y = 0.0f;
		random1 = AiPerlin2(pp1);
		random2 = AiPerlin2(pp2);
		umax = umax + random1 * dUmaxOverDWarp + random2 * dUmaxOverDWeft;
	}
*/

	// Compute u and v.
	// See Chapter 6.
	AtFloat u = xy.y / (l / 2.0f) * umax;
	AtFloat v = xy.x * AI_PI / w;

	// Compute specular contribution.
	AtRGB result = AI_RGB_BLACK;

	AtFloat integrand = 0.0f;
	if (psi != 0.0f)
		integrand = evalStapleIntegrand(u, v, om_i, om_r, _pattern.alpha,
				_pattern.beta, psi, umax, kappa, w, l);
	else
		integrand = evalFilamentIntegrand(u, v, om_i, om_r, _pattern.alpha,
				_pattern.beta, _pattern.ss, umax, kappa, w, l);

	result = AI_RGB_WHITE * integrand;


	// Initialize random number generator based on texture location.
	AtFloat intensityVariation = 1.0f;
	if (_pattern.fineness > 0.0f)
	{
		// Compute random variation and scale specular component.
		// Generate fineness^2 seeds per 1 unit of texture.
		AtUInt32 index1 = (AtUInt32) ((center.x + xy.x) * _pattern.fineness);
		AtUInt32 index2 = (AtUInt32) ((center.y + xy.y) * _pattern.fineness);

		AtFloat xi = sampleTEAFloat(index1, index2, teaIterations);
		intensityVariation = std::min(-logf(xi), (AtFloat) 10.0f);
	}

	// if m_initialization is true, we are running an initialization pass to find the integral of the
	// specular reflectance under uniform diffuse illumination, which we'll use to normalize the specular
	// result during shading evaluation

	if (!_initialization)
		result = ks * (intensityVariation * integrand * _specularNormalization);
	else
		result = AI_RGB_WHITE * (intensityVariation * integrand);


	//result = AI_RGB_WHITE;
	if (type == Yarn::kWarp)
		result *= (_pattern.warpArea + _pattern.weftArea) / _pattern.warpArea;
	else
		result *= (_pattern.warpArea + _pattern.weftArea) / _pattern.weftArea;

	resultDiffuse = kd * fabsf(cosTheta(wo));

	return result * fabsf(cosTheta(wo));

}

AtFloat IrawanBRDF::pdf(const AtVector& wi, const AtVector& wo) const
{
	// check we're on the right side of the hemisphere
	if (cosTheta(wi) <= 0 || cosTheta(wo) <= 0)
		return 0.0f;

	// No good sampling strategy so we're using cosine hemisphere sampling
	return squareToCosineHemispherePdf(wi);
}

AtVector IrawanBRDF::sample(AtFloat r0, AtFloat r1) const
{
	/* Lacking a better sampling method, generate cosine-weighted directions */
	return squareToCosineHemisphere(r0, r1);
}

/** parameters:
 *	u	 to be compared to u(v) in texturing
 *	v	 for filament, we compute u(v)
 *	om_i  incident direction
 *	om_r  exitant direction
 *	ss	filament smoothing parameter
 *  fiber properties
 *	alpha uniform scattering
 *	beta  forward scattering
 *  yarn geometry
 *	psi   fiber twist angle; because this is filament, psi = pi/2
 *	umax  maximum inclination angle
 *	kappa spine curvature parameter
 *  weave pattern
 *	w	 width of segment rectangle
 *	l	 length of segment rectangle
 */
AtFloat IrawanBRDF::evalFilamentIntegrand(AtFloat u, AtFloat v, const AtVector& om_i,
		const AtVector& om_r, AtFloat alpha, AtFloat beta, AtFloat ss,
		AtFloat umax, AtFloat kappa, AtFloat w, AtFloat l) const
{
	// 0 <= ss < 1.0
	if (ss < 0.0f || ss >= 1.0f)
		return 0.0f;

	// w * sin(umax) < l
	if (w * std::sin(umax) >= l)
		return 0.0f;

	// -1 < kappa < inf
	if (kappa < -1.0f)
		return 0.0f;

	// h is the half vector
	AtVector h = AiV3Normalize(om_r + om_i);

	// u_of_v is location of specular reflection.
	AtFloat u_of_v = atan2f(h.y, h.z);

	// Check if u_of_v within the range of valid u values
	if (fabsf(u_of_v) < umax)
	{
		// n is normal to the yarn surface
		// t is tangent of the fibers.
		AtVector n;
		AiV3Create(n, std::sin(v), std::sin(u_of_v) * std::cos(v), std::cos(u_of_v) * std::cos(v));
		n = AiV3Normalize(n);
		AtVector t;
		AiV3Create(t, 0.0f, std::cos(u_of_v), -std::sin(u_of_v));
		t = AiV3Normalize(t);

		// R is radius of curvature.
		AtFloat R = radiusOfCurvature(std::min(std::abs(u_of_v),
			(1-ss)*umax), (1-ss)*umax, kappa, w, l);

		// G is geometry factor.
		AtFloat a = 0.5f * w;
		AtVector om_i_plus_om_r = om_i + om_r,
			   t_cross_h = AiV3Cross(t, h);

		AtFloat Gu = a * (R + a * std::cos(v)) /
			(AiV3Length(om_i_plus_om_r) * std::abs(t_cross_h.x));

		// fc is phase function
		AtFloat fc = alpha + vonMises(-AiV3Dot(om_i, om_r), beta);

		// A is attenuation function without smoothing.
		// As is attenuation function with smoothing.
		AtFloat A = seeliger(AiV3Dot(n, om_i), AiV3Dot(n, om_r), 0, 1);
		AtFloat As;
		if (ss == 0.0f)
			As = A;
		else
			As = A * (1.0f - smoothstep(0, 1, (std::abs(u_of_v)
				- (1.0f - ss) * umax) / (ss * umax)));

		// fs is scattering function.
		AtFloat fs = Gu * fc * As;


		// Domain transform.
		fs = fs * AI_PI * l;

		// Highlight has constant width delta_y on screen.
		AtFloat delta_y = l * _pattern.hWidth;

		// Clamp y_of_v between -(l - delta_y)/2 and (l - delta_y)/2.
		AtFloat y_of_v = u_of_v * 0.5f * l / umax;
		if (y_of_v > 0.5f * (l - delta_y))
			y_of_v = 0.5f * (l - delta_y);
		else if (y_of_v < 0.5f * (delta_y - l))
			y_of_v = 0.5f * (delta_y - l);

		// Check if |y(u(v)) - y(u)| < delta_y/2.
		if (fabsf(y_of_v - u * 0.5f * l / umax) < 0.5f * delta_y)
			return fs / delta_y;
		else return 0.0f;
	}
	return 0.0f;
}

/** parameters:
 *	u	 for staple, we compute v(u)
 *	v	 to be compared to v(u) in texturing
 *	om_i  incident direction
 *	om_r  exitant direction
 *  fiber properties
 *	alpha uniform scattering
 *	beta  forward scattering
 *  yarn geometry
 *	psi   fiber twist angle
 *	umax  maximum inclination angle
 *	kappa spine curvature parameter
 *  weave pattern
 *	w	 width of segment rectangle
 *	l	 length of segment rectangle
 */
AtFloat IrawanBRDF::evalStapleIntegrand(AtFloat u, AtFloat v, const AtVector& om_i,
		const AtVector& om_r, AtFloat alpha, AtFloat beta, AtFloat psi,
		AtFloat umax, AtFloat kappa, AtFloat w, AtFloat l) const
{
	// w * sin(umax) < l
	if (w * std::sin(umax) >= l)
		return 0.0f;

	// -1 < kappa < inf
	if (kappa < -1.0f)
		return 0.0f;

	// h is the half vector
	AtVector h = AiV3Normalize(om_i + om_r);

	// v_of_u is location of specular reflection.
	AtFloat D = (h.y*std::cos(u) - h.z*std::sin(u))
		/ (std::sqrt(h.x * h.x + std::pow(h.y * std::sin(u) + h.z * std::cos(u), (AtFloat) 2.0f)) * std::tan(psi));
	AtFloat v_of_u = std::atan2(-h.y * std::sin(u) - h.z * std::cos(u), h.x) + safe_acos(D);

	// Check if v_of_u within the range of valid v values
	if (std::abs(D) < 1.0f && std::abs(v_of_u) < AI_PI / 2.0f) {
		// n is normal to the yarn surface.
		// t is tangent of the fibers.

		AtVector n;
		AiV3Create(n, std::sin(v_of_u), std::sin(u)
				* std::cos(v_of_u), std::cos(u) * std::cos(v_of_u));
		n = AiV3Normalize(n);

		/*Vector t = normalize(Vector(-std::cos(v_of_u) * std::sin(psi),
				std::cos(u) * std::cos(psi) + std::sin(u) * std::sin(v_of_u) * std::sin(psi),
				-std::sin(u) * std::cos(psi) + std::cos(u) * std::sin(v_of_u) * std::sin(psi))); */

		// R is radius of curvature.
		AtFloat R = radiusOfCurvature(std::abs(u), umax, kappa, w, l);
		R = fabs(R);
		// G is geometry factor.
		AtFloat a = 0.5f * w;
		AtVector om_i_plus_om_r(om_i + om_r);
		AtFloat Gv = a * (R + a * std::cos(v_of_u))
			/ (AiV3Length(om_i_plus_om_r) * AiV3Dot(n, h) * std::abs(std::sin(psi)));

		// fc is phase function.
		AtFloat fc = alpha + vonMises(-AiV3Dot(om_i, om_r), beta);

		// A is attenuation function without smoothing.
		AtFloat A = seeliger(AiV3Dot(n, om_i), AiV3Dot(n, om_r), 0, 1);

		// fs is scattering function.
		AtFloat fs = Gv * fc * A;

		// Domain transform.
		fs = fs * 2.0f * w * umax;

		// Highlight has constant width delta_x on screen.
		AtFloat delta_x = w * _pattern.hWidth;

		// Clamp x_of_u between (w - delta_x)/2 and -(w - delta_x)/2.
		AtFloat x_of_u = v_of_u * w / AI_PI;
		if (x_of_u > 0.5f * (w - delta_x))
			x_of_u = 0.5f * (w - delta_x);
		else if (x_of_u < 0.5f * (delta_x - w))
			x_of_u = 0.5f * (delta_x - w);

		// Check if |x(v(u)) - x(v)| < delta_x/2.
		if (std::abs(x_of_u - v * w / AI_PI) < 0.5f * delta_x)
			return fs / delta_x;
	}
	return 0.0f;
}

AtFloat IrawanBRDF::radiusOfCurvature(AtFloat u, AtFloat umax, AtFloat kappa, AtFloat w, AtFloat l) const
{
	// rhat determines whether the spine is a segment
	// of an ellipse, a parabole, or a hyperbola.
	// See Section 5.3.
	AtFloat rhat = 1.0f + kappa * (1.0f + 1.0f / tanf(umax));

	AtFloat a = 0.5f * w, R = 0;
	if (rhat == 1.0f)
	{ // circle; see Subsection 5.3.1.
		R = 0.5f * l / sinf(umax) - a;
	}
	else if (rhat > 0.0f)
	{

		AtFloat tmax = atanf(rhat * tanf(umax));
		AtFloat bhat = (0.5f * l - a * sinf(umax)) / sinf(tmax);
		AtFloat ahat = bhat / rhat;
		AtFloat t = atanf(rhat * tanf(u));
		R = powf(bhat * bhat * cosf(t) * cosf(t)
		  + ahat * ahat * sinf(t) * sinf(t),(AtFloat) 1.5f) / (ahat * bhat);
	}
	else if (rhat < 0.0f)
	{ // hyperbola; see Subsection 5.3.3.

		AtFloat tmax = -atanhf(rhat * tanf(umax));
		AtFloat bhat = (0.5f * l - a * sinf(umax)) / sinhf(tmax);
		AtFloat ahat = bhat / rhat;
		AtFloat t = -atanhf(rhat * tanf(u));
		R = -powf(bhat * bhat * coshf(t) * coshf(t)
			+ ahat * ahat * sinhf(t) * sinhf(t), (AtFloat) 1.5f) / (ahat * bhat);
	}
	else
	{ // rhat == 0  // parabola; see Subsection 5.3.2.

		AtFloat tmax = tanf(umax);
		AtFloat ahat = (0.5f * l - a * sinf(umax)) / (2 * tmax);
		AtFloat t = tanf(u);
		R = 2 * ahat * powf(1 + t * t, (AtFloat) 1.5f);
	}
	return R;
}

// von Mises Distribution
AtFloat IrawanBRDF::vonMises(AtFloat cos_x, AtFloat b) const
{
	// assumes a = 0, b > 0 is a concentration parameter.

	AtFloat I0, absB = std::abs(b);
	if (std::abs(b) <= 3.75f) {
		AtFloat t = absB / 3.75f;
		t = t * t;
		I0 = 1.0f + t*(3.5156229f + t*(3.0899424f + t*(1.2067492f
				+ t*(0.2659732f + t*(0.0360768f + t*0.0045813f)))));
	} else {
		AtFloat t = 3.75f / absB;
		I0 = expf(absB) / std::sqrt(absB) * (0.39894228f + t*(0.01328592f
			+ t*(0.00225319f + t*(-0.00157565f + t*(0.00916281f + t*(-0.02057706f
			+ t*(0.02635537f + t*(-0.01647633f + t*0.00392377f))))))));
	}

	return expf(b * cos_x) / (2 * AI_PI * I0);
}

/// Attenuation term
AtFloat IrawanBRDF::seeliger(AtFloat cos_th1, AtFloat cos_th2, AtFloat sg_a, AtFloat sg_s) const
{
	AtFloat al = sg_s / (sg_a + sg_s); // albedo
	AtFloat c1 = std::max((AtFloat) 0, cos_th1);
	AtFloat c2 = std::max((AtFloat) 0, cos_th2);
	if (c1 == 0.0f || c2 == 0.0f)
		return 0.0f;
	return al / (4.0f * AI_PI) * c1 * c2 / (c1 + c2);
}

WeavePattern* createDenimPreset()
{
	WeavePattern* wv = new WeavePattern(
			"denim",
			3,			// tileWidth
			6,			// tileHeight
			0.01f,		// alpha
			4.0f,		// beta
			0.0f,		// ss
			0.5f,		// hWidth
			5.0f,		// warpArea
			1.0f,		// weftArea
			3.0f,		// fineness
			50, 50,		// repeatU, repeatV
			0, 0, 0, 0,	// dWarpdWarp, dWarpdWeft, dWeftdWarp, dWeftdWeft
			0			// period
	);

	AtUInt32 patterns[] =
	{
		1, 3, 8,
		1, 3, 5,
		1, 7, 5,
		1, 4, 5,
		6, 4, 5,
		2, 4, 5
	};
	for (int i = 0; i < wv->tileHeight * wv->tileWidth; i++)
		wv->pattern.push_back(patterns[i]);

	AtRGB blue = AiColorCreate(0.1, 0.1, 0.6);
	AtRGB white = AiColorCreate(0.7, 0.7, 0.7);

	wv->yarns.push_back(Yarn(Yarn::kWarp, -30, 12, 0, 1, 5, 0.1667, 0.75, blue, blue));
	wv->yarns.push_back(Yarn(Yarn::kWarp, -30, 12, 0, 1, 5, 0.1667, -0.25, blue, blue));
	wv->yarns.push_back(Yarn(Yarn::kWarp, -30, 12, 0, 1, 5, 0.5, 1.0833, blue, blue));
	wv->yarns.push_back(Yarn(Yarn::kWarp, -30, 12, 0, 1, 5, 0.5, 0.0833, blue, blue));
	wv->yarns.push_back(Yarn(Yarn::kWarp, -30, 12, 0, 1, 5, 0.8333, 0.4167, blue, blue));
	wv->yarns.push_back(Yarn(Yarn::kWeft, -30, 38, 0, 1, 1, 0.1667, 0.25, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWeft, -30, 38, 0, 1, 1, 0.5, 0.5833, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWeft, -30, 38, 0, 1, 1, 0.8333, 0.9167, white, white));
	return wv;
}

WeavePattern* createSilkCharmeusePreset()
{
	WeavePattern* wv = new WeavePattern("silk_charmeuse", 5, 10, 0.02f,7.3f, 0.5f, 0.5f, 9.0f, 1.0f, 3.0f, 50, 50,
											0.0f,0.0f,0.0f,0.0f, 0.0f);

	int patterns[] =
	{
		10, 2,  4,  6,  8,
		1,  2,  4,  6,  8,
		1,  2,  4, 13,  8,
		1,  2,  4,  7,  8,
		1, 11,  4,  7,  8,
		1,  3,  4,  7,  8,
		1,  3,  4,  7, 14,
		1,  3,  4,  7,  9,
		1,  3, 12,  7,  9,
		1,  3,  5,  7,  9
	};

	for (AtUInt32 i = 0; i < wv->tileHeight * wv->tileWidth; i++)
		wv->pattern.push_back(patterns[i]);

	AtRGB white = AiColorCreate(0.7, 0.7, 0.7);

	wv->yarns.push_back(Yarn(Yarn::kWarp, 0, 40, 2, 1, 9, 0.1, 0.45, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWarp, 0, 40, 2, 1, 9, 0.3, 1.05, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWarp, 0, 40, 2, 1, 9, 0.3, 0.05, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWarp, 0, 40, 2, 1, 9, 0.5, 0.65, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWarp, 0, 40, 2, 1, 9, 0.5, -0.35, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWarp, 0, 40, 2, 1, 9, 0.7, 1.25, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWarp, 0, 40, 2, 1, 9, 0.7, 0.25, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWarp, 0, 40, 2, 1, 9, 0.9, 0.85, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWarp, 0, 40, 2, 1, 9, 0.9, -0.15, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWeft, 0, 60, 0, 1, 1, 0.1, 0.95, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWeft, 0, 60, 0, 1, 1, 0.3, 0.55, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWeft, 0, 60, 0, 1, 1, 0.5, 0.15, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWeft, 0, 60, 0, 1, 1, 0.7, 0.75, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWeft, 0, 60, 0, 1, 1, 0.9, 0.35, white, white));

	return wv;
}

WeavePattern* createCottonTwillPreset()
{
	WeavePattern* wv = new WeavePattern("cotton_twill", 4, 8, 0.01f,4.0f, 0.0f, 0.5f, 6.0f, 2.0f, 4.0f, 50, 50,
											0.0f,0.0f,0.0f,0.0f, 0.0f);

	int patterns[] =
	{
		7, 2, 4, 6,
		7, 2, 4, 6,
		1, 8, 4, 6,
		1, 8, 4, 6,
		1, 3, 9, 6,
		1, 3, 9, 6,
		1, 3, 5, 10,
		1, 3, 5, 10
	};

	for (AtUInt32 i = 0; i < wv->tileHeight * wv->tileWidth; i++)
		wv->pattern.push_back(patterns[i]);

	AtRGB white = AiColorCreate(0.7, 0.7, 0.7);

	wv->yarns.push_back(Yarn(Yarn::kWarp, -30, 24, 0, 1, 6, 0.125, 0.375, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWarp, -30, 24, 0, 1, 6, 0.375, 1.125, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWarp, -30, 24, 0, 1, 6, 0.375, 0.125, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWarp, -30, 24, 0, 1, 6, 0.625, 0.875, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWarp, -30, 24, 0, 1, 6, 0.625, -0.125, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWarp, -30, 24, 0, 1, 6, 0.875, 0.625, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWeft, -30, 36, 0, 2, 1, 0.125, 0.875, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWeft, -30, 36, 0, 2, 1, 0.375, 0.625, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWeft, -30, 36, 0, 2, 1, 0.625, 0.375, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWeft, -30, 36, 0, 2, 1, 0.875, 0.125, white, white));

	return wv;
}

WeavePattern* createWoolGarbadinePreset()
{
	WeavePattern* wv = new WeavePattern("wool_gabardine", 6, 9, 0.01f,4.0f, 0.0f, 0.5f, 12.0f, 6.0f, 0.0f, 50,50,
											0.0f,0.0f,0.0f,0.0f, 0.0f);

	int patterns[] =
	{
		1, 1, 2, 2, 7, 7,
		1, 1, 2, 2, 7, 7,
		1, 1, 2, 2, 7, 7,
		1, 1, 6, 6, 4, 4,
		1, 1, 6, 6, 4, 4,
		1, 1, 6, 6, 4, 4,
		5, 5, 3, 3, 4, 4,
		5, 5, 3, 3, 4, 4,
		5, 5, 3, 3, 4, 4
	};

	for (AtUInt32 i = 0; i < wv->tileHeight * wv->tileWidth; i++)
		wv->pattern.push_back(patterns[i]);

	AtRGB white = AiColorCreate(0.7, 0.7, 0.7);

	wv->yarns.push_back(Yarn(Yarn::kWarp, 30, 30, 0, 2, 6, 0.167, 0.667, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWarp, 30, 30, 0, 2, 6, 0.5, 1.0, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWarp, 30, 30, 0, 2, 6, 0.5, 0.0, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWarp, 30, 30, 0, 2, 6, 0.833, 0.333, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWeft, 30, 30, 0, 3, 2, 0.167, 0.167, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWeft, 30, 30, 0, 3, 2, 0.5, 0.5, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWeft, 30, 30, 0, 3, 2, 0.833, 0.833, white, white));

	return wv;
}

WeavePattern* createPolyesterLiningClothPreset()
{
	WeavePattern* wv = new WeavePattern("polyester_lining_cloth", 2, 2, 0.015f,4.0f, 0.5f, 0.5f, 1.0f, 1.0f, 0.0f, 50,50,
											8.0f,8.0f,6.0f,6.0f, 50.0f);

	int patterns[] =
	{
		3, 2,
		1, 4
	};

	for (AtUInt32 i = 0; i < wv->tileHeight * wv->tileWidth; i++)
		wv->pattern.push_back(patterns[i]);

	AtRGB white = AiColorCreate(0.7, 0.7, 0.7);

	wv->yarns.push_back(Yarn(Yarn::kWarp, 0, 22, -0.7, 1, 1, 0.25, 0.25, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWarp, 0, 22, -0.7, 1, 1, 0.75, 0.75, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWeft, 0, 16, -0.7, 1, 1, 0.25, 0.75, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWeft, 0, 16, -0.7, 1, 1, 0.75, 0.25, white, white));

	return wv;
}

WeavePattern* createSilkShantungPreset()
{
	WeavePattern* wv = new WeavePattern("silk_shantung", 6, 8, 0.02f,1.5f, 0.5f, 0.5f, 8.0f, 16.0f, 0.0f, 50,50,
											20.0f,20.0f,10.0f,10.0f, 500.0f);

	int patterns[] =
	{
		3, 3, 3, 3, 2, 2,
		3, 3, 3, 3, 2, 2,
		3, 3, 3, 3, 2, 2,
		3, 3, 3, 3, 2, 2,
		4, 1, 1, 5, 5, 5,
		4, 1, 1, 5, 5, 5,
		4, 1, 1, 5, 5, 5,
		4, 1, 1, 5, 5, 5
	};

	for (AtUInt32 i = 0; i < wv->tileHeight * wv->tileWidth; i++)
		wv->pattern.push_back(patterns[i]);

	AtRGB white = AiColorCreate(0.7, 0.7, 0.7);

	wv->yarns.push_back(Yarn(Yarn::kWarp, 0, 50, -0.5, 2, 4, 0.3333, 0.25, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWarp, 0, 50, -0.5, 2, 4, 0.8333, 0.75, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWeft, 0, 23, -0.3, 4, 4, 0.3333, 0.75, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWeft, 0, 23, -0.3, 4, 4, -0.1667, 0.25, white, white));
	wv->yarns.push_back(Yarn(Yarn::kWeft, 0, 23, -0.3, 4, 4, 0.8333, 0.25, white, white));

	return wv;
}

// preset creation
WeavePattern* createWeavePreset(AtUInt32 preset)
{
	WeavePattern* wv = NULL;

	switch(preset)
	{
	case kSilkCharmeuse: //
		wv = createSilkCharmeusePreset();
		break;
	case kCottonTwill: //
		wv = createCottonTwillPreset();
		break;
	case kWoolGarbadine: //
		wv = createWoolGarbadinePreset();
		break;
	case kPolyesterLiningCloth: //
		wv = createPolyesterLiningClothPreset();
		break;
	case kSilkShantung: //
		wv = createSilkShantungPreset();
		break;
	default: //
		wv = createDenimPreset();
		break;

	}

	return wv;
}
