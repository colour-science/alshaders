#pragma once

#include "alUtil.h"
#include <ai.h>
#include <string>
#include <vector>

// This code is an implementation of Irawan & Marschner, "Specular Reflection from Woven Cloth":
// http://dl.acm.org/citation.cfm?id=2077352
// The implementation here is a port of the algorithm from Mitsuba, by Wenzel Jakob:
// http://www.mitsuba-renderer.org/
// which is itself port of an original Java implementation by Piti Irawan
// The algorithm is reproduced here with the permission of Irawan and Jakob

/// Data structure describing the properties of a single yarn
struct Yarn
{
	enum YarnType
	{
		kWarp = 0,
		kWeft = 1
	};

	/// Type of yarn (warp or weft)
	YarnType type;
	/// Fiber twist angle
	AtFloat psi;
	// Maximum inclination angle
	AtFloat umax;
	/// Spine curvature
	AtFloat kappa;
	/// Width of segment rectangle
	AtFloat width;
	/// Length of segment rectangle
	AtFloat length;
	/*! u coordinate of the yarn segment center,
	 * assumes that the tile covers 0 <= u, v <= 1.
	 * (0, 0) is lower left corner of the weave pattern
	 */
	AtFloat centerU;
	/// v coordinate of the yarn segment center
	AtFloat centerV;
	/// Diffuse color
	AtRGB kd;
	/// Specular color
	AtRGB ks;

	Yarn() : type(kWarp),
		psi(0), umax(0), kappa(0), width(0), length(0),
		centerU(0), centerV(0), kd(AI_RGB_BLACK), ks(AI_RGB_BLACK) { }

	Yarn(YarnType type_, AtFloat psi_, AtFloat umax_, AtFloat kappa_, AtFloat width_, AtFloat length_, AtFloat centerU_,
			AtFloat centerV_, const AtRGB& kd_, const AtRGB& ks_)
	: type(type_), psi(psi_), umax(umax_), kappa(kappa_), width(width_), length(length_), centerU(centerU_),
	  centerV(centerV_), kd(kd_), ks(ks_)
	{}

// weave file processing. leave out for now...
#if 0
	Yarn(Stream *stream) {
		type = (EYarnType) stream->readInt();
		psi = stream->readFloat();
		umax = stream->readFloat();
		kappa = stream->readFloat();
		width = stream->readFloat();
		length = stream->readFloat();
		centerU = stream->readFloat();
		centerV = stream->readFloat();
		kd = Spectrum(stream);
		ks = Spectrum(stream);
	}

	void serialize(Stream *stream) const {
		stream->writeInt(type);
		stream->writeFloat(psi);
		stream->writeFloat(umax);
		stream->writeFloat(kappa);
		stream->writeFloat(width);
		stream->writeFloat(length);
		stream->writeFloat(centerU);
		stream->writeFloat(centerV);
		kd.serialize(stream);
		ks.serialize(stream);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "yarn {" << endl
			<< "  type = " << ((type == EWarp) ? "warp" : "weft") << "," << endl;
		if (psi != 0)
			oss << "  /* Fiber twist angle */" << endl
				<< "  psi = " << psi * 180 / M_PI << "," << endl;
		oss << "  /* Maximum inclination angle */" << endl
			<< "  umax = " << umax * 180 / M_PI << "," << endl;
		if (kappa != 0)
			oss << "  /* Spine curvature */" << endl
				<< "  kappa = " << kappa << "," << endl;
		oss << "  /* Width and length of the segment rectangle */" << endl
			<< "  width = " << width << "," << endl
			<< "  length = " << length << "," << endl
			<< "  /* Yarn segment center in tile space */" << endl
			<< "  centerU = " << centerU << "," << endl
			<< "  centerV = " << centerV << "," << endl
			<< "  /* Diffuse and specular color */" << endl;
		Float r, g, b;
		kd.toLinearRGB(r, g, b);
		oss << "  kd = {" << r << ", " << g << ", " << b << "}," << endl;
		ks.toLinearRGB(r, g, b);
		oss << "  ks = {" << r << ", " << g << ", " << b << "}," << endl
			<< "}";
		return oss.str();
	}
#endif
};

/// Data structure describing a weave pattern
struct WeavePattern
{
	/// Name of the weave pattern
	std::string name;
	/// Uniform scattering parameter
	AtFloat alpha;
	/// Forward scattering parameter
	AtFloat beta;
	/// Filament smoothing
	AtFloat ss;
	/// Highlight width
	AtFloat hWidth;
	/// Combined area taken up by the warp & weft
	AtFloat warpArea, weftArea;

	/// Size of the weave pattern
	AtUInt32 tileWidth, tileHeight;

	/* Noise-related parameters */
	AtFloat dWarpUmaxOverDWarp;
	AtFloat dWarpUmaxOverDWeft;
	AtFloat dWeftUmaxOverDWarp;
	AtFloat dWeftUmaxOverDWeft;
	AtFloat fineness, period;

	/// Detailed weave pattern
	std::vector<AtUInt32> pattern;

	/// List of all yarns referenced in \c pattern
	std::vector<Yarn> yarns;

	/// patter repitition
	AtFloat repeatU, repeatV;

	inline WeavePattern() : name(""), tileWidth(0), tileHeight(0),
		alpha(0), beta(0), ss(0), hWidth(0),
		warpArea(0), weftArea(0),
		dWarpUmaxOverDWarp(0), dWarpUmaxOverDWeft(0),
		dWeftUmaxOverDWarp(0), dWeftUmaxOverDWeft(0),
		fineness(0), period(0), repeatU(1.0f), repeatV(1.0f) { }

	inline WeavePattern(const std::string& name_, AtFloat tileWidth_, AtFloat tileHeight_, AtFloat alpha_, AtFloat beta_,
						AtFloat ss_,  AtFloat hWidth_, AtFloat warpArea_,
							AtFloat weftArea_, AtFloat fineness_, AtFloat repeatU_, AtFloat repeatV_, AtFloat dWarpdWarp,
							AtFloat dWarpdWeft, AtFloat dWeftdWarp, AtFloat dWeftdWeft,
							AtFloat period_)
	: 	name(name_),
		alpha(alpha_), beta(beta_), ss(ss_), hWidth(hWidth_),
		warpArea(warpArea_), weftArea(weftArea_), tileWidth(tileWidth_), tileHeight(tileHeight_),
		dWarpUmaxOverDWarp(dWarpdWarp), dWarpUmaxOverDWeft(dWarpdWeft),
		dWeftUmaxOverDWarp(dWeftdWarp), dWeftUmaxOverDWeft(dWeftdWeft),
		fineness(fineness_), period(period_), repeatU(repeatU_), repeatV(repeatV_) { }

// we'll handle weave pattern file loading later...
#if 0
	WeavePattern(Stream *stream) {
		name = stream->readString();
		alpha = stream->readFloat();
		beta = stream->readFloat();
		ss = stream->readFloat();
		hWidth = stream->readFloat();
		warpArea = stream->readFloat();
		weftArea = stream->readFloat();
		tileWidth = stream->readUInt();
		tileHeight = stream->readUInt();
		dWarpUmaxOverDWarp = stream->readFloat();
		dWarpUmaxOverDWeft = stream->readFloat();
		dWeftUmaxOverDWarp = stream->readFloat();
		dWeftUmaxOverDWeft = stream->readFloat();
		fineness = stream->readFloat();
		period = stream->readFloat();
		pattern.resize(tileWidth * tileHeight);
		stream->readUIntArray(&pattern[0], pattern.size());
		size_t yarnCount = stream->readSize();
		yarns.resize(yarnCount);
		for (size_t i=0; i<yarnCount; ++i)
			yarns[i] = Yarn(stream);
	}

	void serialize(Stream *stream) const {
		stream->writeString(name);
		stream->writeFloat(alpha);
		stream->writeFloat(beta);
		stream->writeFloat(ss);
		stream->writeFloat(hWidth);
		stream->writeFloat(warpArea);
		stream->writeFloat(weftArea);
		stream->writeUInt(tileWidth);
		stream->writeUInt(tileHeight);
		stream->writeFloat(dWarpUmaxOverDWarp);
		stream->writeFloat(dWarpUmaxOverDWeft);
		stream->writeFloat(dWeftUmaxOverDWarp);
		stream->writeFloat(dWeftUmaxOverDWeft);
		stream->writeFloat(fineness);
		stream->writeFloat(period);
		stream->writeUIntArray(&pattern[0], pattern.size());
		stream->writeSize(yarns.size());
		for (size_t i=0; i<yarns.size(); ++i)
			yarns[i].serialize(stream);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "weave {" << endl
			<< "  name = \"" << name << "\"," << endl << endl
			<< "  /* Tile size of the weave pattern */" << endl
			<< "  tileWidth = " << tileWidth << "," << endl
			<< "  tileHeight = " << tileHeight << "," << endl << endl
			<< "  /* Uniform and forward scattering parameters */" << endl
			<< "  alpha = " << alpha << "," << endl
			<< "  beta = " << beta << "," << endl << endl
			<< "  /* Filament smoothing */" << endl
			<< "  ss = " << ss << "," << endl << endl
			<< "  /* Highlight width */" << endl
			<< "  hWidth = " << hWidth << "," << endl << endl
			<< "  /* Combined warp/weft size */" << endl
			<< "  warpArea = " << warpArea << "," << endl
			<< "  weftArea = " << weftArea << "," << endl << endl
			<< "  /* Noise-related parameters */" << endl;
		if (dWarpUmaxOverDWarp != 0)
			oss << "  dWarpUmaxOverDWarp = " << dWarpUmaxOverDWarp * 180 / M_PI << "," << endl;
		if (dWeftUmaxOverDWeft != 0)
			oss << "  dWarpUmaxOverDWeft = " << dWarpUmaxOverDWeft * 180 / M_PI << "," << endl;
		if (dWarpUmaxOverDWarp != 0)
			oss << "  dWeftUmaxOverDWarp = " << dWeftUmaxOverDWarp * 180 / M_PI << "," << endl;
		if (dWeftUmaxOverDWeft != 0)
			oss << "  dWeftUmaxOverDWeft = " << dWeftUmaxOverDWeft * 180 / M_PI << "," << endl;
		if (fineness != 0)
			oss << "  fineness = " << fineness << "," << endl;
		if (period != 0)
			oss << "  period = " << period << "," << endl;
		oss << endl
			<< "  /* Weave pattern description */" << endl
			<< "  pattern {" << endl
			<< "    ";
		for (size_t i=0; i<pattern.size(); ++i) {
			oss << (int) pattern[i];
			if (i+1<pattern.size())
				oss << ", ";
		}
		oss << endl
			<< "  }," << endl
			<< endl
			<< "  /* Listing of all used yarns */" << endl;
		for (size_t i=0; i<yarns.size(); ++i) {
			oss << "  " << indent(yarns[i].toString());
			if (i+1<yarns.size())
				oss << "," << endl;
			oss << endl;
		}

		oss << "}";
		return oss.str();
	}
#endif
};


class IrawanBRDF
{
public:
	IrawanBRDF(const WeavePattern& pattern);

	void setFrame(const AtVector& n, const AtVector& s, const AtVector& t)
	{
		_n = n;
		_s = s;
		_t = t;
	}

	void setRepeat(AtFloat repeatU, AtFloat repeatV)
	{
		_repeatU = repeatU;
		_repeatV = repeatV;
	}

	void setColors(AtRGB warpColor, AtRGB weftColor)
	{
		_warpColor = warpColor;
		_weftColor = weftColor;
	}

	// Transform a world-space direction into the local shading coordinate system, where the normal, _n is the z axis
	// This isn't how we generally do things in Arnold but it will make the rest of the implementation a lot easier...
	inline AtVector worldToLocal(const AtVector& v)
	{
		//return AtVector(AiV3Dot(v, _s), AiV3Dot(v, _t), AiV3Dot(v, _n));
		AtVector l;
		l.x = AiV3Dot(v, _s);
		l.y = AiV3Dot(v, _t);
		l.z = AiV3Dot(v, _n);
		return l;
	}

	// Transform from local shading space to world space. This is necessary so we can get vectors back into Arnold for
	// tracing rays
	inline AtVector localToWorld(const AtVector& v)
	{
		/*
		return AtVector(
				_s.x*v.x + _t.x*v.y + _n.x*v.z,
				_s.y*v.x + _t.y*v.y + _n.y*v.z,
				_s.z*v.x + _t.z*v.y + _n.z*v.z
		);
		*/
		AtVector w;
		w.x = _s.x*v.x + _t.x*v.y + _n.x*v.z;
		w.y = _s.y*v.x + _t.y*v.y + _n.y*v.z;
		w.z = _s.z*v.x + _t.z*v.y + _n.z*v.z;
		return w;
	}

	// Get the cosine of the angle between the given vector (in shading space) and the cloth normal
	inline AtFloat cosTheta(const AtVector& v) const
	{
		return v.z;
	}


	AtVector squareToCosineHemisphere(AtFloat u1, AtFloat u2) const
	{
	   AtVector ret;
	   concentricSampleDisk(u1, u2, ret.x, ret.y);
	   ret.z = sqrtf(std::max(0.0f, 1.0f - ret.x*ret.x - ret.y*ret.y));
	   return ret;
	}

	AtFloat squareToCosineHemispherePdf(const AtVector& d) const
	{
		return AI_ONEOVERPI * cosTheta(d);
	}


	// Evaluate the BRDF for the given pair of directions
	AtRGB eval(const AtVector& wi, const AtVector& wo, AtFloat u, AtFloat v) const;

	// Calculate the PDF for the given pair of directions
	AtFloat pdf(const AtVector& wi, const AtVector& wo) const;

	// Create a sample for this BRDF
	AtVector sample(AtFloat r0, AtFloat r1) const;

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
	AtFloat evalFilamentIntegrand(AtFloat u, AtFloat v, const AtVector& om_i, const AtVector& om_r, AtFloat alpha,
			AtFloat beta, AtFloat ss, AtFloat umax, AtFloat kappa, AtFloat w, AtFloat l) const;

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
	AtFloat evalStapleIntegrand(AtFloat u, AtFloat v, const AtVector& om_i, const AtVector& om_r, AtFloat alpha, AtFloat beta, AtFloat psi,
			AtFloat umax, AtFloat kappa, AtFloat w, AtFloat l) const;

	AtFloat radiusOfCurvature(AtFloat u, AtFloat umax, AtFloat kappa, AtFloat w, AtFloat l) const;

	inline AtFloat atanh(AtFloat arg) const
	{
		return logf((1.0f + arg) / (1.0f - arg)) / 2.0f;
	}

	// von Mises Distribution
	AtFloat vonMises(AtFloat cos_x, AtFloat b) const;

	/// Attenuation term
	AtFloat seeliger(AtFloat cos_th1, AtFloat cos_th2, AtFloat sg_a, AtFloat sg_s) const;

private:
	// weave pattern
	WeavePattern _pattern;
	// local coordinate frame
	// n becomes the z axis
	AtVector _n;
	AtVector _s;
	AtVector _t;
	// are we doing an initialization pass to calculate the specular normalization factor?
	bool _initialization;
	AtFloat _specularNormalization;
	AtFloat _repeatU;
	AtFloat _repeatV;
	AtRGB _warpColor;
	AtRGB _weftColor;
};

// utility function for creating preset weave patterns
enum WeavePreset
{
	kDenim = 0,
	kSilkCharmeuse,
	kCottonTwill,
	kWoolGarbadine,
	kPolyesterLiningCloth,
	kSilkShantung
};
WeavePattern* createWeavePreset(AtUInt32 preset);
