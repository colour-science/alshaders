#pragma once
#include <ai.h>
#include <cstring>
#include <string>


#define SPD_NUM_SAMPLES 81
#define SPD_INCREMENT 5
#define SPD_MIN 380
#define SPD_MAX 780

class SPD
{
public:
	enum Constant
	{
		kBLACK=0,
		kWHITE,
		kD60,
		kD65,

		kDARK_SKIN,
		kLIGHT_SKIN,
		kBLUE_SKY,
		kFOLIAGE,
		kBLUE_FLOWER,
		kBLUISH_GREEN,
		kORANGE,
		kPURPLISH_BLUE,
		kMODERATE_RED,
		kPURPLE,
		kYELLOW_GREEN,
		kORANGE_YELLOW,
		kBLUE,
		kGREEN,
		kRED,
		kYELLOW,
		kMAGENTA,
		kCYAN,
		kWHITE_95,
		kNEUTRAL_8,
		kNEUTRAL_65,
		kNEUTRAL_5,
		kNEUTRAL_35,
		kBLACK_2
	};

	enum Nk
	{
		kGOLD=0,
		kSILVER,
		kCOPPER,
		kTUNGSTEN
	};

	SPD();
	SPD(int c);

	void set(float f);
	void set(int c);
	const float* get(int c) const;

	AtRGB xyz(Constant t=kWHITE) const;

	SPD operator +(const SPD& s) const;
	SPD operator -(const SPD& s) const;
	SPD operator *(const SPD& s) const;
	SPD operator *(float s) const;

	inline static void fresnel(float cos_theta, int nk, SPD& spd)
	{
		const float* nkdata = NULL;
		switch (nk)
		{
		case kGOLD:
			nkdata = SPD::_Au_nk;
			break;
		case kSILVER:
			nkdata = SPD::_Ag_nk;
			break;
		case kCOPPER:
			nkdata = SPD::_Cu_nk;
			break;
		case kTUNGSTEN:
			nkdata = SPD::_W_nk;
			break;
		default:
			break;
		}

		if (!nkdata) 
		{
			spd.set(kBLACK);
			return;
		}

		float mx = 0.0f;
		for (int i=0; i < SPD_NUM_SAMPLES; ++i)
		{
			spd._spd[i] = frcond(cos_theta, nkdata[i*2], nkdata[i*2+1]);
			mx = MAX(mx, spd._spd[i]);
		}

		mx = (1.0f / mx) * 0.99f;

		for (int i=0; i < SPD_NUM_SAMPLES; ++i)
		{
			spd._spd[i] *= mx;
		}
	}

private:

	inline static float frcond(float cosi, const float eta, const float k)
	{
	    float tmp = (eta*eta + k*k) * cosi*cosi;
	    float Rparl2 = (tmp - (2.f * eta * cosi) + 1) /
	                      (tmp + (2.f * eta * cosi) + 1);
	    float tmp_f = eta*eta + k*k;
	    float Rperp2 =
	        (tmp_f - (2.f * eta * cosi) + cosi*cosi) /
	        (tmp_f + (2.f * eta * cosi) + cosi*cosi);
	    return (Rparl2 + Rperp2) / 2.f;
	}

	// Samples from 380nm - 780nm at 5nm increments
	float _spd[SPD_NUM_SAMPLES];
	const static float _cieMatch[SPD_NUM_SAMPLES][3];
	const static float _D65[SPD_NUM_SAMPLES];
	const static float _D60[SPD_NUM_SAMPLES];
	const static float _black[SPD_NUM_SAMPLES];
	const static float _white[SPD_NUM_SAMPLES];
	const static float _bluish_green[SPD_NUM_SAMPLES];
	const static float _yellow[SPD_NUM_SAMPLES];
	const static float _purplish_blue[SPD_NUM_SAMPLES];
	const static float _neutral_8[SPD_NUM_SAMPLES];
	const static float _black_2[SPD_NUM_SAMPLES];
	const static float _magenta[SPD_NUM_SAMPLES];
	const static float _neutral_5[SPD_NUM_SAMPLES];
	const static float _blue_flower[SPD_NUM_SAMPLES];
	const static float _blue[SPD_NUM_SAMPLES];
	const static float _purple[SPD_NUM_SAMPLES];
	const static float _orange[SPD_NUM_SAMPLES];
	const static float _green[SPD_NUM_SAMPLES];
	const static float _red[SPD_NUM_SAMPLES];
	const static float _neutral_35[SPD_NUM_SAMPLES];
	const static float _light_skin[SPD_NUM_SAMPLES];
	const static float _yellow_green[SPD_NUM_SAMPLES];
	const static float _cyan[SPD_NUM_SAMPLES];
	const static float _white_95[SPD_NUM_SAMPLES];
	const static float _moderate_red[SPD_NUM_SAMPLES];
	const static float _blue_sky[SPD_NUM_SAMPLES];
	const static float _orange_yellow[SPD_NUM_SAMPLES];
	const static float _dark_skin[SPD_NUM_SAMPLES];
	const static float _foliage[SPD_NUM_SAMPLES];
	const static float _neutral_65[SPD_NUM_SAMPLES];

	const static float _Au_nk[SPD_NUM_SAMPLES*2];
	const static float _Cu_nk[SPD_NUM_SAMPLES*2];
	const static float _W_nk[SPD_NUM_SAMPLES*2];
	const static float _Ag_nk[SPD_NUM_SAMPLES*2];
	const static float _Os_nk[SPD_NUM_SAMPLES*2];
	const static float _Stsg0_nk[SPD_NUM_SAMPLES*2];
};

struct ColorSpace
{
	ColorSpace(const std::string& name_, const AtPoint2& r, const AtPoint2& g, const AtPoint2& b, const AtPoint2& w):
		name(name_), red(r), green(g), blue(b), white(w), matrix(false)
	{}

	ColorSpace(float m00, float m10, float m20, float m01, float m11, float m21, float m02, float m12, float m22);

	AtRGB xyzToRgb(const AtRGB& xyz);

	std::string name;	/// name of the system, e.g. rec709, ACES etc
	AtPoint2 red;		/// chromaticity of the red primary
	AtPoint2 green;		/// chromaticity of the green primary
	AtPoint2 blue;		/// chromticity of the blue primary
	AtPoint2 white;		/// chromaticity of the white point
	float m[9];
	bool matrix;
};

extern AtPoint2 xyD60;
extern AtPoint2 xyACES;
extern AtPoint2 xyP3;
extern ColorSpace CsRec709;
extern ColorSpace CsRec2020;
extern ColorSpace CsACES;
extern ColorSpace CsACES_D65;
extern ColorSpace CsP3;
extern ColorSpace CsP3_D65;