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

private:
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
};

struct ColorSpace
{
	ColorSpace(const std::string& name_, const AtPoint2& r, const AtPoint2& g, const AtPoint2& b, const AtPoint2& w):
		name(name_), red(r), green(g), blue(b), white(w)
	{}

	AtRGB xyzToRgb(const AtRGB& xyz);

	std::string name;	/// name of the system, e.g. rec709, ACES etc
	AtPoint2 red;		/// chromaticity of the red primary
	AtPoint2 green;		/// chromaticity of the green primary
	AtPoint2 blue;		/// chromticity of the blue primary
	AtPoint2 white;		/// chromaticity of the white point
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