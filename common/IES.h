#pragma once

#include <string>
#include <cmath>
#include <iostream>

class IESData
{
public:
	IESData(const std::string& fn);

	enum Hemisphere
	{
		kBottom=0,
		kTop,
		kBoth
	};

	enum Symmetry
	{
		k0,
		k90,
		k180,
		k360
	};

	inline int index(float phi, float theta)
	{
		theta = fmodf(theta, _thetaMax);
		int y = static_cast<int>(floorf(theta / _thetaMax * _yend));
		int x;
		if (_symmetry==k0)
		{
			x = 0;
		}
		else
		{
			phi = fmodf(phi, _phiMax);
			x = static_cast<int>(floorf(phi / _phiMax * _xend));
		}

		return y*_xend+x;
	}

	inline float lookup(float phi, float theta)
	{
		if (!_valid) return 1.0f;
		return _data[index(phi, theta)];
	}

	inline bool isValid()
	{
		return _valid;
	}

private:
	int _xres;
	int _yres;
	int _xend;
	int _yend;
	float _phiMax;
	float _thetaMax;
	int _numV;
	int _numH;
	float* _anglesV;
	float* _anglesH;
	float* _intensities;
	float* _data;
	Hemisphere _hemisphere;
	Symmetry _symmetry;

	bool _valid;
};
