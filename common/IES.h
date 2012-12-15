#pragma once

#include <string>
#include <cmath>
#include <iostream>

class IESData
{
public:
	IESData(const std::string& fn);
	~IESData();

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
		int y = static_cast<int>(floorf(theta / _thetaMax * float(_yend)));
		int x;
		if (_symmetry==k0)
		{
			x = 0;
		}
		else
		{

			if (_symmetry==k90)
			{
				if (phi > M_PI) phi -= M_PI;
				if (phi > M_PI/2.0f && phi < M_PI)
				{
					phi -= M_PI/2.0f;
					phi = (M_PI/2.0f) - phi;
				}
			}
			else if (_symmetry == k180)
			{
				if (phi > M_PI)
				{
					phi -= M_PI;
					phi = M_PI - phi;
				}
			}

			phi = fmodf(phi, _phiMax);
			x = static_cast<int>(floorf(phi / _phiMax * float(_xend)));
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
	float* _data;
	Hemisphere _hemisphere;
	Symmetry _symmetry;

	bool _valid;
};
