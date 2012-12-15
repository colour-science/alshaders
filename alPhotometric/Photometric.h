#pragma once

#include <string>
#include <cmath>
#include <iostream>

class PhotometricData
{
public:
	/// Creates a new PhotometricData object from the given photometric data file.
	/// Currently only IES files are supported
	PhotometricData(const std::string& fn);
	~PhotometricData();

	/// Return the normalized intensity of the light in the given spherical direction
	inline float lookup(float phi, float theta)
	{
		if (!_valid) return 1.0f;
		return _data[index(phi, theta)];
	}

	/// Is the data valid? If this returns false,
	inline bool isValid()
	{
		return _valid;
	}



private:
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
			// TODO: this is probably the worst way of handling symmetry ever... but it works.
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


	void readIES(std::ifstream& in);
	void readLDT(std::ifstream& in);

};
