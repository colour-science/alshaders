#include "IES.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>

IESData::IESData(const std::string& fn)
: _xres(512), _yres(256), _valid(false), _data(NULL)
{
	float* _anglesV = NULL;
	float* _anglesH = NULL;
	float* _intensities = NULL;
	std::string line;
	std::ifstream in;
	in.open(fn.c_str());
	if (!in.is_open())
	{
		std::cerr << "[IES] Could not open file: " << fn << std::endl;
		return;
	}

	try
	{
		// fast-forward to TILT line
		while (line.substr(0, 4) != "TILT")
		{
			in >> line;
		}
		// fast-forward to the '1' character that marks the start of the bit we're interested in
		while (line != "1")
		{
			in >> line;
		}
	}
	catch (std::exception& e)
	{
		std::cerr << "[IES] Error reading file header " << fn << std::endl;
		return;
	}

	try
	{
		in >> line; // read lumens
		in >> line; // read candela scale
		in >> line; // read num vertical angles
		_numV = boost::lexical_cast<int>(line);
		_anglesV = new float[_numV];
		in >> line; // read num horizontal angles
		_numH = boost::lexical_cast<int>(line);
		_anglesH = new float[_numH];
		_intensities = new float[_numV*_numH];
		in >> line; // read another '1'
		in >> line; // read unit
		in >> line; // read width
		in >> line; // read length
		in >> line; // read height
		in >> line; // read a number
		in >> line; // another number
		in >> line; // ...and another number
		// now read the set of vertical angles
		for (int i=0; i < _numV; ++i)
		{
			in >> line;
			_anglesV[i] = boost::lexical_cast<float>(line);
		}

		// set our internal hemisphere flag
		if (_anglesV[0] == 0.0f)
		{
			if (_anglesV[_numV-1] == 90.0f)
			{
				_hemisphere = kBottom;
			}
			else
			{
				_hemisphere = kBoth;
			}
		}
		else
		{
			_hemisphere = kTop;
		}

		// convert to radians
		for (int i=0; i < _numV; ++i)
		{
			// we don't care which way it's pointing, since we'll be doing orientation ourselves
			if (_hemisphere == kTop)
			{
				_anglesV[i] -= 90.0f;
			}
			_anglesV[i] *= M_PI / 180.0f;
		}

		// now read the set of horizontal angles
		for (int i=0; i < _numH; ++i)
		{
			in >> line;
			_anglesH[i] = boost::lexical_cast<float>(line);
		}

		float lastAngle = _anglesH[_numH-1];
		if (lastAngle == 0.0f)
		{
			_symmetry = k0;
		}
		else if (lastAngle == 90.0f)
		{
			_symmetry = k90;
		}
		else if (lastAngle == 180.0f)
		{
			_symmetry = k180;
		}
		else
		{
			_symmetry = k360;
		}

		// convert to radians
		for (int i=0; i < _numH; ++i)
		{
			_anglesH[i] *= M_PI / 180.0f;
		}

		// now read the candela values
		float maxc = 0.0f;
		int n = _numV*_numH;
		for (int i=0; i < n; ++i)
		{
			in >> line;
			_intensities[i] = boost::lexical_cast<float>(line);
			if (_intensities[i] > maxc) maxc = _intensities[i];
		}
		// now normalize the candela values since we don't really care about physical units, just the shape...
		maxc = 1.0f / maxc;
		for (int i=0; i < n; ++i)
		{
			_intensities[i] *= maxc;
		}

		// rather than dicking around trying to find nearest angles at rendertime we'll precalculate a lookup table
		// figure out what the actual dimension of our array will be
		_hemisphere == kBoth ? _yend = _yres : _yend = _yres/2;
		_hemisphere == kBottom ? _thetaMax = M_PI/2 : _thetaMax = M_PI;
		switch (_symmetry)
		{
		case k0: _xend = 1; _phiMax = 0.0f; break;
		case k90: _xend = _xres/4; _phiMax = M_PI/2; break;
		case k180: _xend = _xres/2; _phiMax = M_PI; break;
		default: _xend = _xres; _phiMax = 2.0f * M_PI; break;
		}

		// allocate an array big enough to hold all our samples. We'll need to take account of symmetry when we index
		// into this during lookups
		_data = new float[_xend*_yend];
		memset(_data, 0, sizeof(float)*_xend*_yend);

		float phi = 0.0f; float phiStep; if (_phiMax==0.0f) phiStep = 0.0f; else phiStep = _phiMax / float(_xend);
		float theta = 0.0f, thetaStep = _thetaMax / float(_yend);

		// iterate over the intensities array and interpolate values from the normalized candela values read from the
		// file
		int lasth = 0, lastv = 0;
		int h, v;
		for (int x = 0; x < _xend;  ++x)
		{
			h = 0;
			while (_anglesH[h] < phi) ++h;
			theta = 0.0f;
			for (int y=0; y < _yend; ++y)
			{
				v = 0;
				while (_anglesV[v] < theta) ++v;
				int idx = y*_xend+x;
				if (v==0)
				{
					_data[idx] = _intensities[h*_numV];
				}
				else
				{
					float a0 = _anglesV[v-1];
					float a1 = _anglesV[v];
					float i0 = _intensities[h*_numV+v-1];
					float i1 = _intensities[h*_numV+v];
					float t = (theta-a0)/(a1-a0);
					_data[idx] = (1-t)*i0 + t*i1;
				}
				theta += thetaStep;
			}
			phi += phiStep;
		}

	}
	catch (std::exception& e)
	{
		std::cerr << "[IES] Error reading IES file " << fn << std::endl;
		in.close();
		delete[] _anglesH;
		delete[] _anglesV;
		delete[] _intensities;
		return;
	}

	// we're ok to access the data
	_valid = true;

	// clean up
	in.close();
	delete[] _anglesH;
	delete[] _anglesV;
	delete[] _intensities;
}

IESData::~IESData()
{
	delete[] _data;
}
