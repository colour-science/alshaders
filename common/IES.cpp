#include "IES.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <boost/lexical_cast.hpp>

IESData::IESData(const std::string& fn)
: _xres(1024), _yres(512), _anglesH(NULL), _anglesV(NULL), _intensities(NULL), _valid(false)
{
	std::string line;
	std::ifstream in;
	in.open(fn.c_str());
	if (!in.is_open())
	{
		std::cerr << "Could not open IES file: " << fn << std::endl;
		return;
	}

	while (line.substr(0, 4) != "TILT")
	{
		in >> line;
		//std::cout << line << std::endl;
	}
	//std::cout << "FOUND TILT" << std::endl;
	while (line != "1")
	{
		in >> line;
		//std::cout << line << std::endl;
	}
	//std::cout << "FOUND 1" << std::endl;

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
				std::cout << "Bottom" << std::endl;
			}
			else
			{
				_hemisphere = kBoth;
				std::cout << "Both" << std::endl;
			}
		}
		else
		{
			_hemisphere = kTop;
			std::cout << "Top" << std::endl;
		}

		// convert to radians
		for (int i=0; i < _numV; ++i)
		{
			_anglesV[i] *= M_PI / 180.0f;
			// we don't care which way it's pointing, since we'll be doing orientation ourselves
			if (_hemisphere == kTop)
			{
				_anglesV[i] -= M_PI / 2.0f;
			}
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

		// convert t radians
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
		// now normalize the candela values since we don't really care...
		for (int i=0; i < n; ++i)
		{
			_intensities[i] /= maxc;
		}

		// rather than dicking around trying to find nearest angles at rendertime we'll precalculate a lookup table
		// at 512x256 this implies an overhead of ~500k per IES file. This is extremely inefficient for compact files
		// but it will be fast...

		float phi = 0.0f, phiStep = (2.0f * M_PI) / float(_xres);
		float theta = 0.0f, thetaStep = M_PI / float(_yres);

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

		_data = new float[_xend*_yend];
		memset(_data, 0, sizeof(float)*_xend*_yend);

		int lasth = 0, lastv = 0;
		int h, v;
		for (int y=0; y < _yend; theta+=thetaStep, y++)
		{
			//v = lastv;
			v = 0;
			while (theta > _anglesV[v])
			{
				++v;
			}
			for (int x=0; x < _xend; phi += phiStep, x++)
			{
				// find horizontal angle that is <= phi
				//h = lasth;
				h = 0;
				while (phi > _anglesH[h])
				{
					++h;
				}
				int i_idx = h * _numV + v;
				int idx = index(phi, theta);
				_data[idx] = _intensities[i_idx];
			}
		}
	}
	catch (std::exception& e)
	{
		std::cerr << "Error reading IES file " << fn << std::endl;
		delete[] _anglesH;
		delete[] _anglesV;
		delete[] _intensities;
		return;
	}

	// we're ok to access the data
	std::cerr << "OK TO READ BITCHES" << std::endl;

	_valid = true;
}
