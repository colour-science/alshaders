#include "Photometric.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>
#include <ai.h>

#define _USE_MATH_DEFINES

PhotometricData::PhotometricData(const std::string& fn)
: _xres(1024), _yres(512), _valid(false), _data(NULL)
{
	std::ifstream in;
	in.open(fn.c_str());

	if (!in.is_open())
	{
		std::cerr << "[Photometric] Could not open file: " << fn << std::endl;
		return;
	}

	if (fn.substr(fn.size()-3, fn.size()) == "ies" ||
		fn.substr(fn.size()-3, fn.size()) == "IES")
	{
		readIES(in);
	}
	else if (fn.substr(fn.size()-3, fn.size()) == "ldt" ||
				fn.substr(fn.size()-3, fn.size()) == "LDT")
	{
		readLDT(in);
	}
	else
	{
		std::cerr << "[Photometric] unknown file extension '." <<  fn.substr(fn.size()-3, fn.size()) << "'" << std::endl;
	}

	in.close();
}

void PhotometricData::readIES(std::ifstream& in)
{
	float* _anglesV = NULL;
	float* _anglesH = NULL;
	float* _intensities = NULL;
	std::string line;

	char* array = new char[128000];
	try
	{
		// fast-forward to TILT line
		while (line.substr(0, 4) != "TILT")
		{
			in.getline(array, 128000);
			line = array;
		}
	}
	catch (std::exception& e)
	{
		std::cerr << "[IES] Error reading file header " << std::endl;
		delete[] array;
		return;
	}

	try
	{
		// the delimiting of the angle and candela values in any given ies file isn't particularly standardized
		// we'll just read a chunk of data and split it up afterwards
		in.read(array, 128000);
		if (in.bad())
		{
			std::cerr << "[IES] bad read data from " << std::endl;
			delete[] array;
			return;
		}
		if (!in.eof())
		{
			std::cerr << "[IES] failed to read whole file " << std::endl;
			delete[] array;
			return;
		}
		array[in.gcount()] = '\0';
		line = array;
		typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
		boost::char_separator<char> cs(", \n\t\r");
		tokenizer toks(line, cs);
		tokenizer::iterator it = toks.begin();


		++it; // read lumens
		++it; // read candela scale
		++it; // read num vertical angles
		_numV = boost::lexical_cast<int>(*it);
		++it; // read num horizontal angles
		_numH = boost::lexical_cast<int>(*it);
		++it; // read another '1'
		++it; // read unit
		++it; // read width
		++it; // read length
		++it; // read height
		++it; // read a number
		++it; // another number
		++it; // ...and another number
		++it; // ...

		_anglesV = new float[_numV];
		_anglesH = new float[_numH];
		_intensities = new float[_numV*_numH];


		// now read the set of vertical angles
		for (int i=0; i < _numV; ++i, ++it)
		{
			_anglesV[i] = boost::lexical_cast<float>(*it);
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
			_anglesV[i] *= AI_PI / 180.0f;
		}

		// now read the set of horizontal angles
		for (int i=0; i < _numH; ++i, ++it)
		{
			_anglesH[i] = boost::lexical_cast<float>(*it);
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
			_anglesH[i] *= AI_PI / 180.0f;
		}

		// now read the candela values
		float maxc = 0.0f;
		int n = _numV*_numH;
		for (int i=0; i < n; ++i, ++it)
		{
			_intensities[i] = boost::lexical_cast<float>(*it);
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
		_hemisphere == kBottom ? _thetaMax = AI_PIOVER2 : _thetaMax = AI_PI;
		switch (_symmetry)
		{
		case k0: _xend = 1; _phiMax = 0.0f; break;
		case k90: _xend = _xres/4; _phiMax = AI_PIOVER2; break;
		case k180: _xend = _xres/2; _phiMax = AI_PI; break;
		default: _xend = _xres; _phiMax = 2.0f * AI_PI; break;
		}

		// allocate an array big enough to hold all our samples. We'll need to take account of symmetry when we index
		// into this during lookups
		_data = new float[_xend*_yend];
		memset(_data, 0, sizeof(float)*_xend*_yend);

		float phi = 0.0f; float phiStep; if (_phiMax==0.0f) phiStep = 0.0f; else phiStep = _phiMax / float(_xend);
		float theta = 0.0f, thetaStep = _thetaMax / float(_yend);

		// iterate over the intensities array and interpolate values from the normalized candela values read from the
		// file into the data array we'll use to lookup aat render time.
		int lasth = 0, lastv = 0;
		int h, v;
		float av0, av1;
		float ah0, ah1;
		float i00, i01, i10, i11, i0, i1;
		int hm1, vm1;
		float tv, th;
		for (int x = 0; x < _xend;  ++x)
		{
			// find the pair of horinzontal angles neighbouring our lookup angle
			h = 0;
			while (_anglesH[h] < phi) ++h;
			if (h == 0)
			{
				hm1 = 0;
				th = 0.0f;
			}
			else
			{
				hm1 = h-1;
				th = (phi-_anglesH[hm1]) / (_anglesH[h]-_anglesH[hm1]);
			}

			theta = 0.0f;
			for (int y=0; y < _yend; ++y)
			{
				int idx = y*_xend+x;

				// find the pair of vertical angles neighbouring our lookup angle
				v = 0;
				while (_anglesV[v] < theta) ++v;
				if (v == 0)
				{
					vm1 = 0;
					tv = 0.0f;
				}
				else
				{
					vm1 = v-1;
					tv = (theta-_anglesV[vm1]) / (_anglesV[v]-_anglesV[vm1]);
				}

				// bilinearly interpolate the final lookup intensity value
				i00 = _intensities[hm1*_numV+vm1];
				i01 = _intensities[hm1*_numV+v];
				i10 = _intensities[h*_numV+vm1];
				i11 = _intensities[h*_numV+v];
				i0 = (1-tv)*i00 + tv*i01;
				i1 = (1-tv)*i10 + tv*i11;

				_data[idx] = (1-th)*i0 + th*i1;

				theta += thetaStep;
			}

			phi += phiStep;
		}
	}
	catch (std::exception& e)
	{
		std::cerr << "[IES] Error reading IES file " << " " << e.what() <<  std::endl;
		delete[] _anglesH;
		delete[] _anglesV;
		delete[] _intensities;
		delete[] array;
		return;
	}

	// we're ok to access the data
	_valid = true;

	// clean up
	delete[] _anglesH;
	delete[] _anglesV;
	delete[] _intensities;
	delete[] array;
}

void PhotometricData::readLDT(std::ifstream& in)
{
	// TODO: implement LDT file reading
}

PhotometricData::~PhotometricData()
{
	delete[] _data;
}
