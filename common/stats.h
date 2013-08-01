#include <iostream>

// stats.h
struct Range
{
	Range(const std::string& nm, bool rod=false)
	: name(nm), min(HUGE), max(-HUGE), total(0.0), n(0.0), reportOnDestruction(rod)
	{}

	~Range()
	{
		report(std::cerr);
	}

	void addSample(double x)
	{
		min = std::min(x, min);
		max = std::max(x, max);
		total += x;
		n++;
	}

	void report(std::ostream& os)
	{
		os << "[" << name << "] " << "min: " << min << " max: " << max << " avg: " << total/n << std::endl;
	}

	std::string name;

	double min;
	double max;
	double total;
	double n;
	bool reportOnDestruction;
};