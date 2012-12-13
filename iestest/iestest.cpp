#include "IES.h"
#include <iostream>

int main(int argc, char** argv)
{
	if (argc != 2)
	{
		std::cout << "\n\tusage: iestest <ies_file.ies>" << std::endl;
		return -1;
	}

	IESData iesdata(argv[1]);
}
