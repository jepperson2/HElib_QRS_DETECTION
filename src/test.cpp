#include <fstream>
#include <sstream>
#include <iostream>
#include <bitset>
#include "helper_functions.h"

using namespace std;

int main(int argc, char **argv)
{
	Conversion conv;
	
	long z = 5821253928;
//	bool x = true;
	string y = conv.long2Str(z);
	cout << y << endl;
    return 0;		
}


