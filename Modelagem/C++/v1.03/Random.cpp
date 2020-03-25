#include "Random.h"

#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <math.h>
#include <functional>
#include <limits>
#include <algorithm>

double Random::NumeroRandomico(double Min,double Max)
{
	double RandFactor = 1.0/(double(RAND_MAX)+1.0);
    return rand() * RandFactor* (Max - Min) + Min;
}
