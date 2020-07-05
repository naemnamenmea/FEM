#include "numeric_math.hpp"

double round(double val, int precision)
{
	double mult = pow(10, precision);
	val *= mult;
	return (val < 0 ? ceil(val - 0.5) : floor(val + 0.5)) / mult;
}
