#include "GaussIntegr.hpp"

double LegendrePoly(int ord_, double x_)  // ord_ >= 0
{
	if (ord_ == 0)
		return 1.;

	if (ord_ == 1)
		return x_;

	return ((2. * ord_ - 1.) * x_ * LegendrePoly(ord_ - 1, x_) -
			LegendrePoly(ord_ - 2, x_) * (ord_ - 1.)) /
		   ord_;
}

double LegendrePolyDiff(int ord_, double x)
{
	return (LegendrePoly(ord_ - 1, x) - x * LegendrePoly(ord_, x)) * ord_ / (1. - x * x);
}

double GaussX(int ord_, int point_num_, int newton_steps_)
{
	double x = -cos(acos(-1.) * (2. * (point_num_ + 1.) - 0.5) / (2. * ord_ + 1.));
	for (int i = 0; i < newton_steps_; ++i)
	{
		x -= LegendrePoly(ord_, x) / LegendrePolyDiff(ord_, x);
	}
	return x;
}

double GaussW(int ord_, double x_)
{
	return 2. / ((1. - x_ * x_) * pow(LegendrePolyDiff(ord_, x_), 2));
}
