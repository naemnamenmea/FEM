// Author: Chekhov Vladimir Valerevich

#include "stdafx.hpp"

#include "GaussIntegr.hpp"

namespace GaussIntegr
{
real_t LegendrePoly(int ord_, real_t x_)  // ord_>=0
{
	return ord_ == 0 ? 1.
					 : ord_ == 1 ? x_
								 : ((2. * ord_ - 1.) * x_ * LegendrePoly(ord_ - 1, x_) -
									LegendrePoly(ord_ - 2, x_) * (ord_ - 1.)) /
									   ord_;
}

real_t LegendrePolyDiff(int ord_, real_t x)
{
	return (LegendrePoly(ord_ - 1, x) - x * LegendrePoly(ord_, x)) * ord_ / (1. - x * x);
}

real_t GaussX(int ord_, int point_num_, int newton_steps_)
{
	// if (ord_==1)  return 0.;
	real_t x = -cos(acos(-1.) * (2. * (point_num_ + 1.) - 0.5) / (2. * ord_ + 1.));
	for (int i(0); i < newton_steps_; ++i) x -= LegendrePoly(ord_, x) / LegendrePolyDiff(ord_, x);
	return x;
}

real_t GaussW(int ord_, real_t x_)
{
	// if (ord_ == 1)  return 2.;
	return 2. / ((1. - x_ * x_) * pow(LegendrePolyDiff(ord_, x_), 2));
}
}  // namespace GaussIntegr
