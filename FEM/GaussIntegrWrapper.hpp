#pragma once
#include "GaussIntegr.hpp"

namespace GaussIntegr
{
template <int Dim, int Ord>
class GaussIntegrWrapper
{
public:
	const auto begin() const
	{
		return integr.begin();
	}

	const auto end() const
	{
		return integr.end();
	}

	template <typename Targ, typename Tret, typename FUNC>
	Tret ByPlus(FUNC f_) const
	{
		return integr.ByPlus<Targ, Tret, FUNC>(f_);
	}

	template <typename Targ, typename Tret, typename FUNC>
	Tret& ByPlus(FUNC f_, Tret& result_) const
	{
		return integr.ByPlusAssgn<Targ, Tret, FUNC>(f_, result_);
	}

	template <typename Tret, typename FUNC>
	Tret& ByPlus_ArrArg(FUNC f_, Tret& result_) const
	{
		return integr.ByPlusAssgn_ArrArg<Tret, FUNC>(f_, result_);
	}

	template <typename Tret, typename FUNC>
	Tret ByPlus_ArrArg(FUNC f_) const
	{
		return integr.ByPlus_ArrArg<Tret, FUNC>(f_);
	}

private:
	fIntegrate<Dim, Ord> integr;
};
}  // namespace GaussIntegr
