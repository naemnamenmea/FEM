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
		return m_integr.begin();
	}

	const auto end() const
	{
		return m_integr.end();
	}

	template <typename Targ, typename Tret, typename FUNC>
	Tret Calculate(FUNC f_) const
	{
		return m_integr.ByPlus<Targ, Tret, FUNC>(f_);
	}

	template <typename Targ, typename Tret, typename FUNC>
	Tret& Calculate(FUNC f_, Tret& result_) const
	{
		return m_integr.ByPlusAssgn<Targ, Tret, FUNC>(f_, result_);
	}

	template <typename Tret, typename FUNC>
	Tret& ByPlus_ArrArg(FUNC f_, Tret& result_) const
	{
		return m_integr.ByPlusAssgn_ArrArg<Tret, FUNC>(f_, result_);
	}

	template <typename Tret, typename FUNC>
	Tret ByPlus_ArrArg(FUNC f_) const
	{
		return m_integr.ByPlus_ArrArg<Tret, FUNC>(f_);
	}

private:
	fIntegrate<Dim, Ord> m_integr;
};
}  // namespace GaussIntegr
