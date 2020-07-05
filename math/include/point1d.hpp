#pragma once

#include <iostream>

template <class T>
class point1d
{
public:
	typedef T ValueType;
	constexpr static unsigned int DIM = 1;

	point1d()
	{
		m_elems[0] = T{};
	}

	point1d(T x)
	{
		m_elems[0] = x;
	}

	const T x() const
	{
		return m_elems[0];
	}

	T& x()
	{
		return m_elems[0];
	}

	operator T*()
	{
		return m_elems;
	}

	operator T const*() const
	{
		return m_elems;
	}

	const point1d<T> Normilized() const
	{
		return point1d<T>(*this).Normilize();
	}

	const point1d<T>& Normilize()
	{
		(*this) /= ~(*this);
		return *this;
	}

	T operator*(const point1d<T>& p) const
	{
		return this->x() * p.x();
	}

	point1d<T> operator-(const point1d<T>& p) const
	{
		return point1d<T>(this->x() - p.x());
	}

	point1d<T> operator+(const point1d<T>& p) const
	{
		return point1d<T>(this->x() + p.x());
	}

	point1d<T> operator/(const double factor) const
	{
		return point1d<T>(this->x() / factor);
	}

	void operator/=(const double factor)
	{
		m_elems[0] /= factor;
	}

	point1d<T> operator*(const double factor) const
	{
		return point1d<T>(this->x() * factor);
	}

	bool operator==(const point1d<T>& p) const
	{
		return this->x() == p.x();
	}

	point1d<T> operator-() const
	{
		return point1d<T>(-this->x());
	}

	void operator+=(const point1d<T>& p)
	{
		this->x() += p.x();
	}

	double operator~() const
	{
		return sqrt((*this) * (*this));
	}

	void operator*=(const double factor)
	{
		this->x() *= factor;
	}

	point1d<T> operator-(const T v) const
	{
		return point1d<T>(this->x() - v);
	}
	point1d<T> operator+(const T v) const
	{
		return point1d<T>(this->x() + v);
	}

	bool operator!=(const point1d<T>& p1) const
	{
		return !(*this == p1);
	}

private:
	T m_elems[1];
};

template class point1d<double>;

template <class T>
std::istream& operator>>(std::istream& is, point1d<T>& p)
{
	is >> p.x();
	return is;
}

template <class T>
std::ostream& operator<<(std::ostream& os, const point1d<T>& p)
{
	os << p.x();
	return os;
}
