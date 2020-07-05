#pragma once

#include <iostream>

template <class T>
class point3d
{
public:
	typedef T ValueType;
	constexpr static unsigned int DIM = 3;

	point3d()
	{
		m_elems[0] = T{};
		m_elems[1] = T{};
		m_elems[2] = T{};
	}

	point3d(T value)
	{
		m_elems[0] = value;
		m_elems[1] = value;
		m_elems[2] = value;
	}

	point3d(T x, T y, T z)
	{
		m_elems[0] = x;
		m_elems[1] = y;
		m_elems[2] = z;
	}

	const T x() const
	{
		return m_elems[0];
	}

	const T y() const
	{
		return m_elems[1];
	}

	const T z() const
	{
		return m_elems[2];
	}

	T& x()
	{
		return m_elems[0];
	}

	T& y()
	{
		return m_elems[1];
	}

	T& z()
	{
		return m_elems[2];
	}

	operator T*()
	{
		return m_elems;
	}

	operator T const*() const
	{
		return m_elems;
	}

	const point3d<T> Normilized() const
	{
		return point3d<T>(*this).Normilize();
	}

	const point3d<T>& Normilize()
	{
		(*this) /= ~(*this);
		return *this;
	}

	T operator*(const point3d<T>& p) const
	{
		return this->x() * p.x() + this->y() * p.y() + this->z() * p.z();
	}

	point3d<T> operator-(const point3d<T>& p) const
	{
		return point3d<T>(this->x() - p.x(), this->y() - p.y(), this->z() - p.z());
	}

	point3d<T> operator+(const point3d<T>& p) const
	{
		return point3d<T>(this->x() + p.x(), this->y() + p.y(), this->z() + p.z());
	}

	point3d<T> operator/(const double factor) const
	{
		return point3d<T>(this->x() / factor, this->y() / factor, this->z() / factor);
	}

	void operator/=(const double factor)
	{
		m_elems[0] /= factor;
		m_elems[1] /= factor;
		m_elems[2] /= factor;
	}

	point3d<T> operator*(const double factor) const
	{
		return point3d<T>(this->x() * factor, this->y() * factor, this->z() * factor);
	}

	bool operator==(const point3d<T>& p) const
	{
		return this->x() == p.x() && this->y() == p.y() && this->z() == p.z();
	}

	point3d<T> operator-() const
	{
		return point3d<T>(-this->x(), -this->y(), -this->z());
	}

	void operator+=(const point3d<T>& p)
	{
		this->x() += p.x();
		this->y() += p.y();
		this->z() += p.z();
	}

	double operator~() const
	{
		return sqrt((*this) * (*this));
	}

	void operator*=(const double factor)
	{
		this->x() *= factor;
		this->y() *= factor;
		this->z() *= factor;
	}

	point3d<T> operator-(const T v) const
	{
		return point3d<T>(this->x() - v, this->y() - v, this->z() - v);
	}
	point3d<T> operator+(const T v) const
	{
		return point3d<T>(this->x() + v, this->y() + v, this->z() + v);
	}

	bool operator!=(const point3d<T>& p1) const
	{
		return !(*this == p1);
	}

	point3d<T> operator%(const point3d<T>& p) const
	{
		return point3d<T>(
			this->y() * p.z() - this->z() * p.y(),
			-(this->x() * p.z() - this->z() * p.x()),
			this->x() * p.y() - this->y() * p.x());
	}

private:
	T m_elems[3];
};

template <class T>
std::istream& operator>>(std::istream& is, point3d<T>& p)
{
	is >> p.x() >> p.y() >> p.z();
	return is;
}

template <class T>
std::ostream& operator<<(std::ostream& os, const point3d<T>& p)
{
	os << p.x() << ' ' << p.y() << ' ' << p.z();
	return os;
}
