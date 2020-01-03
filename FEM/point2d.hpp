#pragma once

#include <iostream>

template <class T>
class point2d
{
 public:
  typedef T ValueType;

  point2d()
  {
    elems[0] = T{};
    elems[1] = T{};
  }

  point2d(T x, T y)
  {
    elems[0] = x;
    elems[1] = y;
  }

  const T x() const { return elems[0]; }

  const T y() const { return elems[1]; }

  T& x() { return elems[0]; }

  T& y() { return elems[1]; }

  inline operator T*() { return elems; }

  inline operator T const*() const { return elems; }

  inline T operator*(const point2d<T>& p) const
  {
    return this->x() * p.x() + this->y() * p.y();
  }

  inline point2d<T> operator-(const point2d<T>& p) const
  {
    return point2d<T>(this->x() - p.x(), this->y() - p.y());
  }

  inline point2d<T> operator+(const point2d<T>& p) const
  {
    return point2d<T>(this->x() + p.x(), this->y() + p.y());
  }

  inline point2d<T> operator/(const double factor) const
  {
    return point2d<T>(this->x() / factor, this->y() / factor);
  }

  inline point2d<T> operator*(const double factor) const
  {
    return point2d<T>(this->x() * factor, this->y() * factor);
  }

  inline bool operator==(const point2d<T>& p) const
  {
    return this->x() == p.x() && this->y() == p.y();
  }

  inline point2d<T> operator-() const
  {
    return point2d<T>(-this->x(), -this->y());
  }

  inline void operator+=(const point2d<T>& p)
  {
    this->x() += p.x();
    this->y() += p.y();
  }

  inline double operator~() const { return sqrt((*this) * (*this)); }

  void operator*=(const double factor)
  {
    this->x() *= factor;
    this->y() *= factor;
  }

  point2d<T> operator-(const T v) const
  {
    return point2d<T>(this->x() - v, this->y() - v);
  }
  point2d<T> operator+(const T v) const
  {
    return point2d<T>(this->x() + v, this->y() + v);
  }

  bool operator!=(const point2d<T>& p1) const { return !(*this == p1); }

 private:
  T elems[2];
};

template class point2d<double>;

template <class T>
std::istream& operator>>(std::istream& is, point2d<T>& p)
{
  is >> p.x() >> p.y();
  return is;
}

template <class T>
std::ostream& operator<<(std::ostream& os, const point2d<T>& p)
{
  os << p.x() << ' ' << p.y();
  return os;
}