#pragma once

#include <iostream>

template <class T>
class point3d
{
 public:
  typedef T ValueType;

  point3d()
  {
    elems[0] = T{};
    elems[1] = T{};
    elems[2] = T{};
  }

  point3d(T x, T y, T z)
  {
    elems[0] = x;
    elems[1] = y;
    elems[2] = z;
  }

  const T x() const { return elems[0]; }

  const T y() const { return elems[1]; }

  const T z() const { return elems[2]; }

  T& x() { return elems[0]; }

  T& y() { return elems[1]; }

  T& z() { return elems[2]; }

  inline operator T*() { return elems; }

  inline operator T const*() const { return elems; }

  inline T operator*(const point3d<T>& p) const
  {
    return this->x() * p.x() + this->y() * p.y() + this->z() * p.z();
  }

  inline point3d<T> operator-(const point3d<T>& p) const
  {
    return point3d<T>(this->x() - p.x(), this->y() - p.y(), this->z() - p.z());
  }

  inline point3d<T> operator+(const point3d<T>& p) const
  {
    return point3d<T>(this->x() + p.x(), this->y() + p.y(), this->z() + p.z());
  }

  inline point3d<T> operator/(const double factor) const
  {
    return point3d<T>(this->x() / factor, this->y() / factor,
                      this->z() / factor);
  }

  inline point3d<T> operator*(const double factor) const
  {
    return point3d<T>(this->x() * factor, this->y() * factor,
                      this->z() * factor);
  }

  inline bool operator==(const point3d<T>& p) const
  {
    return this->x() == p.x() && this->y() == p.y() && this->z() == p.z();
  }

  inline point3d<T> operator-() const
  {
    return point3d<T>(-this->x(), -this->y(), -this->z());
  }

  inline void operator+=(const point3d<T>& p)
  {
    this->x() += p.x();
    this->y() += p.y();
    this->z() += p.z();
  }

  inline double operator~() const { return sqrt((*this) * (*this)); }

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

  bool operator!=(const point3d<T>& p1) const { return !(*this == p1); }

  point3d<T> operator%(const point3d<T>& p) const
  {
    return point3d<T>(this->y() * p.z() - this->z() * p.y(),
                      -(this->x() * p.z() - this->z() * p.x()),
                      this->x() * p.y() - this->y() * p.x());
  }

 private:
  T elems[3];
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
