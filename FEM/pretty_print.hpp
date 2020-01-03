#pragma once

#include <algorithm>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

template <class Container>
string Join(const Container& container, char delimiter);

template <typename T>
ostream& operator<<(ostream& os, const vector<T>& v)
{
  return os << '[' << Join(v, ',') << ']';
}

template <typename K, typename V>
ostream& operator<<(ostream& os, const map<K, V>& m)
{
  return os << '{' << Join(m, ',') << '}';
}

template <typename K, typename V>
ostream& operator<<(ostream& os, const unordered_multimap<K, V>& m)
{
  return os << '{' << Join(m, ',') << '}';
}

template <typename L, typename R>
ostream& operator<<(ostream& os, const pair<L, R>& p)
{
  return os << '(' << p.first << ',' << p.second << ')';
}

template <typename T>
ostream& operator<<(ostream& os, const set<T>& s)
{
  return os << '{' << Join(s, ',') << '}';
}

template <typename T, class Comp>
ostream& operator<<(ostream& os, const multiset<T, Comp>& s)
{
  return os << '{' << Join(s, ',') << '}';
}

template <typename T>
ostream& operator<<(ostream& os, const list<T>& s)
{
  return os << '(' << Join(s, ',') << ')';
}

template <typename Collection>
string Join(const Collection& c, char d)
{
  stringstream ss;

  bool flag = false;
  for (const auto& el : c)
  {
    if (flag)
    {
      ss << d;
    }
    flag = true;
    ss << el;
  }
  return ss.str();
}
