#pragma once

#include <string>
#include <vector>
#include <utility>
#include <algorithm>

#include "NumTypes.h"

/*
  template <typename T> struct fLessName
  {
   bool operator()(std::pair<std::string,T> pr1, std::pair<std::string,T> pr2) const
           {return pr1.first < pr2.first;}
  };
*/

template <typename T> struct fTheSame2nd
{
  fTheSame2nd(T p) : m_second(p) {}
  bool operator()(std::pair<std::string, T> pr) const
  {
    return pr.second == m_second;
  }

  T m_second;
};

struct fTheSame1st
{
  fTheSame1st(const std::string& str) : m_first(str) {}
  template <typename T>
  bool operator()(std::pair<std::string, T> pr) const
  {
    return pr.first == m_first;
  }

  const std::string& m_first;
};

template <typename T>
class tNamesToFuns
{
private:
  typedef T* (*tPtrToFun)();
  std::vector<std::pair<std::string, tPtrToFun> > Array;

public:
  tNamesToFuns(const char* names_[], const tPtrToFun functions_[], cardinal_t dim_) : Array()
  {
    Array.reserve(dim_);
    for (cardinal_t i = 0; i < dim_; ++i)
      Array.push_back(std::make_pair(std::string(names_[i]), functions_[i]));
    //      std::sort(Array.begin(),Array.end(),fLessName<tPtrToFun>());
  }
  const std::string& FindName(tPtrToFun p) const { return std::find_if(Array.begin(), Array.end(), fTheSame2nd<tPtrToFun>(p))->first; }
  T* CallFunction(const std::string& name_) const { return std::find_if(Array.begin(), Array.end(), fTheSame1st(name_))->second(); }
};
