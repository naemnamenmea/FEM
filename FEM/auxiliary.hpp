#pragma once

#include <string_view>

void TrimFront(std::string_view& sv, char symbol = ' ')
{
  while (!sv.empty() && sv.front() == symbol)
  {
    sv.remove_prefix(1);
  }
}

void TrimBack(std::string_view& sv, char symbol = ' ')
{
  while (!sv.empty() && sv.back() == symbol)
  {
    sv.remove_suffix(1);
  }
}

void Trim(std::string_view& sv, char symbol = ' ')
{
  TrimFront(sv, symbol);
  TrimBack(sv, symbol);
}
