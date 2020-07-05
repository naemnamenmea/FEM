#pragma once

#include <string_view>
#include <stdexcept>

std::string_view& TrimFront(std::string_view& str, char delimiter = ' ');

std::string_view& TrimBack(std::string_view& str, char delimiter = ' ');

std::string_view& Trim(std::string_view& str, char delimiter = ' ');
