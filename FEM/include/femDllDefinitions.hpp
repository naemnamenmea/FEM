#pragma once

#ifdef _WIN32
#define FEM_API __declspec(dllexport)
#elif defined(linux)
#define FEM_API __attribute__((visibility("default")))
#endif	// _WIN32
