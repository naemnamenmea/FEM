#pragma once

#include <chrono>
#include <iostream>
#include <sstream>
#include <string>

using namespace std::chrono;

template <typename TOStream = std::ostream>
class LogDuration
{
public:
	explicit LogDuration(const std::string& msg = "", TOStream& oStream = std::cerr)
		: m_message(msg + ": "), m_OStream(oStream), m_start(steady_clock::now())
	{
	}

	~LogDuration()
	{
		auto finish = steady_clock::now();
		auto dur = finish - m_start;
		m_OStream << m_message << duration_cast<milliseconds>(dur).count() << " ms" << std::endl;
	}

private:
	std::string m_message;
	steady_clock::time_point m_start;
	TOStream& m_OStream;
};

#define UNIQ_ID_IMPL(lineno) _a_local_var_##lineno
#define UNIQ_ID(lineno) UNIQ_ID_IMPL(lineno)

#define LOG_DURATION(message) LogDuration UNIQ_ID(__LINE__){message};
#define LOG_DURATION_STREAM(stream, message) LogDuration UNIQ_ID(__LINE__){message, stream};

struct TotalDuration
{
	explicit TotalDuration(const std::string& msg = "") : message(msg + ": "), value(0)
	{
	}
	~TotalDuration()
	{
		std::ostringstream os;
		os << message << duration_cast<milliseconds>(value).count() << " ms" << std::endl;
		std::cerr << os.str();
	}

	std::string message;
	steady_clock::duration value;
};
class AddDuration
{
public:
	explicit AddDuration(steady_clock::duration& dest) : m_addTo(dest), m_start(steady_clock::now())
	{
	}
	explicit AddDuration(TotalDuration& dest) : AddDuration(dest.value)
	{
	}
	~AddDuration()
	{
		m_addTo += steady_clock::now() - m_start;
	}

private:
	steady_clock::duration& m_addTo;
	steady_clock::time_point m_start;
};

#define ADD_DURATION(value) AddDuration UNIQ_ID(__LINE__){value};
