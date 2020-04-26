#pragma once

#include <chrono>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;
using namespace std::chrono;

class LogDuration
{
public:
	explicit LogDuration(const string& msg = "") : m_message(msg + ": "), m_start(steady_clock::now())
	{
	}

	~LogDuration()
	{
		auto finish = steady_clock::now();
		auto dur = finish - m_start;
		cerr << m_message << duration_cast<milliseconds>(dur).count() << " ms" << endl;
	}

private:
	string m_message;
	steady_clock::time_point m_start;
};

#define UNIQ_ID_IMPL(lineno) _a_local_var_##lineno
#define UNIQ_ID(lineno) UNIQ_ID_IMPL(lineno)

#define LOG_DURATION(message) LogDuration UNIQ_ID(__LINE__){message};

struct TotalDuration
{
	string message;
	steady_clock::duration value;
	explicit TotalDuration(const string& msg = "") : message(msg + ": "), value(0)
	{
	}
	~TotalDuration()
	{
		ostringstream os;
		os << message << duration_cast<milliseconds>(value).count() << " ms" << endl;
		cerr << os.str();
	}
};
class AddDuration
{
public:
	explicit AddDuration(steady_clock::duration& dest) : m_add_to(dest), m_start(steady_clock::now())
	{
	}
	explicit AddDuration(TotalDuration& dest) : AddDuration(dest.value)
	{
	}
	~AddDuration()
	{
		m_add_to += steady_clock::now() - m_start;
	}

private:
	steady_clock::duration& m_add_to;
	steady_clock::time_point m_start;
};

#define ADD_DURATION(value) AddDuration UNIQ_ID(__LINE__){value};
