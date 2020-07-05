#pragma once

#include <sstream>
#include <string>
#include <iostream>
#include <tchar.h>

#if defined(linux) && __GNUC__ >= 4
#define THROW_MESSAGE(errorText) \
	throw tMessage(std::string(__func__) + _T(": ") + std::string(errorText))
#else
#define THROW_MESSAGE(errorText) \
	throw tMessage(std::string(_T(__FUNCTION__)) + _T(": ") + std::string(errorText))
#endif

class tMessage
{
public:
	tMessage(const char* newMessage_) : m_data(newMessage_)
	{
	}
	tMessage(const std::string& newMessage_) : m_data(newMessage_)
	{
	}
	tMessage& operator=(const tMessage& o_)
	{
		m_data = o_.m_data;
		return *this;
	}
	tMessage& operator<<(const std::string& str_)
	{
		m_data += str_;
		return *this;
	}
	tMessage& operator<<(const char* str_)
	{
		m_data += str_;
		return *this;
	}
	tMessage& operator<<(char ch_)
	{
		m_buffer[0] = ch_;
		m_buffer[1] = '\0';
		m_data += m_buffer;
		return *this;
	}
	template <typename T>
	tMessage& operator<<(const T&);
	operator const char*() const
	{
		return m_data.c_str();
	}
	operator char*() const
	{
		m_data.copy(m_buffer, m_data.length());
		return m_buffer;
	}

	std::string m_data;

private:
	static char m_buffer[];
};

inline std::ostream& operator<<(std::ostream& out, const tMessage& msg)
{
	out << msg.m_data;
	return out;
}

inline void Assert(bool assertion_, const tMessage& msg_)
{
	if (!assertion_)
		throw msg_;
}

template <typename T>
inline tMessage& tMessage::operator<<(const T& v_)
{
	/*#if !defined __BORLANDC__  ||  __BORLANDC__ >= 0x550
	 std::
	#endif*/
	/* ostrstream(Buffer,4) << v_ << '\0';
	 Data += Buffer;*/
	std::ostringstream os;
	os << v_;
	m_data += os.str();
	return *this;
}
