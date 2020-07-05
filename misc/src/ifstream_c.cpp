#include "ifstream_c.hpp"
#include <cctype>

using namespace std;

ifstream_c& ifstream_c::operator>>(string& str_)
{
	if (m_pStream != nullptr)
	{
		str_.clear();
		char ch;
		for (;;)
		{
			ch = static_cast<char>(fgetc(m_pStream));
			if (ch == EOF)
				return *this;
			if (isspace(ch))
			{
				unget();
				return *this;
			}
			str_ += ch;
		}
	}
	return *this;
}

ifstream_c& ifstream_c::getline(string& str_, char eol_)
{
	if (m_pStream != nullptr)
	{
		str_.clear();
		for (char ch;;)
		{
			ch = static_cast<char>(fgetc(m_pStream));
			if (ch == EOF)
				return *this;
			if (ch == eol_)
			{
				//                 unget();
				return *this;
			}
			str_ += ch;
		}
	}
	return *this;
}
