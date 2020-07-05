#pragma once

#include <cstdio>
#include <string>
#include <iterator>

#pragma warning(push)
#pragma warning(disable : 4996)
class ifstream_c
{
public:
	ifstream_c(const char* filename_ = nullptr)
		: m_pStream(filename_ == nullptr ? (std::FILE*)nullptr : std::fopen(filename_, "r")),
		  m_lastUnformattedRead(0)
	{
	}
	virtual ~ifstream_c()
	{
		if (m_pStream != nullptr)
			std::fclose(m_pStream);
	}
	void open(const char* filename_)
	{
		if (m_pStream != nullptr)
			std::fclose(m_pStream);
		m_pStream = std::fopen(filename_, "r");
	}
	void close()
	{
		if (m_pStream != nullptr)
			std::fclose(m_pStream);
		m_pStream = nullptr;
	}
	bool fail() const
	{
		return m_pStream == nullptr ? false : /*std::*/ ferror(m_pStream);
	}
	bool bad() const
	{
		return m_pStream == nullptr ? false : /*std::*/ ferror(m_pStream);
	}
	bool operator!() const
	{
		return fail();
	}
	void clear()
	{
	}
	bool eof() const
	{
		return m_pStream == nullptr ? true : /*std::*/ feof(m_pStream);
	}
	long tellg()
	{
		return m_pStream == nullptr ? -1l : std::ftell(m_pStream);
	}
	ifstream_c& seekg(long& offset_)
	{
		if (m_pStream != nullptr)
			std::fseek(m_pStream, offset_, SEEK_SET);
		return *this;
	}
	int get()
	{
		return m_pStream == nullptr ? -1 : std::fgetc(m_pStream);
	}
	ifstream_c& get(char& c_)
	{
		if (m_pStream != nullptr)
			c_ = static_cast<char>(std::fgetc(m_pStream));
		return *this;
	}
	ifstream_c& unget()
	{
		if (m_pStream != nullptr)
			std::fseek(m_pStream, -1, SEEK_CUR);
		return *this;
	}
	ifstream_c& ignore(int n_ = 1)
	{
		if (m_pStream != nullptr)
			std::fseek(m_pStream, n_, SEEK_CUR);
		return *this;
	}
	int peek()
	{
		if (m_pStream == nullptr)
			return -1;
		int c = std::fgetc(m_pStream);
		unget();
		return c;
	}
	ifstream_c& read(char* str_, long n_)
	{
		if (m_pStream != nullptr)
		{
			std::fgets(str_, n_, m_pStream);
			m_lastUnformattedRead = strlen(str_);
		}
		return *this;
	}
	size_t gcount() const
	{
		return m_lastUnformattedRead;
	}
	ifstream_c& operator>>(char& x_)
	{
		if (m_pStream != nullptr)
			get(x_);
		return *this;
	}
	ifstream_c& operator>>(int& x_)
	{
		if (m_pStream != nullptr)
			std::fscanf(m_pStream, "%d", &x_);
		return *this;
	}
	ifstream_c& operator>>(long int& x_)
	{
		if (m_pStream != nullptr)
			std::fscanf(m_pStream, "%dl", &x_);
		return *this;
	}
	ifstream_c& operator>>(unsigned int& x_)
	{
		if (m_pStream != nullptr)
			std::fscanf(m_pStream, "%u", &x_);
		return *this;
	}
	ifstream_c& operator>>(unsigned long int& x_)
	{
		if (m_pStream != nullptr)
			std::fscanf(m_pStream, "%ul", &x_);
		return *this;
	}
	ifstream_c& operator>>(float& x_)
	{
		if (m_pStream != nullptr)
			std::fscanf(m_pStream, "%f", &x_);
		return *this;
	}
	ifstream_c& operator>>(double& x_)
	{
		if (m_pStream != nullptr)
			std::fscanf(m_pStream, "%lf", &x_);
		return *this;
	}
	ifstream_c& operator>>(long double& x_)
	{
		if (m_pStream != nullptr)
			std::fscanf(m_pStream, "%Lf", &x_);
		return *this;
	}
	ifstream_c& operator>>(std::string&);
	ifstream_c& getline(std::string&, char eol_);

private:
	std::FILE* m_pStream;
	size_t m_lastUnformattedRead;
};
#pragma warning(pop)

template <typename T>
class ifstream_c_iterator : public std::iterator<std::input_iterator_tag, T>
{
public:
	ifstream_c_iterator() : m_pStream(nullptr)
	{
	}
	ifstream_c_iterator(ifstream_c& ifs_) : m_pStream(&ifs_)
	{
		operator++();
	}

	const T& operator*() const
	{
		return m_value;
	}

	ifstream_c_iterator& operator++()
	{
		m_pStream->operator>>(m_value);
		if (m_pStream->eof())
			m_pStream = nullptr;
		return *this;
	}

	bool operator==(const ifstream_c_iterator& o_) const
	{
		return m_pStream == o_.m_pStream;
	}
	bool operator!=(const ifstream_c_iterator& o_) const
	{
		return !(*this == o_);
	}

private:
	ifstream_c* m_pStream;
	T m_value;
};
