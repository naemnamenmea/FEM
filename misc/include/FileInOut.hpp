#pragma once

#include "ThrowMessage.hpp"
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>

#if /*!defined __GNUG__  ||*/ defined __MINGW32__ || defined __BORLANDC__
//#if defined __BORLANDC__ && defined _MSC_VER
#include "ifstream_c.h"
#endif

#ifdef IFSTREAM_C
#define IFSTREAM ifstream_c
#define GETLINE(STREAM, WORD, EOL) (STREAM).getline((WORD), (EOL))
#define ISTREAM_ITERATOR ifstream_c_iterator
#else
#include <iterator>
#define IFSTREAM std::ifstream
#define GETLINE(STREAM, WORD, EOL) std::getline(STREAM, (WORD), (EOL))
#define ISTREAM_ITERATOR std::istream_iterator
#endif

//#if defined __BORLANDC__  &&  __BORLANDC__ < 0x550  ||  defined __GNUG__
//#define POS_TYPE long
//#else
//  #define POS_TYPE /*std::*/pos_type
//#endif

class tXML_Reader;

class ifstream_XML : public IFSTREAM
{
public:
	virtual ~ifstream_XML();
	ifstream_XML(const char* file_name_, bool = false);
	bool FindChar(char);
	char FindChar1stOf(const char*);
	bool FindWord(const char*);
	const char* FindWord1stOf(const char**, size_t);
	std::string& ReadWord(std::string&);
	const std::string& ReadWord();
	std::string& ReadWord(std::string&, char eol_);
	const std::string& ReadWord(char eol_);
	size_t WordCount(const std::string&);
	size_t WordCountBefore(const char*, const char*);

	bool FindHighLevelTag(const char*);
	void ReadObject(tXML_Reader& reader_, const std::string& tag = std::string());
	size_t CountArray(const char*);
	size_t CountItems(const char*, const char*);
	char ReadAttribs(tXML_Reader& reader_);
	void SkipObject();

private:
	const std::string m_commentBegin, m_commentEnd;
	char m_noCommentTmpFileName[L_tmpnam + 5];
};

class tXML_Reader
{
public:
	virtual ~tXML_Reader()
	{
	}
	virtual void ReadData(const std::string& /*tag_*/, ifstream_XML& /*file_*/)
	{
	}
	virtual void SetValue(const std::string& /*name_*/, const std::string& /*value_*/)
	{
	}
};

class tOutput_XML
{
public:
	tOutput_XML(const char* filename_) : m_outFile(filename_), m_spacesNum(0)
	{
	}
	bool operator!() const
	{
		return !m_outFile;
	}
	void OpenTag(const char* const tag_name_, bool has_attrs)
	{
		for (size_t i = 1; i <= m_spacesNum; ++i) m_outFile.put(' ');
		m_spacesNum += 3;
		m_outFile << '<' << tag_name_;
		if (!has_attrs)
			m_outFile << ">\n";
	}
	void CloseTag(const char* const tag_name_)
	{
		if (m_spacesNum > 2)
			m_spacesNum -= 3;
		else
			m_spacesNum = 0;

		for (size_t i = 1; i <= m_spacesNum; ++i) m_outFile.put(' ');

		m_outFile << "</" << tag_name_ << ">\n";
	}
	void CloseTag()
	{
		m_outFile << "/>\n";
		m_spacesNum -= 3;
	}

	template <typename Attr_type>
	tOutput_XML& AddAtr(const char* const attr_name_, Attr_type Attr_value_, bool islast = false)
	{
		m_outFile << ' ' << attr_name_ << "=\"" << Attr_value_ << '\"';
		if (islast)
			m_outFile << ">\n";	 // CloseTag();
		return *this;
	}

private:
	std::ofstream m_outFile;
	size_t m_spacesNum;
};

inline bool ifstream_XML::FindChar(char ch_)
{
	return std::find(ISTREAM_ITERATOR<char>(*this), ISTREAM_ITERATOR<char>(), ch_) !=
		   ISTREAM_ITERATOR<char>();
}

inline char ifstream_XML::FindChar1stOf(const char* chars_)
{
	return *std::find_first_of(
		ISTREAM_ITERATOR<char>(*this),
		ISTREAM_ITERATOR<char>(),
		chars_,
		chars_ + std::strlen(chars_));
}

inline std::string& ifstream_XML::ReadWord(std::string& word_)
{
	std::streampos beg_pos = tellg();
	(*this) >> word_;
	std::string::size_type pos = word_.find('>');
	if (pos != std::string::npos && pos != word_.size() - 1)
	{
		//    seekg(beg_pos + pos+1);
		seekg(beg_pos);
		ignore(pos + 1);
		word_.resize(pos + 1);
	}
	return word_;
}

inline std::string& ifstream_XML::ReadWord(std::string& word_, char eol_)
{
	GETLINE(*this, word_, eol_);
	return word_;
}

inline size_t ifstream_XML::WordCount(const std::string& word_)
{
	return count(ISTREAM_ITERATOR<std::string>(*this), ISTREAM_ITERATOR<std::string>(), word_);
}

inline size_t ifstream_XML::WordCountBefore(const char* word_, const char* end_word_)
{
	std::streampos beg_pos = tellg();
	size_t result(0);
	const char* words[] = {word_, end_word_};
	for (;; ++result)
		if (FindWord1stOf(words, 2) == end_word_)
			break;
	seekg(beg_pos);
	return result;
}

inline size_t ifstream_XML::CountItems(const char* end_tag_, const char* item_tag_)
{
	std::streampos beg_pos = tellg();
	const size_t items_num(WordCountBefore(item_tag_, end_tag_));
	Assert(
		!bad(),
		tMessage("Error while counting items \"")
			<< item_tag_ << "\" within object \"" << end_tag_ << '"');
	seekg(beg_pos);
	return items_num;
}

inline void ifstream_XML::SkipObject()
{
	const char terms[] = "<>";
	for (int n = 1; n > 0;)
		if (FindChar1stOf(terms) == terms[0])
			++n;
		else
			--n;
}
