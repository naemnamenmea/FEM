#include "FileInOut.hpp"

using namespace std;

ifstream_XML::ifstream_XML(const char* file_name_, bool no_comments_)
	: m_commentBegin("<!--"), m_commentEnd("-->"), m_noCommentTmpFileName()
{
	open(file_name_);

	if (fail())
		return;

	FindWord("<?xml");
	SkipObject();
	if (no_comments_)
	{
		*m_noCommentTmpFileName = '\0';
		return;
	}
	/*char nocomment_fname[L_tmpnam+5]("tmp.xml");*/
#pragma warning(push)
#pragma warning(disable : 4996)
	tmpnam(m_noCommentTmpFileName);
	strcpy(m_noCommentTmpFileName + strlen(m_noCommentTmpFileName), ".xml");
#pragma warning(pop)
	ofstream dest_without_comments(m_noCommentTmpFileName);
	const unsigned short int bufsize(256);
	char buffer[bufsize], *buf_end(buffer + bufsize), *buf_beg(buf_end);
	unsigned short int buf_vol(0);

	for (;;)
	{
		if (buf_vol < m_commentBegin.size())
			if (eof())
			{
				dest_without_comments.write(buf_beg, buf_vol);
				break;
			}
			else
			{
				copy(buf_beg, buf_end, buffer);
				read(buffer + buf_vol, bufsize - buf_vol);
				buf_vol += static_cast<unsigned short>(gcount());
				buf_beg = buffer;
				buf_end = buf_beg + buf_vol;
			}

		if (m_commentBegin.compare(0, m_commentBegin.size(), buf_beg, m_commentBegin.size()) == 0)
		{
			char* p = search(
				buf_beg + m_commentBegin.size(), buf_end, m_commentEnd.begin(), m_commentEnd.end());
			if (p < buf_end)
			{
				buf_beg = p + m_commentEnd.size();
				buf_vol = static_cast<unsigned short>(buf_end - buf_beg);
			}
			else
			{
				for (size_t i = 1; i < m_commentEnd.size(); ++i) unget();
				ISTREAM_ITERATOR<char> it = search(
					ISTREAM_ITERATOR<char>(*this),
					ISTREAM_ITERATOR<char>(),
					m_commentEnd.begin(),
					m_commentEnd.end());
				read(buffer, bufsize);
				buf_vol = static_cast<unsigned short>(gcount());
				buf_beg = buffer;
				buf_end = buf_beg + buf_vol;
			}
		}

		dest_without_comments.write(buf_beg, 1);
		++buf_beg;
		--buf_vol;
	}
	close();
	dest_without_comments.close();
	clear();
	open(m_noCommentTmpFileName);
}

ifstream_XML::~ifstream_XML()
{
	if (strlen(m_noCommentTmpFileName) > 0)
	{
		close();
		remove(m_noCommentTmpFileName);
	}
}

bool ifstream_XML::FindWord(const char* word_)	// isn't inline because of BCB6's gluck
{
#ifndef _MSC_VER
	return search(
			   ISTREAM_ITERATOR<char>(*this),
			   ISTREAM_ITERATOR<char>(),
			   word_,
			   word_ + strlen(word_)) != ISTREAM_ITERATOR<char>();
#else
	/* const char* words[] = {word_};
	 return FindWord1stOf(words,1) == word_; */
	char symbolToCheck;
	for (;;)
	{
		if (!FindChar(word_[0]))
			return false;
		istream::pos_type position = tellg();
		const size_t wordLength = strlen(word_);
		for (int i = 1; i < wordLength; ++i)
		{
			get(symbolToCheck);
			if (eof())
				return false;
			if (symbolToCheck != word_[i])
			{
				seekg(position);
				break;
			}
			if (i == wordLength - 1)
				return true;
		}
	}
#endif
}

const string& ifstream_XML::ReadWord()
{
	static string last_read;
	ReadWord(last_read);
	return last_read;
}

const string& ifstream_XML::ReadWord(char eol_)
{
	static string last_read;
	ReadWord(last_read, eol_);
	return last_read;
}

struct nonzero
{
	bool operator()(size_t n)
	{
		return n != 0;
	}
};

const char* ifstream_XML::FindWord1stOf(const char** words_, const size_t n_)
{
	char* begs = new char[n_ + 1];
	size_t i(0);
	for (; i < n_; ++i) begs[i] = words_[i][0];
	begs[n_] = '\0';

	vector<size_t> counts(n_, 0);

	for (char ch;;)
	{
		ch = FindChar1stOf(begs);
		streampos beg_pos = tellg();
		for (;;)
		{
			for (i = 0; i < n_; ++i)
				if (ch == words_[i][counts[i]])
				{
					if (++counts[i] == strlen(words_[i]))
					{
						delete[] begs;
						return words_[i];
					}
				}
				else
					counts[i] = 0;
			if (find_if(counts.begin(), counts.end(), nonzero()) == counts.end())
			{
				seekg(beg_pos);
				break;
			}
			get(ch);
		}
	}
	delete[] begs;
	tMessage msg("No of words: ");
	for (i = 0; i < n_; ++i) msg << '"' << words_[i] << "\" ";
	throw msg << "- has been found in XML data file";
}

bool ifstream_XML::FindHighLevelTag(const char* tag_)
{
	static string buffer;
	streampos beg_pos = tellg();

	size_t num_of_highlevel_subtags(0);
	for (;; ++num_of_highlevel_subtags)
	{
		FindChar('<');
		ReadWord(buffer);
		if (buffer[0] == '/' || eof())	// end of current level found
			break;

		size_t buf_len(buffer.length());
		if (buffer[buf_len - 1] == '>')
		{
			--buf_len;
			if (buffer[buf_len - 1] == '/')
				--buf_len;
		}
		if (buffer.compare(0, buf_len, tag_, strlen(tag_)) == 0)
			return true;

		// find end of current object
		if (buffer[buffer.length() - 2] != '/')
			if (FindChar1stOf("/>") == '>')
			{
				buffer.insert(0, "</");
				FindWord(buffer.c_str());
			}
	}

	clear();
	seekg(beg_pos);
	for (size_t i(0); i < num_of_highlevel_subtags; ++i)
	{
		FindChar('<');
		if (FindHighLevelTag(tag_))
			return true;
	}
	return false;
}

void ifstream_XML::ReadObject(tXML_Reader& reader_, const string& tag_)
{
	string tag;
	if (tag_.empty())
	{
		Assert(FindChar('<'), "Character \"<\" not found in XML data file");
		ReadWord(tag);
	}
	else
	{
		tag = '<';
		tag += tag_;
		Assert(FindWord(tag.c_str()), tMessage("Tag \"") << tag << "\" not found in XML data file");
		tag.erase(0, 1);
		char ch = static_cast<char>(peek());
		if (ch == '>')
		{
			ignore();
			tag += ch;
		}
	}
	Assert(!bad(), tMessage("Possible unexpected end of file before reading object ") << tag_);
	Assert(
		tag[0] != '/', tMessage("Finishing tag found \"") << tag << "\" (starting one expected)");

	if (tag[tag.length() - 1] == '>')
	{
		Assert(
			tag[tag.length() - 2] != '/',
			tMessage("Object \"") << tag << "\" contains tag only, no data");
		tag.resize(tag.length() - 1);  // pop_back();
	}
	else if (ReadAttribs(reader_) == '/')  //
	{
		ignore();
		return;
	}
	reader_.ReadData(tag, *this);
	Assert(!fail(), tMessage("Error while reading data for object \"") << tag << '"');

	string next_tag;
	for (;;)
	{
		FindChar('<');
		streampos beg_pos = tellg();
		ReadWord(next_tag);
		Assert(
			!bad(),
			tMessage("Possible unexpected end of file while reading next data for \"")
				<< tag << '"');
		if (next_tag[0] == '/')
			break;
		seekg(beg_pos -= 1);
		reader_.ReadData(next_tag, *this);
	}
}

char ifstream_XML::ReadAttribs(tXML_Reader& reader_)
{
	static string attr_name, value;
	char ch;
	do
	{
		GETLINE(*this, attr_name, '=');
		// trimming leading and trailing spaces:
		while (isspace(attr_name[attr_name.size() - 1]))
			attr_name.resize(attr_name.size() - 1);	 // pop_back();
		while (isspace(attr_name[0])) attr_name.erase(0, 1);
		get(ch);  // - read '"'
		GETLINE(*this, value, ch);
		reader_.SetValue(attr_name, value);
		Assert(
			!fail(),
			tMessage("Error while reading attributes. Last attribute = \"")
				<< attr_name << "\", value = \"" << value << '"');
		for (;;)
		{
			ch = static_cast<char>(peek());
			if (ch == ' ' || ch == '?')
				get();
			else
				break;
		}
	} while (ch != '>' && ch != '/');
	return ch;
}

size_t ifstream_XML::CountArray(const char* array_tag_)
{
	if (!FindWord(array_tag_))
		// if (!FindHighLevelTag(array_tag_)) // Почему-то очень долго ищет тег finite_elements для
		// больших моделей, а при 10000 КЭ программа вылетает
		throw tMessage("Can't find tag ") << array_tag_;

	streampos beg_pos = tellg();

	FindChar('<');
	unget();
	string item_tag;
	ReadWord(item_tag);

	if (item_tag[item_tag.length() - 1] !=
		'>')  // because of item_tag may be substring of array_tag_
		item_tag += to_string(get());
	const size_t array_dim(CountItems(array_tag_, item_tag.c_str()) + 1);
	seekg(beg_pos);
	// reader_.ReadArray(array_tag_,array_dim,*this);
	return array_dim;
}
