#pragma once

#define MAX_NUM_OF_PROCESS_DISPLAYS 12
#define CAPTION_LEN 75

#ifdef STRONGCHECK
#include "ThrowMessage.hpp"
#endif

#if (defined __BORLANDC__ && defined __WIN32__ && defined __CONSOLE__) || \
	(defined _MSC_VER && defined _WIN32 /*&& defined __CONSOLE__*/) ||    \
	defined __MINGW32__ /*|| defined __CYGWIN__*/

#include "win_con.h"

#elif defined __GNUG__
#include "unix_con.h"
#else

struct tVisibleUI
{
	void PrintCaption(int, const char*, int, int)
	{
	}
	void PrintCurrent(int, int, int, int)
	{
	}
	void Clear(int nomer_)
	{
	}
};
#endif

class tPDArray;

extern tPDArray PDArray;

typedef int data_t;

class tProcessDisplayCore
{
public:
	void Reset(const char* capt_, data_t min_, data_t max_, data_t step_);
	tProcessDisplayCore& operator++();
	//    tProcessDisplayCore& operator-- ()
	bool Done() const
	{
		return m_curValue >= m_maxValue;
	}

private:
	friend class tPDArray;

	char m_caption[CAPTION_LEN];
	data_t m_minValue, m_maxValue, m_step, m_curValue;
};

class tPDArray
{
public:
	tPDArray() : m_array(), m_pTop(m_array)
	{
	}

	tProcessDisplayCore* Push(const char* capt_, data_t min_, data_t max_, data_t stp_)
	{
#ifdef STRONGCHECK
		Assert(stp_ > 0, "negative step in tPDArray::Push");
#endif
		//        for (tProcessDisplayCore* p=pTop-1; p>=Array && p->Done(); --p)
		//                 Pop(p);
		if (m_pTop >= m_array + MAX_NUM_OF_PROCESS_DISPLAYS)
			return nullptr;
		m_pTop->Reset(capt_, min_, max_, stp_);
		return m_pTop++;
	}

	void Pop(const tProcessDisplayCore* p_)
	{
		if (p_ == nullptr /*&VoidDspl || p_>=pTop*/)
			return;
#ifdef STRONGCHECK
		Assert(m_pTop >= m_array, "Deleting from empty array of process displays");
#endif
		--m_pTop;
		m_userScreen.Clear(m_pTop - m_array);
	}

	void PrintCaption(const tProcessDisplayCore* pdisplay_)
	{
		m_userScreen.Clear(pdisplay_ - m_array);
		m_userScreen.PrintCaption(
			pdisplay_ - m_array,
			pdisplay_->m_caption,
			pdisplay_->m_minValue,
			pdisplay_->m_maxValue);
	}
	void PrintCurrent(const tProcessDisplayCore* pdisplay_)
	{
		m_userScreen.PrintCurrent(
			pdisplay_ - m_array,
			pdisplay_->m_minValue,
			pdisplay_->m_maxValue,
			pdisplay_->m_curValue);
	}

private:
	tProcessDisplayCore m_array[MAX_NUM_OF_PROCESS_DISPLAYS];
	tProcessDisplayCore* m_pTop;

	tVisibleUI m_userScreen;
};

inline void tProcessDisplayCore::Reset(const char* capt_, data_t min_, data_t max_, data_t step_)
{
#pragma warning(push)
#pragma warning(disable : 4996)
	strcpy(m_caption, capt_);
#pragma warning(pop)
	m_minValue = m_curValue = min_;
	m_maxValue = max_;
	m_step = step_;
	PDArray.PrintCaption(this);
}
inline tProcessDisplayCore& tProcessDisplayCore::operator++()
{
	m_curValue += m_step;
	if (m_curValue > m_maxValue)
		m_curValue = m_maxValue;
	PDArray.PrintCurrent(this);
	return *this;
}

class tProcessDisplay
{
public:
	tProcessDisplay(const char* capt_, data_t min_, data_t max_, data_t stp_)
	{
		m_pCore = PDArray.Push(capt_, min_, max_, stp_);
	}
	tProcessDisplay(const char* capt_, data_t min_, data_t max_)
	{
		m_pCore = PDArray.Push(capt_, min_, max_, 1);
	}
	tProcessDisplay(const char* capt_, data_t max_)
	{
		m_pCore = PDArray.Push(capt_, 0, max_, 1);
	}
	~tProcessDisplay()
	{
		PDArray.Pop(m_pCore);
	}

	void Reset(const char* capt_, data_t min_, data_t max_, data_t step_)
	{
		if (!is_void())
			m_pCore->Reset(capt_, min_, max_, step_);
	}
	void Reset(const char* capt_, data_t min_, data_t max_)
	{
		if (!is_void())
			m_pCore->Reset(capt_, min_, max_, 1);
	}
	void Reset(const char* capt_, data_t max_)
	{
		if (!is_void())
			m_pCore->Reset(capt_, 0, max_, 1);
	}

	tProcessDisplay& operator++()
	{
		if (!is_void())
			++(*m_pCore);
		return *this;
	}
	bool Done() const
	{
		return is_void() ? true : m_pCore->Done();
	}

private:
	bool is_void() const
	{
		return m_pCore == nullptr;
	}

	tProcessDisplayCore* m_pCore;
};
