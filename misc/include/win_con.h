#pragma once

#include <windows.h>
#include <algorithm>
#ifndef MAX_SCREEN_WIDTH
#define MAX_SCREEN_WIDTH 300
#endif

class tVisibleUI
{
public:
	tVisibleUI()
		: m_hCon(GetStdHandle(STD_OUTPUT_HANDLE)), m_beginPos(0), m_buffer(), m_tmp(), m_xy()
	{
#ifdef STRONGCHECK
		Assert(m_hCon != INVALID_HANDLE_VALUE, "Can't get handle");
#endif
		CONSOLE_SCREEN_BUFFER_INFO sb_inf;
		GetConsoleScreenBufferInfo(m_hCon, &sb_inf);
		m_endPos = sb_inf.srWindow.Right - sb_inf.srWindow.Left - 5;
		std::fill_n(const_cast<char*>(m_blank), MAX_SCREEN_WIDTH, ' ');
		std::fill_n(const_cast<char*>(m_bar), MAX_SCREEN_WIDTH, '°');
	}

	void PrintCaption(__int64 nomer_, const char* capt_, int min_, int max_)
	{
		m_xy.X = m_beginPos;
		m_xy.Y = static_cast<SHORT>(nomer_ * 2);
		WriteConsoleOutputCharacter(m_hCon, capt_, (DWORD)strlen(capt_), m_xy, &m_tmp);
		++m_xy.Y;
		_itoa_s(min_, m_buffer, 10);  // std::sprintf(buffer,"%i",min_);
		WriteConsoleOutputCharacter(m_hCon, m_buffer, (DWORD)strlen(m_buffer), m_xy, &m_tmp);
		m_xy.X = m_endPos;
		_itoa_s(max_, m_buffer, 10);  // std::sprintf(buffer,"%i",max_);
		WriteConsoleOutputCharacter(m_hCon, m_buffer, (DWORD)strlen(m_buffer), m_xy, &m_tmp);
	}

	void PrintCurrent(__int64 nomer_, int min_, int max_, int val_)
	{
#pragma warning(push)
#pragma warning(disable : 4996)
		const size_t min_len = strlen(itoa(min_, m_buffer, 10));  // std::sprintf(buffer,"%i",min_),
#pragma warning(pop)
		const size_t current_pos = (val_ - min_) * (m_endPos - m_beginPos) / (max_ - min_);
		m_xy.X = static_cast<SHORT>(m_beginPos + min_len);
		m_xy.Y = static_cast<SHORT>(nomer_ * 2 + 1);
		if (current_pos > min_len)
			WriteConsoleOutputCharacter(m_hCon, m_bar, DWORD(current_pos - min_len), m_xy, &m_tmp);

		_itoa_s(val_, m_buffer, 10);  // std::sprintf(buffer,"%i",val_);
		m_xy.X = static_cast<SHORT>(m_beginPos + current_pos);
		WriteConsoleOutputCharacter(m_hCon, m_buffer, (DWORD)strlen(m_buffer), m_xy, &m_tmp);
	}

	void Clear(__int64 nomer_)
	{
		m_xy.X = m_beginPos;
		m_xy.Y = static_cast<SHORT>(nomer_ * 2);
		WriteConsoleOutputCharacter(m_hCon, m_blank, m_endPos - m_beginPos + 5, m_xy, &m_tmp);
		m_xy.Y = static_cast<SHORT>(nomer_ * 2 + 1);
		WriteConsoleOutputCharacter(m_hCon, m_blank, m_endPos - m_beginPos + 5, m_xy, &m_tmp);
	}

private:
	HANDLE m_hCon;
	short int m_beginPos, m_endPos;
	char m_blank[MAX_SCREEN_WIDTH];
	char m_bar[MAX_SCREEN_WIDTH];
	DWORD m_tmp;
	char m_buffer[25];
	COORD m_xy;
};
