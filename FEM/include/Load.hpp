#pragma once

#include "HybridTensors.hpp"
#include "NameList.hpp"
#include "Node.hpp"
#include "FileInOut.hpp"
#include "femDllDefinitions.hpp"

#include <vector>
#include <string>

class tNodalTensor1Column;

class tNodalLoad
{
public:
	static tNodalLoad* NewNodalLoad(const std::string& kindName_)
	{
		return Factory.CallFunction(kindName_);
	}

	virtual ~tNodalLoad()
	{
	}
	virtual tNodalLoad& Link(const tNode&) = 0;
	virtual tNodalLoad& ReadData(std::istream&) = 0;
	virtual const std::string& Kind() const = 0;
	virtual const tNode& Node() const = 0;
	virtual tNodalTensor1& Level(
		tNodalTensor1& /* result_=tNodalTensor1()*/, real_t = 1.) const = 0;

protected:
	static const tNamesToFuns<tNodalLoad> Factory;
};

class tStatictNodalLoad : public tNodalLoad, public tNodalTensor1
{
public:
	static tNodalLoad* NewLoad()
	{
		return new tStatictNodalLoad();
	}

	tNodalLoad& Link(const tNode& node_) override
	{
		tNodalTensor1::Link(node_);
		return *this;
	}
	tNodalLoad& ReadData(std::istream&) override;
	const std::string& Kind() const override
	{
		return KindName;
	}
	const tNode& Node() const override
	{
		return tNodalTensor1::Node();
	}
	tNodalTensor1& Level(tNodalTensor1& /* result_=tNodalTensor1()*/, real_t = 1.) const override;

private:
	static const std::string& KindName;
};

class FEM_API tLoadCase
{
public:
	tLoadCase() : m_caption(), m_pNodalLoads() /*, CurrentLoadLevel(0.)*/
	{
	}
	~tLoadCase();
	tLoadCase(const std::string& name_)
		: m_caption(name_), m_pNodalLoads() /*, CurrentLoadLevel(0.)*/
	{
	}
	tLoadCase& SetName(const std::string& name_)
	{
		m_caption = name_;
		return *this;
	}
	const std::string& Name() const
	{
		return m_caption;
	}
	tNodalLoad& AddNodalForce(const std::string&, const tNode&);
	tNodalTensor1Column& LoadLevel(tNodalTensor1Column& result_, real_t = 1.) const;
	tLoadCase& ImportFromFile(ifstream_XML&, const tFE_model&);

private:
#pragma warning(push)
#pragma warning(disable : 4251)
	std::string m_caption;
	std::vector<tNodalLoad*> m_pNodalLoads;
	//   real_t CurrentLoadLevel;
#pragma warning(pop)
};

struct fHasCase
{
	fHasCase(const char* name_) : Name(name_)
	{
	}
	bool operator()(const tLoadCase& case_) const
	{
		return case_.Name() == Name;
	}

	const char* Name;
};

// class tFE_model;

class FEM_API tLoads
{
public:
	tLoads() : m_loadCases()
	{
	}
	tLoads& ImportFromFile(const char*, const tFE_model&);
	bool IsEmpty() const
	{
		return m_loadCases.empty();
	};
	size_t HowManyCases() const
	{
		return m_loadCases.size();
	}
	tLoadCase& Case(size_t);
	const tLoadCase& Case(size_t) const;
	tLoadCase& Case(const std::string&);
	const tLoadCase& Case(const std::string&) const;
	bool HasCase(const char* name_) const
	{
		return std::find_if(m_loadCases.begin(), m_loadCases.end(), fHasCase(name_)) !=
			   m_loadCases.end();
	}
	bool HasCase(const std::string& name_) const
	{
		return HasCase(name_.c_str());
	}
	void Clear()
	{
		m_loadCases.clear();
	}
	void Dump(std::ostream&, const tFE_model&) const;

private:
#pragma warning(suppress : 4251)
	std::vector<tLoadCase> m_loadCases;
};

inline tLoadCase& tLoads::Case(size_t i_)
{
#ifdef STRONGCHECK
	Assert(i_ > 0 && i_ <= m_loadCases.size(), "invalid index in tLoads::Case");
#endif
	return m_loadCases[--i_];
}

inline const tLoadCase& tLoads::Case(size_t i_) const
{
#ifdef STRONGCHECK
	Assert(i_ > 0 && i_ <= m_loadCases.size(), "invalid index in tLoads::Case const");
#endif
	return m_loadCases[--i_];
}
