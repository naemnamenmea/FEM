#include "app_constants.hpp"
#include "cmdLineParams.hpp"
#include <sstream>

namespace fs = std::filesystem;

char CmdLineParams::KEY_VALUE_SEP = '=';

namespace
{
auto GetParamPointer(const std::string& key)
{
	auto it = cmd_args.find(key);
	if (it == cmd_args.end() && cmd_args_alias.count(key) != 0)
	{
		it = cmd_args.find(cmd_args_alias[key]);
	}

	if (it == cmd_args.end())
	{
		std::stringstream error_msg;
		error_msg << "Unknown parameter: \"" << key << "\"";
		throw std::runtime_error(error_msg.str());
	}

	return it;
}
}  // namespace

void CmdLineParams::SetParams(unsigned int argc, const char* argv[])
{
	for (unsigned int i = 1; i < argc; i += 2)
	{
		std::string paramPair = argv[i];

		size_t sepPos = paramPair.find(KEY_VALUE_SEP);
		bool isKeyOnly = sepPos == paramPair.npos;

		std::string key = paramPair.substr(0, isKeyOnly ? paramPair.size() : sepPos);

		auto it = cmd_args.end();

		it = GetParamPointer(key);

		if (!isKeyOnly)
		{
			std::string value = paramPair.substr(sepPos + 1);

			FEM::CMD_ARGS param = it->second;

			if (param == FEM::CMD_ARGS::DUMMY)
			{
			}
			else if (param == FEM::CMD_ARGS::INPUT_FILEPATH)
			{
				SetInputFilePath(value);
			}
			else if (param == FEM::CMD_ARGS::OUTPUT_FOLDER)
			{
				SetOutputFolderName(value);
			}
			else
			{
			}
		}
	}
}

void CmdLineParams::InitToDefaults()
{
}

CmdLineParams::CmdLineParams()
{
	InitToDefaults();
}

CmdLineParams::CmdLineParams(int argc, const char* argv[])
{
	SetParams(argc, argv);
}

void CmdLineParams::SetInputFilePath(const std::string& path)
{
	m_inputFileName = path;
}

void CmdLineParams::SetOutputFolderName(const std::string& path)
{
	m_outputFolder = path;
}

const fs::path& CmdLineParams::GetInputFilePath() const
{
	return m_inputFileName;
}

fs::path CmdLineParams::GetOutputFilePath() const
{
	fs::path outputFileName = m_inputFileName.filename();
	outputFileName.replace_extension("_out" + m_inputFileName.extension().string());
	return m_outputFolder / outputFileName;
}

const fs::path& CmdLineParams::GetOutputFolder() const
{
	return m_outputFolder;
}
