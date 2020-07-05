#pragma once

#include <string>
#include <filesystem>

namespace fs = std::filesystem;

class CmdLineParams
{
public:
	static char KEY_VALUE_SEP;

	CmdLineParams();

	CmdLineParams(int argc, const char* argv[]);

	void SetInputFilePath(const std::string& path);
	void SetOutputFolderName(const std::string& path);

	fs::path GetOutputFilePath() const;
	const fs::path& GetInputFilePath() const;
	const fs::path& GetOutputFolder() const;

private:
	void SetParams(unsigned int argc, const char* argv[]);
	void InitToDefaults();

	fs::path m_inputFileName;
	fs::path m_outputFolder;
};