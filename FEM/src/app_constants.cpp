#include "stdafx.hpp"

#include "app_constants.hpp"

std::unordered_map<std::string, FEM::CMD_ARGS> cmd_args = {
	{CMD_ARG_FULL_NAME_INPUT_FILE, FEM::CMD_ARGS::INPUT_FILEPATH},
	{CMD_ARG_OUTPUT_FOLDER_NAME, FEM::CMD_ARGS::OUTPUT_FOLDER},
};

std::unordered_map<std::string, std::string> cmd_args_alias = {
	{CMD_ARG_SHORT_NAME_INPUT_FILE, CMD_ARG_FULL_NAME_INPUT_FILE},
};

std::unordered_map<std::string, FEM::IO::TAG> input_tags = {
	{INPUT_TAG_END_OF_INPUT, FEM::IO::TAG::END_OF_INPUT},
	{INPUT_TAG_NODES, FEM::IO::TAG::NODES},
	{INPUT_TAG_MATERIALS, FEM::IO::TAG::MATERIALS},
	{INPUT_TAG_FINITE_ELEMENTS, FEM::IO::TAG::FINITE_ELEMENTS},
	{INPUT_TAG_LOADS, FEM::IO::TAG::LOADS},
};

std::unordered_map<std::string, std::pair<size_t, size_t>> fe_type = {
	{DIM_1_NODES_2, {1, 2}},
	{DIM_2_NODES_4, {2, 4}},
	{DIM_3_NODES_8, {3, 8}},
};

std::ostream& operator<<(std::ostream& os, const FEM::DEGREES_OF_FREEDOM& freedomDegrees)
{
	os << (freedomDegrees == FEM::DEGREES_OF_FREEDOM::FIXED ? "fixed" : "free");
	return os;
}
