#pragma once

#include <ostream>
#include <string>
#include <unordered_map>
#include <utility>

#define CMD_ARG_FULL_NAME_INPUT_FILE "--input_file"
#define CMD_ARG_SHORT_NAME_INPUT_FILE "-ifile"
#define CMD_ARG_OUTPUT_FOLDER_NAME "--result_folder"

#define INPUT_TAG_END_OF_INPUT "[End]"
#define INPUT_TAG_NODES "[Nodes]"
#define INPUT_TAG_MATERIALS "[Materials]"
#define INPUT_TAG_FINITE_ELEMENTS "[Finite Elements]"
#define INPUT_TAG_LOADS "[tLoads]"

#define DIM_1_NODES_2 "Rod2"   // rod with two nodes
#define DIM_2_NODES_4 "Quad4"  // quadrilateral with four nodes
#define DIM_3_NODES_8 "Cube8"  // cube with 8 nodes

namespace FEM
{
namespace IO
{
enum class TAG
{
	END_OF_INPUT,
	NODES,
	MATERIALS,
	FINITE_ELEMENTS,
	LOADS
};
}  // namespace IO

enum class IO_FORMAT
{
	CHEKHOV,
	/*
	  [<tag_name>]
	  <row>
	  ...
	  [End]
	*/
};

enum class DEGREES_OF_FREEDOM
{
	FREE,
	FIXED,
};

enum class MATERIAL
{
	AL_ALLOY,  // Al_Alloy
	STEEL,	   // Steel
};

enum class CMD_ARGS
{
	INPUT_FILEPATH,
	OUTPUT_FOLDER,
	DUMMY
};

}  // namespace FEM

std::ostream& operator<<(std::ostream& os, const FEM::DEGREES_OF_FREEDOM& freedomDegrees);

extern std::unordered_map<std::string, FEM::CMD_ARGS> cmd_args;

extern std::unordered_map<std::string, std::string> cmd_args_alias;

extern std::unordered_map<std::string, FEM::IO::TAG> input_tags;

extern std::unordered_map<std::string, std::pair<size_t, size_t>> fe_type;
