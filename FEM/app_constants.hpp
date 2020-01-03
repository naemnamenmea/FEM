#pragma once

#include <string>
#include <unordered_map>

#define CMD_ARG_FULL_NAME_INPUT_FILE "--data"
#define CMD_ARG_SHORT_NAME_INPUT_FILE "-d"

namespace FEM
{
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
  STEEL,     // Steel
};

enum class FE_TYPE
{
  QUAD4,  // Quad4
  ROD2,   // Rod2
};

enum class CMD_ARGS
{
  DUMMY
};

}  // namespace FEM

std::unordered_map<std::string, FEM::CMD_ARGS> cmd_args = {
    {CMD_ARG_FULL_NAME_INPUT_FILE, FEM::CMD_ARGS::DUMMY},
};

std::unordered_map<std::string, std::string> arg_alias = {
    {CMD_ARG_SHORT_NAME_INPUT_FILE, CMD_ARG_FULL_NAME_INPUT_FILE}};