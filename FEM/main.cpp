#include "stdafx.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include "app_constants.hpp"
#include "auxiliary.hpp"
#include "numeric_math.hpp"
#include "point3d.hpp"
#include "test_runner.hpp"
#include "tests.hpp"

namespace fs = std::filesystem;

class Node
{
};

class Material
{
 public:
  Material() {}
  Material(std::string name) : name(name), density(0) {}

  void SetDensity(double newDensity) { density = newDensity; }

  double GetDensity() const { return density; }

 private:
  std::string name;
  double density;
};

/*
 * 1 dim - 2 nodes
 * 2 dim - 4 nodes
 * 3 dim - 8 nodes
 */
class FiniteElement
{
 public:
 protected:
 private:
  std::vector<std::shared_ptr<Node>> nodePointers;
};

class FiniteElement1d : public FiniteElement
{
};

class FiniteElement2d : public FiniteElement
{
};

class FiniteElement3d : public FiniteElement
{
};

template <typename Point>
class FiniteElementModel
{
 public:
  FiniteElementModel()
  {
    size_t spaceDim = (sizeof(Point) / sizeof(typename Point::ValueType));

    if (spaceDim > 3)
    {
      throw std::runtime_error(spaceDim + " dimension is not supported");
    }

    std::cerr << "You are working in " << spaceDim << "-dimentional space"
              << std::endl;
    std::cerr << "Nodes for FE: " << quick_pow(2, spaceDim) << std::endl;
  }
  void SetParams(unsigned int argc, const char* argv[])
  {
    auto isValidParam = [argc, &argv](unsigned int i) -> bool {
      return i + 1 < argc && (i + 2 >= argc || argv[i + 2][0] == '-');
    };

    for (unsigned int i = 1; i < argc; i += 2)
    {
      // std::cerr << std::setw(2) << i << ": \"" << argv[i] << "\"" <<
      // std::endl;

      auto it = cmd_args.find(argv[i]);
      if (it == cmd_args.end() && arg_alias.count(argv[i]) != 0)
      {
        it = cmd_args.find(arg_alias[argv[i]]);
      }

      std::stringstream error_msg;
      if (it == cmd_args.end())
      {
        error_msg << "Unknown parameter: \"" << argv[i] << "\"";
        throw std::runtime_error(error_msg.str());
      }

      if (!isValidParam(i))
      {
        error_msg << "Parameter \"" << argv[i]
                  << "\" should accept exactly 1 value";
        throw std::runtime_error(error_msg.str());
      }

      FEM::CMD_ARGS param = it->second;

      if (param == FEM::CMD_ARGS::DUMMY)
      {
      }
      else
      {
      }
    }
  }
  void ReadData(fs::path path,
                FEM::IO_FORMAT inputFormat = FEM::IO_FORMAT::CHEKHOV)
  {
    switch (inputFormat)
    {
      case FEM::IO_FORMAT::CHEKHOV:
      {
        std::string inputLine;
        std::string tagName;
        bool isPrevLineTag = false;
        std::ifstream ifs(path);

        if (ifs.fail())
        {
          std::ostringstream oss;
          oss << "Can`t open file " << path;

          throw std::runtime_error(oss.str());
        }

        while (std::getline(ifs, inputLine))
        {
          std::istringstream iss(inputLine);
          iss >> tagName;
          if (tagName == "[End]") break;

          if (tagName == "[Nodes]")
          {
            /*
             * example of input:
             * 1 0. −10. 0. fixed fixed fixed
             */
            uint32_t nodeNumber;
            Point point;
            FEM::DEGREES_OF_FREEDOM freedomDegrees;
          }
          else if (tagName == "[Materials]")
          {
            /*
             * example of input:
             * Al_Alloy 0.0028
             */
            Material material;
            double materialProperty;

            // std::cin >> material >> materialProperty;
          }
          else if (tagName == "[Finite Elements]")
          {
            /*
             * example of input:
             * 1 Quad4 1 2 5 4 Al_Alloy .1
             */
            uint32_t FENumber;
            FEM::FE_TYPE FEType;
            std::vector<uint32_t> nodeNumbers;
            // Material material;
            double parameter;
          }
        }
      }
      break;
      default:
        std::cerr << "Unrecognized format." << std::endl;
    }
  }
  void WriteData(fs::path path,
                 FEM::IO_FORMAT outputFormat = FEM::IO_FORMAT::CHEKHOV) const
  {
  }

 private:
  std::vector<Node> nodePointers;
  std::vector<Material> materialPointers;
  std::vector<std::shared_ptr<FiniteElement>> finiteElements;
};

template class FiniteElementModel<point3d<double>>;

template <class Point>
std::istream& operator>>(std::istream& is, FiniteElementModel<Point>& feModel)
{
  return is;
}

template <class Point>
std::ostream& operator<<(std::ostream& os,
                         const FiniteElementModel<Point>& feModel)
{
  return os;
}

int main(unsigned int argc, const char* argv[])
{
  RunAllTests();

  try
  {
    FiniteElementModel<point3d<double>> FEModel;
    try
    {
      FEModel.SetParams(argc, argv);
    }
    catch (std::exception& e)
    {
      std::cerr << e.what() << std::endl;
      exit(1);
    }

    try
    {
      FEModel.ReadData("data/in/2dim_model_001.txt", FEM::IO_FORMAT::CHEKHOV);
    }
    catch (std::exception& e)
    {
      std::cerr << e.what() << std::endl;
    }

    FEModel.WriteData("data/out/2dim_model_001.txt", FEM::IO_FORMAT::CHEKHOV);
  }
  catch (std::exception& e)
  {
    std::cerr << e.what() << std::endl;
  }
  catch (...)
  {
    std::cerr << "Unknown error occured" << std::endl;
  }

  // std::system("pause");

  return 0;
}