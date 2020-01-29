#include "stdafx.hpp"

#include "fem.hpp"

shared_ptr<FiniteElementBase> CreateFiniteElement(size_t dim, size_t nodes)
{
	shared_ptr<FiniteElementBase> res;
	switch (dim)
	{
		case 1:
			res = std::make_shared<FiniteElement1d>(nodes);
			break;
		case 2:
			res = std::make_shared<FiniteElement2d>(nodes);
			break;
		case 3:
			res = std::make_shared<FiniteElement3d>(nodes);
			break;
		default:
			break;
	}
	return res;
}

FiniteElementBase::FiniteElementBase()
{
}

FiniteElementBase::FiniteElementBase(size_t nodes) : nodeNumbers(nodes)
{
}

void FiniteElementBase::ReadProperties(std::istream& is)
{
}

FiniteElement1d::FiniteElement1d(size_t nodes) : FiniteElementBase(nodes)
{
}

void FiniteElement1d::ReadProperties(std::istream& is)
{
	is >> crossSectionalArea;
}

FiniteElement2d::FiniteElement2d(size_t nodes) : FiniteElementBase(nodes)
{
}

void FiniteElement2d::ReadProperties(std::istream& is)
{
	is >> thickness;
}

FiniteElement3d::FiniteElement3d(size_t nodes) : FiniteElementBase(nodes)
{
}

template <size_t DIM>
FiniteElementModel<DIM>::FiniteElementModel()
{
	if (DIM > 3)
	{
		throw std::runtime_error(DIM + " dimension is not supported");
	}

	// std::cerr << "You are working in " << DIM << "-dimentional space"
	//          << std::endl;
}

template <size_t DIM>
void FiniteElementModel<DIM>::SetParams(unsigned int argc, const char* argv[])
{
	auto isValidParam = [argc, &argv](unsigned int i) -> bool {
		return i + 1 < argc && (i + 2 >= argc || argv[i + 2][0] == '-');
	};

	for (unsigned int i = 1; i < argc; i += 2)
	{
		// std::cerr << std::setw(2) << i << ": \"" << argv[i] << "\"" <<
		// std::endl;

		auto it = cmd_args.find(argv[i]);
		if (it == std::end(cmd_args) && arg_alias.count(argv[i]) != 0)
		{
			it = cmd_args.find(arg_alias[argv[i]]);
		}

		std::stringstream error_msg;
		if (it == std::end(cmd_args))
		{
			error_msg << "Unknown parameter: \"" << argv[i] << "\"";
			throw std::runtime_error(error_msg.str());
		}

		if (!isValidParam(i))
		{
			error_msg << "Parameter \"" << argv[i] << "\" should accept exactly 1 value";
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

template <size_t DIM>
void FiniteElementModel<DIM>::ReadData(fs::path path, FEM::IO_FORMAT inputFormat)
{
	std::ifstream ifs(path);

	if (ifs.fail())
	{
		std::ostringstream oss;
		oss << "Can`t open file " << path;

		throw std::runtime_error(oss.str());
	}

	switch (inputFormat)
	{
		case FEM::IO_FORMAT::CHEKHOV:
		{
			std::string inputLine;
			FEM::IO::TAG curTag = FEM::IO::TAG::END_OF_INPUT;
			size_t lineNumber = 1;
			while (std::getline(ifs, inputLine))
			{
				std::string_view curTagName(inputLine);
				Trim(curTagName);

				auto it = input_tags.find(std::string(curTagName));

				if (it != std::end(input_tags))
				{
					curTag = it->second;
					if (curTag == FEM::IO::TAG::END_OF_INPUT)
					{
						break;
					}
					continue;
				}
				// the invalid tag is considered to be equal to the input data string
				// else {
				//  std::ostringstream os;
				//  os << "Wrong tag name: \"" << curTagName << "\"";
				//  throw std::runtime_error(os.str());
				//}

				if (curTag == FEM::IO::TAG::END_OF_INPUT)
				{
					break;
				}

				std::istringstream iss(inputLine);

				if (curTag == FEM::IO::TAG::NODES)
				{
					size_t nodeNumber;
					Node<DIM> node;
					iss >> nodeNumber >> node;
					nodes.emplace(nodeNumber, node);
				}
				else if (curTag == FEM::IO::TAG::MATERIALS)
				{
					Material material;
					iss >> material;
					materialCollection.insert(material);
				}
				else if (curTag == FEM::IO::TAG::FINITE_ELEMENTS)
				{
					size_t number;
					std::string type;
					iss >> number >> type;

					auto it = fe_type.find(type);
					if (it == std::end(fe_type))
					{
						std::ostringstream os;
						os << "Wrong FE type: \"" << type << "\" parsed from string \"" << iss.str()
						   << "\"";
						throw std::runtime_error(os.str());
					}

					size_t dim = it->second.first;
					size_t nodes = it->second.second;
					shared_ptr<FiniteElementBase> finiteElement = CreateFiniteElement(dim, nodes);

					try
					{
						ReadFiniteElement(iss, finiteElement);
					}
					catch (std::runtime_error& e)
					{
						std::ostringstream err_msg;
						err_msg << fs::absolute(path).string() << ":" << lineNumber << " "
								<< e.what();
						throw std::runtime_error(err_msg.str());
					}

					finiteElements.emplace(number, finiteElement);
				}
				++lineNumber;
			}
		}
		break;
		default:
			std::cerr << "Unrecognized format." << std::endl;
	}
}

/*
 * example of input:
 * 1 Quad4 1 2 5 4 Al_Alloy .1
 */
template <size_t DIM>
void FiniteElementModel<DIM>::ReadGenericFEEntry(
	std::istream& is, std::shared_ptr<FiniteElementBase> finiteElement)
{
	std::vector<size_t> nodep(finiteElement->GetNodeNumbers().size());
	for (size_t i = 0; i < nodep.size(); ++i)
	{
		is >> nodep[i];
		if (nodes.find(nodep[i]) == std::end(nodes))
		{
			std::ostringstream os;
			os << "Node with number \"" << nodep[i] << "\" does not exist";
			throw std::runtime_error(os.str());
		}
	}

	finiteElement->SetNodeNumbers(std::move(nodep));

	std::string materialName;
	is >> materialName;

	auto it = materialCollection.find(materialName);
	if (it == std::end(materialCollection))
	{
		std::ostringstream os;
		os << "Unknown material \"" << materialName << "\" found during creation of FE";
		throw std::runtime_error(os.str());
	}

	finiteElement->SetMaterial(it);
}

template<>
void FiniteElementModel<3>::ReadSpecFEEntry(
	std::istream& is, std::shared_ptr<FiniteElementBase> finiteElement)
{
	finiteElement->ReadProperties(is);
}

template <size_t DIM>
void FiniteElementModel<DIM>::ReadSpecFEEntry(
	std::istream& is, std::shared_ptr<FiniteElementBase> finiteElement)
{
	finiteElement->ReadProperties(is);
}

template class FiniteElementModel<1>;
template class FiniteElementModel<2>;
template class FiniteElementModel<3>;
