#include "stdafx.hpp"

#include "fem.hpp"

template <typename T, size_t DIM>
FiniteElementModel<T, DIM>::FiniteElementModel()
{
#pragma warning(push)
#pragma warning(disable : 4127)
	if (DIM > 3)
#pragma warning(pop)
	{
		throw std::runtime_error(DIM + " dimension is not supported");
	}

	// std::cerr << "You are working in " << DIM << "-dimentional space"
	//          << std::endl;
}

template <typename T, size_t DIM>
void FiniteElementModel<T, DIM>::ReadData(fs::path path, FEM::IO_FORMAT inputFormat)
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

				if (it != input_tags.end())
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
					Node3d<T> node;
					iss >> nodeNumber >> node;
					m_nodes.emplace(nodeNumber, node);
				}
				else if (curTag == FEM::IO::TAG::MATERIALS)
				{
					Material material;
					iss >> material;
					m_materialCollection.insert(material);
				}
				else if (curTag == FEM::IO::TAG::FINITE_ELEMENTS)
				{
					size_t number;
					std::string type;
					iss >> number >> type;

					auto it2 = fe_type.find(type);
					if (it2 == fe_type.end())
					{
						std::ostringstream os;
						os << "Wrong FE type: \"" << type << "\" parsed from string \"" << iss.str()
						   << "\"";
						throw std::runtime_error(os.str());
					}

					size_t dim = it2->second.first;
					size_t nodes = it2->second.second;
					shared_ptr<tFinitElement<T>> finiteElement = CreateFiniteElement<T>(dim, nodes);

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

					m_finiteElements.emplace(number, finiteElement);
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
template <typename T, size_t DIM>
void FiniteElementModel<T, DIM>::ReadGenericFEEntry(
	std::istream& is, std::shared_ptr<tFinitElement<T>> finiteElement)
{
	std::vector<size_t> nodep(finiteElement->GetNodeNumbers().size());
	for (size_t i = 0; i < nodep.size(); ++i)
	{
		is >> nodep[i];
		if (m_nodes.find(nodep[i]) == m_nodes.end())
		{
			std::ostringstream os;
			os << "Node with number \"" << nodep[i] << "\" does not exist";
			throw std::runtime_error(os.str());
		}
	}

	finiteElement->SetNodeNumbers(std::move(nodep));

	std::string materialName;
	is >> materialName;

	auto it = m_materialCollection.find(materialName);
	if (it == m_materialCollection.end())
	{
		std::ostringstream os;
		os << "Unknown material \"" << materialName << "\" found during creation of FE";
		throw std::runtime_error(os.str());
	}

	finiteElement->SetMaterial(it);
}

template <typename T, size_t DIM>
void FiniteElementModel<T, DIM>::ReadSpecFEEntry(
	std::istream& is, std::shared_ptr<tFinitElement<T>> finiteElement)
{
	finiteElement->ReadProperties(is);
}

template class FiniteElementModel<double, 1>;
template class FiniteElementModel<double, 2>;
template class FiniteElementModel<double, 3>;
