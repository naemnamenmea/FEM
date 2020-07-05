#include "ThrowMessage.hpp"
#include "Structure.hpp"
#include "FileInOut.hpp"
#include "cmdLineParams.hpp"
#include <iostream>

int main(int argc, const char* argv[])
{
	if (argc < 2)
	{
		std::cerr << "no input file\n";
		return EXIT_FAILURE;
	}

	try
	{
		CmdLineParams cmdParams(argc, argv);

		//   std::cout << " file name " << argv[1] << std::endl;
		tStructure structure(cmdParams.GetInputFilePath().string().c_str());

		// comment-block begin
		std::cout.precision(15);
		real_t dx = dynamic_cast<const t1D_FE&>(structure.m_FEModel.GetFE(113)).Node(1).Coord()(1) -
					dynamic_cast<const t1D_FE&>(structure.m_FEModel.GetFE(113)).Node(2).Coord()(1),
			   dy = dynamic_cast<const t1D_FE&>(structure.m_FEModel.GetFE(113)).Node(1).Coord()(2) -
					dynamic_cast<const t1D_FE&>(structure.m_FEModel.GetFE(113)).Node(2).Coord()(2);
		std::cout << "Len =\n"
				  << sqrt(dx * dx - dy * dy) << " ===\n"
				  << dynamic_cast<const t1D_FE&>(structure.m_FEModel.GetFE(113)).Length()
				  << std::endl;
		std::cout << structure.m_FEModel.GetFE(113).Volume() << "   - vol" << std::endl;
		dx = dynamic_cast<const t1D_FE&>(structure.m_FEModel.GetFE(168)).Node(1).Coord()(1) -
			 dynamic_cast<const t1D_FE&>(structure.m_FEModel.GetFE(168)).Node(2).Coord()(1),
		dy = dynamic_cast<const t1D_FE&>(structure.m_FEModel.GetFE(168)).Node(1).Coord()(2) -
			 dynamic_cast<const t1D_FE&>(structure.m_FEModel.GetFE(168)).Node(2).Coord()(2);
		std::cout << "Len =\n"
				  << sqrt(dx * dx - dy * dy) << " ===\n"
				  << dynamic_cast<const t1D_FE&>(structure.m_FEModel.GetFE(168)).Length()
				  << std::endl;
		std::cout << structure.m_FEModel.GetFE(168).Volume() << "   - vol" << std::endl;

		real_t x1 = dynamic_cast<const t2D_FE&>(structure.m_FEModel.GetFE(1)).Node(1).Coord()(1),
			   x2 = dynamic_cast<const t2D_FE&>(structure.m_FEModel.GetFE(1)).Node(2).Coord()(1),
			   x3 = dynamic_cast<const t2D_FE&>(structure.m_FEModel.GetFE(1)).Node(3).Coord()(1),
			   x4 = dynamic_cast<const t2D_FE&>(structure.m_FEModel.GetFE(1)).Node(4).Coord()(1),
			   y1 = dynamic_cast<const t2D_FE&>(structure.m_FEModel.GetFE(1)).Node(1).Coord()(2),
			   y2 = dynamic_cast<const t2D_FE&>(structure.m_FEModel.GetFE(1)).Node(2).Coord()(2),
			   y3 = dynamic_cast<const t2D_FE&>(structure.m_FEModel.GetFE(1)).Node(3).Coord()(2),
			   y4 = dynamic_cast<const t2D_FE&>(structure.m_FEModel.GetFE(1)).Node(4).Coord()(2),
			   d1 = sqrt((x3 - x1) * (x3 - x1) + (y3 - y1) * (y3 - y1)),
			   d2 = sqrt((x4 - x2) * (x4 - x2) + (y4 - y2) * (y4 - y2));
		std::cout << "Area =\n"
				  << 0.5 * d1 * d2 << " ===\n"
				  << dynamic_cast<const t2D_FE&>(structure.m_FEModel.GetFE(1)).ComputeArea()
				  << std::endl;
		std::cout << structure.m_FEModel.GetFE(1).Volume() << "   - vol" << std::endl;
		x1 = dynamic_cast<const t2D_FE&>(structure.m_FEModel.GetFE(112)).Node(1).Coord()(1),
		x2 = dynamic_cast<const t2D_FE&>(structure.m_FEModel.GetFE(112)).Node(2).Coord()(1),
		x3 = dynamic_cast<const t2D_FE&>(structure.m_FEModel.GetFE(112)).Node(3).Coord()(1),
		x4 = dynamic_cast<const t2D_FE&>(structure.m_FEModel.GetFE(112)).Node(4).Coord()(1),
		y1 = dynamic_cast<const t2D_FE&>(structure.m_FEModel.GetFE(112)).Node(1).Coord()(2),
		y2 = dynamic_cast<const t2D_FE&>(structure.m_FEModel.GetFE(112)).Node(2).Coord()(2),
		y3 = dynamic_cast<const t2D_FE&>(structure.m_FEModel.GetFE(112)).Node(3).Coord()(2),
		y4 = dynamic_cast<const t2D_FE&>(structure.m_FEModel.GetFE(112)).Node(4).Coord()(2),
		d1 = sqrt((x3 - x1) * (x3 - x1) + (y3 - y1) * (y3 - y1)),
		d2 = sqrt((x4 - x2) * (x4 - x2) + (y4 - y2) * (y4 - y2));
		std::cout << "Area =\n"
				  << 0.5 * d1 * d2 << " ===\n"
				  << dynamic_cast<const t2D_FE&>(structure.m_FEModel.GetFE(112)).ComputeArea()
				  << std::endl;
		std::cout << structure.m_FEModel.GetFE(112).Volume() << "   - vol" << std::endl;
		return 0;
		//     std::cout << structure.FE_model.Volume() << std::endl;
		//     return 0;
		// comment-block end

		std::string dump_file_name(cmdParams.GetInputFilePath().string());
		dump_file_name.replace(dump_file_name.length() - 4, 4, ".dump.txt");
		structure.Dump(dump_file_name.c_str());
		structure.m_analyses.RunCurrent();

		structure.m_analyses.GetCurrent().WriteResultXML();

		// comment-block begin
		// std::cout << "SizeOf Jacobians = "
		//	 << dynamic_cast<const t2D_FE&>(structure.FE_model.FE(1)).Jacobians.size() << std::endl;
		// std::cout << "SizeOf Jacobians ="
		//	 << dynamic_cast<const t2D_FE&>(structure.FE_model.FE(112)).Jacobians.size() <<
		// std::endl; std::cout << "SizeOf Jacobians = "
		//	 << dynamic_cast<const t1D_FE&>(structure.FE_model.FE(113)).Jacobians.size() <<
		// std::endl; std::cout << "SizeOf Jacobians = "
		//	 << dynamic_cast<const t1D_FE&>(structure.FE_model.FE(168)).Jacobians.size() <<
		// std::endl;
		// comment-block end

		std::cout << "FEM executed successfully." << std::endl;
	}
	catch (const tMessage& report_)
	{
		std::cerr << report_ << std::endl;
		return EXIT_FAILURE;
	}
	catch (const std::exception e)
	{
		std::cerr << e.what() << std::endl;
		return EXIT_FAILURE;
	}
	catch (...)
	{
		std::cerr << "Unknown error" << std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
