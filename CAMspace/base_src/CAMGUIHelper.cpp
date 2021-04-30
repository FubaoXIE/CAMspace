#include "CAMGUIHelper.h"
#include <igl/readOBJ.h>
#include <igl/readSTL.h>
namespace CAMGUI {
	void init()
	{
		polyscope::view::setUpDir(polyscope::view::UpDir::ZUp);
		polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
		polyscope::options::alwaysRedraw = false;
		polyscope::options::autoscaleStructures = false;
		polyscope::options::autocenterStructures = false;
		polyscope::options::verbosity = 3;
		polyscope::options::ssaaFactor = 3;
		Eigen::MatrixXd _V = Eigen::MatrixXd::Zero(3, 3);
		Eigen::MatrixXi _F(1, 3);
		_F << 0, 1, 2;
		polyscope::registerSurfaceMesh("default", _V, _F);
	}

	void ImportMesh(std::string filename, Eigen::MatrixXd& meshV, Eigen::MatrixXi& meshF) {
		// Read the mesh
		std::string stl_flag = "stl";
		std::string obj_flag = "obj";
		std::string input_format(filename.end() - 3, filename.end());
		if (input_format == obj_flag)
		{
			igl::readOBJ(filename, meshV, meshF);
		}
		
		if (input_format == stl_flag)
		{
			auto meshN = meshV;
			igl::readSTL(filename, meshV, meshF, meshN);
		}
		
		// Register the mesh with Polyscope
		auto p = polyscope::registerSurfaceMesh("input mesh", meshV, meshF);
		polyscope::view::resetCameraToHomeView();
	}




}