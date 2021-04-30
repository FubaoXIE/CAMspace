#include "Entry.h"

#include <igl/PI.h>
#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/boundary_loop.h>
#include <igl/exact_geodesic.h>
#include <igl/gaussian_curvature.h>
#include <igl/invert_diag.h>
#include <igl/lscm.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/embree/ambient_occlusion.h>
#include <igl/file_dialog_open.h>

using namespace std;
// The mesh, Eigen representation
Eigen::MatrixXd meshV;
Eigen::MatrixXi meshF;

void addCurvatureScalar() {
	using namespace Eigen;
	using namespace std;

	VectorXd K;
	igl::gaussian_curvature(meshV, meshF, K);
	SparseMatrix<double> M, Minv;
	igl::massmatrix(meshV, meshF, igl::MASSMATRIX_TYPE_DEFAULT, M);
	igl::invert_diag(M, Minv);
	K = (Minv * K).eval();

	polyscope::getSurfaceMesh("input mesh")
		->addVertexScalarQuantity("gaussian curvature", K, polyscope::DataType::SYMMETRIC);
}


void computeParameterization() {
	using namespace Eigen;
	using namespace std;

	// Fix two points on the boundary
	VectorXi bnd, b(2, 1);
	igl::boundary_loop(meshF, bnd);

	if (bnd.size() == 0) {
		polyscope::warning("mesh has no boundary, cannot parameterize");
		return;
	}

	b(0) = bnd(0);
	b(1) = bnd(round(bnd.size() / 2));
	MatrixXd bc(2, 2);
	bc << 0, 0, 1, 0;

	// LSCM parametrization
	Eigen::MatrixXd V_uv;
	igl::lscm(meshV, meshF, b, bc, V_uv);

	polyscope::getSurfaceMesh("input mesh")->addVertexParameterizationQuantity("LSCM parameterization", V_uv);
}

void computeNormals() {
	Eigen::MatrixXd N_vertices;
	igl::per_vertex_normals(meshV, meshF, N_vertices);
	//std::vector<Eigen::Vector3d> temp_v, temp_n;
	//for (size_t i = 0; i < meshV.rows(); i++)
	//{
	//	Eigen::Vector3d _v = meshV.row(i);
	//	Eigen::Vector3d _n = N_vertices.row(i);
	//	temp_v.push_back(_v);
	//	temp_n.push_back(_n);
	//}
	//CAMGUI::addVectorsArrow("libIGL vertex normals", temp_v, temp_n);
	auto psMesh = polyscope::getSurfaceMesh("input mesh");
	psMesh->setTransparency(0.5);
	auto n = psMesh->addVertexVectorQuantity("libIGL vertex normals", N_vertices);
	n->setVectorColor(glm::vec3(1, 0, 0));
}




void Basic_entry() {

	// static int numPoints = 2000;
	// static float param = 3.14;
	if (ImGui::Button("Import Mesh")) {
		std::string fname = igl::file_dialog_open();
		if (fname.length() == 0)
			return;
		CAMGUI::ImportMesh(fname, meshV, meshF);
		//CAMGUI::init();
	}

	// Curvature
	if (ImGui::Button("add curvature")) {
		addCurvatureScalar();
	}


	// Normals
	if (ImGui::Button("add normals")) {
		computeNormals();
	}

	// Param
	if (ImGui::Button("add parameterization")) {
		computeParameterization();
	}
}

void my_entry() {
	if (ImGui::Button("MY TEST")) {
		my_test();
	}
}

void my_test()
{
	return;
}

