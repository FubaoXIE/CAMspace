#pragma once
#include "polyscope/messages.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"

#include <iostream>
#include <unordered_set>
#include <utility>

#include <Eigen/Dense>
#include <Eigen/Geometry>


namespace CAMGUI {

	void init();
	void ImportMesh(std::string filename, Eigen::MatrixXd& meshV, Eigen::MatrixXi& meshF);

	inline void addPoint(std::string name, double x, double y, double z) {
		Eigen::Vector3d p(x, y, z);
		std::vector<Eigen::Vector3d> ps = {p};
		polyscope::registerPointCloud(name, ps);
	};

	template <typename P>
	inline void addPoints(std::string name, const std::vector<P>& points) {
		polyscope::registerPointCloud(name, points);
	};

	template <typename P>
	inline void addLines(std::string name, const std::vector<P>& lines) {
		auto p = polyscope::getSurfaceMesh("default");
		std::vector<std::vector<P>> _tmpLines = { lines };
		p->addSurfaceGraphQuantity(name, _tmpLines);
	};

	template <typename P>
	inline void addLines(std::string name, const std::vector<std::vector<P>>& lines) {
		auto p = polyscope::getSurfaceMesh("default");
		p->addSurfaceGraphQuantity(name, lines);
	};

	template <typename P>
	inline void addLine(std::string name, const P& p0, const P& p1) {
		std::vector<P> _tmpLines = { p0, p1 };
		addLines(name, _tmpLines);
	};


	inline void addFace(std::string name, const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2) {
		Eigen::MatrixXd _V = Eigen::MatrixXd::Zero(3, 3);
		Eigen::MatrixXi _F(1, 3);
		_F << 0, 1, 2;
		_V << p0.x(), p0.y(), p0.z(), p1.x(), p1.y(), p1.z(), p2.x(), p2.y(), p2.z();
		polyscope::registerSurfaceMesh(name, _V, _F);
	};

	template <typename P0, typename P1>
	void addVectors(std::string name, const std::vector<P0>& pos, const std::vector<P1>& vec) {

		auto p = polyscope::getSurfaceMesh("default");
		std::vector<std::vector<P0>> _tmpLines;
		std::vector<P0> _v(2);
		if (pos.size() != vec.size()) {
			std::cout << "The length of pos and vec are not equal" << std::endl;
		}
		for (size_t i; i < pos.size(); i++) {
			_v[0] = pos[i];
			_v[1] = pos[i] + vec[i];
			_tmpLines.push_back(_v);
		}
		p->addSurfaceGraphQuantity(name, _tmpLines);
	};

	template <typename P0, typename P1>
	void addVectorsArrow(std::string name, const std::vector<P0>& pos, const std::vector<P1>& vec) {
		auto p = polyscope::registerPointCloud(name, pos);
		p->addVectorQuantity(name, vec);
	};

}; // namespace CAMGUI