#ifndef BASEUTILITIES_H
#define BASEUTILITIES_H

#include "BASEDefines.h"



	//Basic rotation matrix
	//Rotation along x-axis with _theta
	mAfftransformation mRot_x(double _theta);

	//Rotation along y-axis with _theta
	mAfftransformation mRot_y(double _theta);

	//Rotation along z-axis with _theta
	mAfftransformation mRot_z(double _theta);

	//Check the point if it is in the model or not
	bool is_inModel(mMesh* _Mesh, mRay* _r);

	//Calculate half diagonal line length of a BBox
	double BBoxLength(std::vector<mPoint> _pointlist);

	//Uniform point sampling on the unit sphere
	void UniformSamplingOnSphere(size_t K, SphericalBBOX _SBBox, Eigen::Matrix2Xd& sampling_angles);

	// Check the intersection between two arcs on the unit sphere
	bool is_ArcIntersect(mVector a0, mVector a1, mVector b0, mVector b1);

	//mix product
	double mixproduct(mVector a, mVector b, mVector c);
#endif