
#include "BASEUtilities.h"

//Rotation along x-axis with _theta
mAfftransformation mRot_x(double _theta)
{
	return mAfftransformation(1, 0, 0, 0,
		0, cos(_theta), -sin(_theta), 0,
		0, sin(_theta), cos(_theta), 0,
		1);


}

//Rotation along y-axis with _theta
mAfftransformation mRot_y(double _theta)
{
	return mAfftransformation(cos(_theta), 0, sin(_theta), 0,
		0, 1, 0, 0,
		-sin(_theta), 0, cos(_theta), 0,
		1);

}

//Rotation along z-axis with _theta
mAfftransformation mRot_z(double _theta)
{
	return mAfftransformation(cos(_theta), -sin(_theta), 0, 0,
		sin(_theta), cos(_theta), 0, 0,
		0, 0, 1, 0,
		1);
}

	
bool is_inModel(mMesh* _Mesh, mRay* _r)
{
	mAABBtree tree(faces(*_Mesh).first, faces(*_Mesh).second, *_Mesh);
	int num_interp = tree.number_of_intersected_primitives(*_r);
	//std::cout << tree.size() <<std::endl;
	//std::cout << num_interp <<std::endl;

	if (num_interp % 2 == 0)
		return false;
	else
		return true;

}


double BBoxLength(std::vector<mPoint> _pointlist)
{
	int num_p = _pointlist.size();
	Eigen::Matrix3Xf PointMat = Eigen::Matrix3Xf::Zero(3, num_p);

	for (size_t i = 0; i < num_p; i++)
	{
		Eigen::Vector3f _p(_pointlist[i].x(), _pointlist[i].y(), _pointlist[i].z());
		PointMat.col(i) = _p;
	}

	Eigen::Vector3f bboxMin = (PointMat.rowwise().minCoeff()).transpose();
	Eigen::Vector3f bboxMax = (PointMat.rowwise().maxCoeff()).transpose();

	double diagnal_length = (bboxMax - bboxMin).norm() / 2.0;

	return diagnal_length;
}

void UniformSamplingOnSphere(size_t K, SphericalBBOX _SBBox, Eigen::Matrix2Xd& sampling_angles)
{
	sampling_angles = Eigen::Matrix2Xd::Zero(2, K);
	double diff_theta = _SBBox.max_theta - _SBBox.min_theta;
	double diff_phi = _SBBox.max_phi - _SBBox.min_phi;

	double num_up = K * diff_theta*diff_theta;
	double num_down = 2 * diff_phi*sin(diff_theta / 2)*cos(_SBBox.min_theta + diff_theta / 2);
	double num = ceil(sqrt(num_up / num_down));

	double L = diff_phi * sin(diff_theta / 2)*cos(_SBBox.min_theta + diff_theta / 2) / sin(0.5*diff_theta / num);
		
	rsize_t sum = 0;
	size_t index = 0;

	for (size_t i = 0; i < num-1; i++)
	{
		double theta = _SBBox.min_theta + (i + 0.5)*diff_theta / num;
		double k = round(diff_phi*cos(theta)*K / L);
		sum = sum + k;
		if (sum > K)
		{
			sum = sum - k;
			k = K - sum;
		}

		for (size_t j = 0; j < k; j++)
		{
			double phi = _SBBox.min_phi + (j + 0.5)*diff_phi / k;
			Eigen::Vector2d tmp_angle(theta, phi);
			sampling_angles.col(index) = tmp_angle;
			index++;
		}

	}

	double theta = _SBBox.min_theta + (num - 0.5)*diff_theta / num;
	double k = K - sum;
	for (size_t j = 0; j < k; j++)
	{
		double phi = _SBBox.min_phi + (j + 0.5)*diff_phi / k;
		Eigen::Vector2d tmp_angle(theta, phi);
		sampling_angles.col(index) = tmp_angle;
		index++;
	}

}


bool is_ArcIntersect(mVector a0, mVector a1, mVector b0, mVector b1 )
{
	mVector p = CGAL::cross_product(a0, a1);
	mVector q = CGAL::cross_product(b0, b1);

	mVector t = CGAL::cross_product(p, q);
	t = t / sqrt(t.squared_length());

	//double s1 = CGAL::cross_product(a0, p)*t;
	//double s2 = CGAL::cross_product(a1, p)*t;
	//double s3 = CGAL::cross_product(b0, q)*t;
	//double s4 = CGAL::cross_product(b1, q)*t;
	double s1 = mixproduct(t, a0, p);
	double s2 = mixproduct(t, a1, p);
	double s3 = mixproduct(t, b0, q);
	double s4 = mixproduct(t, b1, q);

	if ((((-s1) > 0) && (s2 > 0) && ((-s3) > 0) && (s4 > 0)) || (((-s1) < 0) && (s2 < 0) && ((-s3) < 0) && (s4 < 0)))
		return true;
	else
		return false;

}
//a*(b x c)
double mixproduct(mVector a, mVector b, mVector c)
{
	Eigen::Matrix3d _m;

	_m << a.x(), a.y(), a.z(),
			b.x(), b.y(), b.z(),
			c.x(), c.y(), c.z();
		
	return _m.determinant();
}

