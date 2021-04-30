//Basic defines
#ifndef BASEDEFINES_H
#define BASEDEFINES_H
#define MY_DEBUG

//CGAL related
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/aff_transformation_tags.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/enum.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/jet_smooth_point_set.h>
#include <CGAL/IO/OBJ_reader.h>

//Boost related
#include <boost/foreach.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/undirected_dfs.hpp>

//std
#include <vector>
#include<algorithm>
#include <fstream>
#include <iostream>
#include <chrono>
#include <iomanip>

//Eigen related
#include <Eigen/Dense>
#include <Eigen/Geometry>

//polyscope show function
#include "CAMGUIHelper.h"


#define mPI 3.14159265358979323846
#define ERR 10e-8
#define NozzleRadius 1


typedef CGAL::Simple_cartesian<double>											mKernel;
typedef mKernel::Vector_3														mVector;
typedef mKernel::Vector_3														mNormal;
typedef mKernel::Plane_3														mPlane;
typedef mKernel::Point_3														mPoint;
typedef mKernel::Aff_transformation_3											mAfftransformation;
typedef mKernel::Ray_3															mRay;

typedef CGAL::Surface_mesh<mPoint>												mMesh;
typedef boost::graph_traits<mMesh>::vertex_descriptor							mvertex_descriptor;		//mMesh::Vertex_Index
typedef boost::graph_traits<mMesh>::face_descriptor								mface_descriptor;		//mMesh::Face_Index
typedef boost::graph_traits<mMesh>::edge_descriptor								medge_descriptor;		//mMesh::Edge_Index
typedef boost::graph_traits<mMesh>::halfedge_descriptor							mhalfedge_descriptor;	//mMesh::Halfedge_Index

typedef CGAL::Mean_curvature_flow_skeletonization<mMesh>						mSkeletonization;
typedef mSkeletonization::Skeleton												mSkeleton;
typedef mSkeleton::vertex_descriptor											mSkeleton_vertex;
typedef mSkeleton::edge_descriptor												mSkeleton_edge;
typedef mSkeleton::vertex_iterator												mSkeleton_viterator;

typedef CGAL::AABB_face_graph_triangle_primitive<mMesh>							mPrimitive;
typedef CGAL::AABB_traits<mKernel, mPrimitive>									mTraits;
typedef CGAL::AABB_tree<mTraits>												mAABBtree;
typedef boost::optional<mAABBtree::Intersection_and_primitive_id<mRay>::Type>	mRay_intersection;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
			boost::no_property, 
			boost::property<boost::edge_color_t, boost::default_color_type>>	mGraph;
typedef boost::graph_traits<mGraph>::vertex_descriptor							mGraphvertex;
typedef boost::graph_traits<mGraph>::edge_descriptor							mGraphedge;

typedef std::chrono::high_resolution_clock										mClock;

//Structure for recording the ellipsoid polar and rotation axis
struct SphericalBBOX
{
	double max_theta;
	double min_theta;
	double max_phi;
	double min_phi;
};

#endif