#include "mesh_utils.h"
#include "mesh_io/mesh_io.h"

namespace mesh_utils {
	bool isInVector(const std::vector < size_t >& vec,
		const int& value
	)
	{
		std::vector < size_t >::const_iterator it = std::find(vec.begin(), vec.end(), value);
		if (it != vec.end()) return true;
		else return false;
	}

	bool isInVector(const std::vector < size_t >& vec,
		const size_t& value
	)
	{
		std::vector < size_t >::const_iterator it = std::find(vec.begin(), vec.end(), value);
		if (it != vec.end()) return true;
		else return false;
	}

	bool isInVector(const std::vector < int >& vec,
		const int& value
	)
	{
		std::vector < int >::const_iterator it = std::find(vec.begin(), vec.end(), value);
		if (it != vec.end()) return true;
		else return false;
	}

	int convertVecToMat(global_type::Mesh& mesh)
	{
		if (mesh.vecVertices.size() == 0) return 1;
		if (mesh.vecCells.size() == 0) return 1;
		int vertexNum = mesh.vecVertices[0].size();
		int cellNum = mesh.vecCells[0].size();
		for (size_t i = 0; i < mesh.vecVertices.size(); ++i) {
			if (mesh.vecVertices[i].size() != vertexNum) {
				std::cout << "vertices error" << std::endl;
				return 0;
			}
		}
		for (size_t i = 0; i < mesh.vecCells.size(); ++i) {
			if (mesh.vecCells[i].size() != cellNum) {
				std::cout << "cells error" << std::endl;
				return 0;
			}
		}
		mesh.matVertices.resize(mesh.vecVertices.size(), vertexNum);
		mesh.matCells.resize(mesh.vecCells.size(), cellNum);
		for (size_t i = 0; i < mesh.vecVertices.size(); ++i) {
			for (int j = 0; j < vertexNum; ++j) {
				mesh.matVertices(i, j) = mesh.vecVertices[i][j];
			}
		}
		for (size_t i = 0; i < mesh.vecCells.size(); ++i) {
			for (int j = 0; j < cellNum; ++j) {
				mesh.matCells(i, j) = mesh.vecCells[i][j];
			}
		}
		return 1;
	}

	int convertMatToVec(global_type::Mesh& mesh)
	{
		if (mesh.matVertices.rows() == 0) return 1;
		if (mesh.matCells.rows() == 0) return 1;
		int vertexNum = mesh.matVertices.rows();
		int cellNum = mesh.matCells.rows();
		mesh.vecVertices.resize(vertexNum);
		mesh.vecCells.resize(cellNum);
		for (size_t i = 0; i < mesh.matVertices.rows(); ++i) {
			mesh.vecVertices[i].emplace_back(mesh.matVertices(i, 0));
			mesh.vecVertices[i].emplace_back(mesh.matVertices(i, 1));
			mesh.vecVertices[i].emplace_back(mesh.matVertices(i, 2));
		}
		for (size_t i = 0; i < mesh.matCells.rows(); ++i) {
			for (int j = 0; j < mesh.matCells.cols(); ++j) {
				mesh.vecCells[i].emplace_back(mesh.matCells(i, j));
			}
		}
		return 1;
	}

	int computeTriNormal(const global_type::Mesh& mesh,
		const std::vector < size_t >& boundary_layer_cell,
		global_type::MeshNormal& meshNormal
	)
	{
		Eigen::Vector3d v0;
		Eigen::Vector3d v1;
		Eigen::Vector3d v2;
		Eigen::Vector3d e0;
		Eigen::Vector3d e1;
		Eigen::Vector3d normalizedNormal;
		std::vector < double > singleNormal(3);
		for (size_t i = 0; i < mesh.matCells.rows(); ++i) {
			v0 = mesh.matVertices.row(mesh.matCells(i, 0));
			v1 = mesh.matVertices.row(mesh.matCells(i, 1));
			v2 = mesh.matVertices.row(mesh.matCells(i, 2));
			e0 = v1 - v0;
			e1 = v2 - v0;
			normalizedNormal = e1.cross(e0).normalized();
			for (int j = 0; j < 3; ++j) {
				singleNormal[j] = normalizedNormal[j];
			}
			meshNormal.cellsNormal.emplace_back(singleNormal);
		}
		return 1;
	}

	int constructVerTriMap(const Eigen::MatrixXi& tri,
		std::map < size_t, std::vector < size_t > >& verTriMap
	)
	{
		for (size_t i = 0; i < tri.rows(); ++i) {
			for (int j = 0; j < 3; ++j) {
				verTriMap[tri(i, j)].emplace_back(i);
			}
		}
		return 1;
	}

	int constructVerTriMap(const Eigen::MatrixXi& tri,
		const std::vector < size_t >& boundary_layer_cell,
		std::map < size_t, std::vector < size_t > >& verTriMap
	)
	{
		for (size_t i = 0; i < tri.rows(); ++i) {
			if (!isInVector(boundary_layer_cell, i)) continue;
			for (int j = 0; j < 3; ++j) {
				verTriMap[tri(i, j)].emplace_back(i);
			}
		}
		return 1;
	}

	int constructVerCellMap(const std::vector < std::vector < size_t > >& cell,
		std::unordered_map < size_t, std::vector < size_t > >& verCellMap
	)
	{
		for (size_t i = 0; i < cell.size(); ++i) {
			for (int j = 0; j < cell[i].size(); ++j) {
				verCellMap[cell[i][j]].emplace_back(i);
			}
		}
		return 1;
	}

	double a_jacobian(const Eigen::Vector3d& v0,
		const Eigen::Vector3d& v1,
		const Eigen::Vector3d& v2,
		const Eigen::Vector3d& v3
	)
	{
		Eigen::Matrix3d Jacobian;

		Jacobian.col(0) = (v1 - v0) * .5;
		Jacobian.col(1) = (v2 - v0) * .5;
		Jacobian.col(2) = (v3 - v0) * .5;


		double norm1 = Jacobian.col(0).norm();
		double norm2 = Jacobian.col(1).norm();
		double norm3 = Jacobian.col(2).norm();

		double scaled_jacobian = Jacobian.determinant();
		if (std::abs(norm1) < 1.e-7 || std::abs(norm2) < 1.e-7 || std::abs(norm3) < 1.e-7) {
			std::cout << "Potential Bug, check!" << std::endl; //system("PAUSE");
			return scaled_jacobian;
		}
		scaled_jacobian /= norm1 * norm2 * norm3;
		return scaled_jacobian;
	}

	int tetJacobian(const global_type::Mesh& tetMesh, const std::string& filename)
	{
		Eigen::Vector3d c0;
		Eigen::Vector3d c1;
		Eigen::Vector3d c2;
		Eigen::Vector3d c3;
		int inverse = 0;
		std::vector < int > inverseTet;
		std::vector < std::vector < size_t > > inverseCell;
		for (size_t i = 0; i < tetMesh.matCells.rows(); ++i) {
			c0 = tetMesh.matVertices.row(tetMesh.matCells(i, 0));
			c1 = tetMesh.matVertices.row(tetMesh.matCells(i, 1));
			c2 = tetMesh.matVertices.row(tetMesh.matCells(i, 2));
			c3 = tetMesh.matVertices.row(tetMesh.matCells(i, 3));

			double jacobianVal = a_jacobian(c0, c1, c2, c3);
			if (jacobianVal < 0) {
				++inverse;
				inverseTet.emplace_back(i);
			}
		}
		if (inverse > 0) {
			std::cout << "flipped elements: " << inverse << std::endl;
			for (int i = 0; i < inverseTet.size(); ++i) {
				//std::cout << inverseTet[i] << " " << std::endl;
				inverseCell.emplace_back(tetMesh.vecCells[inverseTet[i]]);
			}
		}
		mesh_io::saveVTK(filename, tetMesh.vecVertices, inverseCell);
		return 1;
	}

	int prismJacobian(const global_type::Mesh& prismMesh, const std::string& filename)
	{
		Eigen::Vector3d c0;
		Eigen::Vector3d c1;
		Eigen::Vector3d c2;
		Eigen::Vector3d c3;
		int inverse = 0;
		std::vector < int > inversePrism;
		std::vector < std::vector < size_t > > inverseCell;
		for (size_t i = 0; i < prismMesh.matCells.rows(); ++i) {
			for (int j = 0; j < 6; ++j) {
				c0 = prismMesh.matVertices.row(prismMesh.matCells(i, global_type::prismSixTetCells[j][0]));
				c1 = prismMesh.matVertices.row(prismMesh.matCells(i, global_type::prismSixTetCells[j][1]));
				c2 = prismMesh.matVertices.row(prismMesh.matCells(i, global_type::prismSixTetCells[j][2]));
				c3 = prismMesh.matVertices.row(prismMesh.matCells(i, global_type::prismSixTetCells[j][3]));
				double jacobianVal = a_jacobian(c0, c1, c2, c3);
				if (jacobianVal < 0) {
					++inverse;
					inversePrism.emplace_back(i);
					break;
				}
			}
		}
		if (inverse > 0) {
			std::cout << "flipped elements: " << inverse << std::endl;
			for (int i = 0; i < inversePrism.size(); ++i) {
				//std::cout << inversePrism[i] << " " << std::endl;
				inverseCell.emplace_back(prismMesh.vecCells[inversePrism[i]]);
			}
		}
		mesh_io::saveVTK(filename, prismMesh.vecVertices, inverseCell);
		return 1;
	}

	int pyramidJacobian(const global_type::Mesh& pyramidMesh, const std::string& filename)
	{
		Eigen::Vector3d c0;
		Eigen::Vector3d c1;
		Eigen::Vector3d c2;
		Eigen::Vector3d c3;
		int inverse = 0;
		std::vector < int > inversePyramid;
		std::vector < std::vector < size_t > > inverseCell;
		for (size_t i = 0; i < pyramidMesh.matCells.rows(); ++i) {
			for (int j = 0; j < 4; ++j) {
				c0 = pyramidMesh.matVertices.row(pyramidMesh.matCells(i, global_type::pyramidTetCells[j][0]));
				c1 = pyramidMesh.matVertices.row(pyramidMesh.matCells(i, global_type::pyramidTetCells[j][1]));
				c2 = pyramidMesh.matVertices.row(pyramidMesh.matCells(i, global_type::pyramidTetCells[j][2]));
				c3 = pyramidMesh.matVertices.row(pyramidMesh.matCells(i, global_type::pyramidTetCells[j][3]));
				double jacobianVal = a_jacobian(c0, c1, c2, c3);
				if (jacobianVal < 0) {
					++inverse;
					inversePyramid.emplace_back(i);
					break;
				}
			}
		}
		if (inverse > 0) {
			std::cout << "flipped elements: " << inverse << std::endl;
			for (int i = 0; i < inversePyramid.size(); ++i) {
				//std::cout << inversePyramid[i] << " " << std::endl;
				inverseCell.emplace_back(pyramidMesh.vecCells[inversePyramid[i]]);
			}
		}
		mesh_io::saveVTK(filename, pyramidMesh.vecVertices, inverseCell);
		return 1;
	}

	int Jacobian(const global_type::Mesh& hybridMesh, const std::string& filename)
	{
		global_type::Mesh* prism = new global_type::Mesh;
		global_type::Mesh* tet = new global_type::Mesh;
		global_type::Mesh* pyramid = new global_type::Mesh;
		tet->vecVertices = hybridMesh.vecVertices;
		for (size_t i = 0; i < hybridMesh.vecCells.size(); ++i) {
			if (hybridMesh.vecCells[i].size() == 4) {
				tet->vecCells.emplace_back(hybridMesh.vecCells[i]);
			}
		}
		convertVecToMat(*tet);
		std::string tet_filename = filename.substr(0, filename.rfind("."));
		tetJacobian(*tet, tet_filename + "_tet.vtk");
		delete tet;

		prism->vecVertices = hybridMesh.vecVertices;
		for (size_t i = 0; i < hybridMesh.vecCells.size(); ++i) {
			if (hybridMesh.vecCells[i].size() == 6) {
				prism->vecCells.emplace_back(hybridMesh.vecCells[i]);
			}
		}
		convertVecToMat(*prism);
		std::string prism_filename = filename.substr(0, filename.rfind("."));
		prismJacobian(*prism, prism_filename + "_prism.vtk");
		delete prism;

		pyramid->vecVertices = hybridMesh.vecVertices;
		for (size_t i = 0; i < hybridMesh.vecCells.size(); ++i) {
			if (hybridMesh.vecCells[i].size() == 5) {
				pyramid->vecCells.emplace_back(hybridMesh.vecCells[i]);
			}
		}
		convertVecToMat(*pyramid);
		std::string pyramid_filename = filename.substr(0, filename.rfind("."));
		pyramidJacobian(*pyramid, pyramid_filename + "_pyramid.vtk");
		delete pyramid;
		return 1;
	}

	int scaled_tet_jacobian(const global_type::Mesh& tetMesh)
	{
		double ave_jacobian = 0;
		double min_jacobian = 1;
		for (size_t i = 0; i < tetMesh.matCells.rows(); ++i) {
			Eigen::Vector3d c0 = tetMesh.matVertices.row(tetMesh.matCells(i, 0));
			Eigen::Vector3d c1 = tetMesh.matVertices.row(tetMesh.matCells(i, 1));
			Eigen::Vector3d c2 = tetMesh.matVertices.row(tetMesh.matCells(i, 2));
			Eigen::Vector3d c3 = tetMesh.matVertices.row(tetMesh.matCells(i, 3));

			double jacobianVal = a_jacobian(c0, c1, c2, c3);
			ave_jacobian += jacobianVal;
			if (min_jacobian > jacobianVal) min_jacobian = jacobianVal;
		}
		ave_jacobian /= tetMesh.matCells.rows();
		std::cout << "tet MSJ : " << min_jacobian << std::endl;
		std::cout << "tet ASJ : " << ave_jacobian << std::endl;
		return 1;
	}

	int scaled_prism_jacobian(const global_type::Mesh& prismMesh)
	{
		double ave_jacobian = 0;
		double min_jacobian = 1;
		for (size_t i = 0; i < prismMesh.matCells.rows(); ++i) {
			double prism_minJ = 1;
			for (int j = 0; j < 6; ++j) {
				Eigen::Vector3d c0 = prismMesh.matVertices.row(prismMesh.matCells(i, global_type::prismSixTetCells[j][0]));
				Eigen::Vector3d c1 = prismMesh.matVertices.row(prismMesh.matCells(i, global_type::prismSixTetCells[j][1]));
				Eigen::Vector3d c2 = prismMesh.matVertices.row(prismMesh.matCells(i, global_type::prismSixTetCells[j][2]));
				Eigen::Vector3d c3 = prismMesh.matVertices.row(prismMesh.matCells(i, global_type::prismSixTetCells[j][3]));
				double jacobianVal = a_jacobian(c0, c1, c2, c3);
				if (prism_minJ > jacobianVal) prism_minJ = jacobianVal;
			}
			ave_jacobian += prism_minJ;
			if (min_jacobian > prism_minJ) min_jacobian = prism_minJ;
		}
		ave_jacobian /= prismMesh.matCells.rows();
		std::cout << "prism MSJ : " << min_jacobian << std::endl;
		std::cout << "prism ASJ : " << ave_jacobian << std::endl;
		return 1;
	}

	int scaled_jacobian(const global_type::Mesh& hybridMesh)
	{
		global_type::Mesh* tet = new global_type::Mesh;
		tet->vecVertices = hybridMesh.vecVertices;
		for (size_t i = 0; i < hybridMesh.vecCells.size(); ++i) {
			if (hybridMesh.vecCells[i].size() == 4) {
				tet->vecCells.emplace_back(hybridMesh.vecCells[i]);
			}
		}
		convertVecToMat(*tet);
		scaled_tet_jacobian(*tet);
		delete tet;

		global_type::Mesh* prism = new global_type::Mesh;
		prism->vecVertices = hybridMesh.vecVertices;
		for (size_t i = 0; i < hybridMesh.vecCells.size(); ++i) {
			if (hybridMesh.vecCells[i].size() == 6) {
				prism->vecCells.emplace_back(hybridMesh.vecCells[i]);
			}
		}
		convertVecToMat(*prism);
		scaled_prism_jacobian(*prism);
		delete prism;
		return 1;
	}

	int computeForwardDistance(const global_type::Mesh& mesh,
		const global_type::MeshNormal& meshNormal,
		global_type::Parameter& param)
	{
		typedef CGAL::Simple_cartesian<double> K;
		typedef K::Triangle_3 Triangle;
		typedef K::Point_3 Point;
		typedef std::list<Triangle>::iterator TriIterator;
		typedef CGAL::AABB_triangle_primitive<K, TriIterator> TriPrimitive;
		typedef CGAL::AABB_traits<K, TriPrimitive> TriTraits;
		typedef CGAL::AABB_tree<TriTraits> TriTree;
		typedef K::Segment_3 Segment;
		typedef K::Ray_3 Ray;
		typedef K::Vector_3 Vector;
		typedef boost::optional<TriTree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;
		std::vector < double > idealHeight(mesh.vecVertices.size());
		for (size_t i = 0; i < idealHeight.size(); ++i) {
			if (!isInVector(meshNormal.boundary_vertex, i)) {
				idealHeight[i] = -1;
			}
			else {
				idealHeight[i] = param.userHeight;
			}
		}
		std::list < Triangle > triangles;
		std::vector < Point > ps(3);
		for (size_t i = 0; i < mesh.vecCells.size(); ++i) {
			for (int j = 0; j < 3; ++j) {
				ps[j] = { mesh.vecVertices[mesh.vecCells[i][j]][0], mesh.vecVertices[mesh.vecCells[i][j]][1], mesh.vecVertices[mesh.vecCells[i][j]][2] };
			}
			triangles.push_back({ ps[0],ps[1],ps[2] });
		}
		TriTree triTree(triangles.begin(), triangles.end());
		triTree.accelerate_distance_queries();
		if (mesh.vecVertices.size() != mesh.matVertices.rows()) {
			std::cout << "error: vertices vector sizes not equal to matrix rows" << std::endl;
			return 0;
		}
		if (mesh.vecCells.size() != mesh.matCells.rows()) {
			std::cout << "error: cells vector sizes not equal to matrix rows" << std::endl;
			return 0;
		}
		Eigen::Vector3d dir;
		Eigen::Vector3d p;
		Eigen::Vector3d q;
		for (size_t i = 0; i < mesh.matVertices.rows(); ++i) {
			if (!isInVector(meshNormal.boundary_vertex, i)) continue;

			if (meshNormal.verticesNormalizedNormal[i][0] == 0 &&
				meshNormal.verticesNormalizedNormal[i][1] == 0 &&
				meshNormal.verticesNormalizedNormal[i][2] == 0) continue;

			for (int j = 0; j < 3; ++j) {
				dir(j) = meshNormal.verticesNormalizedNormal[i][j];
			}
			for (int j = 0; j < 3; ++j) {
				p(j) = mesh.matVertices(i, j) + 0.01 * dir(j);
				q(j) = mesh.matVertices(i, j) + 3 * idealHeight[i] * dir(j);
			}
			Segment segmentQuery({ p(0),p(1),p(2) }, { q(0),q(1),q(2) });
			if (triTree.do_intersect(segmentQuery)) {
				Ray ray(Point(p(0), p(1), p(2)), Vector(dir(0), dir(1), dir(2)));
				Ray_intersection hit = triTree.first_intersection(ray);
				if (hit) {
					const TriTree::Primitive_id& primitive_id = boost::get < TriTree::Primitive_id >(hit->second);
					Eigen::Matrix3d mat;
					for (int j = 0; j < 3; ++j) {
						mat.col(j) << primitive_id->vertex(j).x(), primitive_id->vertex(j).y(), primitive_id->vertex(j).z();
					}
					int idx;
					(mat.colwise() - p).colwise().squaredNorm().minCoeff(&idx);
					idealHeight[i] = (p - mat.col(idx)).norm() / 3;
				}
			}
		}
		param.idealHeight = idealHeight;
		return 1;
	}

	void compute_volume(Eigen::MatrixXd& L, Eigen::VectorXd& vol)
	{
		const int m = L.rows();
		vol.resize(m, 1);
		for (int t = 0; t < m; t++)
		{
			const double u = L(t, 0);
			const double v = L(t, 1);
			const double w = L(t, 2);
			const double U = L(t, 3);
			const double V = L(t, 4);
			const double W = L(t, 5);
			const double X = (w - U + v) * (U + v + w);
			const double x = (U - v + w) * (v - w + U);
			const double Y = (u - V + w) * (V + w + u);
			const double y = (V - w + u) * (w - u + V);
			const double Z = (v - W + u) * (W + u + v);
			const double z = (W - u + v) * (u - v + W);
			const double a = sqrt(x * Y * Z);
			const double b = sqrt(y * Z * X);
			const double c = sqrt(z * X * Y);
			const double d = sqrt(x * y * z);
			vol(t) = sqrt(
				(-a + b + c + d) *
				(a - b + c + d) *
				(a + b - c + d) *
				(a + b + c - d)) /
				(192. * u * v * w);
		}
	}

	int compute_edge_length(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& L)
	{
		const int m = F.rows();
		switch (F.cols())
		{
		case 2:
		{
			L.resize(F.rows(), 1);
			for (int i = 0; i < F.rows(); ++i) {
				L(i, 0) = (V.row(F(i, 1)) - V.row(F(i, 0))).squaredNorm();
			}
			break;
		}
		case 3:
		{
			L.resize(m, 3);
			for (int i = 0; i < F.rows(); ++i) {
				L(i, 0) = (V.row(F(i, 1)) - V.row(F(i, 2))).squaredNorm();
				L(i, 1) = (V.row(F(i, 2)) - V.row(F(i, 0))).squaredNorm();
				L(i, 2) = (V.row(F(i, 0)) - V.row(F(i, 1))).squaredNorm();
			}
			break;
		}
		case 4:
		{
			L.resize(m, 6);
			for (int i = 0; i < F.rows(); ++i) {
				L(i, 0) = (V.row(F(i, 3)) - V.row(F(i, 0))).squaredNorm();
				L(i, 1) = (V.row(F(i, 3)) - V.row(F(i, 1))).squaredNorm();
				L(i, 2) = (V.row(F(i, 3)) - V.row(F(i, 2))).squaredNorm();
				L(i, 3) = (V.row(F(i, 1)) - V.row(F(i, 2))).squaredNorm();
				L(i, 4) = (V.row(F(i, 2)) - V.row(F(i, 0))).squaredNorm();
				L(i, 5) = (V.row(F(i, 0)) - V.row(F(i, 1))).squaredNorm();
			}
			break;
		}
		}
		L = L.array().sqrt().eval();
		return 1;
	}

	int prism_quality(const global_type::Mesh& prism_mesh)
	{
		Eigen::MatrixXd prism_L;
		Eigen::VectorXd prism_vol;
		prism_L.resize(prism_mesh.matCells.rows(), 9);
		prism_vol.resize(prism_mesh.matCells.rows());
		for (size_t i = 0; i < prism_mesh.matCells.rows(); ++i) {
			Eigen::VectorXd vol1;
			Eigen::VectorXd vol2;
			Eigen::VectorXd vol3;
			Eigen::MatrixXd L1;
			Eigen::MatrixXd L2;
			Eigen::MatrixXd L3;
			Eigen::MatrixXd tet_v_1; Eigen::MatrixXi tet_c_1;
			Eigen::MatrixXd tet_v_2; Eigen::MatrixXi tet_c_2;
			Eigen::MatrixXd tet_v_3; Eigen::MatrixXi tet_c_3;
			tet_v_1.resize(4, 3); tet_c_1.resize(1, 4);
			tet_v_2.resize(4, 3); tet_c_2.resize(1, 4);
			tet_v_3.resize(4, 3); tet_c_3.resize(1, 4);
			std::vector < size_t > c_prism(6);
			for (int j = 0; j < 6; ++j) {
				c_prism[j] = prism_mesh.matCells(i, j);
			}
			tet_v_1.row(0) = prism_mesh.matVertices.row(c_prism[0]);
			tet_v_1.row(1) = prism_mesh.matVertices.row(c_prism[3]);
			tet_v_1.row(2) = prism_mesh.matVertices.row(c_prism[5]);
			tet_v_1.row(3) = prism_mesh.matVertices.row(c_prism[4]);
			tet_c_1(0, 0) = 0; tet_c_1(0, 1) = 1; tet_c_1(0, 2) = 2; tet_c_1(0, 3) = 3;

			tet_v_2.row(0) = prism_mesh.matVertices.row(c_prism[1]);
			tet_v_2.row(1) = prism_mesh.matVertices.row(c_prism[0]);
			tet_v_2.row(2) = prism_mesh.matVertices.row(c_prism[5]);
			tet_v_2.row(3) = prism_mesh.matVertices.row(c_prism[4]);
			tet_c_2(0, 0) = 0; tet_c_2(0, 1) = 1; tet_c_2(0, 2) = 2; tet_c_2(0, 3) = 3;

			tet_v_3.row(0) = prism_mesh.matVertices.row(c_prism[2]);
			tet_v_3.row(1) = prism_mesh.matVertices.row(c_prism[0]);
			tet_v_3.row(2) = prism_mesh.matVertices.row(c_prism[5]);
			tet_v_3.row(3) = prism_mesh.matVertices.row(c_prism[1]);
			tet_c_3(0, 0) = 0; tet_c_3(0, 1) = 1; tet_c_3(0, 2) = 2; tet_c_3(0, 3) = 3;

			compute_edge_length(tet_v_1, tet_c_1, L1); compute_volume(L1, vol1);
			compute_edge_length(tet_v_2, tet_c_2, L2); compute_volume(L2, vol2);
			compute_edge_length(tet_v_3, tet_c_3, L3); compute_volume(L3, vol3);
			prism_vol[i] = vol1[0] + vol2[0] + vol3[0];

			double l1 = (prism_mesh.matVertices.row(c_prism[0]) - prism_mesh.matVertices.row(c_prism[1])).norm();
			double l2 = (prism_mesh.matVertices.row(c_prism[1]) - prism_mesh.matVertices.row(c_prism[2])).norm();
			double l3 = (prism_mesh.matVertices.row(c_prism[2]) - prism_mesh.matVertices.row(c_prism[0])).norm();
			double l4 = (prism_mesh.matVertices.row(c_prism[0]) - prism_mesh.matVertices.row(c_prism[3])).norm();
			double l5 = (prism_mesh.matVertices.row(c_prism[1]) - prism_mesh.matVertices.row(c_prism[4])).norm();
			double l6 = (prism_mesh.matVertices.row(c_prism[2]) - prism_mesh.matVertices.row(c_prism[5])).norm();
			double l7 = (prism_mesh.matVertices.row(c_prism[3]) - prism_mesh.matVertices.row(c_prism[4])).norm();
			double l8 = (prism_mesh.matVertices.row(c_prism[4]) - prism_mesh.matVertices.row(c_prism[5])).norm();
			double l9 = (prism_mesh.matVertices.row(c_prism[5]) - prism_mesh.matVertices.row(c_prism[3])).norm();
			prism_L(i, 0) = l1; prism_L(i, 1) = l2; prism_L(i, 2) = l3; prism_L(i, 3) = l4; prism_L(i, 4) = l5; prism_L(i, 5) = l6;
			prism_L(i, 6) = l7; prism_L(i, 7) = l8; prism_L(i, 8) = l9;
		}

		Eigen::VectorXd prism_quality(prism_vol.size());
		std::vector < size_t > bad_prism_idx;
		for (size_t i = 0; i < prism_quality.size(); ++i) {
			double numerator = 62.3538 * prism_vol[i];
			double edge_length_sum = prism_L(i, 0) * prism_L(i, 0) + prism_L(i, 1) * prism_L(i, 1) + prism_L(i, 2) * prism_L(i, 2) +
				prism_L(i, 3) * prism_L(i, 3) + prism_L(i, 4) * prism_L(i, 4) + prism_L(i, 5) * prism_L(i, 5) + prism_L(i, 6) * prism_L(i, 6) +
				prism_L(i, 7) * prism_L(i, 7) + prism_L(i, 8) * prism_L(i, 8);
			double denominator = pow(edge_length_sum, 1.5);
			prism_quality[i] = numerator / denominator;
			if (prism_quality[i] < 0.1) {
				bad_prism_idx.emplace_back(i);
			}
		}

		double prism_min = prism_quality[0];
		double prism_ave = 0;
		for (size_t i = 0; i < prism_quality.size(); ++i) {
			if (prism_quality[i] < prism_min) prism_min = prism_quality[i];
			prism_ave += prism_quality[i];
		}
		prism_ave /= prism_quality.size();
		std::cout << "prism min quality " << prism_min << std::endl;
		std::cout << "prism ave quality " << prism_ave << std::endl;
		std::vector < std::vector < double > > bad_prism_ver = prism_mesh.vecVertices;
		std::vector < std::vector < size_t > > bad_prism_cell(bad_prism_idx.size());
		for (size_t i = 0; i < bad_prism_cell.size(); ++i) {
			bad_prism_cell[i] = prism_mesh.vecCells[bad_prism_idx[i]];
		}
		mesh_io::saveVTK("data/bad_prism.vtk", bad_prism_ver, bad_prism_cell);
		return 1;
	}

	int tet_quality(const global_type::Mesh& tet_mesh)
	{
		Eigen::MatrixXd tet_L;
		Eigen::MatrixXd tet_v_L;
		Eigen::MatrixXi tet_c_L;
		tet_v_L.resize(tet_mesh.matVertices.rows(), 3);
		tet_c_L.resize(tet_mesh.matCells.rows(), 4);
		for (size_t i = 0; i < tet_v_L.rows(); ++i) {
			tet_v_L.row(i) = tet_mesh.matVertices.row(i);
		}
		for (size_t i = 0; i < tet_c_L.rows(); ++i) {
			tet_c_L.row(i) = tet_mesh.matCells.row(i);
		}
		compute_edge_length(tet_v_L, tet_c_L, tet_L);
		Eigen::VectorXd tet_vol;
		compute_volume(tet_L, tet_vol);

		Eigen::VectorXd tet_quality(tet_vol.size());
		for (size_t i = 0; i < tet_quality.size(); ++i) {
			double numerator = 124.707 * tet_vol[i];
			double edge_length_sum = tet_L(i, 0) * tet_L(i, 0) + tet_L(i, 1) * tet_L(i, 1) + tet_L(i, 2) * tet_L(i, 2) +
				tet_L(i, 3) * tet_L(i, 3) + tet_L(i, 4) * tet_L(i, 4) + tet_L(i, 5) * tet_L(i, 5);
			double denominator = sqrt(pow(edge_length_sum, 3));
			tet_quality[i] = numerator / denominator;
		}
		double tet_min = tet_quality[0];
		double tet_ave = 0;
		for (size_t i = 0; i < tet_quality.size(); ++i) {
			if (tet_quality[i] < tet_min) tet_min = tet_quality[i];
			tet_ave += tet_quality[i];
		}
		tet_ave /= tet_quality.size();
		std::cout << "tet min quality " << tet_min << std::endl;
		std::cout << "tet ave quality " << tet_ave << std::endl;
		return 1;
	}

	int evaluate_quality(const global_type::Mesh& hybridMesh)
	{
		global_type::Mesh* tet = new global_type::Mesh;
		tet->vecVertices = hybridMesh.vecVertices;
		for (size_t i = 0; i < hybridMesh.vecCells.size(); ++i) {
			if (hybridMesh.vecCells[i].size() == 4) {
				tet->vecCells.emplace_back(hybridMesh.vecCells[i]);
			}
		}
		convertVecToMat(*tet);
		tet_quality(*tet);
		delete tet;

		global_type::Mesh* prism = new global_type::Mesh;
		prism->vecVertices = hybridMesh.vecVertices;
		for (size_t i = 0; i < hybridMesh.vecCells.size(); ++i) {
			if (hybridMesh.vecCells[i].size() == 6) {
				prism->vecCells.emplace_back(hybridMesh.vecCells[i]);
			}
		}
		convertVecToMat(*prism);
		prism_quality(*prism);
		delete prism;
		return 1;
	}
}
