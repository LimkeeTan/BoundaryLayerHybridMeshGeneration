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
		std::unordered_map < size_t, std::vector < size_t > >& verTriMap
	)
	{
		for (size_t i = 0; i < tri.rows(); ++i) {
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
		scaled_jacobian;
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
				std::cout << inverseTet[i] << " " << std::endl;
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
				std::cout << inversePrism[i] << " " << std::endl;
				inverseCell.emplace_back(prismMesh.vecCells[inversePrism[i]]);
			}
		}
		mesh_io::saveVTK(filename, prismMesh.vecVertices, inverseCell);
		return 1;
	}

	int Jacobian(const global_type::Mesh& hybridMesh, const std::string& filename)
	{
		global_type::Mesh prism;
		global_type::Mesh tet;
		tet.vecVertices = hybridMesh.vecVertices;
		for (size_t i = 0; i < hybridMesh.vecCells.size(); ++i) {
			if (hybridMesh.vecCells[i].size() == 4) {
				tet.vecCells.emplace_back(hybridMesh.vecCells[i]);
			}
		}
		convertVecToMat(tet);
		prism.vecVertices = hybridMesh.vecVertices;
		for (size_t i = 0; i < hybridMesh.vecCells.size(); ++i) {
			if (hybridMesh.vecCells[i].size() == 6) {
				prism.vecCells.emplace_back(hybridMesh.vecCells[i]);
			}
		}
		convertVecToMat(prism);
		std::string tet_filename = filename.substr(0, filename.rfind("."));
		std::string prism_filename = filename.substr(0, filename.rfind("."));
		tetJacobian(tet, tet_filename + "_tet.vtk");
		prismJacobian(prism, prism_filename + "_prism.vtk");
		return 1;
	}
}
