#include "mesh_utils.h"

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

	int computeTriNormal(const global_type::Mesh& mesh,
		global_type::MeshNormal& meshNormal
	)
	{
		Eigen::Vector3d v0;
		Eigen::Vector3d v1;
		Eigen::Vector3d v2;
		Eigen::Vector3d e0;
		Eigen::Vector3d e1;
		Eigen::Vector3d normal;
		Eigen::Vector3d normalizedNormal;
		std::vector < double > singleNormal(3);
		for (size_t i = 0; i < mesh.matCells.rows(); ++i) {
			v0 = mesh.matVertices.row(mesh.matCells(i, 0));
			v1 = mesh.matVertices.row(mesh.matCells(i, 1));
			v2 = mesh.matVertices.row(mesh.matCells(i, 2));
			e0 = v1 - v0;
			e1 = v2 - v0;
			normal = e0.cross(e1);
			normalizedNormal = normal.normalized();
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
}