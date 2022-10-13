#include "optimize_mesh.h"
#include "../mesh_utils.h"
#include "m_slim.h"

namespace optimize_mesh {
	const Eigen::Vector3d& OptimizeMesh::computeTriNorm(const size_t& idx) const
	{
		Eigen::Vector3d v0;
		Eigen::Vector3d v1;
		Eigen::Vector3d v2;
		Eigen::Vector3d e0;
		Eigen::Vector3d e1;
		for (int i = 0; i < 3; ++i) {
			v0(i) = m_mesh->vecVertices[m_mesh->vecCells[idx][0]][i];
			v1(i) = m_mesh->vecVertices[m_mesh->vecCells[idx][1]][i];
			v2(i) = m_mesh->vecVertices[m_mesh->vecCells[idx][2]][i];
		}
		e0 = v1 - v0;
		e1 = v2 - v0;
		return e1.cross(e0).normalized();
	}

	int OptimizeMesh::constructWholeTet(global_type::Mesh& tetMesh, 
		std::vector < Eigen::MatrixXd >& targetPrismTet
	)
	{
		tetMesh.matVertices.resize(m_mesh->vecVertices.size(), 3);
		for (size_t i = 0; i < m_mesh->vecVertices.size(); ++i) {
			for (int j = 0; j < 3; ++j) {
				tetMesh.matVertices(i, j) = m_mesh->vecVertices[i][j];
			}
		}
		std::vector < std::vector < size_t > > tmpTetCell;
		std::vector < size_t > singleTetCell(4);
		size_t triangleNums = 0;
		for (size_t i = 0; i < m_mesh->vecCells.size(); ++i) {
			//triangle
			if (m_mesh->vecCells[i].size() == 3) {
				++triangleNums;
			}

			//prism
			if (m_mesh->vecCells[i].size() == 6) {
				for (int j = 0; j < 6; ++j) {
					for (int k = 0; k < 4; ++k) {
						singleTetCell[k] = m_mesh->vecCells[i][global_type::prismTetCells[j][k]];
					}
					tmpTetCell.emplace_back(singleTetCell);
				}
			}
			//pyramid
			if (m_mesh->vecCells[i].size() == 5) {
				for (int j = 0; j < 4; ++j) {
					for (int k = 0; k < 4; ++k) {
						singleTetCell[k] = m_mesh->vecCells[i][global_type::pyramidTetCells[j][k]];
					}
					tmpTetCell.emplace_back(singleTetCell);
				}
			}
			//tet
			if (m_mesh->vecCells[i].size() == 4) {
				singleTetCell = m_mesh->vecCells[i];
				tmpTetCell.emplace_back(singleTetCell);
			}
		}
		tetMesh.matCells.resize(tmpTetCell.size(), 4);
		for (size_t i = 0; i < tetMesh.matCells.rows(); ++i) {
			for (int j = 0; j < 4; ++j) {
				tetMesh.matCells(i, j) = tmpTetCell[i][j];
			}
		}
		targetPrismTet.resize(tetMesh.matCells.rows());
		for (size_t i = 0; i < targetPrismTet.size(); ++i) {
			targetPrismTet[i].resize(4, 3);
		}
		Eigen::Vector3d triangleNorm;
		for (size_t i = 0; i < triangleNums; ++i) {
			triangleNorm = computeTriNorm(i);
		}
		return 1;
	}

	int OptimizeMesh::optimize()
	{
		global_type::Mesh tetMesh;
		Eigen::MatrixXd initTetVer;
		std::vector < Eigen::MatrixXd > targetPrismTet;
		constructWholeTet(tetMesh, targetPrismTet);
		initTetVer = tetMesh.matVertices;
		tetMesh.boundaryVerNums = m_mesh->boundaryVerNums;
		std::cout << "mesh optimization..." << std::endl;
		slim_opt::slimOptimization(tetMesh, initTetVer);
		for (size_t i = 0; i < tetMesh.matVertices.rows(); ++i) {
			m_mesh->vecVertices[i][0] = tetMesh.matVertices(i, 0);
			m_mesh->vecVertices[i][1] = tetMesh.matVertices(i, 1);
			m_mesh->vecVertices[i][2] = tetMesh.matVertices(i, 2);
		}
		mesh_utils::Jacobian(*m_mesh, "opt_inverse.vtk");
		return 1;
	}
}