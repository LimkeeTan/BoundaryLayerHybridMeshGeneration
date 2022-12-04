#include "optimize_mesh.h"
#include "../mesh_utils.h"
#include "m_slim.h"

namespace optimize_mesh {
	int OptimizeMesh::constructTargetTet(std::vector < Eigen::MatrixXd >& target)
	{
		Eigen::MatrixXd singleTet;
		Eigen::MatrixXd tagTet;
		std::vector < size_t > singleTetCell(4);
		singleTet.resize(4, 3);
		tagTet.resize(0, 0);
		for (size_t i = 0; i < m_mesh->vecCells.size(); ++i) {
			//prism
			if (m_mesh->vecCells[i].size() == 6) {
				for (int j = 0; j < 6; ++j) {
					for (int k = 0; k < 4; ++k) {
						singleTetCell[k] = m_mesh->vecCells[i][global_type::prismSixTetCells[j][k]];
						singleTet(k, 0) = m_mesh->vecVertices[singleTetCell[k]][0];
						singleTet(k, 1) = m_mesh->vecVertices[singleTetCell[k]][1];
						singleTet(k, 2) = m_mesh->vecVertices[singleTetCell[k]][2];
					}
					target.emplace_back(singleTet);
				}
			}
			//pyramid
			if (m_mesh->vecCells[i].size() == 5) {
				for (int j = 0; j < 4; ++j) {
					for (int k = 0; k < 4; ++k) {
						singleTetCell[k] = m_mesh->vecCells[i][global_type::pyramidTetCells[j][k]];
						singleTet(k, 0) = m_mesh->vecVertices[singleTetCell[k]][0];
						singleTet(k, 1) = m_mesh->vecVertices[singleTetCell[k]][1];
						singleTet(k, 2) = m_mesh->vecVertices[singleTetCell[k]][2];
					}
					target.emplace_back(singleTet);
				}
			}
			//tet
			if (m_mesh->vecCells[i].size() == 4) {
				singleTetCell = m_mesh->vecCells[i];
				for (int k = 0; k < 4; ++k) {
					singleTet(k, 0) = m_mesh->vecVertices[singleTetCell[k]][0];
					singleTet(k, 1) = m_mesh->vecVertices[singleTetCell[k]][1];
					singleTet(k, 2) = m_mesh->vecVertices[singleTetCell[k]][2];
				}
				target.emplace_back(singleTet);
			}
		}
		return 1;
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
		for (size_t i = 0; i < m_mesh->vecCells.size(); ++i) {
			//prism
			if (m_mesh->vecCells[i].size() == 6) {
				for (int j = 0; j < 6; ++j) {
					for (int k = 0; k < 4; ++k) {
						singleTetCell[k] = m_mesh->vecCells[i][global_type::prismSixTetCells[j][k]];
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
		targetPrismTet.reserve(tetMesh.matCells.rows());
		constructTargetTet(targetPrismTet);
		return 1;
	}

	int OptimizeMesh::optimize()
	{
		for (int i = 0; i < 2; ++i) {
			global_type::Mesh tetMesh;
			Eigen::MatrixXd initTetVer;
			std::vector < Eigen::MatrixXd > targetPrismTet;
			constructWholeTet(tetMesh, targetPrismTet);
			initTetVer = tetMesh.matVertices;
			tetMesh.boundaryVerNums = m_mesh->boundaryVerNums;
			std::cout << "mesh optimization..." << std::endl;
			slim_opt::slimOptimization(m_param, *m_mesh, tetMesh, initTetVer, targetPrismTet);
			for (size_t j = 0; j < tetMesh.matVertices.rows(); ++j) {
				m_mesh->vecVertices[j][0] = tetMesh.matVertices(j, 0);
				m_mesh->vecVertices[j][1] = tetMesh.matVertices(j, 1);
				m_mesh->vecVertices[j][2] = tetMesh.matVertices(j, 2);
			}
			mesh_utils::Jacobian(*m_mesh, "opt_inverse.vtk");
		}
		return 1;
	}
}