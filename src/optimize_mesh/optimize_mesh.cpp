#include "optimize_mesh.h"
#include "../mesh_utils.h"
#include "m_slim.h"

namespace optimize_mesh {
	int OptimizeMesh::constructWholeTet(Eigen::MatrixXd& tetVer,
		Eigen::MatrixXi& tetCell
	)
	{
		tetVer.resize(m_mesh->vecVertices.size(), 3);
		for (size_t i = 0; i < m_mesh->vecVertices.size(); ++i) {
			for (int j = 0; j < 3; ++j) {
				tetVer(i, j) = m_mesh->vecVertices[i][j];
			}
		}
		std::vector < std::vector < size_t > > tmpTetCell;
		std::vector < size_t > singleTetCell(4);
		for (size_t i = 0; i < m_mesh->vecCells.size(); ++i) {
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
		tetCell.resize(tmpTetCell.size(), 4);
		for (size_t i = 0; i < tetCell.rows(); ++i) {
			for (int j = 0; j < 4; ++j) {
				tetCell(i, j) = tmpTetCell[i][j];
			}
		}
		return 1;
	}

	int OptimizeMesh::optimize()
	{
		Eigen::MatrixXd tetVer;
		Eigen::MatrixXi tetCell;
		Eigen::MatrixXd initTetVer;
		constructWholeTet(tetVer, tetCell);
		initTetVer = tetVer;
		std::cout << "mesh optimization..." << std::endl;
		slim_opt::slimOptimization(tetVer, tetCell, initTetVer);
		for (size_t i = 0; i < tetVer.rows(); ++i) {
			m_mesh->vecVertices[i][0] = tetVer(i, 0);
			m_mesh->vecVertices[i][1] = tetVer(i, 1);
			m_mesh->vecVertices[i][2] = tetVer(i, 2);
		}
		return 1;
	}
}