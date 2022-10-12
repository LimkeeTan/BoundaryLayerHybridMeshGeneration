#include "optimize_mesh.h"
#include "../mesh_utils.h"
#include "m_slim.h"

namespace optimize_mesh {
	int OptimizeMesh::constructWholeTet(global_type::Mesh& tetMesh)
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
		return 1;
	}

	int OptimizeMesh::optimize()
	{
		global_type::Mesh tetMesh;
		Eigen::MatrixXd initTetVer;
		constructWholeTet(tetMesh);
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