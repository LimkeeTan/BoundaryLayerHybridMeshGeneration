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
				target.emplace_back(tagTet);
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
		for (int i = 0; i < 1; ++i) {
			global_type::Mesh tetMesh;
			Eigen::MatrixXd initTetVer;
			std::vector < Eigen::MatrixXd > targetPrismTet;
			constructWholeTet(tetMesh, targetPrismTet);
			initTetVer = tetMesh.matVertices;
			tetMesh.boundaryVerNums = m_mesh->boundaryVerNums;
			std::cout << "mesh optimization..." << std::endl;
			slim_opt::slimOptimization(m_param, *m_mesh, tetMesh, initTetVer, targetPrismTet);
			for (size_t i = 0; i < tetMesh.matVertices.rows(); ++i) {
				m_mesh->vecVertices[i][0] = tetMesh.matVertices(i, 0);
				m_mesh->vecVertices[i][1] = tetMesh.matVertices(i, 1);
				m_mesh->vecVertices[i][2] = tetMesh.matVertices(i, 2);
			}
			mesh_utils::Jacobian(*m_mesh, "opt_inverse.vtk");
		}
		Eigen::Vector3d e0; Eigen::Vector3d e1; Eigen::Vector3d e2;
		double h0; double h1; double h2;
		int layerNum = m_param.layerNumber;
		if(layerNum > 1) {
			double increaseRatio = m_param.increaseRatio;
			double sum = 0;
			for (int i = 0; i < layerNum; ++i) {
				sum += std::pow(increaseRatio, i);
			}
			std::vector < double  >eps0(layerNum);
			std::vector < double > eps1(layerNum);
			std::vector < double > eps2(layerNum);
			double addEps0; double addEps1; double addEps2;
			size_t top0; size_t top1; size_t top2;
			std::unordered_map < size_t, std::vector < size_t > > v_map;
			std::vector < double > marchVertex(3);
			std::vector < double > initVertex(3);
			size_t cellsNum = m_mesh->boundaryCellsNums;
			for (size_t i = 0; i < cellsNum; ++i) {
				//prism
				if (m_mesh->vecCells[i].size() == 6) {
					if (!v_map.count(m_mesh->vecCells[i][0])) {
						for (int j = 0; j < 3; ++j) {
							e0(j) = m_mesh->vecVertices[m_mesh->vecCells[i][3]][j] - m_mesh->vecVertices[m_mesh->vecCells[i][0]][j];
						}
						h0 = e0.norm(); eps0[0] = h0 / sum; addEps0 = eps0[0];
						for (int j = 1; j < layerNum; ++j) {
							eps0[j] = eps0[0] * std::pow(increaseRatio, j);
							eps0[j] += addEps0;
							addEps0 = eps0[j];
						}
						e0.normalize();
						initVertex = m_mesh->vecVertices[m_mesh->vecCells[i][0]];
						for (int j = 0; j < layerNum - 1; ++j) {
							for (int k = 0; k < 3; ++k) marchVertex[k] = initVertex[k] + eps0[j] * e0(k);
							m_mesh->vecVertices.emplace_back(marchVertex);
							v_map[m_mesh->vecCells[i][0]].emplace_back(m_mesh->vecVertices.size() - 1);
						}
					}

					if (!v_map.count(m_mesh->vecCells[i][1])) {
						for (int j = 0; j < 3; ++j) {
							e1(j) = m_mesh->vecVertices[m_mesh->vecCells[i][4]][j] - m_mesh->vecVertices[m_mesh->vecCells[i][1]][j];
						}
						h1 = e1.norm(); eps1[0] = h1 / sum; addEps1 = eps1[0];
						for (int j = 1; j < layerNum; ++j) {
							eps1[j] = eps1[0] * std::pow(increaseRatio, j);
							eps1[j] += addEps1;
							addEps1 = eps1[j];
						}
						e1.normalize();
						initVertex = m_mesh->vecVertices[m_mesh->vecCells[i][1]];
						for (int j = 0; j < layerNum - 1; ++j) {
							for (int k = 0; k < 3; ++k) marchVertex[k] = initVertex[k] + eps1[j] * e1(k);
							m_mesh->vecVertices.emplace_back(marchVertex);
							v_map[m_mesh->vecCells[i][1]].emplace_back(m_mesh->vecVertices.size() - 1);
						}
					}

					if (!v_map.count(m_mesh->vecCells[i][2])) {
						for (int j = 0; j < 3; ++j) {
							e2(j) = m_mesh->vecVertices[m_mesh->vecCells[i][5]][j] - m_mesh->vecVertices[m_mesh->vecCells[i][2]][j];
						}
						h2 = e2.norm(); eps2[0] = h2 / sum; addEps2 = eps2[0];
						for (int j = 1; j < layerNum; ++j) {
							eps2[j] = eps2[0] * std::pow(increaseRatio, j);
							eps2[j] += addEps2;
							addEps2 = eps2[j];
						}
						e2.normalize();
						initVertex = m_mesh->vecVertices[m_mesh->vecCells[i][2]];
						for (int j = 0; j < layerNum - 1; ++j) {
							for (int k = 0; k < 3; ++k) marchVertex[k] = initVertex[k] + eps2[j] * e2(k);
							m_mesh->vecVertices.emplace_back(marchVertex);
							v_map[m_mesh->vecCells[i][2]].emplace_back(m_mesh->vecVertices.size() - 1);
						}
					}
					top0 = m_mesh->vecCells[i][3]; top1 = m_mesh->vecCells[i][4]; top2 = m_mesh->vecCells[i][5];
					m_mesh->vecCells[i][3] = v_map[m_mesh->vecCells[i][0]][0];
					m_mesh->vecCells[i][4] = v_map[m_mesh->vecCells[i][1]][0];
					m_mesh->vecCells[i][5] = v_map[m_mesh->vecCells[i][2]][0];
					std::vector < size_t > middleCell(6);
					for (int j = 1; j < layerNum - 1; ++j) {
						for (int k = 0; k < 3; ++k) {
							middleCell[k] = v_map[m_mesh->vecCells[i][k]][j - 1];
						}
						for (int k = 0; k < 3; ++k) {
							middleCell[k + 3] = v_map[m_mesh->vecCells[i][k]][j];
						}
						m_mesh->vecCells.emplace_back(middleCell);
					}
					std::vector < size_t > topCell(6);
					topCell[3] = top0; topCell[4] = top1; topCell[5] = top2;
					topCell[0] = v_map[m_mesh->vecCells[i][0]][layerNum - 2];
					topCell[1] = v_map[m_mesh->vecCells[i][1]][layerNum - 2];
					topCell[2] = v_map[m_mesh->vecCells[i][2]][layerNum - 2];
					m_mesh->vecCells.emplace_back(topCell);
				}

				//pyramid
				if (m_mesh->vecCells[i].size() == 5) {
					if (!v_map.count(m_mesh->vecCells[i][0])) {
						for (int j = 0; j < 3; ++j) {
							e0(j) = m_mesh->vecVertices[m_mesh->vecCells[i][3]][j] - m_mesh->vecVertices[m_mesh->vecCells[i][0]][j];
						}
						h0 = e0.norm(); eps0[0] = h0 / sum; addEps0 = eps0[0];
						for (int j = 1; j < layerNum; ++j) {
							eps0[j] = eps0[0] * std::pow(increaseRatio, j);
							eps0[j] += addEps0;
							addEps0 = eps0[j];
						}
						e0.normalize();
						initVertex = m_mesh->vecVertices[m_mesh->vecCells[i][0]];
						for (int j = 0; j < layerNum - 1; ++j) {
							for (int k = 0; k < 3; ++k) marchVertex[k] = initVertex[k] + eps0[j] * e0(k);
							m_mesh->vecVertices.emplace_back(marchVertex);
							v_map[m_mesh->vecCells[i][0]].emplace_back(m_mesh->vecVertices.size() - 1);
						}
					}

					if (!v_map.count(m_mesh->vecCells[i][1])) {
						for (int j = 0; j < 3; ++j) {
							e1(j) = m_mesh->vecVertices[m_mesh->vecCells[i][2]][j] - m_mesh->vecVertices[m_mesh->vecCells[i][1]][j];
						}
						h1 = e1.norm(); eps1[0] = h1 / sum; addEps1 = eps1[0];
						for (int j = 1; j < layerNum; ++j) {
							eps1[j] = eps1[0] * std::pow(increaseRatio, j);
							eps1[j] += addEps1;
							addEps1 = eps1[j];
						}
						e1.normalized();
						initVertex = m_mesh->vecVertices[m_mesh->vecCells[i][1]];
						for (int j = 0; j < layerNum - 1; ++j) {
							for (int k = 0; k < 3; ++k) marchVertex[k] = initVertex[k] + eps1[j] * e1(k);
							m_mesh->vecVertices.emplace_back(marchVertex);
							v_map[m_mesh->vecCells[i][1]].emplace_back(m_mesh->vecVertices.size() - 1);
						}
					}
					top0 = m_mesh->vecCells[i][2]; top1 = m_mesh->vecCells[i][3];
					m_mesh->vecCells[i][2] = v_map[m_mesh->vecCells[i][1]][0];
					m_mesh->vecCells[i][3] = v_map[m_mesh->vecCells[i][0]][0];
					std::vector < size_t > middleCell(5);
					for (int j = 1; j < layerNum - 1; ++j) {
						middleCell[0] = v_map[m_mesh->vecCells[i][0]][j - 1];
						middleCell[1] = v_map[m_mesh->vecCells[i][1]][j - 1];
						middleCell[2] = v_map[m_mesh->vecCells[i][1]][j];
						middleCell[3] = v_map[m_mesh->vecCells[i][0]][j];
						middleCell[4] = m_mesh->vecCells[i][4];
						m_mesh->vecCells.emplace_back(middleCell);
					}
					std::vector < size_t > topCell(5);
					topCell[2] = top0; topCell[3] = top1;
					topCell[0] = v_map[m_mesh->vecCells[i][0]][layerNum - 2];
					topCell[1] = v_map[m_mesh->vecCells[i][1]][layerNum - 2];
					topCell[4] = m_mesh->vecCells[i][4];
					m_mesh->vecCells.emplace_back(topCell);
				}

				//tet
				if (m_mesh->vecCells[i].size() == 4) {
					if (!v_map.count(m_mesh->vecCells[i][0])) {
						for (int j = 0; j < 3; ++j) {
							e0(j) = m_mesh->vecVertices[m_mesh->vecCells[i][3]][j] - m_mesh->vecVertices[m_mesh->vecCells[i][0]][j];
						}
						h0 = e0.norm(); eps0[0] = h0 / sum; addEps0 = eps0[0];
						for (int j = 1; j < layerNum; ++j) {
							eps0[j] = eps0[0] * std::pow(increaseRatio, j);
							eps0[j] += addEps0;
							addEps0 = eps0[j];
						}
						e0.normalize();
						initVertex = m_mesh->vecVertices[m_mesh->vecCells[i][0]];
						for (int j = 0; j < layerNum - 1; ++j) {
							for (int k = 0; k < 3; ++k) marchVertex[k] = initVertex[k] + eps0[j] * e0(k);
							m_mesh->vecVertices.emplace_back(marchVertex);
							v_map[m_mesh->vecCells[i][0]].emplace_back(m_mesh->vecVertices.size() - 1);
						}
					}
					top0 = m_mesh->vecCells[i][3];
					m_mesh->vecCells[i][3] = v_map[m_mesh->vecCells[i][0]][0];
					std::vector < size_t > middleCell(4);
					for (int j = 1; j < layerNum - 1; ++j) {
						middleCell[0] = v_map[m_mesh->vecCells[i][0]][j - 1];
						middleCell[1] = m_mesh->vecCells[i][1];
						middleCell[2] = m_mesh->vecCells[i][2];
						middleCell[3] = v_map[m_mesh->vecCells[i][0]][j];
						m_mesh->vecCells.emplace_back(middleCell);
					}
					std::vector < size_t > topCell(4);
					topCell[3] = top0;
					topCell[0] = v_map[m_mesh->vecCells[i][0]][layerNum - 2];
					topCell[1] = m_mesh->vecCells[i][1];
					topCell[2] = m_mesh->vecCells[i][2];
					m_mesh->vecCells.emplace_back(topCell);
				}
			}
		}
		return 1;
	}
}