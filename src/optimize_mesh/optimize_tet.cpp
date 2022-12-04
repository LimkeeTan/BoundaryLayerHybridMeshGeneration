#include "optimize_mesh.h"
#include "../mesh_utils.h"
#include "m_slim.h"

namespace optimize_tet {
	void extract_tet_mesh(
		const global_type::Mesh& hybrid_mesh, 
		global_type::Mesh& tet_mesh, 
		std::unordered_map < size_t, size_t >& tet_ver_hybrid_ver_map,
		std::unordered_map < size_t, size_t >& tet_cell_hybrid_cell_map
	)
	{
		std::vector < std::vector < double > > tet_vertex_vec;
		std::vector < std::vector < size_t > > tet_cell_vec;
		std::unordered_map < size_t, size_t > ver_ver_map;
		for (size_t i = 0; i < hybrid_mesh.vecCells.size(); ++i) {
			if (hybrid_mesh.vecCells[i].size() == 4) {
				std::vector < size_t > single_tet(4);
				for (int j = 0; j < 4; ++j) {
					if (!ver_ver_map.count(hybrid_mesh.vecCells[i][j])) {
						std::vector < double > vertex = hybrid_mesh.vecVertices[hybrid_mesh.vecCells[i][j]];
						tet_vertex_vec.emplace_back(vertex);
						ver_ver_map[hybrid_mesh.vecCells[i][j]] = tet_vertex_vec.size() - 1;
						single_tet[j] = tet_vertex_vec.size() - 1;
						tet_ver_hybrid_ver_map[tet_vertex_vec.size() - 1] = hybrid_mesh.vecCells[i][j];
					}
					else {
						single_tet[j] = ver_ver_map[hybrid_mesh.vecCells[i][j]];
					}
				}
				tet_cell_vec.emplace_back(single_tet);
				tet_cell_hybrid_cell_map[tet_cell_vec.size() - 1] = i;
			}
		}
		tet_mesh.vecVertices = tet_vertex_vec;
		tet_mesh.vecCells = tet_cell_vec;
	}

	void construct_target_tet(std::vector < Eigen::MatrixXd >& target_tet)
	{
		Eigen::MatrixXd tagTet;
		tagTet.resize(0, 0);
		for (size_t i = 0; i < target_tet.size(); ++i) {
			target_tet[i] = tagTet;
		}
	}

	int optimize_tet(const global_type::Parameter& param, global_type::Mesh& hybrid_mesh)
	{
		global_type::Mesh tet;
		std::unordered_map < size_t, size_t > tet_cell_hybrid_cell_map;
		tet.vecVertices = hybrid_mesh.vecVertices;
		for (size_t i = 0; i < hybrid_mesh.vecCells.size(); ++i) {
			if (hybrid_mesh.vecCells[i].size() == 4) {
				tet.vecCells.emplace_back(hybrid_mesh.vecCells[i]);
				tet_cell_hybrid_cell_map[tet.vecCells.size() - 1] = i;
			}
		}
		mesh_utils::convertVecToMat(tet);
		Eigen::MatrixXd init_tet_ver = tet.matVertices;
		std::vector < Eigen::MatrixXd > target_tet;
		target_tet.resize(tet.matCells.rows());
		construct_target_tet(target_tet);
		tet.boundaryVerNums = hybrid_mesh.boundaryVerNums;
		std::unordered_map < size_t, Eigen::Vector3d > v_boundary_map;
		slim_opt::tet_optimization(hybrid_mesh, tet, init_tet_ver, target_tet, v_boundary_map);
		for (size_t i = 0; i < tet.matVertices.rows(); ++i) {
			if (v_boundary_map.count(i)) {
				hybrid_mesh.vecVertices[i][0] = v_boundary_map[i][0];
				hybrid_mesh.vecVertices[i][1] = v_boundary_map[i][1];
				hybrid_mesh.vecVertices[i][2] = v_boundary_map[i][2];
				continue;
			}
			hybrid_mesh.vecVertices[i][0] = tet.matVertices(i, 0);
			hybrid_mesh.vecVertices[i][1] = tet.matVertices(i, 1);
			hybrid_mesh.vecVertices[i][2] = tet.matVertices(i, 2);
		}
		Eigen::Vector3d e0; Eigen::Vector3d e1; Eigen::Vector3d e2;
		double h0; double h1; double h2;
		int layerNum = param.layerNumber;
		if(layerNum > 1) {
			double increaseRatio = param.increaseRatio;
			double sum = 0;
			for (int i = 0; i < layerNum; ++i) {
				sum += std::pow(increaseRatio, i);
			}
			std::vector < double > eps0(layerNum);
			std::vector < double > eps1(layerNum);
			std::vector < double > eps2(layerNum);
			double addEps0; double addEps1; double addEps2;
			size_t top0; size_t top1; size_t top2;
			std::unordered_map < size_t, std::vector < size_t > > v_map;
			std::vector < double > marchVertex(3);
			std::vector < double > initVertex(3);
			size_t cellsNum = hybrid_mesh.boundaryCellsNums;
			std::unordered_map < size_t, size_t > v_layer_map;
			for (size_t i = 0; i < cellsNum; ++i) {
				//prism
				if (hybrid_mesh.vecCells[i].size() == 6) {
					if (!v_layer_map.count(hybrid_mesh.vecCells[i][0])) {
						for (int j = 0; j < 3; ++j) {
							e0(j) = hybrid_mesh.vecVertices[hybrid_mesh.vecCells[i][3]][j] - hybrid_mesh.vecVertices[hybrid_mesh.vecCells[i][0]][j];
						}
						h0 = e0.norm();
						if (h0 >= param.userHeight) v_layer_map[hybrid_mesh.vecCells[i][0]] = param.userHeight;
						if (h0 < param.userHeight) {
							double every_layer_height = param.userHeight / layerNum;
							int new_layer = int((h0 / every_layer_height) + global_type::LAYER_THRESHOLD);
							if (new_layer == 0) new_layer = 1;
							v_layer_map[hybrid_mesh.vecCells[i][0]] = new_layer;
						}
					}

					if (!v_layer_map.count(hybrid_mesh.vecCells[i][1])) {
						for (int j = 0; j < 3; ++j) {
							e1(j) = hybrid_mesh.vecVertices[hybrid_mesh.vecCells[i][4]][j] - hybrid_mesh.vecVertices[hybrid_mesh.vecCells[i][1]][j];
						}
						h1 = e1.norm();
						if (h1 >= param.userHeight) v_layer_map[hybrid_mesh.vecCells[i][1]] = param.userHeight;
						if (h1 < param.userHeight) {
							double every_layer_height = param.userHeight / layerNum;
							int new_layer = int((h1 / every_layer_height) + global_type::LAYER_THRESHOLD);
							if (new_layer == 0) new_layer = 1;
							v_layer_map[hybrid_mesh.vecCells[i][1]] = new_layer;
						}
					}

					if (!v_layer_map.count(hybrid_mesh.vecCells[i][2])) {
						for (int j = 0; j < 3; ++j) {
							e2(j) = hybrid_mesh.vecVertices[hybrid_mesh.vecCells[i][5]][j] - hybrid_mesh.vecVertices[hybrid_mesh.vecCells[i][2]][j];
						}
						h2 = e2.norm();
						if (h2 >= param.userHeight) v_layer_map[hybrid_mesh.vecCells[i][2]] = param.userHeight;
						if (h2 < param.userHeight) {
							double every_layer_height = param.userHeight / layerNum;
							int new_layer = int((h2 / every_layer_height) + global_type::LAYER_THRESHOLD);
							if (new_layer == 0) new_layer = 1;
							v_layer_map[hybrid_mesh.vecCells[i][2]] = new_layer;
						}
					}
				}
			}

			for (size_t i = 0; i < cellsNum; ++i) {
				//prism
				if (hybrid_mesh.vecCells[i].size() == 6) {
					if (!v_map.count(hybrid_mesh.vecCells[i][0])) {
						for (int j = 0; j < 3; ++j) {
							e0(j) = hybrid_mesh.vecVertices[hybrid_mesh.vecCells[i][3]][j] - hybrid_mesh.vecVertices[hybrid_mesh.vecCells[i][0]][j];
						}
						h0 = e0.norm(); eps0[0] = h0 / sum; addEps0 = eps0[0];
						for (int j = 1; j < layerNum; ++j) {
							eps0[j] = eps0[0] * std::pow(increaseRatio, j);
							eps0[j] += addEps0;
							addEps0 = eps0[j];
						}
						e0.normalize();
						initVertex = hybrid_mesh.vecVertices[hybrid_mesh.vecCells[i][0]];
						for (int j = 0; j < layerNum - 1; ++j) {
							for (int k = 0; k < 3; ++k) marchVertex[k] = initVertex[k] + eps0[j] * e0(k);
							hybrid_mesh.vecVertices.emplace_back(marchVertex);
							v_map[hybrid_mesh.vecCells[i][0]].emplace_back(hybrid_mesh.vecVertices.size() - 1);
						}
					}

					if (!v_map.count(hybrid_mesh.vecCells[i][1])) {
						for (int j = 0; j < 3; ++j) {
							e1(j) = hybrid_mesh.vecVertices[hybrid_mesh.vecCells[i][4]][j] - hybrid_mesh.vecVertices[hybrid_mesh.vecCells[i][1]][j];
						}
						h1 = e1.norm(); eps1[0] = h1 / sum; addEps1 = eps1[0];
						for (int j = 1; j < layerNum; ++j) {
							eps1[j] = eps1[0] * std::pow(increaseRatio, j);
							eps1[j] += addEps1;
							addEps1 = eps1[j];
						}
						e1.normalize();
						initVertex = hybrid_mesh.vecVertices[hybrid_mesh.vecCells[i][1]];
						for (int j = 0; j < layerNum - 1; ++j) {
							for (int k = 0; k < 3; ++k) marchVertex[k] = initVertex[k] + eps1[j] * e1(k);
							hybrid_mesh.vecVertices.emplace_back(marchVertex);
							v_map[hybrid_mesh.vecCells[i][1]].emplace_back(hybrid_mesh.vecVertices.size() - 1);
						}
					}

					if (!v_map.count(hybrid_mesh.vecCells[i][2])) {
						for (int j = 0; j < 3; ++j) {
							e2(j) = hybrid_mesh.vecVertices[hybrid_mesh.vecCells[i][5]][j] - hybrid_mesh.vecVertices[hybrid_mesh.vecCells[i][2]][j];
						}
						h2 = e2.norm(); eps2[0] = h2 / sum; addEps2 = eps2[0];
						for (int j = 1; j < layerNum; ++j) {
							eps2[j] = eps2[0] * std::pow(increaseRatio, j);
							eps2[j] += addEps2;
							addEps2 = eps2[j];
						}
						e2.normalize();
						initVertex = hybrid_mesh.vecVertices[hybrid_mesh.vecCells[i][2]];
						for (int j = 0; j < layerNum - 1; ++j) {
							for (int k = 0; k < 3; ++k) marchVertex[k] = initVertex[k] + eps2[j] * e2(k);
							hybrid_mesh.vecVertices.emplace_back(marchVertex);
							v_map[hybrid_mesh.vecCells[i][2]].emplace_back(hybrid_mesh.vecVertices.size() - 1);
						}
					}
					top0 = hybrid_mesh.vecCells[i][3]; top1 = hybrid_mesh.vecCells[i][4]; top2 = hybrid_mesh.vecCells[i][5];
					hybrid_mesh.vecCells[i][3] = v_map[hybrid_mesh.vecCells[i][0]][0];
					hybrid_mesh.vecCells[i][4] = v_map[hybrid_mesh.vecCells[i][1]][0];
					hybrid_mesh.vecCells[i][5] = v_map[hybrid_mesh.vecCells[i][2]][0];
					std::vector < size_t > middleCell(6);
					for (int j = 1; j < layerNum - 1; ++j) {
						for (int k = 0; k < 3; ++k) {
							middleCell[k] = v_map[hybrid_mesh.vecCells[i][k]][j - 1];
						}
						for (int k = 0; k < 3; ++k) {
							middleCell[k + 3] = v_map[hybrid_mesh.vecCells[i][k]][j];
						}
						hybrid_mesh.vecCells.emplace_back(middleCell);
					}
					std::vector < size_t > topCell(6);
					topCell[3] = top0; topCell[4] = top1; topCell[5] = top2;
					topCell[0] = v_map[hybrid_mesh.vecCells[i][0]][layerNum - 2];
					topCell[1] = v_map[hybrid_mesh.vecCells[i][1]][layerNum - 2];
					topCell[2] = v_map[hybrid_mesh.vecCells[i][2]][layerNum - 2];
					hybrid_mesh.vecCells.emplace_back(topCell);
				}

				////pyramid
				//if (hybrid_mesh.vecCells[i].size() == 5) {
				//	if (!v_map.count(hybrid_mesh.vecCells[i][0])) {
				//		for (int j = 0; j < 3; ++j) {
				//			e0(j) = hybrid_mesh.vecVertices[hybrid_mesh.vecCells[i][3]][j] - hybrid_mesh.vecVertices[hybrid_mesh.vecCells[i][0]][j];
				//		}
				//		h0 = e0.norm(); eps0[0] = h0 / sum; addEps0 = eps0[0];
				//		for (int j = 1; j < layerNum; ++j) {
				//			eps0[j] = eps0[0] * std::pow(increaseRatio, j);
				//			eps0[j] += addEps0;
				//			addEps0 = eps0[j];
				//		}
				//		e0.normalize();
				//		initVertex = hybrid_mesh.vecVertices[hybrid_mesh.vecCells[i][0]];
				//		for (int j = 0; j < layerNum - 1; ++j) {
				//			for (int k = 0; k < 3; ++k) marchVertex[k] = initVertex[k] + eps0[j] * e0(k);
				//			hybrid_mesh.vecVertices.emplace_back(marchVertex);
				//			v_map[hybrid_mesh.vecCells[i][0]].emplace_back(hybrid_mesh.vecVertices.size() - 1);
				//		}
				//	}

				//	if (!v_map.count(hybrid_mesh.vecCells[i][1])) {
				//		for (int j = 0; j < 3; ++j) {
				//			e1(j) = hybrid_mesh.vecVertices[hybrid_mesh.vecCells[i][2]][j] - hybrid_mesh.vecVertices[hybrid_mesh.vecCells[i][1]][j];
				//		}
				//		h1 = e1.norm(); eps1[0] = h1 / sum; addEps1 = eps1[0];
				//		for (int j = 1; j < layerNum; ++j) {
				//			eps1[j] = eps1[0] * std::pow(increaseRatio, j);
				//			eps1[j] += addEps1;
				//			addEps1 = eps1[j];
				//		}
				//		e1.normalized();
				//		initVertex = hybrid_mesh.vecVertices[hybrid_mesh.vecCells[i][1]];
				//		for (int j = 0; j < layerNum - 1; ++j) {
				//			for (int k = 0; k < 3; ++k) marchVertex[k] = initVertex[k] + eps1[j] * e1(k);
				//			hybrid_mesh.vecVertices.emplace_back(marchVertex);
				//			v_map[hybrid_mesh.vecCells[i][1]].emplace_back(hybrid_mesh.vecVertices.size() - 1);
				//		}
				//	}
				//	top0 = hybrid_mesh.vecCells[i][2]; top1 = hybrid_mesh.vecCells[i][3];
				//	hybrid_mesh.vecCells[i][2] = v_map[hybrid_mesh.vecCells[i][1]][0];
				//	hybrid_mesh.vecCells[i][3] = v_map[hybrid_mesh.vecCells[i][0]][0];
				//	std::vector < size_t > middleCell(5);
				//	for (int j = 1; j < layerNum - 1; ++j) {
				//		middleCell[0] = v_map[hybrid_mesh.vecCells[i][0]][j - 1];
				//		middleCell[1] = v_map[hybrid_mesh.vecCells[i][1]][j - 1];
				//		middleCell[2] = v_map[hybrid_mesh.vecCells[i][1]][j];
				//		middleCell[3] = v_map[hybrid_mesh.vecCells[i][0]][j];
				//		middleCell[4] = hybrid_mesh.vecCells[i][4];
				//		hybrid_mesh.vecCells.emplace_back(middleCell);
				//	}
				//	std::vector < size_t > topCell(5);
				//	topCell[2] = top0; topCell[3] = top1;
				//	topCell[0] = v_map[hybrid_mesh.vecCells[i][0]][layerNum - 2];
				//	topCell[1] = v_map[hybrid_mesh.vecCells[i][1]][layerNum - 2];
				//	topCell[4] = hybrid_mesh.vecCells[i][4];
				//	hybrid_mesh.vecCells.emplace_back(topCell);
				//}
			}
		}
		//mesh_io::saveVTK("data/wanxiangjie_tet.vtk", tet.vecVertices, tet.vecCells);
		return 1;
	}
}