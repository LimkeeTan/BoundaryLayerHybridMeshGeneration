#include "topo_hybrid.h"

namespace topo_hybrid {
	void computeTriNormal(const std::vector< Eigen::Vector3d >& vertices,
		Eigen::Vector3d& triNormal
	)
	{
		Eigen::Vector3d E0;
		Eigen::Vector3d E1;
		E0 = vertices[1] - vertices[0];
		E1 = vertices[2] - vertices[0];
		triNormal = E0.cross(E1).normalized();
	}

	int validTetCheck(const global_type::Mesh& tetTopo,
		const global_type::MeshNormal& meshNormal,
		const std::unordered_map < size_t, std::vector < size_t > >& vertMap,
		const global_type::Parameter& param,
		global_type::Mesh& prismTopo
	)
	{
		std::unordered_map < size_t, std::vector < size_t > > verCellMap;
		mesh_utils::constructVerCellMap(tetTopo.vecCells, verCellMap);
		std::vector < size_t > relateCell;
		std::vector < size_t > singleCell;
		Eigen::Vector3d triNormal;
		Eigen::Vector3d queryVer;
		Eigen::Vector3d triVer;
		std::vector < Eigen::Vector3d > threeVer(3);
		int triFace[3]{};
		int idx;
		double distance;
		double projection;
		double height = param.initHeight;
		std::vector < double > eps(param.layerNumber);
		double sum = 0;
		double addEps;
		for (int i = 0; i < param.layerNumber; ++i) {
			sum += std::pow(param.increaseRatio, i);
		}
		for (size_t i = 0; i < vertMap.size(); ++i) {
			if (vertMap.at(i).size() == 1) {
				continue;
			}
			for (int j = 0; j < 3; ++j) {
				queryVer[j] = prismTopo.vecVertices[vertMap.at(i).back()][j];
				triVer[j] = prismTopo.vecVertices[i][j];
			}
			relateCell = verCellMap[i];
			for (int j = 0; j < relateCell.size(); ++j) {
				singleCell = tetTopo.vecCells[relateCell[j]];
				for (int k = 0; k < 4; ++k) {
					if (singleCell[k] == i) {
						idx = k;
						break;
					}
				}
				for (int k = 0; k < 3; ++k) {
					triFace[k] = global_type::tetFaces[global_type::tetVerOppositeFace[idx]][k];
					if (vertMap.count(singleCell[triFace[k]]) && vertMap.at(singleCell[triFace[k]]).size() > 1) {
						for (int p = 0; p < 3; ++p) {
							threeVer[k][p] = prismTopo.vecVertices[vertMap.at(singleCell[triFace[k]]).back()][p];
						}
					}
					else {
						for (int p = 0; p < 3; ++p) {
							threeVer[k][p] = tetTopo.vecVertices[singleCell[triFace[k]]][p];
						}
					}
				}
				computeTriNormal(threeVer, triNormal);
				projection = (queryVer - threeVer[0]).dot(triNormal);
				distance = (triNormal * projection).norm();
				while (distance < 0.1 || projection < 0) {
					if ((queryVer - triVer).norm() < 0.001 && projection < 0) {
						std::cout << "inverse tet, tet generation failed" << std::endl;
						return 0;
					}
					if ((queryVer - triVer).norm() < 0.001 && projection > 0) {
						break;
					}
					std::cout << "invalid tet, should reduce marching distance" << std::endl;
					std::cout << distance << std::endl;
					height *= 0.8;
					eps[0] = height / sum;
					addEps = eps[0];
					for (int k = 1; k < param.layerNumber; ++k) {
						eps[k] = eps[0] * std::pow(param.increaseRatio, k);
						eps[k] += addEps;
						addEps = eps[k];
					}
					for (int k = 1; k < vertMap.at(i).size(); ++k) {
						for (int p = 0; p < 3; ++p) {
							prismTopo.vecVertices[vertMap.at(i)[k]][p] = prismTopo.vecVertices[vertMap.at(i)[0]][p] + 
								eps[k - 1] * meshNormal.verticesNormalizedNormal[i][p];
						}
					}
					for (int k = 0; k < 3; ++k) {
						queryVer[k] = prismTopo.vecVertices[vertMap.at(i).back()][k];
					}
					projection = (queryVer - threeVer[0]).dot(triNormal);
					distance = (triNormal * projection).norm();
				}
			}
		}
		return 1;
	}

	int topoOperation(const global_type::Mesh& tetTopo,
		const global_type::MeshNormal& meshNormal,
		global_type::Mesh& prismTopo,
		std::unordered_map < size_t, std::vector < size_t > >& vertMap,
		global_type::Parameter& param,
		global_type::Mesh& hybridMesh
	)
	{
		// valid check
		if (!validTetCheck(tetTopo, meshNormal, vertMap, param, prismTopo)) {
			std::cout << "failed to check tet" << std::endl;
			return 0;
		}

		// combine to obtain hybridMesh
		size_t triVerNum = hybridMesh.vecVertices.size();
		hybridMesh.vecVertices = prismTopo.vecVertices;
		hybridMesh.vecCells = prismTopo.vecCells;
		for (size_t i = triVerNum; i < tetTopo.vecVertices.size(); ++i) {
			hybridMesh.vecVertices.emplace_back(tetTopo.vecVertices[i]);
			vertMap[i].emplace_back(hybridMesh.vecVertices.size() - 1);
		}
		std::vector < size_t > singleTet;
		for (size_t i = 0; i < tetTopo.vecCells.size(); ++i) {
			singleTet = tetTopo.vecCells[i];
			for (int j = 0; j < 4; ++j) {
				singleTet[j] = vertMap.at(singleTet[j]).back();
			}
			hybridMesh.vecCells.emplace_back(singleTet);
		}
		return 1;
	}
}
