#include "topo_hybrid.h"

namespace topo_hybrid {
	void computeTriNormal(const std::vector< std::vector < double > >& vertices,
		std::vector < double >& triNormal
	)
	{
		std::vector < Eigen::Vector3d > V(vertices.size());
		Eigen::Vector3d E0;
		Eigen::Vector3d E1;
		Eigen::Vector3d Normal;
		for (int i = 0; i < vertices.size(); ++i) {
			for (int j = 0; j < 3; ++j) {
				V[i][j] = vertices[i][j];
			}
		}
		E0 = V[1] - V[0]; E1 = V[2] - V[0];
		Normal = E0.cross(E1).normalized();
		for (int i = 0; i < 3; ++i) {
			triNormal[i] = Normal[i];
		}
	}

	int validTetCheck(const global_type::Mesh& tetTopo,
		const global_type::MeshNormal& meshNormal,
		const std::unordered_map < size_t, std::vector < size_t > >& vertMap,
		global_type::Mesh& prismTopo,
		global_type::Parameter& param
	)
	{
		std::unordered_map < size_t, std::vector < size_t > > verCellMap;
		mesh_utils::constructVerCellMap(tetTopo.vecCells, verCellMap);
		std::vector < size_t > relateCell;
		std::vector < size_t > singleCell;
		std::vector < double > triNormal(3);
		std::vector < std::vector < double > > threeVer(3);
		int triFace[3];
		int idx;
		for (size_t i = 0; i < vertMap.size(); ++i) {
			relateCell = verCellMap[i];
			for (int j = 0; j < relateCell.size(); ++j) {
				singleCell = tetTopo.vecCells[relateCell[j]];
				for (int k = 0; k < 4; ++k) {
					if (singleCell[k] == i) {
						idx = k;
						break;
					}
				}
				switch (idx) {
				case 0:
					for (int k = 0; k < 3; ++k) {
						triFace[k] = global_type::tetFaces[2][k];
						threeVer[k] = tetTopo.vecVertices[singleCell[triFace[k]]];
					}
					break;
				case 1:
					for (int k = 0; k < 3; ++k) {
						triFace[k] = global_type::tetFaces[1][k];
						threeVer[k] = tetTopo.vecVertices[singleCell[triFace[k]]];
					}
					break;
				case 2:
					for (int k = 0; k < 3; ++k) {
						triFace[k] = global_type::tetFaces[3][k];
						threeVer[k] = tetTopo.vecVertices[singleCell[triFace[k]]];
					}
					break;
				case 3:
					for (int k = 0; k < 3; ++k) {
						triFace[k] = global_type::tetFaces[0][k];
						threeVer[k] = tetTopo.vecVertices[singleCell[triFace[k]]];
					}
					break;
				default:
					std::cout << "unknow index of tet" << std::endl;
					return 0;
					break;
				}
				computeTriNormal(threeVer, triNormal);
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
		if (!validTetCheck(tetTopo, meshNormal, vertMap, prismTopo, param)) {
			std::cout << "failed to check tet" << std::endl;
			return 0;
		}

		// combine to generate hybridMesh
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
