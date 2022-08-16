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
		global_type::Mesh& prismTopo,
		global_type::Parameter& param
	)
	{
		std::unordered_map < size_t, std::vector < size_t > > verCellMap;
		mesh_utils::constructVerCellMap(tetTopo.vecCells, verCellMap);
		std::vector < size_t > relateCell;
		std::vector < size_t > singleCell;
		Eigen::Vector3d triNormal;
		Eigen::Vector3d queryVer;
		std::vector < Eigen::Vector3d > threeVer(3);
		int triFace[3]{};
		int idx;
		double distance;
		for (size_t i = 0; i < vertMap.size(); ++i) {
			if (vertMap.at(i).size() == 1) {
				continue;
			}
			for (int j = 0; j < 3; ++j) {
				queryVer[j] = prismTopo.vecVertices[vertMap.at(i).back()][j];
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
				distance = (triNormal * (queryVer - threeVer[0]).dot(triNormal)).norm();
				while (distance < 0.1) {
					std::cout << "flatten tet, should reduce marching distance" << std::endl;
					std::cout << relateCell[j] << " " << distance << std::endl;
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
		if (!validTetCheck(tetTopo, meshNormal, vertMap, prismTopo, param)) {
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
