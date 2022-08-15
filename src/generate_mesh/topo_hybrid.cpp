#include "topo_hybrid.h"

namespace topo_hybrid {
	void validTetCheck(const global_type::Mesh& tetTopo,
		const global_type::MeshNormal& meshNormal,
		const std::unordered_map < size_t, std::vector < size_t > >& vertMap,
		global_type::Mesh& prismTopo,
		global_type::Parameter& param
	)
	{
		std::unordered_map < size_t, std::vector < size_t > > verCellMap;
		mesh_utils::constructVerCellMap(tetTopo.vecCells, verCellMap);

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
		validTetCheck(tetTopo, meshNormal, vertMap, prismTopo, param);

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