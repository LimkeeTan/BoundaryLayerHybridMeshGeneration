#ifndef TOPO_HYBRID_H_
#define TOPO_HYBRID_H_
#include "../global_type.h"
#include "../mesh_utils.h"

namespace topo_hybrid {
	int topoOperation(const global_type::Mesh& tetTopo,
		const global_type::MeshNormal& meshNormal,
		global_type::Mesh& prismTopo,
		std::unordered_map < size_t, std::vector < size_t > >& vertMap,
		global_type::Parameter& param,
		global_type::Mesh& hybridMesh
	);
}

#endif