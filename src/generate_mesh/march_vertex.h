#ifndef MARCH_VERTEX_H_
#define MARCH_VERTEX_H_
#include "../global_type.h"
#include "../mesh_utils.h"
#include "../mesh_io/mesh_io.h"

namespace march_vertex {
	int computeMarchVertex(const global_type::Mesh& mesh,
		const global_type::MeshNormal& meshNormal,
		global_type::Parameter& param,
		global_type::Mesh& prismTopo,
		std::unordered_map < size_t, std::vector < size_t > >& vertMap
	);
}

#endif