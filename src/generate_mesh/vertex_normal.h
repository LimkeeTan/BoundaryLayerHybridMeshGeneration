#ifndef VERTEX_NORMAL_H_
#define VERTEX_NORMAL_H_
#include "../mesh_utils.h"

namespace vertex_normal {
	int computeVertexNormal(const global_type::Mesh& mesh,
		global_type::MeshNormal& meshNormal
	);
}

#endif
