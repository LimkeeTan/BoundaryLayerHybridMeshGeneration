#ifndef VERTEX_NORMAL_H_
#define VERTEX_NORMAL_H_
#include "../mesh_utils.h"
#include "../osqp_eigen/OsqpEigen.h"
#include "../mesh_io/mesh_io.h"

namespace vertex_normal {
	int computeVertexNormal(const global_type::Mesh& mesh,
		const global_type::Parameter& param,
		global_type::MeshNormal& meshNormal
	);
}

#endif