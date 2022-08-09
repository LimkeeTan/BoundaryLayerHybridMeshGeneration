#include "generate_mesh.h"

namespace generate_mesh {
	int GenerateMesh::generate()
	{
		if (!mesh_utils::convertVecToMat(*m_mesh)) {
			std::cout << "failed to convert vec to mat" << std::endl;
			return 0;
		}
		if (!vertex_normal::computeVertexNormal(*m_mesh, m_meshNormal)) {
			std::cout << "failed to compute vertex normal" << std::endl;
			return 0;
		}
		return 1;
	}
} 
