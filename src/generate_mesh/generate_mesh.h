#ifndef GENERATE_MESH_H_
#define GENERATE_MESH_H_
#include "../mesh_utils.h"
#include "vertex_normal.h"
#include "march_vertex.h"

namespace generate_mesh {
	class GenerateMesh
	{
	public:
		GenerateMesh(global_type::Mesh* mesh) :m_mesh(mesh) {}

		~GenerateMesh() {}

		GenerateMesh(const GenerateMesh& generateMesh) :m_mesh(generateMesh.m_mesh) {}

		int generate();
	private:
		global_type::Mesh* m_mesh;
		global_type::MeshNormal m_meshNormal;
		global_type::Mesh m_prismTopo;
	};
}

#endif