#ifndef GENERATE_MESH_H_
#define GENERATE_MESH_H_
#include "../mesh_utils.h"
#include "vertex_normal.h"
#include "march_vertex.h"
#include "tetrahedralize.h"

namespace generate_mesh {
	class GenerateMesh
	{
	public:
		GenerateMesh(const global_type::Parameter& param,
			global_type::Mesh* mesh) :
			m_param(param),
			m_mesh(mesh)
		{}

		~GenerateMesh() {}

		GenerateMesh(const GenerateMesh& generateMesh) :
			m_mesh(generateMesh.m_mesh),
			m_param(generateMesh.m_param)
		{}

		int generate();
	private:
		global_type::Parameter m_param;
		global_type::Mesh* m_mesh;
		global_type::MeshNormal m_meshNormal;
		global_type::Mesh m_prismTopo;
		global_type::Mesh m_tetTopo;
	};
}

#endif