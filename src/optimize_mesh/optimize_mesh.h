#ifndef OPTIMIZE_MESH_H_
#define OPTIMIZE_MESH_H_
#include "../global_type.h"

namespace optimize_mesh {
	class OptimizeMesh
	{
	public:
		OptimizeMesh(const global_type::Parameter& param,
			global_type::Mesh* mesh) :
			m_mesh(mesh),
			m_param(param)
		{}

		~OptimizeMesh() {}

		OptimizeMesh(const OptimizeMesh& optimizeMesh) :m_mesh(optimizeMesh.m_mesh) {}

		int optimize();
	private:
		global_type::Mesh* m_mesh;
		global_type::Parameter m_param;
		int constructWholeTet(global_type::Mesh& tetMesh);
	};
}

#endif