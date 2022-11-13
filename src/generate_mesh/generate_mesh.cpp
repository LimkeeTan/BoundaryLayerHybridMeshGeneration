#include "generate_mesh.h"

namespace generate_mesh {
	int GenerateMesh::generate()
	{
		m_mesh.boundaryVerNums = m_mesh.vecVertices.size();
		if (!mesh_utils::convertVecToMat(m_mesh)) {
			std::cout << "failed to convert vec to mat" << std::endl;
			return 0;
		}
		if (!vertex_normal::computeVertexNormal(m_mesh, m_param, m_meshNormal)) {
			std::cout << "failed to compute vertex normal" << std::endl;
			return 0;
		}
		if (!mesh_utils::computeForwardDistance(m_mesh, m_meshNormal, m_param));
		std::unordered_map < size_t, std::vector < size_t > > vertMap;
		if (!march_vertex::computeMarchVertex(m_mesh, m_meshNormal, m_param, m_prismTopo, vertMap)) {
			std::cout << "failed to generate prism topo" << std::endl;
			return 0;
		}
		std::vector < Eigen::Vector3d > holePoints;
		std::vector<size_t> intersectCell;
		if (!tetgen::tetrahedralization(holePoints, m_mesh.vecVertices, m_mesh.vecCells, m_tetTopo.matVertices, m_tetTopo.matCells, intersectCell)) {
			std::cout << "falied to generate tetrahedron topo" << std::endl;
			return 0;
		}
		if (!mesh_utils::convertMatToVec(m_tetTopo)) {
			std::cout << "failed to convert vec to mat" << std::endl;
			return 0;
		}
		if (!topo_hybrid::topoOperation(m_tetTopo, m_meshNormal, m_prismTopo, vertMap, m_param, m_mesh)) {
			std::cout << "failed to apply topo operation" << std::endl;
			return 0;
		}
		//Test
		mesh_io::saveVTK("data/daijinzhijia_hybrid.vtk", m_mesh.vecVertices, m_mesh.vecCells);
		return 1;
	}
}