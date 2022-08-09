#include "mesh_io/mesh_io.h"
#include "generate_mesh/generate_mesh.h"
#include "optimize_mesh/optimize_mesh.h"

int main(int argc, char** argv)
{
	global_type::Mesh mesh;
	std::string inputMeshFile = "data/wanxiangjie.obj";
	std::string outputMeshFile = "data/wanxiangjie.vtk";

	if (!mesh_io::readTriOBJ(inputMeshFile, mesh.vecVertices, mesh.vecCells)) {
		std::cout << "failed to read obj mesh" << std::endl;
		return 0;
	}

	generate_mesh::GenerateMesh generateMesh(&mesh);
	if (!generateMesh.generate()) {
		std::cout << "failed to generate mesh" << std::endl;
		return 0;
	}

	optimize_mesh::OptimizeMesh optimizeMesh(&mesh);
	if (!optimizeMesh.optimize()) {
		std::cout << "failed to optimize mesh" << std::endl;
		return 0;
	}

	if (!mesh_io::saveVTK(outputMeshFile, mesh.vecVertices, mesh.vecCells)) {
		std::cout << "failed to save vtk mesh" << std::endl;
		return 0;
	}
	return 0;
}
