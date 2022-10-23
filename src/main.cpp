#include "mesh_io/mesh_io.h"
#include "generate_mesh/generate_mesh.h"
#include "optimize_mesh/optimize_mesh.h"

int main(int argc, char** argv)
{
	global_type::Mesh mesh;
	std::string inputMeshFile = argv[1];
	std::string outputMeshFile = argv[2];
	global_type::Parameter param;
	param.layerNumber = 1;
	param.initHeight = 0.2;
	param.idealHeight = 10.0;
	param.increaseRatio = 1.2;

	if (!mesh_io::readTriOBJ(inputMeshFile, mesh.vecVertices, mesh.vecCells)) {
		std::cout << "failed to read obj mesh" << std::endl;
		return 0;
	}

	generate_mesh::GenerateMesh generateMesh(param, &mesh);
	if (!generateMesh.generate()) {
		std::cout << "failed to generate mesh" << std::endl;
		return 0;
	}

	optimize_mesh::OptimizeMesh optimizeMesh(param, &mesh);
	if (!optimizeMesh.optimize()) {
		std::cout << "failed to optimize mesh" << std::endl;
		return 0;
	}

	//check
	std::vector < std::vector < size_t > > triCells;
	std::vector < size_t > tri(3);
	for (size_t i = 0; i < mesh.vecCells.size(); ++i) {
		if (mesh.vecCells[i].size() == 6) {
			tri[0] = mesh.vecCells[i][3];
			tri[1] = mesh.vecCells[i][4];
			tri[2] = mesh.vecCells[i][5];
			triCells.emplace_back(tri);
		}
	}
	mesh_io::saveTriOBJ("data/topTri.obj", mesh.vecVertices, triCells);

	if (!mesh_io::saveVTK(outputMeshFile, mesh.vecVertices, mesh.vecCells)) {
		std::cout << "failed to save vtk mesh" << std::endl;
		return 0;
	}
	return 0;
}