#ifndef MESH_UTILS_H_
#define MESH_UTILS_H_
#include "global_type.h"

namespace mesh_utils {
	bool isInVector(const std::vector < size_t >& vec,
		const int& value
	);

	bool isInVector(const std::vector < size_t >& vec,
		const size_t& value
	);

	int convertVecToMat(global_type::Mesh& mesh);

	int convertMatToVec(global_type::Mesh& mesh);

	int computeTriNormal(const global_type::Mesh& mesh,
		global_type::MeshNormal& meshNormal
	);

	int constructVerTriMap(const Eigen::MatrixXi& tri,
		std::unordered_map < size_t, std::vector < size_t > >& verTriMap
	);

	int constructVerCellMap(const std::vector < std::vector < size_t > >& cell,
		std::unordered_map < size_t, std::vector < size_t > >& verCellMap
	);

	int tetJacobian(const global_type::Mesh& tetMesh, const std::string& filename);

	int prismJacobian(const global_type::Mesh & prismMesh, const std::string& filename);

	int Jacobian(const global_type::Mesh& hybridMesh, const std::string& filename);
}

#endif