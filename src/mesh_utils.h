#ifndef MESH_UTILS_H_
#define MESH_UTILS_H_
#include "global_type.h"
#include "CGAL/Simple_cartesian.h"
#include "CGAL/AABB_tree.h"
#include "CGAL/AABB_traits.h"
#include "CGAL/AABB_face_graph_triangle_primitive.h"
#include "CGAL/AABB_triangle_primitive.h"
#include "CGAL/AABB_segment_primitive.h"

namespace mesh_utils {
	bool isInVector(const std::vector < size_t >& vec,
		const int& value
	);

	bool isInVector(const std::vector < size_t >& vec,
		const size_t& value
	);

	bool isInVector(const std::vector < int >& vec,
		const int& value
	);

	int convertVecToMat(global_type::Mesh& mesh);

	int convertMatToVec(global_type::Mesh& mesh);

	int computeTriNormal(const global_type::Mesh& mesh,
		const std::vector < size_t >& boundary_layer_cell,
		global_type::MeshNormal& meshNormal
	);

	int constructVerTriMap(const Eigen::MatrixXi& tri,
		std::map < size_t, std::vector < size_t > >& verTriMap
	);

	int constructVerTriMap(const Eigen::MatrixXi& tri,
		const std::vector < size_t >& boundary_layer_cell,
		std::map < size_t, std::vector < size_t > >& verTriMap
	);

	int constructVerCellMap(const std::vector < std::vector < size_t > >& cell,
		std::unordered_map < size_t, std::vector < size_t > >& verCellMap
	);

	int tetJacobian(const global_type::Mesh& tetMesh, const std::string& filename);

	int prismJacobian(const global_type::Mesh & prismMesh, const std::string& filename);

	int Jacobian(const global_type::Mesh& hybridMesh, const std::string& filename);

	int scaled_tet_jacobian(const global_type::Mesh& tetMesh);

	int scaled_prism_jacobian(const global_type::Mesh& prismMesh);

	int scaled_jacobian(const global_type::Mesh& hybridMesh);

	int computeForwardDistance(const global_type::Mesh& mesh,
		const global_type::MeshNormal& meshNormal,
		global_type::Parameter& param);

	int prism_quality(const global_type::Mesh& prism_mesh);

	int tet_quality(const global_type::Mesh& tet_mesh);

	int evaluate_quality(const global_type::Mesh& hybridMesh);
}

#endif