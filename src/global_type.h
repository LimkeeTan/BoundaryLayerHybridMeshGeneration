#ifndef GLOBAL_TYPE_H_
#define GLOBAL_TYPE_H_
#include <vector>
#include <Eigen/Eigen>
#include <iostream>
#include <unordered_map>

namespace global_type {
	struct Mesh {
		std::vector < std::vector < double > > vecVertices;
		std::vector < std::vector < size_t > > vecCells;
		Eigen::MatrixXd matVertices;
		Eigen::MatrixXi matCells;
	};

	struct MeshNormal {
		std::vector < std::vector < double > > verticesNormal;
		std::vector < std::vector < double > > cellsNormal;
		std::vector < std::vector < double > > verticesNormalizedNormal;
	};

	struct Parameter {
		int layerNumber{};
		double initHeight{};
		double increaseRatio{};
		std::vector < double > eps;
	};

	const int tetFaces[4][3] = 
	{
		{ 0, 1, 2 },
		{ 2, 3, 0 },
		{ 2, 1, 3 },
		{ 1, 0, 3 }
	};

	const int tetVerOppositeFace[4] = { 2,1,3,0 };
}

#endif
