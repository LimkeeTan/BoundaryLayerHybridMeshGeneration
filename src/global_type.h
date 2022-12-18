#ifndef GLOBAL_TYPE_H_
#define GLOBAL_TYPE_H_
#include <vector>
#include <Eigen/Eigen>
#include <iostream>
#include <unordered_map>

namespace global_type {
	const double LAYER_THRESHOLD = 0.5;

	struct Mesh {
		std::vector < std::vector < double > > vecVertices;
		std::vector < std::vector < size_t > > vecCells;
		Eigen::MatrixXd matVertices;
		Eigen::MatrixXi matCells;
		size_t boundaryVerNums;
		size_t boundaryCellsNums;
	};

	struct MeshNormal {
		std::vector < size_t > boundary_vertex;
		std::vector < std::vector < double > > cellsNormal;
		std::vector < std::vector < double > > verticesNormalizedNormal;
	};

	struct Parameter {
		int initLayerNumber{};
		int layerNumber{};
		double initHeight{};
		double userHeight{};
		std::vector < double > idealHeight;
		double increaseRatio{};
		std::vector < double > eps;
		std::vector < int > boundary_layer_label;
		std::vector < size_t > boundary_layer_cell;
	};

	const int tetFaces[4][3] =
	{
		{ 0, 1, 2 },
		{ 2, 3, 0 },
		{ 2, 1, 3 },
		{ 1, 0, 3 }
	};

	const int tetVerOppositeFace[4] = { 2,1,3,0 };

	const int prismSixTetCells[6][4] =
	{
		{1, 3, 2, 0},
		{2, 4, 0, 1},
		{0, 4, 1, 2},
		{5, 0, 4, 3},
		{3, 1, 5, 4},
		{4, 2, 3, 5}
	};

	const int prismThreeTetCells[3][4] =
	{
		{0, 1, 2, 3},
		{2, 3, 4, 5},
		{1, 2, 3, 4}
	};

	const int pyramidTetCells[4][4] =
	{
		{1, 4, 3, 0},
		{1, 3, 4, 2},
		{2, 4, 0, 1},
		{4, 2, 0, 3}
	};
}

#endif