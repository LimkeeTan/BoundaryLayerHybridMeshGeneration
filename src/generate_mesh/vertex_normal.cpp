#include "vertex_normal.h"

namespace vertex_normal {
	int constructTriNormalMat(const std::vector < std::vector < double > >& triNormal,
		const std::unordered_map < size_t, std::vector < size_t > >& verTriMap,
		Eigen::SparseMatrix < double >& C
	)
	{
		std::vector < Eigen::Triplet < double > > triplet;
		size_t rows = 0;
		size_t cols = verTriMap.size() * 3;
		int faceIdx;
		for (size_t i = 0; i < verTriMap.size(); ++i) {
			size_t iOffset = i * 3;
			for (int j = 0; j < verTriMap.at(i).size(); ++j) {
				faceIdx = verTriMap.at(i)[j];
				triplet.emplace_back(rows, iOffset, triNormal[faceIdx][0]);
				triplet.emplace_back(rows, iOffset + 1, triNormal[faceIdx][1]);
				triplet.emplace_back(rows, iOffset + 2, triNormal[faceIdx][2]);
				++rows;
			}
		}
		C.resize(rows, cols);
		C.setFromTriplets(triplet.begin(), triplet.end());
	}

	void optimization(const Eigen::SparseMatrix<double>& C, const size_t& sizeOfX, Eigen::VectorXd& x)
	{
		Eigen::SparseMatrix<double> H;
		Eigen::VectorXd f;
		Eigen::VectorXd lowerBound;
		Eigen::VectorXd upperBound;

		H.resize(sizeOfX, sizeOfX);
		std::vector<Eigen::Triplet<double>> triplet;
		for (size_t i = 0; i < sizeOfX; ++i)
		{
			triplet.emplace_back(i, i, 2);
		}
		H.setFromTriplets(triplet.begin(), triplet.end());
		f.resize(sizeOfX);
		f.setZero();
		lowerBound.resize(C.rows());
		for (size_t i = 0; i < lowerBound.size(); ++i)
		{
			lowerBound[i] = 1;
		}
		upperBound.resize(C.rows());
		for (size_t i = 0; i < upperBound.size(); ++i)
		{
			upperBound[i] = OsqpEigen::INFTY;
		}
		int NumberOfVariables = C.cols();
		int NumberOfConstraints = C.rows();
		OsqpEigen::Solver solver;
		solver.settings()->setWarmStart(true);
		solver.data()->setNumberOfVariables(NumberOfVariables);
		solver.data()->setNumberOfConstraints(NumberOfConstraints);
		if (!solver.data()->setHessianMatrix(H))
		{
			std::cout << "setHessianMatrix failed\n";
			return;
		}
		if (!solver.data()->setGradient(f))
		{
			std::cout << "setGradient failed\n";
			return;
		}
		if (!solver.data()->setLinearConstraintsMatrix(C))
		{
			std::cout << "setLinearConstraintsMatrix failed\n";
			return;
		}
		if (!solver.data()->setLowerBound(lowerBound))
		{
			std::cout << "setLowerBound failed\n";
			return;
		}
		if (!solver.data()->setUpperBound(upperBound))
		{
			std::cout << "setUpperBound failed\n";
			return;
		}
		// instantiate the solver

		if (!solver.initSolver())
		{
			std::cout << "initSolver failed\n";
			return;
		}

		// solve the QP problem

		if (!solver.solve())
		{
			std::cout << "solve failed\n";
			return;
		}
		clock_t time_start = clock();
		clock_t time_end = clock();
		time_start = clock();
		x = solver.getSolution();
		time_end = clock();
		std::cout << "time use:" << 1000 * (time_end - time_start) / (double)CLOCKS_PER_SEC << "ms" << std::endl;
	}

	int computeVertexNormal(const global_type::Mesh& mesh,
		global_type::MeshNormal& meshNormal
	)
	{
		mesh_utils::computeTriNormal(mesh, meshNormal);
		std::unordered_map < size_t, std::vector < size_t > > verTriMap;
		mesh_utils::constructVerTriMap(mesh.matCells, verTriMap);
		Eigen::SparseMatrix < double > C;
		size_t sizeOfX = mesh.matVertices.rows() * 3;
		constructTriNormalMat(meshNormal.cellsNormal, verTriMap, C);
		Eigen::VectorXd x;
		optimization(C, sizeOfX, x);
		meshNormal.verticesNormal.resize(mesh.matVertices.rows());
		Eigen::Vector3d normal;
		Eigen::Vector3d normalizedNormal;
		std::vector < double > singleNormal(3);
		for (size_t i = 0; i < meshNormal.verticesNormal.size(); ++i)
		{
			normal[0] = x[i * 3];
			normal[1] = x[i * 3 + 1];
			normal[2] = x[i * 3 + 2];
			normalizedNormal = normal.normalized();
			for (int j = 0; j < 3; ++j) {
				singleNormal[j] = normalizedNormal[j];
			}
			meshNormal.verticesNormal[i] = singleNormal;
		}
		return 1;
	}
}