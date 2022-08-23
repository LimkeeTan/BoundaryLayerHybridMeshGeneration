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

	void globalOptimization(const Eigen::SparseMatrix < double >& C, const size_t& sizeOfX, Eigen::VectorXd& x)
	{
		Eigen::SparseMatrix < double > H;
		Eigen::VectorXd f;
		Eigen::VectorXd lowerBound;
		Eigen::VectorXd upperBound;

		H.resize(sizeOfX, sizeOfX);
		std::vector < Eigen::Triplet < double > > triplet;
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
		x = solver.getSolution();
	}

	void localOptimization(const std::unordered_map < size_t, std::vector < size_t > >& verTriMap,
		global_type::MeshNormal& meshNormal
	)
	{
		meshNormal.verticesNormal.resize(verTriMap.size());
		meshNormal.verticesNormalizedNormal.resize(verTriMap.size());
		Eigen::SparseMatrix < double > H;
		Eigen::VectorXd f;
		H.resize(3, 3);
		H.insert(0, 0) = 2;
		H.insert(1, 0) = 0;
		H.insert(2, 0) = 0;
		H.insert(0, 1) = 0;
		H.insert(1, 1) = 2;
		H.insert(2, 1) = 0;
		H.insert(0, 2) = 0;
		H.insert(1, 2) = 0;
		H.insert(2, 2) = 2;
		f.resize(3);
		f.setZero();
		for (size_t i = 0; i < verTriMap.size(); ++i) {
			Eigen::SparseMatrix< double > C;
			Eigen::VectorXd x;
			C.resize(verTriMap.at(i).size(), 3);
			for (int j = 0; j < verTriMap.at(i).size(); ++j) {
				for (int k = 0; k < 3; ++k) {
					C.insert(j, k) = meshNormal.cellsNormal[verTriMap.at(i)[j]][k];
				}
			}
			Eigen::VectorXd lowerBound;
			Eigen::VectorXd upperBound;
			lowerBound.resize(C.rows());
			upperBound.resize(C.rows());
			for (int j = 0; j < lowerBound.size(); ++j) {
				lowerBound[j] = 1;
				upperBound[j] = OsqpEigen::INFTY;
			}
			OsqpEigen::Solver solver;
			solver.settings()->setWarmStart(true);
			solver.data()->setNumberOfVariables(C.cols());
			solver.data()->setNumberOfConstraints(C.rows());
			solver.settings()->setVerbosity(false);
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
			Eigen::Vector3d normalizedNormal;
			if (!solver.solve())
			{
				std::cout << "no solution\n";
				for (int j = 0; j < 3; ++j) {
					meshNormal.verticesNormal[i].emplace_back(0);
					meshNormal.verticesNormalizedNormal[i].emplace_back(0);
				}
			}
			else {
				x = solver.getSolution();
				normalizedNormal = x.normalized();
				for (int j = 0; j < 3; ++j) {
					meshNormal.verticesNormal[i].emplace_back(x[j]);
					meshNormal.verticesNormalizedNormal[i].emplace_back(normalizedNormal[j]);
				}
			}
		}
	}

	bool discardNormal(const std::vector < double >& normal)
	{
		if (std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]) > 100) return true;
		else return false;
	}

	int computeVertexNormal(const global_type::Mesh& mesh,
		global_type::MeshNormal& meshNormal
	)
	{
		mesh_utils::computeTriNormal(mesh, meshNormal);
		std::unordered_map < size_t, std::vector < size_t > > verTriMap;
		mesh_utils::constructVerTriMap(mesh.matCells, verTriMap);
		localOptimization(verTriMap, meshNormal);

		//mesh_utils::computeTriNormal(mesh, meshNormal);
		//std::unordered_map < size_t, std::vector < size_t > > verTriMap;
		//mesh_utils::constructVerTriMap(mesh.matCells, verTriMap);
		//Eigen::SparseMatrix < double > C;
		//size_t sizeOfX = mesh.matVertices.rows() * 3;
		//constructTriNormalMat(meshNormal.cellsNormal, verTriMap, C);
		//Eigen::VectorXd x;
		//globalOptimization(C, sizeOfX, x);
		//meshNormal.verticesNormal.resize(mesh.matVertices.rows());
		//meshNormal.verticesNormalizedNormal.resize(mesh.matVertices.rows());
		//Eigen::Vector3d normal;
		//std::vector < double > singleNormal(3);
		//for (size_t i = 0; i < meshNormal.verticesNormal.size(); ++i)
		//{
		//	normal[0] = x[i * 3];
		//	normal[1] = x[i * 3 + 1];
		//	normal[2] = x[i * 3 + 2];
		//	for (int j = 0; j < 3; ++j) {
		//		singleNormal[j] = normal[j];
		//	}
		//	if (discardNormal(singleNormal)) {
		//		for (int j = 0; j < 3; ++j) singleNormal[j] = 0;
		//		meshNormal.verticesNormal[i] = singleNormal;
		//		meshNormal.verticesNormalizedNormal[i] = singleNormal;
		//	}
		//	else {
		//		meshNormal.verticesNormal[i] = singleNormal;
		//		normal.normalize();
		//		for (int j = 0; j < 3; ++j) singleNormal[j] = normal[j];
		//		meshNormal.verticesNormalizedNormal[i] = singleNormal;
		//	}
		//}

		//Test
		//mesh_io::saveVTK("data/wanxiangjie_global_normal.vtk", mesh.matVertices, mesh.matCells, meshNormal.verticesNormal);
		//mesh_io::saveNormalFile("data/wanxiangjie_global_normal.txt", meshNormal.verticesNormal);
		return 1;
	}
}
