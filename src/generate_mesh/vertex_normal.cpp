#include "vertex_normal.h"

namespace vertex_normal {
	void localOptimization(const std::map < size_t, std::vector < size_t > >& verTriMap,
		const size_t& vertices_num,
		global_type::MeshNormal& meshNormal
	)
	{
		meshNormal.verticesNormalizedNormal.resize(vertices_num);
		for (size_t i = 0; i < meshNormal.verticesNormalizedNormal.size(); ++i) {
			meshNormal.verticesNormalizedNormal[i].resize(3, 0);
		}
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
		Eigen::VectorXd lowerBound;
		Eigen::VectorXd upperBound;
		std::map < size_t, std::vector < size_t > >::const_iterator c_it = verTriMap.begin();
		while (c_it != verTriMap.end()) {
			Eigen::SparseMatrix< double > C;
			Eigen::VectorXd x;
			C.resize(c_it->second.size(), 3);
			for (size_t j = 0; j < c_it->second.size(); ++j) {
				for (int k = 0; k < 3; ++k) {
					C.insert(j, k) = meshNormal.cellsNormal[c_it->second[j]][k];
				}
			}
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
			if (!solver.solve())
			{
				std::cout << "no solution\n";
			}
			else {
				x = solver.getSolution().normalized();
				for (int j = 0; j < 3; ++j) {
					meshNormal.verticesNormalizedNormal[c_it->first][j] = x[j];
				}
			}
			++c_it;
		}
	}

	int computeVertexNormal(const global_type::Mesh& mesh,
		const global_type::Parameter& param,
		global_type::MeshNormal& meshNormal
	)
	{
		mesh_utils::computeTriNormal(mesh, param.boundary_layer_cell, meshNormal);
		std::map < size_t, std::vector < size_t > > verTriMap;
		mesh_utils::constructVerTriMap(mesh.matCells, param.boundary_layer_cell, verTriMap);
		std::map < size_t, std::vector < size_t > >::iterator it = verTriMap.begin();
		while (it != verTriMap.end()) {
			meshNormal.boundary_vertex.emplace_back(it->first);
			++it;
		}
		localOptimization(verTriMap, mesh.vecVertices.size(), meshNormal);
		//Test
		//mesh_io::saveVTK("data/new_opt_0_1_2_normal.vtk", mesh.matVertices, mesh.matCells, meshNormal.verticesNormalizedNormal);
		return 1;
	}
}
