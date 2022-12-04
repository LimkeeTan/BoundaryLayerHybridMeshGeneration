#include "march_vertex.h"

namespace march_vertex {
	bool noNormal(const std::vector < double >& normal)
	{
		if (normal[0] == 0 && normal[1] == 0 && normal[2] == 0) return true;
		else return false;
	}

	int generateVertices(const std::vector < std::vector < double > >& vertices,
		const std::vector < size_t >& boundary_vertex,
		const std::vector < std::vector < double > >& normalizedNorm,
		const int& layerNum,
		const std::vector < double >& eps,
		global_type::Mesh& prismTopo,
		std::unordered_map < size_t, std::vector < size_t > >& vertMap
	)
	{
		std::vector < double > initVertex(3);
		std::vector < double > marchVertex(3);
		for (size_t i = 0; i < vertices.size(); ++i) {
			if (!mesh_utils::isInVector(boundary_vertex, i)) continue;
			if (noNormal(normalizedNorm[i])) continue;
			initVertex = vertices[i];
			vertMap[i].emplace_back(i);
			for (int j = 0; j < layerNum; ++j) {
				for (int k = 0; k < 3; ++k) {
					marchVertex[k] = initVertex[k] + eps[j] * normalizedNorm[i][k];
				}
				prismTopo.vecVertices.emplace_back(marchVertex);
				vertMap[i].emplace_back(prismTopo.vecVertices.size() - 1);
			}
		}
		return 1;
	}

	int generateCell(const std::vector < std::vector < size_t > >& cells,
		const std::unordered_map < size_t, std::vector < size_t > >& vertMap,
		const std::vector < size_t >& boundary_layer_cell,
		global_type::Mesh& prismTopo
	)
	{
		size_t bottomVer[3];
		size_t topVer[3];
		std::vector < size_t > prismCell(6);
		std::vector < size_t > pyramidCell(5);
		std::vector < size_t > tetraCell(4);
		for (size_t i = 0; i < cells.size(); ++i) {
			if (!mesh_utils::isInVector(boundary_layer_cell, i)) {
				if (vertMap.count(cells[i][0])) {
					prismTopo.vecCells[i][0] = vertMap.at(cells[i][0]).back();
				}
				if (vertMap.count(cells[i][1])) {
					prismTopo.vecCells[i][1] = vertMap.at(cells[i][1]).back();
				}
				if (vertMap.count(cells[i][2])) {
					prismTopo.vecCells[i][2] = vertMap.at(cells[i][2]).back();
				}
				continue;
			}
			//prism
			if (vertMap.at(cells[i][0]).size() > 1 && vertMap.at(cells[i][1]).size() > 1 && vertMap.at(cells[i][2]).size() > 1) {
				for (int j = 0; j < vertMap.at(cells[i][0]).size() - 1; ++j) {
					prismCell[0] = vertMap.at(cells[i][0])[j];
					prismCell[1] = vertMap.at(cells[i][2])[j];
					prismCell[2] = vertMap.at(cells[i][1])[j];
					prismCell[3] = vertMap.at(cells[i][0])[j + 1];
					prismCell[4] = vertMap.at(cells[i][2])[j + 1];
					prismCell[5] = vertMap.at(cells[i][1])[j + 1];
					prismTopo.vecCells.emplace_back(prismCell);
				}
			}
			//pyramid
			if (vertMap.at(cells[i][0]).size() == 1 && vertMap.at(cells[i][1]).size() > 1 && vertMap.at(cells[i][2]).size() > 1) {
				for (int j = 0; j < vertMap.at(cells[i][1]).size() - 1; ++j) {
					pyramidCell[0] = vertMap.at(cells[i][1])[j];
					pyramidCell[1] = vertMap.at(cells[i][2])[j];
					pyramidCell[2] = vertMap.at(cells[i][2])[j + 1];
					pyramidCell[3] = vertMap.at(cells[i][1])[j + 1];
					pyramidCell[4] = vertMap.at(cells[i][0])[0];
					prismTopo.vecCells.emplace_back(pyramidCell);
				}
			}
			if (vertMap.at(cells[i][0]).size() > 1 && vertMap.at(cells[i][1]).size() == 1 && vertMap.at(cells[i][2]).size() > 1) {
				for (int j = 0; j < vertMap.at(cells[i][0]).size() - 1; ++j) {
					pyramidCell[0] = vertMap.at(cells[i][2])[j];
					pyramidCell[1] = vertMap.at(cells[i][0])[j];
					pyramidCell[2] = vertMap.at(cells[i][0])[j + 1];
					pyramidCell[3] = vertMap.at(cells[i][2])[j + 1];
					pyramidCell[4] = vertMap.at(cells[i][1])[0];
					prismTopo.vecCells.emplace_back(pyramidCell);
				}
			}
			if (vertMap.at(cells[i][0]).size() > 1 && vertMap.at(cells[i][1]).size() > 1 && vertMap.at(cells[i][2]).size() == 1) {
				for (int j = 0; j < vertMap.at(cells[i][0]).size() - 1; ++j) {
					pyramidCell[0] = vertMap.at(cells[i][0])[j];
					pyramidCell[1] = vertMap.at(cells[i][1])[j];
					pyramidCell[2] = vertMap.at(cells[i][1])[j + 1];
					pyramidCell[3] = vertMap.at(cells[i][0])[j + 1];
					pyramidCell[4] = vertMap.at(cells[i][2])[0];
					prismTopo.vecCells.emplace_back(pyramidCell);
				}
			}
			//tetrahedron
			if (vertMap.at(cells[i][0]).size() > 1 && vertMap.at(cells[i][1]).size() == 1 && vertMap.at(cells[i][2]).size() == 1) {
				for (int j = 0; j < vertMap.at(cells[i][0]).size() - 1; ++j) {
					tetraCell[0] = vertMap.at(cells[i][0])[j];
					tetraCell[1] = vertMap.at(cells[i][2])[0];
					tetraCell[2] = vertMap.at(cells[i][1])[0];
					tetraCell[3] = vertMap.at(cells[i][0])[j + 1];
					prismTopo.vecCells.emplace_back(tetraCell);
				}
			}
			if (vertMap.at(cells[i][0]).size() == 1 && vertMap.at(cells[i][1]).size() > 1 && vertMap.at(cells[i][2]).size() == 1) {
				for (int j = 0; j < vertMap.at(cells[i][1]).size() - 1; ++j) {
					tetraCell[0] = vertMap.at(cells[i][1])[j];
					tetraCell[1] = vertMap.at(cells[i][0])[0];
					tetraCell[2] = vertMap.at(cells[i][2])[0];
					tetraCell[3] = vertMap.at(cells[i][1])[j + 1];
					prismTopo.vecCells.emplace_back(tetraCell);
				}
			}
			if (vertMap.at(cells[i][0]).size() == 1 && vertMap.at(cells[i][1]).size() == 1 && vertMap.at(cells[i][2]).size() > 1) {
				for (int j = 0; j < vertMap.at(cells[i][2]).size() - 1; ++j) {
					tetraCell[0] = vertMap.at(cells[i][2])[j];
					tetraCell[1] = vertMap.at(cells[i][1])[0];
					tetraCell[2] = vertMap.at(cells[i][0])[0];
					tetraCell[3] = vertMap.at(cells[i][2])[j + 1];
					prismTopo.vecCells.emplace_back(tetraCell);
				}
			}
		}
		return 1;
	}

	int computeMarchVertex(const global_type::Mesh& mesh,
		const global_type::MeshNormal& meshNormal,
		global_type::Parameter& param,
		global_type::Mesh& prismTopo,
		std::unordered_map < size_t, std::vector < size_t > >& vertMap
	)
	{
		double initHeight = param.initHeight;
		double increaseRatio = param.increaseRatio;
		int layerNum = param.initLayerNumber;
		std::vector < double > eps(layerNum);
		double sum = 0;
		for (int i = 0; i < layerNum; ++i) {
			sum += std::pow(increaseRatio, i);
		}
		eps[0] = initHeight / sum;
		double addEps = eps[0];
		for (int i = 1; i < layerNum; ++i) {
			eps[i] = eps[0] * std::pow(increaseRatio, i);
			eps[i] += addEps;
			addEps = eps[i];
		}
		param.eps = eps;

		prismTopo.vecVertices = mesh.vecVertices;
		prismTopo.matVertices = mesh.matVertices;
		prismTopo.vecCells = mesh.vecCells;
		prismTopo.matCells = mesh.matCells;
		generateVertices(mesh.vecVertices, meshNormal.boundary_vertex, meshNormal.verticesNormalizedNormal, layerNum, eps, prismTopo, vertMap);
		generateCell(mesh.vecCells, vertMap, param.boundary_layer_cell, prismTopo);

		//Test
		//mesh_io::saveVTK("data/new_opt_0_1_2_prism.vtk", prismTopo.vecVertices, prismTopo.vecCells);
		return 1;
	}
}
