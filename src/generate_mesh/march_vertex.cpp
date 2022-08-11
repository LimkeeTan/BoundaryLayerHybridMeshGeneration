#include "march_vertex.h"

namespace march_vertex {
	int generateCell(const size_t& B0,
		const size_t& T0,
		const size_t& B1,
		const size_t& T1,
		const size_t& B2,
		const size_t& T2,
		std::vector < size_t >& singleCell
	)
	{
		//prism
		if (T0 != B0 && T1 != B1 && T2 != B2) {
			singleCell.emplace_back(B0);
			singleCell.emplace_back(B1);
			singleCell.emplace_back(B2);
			singleCell.emplace_back(T0);
			singleCell.emplace_back(T1);
			singleCell.emplace_back(T2);
		}

		//pyramid
		if (T0 == B0 && T1 != B1 && T2 != B2) {
			singleCell.emplace_back(B1);
			singleCell.emplace_back(T1);
			singleCell.emplace_back(T2);
			singleCell.emplace_back(B2);
			singleCell.emplace_back(B0);
		}
		if (T1 == B1 && T0 != B0 && T2 != B2) {
			singleCell.emplace_back(B2);
			singleCell.emplace_back(T2);
			singleCell.emplace_back(T0);
			singleCell.emplace_back(B0);
			singleCell.emplace_back(B1);
		}
		if (T2 == B2 && T0 != B0 && T1 != B1) {
			singleCell.emplace_back(B0);
			singleCell.emplace_back(T0);
			singleCell.emplace_back(T1);
			singleCell.emplace_back(B1);
			singleCell.emplace_back(B2);
		}

		//tetrahedral
		if (T0 != B0 && T1 == B1 && T2 == B2) {
			singleCell.emplace_back(B0);
			singleCell.emplace_back(B1);
			singleCell.emplace_back(B2);
			singleCell.emplace_back(T0);
		}
		if (T1 != B1 && T0 == B0 && T2 == B2) {
			singleCell.emplace_back(B0);
			singleCell.emplace_back(B1);
			singleCell.emplace_back(B2);
			singleCell.emplace_back(T1);
		}
		if (T2 != B2 && T0 == B0 && T1 == B1) {
			singleCell.emplace_back(B0);
			singleCell.emplace_back(B1);
			singleCell.emplace_back(B2);
			singleCell.emplace_back(T2);
		}

		//triangle
		if (T0 == B0 && T1 == B1 && T2 == B2) {
			singleCell.emplace_back(B0);
			singleCell.emplace_back(B1);
			singleCell.emplace_back(B2);
		}
		return 1;
	}

	int computeMarchVertex(const global_type::Mesh& mesh,
		const global_type::MeshNormal& meshNormal,
		global_type::Mesh& prismTopo
	)
	{
		double eps = 1;
		prismTopo.vecVertices = mesh.vecVertices;
		prismTopo.matVertices = mesh.matVertices;
		std::unordered_map < size_t, size_t > vertMap;
		std::vector < size_t > singleTri;
		size_t B0, T0, B1, T1, B2, T2;
		std::vector < double > b0(3), t0(3), b1(3), t1(3), b2(3), t2(3);
		for (size_t i = 0; i < mesh.vecCells.size(); ++i) {
			singleTri = mesh.vecCells[i];
			std::vector < size_t > singleCell;
			if (!vertMap.count(singleTri[0])) {
				b0 = mesh.vecVertices[singleTri[0]];
				if (meshNormal.verticesNormal[singleTri[0]][0] == 0 &&
					meshNormal.verticesNormal[singleTri[0]][1] == 0 &&
					meshNormal.verticesNormal[singleTri[0]][2] == 0)
				{
					B0 = singleTri[0];
					T0 = B0;
					vertMap.insert(std::pair < size_t, size_t >(singleTri[0], singleTri[0]));
				}
				else 
				{
					t0[0] = b0[0] + eps * meshNormal.verticesNormalizedNormal[singleTri[0]][0];
					t0[1] = b0[1] + eps * meshNormal.verticesNormalizedNormal[singleTri[0]][1];
					t0[2] = b0[2] + eps * meshNormal.verticesNormalizedNormal[singleTri[0]][2];
					prismTopo.vecVertices.emplace_back(t0);
					B0 = singleTri[0];
					T0 = prismTopo.vecVertices.size() - 1;
					vertMap.insert(std::pair < size_t, size_t >(singleTri[0], prismTopo.vecVertices.size() - 1));
				}
			}
			else {
				B0 = singleTri[0];
				T0 = vertMap[singleTri[0]];
			}

			if (!vertMap.count(singleTri[1])) {
				b1 = mesh.vecVertices[singleTri[1]];
				if (meshNormal.verticesNormal[singleTri[1]][0] == 0 &&
					meshNormal.verticesNormal[singleTri[1]][1] == 0 &&
					meshNormal.verticesNormal[singleTri[1]][2] == 0)
				{
					B1 = singleTri[1];
					T1 = B1;
					vertMap.insert(std::pair < size_t, size_t >(singleTri[1], singleTri[1]));
				}
				else
				{
					t1[0] = b1[0] + eps * meshNormal.verticesNormalizedNormal[singleTri[1]][0];
					t1[1] = b1[1] + eps * meshNormal.verticesNormalizedNormal[singleTri[1]][1];
					t1[2] = b1[2] + eps * meshNormal.verticesNormalizedNormal[singleTri[1]][2];
					prismTopo.vecVertices.emplace_back(t1);
					B1 = singleTri[1];
					T1 = prismTopo.vecVertices.size() - 1;
					vertMap.insert(std::pair < size_t, size_t >(singleTri[1], prismTopo.vecVertices.size() - 1));
				}
			}
			else {
				B1 = singleTri[1];
				T1 = vertMap[singleTri[1]];
			}

			if (!vertMap.count(singleTri[2])) {
				b2 = mesh.vecVertices[singleTri[2]];
				if (meshNormal.verticesNormal[singleTri[2]][0] == 0 &&
					meshNormal.verticesNormal[singleTri[2]][1] == 0 &&
					meshNormal.verticesNormal[singleTri[2]][2] == 0)
				{
					B2 = singleTri[2];
					T2 = B2;
					vertMap.insert(std::pair < size_t, size_t >(singleTri[2], singleTri[2]));
				}
				else 
				{
					t2[0] = b2[0] + eps * meshNormal.verticesNormalizedNormal[singleTri[2]][0];
					t2[1] = b2[1] + eps * meshNormal.verticesNormalizedNormal[singleTri[2]][1];
					t2[2] = b2[2] + eps * meshNormal.verticesNormalizedNormal[singleTri[2]][2];
					prismTopo.vecVertices.emplace_back(t2);
					B2 = singleTri[2];
					T2 = prismTopo.vecVertices.size() - 1;
					vertMap.insert(std::pair < size_t, size_t >(singleTri[2], prismTopo.vecVertices.size() - 1));
				}
			}
			else {
				B2 = singleTri[2];
				T2 = vertMap[singleTri[2]];
			}
			generateCell(B0, T0, B1, T1, B2, T2, singleCell);
			prismTopo.vecCells.emplace_back(singleCell);
		}

		//Test
		mesh_io::saveVTK("data/wanxiangjie_prism.vtk", prismTopo.vecVertices, prismTopo.vecCells);
		return 1;
	}
}