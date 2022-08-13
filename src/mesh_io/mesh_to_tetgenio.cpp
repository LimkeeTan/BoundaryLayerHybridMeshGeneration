#include "mesh_to_tetgenio.h"

#include <cassert>
#include <Eigen/Core>
#include "../generate_mesh/tetgen.h"


bool tetgen::mesh_to_tetgenio(
    const std::vector<Eigen::Vector3d>& hole_points,
    const std::vector<std::vector<REAL > >& V,
    const std::vector<std::vector<int> >& F,
    tetgenio& in)
{
    using namespace std;
    // all indices start from 0
    in.firstnumber = 0;

    in.numberofpoints = V.size();
    in.pointlist = new REAL[in.numberofpoints * 3];
    // loop over points
    for (int i = 0; i < (int)V.size(); i++)
    {
        assert(V[i].size() == 3);
        in.pointlist[i * 3 + 0] = V[i][0];
        in.pointlist[i * 3 + 1] = V[i][1];
        in.pointlist[i * 3 + 2] = V[i][2];
    }

    in.numberoffacets = F.size();
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];

    // loop over face
    for (int i = 0; i < (int)F.size(); i++)
    {
        in.facetmarkerlist[i] = i;
        tetgenio::facet* f = &in.facetlist[i];
        f->numberofpolygons = 1;
        f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
        f->numberofholes = 0;
        f->holelist = NULL;
        tetgenio::polygon* p = &f->polygonlist[0];
        p->numberofvertices = F[i].size();
        p->vertexlist = new int[p->numberofvertices];
        // loop around face
        for (int j = 0; j < (int)F[i].size(); j++)
        {
            p->vertexlist[j] = F[i][j];
        }
    }

    // loop over hole 
    //in.numberofholes; in.holelist;
    if (hole_points.size() > 0) {
        in.numberofholes = hole_points.size();
        in.holelist = new REAL[in.numberofholes * 3];
        for (int i = 0; i < 3 * in.numberofholes; i += 3) {
            //holelist[i] = ?
            in.holelist[i] = hole_points[i / 3][0];
            in.holelist[i + 1] = hole_points[i / 3][1];
            in.holelist[i + 2] = hole_points[i / 3][2];
        }
    }
    return true;
}
