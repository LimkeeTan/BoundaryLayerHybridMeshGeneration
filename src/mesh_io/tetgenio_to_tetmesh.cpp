#include "tetgenio_to_tetmesh.h"
#include <iostream>
#include <unordered_map>
#include <Eigen/Core>

bool tetgen::tetgenio_to_tetmesh(
    const tetgenio& out,
    std::vector<std::vector<REAL > >& V,
    std::vector<std::vector<int> >& T)
{
    using namespace std;
    // process points
    if (out.pointlist == NULL)
    {
        cerr << "^tetgenio_to_tetmesh Error: point list is NULL\n" << endl;
        return false;
    }
    V.resize(out.numberofpoints, vector<REAL>(3));
    // loop over points
    for (int i = 0; i < out.numberofpoints; i++)
    {
        V[i][0] = out.pointlist[i * 3 + 0];
        V[i][1] = out.pointlist[i * 3 + 1];
        V[i][2] = out.pointlist[i * 3 + 2];
    }


    // process tets
    if (out.tetrahedronlist == NULL)
    {
        cerr << "^tetgenio_to_tetmesh Error: tet list is NULL\n" << endl;
        return false;
    }

    // When would this not be 4?
    assert(out.numberofcorners == 4);
    T.resize(out.numberoftetrahedra, vector<int>(out.numberofcorners));
    int min_index = 1e7;
    int max_index = -1e7;
    // loop over tetrahedra
    for (int i = 0; i < out.numberoftetrahedra; i++)
    {
        for (int j = 0; j < out.numberofcorners; j++)
        {
            int index = out.tetrahedronlist[i * out.numberofcorners + j];
            T[i][j] = index;
            min_index = (min_index > index ? index : min_index);
            max_index = (max_index < index ? index : max_index);
        }
    }
    assert(min_index >= 0);
    assert(max_index >= 0);
    assert(max_index < (int)V.size());
    return true;
}
