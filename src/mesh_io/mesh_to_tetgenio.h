#pragma once
#include <vector>
#include <Eigen/Dense>
#include "../generate_mesh/tetgen.h"

namespace tetgen
{
    bool mesh_to_tetgenio(
        const std::vector<Eigen::Vector3d>& hole_points,
        const std::vector<std::vector<REAL > >& V,
        const std::vector<std::vector<int> >& F,
        tetgenio& in);
}
