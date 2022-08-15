#pragma once
#include <vector>
#include <string>
#include <Eigen/Core>
#ifndef TETLIBRARY
#define TETLIBRARY 
#endif
#include "tetgen.h"
namespace tetgen
{
    int tetrahedralization(const std::vector< Eigen::Vector3d >& hole_points,
        const std::vector< std::vector < double > >& tri_v, const std::vector< std::vector < size_t > >& tri_c,
        Eigen::MatrixXd& TV, Eigen::MatrixXi& TT, std::vector< size_t >& intersectF);
}

