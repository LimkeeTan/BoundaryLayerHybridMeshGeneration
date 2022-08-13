#pragma once
// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.

#include <vector>
#include <string>
#include <Eigen/Core>
#ifndef TETLIBRARY
#define TETLIBRARY 
#endif
#include "tetgen.h"
namespace tetgen
{
    //// Mesh the interior of a surface mesh (V,F) using tetgen
    ////
    //// Inputs:
    ////   V  #V by 3 vertex position list
    ////   F  #F list of polygon face indices into V (0-indexed)
    ////   switches  string of tetgen options (See tetgen documentation) e.g.
    ////     "pq1.414a0.01" tries to mesh the interior of a given surface with
    ////       quality and area constraints
    ////     "" will mesh the convex hull constrained to pass through V (ignores F)
    //// Outputs:
    ////   TV  #V by 3 vertex position list
    ////   TT  #T by 4 list of tet face indices
    ////   TF  #F by 3 list of triangle face indices
    //// Returns status:
    ////   0 success
    ////   1 tetgen threw exception
    ////   2 tetgen did not crash but could not create any tets (probably there are
    ////     holes, duplicate faces etc.)
    ////   -1 other error
    //int tetrahedralize(
    //    const std::vector<Eigen::Vector3d>& hole_points,
    //    const std::vector<std::vector<REAL > >& V,
    //    const std::vector<std::vector<int> >& F,
    //    const std::string switches,
    //    std::vector<std::vector<REAL > >& TV,
    //    std::vector<std::vector<int > >& TT,
    //    std::vector<size_t>& intersectF);

    //// Wrapper with Eigen types
    //// Templates:
    ////   DerivedV  real-value: i.e. from MatrixXd
    ////   DerivedF  integer-value: i.e. from MatrixXi
    //int tetrahedralize(
    //    const std::vector<Eigen::Vector3d>& hole_points,
    //    const Eigen::MatrixXd& V,
    //    const Eigen::MatrixXi& F,
    //    const std::string switches,
    //    Eigen::MatrixXd& TV,
    //    Eigen::MatrixXi& TT,
    //    std::vector<size_t>& intersectF);

    int tetrahedralization(const std::vector< Eigen::Vector3d >& hole_points,
        const std::vector< std::vector < double > >& tri_v, const std::vector< std::vector < size_t > >& tri_c,
        Eigen::MatrixXd& TV, Eigen::MatrixXi& TT, std::vector< size_t >& intersectF);
}

