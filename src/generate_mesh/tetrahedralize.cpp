#include "tetrahedralize.h"
#include <cassert>
#include <iostream>
#include "tetgen.h"
#include "../mesh_io/mesh_to_tetgenio.h"
#include "../mesh_io/tetgenio_to_tetmesh.h"


void matrix_to_list(
    const Eigen::MatrixXd& M,
    std::vector<std::vector< double > >& V)
{
    V.resize(M.rows(), std::vector< double >(M.cols()));
    for (int i = 0; i < M.rows(); i++)
    {
        for (int j = 0; j < M.cols(); j++)
        {
            V[i][j] = M(i, j);
        }
    }
}

void matrix_to_list(
    const Eigen::MatrixXi& M,
    std::vector<std::vector< int > >& V)
{
    V.resize(M.rows(), std::vector< int >(M.cols()));
    for (int i = 0; i < M.rows(); i++)
    {
        for (int j = 0; j < M.cols(); j++)
        {
            V[i][j] = M(i, j);
        }
    }
}

void list_to_matrix(const std::vector< std::vector< double > >& V,
    Eigen::MatrixXd& M)
{
    if (V.size() == 0) return;
    int m = V.size();
    int n = V[0].size();
    M.resize(m, n);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            M(i, j) = V[i][j];
        }
    }
}

void list_to_matrix(const std::vector< std::vector< int > >& V,
    Eigen::MatrixXi& M)
{
    if (V.size() == 0) return;
    int m = V.size();
    int n = V[0].size();
    M.resize(m, n);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            M(i, j) = V[i][j];
        }
    }
}

int tetrahedralize(
    const std::vector<Eigen::Vector3d>& hole_points,
    const std::vector<std::vector<REAL > >& V,
    const std::vector<std::vector<int> >& F,
    const std::string switches,
    std::vector<std::vector<REAL > >& TV,
    std::vector<std::vector<int > >& TT,
    std::vector<size_t>& intersectF)
{
    using namespace std;
    tetgenio in, out;
    bool success;
    success = tetgen::mesh_to_tetgenio(hole_points, V, F, in);
    if (!success)
    {
        return -1;
    }
    try
    {
        char* cswitches = new char[switches.size() + 1];
        std::strcpy(cswitches, switches.c_str());
        ::tetrahedralize(intersectF, cswitches, &in, &out);
        delete[] cswitches;
    }
    catch (int e)
    {
        cerr << "^" << __FUNCTION__ << ": TETGEN CRASHED... KABOOOM!!!" << endl;
        return 1;
    }
    if (out.numberoftetrahedra == 0)
    {
        cerr << "^" << __FUNCTION__ << ": Tetgen failed to create tets" << endl;
        return 2;
    }
    success = tetgen::tetgenio_to_tetmesh(out, TV, TT);
    if (!success)
    {
        return -1;
    }
    return 0;
}
int tetrahedralize(
    const std::vector<Eigen::Vector3d>& hole_points,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const std::string switches,
    Eigen::MatrixXd& TV,
    Eigen::MatrixXi& TT,
    std::vector<size_t>& intersectF)
{
    using namespace std;
    vector<vector<REAL> > vV, vTV;
    vector<vector<int> > vF, vTT;
    matrix_to_list(V, vV);
    matrix_to_list(F, vF);
    int e = tetrahedralize(hole_points, vV, vF, switches, vTV, vTT, intersectF);
    if (e == 0)
    {
        list_to_matrix(vTV, TV);
        list_to_matrix(vTT, TT);
    }
    return e;
}

int tetgen::tetrahedralization(const std::vector<Eigen::Vector3d>& hole_points,
    const std::vector< std::vector < double > >& tri_v, const std::vector< std::vector < size_t > >& tri_c,
    Eigen::MatrixXd& TV, Eigen::MatrixXi& TT, std::vector<size_t>& intersectF)
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXd B;
    V.resize(tri_v.size(), 3);
    F.resize(tri_c.size(), 3);
    for (size_t i = 0; i < V.rows(); ++i)
    {
        V(i, 0) = tri_v[i][0];
        V(i, 1) = tri_v[i][1];
        V(i, 2) = tri_v[i][2];
    }
    for (size_t i = 0; i < F.rows(); ++i)
    {
        F(i, 0) = tri_c[i][0];
        F(i, 1) = tri_c[i][1];
        F(i, 2) = tri_c[i][2];
    }
    std::string param = "pq1.1Y";
    tetrahedralize(hole_points, V, F, param, TV, TT, intersectF);
    return 1;
}