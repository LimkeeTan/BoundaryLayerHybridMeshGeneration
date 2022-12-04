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
    std::string param = "pq1.2Y";
    Eigen::MatrixXd TTV;
    Eigen::MatrixXi TTT;
    tetrahedralize(hole_points, V, F, param, TTV, TTT, intersectF);
    Eigen::MatrixXd tet_L;
    Eigen::MatrixXd tet_v_L;
    Eigen::MatrixXi tet_c_L;
    tet_v_L.resize(TTV.rows(), 3);
    tet_c_L.resize(TTT.rows(), 4);
    for (size_t j = 0; j < tet_v_L.rows(); ++j) {
        tet_v_L.row(j) = TTV.row(j);
    }
    for (size_t j = 0; j < tet_c_L.rows(); ++j) {
        Eigen::Vector4i cc;
        cc[0] = TTT(j, 0);
        cc[1] = TTT(j, 1);
        cc[2] = TTT(j, 2);
        cc[3] = TTT(j, 3);
        tet_c_L.row(j) = cc;
    }
    tet_L.resize(tet_c_L.rows(), 6);
    for (int i = 0; i < tet_c_L.rows(); ++i) {
        tet_L(i, 0) = (tet_v_L.row(tet_c_L(i, 3)) - tet_v_L.row(tet_c_L(i, 0))).squaredNorm();
        tet_L(i, 1) = (tet_v_L.row(tet_c_L(i, 3)) - tet_v_L.row(tet_c_L(i, 1))).squaredNorm();
        tet_L(i, 2) = (tet_v_L.row(tet_c_L(i, 3)) - tet_v_L.row(tet_c_L(i, 2))).squaredNorm();
        tet_L(i, 3) = (tet_v_L.row(tet_c_L(i, 1)) - tet_v_L.row(tet_c_L(i, 2))).squaredNorm();
        tet_L(i, 4) = (tet_v_L.row(tet_c_L(i, 2)) - tet_v_L.row(tet_c_L(i, 0))).squaredNorm();
        tet_L(i, 5) = (tet_v_L.row(tet_c_L(i, 0)) - tet_v_L.row(tet_c_L(i, 1))).squaredNorm();
    }
    tet_L = tet_L.array().sqrt().eval();
    Eigen::VectorXd tet_vol;
    tet_vol.resize(tet_L.rows(), 1);
    for (int t = 0; t < tet_L.rows(); ++t) {
        const double u = tet_L(t, 0);
        const double v = tet_L(t, 1);
        const double w = tet_L(t, 2);
        const double U = tet_L(t, 3);
        const double V = tet_L(t, 4);
        const double W = tet_L(t, 5);
        const double X = (w - U + v) * (U + v + w);
        const double x = (U - v + w) * (v - w + U);
        const double Y = (u - V + w) * (V + w + u);
        const double y = (V - w + u) * (w - u + V);
        const double Z = (v - W + u) * (W + u + v);
        const double z = (W - u + v) * (u - v + W);
        const double a = sqrt(x * Y * Z);
        const double b = sqrt(y * Z * X);
        const double c = sqrt(z * X * Y);
        const double d = sqrt(x * y * z);
        tet_vol(t) = sqrt(
            (-a + b + c + d) *
            (a - b + c + d) *
            (a + b - c + d) *
            (a + b + c - d)) /
            (192. * u * v * w);
    }
    double aveOfTetVol = 0;
    for (int j = 0; j < tet_vol.size(); ++j) {
        aveOfTetVol += tet_vol[j];
    }
    aveOfTetVol /= tet_vol.size();
    param = "pq1.2Ya" + std::to_string(aveOfTetVol);
    tetrahedralize(hole_points, V, F, param, TV, TT, intersectF);
    return 1;
}
