#ifndef M_SLIM_H_
#define M_SLIM_H_
#include <pmp/algorithms/Smoothing.h>
#include <pmp/io/io.h>
#include "Eigen/Eigen"
#include "igl/grad.h"
#include "igl/slim.h"
#include "../global_type.h"

namespace slim_opt {
	enum SLIM_ENERGY
	{
		ARAP,
		LOG_ARAP,
		SYMMETRIC_DIRICHLET,
		CONFORMAL,
		EXP_CONFORMAL,
		EXP_SYMMETRIC_DIRICHLET
	};
	struct SLIMData
	{
		// Input
		Eigen::MatrixXd V; // #V by 3 list of mesh vertex positions
		Eigen::MatrixXi F; // #F by 3/3 list of mesh faces (triangles/tets)
		SLIM_ENERGY slim_energy;

		double weight_opt = 1;
		// Optional Input
		// soft constraints
		Eigen::VectorXi b;
		Eigen::MatrixXd bc;
		double soft_const_p = 1;

		double exp_factor; // used for exponential energies, ignored otherwise
		bool mesh_improvement_3d = true; // only supported for 3d

		double energy_quality = 0;
		double energy_soft = 0;
		double lambda_s = 1;
		// Output
		Eigen::MatrixXd V_o; // #V by dim list of mesh vertex positions (dim = 2 for parametrization, 3 otherwise)
		double energy; // objective value

		// INTERNAL
		Eigen::VectorXd M;
		double mesh_area = 1;
		double avg_edge_length;
		int v_num;
		int f_num;
		double proximal_p;

		Eigen::VectorXd WGL_M;
		Eigen::VectorXd rhs;
		Eigen::MatrixXd Ri, Ji;
		Eigen::VectorXd W_11; Eigen::VectorXd W_12; Eigen::VectorXd W_13;
		Eigen::VectorXd W_21; Eigen::VectorXd W_22; Eigen::VectorXd W_23;
		Eigen::VectorXd W_31; Eigen::VectorXd W_32; Eigen::VectorXd W_33;
		Eigen::MatrixXd W;
		Eigen::SparseMatrix<double> Dx, Dy, Dz;
		int f_n, v_n;
		bool first_solve;
		bool has_pre_calc = false;
		int dim;
	};

	int slimOptimization(const global_type::Parameter& param,
		const global_type::Mesh& hybridMesh,
		global_type::Mesh& tetMesh,
		Eigen::MatrixXd& initTetVer,
		std::vector < Eigen::MatrixXd >& targetPrismTet
	);

	void tet_optimization(const global_type::Mesh& hybrid_mesh,
		global_type::Mesh& tet_mesh,
		Eigen::MatrixXd& init_tet_ver,
		std::vector < Eigen::MatrixXd >& target_tet,
		std::unordered_map < size_t, Eigen::Vector3d >& v_boundary_map);
}

#endif