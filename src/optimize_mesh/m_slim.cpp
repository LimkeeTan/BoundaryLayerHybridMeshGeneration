#include "m_slim.h"

namespace slim_opt {
	int grad_tet_ref(const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& T,
		const std::vector < Eigen::MatrixXd >& RF,
		Eigen::SparseMatrix < double >& G
	)
	{
		using namespace Eigen;
		int n = V.rows();
		int m = T.rows();
		MatrixXi F(4 * m, 3);
		for (int i = 0; i < m; i++) {
			F.row(0 * m + i) << T(i, 0), T(i, 1), T(i, 2);
			F.row(1 * m + i) << T(i, 0), T(i, 2), T(i, 3);
			F.row(2 * m + i) << T(i, 0), T(i, 3), T(i, 1);
			F.row(3 * m + i) << T(i, 1), T(i, 3), T(i, 2);
		}
		VectorXd vol; 
		igl::volume(V, T, vol);
		VectorXd A(F.rows());
		MatrixXd N(F.rows(), 3);
		const int tet_faces[4][3] = {
			{ 0, 1, 2 },
			{ 0, 2, 3 },
			{ 0, 3, 1 },
			{ 1, 3, 2 }
		};
		double coeff = sqrt(2);
		double coeff_2 = sqrt(3);
		for (int i = 0; i < m; i++) {
			if (RF[i].rows() == 4) {
				Vector3d e01, e02, e03;
				e01 = RF[i].row(0) - RF[i].row(1);
				e02 = RF[i].row(0) - RF[i].row(2);
				e03 = RF[i].row(0) - RF[i].row(3);
				double volume_i = std::abs(e01.cross(e02).dot(e03)) / 6;
				double scale_i = vol(i) / volume_i;
				scale_i = std::cbrt(scale_i);

				for (int k = 0; k < 4; k++) {
					Vector3d e0, e1;
					e0 = (RF[i].row(tet_faces[k][0]) - RF[i].row(tet_faces[k][1])) * scale_i;
					e1 = (RF[i].row(tet_faces[k][0]) - RF[i].row(tet_faces[k][2])) * scale_i;

					double area = e0.cross(e1).norm() / 2;
					A(k * m + i) = area;
					e0.normalize(); e1.normalize();
					Vector3d normal = e0.cross(e1);
					N.row(k * m + i) = normal;
					N.row(k * m + i) /= N.row(k * m + i).norm();
				}
			}
			else {
				double currVol = std::cbrt(3 * vol(i));

				N.row(i) << 0, 0, 1;
				double a_0 = coeff * currVol;
				A(i) = (pow(a_0, 2) * coeff_2) / 4;

				N.row(m + i) << 0.8165, -0.4714, -0.3333;
				double a_1 = coeff * currVol;
				A(m + i) = (pow(a_1, 2) * coeff_2) / 4;

				N.row(2 * m + i) << 0, 0.9428, -0.3333;
				double a_2 = coeff * currVol;
				A(2 * m + i) = (pow(a_2, 2) * coeff_2) / 4;

				N.row(3 * m + i) << -0.8165, -0.4714, -0.3333;
				double a_3 = coeff * currVol;
				A(3 * m + i) = (pow(a_3, 2) * coeff_2) / 4;
			}
		}
		std::vector<Triplet<double> > G_t;
		for (int i = 0; i < 4 * m; i++) {
			int T_j; // j indexes : repmat([T(:,4);T(:,2);T(:,3);T(:,1)],3,1)
			switch (i / m) {
			case 0:
				T_j = 3;
				break;
			case 1:
				T_j = 1;
				break;
			case 2:
				T_j = 2;
				break;
			case 3:
				T_j = 0;
				break;
			}
			int i_idx = i % m;
			int j_idx = T(i_idx, T_j);

			double val_before_n = A(i) / (3 * vol(i_idx));
			G_t.push_back(Triplet<double>(0 * m + i_idx, j_idx, val_before_n * N(i, 0)));
			G_t.push_back(Triplet<double>(1 * m + i_idx, j_idx, val_before_n * N(i, 1)));
			G_t.push_back(Triplet<double>(2 * m + i_idx, j_idx, val_before_n * N(i, 2)));
		}
		G.resize(3 * m, n);
		G.setFromTriplets(G_t.begin(), G_t.end());
		return 1;
	}

	int pre_calc(SLIMData& s, std::vector < Eigen::MatrixXd >& targetPrismTet)
	{
		if (!s.has_pre_calc) {
			s.v_n = s.v_num;
			s.f_n = s.f_num;
			s.dim = 3;
			Eigen::SparseMatrix<double> G;
			grad_tet_ref(s.V, s.F, targetPrismTet, G);
			s.Dx = G.block(0, 0, s.F.rows(), s.V.rows());
			s.Dy = G.block(s.F.rows(), 0, s.F.rows(), s.V.rows());
			s.Dz = G.block(2 * s.F.rows(), 0, s.F.rows(), s.V.rows());
			s.W_11.resize(s.f_n);
			s.W_12.resize(s.f_n);
			s.W_13.resize(s.f_n);
			s.W_21.resize(s.f_n);
			s.W_22.resize(s.f_n);
			s.W_23.resize(s.f_n);
			s.W_31.resize(s.f_n);
			s.W_32.resize(s.f_n);
			s.W_33.resize(s.f_n);
			s.Dx.makeCompressed();
			s.Dy.makeCompressed();
			s.Dz.makeCompressed();
			s.Ri.resize(s.f_n, s.dim * s.dim);
			s.Ji.resize(s.f_n, s.dim * s.dim);
			s.rhs.resize(s.dim * s.v_num);
			// flattened weight matrix
			s.WGL_M.resize(s.dim * s.dim * s.f_n);
			for (int i = 0; i < s.dim * s.dim; i++)
				for (int j = 0; j < s.f_n; j++)
					s.WGL_M(i * s.f_n + j) = s.M(j);

			s.first_solve = true;
			s.has_pre_calc = true;
		}
		return 1;
	}

	int compute_jacobians(SLIMData& s, const Eigen::MatrixXd& uv)
	{
		if (s.F.cols() == 3)
		{
			// Ji=[D1*u,D2*u,D1*v,D2*v];
			s.Ji.col(0) = s.Dx * uv.col(0);
			s.Ji.col(1) = s.Dy * uv.col(0);
			s.Ji.col(2) = s.Dx * uv.col(1);
			s.Ji.col(3) = s.Dy * uv.col(1);
		}
		else /*tet mesh*/ {
			// Ji=[D1*u,D2*u,D3*u, D1*v,D2*v, D3*v, D1*w,D2*w,D3*w];
			s.Ji.col(0) = s.Dx * uv.col(0);
			s.Ji.col(1) = s.Dy * uv.col(0);
			s.Ji.col(2) = s.Dz * uv.col(0);
			s.Ji.col(3) = s.Dx * uv.col(1);
			s.Ji.col(4) = s.Dy * uv.col(1);
			s.Ji.col(5) = s.Dz * uv.col(1);
			s.Ji.col(6) = s.Dx * uv.col(2);
			s.Ji.col(7) = s.Dy * uv.col(2);
			s.Ji.col(8) = s.Dz * uv.col(2);
		}
		return 1;
	}

	double compute_energy_with_jacobians(SLIMData& s,
		const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::MatrixXd& Ji,
		Eigen::MatrixXd& uv,
		Eigen::VectorXd& areas)
	{

		double energy = 0;
		if (s.dim == 2)
		{
			Eigen::Matrix<double, 2, 2> ji;
			for (int i = 0; i < s.f_n; i++)
			{
				ji(0, 0) = Ji(i, 0);
				ji(0, 1) = Ji(i, 1);
				ji(1, 0) = Ji(i, 2);
				ji(1, 1) = Ji(i, 3);

				typedef Eigen::Matrix<double, 2, 2> Mat2;
				typedef Eigen::Matrix<double, 2, 1> Vec2;
				Mat2 ri, ti, ui, vi;
				Vec2 sing;
				igl::polar_svd(ji, ri, ti, ui, sing, vi);
				double s1 = sing(0);
				double s2 = sing(1);

				switch (s.slim_energy)
				{
				case SLIM_ENERGY::ARAP:
				{
					energy += areas(i) * (pow(s1 - 1, 2) + pow(s2 - 1, 2));
					break;
				}
				case SLIM_ENERGY::SYMMETRIC_DIRICHLET:
				{
					energy += areas(i) * (pow(s1, 2) + pow(s1, -2) + pow(s2, 2) + pow(s2, -2));
					break;
				}
				case SLIM_ENERGY::EXP_SYMMETRIC_DIRICHLET:
				{
					energy += areas(i) * exp(s.exp_factor * (pow(s1, 2) + pow(s1, -2) + pow(s2, 2) + pow(s2, -2)));
					break;
				}
				case SLIM_ENERGY::LOG_ARAP:
				{
					energy += areas(i) * (pow(log(s1), 2) + pow(log(s2), 2));
					break;
				}
				case SLIM_ENERGY::CONFORMAL:
				{
					energy += areas(i) * ((pow(s1, 2) + pow(s2, 2)) / (2 * s1 * s2));
					break;
				}
				case SLIM_ENERGY::EXP_CONFORMAL:
				{
					energy += areas(i) * exp(s.exp_factor * ((pow(s1, 2) + pow(s2, 2)) / (2 * s1 * s2)));
					break;
				}

				}

			}
		}
		else
		{
			Eigen::Matrix<double, 3, 3> ji;
			for (int i = 0; i < s.f_n; i++)
			{
				ji(0, 0) = Ji(i, 0);
				ji(0, 1) = Ji(i, 1);
				ji(0, 2) = Ji(i, 2);
				ji(1, 0) = Ji(i, 3);
				ji(1, 1) = Ji(i, 4);
				ji(1, 2) = Ji(i, 5);
				ji(2, 0) = Ji(i, 6);
				ji(2, 1) = Ji(i, 7);
				ji(2, 2) = Ji(i, 8);

				typedef Eigen::Matrix<double, 3, 3> Mat3;
				typedef Eigen::Matrix<double, 3, 1> Vec3;
				Mat3 ri, ti, ui, vi;
				Vec3 sing;
				igl::polar_svd(ji, ri, ti, ui, sing, vi);
				double s1 = sing(0);
				double s2 = sing(1);
				double s3 = sing(2);

				switch (s.slim_energy)
				{
				case SLIM_ENERGY::ARAP:
				{
					energy += areas(i) * (pow(s1 - 1, 2) + pow(s2 - 1, 2) + pow(s3 - 1, 2));
					break;
				}
				case SLIM_ENERGY::SYMMETRIC_DIRICHLET:
				{
					energy += areas(i) * (pow(s1, 2) + pow(s1, -2) + pow(s2, 2) + pow(s2, -2) + pow(s3, 2) + pow(s3, -2));
					break;
				}
				case SLIM_ENERGY::EXP_SYMMETRIC_DIRICHLET:
				{
					energy += areas(i) * exp(s.exp_factor *
						(pow(s1, 2) + pow(s1, -2) + pow(s2, 2) + pow(s2, -2) + pow(s3, 2) + pow(s3, -2)));
					break;
				}
				case SLIM_ENERGY::LOG_ARAP:
				{
					energy += areas(i) * (pow(log(s1), 2) + pow(log(std::abs(s2)), 2) + pow(log(std::abs(s3)), 2));
					break;
				}
				case SLIM_ENERGY::CONFORMAL:
				{
					energy += areas(i) * ((pow(s1, 2) + pow(s2, 2) + pow(s3, 2)) / (3 * pow(s1 * s2 * s3, 2. / 3.)));
					break;
				}
				case SLIM_ENERGY::EXP_CONFORMAL:
				{
					energy += areas(i) * exp((pow(s1, 2) + pow(s2, 2) + pow(s3, 2)) / (3 * pow(s1 * s2 * s3, 2. / 3.)));
					break;
				}
				}
			}
		}

		return energy;
	}

	double compute_soft_const_energy(const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::MatrixXd& V_o,
		SLIMData& s
	)
	{
		double e = 0;
		for (int i = 0; i < s.b.rows(); i++)
		{
			e += s.soft_const_p * (s.bc.row(i) - V_o.row(s.b(i))).squaredNorm();
		}
		return e;
	}

	double compute_total_energy(SLIMData& s, Eigen::MatrixXd& V_new)
	{
		compute_jacobians(s, V_new);
		double e = compute_energy_with_jacobians(s, s.V, s.F, s.Ji, V_new, s.M);
		s.energy_quality = e;
		double esoft = compute_soft_const_energy(s.V, s.F, V_new, s);
		s.energy_soft = esoft;
		e += esoft;
		return e;
	}

	int precompute(Eigen::MatrixXd& tetVer,
		Eigen::MatrixXi& tetCell,
		Eigen::MatrixXd& initTetVer,
		SLIMData& data,
		std::vector < Eigen::MatrixXd >& targetPrismTet
	)
	{
		data.V = tetVer;
		data.F = tetCell;
		data.V_o = initTetVer;

		data.v_num = tetVer.rows();
		data.f_num = tetCell.rows();
		data.slim_energy = CONFORMAL;
		data.mesh_improvement_3d = true;
		data.M.resize(tetCell.rows());
		data.weight_opt = 10;
		data.M.setConstant(data.weight_opt);
		data.mesh_area = data.M.sum();
		data.exp_factor = 1.0;

		data.soft_const_p = 1e4;

		data.proximal_p = 0.0001;

		pre_calc(data, targetPrismTet);
		data.energy = compute_total_energy(data, data.V_o) / data.mesh_area;
		std::cout << "total energy: " << data.energy << std::endl;
		std::cout << "ED: " << data.energy_quality / data.mesh_area << std::endl;
		std::cout << "ES: " << data.energy_soft / data.mesh_area << std::endl;
		return 1;
	}

	void update_weights_and_closest_rotations(SLIMData& s,
		const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		Eigen::MatrixXd& uv)
	{
		compute_jacobians(s, uv);

		const double eps = 1e-8;
		double exp_f = s.exp_factor;

		if (s.dim == 2)
		{
			for (int i = 0; i < s.Ji.rows(); ++i)
			{
				typedef Eigen::Matrix<double, 2, 2> Mat2;
				typedef Eigen::Matrix<double, 2, 1> Vec2;
				Mat2 ji, ri, ti, ui, vi;
				Vec2 sing;
				Vec2 closest_sing_vec;
				Mat2 mat_W;
				Vec2 m_sing_new;
				double s1, s2;

				ji(0, 0) = s.Ji(i, 0);
				ji(0, 1) = s.Ji(i, 1);
				ji(1, 0) = s.Ji(i, 2);
				ji(1, 1) = s.Ji(i, 3);

				igl::polar_svd(ji, ri, ti, ui, sing, vi);

				s1 = sing(0);
				s2 = sing(1);

				// Update Weights according to energy
				switch (s.slim_energy)
				{
				case SLIM_ENERGY::ARAP:
				{
					m_sing_new << 1, 1;
					break;
				}
				case SLIM_ENERGY::SYMMETRIC_DIRICHLET:
				{
					double s1_g = 2 * (s1 - pow(s1, -3));
					double s2_g = 2 * (s2 - pow(s2, -3));
					m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1)));
					break;
				}
				case SLIM_ENERGY::LOG_ARAP:
				{
					double s1_g = 2 * (log(s1) / s1);
					double s2_g = 2 * (log(s2) / s2);
					m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1)));
					break;
				}
				case SLIM_ENERGY::CONFORMAL:
				{
					double s1_g = 1 / (2 * s2) - s2 / (2 * pow(s1, 2));
					double s2_g = 1 / (2 * s1) - s1 / (2 * pow(s2, 2));

					double geo_avg = sqrt(s1 * s2);
					double s1_min = geo_avg;
					double s2_min = geo_avg;

					m_sing_new << sqrt(s1_g / (2 * (s1 - s1_min))), sqrt(s2_g / (2 * (s2 - s2_min)));

					// change local step
					closest_sing_vec << s1_min, s2_min;
					ri = ui * closest_sing_vec.asDiagonal() * vi.transpose();
					break;
				}
				case SLIM_ENERGY::EXP_CONFORMAL:
				{
					double s1_g = 2 * (s1 - pow(s1, -3));
					double s2_g = 2 * (s2 - pow(s2, -3));

					double geo_avg = sqrt(s1 * s2);
					double s1_min = geo_avg;
					double s2_min = geo_avg;

					double in_exp = exp_f * ((pow(s1, 2) + pow(s2, 2)) / (2 * s1 * s2));
					double exp_thing = exp(in_exp);

					s1_g *= exp_thing * exp_f;
					s2_g *= exp_thing * exp_f;

					m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1)));
					break;
				}
				case SLIM_ENERGY::EXP_SYMMETRIC_DIRICHLET:
				{
					double s1_g = 2 * (s1 - pow(s1, -3));
					double s2_g = 2 * (s2 - pow(s2, -3));

					double in_exp = exp_f * (pow(s1, 2) + pow(s1, -2) + pow(s2, 2) + pow(s2, -2));
					double exp_thing = exp(in_exp);

					s1_g *= exp_thing * exp_f;
					s2_g *= exp_thing * exp_f;

					m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1)));
					break;
				}
				}

				if (std::abs(s1 - 1) < eps) m_sing_new(0) = 1;
				if (std::abs(s2 - 1) < eps) m_sing_new(1) = 1;
				mat_W = ui * m_sing_new.asDiagonal() * ui.transpose();

				s.W_11(i) = mat_W(0, 0);
				s.W_12(i) = mat_W(0, 1);
				s.W_21(i) = mat_W(1, 0);
				s.W_22(i) = mat_W(1, 1);

				// 2) Update local step (doesn't have to be a rotation, for instance in case of conformal energy)
				s.Ri(i, 0) = ri(0, 0);
				s.Ri(i, 1) = ri(1, 0);
				s.Ri(i, 2) = ri(0, 1);
				s.Ri(i, 3) = ri(1, 1);
			}
		}
		else
		{
			typedef Eigen::Matrix<double, 3, 1> Vec3;
			typedef Eigen::Matrix<double, 3, 3> Mat3;
			Mat3 ji;
			Vec3 m_sing_new;
			Vec3 closest_sing_vec;
			const double sqrt_2 = sqrt(2);
			for (int i = 0; i < s.Ji.rows(); ++i)
			{
				ji(0, 0) = s.Ji(i, 0);
				ji(0, 1) = s.Ji(i, 1);
				ji(0, 2) = s.Ji(i, 2);
				ji(1, 0) = s.Ji(i, 3);
				ji(1, 1) = s.Ji(i, 4);
				ji(1, 2) = s.Ji(i, 5);
				ji(2, 0) = s.Ji(i, 6);
				ji(2, 1) = s.Ji(i, 7);
				ji(2, 2) = s.Ji(i, 8);

				Mat3 ri, ti, ui, vi;
				Vec3 sing;
				igl::polar_svd(ji, ri, ti, ui, sing, vi);

				double s1 = sing(0);
				double s2 = sing(1);
				double s3 = sing(2);

				// 1) Update Weights
				switch (s.slim_energy)
				{
				case SLIM_ENERGY::ARAP:
				{
					m_sing_new << 1, 1, 1;
					break;
				}
				case SLIM_ENERGY::LOG_ARAP:
				{
					double s1_g = 2 * (log(s1) / s1);
					double s2_g = 2 * (log(s2) / s2);
					double s3_g = 2 * (log(s3) / s3);
					m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1))), sqrt(s3_g / (2 * (s3 - 1)));
					break;
				}
				case SLIM_ENERGY::SYMMETRIC_DIRICHLET:
				{
					double s1_g = 2 * (s1 - pow(s1, -3));
					double s2_g = 2 * (s2 - pow(s2, -3));
					double s3_g = 2 * (s3 - pow(s3, -3));
					m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1))), sqrt(s3_g / (2 * (s3 - 1)));
					break;
				}
				case SLIM_ENERGY::EXP_SYMMETRIC_DIRICHLET:
				{
					double s1_g = 2 * (s1 - pow(s1, -3));
					double s2_g = 2 * (s2 - pow(s2, -3));
					double s3_g = 2 * (s3 - pow(s3, -3));
					m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1))), sqrt(s3_g / (2 * (s3 - 1)));

					double in_exp = exp_f * (pow(s1, 2) + pow(s1, -2) + pow(s2, 2) + pow(s2, -2) + pow(s3, 2) + pow(s3, -2));
					double exp_thing = exp(in_exp);

					s1_g *= exp_thing * exp_f;
					s2_g *= exp_thing * exp_f;
					s3_g *= exp_thing * exp_f;

					m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1))), sqrt(s3_g / (2 * (s3 - 1)));

					break;
				}
				case SLIM_ENERGY::CONFORMAL:
				{
					double common_div = 9 * (pow(s1 * s2 * s3, 5. / 3.));

					double s1_g = (-2 * s2 * s3 * (pow(s2, 2) + pow(s3, 2) - 2 * pow(s1, 2))) / common_div;
					double s2_g = (-2 * s1 * s3 * (pow(s1, 2) + pow(s3, 2) - 2 * pow(s2, 2))) / common_div;
					double s3_g = (-2 * s1 * s2 * (pow(s1, 2) + pow(s2, 2) - 2 * pow(s3, 2))) / common_div;

					double closest_s = sqrt(pow(s1, 2) + pow(s3, 2)) / sqrt_2;
					double s1_min = closest_s;
					double s2_min = closest_s;
					double s3_min = closest_s;

					m_sing_new << sqrt(s1_g / (2 * (s1 - s1_min))), sqrt(s2_g / (2 * (s2 - s2_min))), sqrt(
						s3_g / (2 * (s3 - s3_min)));

					// change local step
					closest_sing_vec << s1_min, s2_min, s3_min;
					ri = ui * closest_sing_vec.asDiagonal() * vi.transpose();
					break;
				}
				case SLIM_ENERGY::EXP_CONFORMAL:
				{
					// E_conf = (s1^2 + s2^2 + s3^2)/(3*(s1*s2*s3)^(2/3) )
					// dE_conf/ds1 = (-2*(s2*s3)*(s2^2+s3^2 -2*s1^2) ) / (9*(s1*s2*s3)^(5/3))
					// Argmin E_conf(s1): s1 = sqrt(s1^2+s2^2)/sqrt(2)
					double common_div = 9 * (pow(s1 * s2 * s3, 5. / 3.));

					double s1_g = (-2 * s2 * s3 * (pow(s2, 2) + pow(s3, 2) - 2 * pow(s1, 2))) / common_div;
					double s2_g = (-2 * s1 * s3 * (pow(s1, 2) + pow(s3, 2) - 2 * pow(s2, 2))) / common_div;
					double s3_g = (-2 * s1 * s2 * (pow(s1, 2) + pow(s2, 2) - 2 * pow(s3, 2))) / common_div;

					double in_exp = exp_f * ((pow(s1, 2) + pow(s2, 2) + pow(s3, 2)) / (3 * pow((s1 * s2 * s3), 2. / 3)));;
					double exp_thing = exp(in_exp);

					double closest_s = sqrt(pow(s1, 2) + pow(s3, 2)) / sqrt_2;
					double s1_min = closest_s;
					double s2_min = closest_s;
					double s3_min = closest_s;

					s1_g *= exp_thing * exp_f;
					s2_g *= exp_thing * exp_f;
					s3_g *= exp_thing * exp_f;

					m_sing_new << sqrt(s1_g / (2 * (s1 - s1_min))), sqrt(s2_g / (2 * (s2 - s2_min))), sqrt(
						s3_g / (2 * (s3 - s3_min)));

					// change local step
					closest_sing_vec << s1_min, s2_min, s3_min;
					ri = ui * closest_sing_vec.asDiagonal() * vi.transpose();
				}
				}
				if (std::abs(s1 - 1) < eps) m_sing_new(0) = 1;
				if (std::abs(s2 - 1) < eps) m_sing_new(1) = 1;
				if (std::abs(s3 - 1) < eps) m_sing_new(2) = 1;
				Mat3 mat_W;
				mat_W = ui * m_sing_new.asDiagonal() * ui.transpose();

				s.W_11(i) = mat_W(0, 0);
				s.W_12(i) = mat_W(0, 1);
				s.W_13(i) = mat_W(0, 2);
				s.W_21(i) = mat_W(1, 0);
				s.W_22(i) = mat_W(1, 1);
				s.W_23(i) = mat_W(1, 2);
				s.W_31(i) = mat_W(2, 0);
				s.W_32(i) = mat_W(2, 1);
				s.W_33(i) = mat_W(2, 2);

				// 2) Update closest rotations (not rotations in case of conformal energy)
				s.Ri(i, 0) = ri(0, 0);
				s.Ri(i, 1) = ri(1, 0);
				s.Ri(i, 2) = ri(2, 0);
				s.Ri(i, 3) = ri(0, 1);
				s.Ri(i, 4) = ri(1, 1);
				s.Ri(i, 5) = ri(2, 1);
				s.Ri(i, 6) = ri(0, 2);
				s.Ri(i, 7) = ri(1, 2);
				s.Ri(i, 8) = ri(2, 2);
			} // for loop end
		} // if dim end
	}

	void buildA(SLIMData& s, Eigen::SparseMatrix<double>& A)
	{
		// formula (35) in paper
		std::vector<Eigen::Triplet<double> > IJV;
		if (s.dim == 2)
		{
			IJV.reserve(4 * (s.Dx.outerSize() + s.Dy.outerSize()));
			for (int k = 0; k < s.Dx.outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(s.Dx, k); it; ++it)
				{
					int dx_r = it.row();
					int dx_c = it.col();
					double val = it.value();

					IJV.push_back(Eigen::Triplet<double>(dx_r, dx_c, val * s.W_11(dx_r)));
					IJV.push_back(Eigen::Triplet<double>(dx_r, s.v_n + dx_c, val * s.W_12(dx_r)));

					IJV.push_back(Eigen::Triplet<double>(2 * s.f_n + dx_r, dx_c, val * s.W_21(dx_r)));
					IJV.push_back(Eigen::Triplet<double>(2 * s.f_n + dx_r, s.v_n + dx_c, val * s.W_22(dx_r)));
				}
			}

			for (int k = 0; k < s.Dy.outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(s.Dy, k); it; ++it)
				{
					int dy_r = it.row();
					int dy_c = it.col();
					double val = it.value();

					IJV.push_back(Eigen::Triplet<double>(s.f_n + dy_r, dy_c, val * s.W_11(dy_r)));
					IJV.push_back(Eigen::Triplet<double>(s.f_n + dy_r, s.v_n + dy_c, val * s.W_12(dy_r)));

					IJV.push_back(Eigen::Triplet<double>(3 * s.f_n + dy_r, dy_c, val * s.W_21(dy_r)));
					IJV.push_back(Eigen::Triplet<double>(3 * s.f_n + dy_r, s.v_n + dy_c, val * s.W_22(dy_r)));
				}
			}
		}
		else
		{
			IJV.reserve(9 * (s.Dx.outerSize() + s.Dy.outerSize() + s.Dz.outerSize()));
			for (int k = 0; k < s.Dx.outerSize(); k++)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(s.Dx, k); it; ++it)
				{
					int dx_r = it.row();
					int dx_c = it.col();
					double val = it.value();

					IJV.push_back(Eigen::Triplet<double>(dx_r, dx_c, val * s.W_11(dx_r)));
					IJV.push_back(Eigen::Triplet<double>(dx_r, s.v_n + dx_c, val * s.W_12(dx_r)));
					IJV.push_back(Eigen::Triplet<double>(dx_r, 2 * s.v_n + dx_c, val * s.W_13(dx_r)));

					IJV.push_back(Eigen::Triplet<double>(3 * s.f_n + dx_r, dx_c, val * s.W_21(dx_r)));
					IJV.push_back(Eigen::Triplet<double>(3 * s.f_n + dx_r, s.v_n + dx_c, val * s.W_22(dx_r)));
					IJV.push_back(Eigen::Triplet<double>(3 * s.f_n + dx_r, 2 * s.v_n + dx_c, val * s.W_23(dx_r)));

					IJV.push_back(Eigen::Triplet<double>(6 * s.f_n + dx_r, dx_c, val * s.W_31(dx_r)));
					IJV.push_back(Eigen::Triplet<double>(6 * s.f_n + dx_r, s.v_n + dx_c, val * s.W_32(dx_r)));
					IJV.push_back(Eigen::Triplet<double>(6 * s.f_n + dx_r, 2 * s.v_n + dx_c, val * s.W_33(dx_r)));
				}
			}

			for (int k = 0; k < s.Dy.outerSize(); k++)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(s.Dy, k); it; ++it)
				{
					int dy_r = it.row();
					int dy_c = it.col();
					double val = it.value();

					IJV.push_back(Eigen::Triplet<double>(s.f_n + dy_r, dy_c, val * s.W_11(dy_r)));
					IJV.push_back(Eigen::Triplet<double>(s.f_n + dy_r, s.v_n + dy_c, val * s.W_12(dy_r)));
					IJV.push_back(Eigen::Triplet<double>(s.f_n + dy_r, 2 * s.v_n + dy_c, val * s.W_13(dy_r)));

					IJV.push_back(Eigen::Triplet<double>(4 * s.f_n + dy_r, dy_c, val * s.W_21(dy_r)));
					IJV.push_back(Eigen::Triplet<double>(4 * s.f_n + dy_r, s.v_n + dy_c, val * s.W_22(dy_r)));
					IJV.push_back(Eigen::Triplet<double>(4 * s.f_n + dy_r, 2 * s.v_n + dy_c, val * s.W_23(dy_r)));

					IJV.push_back(Eigen::Triplet<double>(7 * s.f_n + dy_r, dy_c, val * s.W_31(dy_r)));
					IJV.push_back(Eigen::Triplet<double>(7 * s.f_n + dy_r, s.v_n + dy_c, val * s.W_32(dy_r)));
					IJV.push_back(Eigen::Triplet<double>(7 * s.f_n + dy_r, 2 * s.v_n + dy_c, val * s.W_33(dy_r)));
				}
			}

			for (int k = 0; k < s.Dz.outerSize(); k++)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(s.Dz, k); it; ++it)
				{
					int dz_r = it.row();
					int dz_c = it.col();
					double val = it.value();

					IJV.push_back(Eigen::Triplet<double>(2 * s.f_n + dz_r, dz_c, val * s.W_11(dz_r)));
					IJV.push_back(Eigen::Triplet<double>(2 * s.f_n + dz_r, s.v_n + dz_c, val * s.W_12(dz_r)));
					IJV.push_back(Eigen::Triplet<double>(2 * s.f_n + dz_r, 2 * s.v_n + dz_c, val * s.W_13(dz_r)));

					IJV.push_back(Eigen::Triplet<double>(5 * s.f_n + dz_r, dz_c, val * s.W_21(dz_r)));
					IJV.push_back(Eigen::Triplet<double>(5 * s.f_n + dz_r, s.v_n + dz_c, val * s.W_22(dz_r)));
					IJV.push_back(Eigen::Triplet<double>(5 * s.f_n + dz_r, 2 * s.v_n + dz_c, val * s.W_23(dz_r)));

					IJV.push_back(Eigen::Triplet<double>(8 * s.f_n + dz_r, dz_c, val * s.W_31(dz_r)));
					IJV.push_back(Eigen::Triplet<double>(8 * s.f_n + dz_r, s.v_n + dz_c, val * s.W_32(dz_r)));
					IJV.push_back(Eigen::Triplet<double>(8 * s.f_n + dz_r, 2 * s.v_n + dz_c, val * s.W_33(dz_r)));
				}
			}
		}
		A.setFromTriplets(IJV.begin(), IJV.end());
	}

	void buildRhs(SLIMData& s, const Eigen::SparseMatrix<double>& At)
	{
		Eigen::VectorXd f_rhs(s.dim * s.dim * s.f_n);
		f_rhs.setZero();
		if (s.dim == 2)
		{
			for (int i = 0; i < s.f_n; i++)
			{
				f_rhs(i + 0 * s.f_n) = s.W_11(i) * s.Ri(i, 0) + s.W_12(i) * s.Ri(i, 1);
				f_rhs(i + 1 * s.f_n) = s.W_11(i) * s.Ri(i, 2) + s.W_12(i) * s.Ri(i, 3);
				f_rhs(i + 2 * s.f_n) = s.W_21(i) * s.Ri(i, 0) + s.W_22(i) * s.Ri(i, 1);
				f_rhs(i + 3 * s.f_n) = s.W_21(i) * s.Ri(i, 2) + s.W_22(i) * s.Ri(i, 3);
			}
		}
		else
		{
			for (int i = 0; i < s.f_n; i++)
			{
				f_rhs(i + 0 * s.f_n) = s.W_11(i) * s.Ri(i, 0) + s.W_12(i) * s.Ri(i, 1) + s.W_13(i) * s.Ri(i, 2);
				f_rhs(i + 1 * s.f_n) = s.W_11(i) * s.Ri(i, 3) + s.W_12(i) * s.Ri(i, 4) + s.W_13(i) * s.Ri(i, 5);
				f_rhs(i + 2 * s.f_n) = s.W_11(i) * s.Ri(i, 6) + s.W_12(i) * s.Ri(i, 7) + s.W_13(i) * s.Ri(i, 8);
				f_rhs(i + 3 * s.f_n) = s.W_21(i) * s.Ri(i, 0) + s.W_22(i) * s.Ri(i, 1) + s.W_23(i) * s.Ri(i, 2);
				f_rhs(i + 4 * s.f_n) = s.W_21(i) * s.Ri(i, 3) + s.W_22(i) * s.Ri(i, 4) + s.W_23(i) * s.Ri(i, 5);
				f_rhs(i + 5 * s.f_n) = s.W_21(i) * s.Ri(i, 6) + s.W_22(i) * s.Ri(i, 7) + s.W_23(i) * s.Ri(i, 8);
				f_rhs(i + 6 * s.f_n) = s.W_31(i) * s.Ri(i, 0) + s.W_32(i) * s.Ri(i, 1) + s.W_33(i) * s.Ri(i, 2);
				f_rhs(i + 7 * s.f_n) = s.W_31(i) * s.Ri(i, 3) + s.W_32(i) * s.Ri(i, 4) + s.W_33(i) * s.Ri(i, 5);
				f_rhs(i + 8 * s.f_n) = s.W_31(i) * s.Ri(i, 6) + s.W_32(i) * s.Ri(i, 7) + s.W_33(i) * s.Ri(i, 8);
			}
		}
		Eigen::VectorXd uv_flat(s.dim * s.v_n);
		for (int i = 0; i < s.dim; i++)
			for (int j = 0; j < s.v_n; j++)
				uv_flat(s.v_n * i + j) = s.V_o(j, i);

		s.rhs = (At * s.WGL_M.asDiagonal() * f_rhs + s.proximal_p * uv_flat);
	}

	int add_soft_constraints(SLIMData& s, Eigen::SparseMatrix<double>& L)
	{
		int v_n = s.v_num;
		for (int d = 0; d < s.dim; d++)
		{
			for (int i = 0; i < s.b.rows(); i++)
			{
				int v_idx = s.b(i);
				s.rhs(d * v_n + v_idx) += s.soft_const_p * s.bc(i, d); // rhs
				L.coeffRef(d * v_n + v_idx, d * v_n + v_idx) += s.soft_const_p; // diagonal of matrix
			}
		}
		return 1;
	}

	void build_linear_system(SLIMData& s, Eigen::SparseMatrix<double>& L)
	{
		// formula (35) in paper
		Eigen::SparseMatrix<double> A(s.dim * s.dim * s.f_n, s.dim * s.v_n);
		buildA(s, A);

		Eigen::SparseMatrix<double> At = A.transpose();
		At.makeCompressed();

		Eigen::SparseMatrix<double> id_m(At.rows(), At.rows());
		id_m.setIdentity();

		// add proximal penalty
		L = At * s.WGL_M.asDiagonal() * A + s.proximal_p * id_m; //add also a proximal term
		L.makeCompressed();

		buildRhs(s, At);
		Eigen::SparseMatrix<double> OldL = L;

		//cout << L << endl;
		//cout << s.rhs << endl;
		//if (stitching)
		//{
		//	add_soft_constraints(s, L);
		//}
		add_soft_constraints(s, L);
		L.makeCompressed();
	}

	void solve_weighted_arap(SLIMData& s,
		const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		Eigen::MatrixXd& uv,
		Eigen::VectorXi& soft_b_p,
		Eigen::MatrixXd& soft_bc_p)
	{
		using namespace Eigen;

		Eigen::SparseMatrix<double> L;
		build_linear_system(s, L);

		//bool cholmod_definition = false;

		// solve
		Eigen::VectorXd Uc;

		if (s.dim == 2)
		{
			SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
			Uc = solver.compute(L).solve(s.rhs);
		}
		else
		{ // seems like CG performs much worse for 2D and way better for 3D
			uint32_t Num_variables = uv.rows() * s.dim;

			Eigen::VectorXd guess(uv.rows() * s.dim);
			for (int i = 0; i < s.v_n; i++) for (int j = 0; j < s.dim; j++) guess(uv.rows() * j + i) = uv(i, j); // flatten vector
			ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Upper> solver;
			solver.setTolerance(1e-8);
			solver.setMaxIterations(10);
			Uc = solver.compute(L).solveWithGuess(s.rhs, guess);
		}

		for (int i = 0; i < s.dim; i++)
		{
			uv.col(i) = Uc.block(i * s.v_n, 0, s.v_n, 1);
		}
	}

	Eigen::MatrixXd slimSolve(SLIMData& data, int iter_num)
	{
		for (int i = 0; i < iter_num; ++i)
		{
			Eigen::MatrixXd dest_res;
			dest_res = data.V_o;

			// Solve Weighted Proxy
			update_weights_and_closest_rotations(data, data.V, data.F, dest_res);
			solve_weighted_arap(data, data.V, data.F, dest_res, data.b, data.bc);

			double old_energy = data.energy;

			std::function<double(Eigen::MatrixXd&)> compute_energy = [&](
				Eigen::MatrixXd& aaa) { return compute_total_energy(data, aaa); };

			data.energy = igl::flip_avoiding_line_search(data.F, data.V_o, dest_res, compute_energy,
				data.energy * data.mesh_area) / data.mesh_area;
			std::cout << "total energy: " << data.energy << std::endl;
			std::cout << "ED: " << data.energy_quality / data.mesh_area << std::endl;
			std::cout << "ES: " << data.energy_soft / data.mesh_area << std::endl;
			std::cout << std::endl;
		}
		return data.V_o;
	}

	int constructSurfaceMesh(const std::vector < std::vector < double > >& vertices,
		const std::vector < std::vector < size_t > >& faces,
		pmp::SurfaceMesh& surfaceMesh
	) 
	{
		for (size_t i = 0; i < vertices.size(); ++i) {
			pmp::Point p(vertices[i][0], vertices[i][1], vertices[i][2]);
			surfaceMesh.add_vertex(p);
		}
		std::vector < pmp::Vertex > vert(3);
		for (size_t i = 0; i < faces.size(); ++i) {
			pmp::Vertex v0(faces[i][0]);
			pmp::Vertex v1(faces[i][1]);
			pmp::Vertex v2(faces[i][2]);
			vert[0] = v0;
			vert[1] = v1;
			vert[2] = v2;
			surfaceMesh.add_face(vert);
		}
		return 1;
	}

	int getBoundaryVerConstraints(const global_type::Mesh& tetMesh,
		SLIMData& data,
		std::vector < Eigen::VectorXi >& bs,
		std::vector < Eigen::MatrixXd >& bcs,
		std::unordered_map < size_t, Eigen::Vector3d >& v_boundary_map
	)
	{
		for (size_t i = 0; i < tetMesh.boundaryVerNums; ++i) {
			v_boundary_map.insert(std::make_pair(i, tetMesh.matVertices.row(i)));
		}
		return 1;
	}

	int getLaplaceConstraints(const global_type::Mesh& hybridMesh,
		const global_type::Mesh& tetMesh,
		const std::unordered_map < size_t, Eigen::Vector3d >& v_boundary_map,
		SLIMData& data,
		std::vector < Eigen::VectorXi >& bs,
		std::vector < Eigen::MatrixXd >& bcs,
		std::unordered_map < size_t, Eigen::Vector3d >& v_laplace_map
	)
	{
		std::vector < std::vector < double > > copyVertices = hybridMesh.vecVertices;
		std::vector < std::vector < size_t > > topTriangle;
		std::vector < size_t > singleTriangle(3);
		for (size_t i = 0; i < hybridMesh.vecCells.size(); ++i) {
			if (hybridMesh.vecCells[i].size() == 6) {
				singleTriangle[0] = hybridMesh.vecCells[i][3];
				singleTriangle[1] = hybridMesh.vecCells[i][4];
				singleTriangle[2] = hybridMesh.vecCells[i][5];
				topTriangle.emplace_back(singleTriangle);
			}
		}
		pmp::SurfaceMesh mesh;
		constructSurfaceMesh(copyVertices, topTriangle, mesh);
		pmp::Smoothing smooth(mesh);
		smooth.explicit_smoothing(1, false);
		pmp::write(mesh, "data/smooth.obj");
		Eigen::Vector3d secondVert;
		for (size_t i = 0; i < topTriangle.size(); ++i) {
			for (int j = 0; j < 3; ++j) {
				if (!v_boundary_map.count(topTriangle[i][j]) && !v_laplace_map.count(topTriangle[i][j])) {
					secondVert(0) = mesh.position(pmp::Vertex(topTriangle[i][j]))[0];
					secondVert(1) = mesh.position(pmp::Vertex(topTriangle[i][j]))[1];
					secondVert(2) = mesh.position(pmp::Vertex(topTriangle[i][j]))[2];
					v_laplace_map.insert(std::make_pair(topTriangle[i][j], secondVert));
				}
			}
		}
		return 1;
	}

	int getHeightConstraints(const global_type::Parameter& param,
		const global_type::Mesh& hybridMesh,
		const global_type::Mesh& tetMesh,
		const std::unordered_map < size_t, Eigen::Vector3d >& v_boundary_map,
		SLIMData& data,
		std::vector < Eigen::VectorXi >& bs,
		std::vector < Eigen::MatrixXd >& bcs,
		std::unordered_map < size_t, Eigen::Vector3d >& v_height_map
	)
	{
		Eigen::VectorXi b;
		Eigen::MatrixXd bc;
		Eigen::Vector3d topPoint;
		Eigen::Vector3d bottomPoint;
		Eigen::Vector3d normal;
		for (size_t i = 0; i < hybridMesh.boundaryCellsNums; ++i) {
			//prism
			if (hybridMesh.vecCells[i].size() == 6) {
				for (int j = 0; j < 3; ++j) {
					if (!v_boundary_map.count(hybridMesh.vecCells[i][j + 3]) && !v_height_map.count(hybridMesh.vecCells[i][j + 3])) {
						for (int k = 0; k < 3; ++k) {
							topPoint(k) = hybridMesh.vecVertices[hybridMesh.vecCells[i][j + 3]][k];
							bottomPoint(k) = hybridMesh.vecVertices[hybridMesh.vecCells[i][j]][k];
						}
						normal = (topPoint - bottomPoint).normalized();
						normal *= param.idealHeight[hybridMesh.vecCells[i][j]];
						normal += bottomPoint;
						v_height_map.insert(std::make_pair(hybridMesh.vecCells[i][j + 3], normal));
					}
				}
			}
			//pyramid
			if (hybridMesh.vecCells[i].size() == 5) {
				if (!v_boundary_map.count(hybridMesh.vecCells[i][2]) && !v_height_map.count(hybridMesh.vecCells[i][2])) {
					for (int j = 0; j < 3; ++j) {
						topPoint(j) = hybridMesh.vecVertices[hybridMesh.vecCells[i][2]][j];
						bottomPoint(j) = hybridMesh.vecVertices[hybridMesh.vecCells[i][1]][j];
					}
					normal = (topPoint - bottomPoint).normalized();
					normal *= param.idealHeight[hybridMesh.vecCells[i][1]];
					normal += bottomPoint;
					v_height_map.insert(std::make_pair(hybridMesh.vecCells[i][2], normal));
				}
				if (!v_boundary_map.count(hybridMesh.vecCells[i][3]) && !v_height_map.count(hybridMesh.vecCells[i][3])) {
					for (int j = 0; j < 3; ++j) {
						topPoint(j) = hybridMesh.vecVertices[hybridMesh.vecCells[i][3]][j];
						bottomPoint(j) = hybridMesh.vecVertices[hybridMesh.vecCells[i][0]][j];
					}
					normal = (topPoint - bottomPoint).normalized();
					normal *= param.idealHeight[hybridMesh.vecCells[i][0]];
					normal += bottomPoint;
					v_height_map.insert(std::make_pair(hybridMesh.vecCells[i][3], normal));
				}
			}
			//tet
			if (hybridMesh.vecCells[i].size() == 4) {
				if (!v_boundary_map.count(hybridMesh.vecCells[i][3]) && !v_height_map.count(hybridMesh.vecCells[i][3])) {
					for (int j = 0; j < 3; ++j) {
						topPoint(j) = hybridMesh.vecVertices[hybridMesh.vecCells[i][3]][j];
						bottomPoint(j) = hybridMesh.vecVertices[hybridMesh.vecCells[i][0]][j];
					}
					normal = (topPoint - bottomPoint).normalized();
					normal *= param.idealHeight[hybridMesh.vecCells[i][0]];
					normal += bottomPoint;
					v_height_map.insert(std::make_pair(hybridMesh.vecCells[i][3], normal));
				}
			}
		}
		return 1;
	}

	int getSoftConstraints(const global_type::Parameter& param,
		const global_type::Mesh& hybridMesh,
		const global_type::Mesh& tetMesh,
		SLIMData& data
	)
	{
		std::vector < Eigen::VectorXi > bs;
		std::vector < Eigen::MatrixXd > bcs;
		std::unordered_map < size_t, Eigen::Vector3d > v_boundary_map;
		std::unordered_map < size_t, Eigen::Vector3d > v_height_map;
		std::unordered_map < size_t, Eigen::Vector3d > v_laplace_map;
		getBoundaryVerConstraints(tetMesh, data, bs, bcs, v_boundary_map);
		getHeightConstraints(param, hybridMesh, tetMesh, v_boundary_map, data, bs, bcs, v_height_map);
		getLaplaceConstraints(hybridMesh, tetMesh, v_boundary_map, data, bs, bcs, v_laplace_map);
		data.b.resize(v_boundary_map.size() + v_height_map.size() + v_laplace_map.size(), 1);
		data.bc.resize(v_boundary_map.size() + v_height_map.size() + v_laplace_map.size(), 3);
		size_t count = 0;
		std::unordered_map < size_t, Eigen::Vector3d >::iterator it = v_boundary_map.begin();
		while (it != v_boundary_map.end()) {
			data.b(count) = it->first;
			data.bc(count, 0) = it->second(0);
			data.bc(count, 1) = it->second(1);
			data.bc(count, 2) = it->second(2);
			++count;
			++it;
		}
		it = v_height_map.begin();
		while (it != v_height_map.end()) {
			data.b(count) = it->first;
			data.bc(count, 0) = it->second(0);
			data.bc(count, 1) = it->second(1);
			data.bc(count, 2) = it->second(2);
			++count;
			++it;
		}
		it = v_laplace_map.begin();
		while (it != v_laplace_map.end()) {
			data.b(count) = it->first;
			data.bc(count, 0) = it->second(0);
			data.bc(count, 1) = it->second(1);
			data.bc(count, 2) = it->second(2);
			++count;
			++it;
		}
		return 1;
	}

	int slimOptimization(const global_type::Parameter& param,
		const global_type::Mesh& hybridMesh,
		global_type::Mesh& tetMesh,
		Eigen::MatrixXd& initTetVer,
		std::vector < Eigen::MatrixXd >& targetPrismTet
	)
	{
		SLIMData data;
		int iter = 5;
		getSoftConstraints(param, hybridMesh, tetMesh, data);
		std::cout << "Precompute..." << std::endl;
		Eigen::MatrixXd tetVer = tetMesh.matVertices;
		Eigen::MatrixXi tetCell = tetMesh.matCells;
		precompute(tetVer, tetCell, initTetVer, data, targetPrismTet);
		std::cout << "Solve..." << std::endl;
		data.V_o = slimSolve(data, iter);
		initTetVer = data.V_o;
		tetMesh.matVertices = data.V_o;
		return 1;
	}
}