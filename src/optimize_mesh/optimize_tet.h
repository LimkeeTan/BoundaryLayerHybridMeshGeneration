#ifndef OPTIMIZE_TET_H_
#define OPTIMIZE_TET_H_
#include "../global_type.h"

namespace optimize_tet {
	int optimize_tet(const global_type::Parameter& param, global_type::Mesh& hybrid_mesh);
}

#endif