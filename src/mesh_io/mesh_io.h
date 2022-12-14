#ifndef MESH_IO_H_
#define MESH_IO_H_
#include <Eigen/Eigen>
#include <iostream>
#include <fstream>

namespace mesh_io {
	int readVTK(const std::string& filename,
		std::vector < std::vector < double > >& vertices,
		std::vector < std::vector < size_t > >& cells
	);

	int readVTK(const std::string& filename,
		std::vector < std::vector < double > >& vertices,
		std::vector < std::vector < size_t > >& cells,
		std::vector < int >& boundary_layer_label,
		std::vector < size_t >& boundary_layer_cells
	);

	int saveVTK(const std::string& filename,
		const std::vector < std::vector < double > >& vertices,
		const std::vector < std::vector < size_t > >& cells
	);

	int saveVTK(const std::string& filename,
		const Eigen::MatrixXd& vertices,
		const Eigen::MatrixXi& cells,
		const std::vector < std::vector < double > >& vertexNormal
	);

	int readTriOBJ(const std::string& filename,
		std::vector < std::vector < double > >& vertices,
		std::vector < std::vector < size_t > >& cells
	);

	int saveTriOBJ(const std::string& filename,
		const std::vector < std::vector < double > >& vertices,
		const std::vector < std::vector < size_t > >& cells
	);

	int saveNormalFile(const std::string& filename,
		const std::vector < std::vector < double > >& VertexNormal
	);
}

#endif