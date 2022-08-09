#pragma once
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

	int saveVTK(const std::string& filename,
		const std::vector < std::vector < double > >& vertices,
		const std::vector < std::vector < size_t > >& cells
	);

	int readTriOBJ(const std::string& filename,
		std::vector < std::vector < double > >& vertices,
		std::vector < std::vector < size_t > >& cells
	);

	int saveTriOBJ(const std::string& filename,
		const std::vector < std::vector < double > >& vertices,
		const std::vector < std::vector < size_t > >& cells
	);
}

#endif