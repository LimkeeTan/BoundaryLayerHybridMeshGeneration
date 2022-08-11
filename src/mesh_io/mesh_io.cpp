#include "mesh_io.h"

namespace mesh_io {
	int readVTK(const std::string& filename,
		std::vector < std::vector < double > >& vertices,
		std::vector < std::vector < size_t > >& cells
	)
	{
		vertices.clear();
		cells.clear();
		std::ifstream is(filename);
		if (!is.is_open()) {
			std::cout << "Unable to open file";
			return 0;
		}
		std::string line;
		char str[50];
		size_t verticesNum;
		size_t cellsNum;
		size_t cellsLists;
		int cellType;
		double cellVertices[8];
		std::vector < size_t > cell;
		while (std::getline(is, line)) {
			std::istringstream iss(line);
			if (line.find("POINTS") != std::string::npos) {
				iss >> str >> verticesNum >> str;
				vertices.resize(verticesNum);
				break;
			}
		}
		double x, y, z;
		for (size_t i = 0; i < verticesNum; ++i) {
			std::getline(is, line);
			std::istringstream iss(line);
			iss >> x >> y >> z;
			vertices[i].emplace_back(x);
			vertices[i].emplace_back(y);
			vertices[i].emplace_back(z);
		}
		while (std::getline(is, line)) {
			std::istringstream iss(line);
			if (line.find("CELLS") != std::string::npos) {
				iss >> str >> cellsNum >> cellsLists;
				cells.resize(cellsNum);
				break;
			}
		}
		for (size_t i = 0; i < cellsNum; ++i) {
			std::getline(is, line);
			std::istringstream iss(line);
			iss >> cellType;
			switch (cellType) {
			case 1:
				iss >> cellVertices[0];
				cells[i].emplace_back(cellVertices[0]);
				break;
			case 2:
				iss >> cellVertices[0] >> cellVertices[1];
				for (int j = 0; j < 2; ++j) {
					cells[i].emplace_back(cellVertices[j]);
				}
				break;
			case 3:
				iss >> cellVertices[0] >> cellVertices[1]
					>> cellVertices[2];
				for (int j = 0; j < 3; ++j) {
					cells[i].emplace_back(cellVertices[j]);
				}
				break;
			case 4:
				iss >> cellVertices[0] >> cellVertices[1]
					>> cellVertices[2] >> cellVertices[3];
				for (int j = 0; j < 4; ++j) {
					cells[i].emplace_back(cellVertices[j]);
				}
				break;
			case 5:
				iss >> cellVertices[0] >> cellVertices[1]
					>> cellVertices[2] >> cellVertices[3]
					>> cellVertices[4];
				for (int j = 0; j < 5; ++j) {
					cells[i].emplace_back(cellVertices[j]);
				}
				break;
			case 6:
				iss >> cellVertices[0] >> cellVertices[1]
					>> cellVertices[2] >> cellVertices[3]
					>> cellVertices[4] >> cellVertices[5];
				for (int j = 0; j < 6; ++j) {
					cells[i].emplace_back(cellVertices[j]);
				}
				break;
			case 7:
				iss >> cellVertices[0] >> cellVertices[1]
					>> cellVertices[2] >> cellVertices[3]
					>> cellVertices[4] >> cellVertices[5]
					>> cellVertices[6];
				for (int j = 0; j < 7; ++j) {
					cells[i].emplace_back(cellVertices[j]);
				}
				break;
			case 8:
				iss >> cellVertices[0] >> cellVertices[1]
					>> cellVertices[2] >> cellVertices[3]
					>> cellVertices[4] >> cellVertices[5]
					>> cellVertices[6] >> cellVertices[7];
				for (int j = 0; j < 8; ++j) {
					cells[i].emplace_back(cellVertices[j]);
				}
				break;
			default:
				std::cout << "Error : Unsupport Type\n";
				return 0;
				break;
			}
		}
		is.close();
		return 1;
	}

	int saveVTK(const std::string& filename,
		const std::vector < std::vector < double > >& vertices,
		const std::vector < std::vector < size_t > >& cells
	)
	{
		std::ofstream os;
		os.open(filename);
		os << "# vtk DataFile Version 3.0\nUnstructured mesh\nASCII\nDATASET UNSTRUCTURED_GRID\n";

		os << "POINTS " << vertices.size() << " double\n";
		for (size_t i = 0; i < vertices.size(); ++i) {
			os << vertices[i][0] << " " << vertices[i][1] << " " << vertices[i][2] << "\n";
		}
		os << "CELLS " << cells.size() << " ";
		size_t cellListsNum = 0;
		for (size_t i = 0; i < cells.size(); ++i) {
			cellListsNum += (cells[i].size() + 1);
		}
		os << cellListsNum << "\n";
		std::vector < int > vtkTypes;
		for (size_t i = 0; i < cells.size(); ++i) {
			switch (cells[i].size()) {
			case 1:
				os << "1 " << cells[i][0] << "\n";
				vtkTypes.emplace_back(1);
				break;
			case 2:
				os << "2 " << cells[i][0] << " " << cells[i][1] << "\n";
				vtkTypes.emplace_back(3);
				break;
			case 3:
				os << "3 " << cells[i][0] << " " << cells[i][1] << " "
					<< cells[i][2] << "\n";
				vtkTypes.emplace_back(5);
				break;
			case 4:
				os << "4 " << cells[i][0] << " " << cells[i][1] << " "
					<< cells[i][2] << " " << cells[i][3] << "\n";
				vtkTypes.emplace_back(10);
				break;
			case 5:
				os << "5 " << cells[i][0] << " " << cells[i][1] << " "
					<< cells[i][2] << " " << cells[i][3] << " "
					<< cells[i][4] << "\n";
				vtkTypes.emplace_back(14);
				break;
			case 6:
				os << "6 " << cells[i][0] << " " << cells[i][1] << " "
					<< cells[i][2] << " " << cells[i][3] << " "
					<< cells[i][4] << " " << cells[i][5] << "\n";
				vtkTypes.emplace_back(13);
				break;
			case 8:
				os << "8 " << cells[i][0] << " " << cells[i][1] << " "
					<< cells[i][2] << " " << cells[i][3] << " "
					<< cells[i][4] << " " << cells[i][5] << " "
					<< cells[i][6] << " " << cells[i][7] << "\n";
				vtkTypes.emplace_back(12);
				break;
			default:
				std::cout << "Error : Unsupport Type\n";
				return 0;
				break;
			}
		}
		os << "CELL_TYPES " << cells.size() << "\n";
		for (size_t i = 0; i < cells.size(); ++i)
		{
			os << vtkTypes[i] << "\n";
		}
		os.close();
		return 1;
	}

	int saveVTK(const std::string& filename,
		const Eigen::MatrixXd& vertices,
		const Eigen::MatrixXi& cells,
		const std::vector < std::vector < double > >& vertexNormal
	)
	{
		std::ofstream os;
		os.open(filename);
		os << "# vtk DataFile Version 3.0\nVolume mesh\nASCII\n\nDATASET UNSTRUCTURED_GRID\n";

		os << "POINTS " << vertices.rows() + vertexNormal.size() << " double\n";
		for (size_t i = 0; i < vertices.rows(); ++i) {
			os << vertices(i, 0) << " " << vertices(i, 1) << " " << vertices(i, 2) << "\n";
		}
		std::vector < std::vector < double > > normalPoints;
		normalPoints.resize(vertexNormal.size());
		for (size_t i = 0; i < vertexNormal.size(); ++i)
		{
			std::vector < double > oriVertex(3);
			oriVertex[0] = vertices(i, 0);
			oriVertex[1] = vertices(i, 1);
			oriVertex[2] = vertices(i, 2);
			std::vector < double > modifyVertex(3);
			modifyVertex[0] = oriVertex[0] + vertexNormal[i][0];
			modifyVertex[1] = oriVertex[1] + vertexNormal[i][1];
			modifyVertex[2] = oriVertex[2] + vertexNormal[i][2];
			normalPoints[i] = modifyVertex;
		}
		for (size_t i = 0; i < vertexNormal.size(); ++i)
		{
			os << normalPoints[i][0] << " " << normalPoints[i][1] << " " << normalPoints[i][2] << "\n";
		}
		os << "CELLS " << cells.rows() + vertexNormal.size() << " " << cells.rows() * 4 + vertexNormal.size() * 3 << "\n";
		for (size_t i = 0; i < cells.rows(); ++i) {
			os << 3 << " " << cells(i, 0) << " " << cells(i, 1) << " " << cells(i, 2) << "\n";
		}
		for (size_t i = 0; i < vertexNormal.size(); ++i)
		{
			os << 2 << " " << i << " " << i + vertices.rows() << "\n";
		}
		os << "CELL_TYPES " << cells.rows() + vertexNormal.size() << "\n";
		for (size_t i = 0; i < cells.rows(); ++i) {
			os << 5 << "\n";
		}
		for (size_t i = 0; i < vertexNormal.size(); ++i)
		{
			os << 3 << "\n";
		}
		os << "POINT_DATA " << vertices.rows() + vertexNormal.size() << "\n";
		os << "SCALARS fixed int" << "\n";
		os << "LOOKUP_TABLE default" << "\n";
		for (size_t i = 0; i < vertices.rows(); ++i) {
			os << 0 << "\n";
		}
		for (size_t i = 0; i < vertexNormal.size(); ++i)
		{
			os << 1 << "\n";
		}
		os.close();
		return 1;
	}

	int readTriOBJ(const std::string& filename,
		std::vector < std::vector < double > >& vertices,
		std::vector < std::vector < size_t > >& cells
	)
	{
		vertices.clear();
		cells.clear();
		std::ifstream is(filename);
		if (!is) {
			std::cout << "ERROR : FILE IS NOT EXIST!" << "\n";
			abort();
		}
		std::vector<double>  vs;
		std::vector<size_t>  fs;
		std::string line, pair[3];
		double  node[3];
		int  tri;
		while (!is.eof()) {
			std::getline(is, line);
			if (line.empty() || 13 == line[0])
				continue;
			std::istringstream instream(line);
			std::string word;
			instream >> word;
			if ("v" == word || "V" == word) {
				instream >> node[0] >> node[1] >> node[2];
				for (size_t j = 0; j < 3; ++j) {
					vs.emplace_back(node[j]);
				}
			}
			else if ('f' == word[0] || 'F' == word[0]) {
				instream >> pair[0] >> pair[1] >> pair[2];
				for (size_t j = 0; j < 3; ++j) {
					tri = strtoul(pair[j].c_str(), NULL, 10) - 1;
					fs.emplace_back(tri);
				}
			}
		}
		for (size_t i = 0; i < vs.size(); i += 3) {
			std::vector < double > cv(3);
			cv[0] = vs[i];
			cv[1] = vs[i + 1];
			cv[2] = vs[i + 2];
			vertices.emplace_back(cv);
		}
		for (size_t i = 0; i < fs.size(); i += 3) {
			std::vector < size_t > cf(3);
			cf[0] = fs[i];
			cf[1] = fs[i + 1];
			cf[2] = fs[i + 2];
			cells.emplace_back(cf);
		}
		is.close();
		return 1;
	}

	int saveTriOBJ(const std::string& filename,
		const std::vector < std::vector < double > >& vertices,
		const std::vector < std::vector < size_t > >& cells
	)
	{
		std::ofstream os(filename);
		for (size_t i = 0; i < vertices.size(); ++i) {
			os << "v " << vertices[i][0] << " " << vertices[i][1] << " " << vertices[i][2] << "\n";
		}
		for (size_t i = 0; i < cells.size(); ++i) {
			os << "f " << cells[i][0] + 1 << " " << cells[i][1] + 1 << " " << cells[i][2] + 1 << "\n";
		}
		os.close();
		return 1;
	}

	int saveNormalFile(const std::string& filename,
		const std::vector < std::vector < double > >& VertexNormal
	)
	{
		std::ofstream os(filename);
		if (!os)
			return 0;
		for (size_t i = 0; i < VertexNormal.size(); ++i)
		{
			os << VertexNormal[i][0] << " " << VertexNormal[i][1] << " " << VertexNormal[i][2] << "\n";
		}
		os.close();
		return 1;
	}
}