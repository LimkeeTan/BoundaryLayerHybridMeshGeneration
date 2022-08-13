#pragma once
#include <vector>

#include "../generate_mesh/tetgen.h"

namespace tetgen
{
    bool tetgenio_to_tetmesh(
        const tetgenio& out,
        std::vector<std::vector<REAL > >& V,
        std::vector<std::vector<int> >& T);
}
