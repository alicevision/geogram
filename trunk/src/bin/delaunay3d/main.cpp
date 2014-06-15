#include <geogram/delaunay/delaunay.h>
#include <geogram/numerics/predicates.h>
#include <iostream>
#include <fstream>

namespace {
    using namespace GEO;

    /**
     * \brief Computes the 3d Delaunay triangulation of
     *  a set of points obtained from an input stream and
     *  saves the result into an output stream.
     * \param[in] in the input stream, with one point per
     *  line, as x,y,z coordinates separated with spaces.
     * \param[out] out the output stream, in TET file format.
     */
    void delaunay(std::istream& in, std::ostream& out) {
        vector<double> points;
        std::cerr << "Loading" << std::endl;
        while(in) {
            double x;
            in >> x;
            points.push_back(x);
        }

        std::cerr << "Delaunay" << std::endl;
        Delaunay3d del;
        del.set_vertices(index_t(points.size()/3), &points[0]);

        std::cerr << "Saving" << std::endl;
        out << del.nb_vertices() << " vertices" << std::endl;
        out << del.nb_cells() << " cells" << std::endl;
        for(index_t v=0; v<del.nb_vertices(); ++v) {
            out << del.vertex_ptr(v)[0] << " " 
                << del.vertex_ptr(v)[1] << " " 
                << del.vertex_ptr(v)[2] 
                << std::endl;
        }
        for(index_t t=0; t<del.nb_cells(); ++t) {
            out << "4 " 
                << del.cell_vertex(t,0) << " "
                << del.cell_vertex(t,1) << " "
                << del.cell_vertex(t,2) << " "
                << del.cell_vertex(t,3) 
                << std::endl;
        }
        std::cerr << "Predicates statistics:" << std::endl;
        PCK::show_stats();
    }
}


int main(int argc, char** argv) {
    switch(argc) {
    case 1: {
        delaunay(std::cin, std::cout);
    } break;
    case 2: {
        std::ifstream in(argv[1]);
        delaunay(in,std::cout);
    } break;
    case 3: {
        std::ifstream in(argv[1]);
        std::ofstream out(argv[2]);
        delaunay(in,out);
    } break;
    default: {
        std::cerr << argv[0] << ": invalid number of arguments"
                  << std::endl;
        return -1;
    }
    }
    return 0;
}
