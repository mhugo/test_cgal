#include <iostream>

#include <vector>
#include <string>
#include <iterator>
#include <fstream>

#include <CGAL/Cartesian.h>
#include <CGAL/Origin.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Spherical_kernel_3.h>
#include <CGAL/Algebraic_kernel_for_spheres_2_3.h>

using namespace std;
using namespace CGAL;

typedef CGAL::Cartesian<double> CK;
typedef CGAL::Spherical_kernel_3<CK, CGAL::Algebraic_kernel_for_spheres_2_3<double> > K;

int main( int argc, char *argv[] )
{
    ifstream data_file( "../data/m0.off" );
    Polyhedron_3<K> poly;
    data_file >> poly;

    // Convert facet to list of Triangle_3
    Polyhedron_3<K>::Facet_iterator it;
    for ( it = poly.facets_begin(); it != poly.facets_end(); it++ )
    {
	assert( it->size() == 3 );
	Point_3<K> tri[3];
	Polyhedron_3<K>::Halfedge_around_facet_circulator eit;
	size_t s = 0;
	for ( eit = it->facet_begin(); s < circulator_size( eit ); s++ )
	{
	    Point_3<K> p = eit->vertex()->point();
	    tri[s] = p;
	}

	Triangle_3<K> triangle( tri[0], tri[1], tri[2] );
    }
    return 0;
}
