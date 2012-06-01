//
// g++ test.cc -o test -lCGAL -frounding-math `pkg-config QtGui --cflags` -lQtGui -lQtCore
#include <iostream>

#include <vector>
#include <string>
#include <iterator>
#include <fstream>

#include <QImage>

#include <CGAL/Cartesian.h>
#include <CGAL/Origin.h>
#include <CGAL/Sphere_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Spherical_kernel_3.h>
#include <CGAL/Exact_spherical_kernel_3.h>
#include <CGAL/Algebraic_kernel_for_spheres_2_3.h>

using namespace std;

using CGAL::ORIGIN;
using CGAL::Object;
using CGAL::Ray_3;
using CGAL::Point_2;
using CGAL::Point_3;
using CGAL::Vector_3;
using CGAL::Triangle_3;
using CGAL::Line_3;
using CGAL::Segment_3;
using CGAL::Sphere_3;
using CGAL::Polyhedron_3;
using CGAL::Plane_3;
using CGAL::to_double;

typedef CGAL::Cartesian<double> CK;
typedef CGAL::Spherical_kernel_3<CK, CGAL::Algebraic_kernel_for_spheres_2_3<double> > K;
//typedef CGAL::Exact_spherical_kernel_3 K;

struct Color
{
    unsigned char r,g,b;
};

struct SceneObject
{
    CGAL::Object object;
    Color color;
    bool is_reflective;
};

struct SunLight
{
    Point_3<K> location;
};

struct Scene
{
    typedef std::list<SceneObject> Objects;
    Objects objects;
    SunLight light;
};

typedef std::list< CGAL::Triangle_3<K> > TIN;

class ProjectionScreen
{
protected:
    double ambient_factor;
public:
    ProjectionScreen( const Point_2<K>& cam_ul, const Point_2<K>& cam_lr, double cam_z, size_t width, size_t height ) :
	cam_ul_(cam_ul), cam_lr_(cam_lr), cam_z_(cam_z), width_(width), height_(height), ambient_factor(0.0)
    {
    }

    void set_ambient( double f )
    {
	ambient_factor = f;
    }

    Point_3<K> to_scene( size_t x, size_t y )
    {
	return Point_3<K>( x / double(width_) * (cam_lr_.x() - cam_ul_.x()) + cam_ul_.x(),
			y / double(height_) * (cam_lr_.y() - cam_ul_.y()) + cam_ul_.y(),
			cam_z_ );
			
    }

    //
    // Returns the first intersection point between a ray and an object
    // Returns false if no intersection exists
    bool first_intersection( Ray_3<K> ray, CGAL::Object obj, Point_3<K>& ipoint, Vector_3<K>& normal )
    {
	Sphere_3<K> sp;
	Plane_3<K> pl;
	TIN tin;
	bool intersects = false;

	// intersection with a sphere
	if ( assign( sp, obj ) )
	{
	    std::vector<CGAL::Object> olist;
	    intersection( ray.supporting_line(), sp, std::back_inserter(olist) );
	    // Pass if no intersection
	    if ( olist.size() > 0 )
	    {
		typedef std::pair<CGAL::Circular_arc_point_3<K>, unsigned> OT;
		OT o;

		std::vector<CGAL::Object>::const_iterator pit;
		for ( pit = olist.begin(); pit != olist.end(); pit++ )
		{
		    if ( !assign( o, *pit ) )
			continue;

		    // Here we introduce some approximations
		    ipoint = Point_3<K>( to_double(o.first.x()), to_double(o.first.y()), to_double(o.first.z()) );

		    // do not consider the ray source as an intersection point
		    bool is_the_source = squared_distance( ipoint, ray.source() ) < 0.0001;
		    if ( is_the_source )
			continue;
		    
		    Vector_3<K>v = ray.source() - ipoint;
		    double d = ray.to_vector() * v;
		    // If the intersection point is outside the ray, pass
		    // <-> dot product of the ray vs the vector from the ray source to the point is not null
		    if ( d > 0.001 )
			continue;
		    
		    // Only consider the first intersection found (the closest from the camera)
		    break;
		}

		if ( pit != olist.end() )
		{
		    // we have an intersection point and it is not the ray origin
		    intersects = true;
		    normal = ipoint - sp.center();
		    normal = normal / sqrt(normal.squared_length());
		}
	    }
	}
	// intersection with a plane
	else if ( assign( pl, obj ) )
	{
	    CGAL::Object p = intersection( pl, ray );
	    if ( assign( ipoint, p ) )
	    {
		intersects = true;
		normal = pl.orthogonal_vector();
		normal = normal / sqrt(normal.squared_length());
	    }
	}
	// intersection with a mesh
	else if ( assign (tin, obj) )
	{
	    for ( TIN::const_iterator it = tin.begin(); it != tin.end(); it++ )
	    {
		Point_3<K> p;
		Segment_3<K> seg;

		CGAL::Object iobj = intersection( ray, *it );
		if ( assign( p, iobj ) )
		{
		    intersects = true;
		    normal = it->supporting_plane().orthogonal_vector();
		    normal = normal / sqrt( normal.squared_length() );
		}
		else if ( assign ( seg, iobj ) )
		{
		    // we have a colinear segment with our ray, return the closer point to the ray source
		    // TODO
		    intersects = false;
		}
	    }
	}
	return intersects;
    }

    void compute_color( const Scene& scene, const Ray_3<K>& ray, uchar& red, uchar& green, uchar& blue )
    {
	static int recurse_level = 0;
	if (recurse_level > 3)
	{
	    return;
	}
	recurse_level++;

	// maps a depth to an (r,g,b) tuple. Since it is a map, it is ordered by depth
	std::map<double, unsigned long> colors;
	for ( Scene::Objects::const_iterator it = scene.objects.begin(); it != scene.objects.end(); it++ )
	{
	    CGAL::Object obj = it->object;
	    Point_3<K> ipoint;
	    Vector_3<K> normal;
	    bool intersects = first_intersection( ray, obj, ipoint, normal );
	 
	    if ( !intersects )
		continue;

	    // trace a ray from the intersection point to the light source
	    const Point_3<K>& light = scene.light.location;
	    Ray_3<K> to_light( ipoint, light );

	    bool reach_light = true;
	    for ( Scene::Objects::const_iterator it2 = scene.objects.begin(); it2 != scene.objects.end(); it2++ )
	    {
		Point_3<K> sun_point;
		Vector_3<K> unused;
		intersects = first_intersection( to_light, it2->object, sun_point, unused );
		if ( intersects )
		{
		    reach_light = false;
		    break;
		}
	    }

	    double shade = 0.0;
	    // reflective color
	    unsigned char rr, rg, rb;
	    rr = rg = rb = 0;
	    if ( reach_light )
	    {
		Vector_3<K> v = to_light.to_vector();
		v = v / sqrt(v.squared_length());
		shade = to_double( v * normal );
		if ( shade < 0.0 )
		    shade = 0.0;

		// recurse on reflections, if needed
		if ( it->is_reflective )
		{
		    Vector_3<K> v = - to_light.to_vector();
		    double c = - (normal * v);
		    Vector_3<K> nv = v + (2 * c * normal);
		    Ray_3<K> reflective_ray( ipoint, nv.direction() );
		    // recursive call
		    compute_color( scene, reflective_ray, rr, rg, rb );
		}
	    }

	    double factor = ambient_factor + (1 - ambient_factor) * shade;
	    double reflective_factor = 0.4;
	    
	    unsigned char r = (unsigned char)(it->color.r * factor * (1 - reflective_factor) + reflective_factor * rr);
	    unsigned char g = (unsigned char)(it->color.g * factor * (1 - reflective_factor) + reflective_factor * rg);
	    unsigned char b = (unsigned char)(it->color.b * factor * (1 - reflective_factor) + reflective_factor * rb);
	    unsigned long color = b | (g << 8) | (r << 16 );

	    colors[ to_double(ipoint.z()) ] = color;
	}

	if ( colors.size() == 0 )
	{
	    red = green = blue = 0;
	}
	else
	{
	    // get the first color (the closest)
	    unsigned long color = colors.begin()->second;
	    red = (color >> 16) & 0xff;
	    green = (color >> 8) & 0xff;
	    blue = color & 0xff;
	}
	recurse_level--;
    }

    void render( const Scene& scene, const std::string& out_file )
    {
	unsigned char *data = new unsigned char[ width_ * height_ * 3 ];

	for ( size_t y = 0; y < height_; y++ )
	{
	    for ( size_t x = 0; x < width_; x++ )
	    {
		uchar r, g, b;
		Point_3<K> pt = to_scene( x, y );
		Point_3<K> po( pt );
		po = po + Vector_3<K>( 0.0, 0.0, - cam_z_ );
		// The ray of returning light to the camera
		Ray_3<K> ray( pt, po );
		compute_color( scene, ray, r, g, b );

		unsigned long p = (y * width_ + x) * 3;
		data[p] = r;
		data[p+1] = g;
		data[p+2] = b;
	    }
	    cout << y << " / " << height_ << endl;
	}
	QImage img( data, width_, height_, width_ * 3, QImage::Format_RGB888 );
	img.save( out_file.c_str(), "PNG" );
    }
protected:
    Point_2<K> cam_ul_, cam_lr_;
    double cam_z_;
    size_t width_, height_;
};

int main( int argc, char *argv[] )
{
    // 256 x 256
    size_t resolution = 256;
    if ( argc > 1 )
    {
	for ( int i = 1; i < argc; i++ )
	{
	    string arg = argv[i];
	    // resolution
	    if ( arg == "-r" )
	    {
		if ( argc > i + 1 )
		{
		    sscanf( argv[++i], "%lu", &resolution );
		}
	    }
	}
    }
    Sphere_3<K> sp( ORIGIN, 3 );

    Point_2<K> camera_ul( -2.0, 2.0 );
    Point_2<K> camera_lr( 2.0, -2.0 );

    Scene scene;
    SceneObject obj;
    // obj.object = make_object(Sphere_3<K>(Point_3<K>( 0.0, -0.2, 0.0 ), 0.5));
    // obj.color.r = 255;
    // obj.color.g = 255;
    // obj.color.b = 255;
    // obj.is_reflective = true;
    // scene.objects.push_back( obj );
    
    obj.object = make_object(Sphere_3<K>(ORIGIN + Vector_3<K>( 1.0, 0.5, -0.4 ), 0.2 ));
    obj.color.r = 255;
    obj.color.g = 0;
    obj.color.b = 0;
    obj.is_reflective = false;
    scene.objects.push_back( obj );
    
    obj.object = make_object(Sphere_3<K>(Point_3<K>( -1.0, 0.8, 0.3 ), 0.2 ));
    obj.color.r = 0;
    obj.color.g = 120;
    obj.color.b = 255;
    obj.is_reflective = false;
    scene.objects.push_back( obj );
    
    obj.object = make_object(Sphere_3<K>(Point_3<K>( 0.0, 1.2, 0.1 ), 0.3 ));
    obj.color.r = 0;
    obj.color.g = 120;
    obj.color.b = 255;
    obj.is_reflective = false;
    scene.objects.push_back( obj );
    
    // obj.object = make_object(Plane_3<K>(Point_3<K>( 0.0, -0.5, 0.0) , Vector_3<K>( 0.0, 1.0, -0.5 )));
    // obj.color.r = 255;
    // obj.color.g = 0;
    // obj.color.b = 255;
    // scene.objects.push_back( obj );
    
    // TIN loading
    TIN tin;
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
	    for ( eit = it->facet_begin(); s < circulator_size( eit ); s++, eit++ )
	    {
		Point_3<K> p = eit->vertex()->point();
		tri[s] = Point_3<K>( p.x() * 3 - 1.5, p.y() * 3 - 1.5, p.z() * 3 - 1.5);
	    }
	    
	    Triangle_3<K> triangle( tri[0], tri[1], tri[2] );
	    if ( triangle.is_degenerate() )
	    {
		cout << triangle << endl;
		cout << "DEGENERATED triangle !" << endl;
		continue;
	    }
	    tin.push_back( triangle );
	}

	obj.object = make_object( tin );
	obj.color.r = 255;
	obj.color.g = 255;
	obj.color.b = 255;
	obj.is_reflective = true;
	scene.objects.push_back( obj );
    }

    scene.light.location = Point_3<K>( -0.5, 1.0, -1.5 );
    
    ProjectionScreen pscreen( camera_ul, camera_lr, -5.0, resolution, resolution );
    pscreen.set_ambient( 0.05 );
    pscreen.render( scene, "out.png" );
    return 0;
}
