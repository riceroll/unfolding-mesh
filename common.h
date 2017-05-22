//------------------------------------------------------------------------------
//	
//	common.h
//	
//------------------------------------------------------------------------------
#ifndef __COMMON_H_
#define __COMMON_H_

//#include <CGAL/basic.h>
//#include <CGAL/Filtered_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Origin.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/enum.h>

#include <boost/graph/adjacency_list.hpp>
// #include <boost/graph/graph_traits.hpp>


using namespace std;
using namespace boost;

//------------------------------------------------------------------------------
//	Macro definitions
//------------------------------------------------------------------------------
#ifndef TRUE
#define TRUE			(1)
#endif	// TRUE
#ifndef FALSE
#define FALSE			(0)
#endif	// FALSE

#define SQRT2			(1.414213562)

#ifndef INFINITY
#define	INFINITY		(1.0e+8)	
#endif	// INFINITY
#ifndef EPS
#define	EPS			(1.0e-8)
#endif	// EPS
#ifndef LARGEINT			
#define LARGEINT		(32767)
// #define LARGEINT		(2147483647)
#endif	// LARGEINT

#define MAX2( x, y )	(( (x)<(y) )? (y) : (x) )
#define MIN2( x, y )	(( (x)>(y) )? (y) : (x) )

#define SQUARE( x )	((x)*(x))

#define NO_INDEX	(-1)

#define DEFAULT_LABEL		(0)

#define DEFAULT_COLORS		(8)

// Labels for vertices
#define PEAK_VERTEX		(1)
#define SADDLE_VERTEX		(2)

// Number of cuts
#define PEAK_CUTS		(1)
#define SADDLE_CUTS		(2)

// Labels for edges
#define FIXED_EDGE		(1)
#define INVALID_EDGE		(2)

// Matches for edges
#define MATCHED_EDGE		(1)
#define UNMATCHED_EDGE		(2)

// default size of the paper
//#define DEFAULT_PAPER_WIDTH	(0.50)
//#define DEFAULT_PAPER_HEIGHT	(0.707106781)
#define DEFAULT_PAPER_WIDTH	(1.000000000)
#define DEFAULT_PAPER_HEIGHT	(1.414213562)
// #define DEFAULT_PAPER_WIDTH	(2.0500000000)
// #define DEFAULT_PAPER_HEIGHT	(2.05*1.414213562)
//#define DEFAULT_PAPER_WIDTH	(2.500000000)
//#define DEFAULT_PAPER_HEIGHT	(2.5*1.414213562)
//#define EXTEND_PAPER_WIDTH	(0.5)
//#define EXTEND_PAPER_HEIGHT	(0.707106781)
//#define EXTEND_PAPER_WIDTH	(0.1)
//#define EXTEND_PAPER_HEIGHT	(0.1414213562)
#define EXTEND_PAPER_WIDTH	(0.0050)
#define EXTEND_PAPER_HEIGHT	(0.00707106781)
//#define EXTEND_PAPER_WIDTH	(0.000)
//#define EXTEND_PAPER_HEIGHT	(0.000)
//#define EXTEND_PAPER_WIDTH	(1.0)
//#define EXTEND_PAPER_HEIGHT	(1.414213562)
#define MAXIMUM_PAPER_WIDTH	(20.00000000)
#define MAXIMUM_PAPER_HEIGHT	(28.28427124)

#define MARGIN_RATIO		(0.015)
// #define MARGIN_RATIO		(0.030)

// BUFFER SIZE
#define BUFFER_SIZE		(256)

// number of triangle sides
#define NUM_SIDES		(3)

#define SHRINKAGE_RATIO		(0.999)

//#define DEFAULT_DEVIATION	(0.0)// 50%
// #define DEFAULT_DEVIATION	(0.50)
#define DEFAULT_DEVIATION	(0.67)// 75%
//#define DEFAULT_DEVIATION	(1.00)

#define DETAILED_CHECK_LIMIT	(2)

#define LABEL_INSIDE_STRIP


// weight assignement for dual edges
typedef enum {
    UNIFORM_WEIGHT, RANDOM_WEIGHT, 
    MINIMUM_PERIMETER, MAXIMUM_PERIMETER, 
    FLAT_SPANNING, COIL_SPANNING,
    MINIMUM_DIHEDRAL, MAXIMUM_DIHEDRAL,
    MINIMUM_BLENDING, MAXIMUM_BLENDING,
    MINIMUM_CURVATURE, MAXIMUM_CURVATURE,
    MINIMUM_CONCAVITY, MAXIMUM_CONCAVITY,
    MINIMUM_HYPERBOLIC_AVE,
    MINIMUM_HYPERBOLIC_MAX,
    MINIMUM_HYPERBOLIC_MIN,
    SKELETON_FLAT, SKELETON_COIL
} Weight_type;


// maximum number of seam search
//#define MAX_SEAMS		(500)
#define MAX_SEAMS		(1000000)
//#define MAX_RECORDS		(500)
#define MAX_RECORDS		(1000000)

// maximum number of trials for resplit&remerge
#define MAX_RETRIALS		(64)
// #define MAX_RETRIALS		(32)

// angle for patch stlipping
#define INNER_PROD_DEGREE90	(0.0)
#define INNER_PROD_DEGREE60	(0.5)
#define INNER_PROD_ALL		(-1.001)

// Minimum spanning tree
#define KRUSKAL_EDGE (100)



#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif	// M_PI

//------------------------------------------------------------------------------
//	Additional templates from CGAL
//------------------------------------------------------------------------------
// Definition of My vertex
template < class Refs, class Traits, class Point >
class My_vertex : public CGAL::HalfedgeDS_vertex_max_base_with_id< Refs, Point, std::size_t > {
private:
    unsigned char	_label;
    int			_nCuts;
public:
    My_vertex() {
	_label		= DEFAULT_LABEL;
	_nCuts		= 0;
    };
    My_vertex( const Point & p ) : CGAL::HalfedgeDS_vertex_max_base_with_id< Refs, Point, std::size_t >( p ) {
	_label		= DEFAULT_LABEL;
	_nCuts		= 0;
    };
    My_vertex( const My_vertex & obj ) : CGAL::HalfedgeDS_vertex_max_base_with_id< Refs, Point, std::size_t >( obj ) {
	_label		= obj._label;
	_nCuts		= obj._nCuts;
    };
    const unsigned char & label( void ) const	{ return _label; }
    unsigned char & label( void )		{ return _label; } 
    const int & nCuts( void ) const		{ return _nCuts; }
    int & nCuts( void )				{ return _nCuts; } 

    My_vertex & operator = ( const My_vertex & obj ) {
	( CGAL::HalfedgeDS_vertex_max_base_with_id< Refs, Point, std::size_t > & )*this = obj;
	if ( this != &obj ) {
	    _label	= obj._label;
	    _nCuts	= obj._nCuts;
	}
	return *this;
    }
};

// Definition of My half edge
template < class Refs, class Traits, class Point >
class My_halfedge : public CGAL::HalfedgeDS_halfedge_base< Refs > {
private:
    int			_id;		// edge IDs (a half edge and its
					// opposite share the same edge ID.)
    unsigned char	_label;		// for labels of the edge
					// *** If I remember correctly, ...
					// FIXED_EDGE:	 edge is fixed as a cut
					// edge (perfect matching edge)
					// INVALID_EDGE: edge that cannot be
					// changed into a perfect matching edge
    unsigned char	_match;		// perfect matching and MST edge flag
    Point		_mid;		// coordinates of the edge midpoint for
					// drawing the dual graph
    double		_weight;	// weight value (currently unused)
    int			_cycle;		// patch ID
    bool		_connect;	// flag for face connection in unfoleded patches
    bool		_fixed;		// flag for fixing the weight value
    bool		_visit;		// for remerge&resplit process to avoid
					// selecting the same edges as the boundary
    int			_path;		// Boundary path ID for topological surgery
    bool		_orient;	// True if its orientation is positive

public:
    My_halfedge() : CGAL::HalfedgeDS_halfedge_base< Refs >() {
	_id		= NO_INDEX;
	_label		= DEFAULT_LABEL;
	_match		= DEFAULT_LABEL;
	_weight		= 0.0;
	_cycle		= NO_INDEX;
	_connect	= true;
	_fixed		= false;
	_visit		= false;
	_path		= NO_INDEX;
	_orient		= true;
    }
    My_halfedge( const My_halfedge & obj ) : CGAL::HalfedgeDS_halfedge_base< Refs >( obj ) {
	_id		= obj._id;
	_label		= obj._label;
	_match		= obj._match;
	_mid		= obj._mid;
	_weight		= obj._weight;
	_cycle		= obj._cycle;
	_connect	= obj._connect;
	_fixed		= obj._fixed;
	_visit		= obj._visit;
	_path		= obj._path;
	_orient		= obj._orient;
    }

    const int & id( void ) const		{ return _id; }
    int & id( void )				{ return _id; } 

    const unsigned char & label( void ) const	{ return _label; }
    unsigned char & label( void )		{ return _label; } 

    const unsigned char & match( void ) const	{ return _match; }
    unsigned char & match( void )		{ return _match; } 

    const Point & mid( void ) const		{ return _mid; }
    Point & mid( void )				{ return _mid; }

    const double & weight( void ) const		{ return _weight; }
    double & weight( void )			{ return _weight; } 

    const int & cycle( void ) const		{ return _cycle; }
    int & cycle( void )				{ return _cycle; } 

    const bool & connect( void ) const		{ return _connect; }
    bool & connect( void )			{ return _connect; }

    const bool & fixed( void ) const		{ return _fixed; }
    bool & fixed( void )			{ return _fixed; }

    const bool & visit( void ) const		{ return _visit; }
    bool & visit( void )			{ return _visit; }

    const int & path( void ) const		{ return _path; }
    int & path( void )				{ return _path; }

    const bool & orient( void ) const		{ return _orient; }
    bool & orient( void )			{ return _orient; }

    My_halfedge & operator = ( const My_halfedge & obj ) {
	( CGAL::HalfedgeDS_halfedge_base< Refs > & )*this = obj;
	if ( this != &obj ) {
	    _id		= obj._id;
	    _label	= obj._label;
	    _match	= obj._match;
	    _mid	= obj._mid;
	    _weight	= obj._weight;
	    _cycle	= obj._cycle;
	    _connect	= obj._connect;
	    _fixed	= obj._fixed;
	    _visit	= obj._visit;
	    _path	= obj._path;
	    _orient	= obj._orient;
	}
	return *this;
    }
};

// Definition of My face 
template < class Refs, class Traits, class Point, class Plane, class Triangle >
class My_face : public CGAL::HalfedgeDS_face_base< Refs, CGAL::Tag_true, Plane > {
private:
    int			_id;
    // unsigned char	_label;
    Point		_center;
    int			_piece;
    Triangle		_triangle;
public:
    My_face() : CGAL::HalfedgeDS_face_base< Refs, CGAL::Tag_true, Plane >() {
	_id		= NO_INDEX;
	// _label		= DEFAULT_LABEL;
	_piece		= NO_INDEX;
    }
    My_face( const My_face & obj ) : CGAL::HalfedgeDS_face_base< Refs, CGAL::Tag_true, Plane >( obj ) {
	_id		= obj._id;
	// _label		= obj._label;
	_center		= obj._center;
	_piece		= obj._piece;
	_triangle	= obj._triangle;
    }
    const int & id( void ) const		{ return _id; }
    int & id( void )				{ return _id; }
    // const unsigned char & label( void ) const	{ return _label; }
    // unsigned char & label( void )		{ return _label; } 
    const Point & center( void ) const		{ return _center; }
    Point & center( void )			{ return _center; }
    const int & piece( void ) const		{ return _piece; }
    int & piece( void )				{ return _piece; }
    const Triangle & triangle( void ) const	{ return _triangle; }
    Triangle & triangle( void )			{ return _triangle; }

    My_face & operator = ( const My_face & obj ) {
	( CGAL::HalfedgeDS_face_base< Refs, CGAL::Tag_true, Plane > & )*this = obj;
	if ( this != &obj ) {
	    _id		= obj._id;
	    // _label	= obj._label;
	    _center	= obj._center;
	    _piece	= obj._piece;
	    _triangle	= obj._triangle;
	}
	return *this;
    }
};

struct My_items : public CGAL::Polyhedron_items_with_id_3 {
    
    template < class Refs, class Traits >

    struct Vertex_wrapper {
        typedef typename Traits::Point_3 Point;
	typedef My_vertex< Refs, Traits, Point > Vertex;
    };

    template < class Refs, class Traits >
    struct Halfedge_wrapper {
        typedef typename Traits::Point_3 Point;
	typedef My_halfedge< Refs, Traits, Point > Halfedge;
    };

    template < class Refs, class Traits >
    struct Face_wrapper {
        typedef typename Traits::Point_3 Point;
	typedef typename Traits::Plane_3 Plane;
	typedef typename Traits::Triangle_2 Triangle;
	typedef My_face< Refs, Traits, Point, Plane, Triangle > Face;
    };

};


class Attribute {
private:
    // Depends only on 3D shape
    double		_totalArea;
    // Depends on 2D unfolded patches
    int			_nRuns;
    int			_nColors;
    int			_nStars;
    int			_nHyperbolics;
    int			_maxLength;
    double		_sumWeights;
    // Dual edge distance
    double		_aveDist;
    double		_maxDist;
public:
    Attribute() {
	_totalArea	= 0.0;
	//
	_nRuns		= 0;
	_nColors	= 0;
	_nStars		= 0;
	_nHyperbolics	= 0;
	_maxLength	= 0;
	_sumWeights	= 0.0;
	//
	_aveDist	= 0.0;
	_maxDist	= 0.0;
    }
    Attribute( const Attribute & obj ) {
	_totalArea	= obj._totalArea;
	//
	_nRuns		= obj._nRuns;
	_nColors	= obj._nColors;
	_nStars		= obj._nStars;
	_nHyperbolics	= obj._nHyperbolics;
	_maxLength	= obj._maxLength;
	_sumWeights	= obj._sumWeights;
	//
	_aveDist	= obj._aveDist;
	_maxDist	= obj._maxDist;
    }

    void init( void ) {
	_totalArea	= 0.0;
	//
	_nRuns		= 0;
	_nColors	= 0;
	_nStars		= 0;
	_nHyperbolics	= 0;
	_maxLength	= 0;
	_sumWeights	= 0.0;
	//
	_aveDist	= 0.0;
	_maxDist	= 0.0;
    }

    const double & totalArea( void ) const	{ return _totalArea; }
    double & totalArea( void )			{ return _totalArea; }

    const int & nRuns( void ) const		{ return _nRuns; }
    int & nRuns( void )				{ return _nRuns; }

    const int & nColors( void ) const		{ return _nColors; }
    int & nColors( void )			{ return _nColors; }

    const int & nStars( void ) const		{ return _nStars; }
    int & nStars( void )			{ return _nStars; }

    const int & nHyperbolics( void ) const	{ return _nHyperbolics; }
    int & nHyperbolics( void )			{ return _nHyperbolics; }

    const int & maxLength( void ) const		{ return _maxLength; }
    int & maxLength( void )			{ return _maxLength; }

    const double & sumWeights( void ) const	{ return _sumWeights; }
    double & sumWeights( void ) 		{ return _sumWeights; }

    const double & aveDist( void ) const	{ return _aveDist; }
    double & aveDist( void )			{ return _aveDist; }

    const double & maxDist( void ) const	{ return _maxDist; }
    double & maxDist( void )			{ return _maxDist; }

    Attribute & operator = ( const Attribute & obj ) {
	if ( this != &obj ) {
	    _totalArea		= obj._totalArea;
	    //
	    _nRuns		= obj._nRuns;
	    _nColors		= obj._nColors;
	    _nStars		= obj._nStars;
	    _nHyperbolics	= obj._nHyperbolics;
	    _maxLength		= obj._maxLength;
	    _sumWeights		= obj._sumWeights;
	    //
	    _aveDist		= obj._aveDist;
	    _maxDist		= obj._maxDist;
	}
	return *this;
    }
};


//------------------------------------------------------------------------------
//	CGAL type definitions
//------------------------------------------------------------------------------
typedef CGAL::Cartesian< double >			Kernel;
typedef Kernel::Vector_3				Vector3;
typedef Kernel::Point_3					Point3;
typedef CGAL::Polyhedron_3< Kernel, My_items >		Polyhedron;

typedef Polyhedron::Vertex				Vertex;
typedef Polyhedron::Vertex_handle			Vertex_handle;
typedef Polyhedron::Vertex_iterator			Vertex_iterator;

typedef Polyhedron::Halfedge_handle			Halfedge_handle;
typedef Polyhedron::Halfedge_iterator			Halfedge_iterator;
typedef Polyhedron::Halfedge_around_vertex_circulator	Halfedge_vertex_circulator;
typedef Polyhedron::Halfedge_around_facet_circulator	Halfedge_facet_circulator;

typedef Polyhedron::Facet				Facet;
typedef Polyhedron::Facet_handle			Facet_handle;
typedef Polyhedron::Facet_iterator			Facet_iterator;



typedef Kernel::Vector_2				Vector2;
typedef Kernel::Point_2					Point2;


typedef Kernel::Triangle_2				Triangle2;
typedef Kernel::Triangle_3				Triangle3;
typedef Kernel::Line_2					Line2;
typedef Kernel::Line_3					Line3;
typedef Kernel::Segment_2				Segment2;
typedef Kernel::Segment_3				Segment3;
typedef CGAL::Aff_transformation_2< Kernel >		Transformation2;
typedef CGAL::Aff_transformation_3< Kernel >		Transformation3;
typedef CGAL::Bbox_2					Bbox2;

typedef CGAL::Object					Object;
//typedef Kernel::Intersect_3				Intersect3;


//------------------------------------------------------------------------------
//	BGL type definitions
//------------------------------------------------------------------------------
// typedef adjacency_list< vecS, vecS, undirectedS >	Graph;
typedef adjacency_list< vecS, vecS, undirectedS,
			no_property,
			property< edge_weight_t, double > > Graph;
typedef graph_traits< Graph > GraphTraits;
typedef graph_traits< Graph >::vertex_descriptor VertexDescriptor;
typedef graph_traits< Graph >::edge_descriptor EdgeDescriptor;
typedef pair< int, int > Edge;
typedef graph_traits< Graph >::vertex_iterator VertexIterator;
typedef graph_traits< Graph >::edge_iterator EdgeIterator;
typedef property_map< Graph, vertex_index_t >::type VertexIDMap;
typedef property_map< Graph, edge_weight_t >::type EdgeIDMap;


//------------------------------------------------------------------------------
//	local type definitions
//------------------------------------------------------------------------------

struct Outline {
private:    
    Halfedge_handle		_start;
    vector< Halfedge_handle >	_seam;
    vector< Halfedge_handle >	_tape;
public:
    Outline() {
	_seam.clear();
	_tape.clear();
    }
    Outline( const Outline & obj ) {
	_start	= obj._start;
	_seam	= obj._seam;
	_tape	= obj._tape;
    }
    const Halfedge_handle & start( void ) const			{ return _start; }
    Halfedge_handle & start( void )				{ return _start; }
    const vector< Halfedge_handle > & seam( void ) const	{ return _seam; }
    vector< Halfedge_handle > & seam( void )			{ return _seam; }
    const vector< Halfedge_handle > & tape( void ) const	{ return _tape; }
    vector< Halfedge_handle > & tape( void )			{ return _tape; }
};








//------------------------------------------------------------------------------
//	Functions in dual.cpp
//------------------------------------------------------------------------------
extern void printHalfedge	( const Halfedge_handle & hh );
extern void normalizeMesh	( Polyhedron & poly, vector< Segment3 > & line );
extern void initAttrs		( Polyhedron & poly );
extern void renewAttrs		( Polyhedron & poly );
extern void initWeights		( Polyhedron & poly );
extern void computeDual		( Polyhedron & poly, Graph & dual,
				  vector< Point3 > & bary, vector< Point3 > & mid ); 
extern Vector3 principalAxis	( Polyhedron & poly );
extern void rotateMesh		( Polyhedron & poly, char axis, double radius );
extern void alignMesh		( Polyhedron & poly );
extern void smoothMesh		( Polyhedron & poly, unsigned int nTimes = 1 );

//------------------------------------------------------------------------------
//	Functions in feature.cpp
//------------------------------------------------------------------------------
extern double	centralAngle	( Vertex_handle	& vh );
extern bool	isSaddle	( Vertex_handle & vh );

extern double	degreeOfValley	( const Halfedge_handle & hh );
extern double	evalConvexity	( const Halfedge_handle & hh );
extern bool	isConvex	( const Halfedge_handle & hh );

extern void	assignWeights	( Polyhedron & poly, Graph & dual, Weight_type sw );

extern void	weightStatistics( Polyhedron & poly, double & ave, double & stdev );

extern double	totalArea	( Polyhedron & poly );

extern void	labelByCurvature( Polyhedron & poly );
extern bool	fixDual		( Graph & dual, Halfedge_handle & hh );
extern bool	fixWeight	( Graph & dual, Halfedge_handle & hh, double value );
extern void	applyConstraints( Polyhedron & poly, Graph & dual );

//------------------------------------------------------------------------------
//	Functions in matching.cpp
//------------------------------------------------------------------------------
extern void	singleFaces	( Polyhedron & poly,
				  vector< vector< Halfedge_handle > > & spine );
extern bool	isPerfect	( Graph & graph );
extern void	perfectMatching	( Graph & dual, Polyhedron & poly,
				  vector< vector< Halfedge_handle > > & spine );
extern void	weightedMatching( Graph & dual, Polyhedron & poly,
				  vector< vector< Halfedge_handle > > & spine );



//------------------------------------------------------------------------------
//	Functions in craft.cpp
//------------------------------------------------------------------------------
// geometrical
extern Point2	centerTriangle	( const Triangle2 & t );
extern bool	isIntersected	( const Triangle2 & triA, const Triangle2 & triB );
extern bool	isIntersected	( const Triangle2 & tri, const vector< Triangle2 > & strip );

// topological
extern Halfedge_handle
		sharedHalfedge	( const Facet_handle & fhI, const Facet_handle & fhO );
extern void	countNCuts	( Polyhedron & poly );
extern int	extractBoundary	( Polyhedron & poly, vector< vector< Halfedge_handle > > & bound );
extern void	updateAttr	( Polyhedron & poly, 
				  vector< vector< Halfedge_handle > > & bound,
				  int & nRuns, Attribute & attr );
extern void	traverseBoundary( Polyhedron & poly, Attribute & attr, vector< vector< Halfedge_handle > > & bound );
extern void	traverseBoundary( Polyhedron & poly, Attribute & attr );
extern void	labelBoundary	( Polyhedron & poly, Attribute & attr );

// function for traversing dual meshes
extern void	findIntersections
				( vector< Facet_handle > & vfh,
				  vector< pair< Facet_handle, Facet_handle > > & endfaces ) ;
extern bool	findDualPath	( const Facet_handle & fhO, // origin face of the path
				  const Facet_handle & fhD, // destination face of the path
				  Polyhedron & poly,
				  vector< Facet_handle > & path );
extern unsigned int
		dualLength	( const Facet_handle & fhO, // origin face of the path
				  const Facet_handle & fhD, // destination face of the path
				  Polyhedron & poly );

// boolean function
extern void	minimumCover	( Polyhedron & poly,
				  vector< pair< Facet_handle, Facet_handle > > & endfaces,
				  vector< vector< Facet_handle > > & pathlist,
				  vector< Halfedge_handle > & cut );
extern bool	constrainedMinimumCover
				( Polyhedron & poly,
				  vector< pair< Facet_handle, Facet_handle > > & endfaces,
				  vector< vector< Facet_handle > > & pathlist,
				  vector< Halfedge_handle > & cut );
extern void	cutStrips	( vector< vector< Facet_handle > > & patch,
				  Halfedge_handle & cut );
// boolean function
extern void	resolveOverlaps	( Polyhedron & poly,
				  vector< vector< Facet_handle > > & patch );
extern bool	resplitOverlaps	( Polyhedron & poly,
				  vector< vector< Facet_handle > > & patch );
extern bool	resplitOverlaps	( Polyhedron & poly,
				  vector< vector< Facet_handle > > & patch,
				  vector< Facet_handle > & setA, vector< Facet_handle > & setB );




//------------------------------------------------------------------------------
//	Functions in project.cpp
//------------------------------------------------------------------------------
extern void	rewindTriangle	( Halfedge_handle & oldbase, const Halfedge_handle & newbase );
extern void	resetCycleIDs	( vector< vector< Facet_handle > > & patch );
extern void	boundCycles	( vector< int > & cycleID,
				  vector< vector< Facet_handle > > & patch,
				  vector< Bbox2 > & bound );
extern void	layoutSingle	( vector< int > & cycleID,
				  vector< vector< Facet_handle > > & patch,
				  Bbox2 & sheet, vector< Bbox2 > & bound );
extern void	layoutMultiple	( vector< int > & cycleID,
				  vector< vector< Facet_handle > > & patch,
				  Bbox2 & sheet, vector< Bbox2 > & bound );
extern void	arrangeSingle	( vector< vector< Facet_handle > > & patch,
				  Bbox2 & sheet, vector< Bbox2 > & bound );
extern void	arrangeMultiple	( vector< vector< Facet_handle > > & patch,
				  Bbox2 & sheet, vector< Bbox2 > & bound );


extern void	projectCycles	( vector< vector< Halfedge_handle > > & cycle,
				  vector< vector< Facet_handle > > & patch,
				  Bbox2 & sheet, vector< Bbox2 > & bound,
				  double limitCos );

//------------------------------------------------------------------------------
//	Functions in stitch.cpp
//------------------------------------------------------------------------------
extern bool	isLocally	( const Halfedge_handle & hh );
extern void	transformStrip	( Halfedge_handle & bridge,
				  vector< vector< Facet_handle > > & patch,
				  vector< Triangle2 > & transformed );
extern bool	isGlobally	( vector< Facet_handle > & fixed,
				  vector< Triangle2 > & transformed );
extern bool	doEmbed		( Halfedge_handle & bridge,
				  vector< vector< Facet_handle > > & patch,
				  vector< Triangle2 > & transformed );
extern bool	canEmbed	( Halfedge_handle & bridge,
				  vector< vector< Facet_handle > > & patch,
				  vector< Triangle2 > & transformed );
extern void	stitch		( Halfedge_handle & bridge,
				  vector< vector< Facet_handle > > & patch,
				  vector< Triangle2 > & transformed );
extern bool	mergeByOne	( Polyhedron & poly, 
				  Bbox2 & sheet,
				  vector< vector< Facet_handle > > & patch,
				  vector< Bbox2 > & bound );
extern void	greedy		( Polyhedron & poly, 
				  Bbox2 & sheet,
				  vector< vector< Facet_handle > > & patch,
				  vector< Bbox2 > & bound );
extern void	reconnectByOne	( Polyhedron & poly, 
				  Bbox2 & sheet,
				  vector< vector< Facet_handle > > & patch,
				  vector< Bbox2 > & bound );
extern void	reconnect	( Polyhedron & poly, 
				  Bbox2 & sheet,
				  vector< vector< Facet_handle > > & patch,
				  vector< Bbox2 > & bound,
				  Attribute & attr,
				  void (*preview)( void ) );

//------------------------------------------------------------------------------
//	Functions in search.cpp
//------------------------------------------------------------------------------
extern void	divide		( Polyhedron & poly,
				  Bbox2 & sheet,
				  vector< vector< Facet_handle > > & patch,
				  vector< Bbox2 > & bound );

//------------------------------------------------------------------------------
//	Functions in gene.cpp
//------------------------------------------------------------------------------
extern int	optimize	( Polyhedron & poly, 
				  Bbox2 & sheet,
				  vector< vector< Facet_handle > > & patch,
				  vector< Bbox2 > & bound,
				  Attribute & attr,
				  void (*preview)( void ) );
extern int	exhaustive	( Polyhedron & poly, 
				  Bbox2 & sheet,
				  vector< vector< Facet_handle > > & patch,
				  vector< Bbox2 > & bound,
				  Attribute & attr,
				  void (*preview)( void ) );



//------------------------------------------------------------------------------
//	Functions in mst.cpp
//------------------------------------------------------------------------------
// Bbox2 & sheet, vector< Bbox2 > & bound, double limit );

extern void	MSTSplit	( Graph & tree, Polyhedron & poly,
				  vector< vector< Halfedge_handle > > & spine,
				  vector< Halfedge_handle > & joint );
extern void	MSTMerge	( Polyhedron & poly,
				  vector< Halfedge_handle > & joint,
				  vector< vector< Facet_handle > > & patch,
				  Bbox2 & sheet, vector< Bbox2 > & bound );
extern void	MSTDissolve	( Polyhedron & poly,
				  // vector< vector< Halfedge_handle > > & spine,
				  vector< vector< Facet_handle > > & patch,
				  Bbox2 & sheet, vector< Bbox2 > & bound );

extern void	MSTStrips	( Graph & tree, Polyhedron & poly,
				  vector< vector< Halfedge_handle > > & spine,
				  vector< vector< Facet_handle > > & patch,
				  double limitCos,
				  Bbox2 & sheet, vector< Bbox2 > & bound );
extern void	MSTPieces	( Graph & tree, Polyhedron & poly,
				  vector< vector< Halfedge_handle > > & spine,
				  vector< vector< Facet_handle > > & patch,
				  double limitCos,
				  Bbox2 & sheet, vector< Bbox2 > & bound );
extern void	minimumSpanningTree
				( Graph & tree, Polyhedron & poly,
				  vector< vector< Facet_handle > > & patch,
				  Bbox2 & sheet, vector< Bbox2 > & bound );
extern void	spanningTrees	( Graph & tree, Polyhedron & poly,
				  vector< vector< Facet_handle > > & patch,
				  Bbox2 & sheet, vector< Bbox2 > & bound );
#endif	// __COMMON_H_
