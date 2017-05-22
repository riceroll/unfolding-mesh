//------------------------------------------------------------------------------
//
//	feature.cpp
//
//------------------------------------------------------------------------------
#include <iostream>
#include <stack>
#include <queue>
#include <cmath>
#include <cassert>

using namespace std;

#include "common.h"
#include "ui.h"

//------------------------------------------------------------------------------
//	Functions for vertices
//------------------------------------------------------------------------------
double centralAngle( Vertex_handle & vh )
{
#ifdef OLD
    double oldsum = 0.0;
    Halfedge_vertex_circulator hvc = vh->vertex_begin();
    do {
	Point3 pointL = hvc->opposite()->vertex()->point();
	Point3 pointC = vh->point();
	Point3 pointR = hvc->next()->vertex()->point();
	// cerr << "pointL = " << hvc->opposite()->vertex()->id() << " : " << pointL << endl;
	// cerr << "pointC = " << vh->id() << " : " << pointC << endl;
	// cerr << "pointR = " << hvc->next()->vertex()->id() << " : " << pointR << endl;
	Vector3 wingL = pointL - pointC;
	Vector3 wingR = pointR - pointC;
	double lengthL = sqrt( wingL.squared_length() );
	double lengthR = sqrt( wingR.squared_length() );
	double innerProd = ( wingL * wingR ) / ( lengthL*lengthR );
	oldsum += acos( innerProd );
    } while ( ++hvc != vh->vertex_begin() );
#endif	// OLD

    double sum_in_radian = 0.0;
    Halfedge_vertex_circulator hvc = vh->vertex_begin();
    do {
	Point3 pointR = hvc->opposite()->vertex()->point();
	Point3 pointC = vh->point();
	Point3 pointL = (--hvc)->opposite()->vertex()->point();
#ifdef MYDEBUG
	cerr << "pointR = " << pointR << endl;
	cerr << "pointC = " << pointC << endl;
	cerr << "pointL = " << pointL << endl;
#endif	// MYDEBUG
	Vector3 wingL = pointL - pointC;
	Vector3 wingR = pointR - pointC;
	double lengthL = sqrt( wingL.squared_length() );
	double lengthR = sqrt( wingR.squared_length() );
	double innerProd = ( wingL * wingR ) / ( lengthL*lengthR );
	sum_in_radian += acos( innerProd );
    } while ( hvc != vh->vertex_begin() );

    return ( sum_in_radian );
}


bool isPeak( Vertex_handle & vh )
{
    const double peakAngle = M_PI;
    double radian = centralAngle( vh );
    if ( radian < peakAngle ) {
	return TRUE;
    }
    else {
	return FALSE;
    }
}


bool isSaddle( Vertex_handle & vh )
{
    const double eps = 1.0e-4;
    const double flatAngle = 2.0 * M_PI + eps;
    double radian = centralAngle( vh );
    if ( radian > flatAngle ) {
	return TRUE;
    }
    else {
	return FALSE;
    }
}


//------------------------------------------------------------------------------
//	Functions for half edges
//------------------------------------------------------------------------------
double evalConvexity( const Halfedge_handle & hh )
{
#ifdef OLD
    Halfedge_handle hhL, hhR, hhC, hhO;
    hhC = hh;
    hhO = hhC->opposite();

    Vector3 vecC = hhC->opposite()->vertex()->point() - hhC->mid();
    Vector3 vecL = hhC->facet()->center() - hhC->mid();
    Vector3 vecR = hhO->facet()->center() - hhO->mid();
    Vector3 nmL = cross_product( vecC, vecL );
    Vector3 nmR = cross_product( vecR, vecC );
    nmL = nmL / sqrt( nmL.squared_length() );
    nmR = nmR / sqrt( nmR.squared_length() );
    Vector3 orient = cross_product( nmR, nmL );
    double signdir = orient * vecC;
    double innerProd = nmL * nmR;
    double sineAngle = sqrt( 1.0 - innerProd*innerProd );

    if ( signdir < 0 ) sineAngle *= ( -1.0 );

    return sineAngle;
#endif	// OLD

    Halfedge_handle hhF, hhB;
    hhF = hh;
    hhB = hhF->opposite();

    Point3  mid  = hhF->mid();
    Vector3 vecC = hhF->vertex()->point() - mid;
    Vector3 vecL = hhF->facet()->center() - mid;
    Vector3 vecR = hhB->facet()->center() - mid;
    Vector3 nmL = cross_product( vecC, vecL );
    Vector3 nmR = cross_product( vecR, vecC );
    nmL = nmL / sqrt( nmL.squared_length() );
    nmR = nmR / sqrt( nmR.squared_length() );
    Vector3 orient = cross_product( nmL, nmR );
    double signdir = orient * vecC;
    double innerProd = nmL * nmR;
    double sharpness = 0.5 * ( 1.0 - innerProd );

    if ( signdir < 0 ) sharpness *= ( -1.0 );


#ifdef DEBUG
    cerr << " nmL = " << nmL << endl;
    cerr << " nmR = " << nmR << endl;
    cerr << " innerProd = " << innerProd << endl;
    cerr << " sineAngle = " << sineAngle << endl;
#endif	// DEBUG
    
    return sharpness;

    // return ( fabs( sineAngle ) - 0.35 );
}


bool isConvex( const Halfedge_handle & hh )
{
    double innerProd = evalConvexity( hh );

#ifdef DEBUG
    const int idA = 30;
    const int idB = 53;
    if ( ( ( hhC->vertex()->id() == idA ) && 
	   ( hhC->opposite()->vertex()->id() == idB ) ) ||
	 ( ( hhC->vertex()->id() == idB ) && 
	   ( hhC->opposite()->vertex()->id() == idA ) ) ) {

	cerr << "##############################" << endl;
	
	cerr << " hhC : " << hhC->vertex()->id() << " -- " << hhC->opposite()->vertex()->id() << " ** ";
	cerr << " hhC : " << hhC->vertex()->point() << " == " << hhC->opposite()->vertex()->point() << endl;
	
	cerr << " vecC = " << vecC << endl;
	cerr << " vecL = " << vecL << endl;
	cerr << " vecR = " << vecR << endl;
	cerr << " nmL = " << nmL << endl;
	cerr << " nmR = " << nmR << endl;
	cerr << " innerProd = " << innerProd << endl;
    }
#endif	// DEBUG

    if ( innerProd >= 0.0 ) return true;
    else return false;
}


//------------------------------------------------------------------------------
//	Edge weight assignement for the MST computation
//------------------------------------------------------------------------------
void uniformWeights( Polyhedron & poly, Graph & dual )
{
    unsigned int nEdges = num_edges( dual );
    vector< double > weight( nEdges );
	
    // Initialization of edge connections
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	hi->connect() = false;
    }

    // uniform weights
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	weight[ eID ] = 1.0;
	// Special handling
	if ( hi->fixed() ) weight[ eID ] = hi->weight();
	else hi->weight() = hi->opposite()->weight() = weight[ eID ];
    }

    EdgeIDMap weightMap = get( edge_weight, dual );
	
    boost::graph_traits< Graph >::edge_iterator ei, ei_end;
    int count = 0;
    for( tie( ei, ei_end ) = edges( dual ); ei != ei_end; ++ei ) {
	weightMap[ *ei ] = weight[ count++ ];
    }
}


void randomWeights( Polyhedron & poly, Graph & dual )
{
    unsigned int nEdges = num_edges( dual );
    vector< double > weight( nEdges );
	
    // Initialization of edge connections
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	hi->connect() = false;
    }

    // random weights
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	weight[ eID ] = ( double )( rand()%20+1 )/ 20.0;
	// Special handling
	if ( hi->fixed() ) weight[ eID ] = hi->weight();
	else hi->weight() = hi->opposite()->weight() = weight[ eID ];
    }

    EdgeIDMap weightMap = get( edge_weight, dual );
	
    boost::graph_traits< Graph >::edge_iterator ei, ei_end;
    int count = 0;
    for( tie( ei, ei_end ) = edges( dual ); ei != ei_end; ++ei ) {
	//std::cout << *ei << ", ";
	weightMap[ *ei ] = weight[ count++ ];
	// cerr << " i = " << *ei << " weight = " << weightMap[ *ei ] << endl;
    }
}


//------------------------------------------------------------------------------
//	
void minimumPerimeter( Polyhedron & poly, Graph & dual )
{
    unsigned int nEdges = num_edges( dual );
    vector< double > weight( nEdges );
    vector< double > length( nEdges );
    double minLength = 1000000.0, maxLength = -1000000.0;
	
    // Initialization of edge connections
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	hi->connect() = false;
    }

    // Perimeter heuristic weights
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	Vector3 vec = hi->vertex()->point() - hi->opposite()->vertex()->point();
	length[ eID ] = sqrt( vec.squared_length() );
	if ( length[ eID ] < minLength ) minLength = length[ eID ];
	if ( length[ eID ] > maxLength ) maxLength = length[ eID ];
    }
    cerr << " minLength = " << minLength << " maxLength = " << maxLength << endl;
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	assert( length[ eID ] > 0.0 );
	weight[ eID ] = 2.0 * ( maxLength - length[ eID ] ) / ( maxLength - minLength ) - 1.0;
	// Special handling
	if ( hi->fixed() ) weight[ eID ] = hi->weight();
	else hi->weight() = hi->opposite()->weight() = weight[ eID ];
    }

    EdgeIDMap weightMap = get( edge_weight, dual );
	
    boost::graph_traits< Graph >::edge_iterator ei, ei_end;
    int count = 0;
    for( tie( ei, ei_end ) = edges( dual ); ei != ei_end; ++ei ) {
	//std::cout << *ei << ", ";
	weightMap[ *ei ] = weight[ count++ ];
	// cerr << " i = " << *ei << " weight = " << weightMap[ *ei ] << endl;
    }
}


void maximumPerimeter( Polyhedron & poly, Graph & dual )
{
    unsigned int nEdges = num_edges( dual );
    vector< double > weight( nEdges );
    vector< double > length( nEdges );
    double minLength = 1000000.0, maxLength = -1000000.0;
	
    // Initialization of edge connections
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	hi->connect() = false;
    }

    // Perimeter heuristic weights
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	Vector3 vec = hi->vertex()->point() - hi->opposite()->vertex()->point();
	length[ eID ] = sqrt( vec.squared_length() );
	if ( length[ eID ] < minLength ) minLength = length[ eID ];
	if ( length[ eID ] > maxLength ) maxLength = length[ eID ];
    }
    cerr << " minLength = " << minLength << " maxLength = " << maxLength << endl;
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	assert( length[ eID ] > 0.0 );
	weight[ eID ] = 2.0 * ( length[ eID ] - minLength ) / ( maxLength - minLength ) - 1.0;
	// Special handling
	if ( hi->fixed() ) weight[ eID ] = hi->weight();
	else hi->weight() = hi->opposite()->weight() = weight[ eID ];
    }

    EdgeIDMap weightMap = get( edge_weight, dual );
	
    boost::graph_traits< Graph >::edge_iterator ei, ei_end;
    int count = 0;
    for( tie( ei, ei_end ) = edges( dual ); ei != ei_end; ++ei ) {
	//std::cout << *ei << ", ";
	weightMap[ *ei ] = weight[ count++ ];
	// cerr << " i = " << *ei << " weight = " << weightMap[ *ei ] << endl;
    }
}


//------------------------------------------------------------------------------
//	
void flatSpanning( Polyhedron & poly, Graph & dual )
{
    unsigned int nEdges = num_edges( dual );
    vector< double > weight( nEdges );
    static Vector3 orient( 0.0, 0.0, 1.0 );
	
    // Initialization of edge connections
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	hi->connect() = false;
    }

    // flat spanning heuristic weights
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	Vector3 vec = hi->vertex()->point() - hi->opposite()->vertex()->point();
	weight[ eID ] = 2.0 * fabs( orient * vec )/sqrt( vec.squared_length() ) - 1.0;
	// Special handling
	if ( hi->fixed() ) weight[ eID ] = hi->weight();
	else hi->weight() = hi->opposite()->weight() = weight[ eID ];
    }

    for ( unsigned int i = 0; i < nEdges; ++i ) {
	assert( ( 0.0 <= weight[ i ] ) && ( weight[ i ] <= 1.0 ) );
    }

    EdgeIDMap weightMap = get( edge_weight, dual );
	
    boost::graph_traits< Graph >::edge_iterator ei, ei_end;
    int count = 0;
    for( tie( ei, ei_end ) = edges( dual ); ei != ei_end; ++ei ) {
	//std::cout << *ei << ", ";
	weightMap[ *ei ] = weight[ count++ ];
	// cerr << " i = " << *ei << " weight = " << weightMap[ *ei ] << endl;
    }
}


void coilSpanning( Polyhedron & poly, Graph & dual )
{
    unsigned int nEdges = num_edges( dual );
    vector< double > weight( nEdges );
    static Vector3 orient( 0.0, 0.0, 1.0 );
	
    // Initialization of edge connections
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	hi->connect() = false;
    }

    // flat spanning heuristic weights
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	Vector3 vec = hi->vertex()->point() - hi->opposite()->vertex()->point();
	weight[ eID ] = 1.0 - 2.0 * fabs( orient * vec )/sqrt( vec.squared_length() );
	// Special handling
	if ( hi->fixed() ) weight[ eID ] = hi->weight();
	else hi->weight() = hi->opposite()->weight() = weight[ eID ];
    }

    for ( unsigned int i = 0; i < nEdges; ++i ) {
	assert( ( 0.0 <= weight[ i ] ) && ( weight[ i ] <= 1.0 ) );
    }

    EdgeIDMap weightMap = get( edge_weight, dual );
	
    boost::graph_traits< Graph >::edge_iterator ei, ei_end;
    int count = 0;
    for( tie( ei, ei_end ) = edges( dual ); ei != ei_end; ++ei ) {
	//std::cout << *ei << ", ";
	weightMap[ *ei ] = weight[ count++ ];
	// cerr << " i = " << *ei << " weight = " << weightMap[ *ei ] << endl;
    }
}


//------------------------------------------------------------------------------
//	
// Flat regions are more likely to be split  
void minimumDihedral( Polyhedron & poly, Graph & dual )
{
    unsigned int nEdges = num_edges( dual );
    vector< double > weight( nEdges );
    static Vector3 orient( 0.0, 0.0, 1.0 );
	
    // Initialization of edge connections
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	hi->connect() = false;
    }

    // flat spanning heuristic weights
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	Point3	midE	= hi->mid();
	Point3	cenFL	= hi->facet()->center();
	Point3	cenFR	= hi->opposite()->facet()->center();
	Vector3 vecL	= cenFL - midE;
	Vector3 vecR	= cenFR - midE;
	double	iprod	= ( vecL * vecR )/sqrt( vecL.squared_length() * vecR.squared_length() );
	weight[ eID ]	= 1.0 - iprod;
	// Special handling
	if ( hi->fixed() ) weight[ eID ] = hi->weight();
	else hi->weight() = hi->opposite()->weight() = weight[ eID ];
    }

    for ( unsigned int i = 0; i < nEdges; ++i ) {
	assert( ( 0.0 <= weight[ i ] ) && ( weight[ i ] <= 1.0 ) );
    }

    EdgeIDMap weightMap = get( edge_weight, dual );
	
    boost::graph_traits< Graph >::edge_iterator ei, ei_end;
    int count = 0;
    for( tie( ei, ei_end ) = edges( dual ); ei != ei_end; ++ei ) {
	//std::cout << *ei << ", ";
	weightMap[ *ei ] = weight[ count++ ];
	// cerr << " i = " << *ei << " weight = " << weightMap[ *ei ] << endl;
    }
}

// Sharp regions are more likely to be split  
void maximumDihedral( Polyhedron & poly, Graph & dual )
{
    unsigned int nEdges = num_edges( dual );
    vector< double > weight( nEdges );
    static Vector3 orient( 0.0, 0.0, 1.0 );
	
    // Initialization of edge connections
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	hi->connect() = false;
    }

    // flat spanning heuristic weights
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	Point3	midE	= hi->mid();
	Point3	cenFL	= hi->facet()->center();
	Point3	cenFR	= hi->opposite()->facet()->center();
	Vector3 vecL	= cenFL - midE;
	Vector3 vecR	= cenFR - midE;
	double	iprod	= ( vecL * vecR )/sqrt( vecL.squared_length() * vecR.squared_length() );
	weight[ eID ]	= iprod + 1.0;
	// Special handling
	if ( hi->fixed() ) weight[ eID ] = hi->weight();
	else hi->weight() = hi->opposite()->weight() = weight[ eID ];
    }

    for ( unsigned int i = 0; i < nEdges; ++i ) {
	assert( ( 0.0 <= weight[ i ] ) && ( weight[ i ] <= 1.0 ) );
    }

    EdgeIDMap weightMap = get( edge_weight, dual );
	
    boost::graph_traits< Graph >::edge_iterator ei, ei_end;
    int count = 0;
    for( tie( ei, ei_end ) = edges( dual ); ei != ei_end; ++ei ) {
	//std::cout << *ei << ", ";
	weightMap[ *ei ] = weight[ count++ ];
	// cerr << " i = " << *ei << " weight = " << weightMap[ *ei ] << endl;
    }
}


//------------------------------------------------------------------------------
//	
void minimumBlending( Polyhedron & poly, Graph & dual )
{
    unsigned int nEdges = num_edges( dual );
    vector< double > weight( nEdges );
    vector< double > length( nEdges );
    double minLength = 1000000.0, maxLength = -1000000.0;
    vector< double > prod( nEdges );
    static Vector3 orient( 0.0, 0.0, 1.0 );
    const double weightP = 0.5, weightF = 0.5;
	
    // Initialization of edge connections
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	hi->connect() = false;
    }

    // Perimeter heuristic weights
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	Vector3 vec = hi->vertex()->point() - hi->opposite()->vertex()->point();
	length[ eID ] = sqrt( vec.squared_length() );
	if ( length[ eID ] < minLength ) minLength = length[ eID ];
	if ( length[ eID ] > maxLength ) maxLength = length[ eID ];
	prod[ eID ] = fabs( orient * vec )/sqrt( vec.squared_length() );
    }
    cerr << " minLength = " << minLength << " maxLength = " << maxLength << endl;

    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	assert( length[ eID ] > 0.0 );
	weight[ eID ] = 
	    weightP * ( maxLength - length[ eID ] ) / ( maxLength - minLength ) +
	    weightF * prod[ eID ];
	weight[ eID ] = 2.0 * weight[ eID ] - 1.0;
	// Special handling
	if ( hi->fixed() ) weight[ eID ] = hi->weight();
	else hi->weight() = hi->opposite()->weight() = weight[ eID ];
    }

    EdgeIDMap weightMap = get( edge_weight, dual );
	
    boost::graph_traits< Graph >::edge_iterator ei, ei_end;
    int count = 0;
    for( tie( ei, ei_end ) = edges( dual ); ei != ei_end; ++ei ) {
	//std::cout << *ei << ", ";
	weightMap[ *ei ] = weight[ count++ ];
	// cerr << " i = " << *ei << " weight = " << weightMap[ *ei ] << endl;
    }
}


void maximumBlending( Polyhedron & poly, Graph & dual )
{
    unsigned int nEdges = num_edges( dual );
    vector< double > weight( nEdges );
    vector< double > length( nEdges );
    double minLength = 1000000.0, maxLength = -1000000.0;
    vector< double > prod( nEdges );
    static Vector3 orient( 0.0, 0.0, 1.0 );
    const double weightP = 0.5, weightF = 0.5;
	
    // Initialization of edge connections
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	hi->connect() = false;
    }

    // Perimeter heuristic weights
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	Vector3 vec = hi->vertex()->point() - hi->opposite()->vertex()->point();
	length[ eID ] = sqrt( vec.squared_length() );
	if ( length[ eID ] < minLength ) minLength = length[ eID ];
	if ( length[ eID ] > maxLength ) maxLength = length[ eID ];
	prod[ eID ] = fabs( orient * vec )/sqrt( vec.squared_length() );
    }
    cerr << " minLength = " << minLength << " maxLength = " << maxLength << endl;

    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	assert( length[ eID ] > 0.0 );
	weight[ eID ] = 
	    weightP * ( length[ eID ] - minLength ) / ( maxLength - minLength ) +
	    weightF * prod[ eID ];
	weight[ eID ] = 2.0 * weight[ eID ] - 1.0;
	// Special handling
	if ( hi->fixed() ) weight[ eID ] = hi->weight();
	else hi->weight() = hi->opposite()->weight() = weight[ eID ];
    }

    EdgeIDMap weightMap = get( edge_weight, dual );
	
    boost::graph_traits< Graph >::edge_iterator ei, ei_end;
    int count = 0;
    for( tie( ei, ei_end ) = edges( dual ); ei != ei_end; ++ei ) {
	//std::cout << *ei << ", ";
	weightMap[ *ei ] = weight[ count++ ];
	// cerr << " i = " << *ei << " weight = " << weightMap[ *ei ] << endl;
    }
}


//------------------------------------------------------------------------------
//	
void minimumCurvature( Polyhedron & poly, Graph & dual )
{
    unsigned int nEdges = num_edges( dual );
    unsigned int nVertices = num_vertices( dual );
    // unsigned int nEdges = poly.size_of_halfedges() / 2;
    // unsigned int nVertices = poly.size_of_vertices();
    vector< double > angle( nVertices );
    vector< double > rep( nEdges );
    vector< double > weight( nEdges );
    double minRep = 100000.0, maxRep = 0.0;
	
    // Initialization of edge connections
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	hi->connect() = false;
    }

    // Curvature heuristic weights
    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	unsigned int vID = vi->id();
	assert( vID < nVertices );
	angle[ vID ] = centralAngle( vi );
    }
    
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	assert( eID < nEdges );
	unsigned int vIDA = hi->vertex()->id();
	unsigned int vIDB = hi->opposite()->vertex()->id();
	assert( vIDA < nVertices );
	assert( vIDB < nVertices );
	rep[ eID ] = 0.5 * ( angle[ vIDA ] + angle[ vIDB ] );
	// rep[ eID ] = MAX2( angle[ vIDA ], angle[ vIDB ] );
	if ( rep[ eID ] < minRep ) minRep = rep[ eID ];
	if ( rep[ eID ] > maxRep ) maxRep = rep[ eID ];
    }
    cerr << " minRep = " << minRep << " maxRep = " << maxRep << endl;
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	assert( rep[ eID ] >= 0.0 );
	assert( eID < nEdges );
	weight[ eID ] = 2.0 * ( rep[ eID ] - minRep ) / ( maxRep - minRep ) - 1.0;
	// Special handling
	if ( hi->fixed() ) weight[ eID ] = hi->weight();
	else hi->weight() = hi->opposite()->weight() = weight[ eID ];
    }

    EdgeIDMap weightMap = get( edge_weight, dual );
	
    boost::graph_traits< Graph >::edge_iterator ei, ei_end;
    int count = 0;
    for( tie( ei, ei_end ) = edges( dual ); ei != ei_end; ++ei ) {
	weightMap[ *ei ] = weight[ count++ ];
    }

    angle.clear();
    rep.clear();
    weight.clear();
}


void maximumCurvature( Polyhedron & poly, Graph & dual )
{
    unsigned int nEdges = num_edges( dual );
    unsigned int nVertices = num_vertices( dual );
    // unsigned int nEdges = poly.size_of_halfedges() / 2;
    // unsigned int nVertices = poly.size_of_vertices();
    vector< double > angle( nVertices );
    vector< double > rep( nEdges );
    vector< double > weight( nEdges );
    double minRep = 100000.0, maxRep = 0.0;
	
    // Initialization of edge connections
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	hi->connect() = false;
    }

    // Curvature heuristic weights
    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	unsigned int vID = vi->id();
	assert( vID < nVertices );
	angle[ vID ] = centralAngle( vi );
    }
    
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	assert( eID < nEdges );
	unsigned int vIDA = hi->vertex()->id();
	unsigned int vIDB = hi->opposite()->vertex()->id();
	assert( vIDA < nVertices );
	assert( vIDB < nVertices );
	// rep[ eID ] = 0.5 * ( angle[ vIDA ] + angle[ vIDB ] );
	rep[ eID ] = MAX2( angle[ vIDA ], angle[ vIDB ] );
	if ( rep[ eID ] < minRep ) minRep = rep[ eID ];
	if ( rep[ eID ] > maxRep ) maxRep = rep[ eID ];
    }
    cerr << " minRep = " << minRep << " maxRep = " << maxRep << endl;
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	assert( rep[ eID ] >= 0.0 );
	assert( eID < nEdges );
	weight[ eID ] = 2.0 * ( maxRep - rep[ eID ] ) / ( maxRep - minRep ) - 1.0;
	// Special handling
	if ( hi->fixed() ) weight[ eID ] = hi->weight();
	else hi->weight() = hi->opposite()->weight() = weight[ eID ];
    }

    EdgeIDMap weightMap = get( edge_weight, dual );
	
    boost::graph_traits< Graph >::edge_iterator ei, ei_end;
    int count = 0;
    for( tie( ei, ei_end ) = edges( dual ); ei != ei_end; ++ei ) {
	weightMap[ *ei ] = weight[ count++ ];
    }

    angle.clear();
    rep.clear();
    weight.clear();
}


//------------------------------------------------------------------------------
//	
void minimumConcavity( Polyhedron & poly, Graph & dual )
{
    unsigned int nEdges = num_edges( dual );
    // unsigned int nVertices = num_vertices( dual );
    // unsigned int nEdges = poly.size_of_halfedges() / 2;
    // unsigned int nVertices = poly.size_of_vertices();
    vector< double > weight( nEdges );
	
    // Initialization of edge connections
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	hi->connect() = false;
    }

    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	assert( eID < nEdges );
	weight[ eID ] = evalConvexity( hi );
    }

    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	assert( eID < nEdges );
	// Special handling
	if ( hi->fixed() ) weight[ eID ] = hi->weight();
	else hi->weight() = hi->opposite()->weight() = weight[ eID ];
    }

    EdgeIDMap weightMap = get( edge_weight, dual );
	
    boost::graph_traits< Graph >::edge_iterator ei, ei_end;
    int count = 0;
    for( tie( ei, ei_end ) = edges( dual ); ei != ei_end; ++ei ) {
	weightMap[ *ei ] = weight[ count++ ];
    }

    weight.clear();
}


void maximumConcavity( Polyhedron & poly, Graph & dual )
{
    unsigned int nEdges = num_edges( dual );
    // unsigned int nVertices = num_vertices( dual );
    // unsigned int nEdges = poly.size_of_halfedges() / 2;
    // unsigned int nVertices = poly.size_of_vertices();
    vector< double > weight( nEdges );
	
    // Initialization of edge connections
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	hi->connect() = false;
    }

    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	assert( eID < nEdges );
	weight[ eID ] = evalConvexity( hi );
    }

    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	assert( eID < nEdges );
	// Special handling
	if ( hi->fixed() ) weight[ eID ] = hi->weight();
	else hi->weight() = hi->opposite()->weight() = weight[ eID ];
    }

    EdgeIDMap weightMap = get( edge_weight, dual );
	
    boost::graph_traits< Graph >::edge_iterator ei, ei_end;
    int count = 0;
    for( tie( ei, ei_end ) = edges( dual ); ei != ei_end; ++ei ) {
	weightMap[ *ei ] = weight[ count++ ];
    }

    weight.clear();
}


//------------------------------------------------------------------------------
//	
void minimumHyperbolicAve( Polyhedron & poly, Graph & dual )
{
    unsigned int nEdges = num_edges( dual );
    unsigned int nVertices = num_vertices( dual );
    // unsigned int nEdges = poly.size_of_halfedges() / 2;
    // unsigned int nVertices = poly.size_of_vertices();
    vector< double > angle( nVertices );
    vector< double > rep( nEdges );
    vector< double > weight( nEdges );
    double minRep = 100000.0, maxRep = 0.0;
	
    // Initialization of edge connections
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	hi->connect() = false;
    }

    // Curvature heuristic weights
    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	unsigned int vID = vi->id();
	assert( vID < nVertices );
	angle[ vID ] = centralAngle( vi );
    }
    
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	assert( eID < nEdges );
	unsigned int vIDA = hi->vertex()->id();
	unsigned int vIDB = hi->opposite()->vertex()->id();
	assert( vIDA < nVertices );
	assert( vIDB < nVertices );
	rep[ eID ] = 0.5 * ( angle[ vIDA ] + angle[ vIDB ] );
	// rep[ eID ] = MAX2( angle[ vIDA ], angle[ vIDB ] );
	if ( rep[ eID ] < minRep ) minRep = rep[ eID ];
	if ( rep[ eID ] > maxRep ) maxRep = rep[ eID ];
    }
    cerr << " minRep = " << minRep << " maxRep = " << maxRep << endl;
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	assert( rep[ eID ] >= 0.0 );
	assert( eID < nEdges );
	weight[ eID ] = 2.0 * ( rep[ eID ] - minRep ) / ( maxRep - minRep ) - 1.0;
	// Special handling
	if ( hi->fixed() ) weight[ eID ] = hi->weight();
	else hi->weight() = hi->opposite()->weight() = weight[ eID ];
    }

    EdgeIDMap weightMap = get( edge_weight, dual );
	
    boost::graph_traits< Graph >::edge_iterator ei, ei_end;
    int count = 0;
    for( tie( ei, ei_end ) = edges( dual ); ei != ei_end; ++ei ) {
	weightMap[ *ei ] = weight[ count++ ];
    }

    angle.clear();
    rep.clear();
    weight.clear();
}


void minimumHyperbolicMax( Polyhedron & poly, Graph & dual )
{
    unsigned int nEdges = num_edges( dual );
    unsigned int nVertices = num_vertices( dual );
    // unsigned int nEdges = poly.size_of_halfedges() / 2;
    // unsigned int nVertices = poly.size_of_vertices();
    vector< double > angle( nVertices );
    vector< double > rep( nEdges );
    vector< double > weight( nEdges );
    double minRep = 100000.0, maxRep = 0.0;
	
    // Initialization of edge connections
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	hi->connect() = false;
    }

    // Curvature heuristic weights
    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	unsigned int vID = vi->id();
	assert( vID < nVertices );
	angle[ vID ] = centralAngle( vi );
    }
    
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	assert( eID < nEdges );
	unsigned int vIDA = hi->vertex()->id();
	unsigned int vIDB = hi->opposite()->vertex()->id();
	assert( vIDA < nVertices );
	assert( vIDB < nVertices );
	// rep[ eID ] = 0.5 * ( angle[ vIDA ] + angle[ vIDB ] );
	rep[ eID ] = MAX2( angle[ vIDA ], angle[ vIDB ] );
	if ( rep[ eID ] < minRep ) minRep = rep[ eID ];
	if ( rep[ eID ] > maxRep ) maxRep = rep[ eID ];
    }
    cerr << " minRep = " << minRep << " maxRep = " << maxRep << endl;
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	assert( rep[ eID ] >= 0.0 );
	assert( eID < nEdges );
	weight[ eID ] = 2.0 * ( rep[ eID ] - minRep ) / ( maxRep - minRep ) - 1.0;
	// Special handling
	if ( hi->fixed() ) weight[ eID ] = hi->weight();
	else hi->weight() = hi->opposite()->weight() = weight[ eID ];
    }

    EdgeIDMap weightMap = get( edge_weight, dual );
	
    boost::graph_traits< Graph >::edge_iterator ei, ei_end;
    int count = 0;
    for( tie( ei, ei_end ) = edges( dual ); ei != ei_end; ++ei ) {
	weightMap[ *ei ] = weight[ count++ ];
    }

    angle.clear();
    rep.clear();
    weight.clear();
}


void minimumHyperbolicMin( Polyhedron & poly, Graph & dual )
{
    unsigned int nEdges = num_edges( dual );
    unsigned int nVertices = num_vertices( dual );
    // unsigned int nEdges = poly.size_of_halfedges() / 2;
    // unsigned int nVertices = poly.size_of_vertices();
    vector< double > angle( nVertices );
    vector< double > rep( nEdges );
    vector< double > weight( nEdges );
    double minRep = 100000.0, maxRep = 0.0;
	
    // Initialization of edge connections
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	hi->connect() = false;
    }

    // Curvature heuristic weights
    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	unsigned int vID = vi->id();
	assert( vID < nVertices );
	angle[ vID ] = centralAngle( vi );
    }
    
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	assert( eID < nEdges );
	unsigned int vIDA = hi->vertex()->id();
	unsigned int vIDB = hi->opposite()->vertex()->id();
	assert( vIDA < nVertices );
	assert( vIDB < nVertices );
	// rep[ eID ] = 0.5 * ( angle[ vIDA ] + angle[ vIDB ] );
	rep[ eID ] = MIN2( angle[ vIDA ], angle[ vIDB ] );
	if ( rep[ eID ] < minRep ) minRep = rep[ eID ];
	if ( rep[ eID ] > maxRep ) maxRep = rep[ eID ];
    }
    cerr << " minRep = " << minRep << " maxRep = " << maxRep << endl;
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	assert( rep[ eID ] >= 0.0 );
	assert( eID < nEdges );
	weight[ eID ] = 2.0 * ( rep[ eID ] - minRep ) / ( maxRep - minRep ) - 1.0;
	// Special handling
	if ( hi->fixed() ) weight[ eID ] = hi->weight();
	else hi->weight() = hi->opposite()->weight() = weight[ eID ];
    }

    EdgeIDMap weightMap = get( edge_weight, dual );
	
    boost::graph_traits< Graph >::edge_iterator ei, ei_end;
    int count = 0;
    for( tie( ei, ei_end ) = edges( dual ); ei != ei_end; ++ei ) {
	weightMap[ *ei ] = weight[ count++ ];
    }

    angle.clear();
    rep.clear();
    weight.clear();
}

//------------------------------------------------------------------------------
//	
Segment3 closestBone( Point3 & key )
{
    double minDist = 1.0e+10;
    Segment3 closest;
    int index = NO_INDEX;

    for ( unsigned k = 0; k < refer.size(); ++k ) {
	Vector3 offset = refer[ k ] - key;
	double curDist = offset.squared_length();

	if ( curDist < minDist ) {
	    minDist = curDist;
	    closest = bone[ k ];
	    index = k;
	}
    }

    cerr << "Closest index = " << index << endl;
    return closest;
}


void skeletonFlat( Polyhedron & poly, Graph & dual )
{
    unsigned int nEdges = num_edges( dual );
    vector< double > weight( nEdges );
    static Vector3 orient( 0.0, 0.0, 1.0 );
	
    for ( unsigned k = 0; k < bone.size(); ++k ) {
	Vector3 toSrc = bone[ k ].source() - CGAL::ORIGIN;
	Vector3 toTar = bone[ k ].target() - CGAL::ORIGIN;
	Point3 mid = ( CGAL::ORIGIN + 0.5 * toSrc + 0.5 * toTar );
	refer.push_back( mid );
    }

    // Initialization of edge connections
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	hi->connect() = false;
    }

    // flat spanning heuristic weights
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	Segment3 close = closestBone( hi->mid() );
	Vector3 orient = close.to_vector();
	orient = orient / sqrt( orient.squared_length() );
	Vector3 vec = hi->vertex()->point() - hi->opposite()->vertex()->point();
	weight[ eID ] = fabs( orient * vec )/sqrt( vec.squared_length() );
	weight[ eID ] = 2.0 * weight[ eID ] - 1.0;
	// Special handling
	if ( hi->fixed() ) weight[ eID ] = hi->weight();
	else hi->weight() = hi->opposite()->weight() = weight[ eID ];
    }

    for ( unsigned int i = 0; i < nEdges; ++i ) {
	if ( ( 0.0 > weight[ i ] ) && ( weight[ i ] > 1.0 ) ) {
	    cerr << " i = " << i << " weight = " << weight[ i ] << endl;
	    assert( ( 0.0 <= weight[ i ] ) && ( weight[ i ] <= 1.0 ) );
	}
    }

    EdgeIDMap weightMap = get( edge_weight, dual );
	
    boost::graph_traits< Graph >::edge_iterator ei, ei_end;
    int count = 0;
    for( tie( ei, ei_end ) = edges( dual ); ei != ei_end; ++ei ) {
	//std::cout << *ei << ", ";
	weightMap[ *ei ] = weight[ count++ ];
	// cerr << " i = " << *ei << " weight = " << weightMap[ *ei ] << endl;
    }
}


void skeletonCoil( Polyhedron & poly, Graph & dual )
{
    unsigned int nEdges = num_edges( dual );
    vector< double > weight( nEdges );
    static Vector3 orient( 0.0, 0.0, 1.0 );
	
    for ( unsigned k = 0; k < bone.size(); ++k ) {
	Vector3 toSrc = bone[ k ].source() - CGAL::ORIGIN;
	Vector3 toTar = bone[ k ].target() - CGAL::ORIGIN;
	Point3 mid = ( CGAL::ORIGIN + 0.5 * toSrc + 0.5 * toTar );
	refer.push_back( mid );
    }

    // Initialization of edge connections
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	hi->connect() = false;
    }

    // flat spanning heuristic weights
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int eID = hi->id();
	Segment3 close = closestBone( hi->mid() );
	Vector3 orient = close.to_vector();
	orient = orient / sqrt( orient.squared_length() );
	Vector3 vec = hi->vertex()->point() - hi->opposite()->vertex()->point();
	weight[ eID ] = 1.0 - fabs( orient * vec )/sqrt( vec.squared_length() );
	weight[ eID ] = 2.0 * weight[ eID ] - 1.0;
	// Special handling
	if ( hi->fixed() ) weight[ eID ] = hi->weight();
	else hi->weight() = hi->opposite()->weight() = weight[ eID ];
    }

    for ( unsigned int i = 0; i < nEdges; ++i ) {
	if ( ( 0.0 > weight[ i ] ) && ( weight[ i ] > 1.0 ) ) {
	    // cerr << " i = " << i << " weight = " << weight[ i ] << endl;
	    assert( ( 0.0 <= weight[ i ] ) && ( weight[ i ] <= 1.0 ) );
	}
    }

    EdgeIDMap weightMap = get( edge_weight, dual );
	
    boost::graph_traits< Graph >::edge_iterator ei, ei_end;
    int count = 0;
    for( tie( ei, ei_end ) = edges( dual ); ei != ei_end; ++ei ) {
	//std::cout << *ei << ", ";
	weightMap[ *ei ] = weight[ count++ ];
	// cerr << " i = " << *ei << " weight = " << weightMap[ *ei ] << endl;
    }
}


void assignWeights( Polyhedron & poly, Graph & dual, Weight_type sw )
{
    switch ( sw ) {
      case UNIFORM_WEIGHT:
	  uniformWeights	( poly, dual );
	  break;
      case RANDOM_WEIGHT:
	  randomWeights		( poly, dual );
	  break;
      case MINIMUM_PERIMETER:
	  minimumPerimeter	( poly, dual );
	  break;
      case MAXIMUM_PERIMETER:
	  maximumPerimeter	( poly, dual );
	  break;
      case FLAT_SPANNING:
	  flatSpanning		( poly, dual );
	  break;
      case COIL_SPANNING:
	  coilSpanning		( poly, dual );
	  break;
      case MINIMUM_DIHEDRAL:
	  minimumDihedral	( poly, dual );
	  break;
      case MAXIMUM_DIHEDRAL:
	  maximumDihedral	( poly, dual );
	  break;
      case MINIMUM_BLENDING:
	  minimumBlending	( poly, dual );
	  break;
      case MAXIMUM_BLENDING:
	  maximumBlending	( poly, dual );
	  break;
      case MINIMUM_CURVATURE:
	  minimumCurvature	( poly, dual );
	  break;
      case MAXIMUM_CURVATURE:
	  // maximumCurvature	( poly, dual );
	  break;
      case MINIMUM_CONCAVITY:
	  minimumConcavity	( poly, dual );
	  break;
      case MAXIMUM_CONCAVITY:
	  maximumConcavity	( poly, dual );
	  break;
      case MINIMUM_HYPERBOLIC_AVE:
	  minimumHyperbolicAve	( poly, dual );
	  break;
      case MINIMUM_HYPERBOLIC_MAX:
	  minimumHyperbolicMax	( poly, dual );
	  break;
      case MINIMUM_HYPERBOLIC_MIN:
	  minimumHyperbolicMin	( poly, dual );
	  break;
      case SKELETON_FLAT:
	  skeletonFlat		( poly, dual );
	  break;
      case SKELETON_COIL:
	  skeletonCoil		( poly, dual );
	  break;
      default:
	  cerr << "Illegal weight type!!" << endl;
	  break;
    }
}

//------------------------------------------------------------------------------
//	Check the statistics of edge weights
//------------------------------------------------------------------------------
// Calculate the average and standard deviation of the edge weights
void weightStatistics( Polyhedron & poly, double & ave, double & stdev )
{
    double sumWeights = 0.0;
    double sumVariance = 0.0;

    // Average
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	if ( hi->opposite()->vertex()->id() < hi->vertex()->id() ) {
	    sumWeights += hi->weight();
	}
    }
    ave = 2.0 * sumWeights / ( double )poly.size_of_halfedges();
    // Standard deviation
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	if ( hi->opposite()->vertex()->id() < hi->vertex()->id() ) {
	    sumVariance += SQUARE( ave - hi->weight() );
	}
    }
    stdev = sqrt( 2.0 * sumVariance / ( double )poly.size_of_halfedges() ); 

    cerr << " average = " << ave << " standard deviation = " << stdev << endl;

    return;
}


//------------------------------------------------------------------------------
//	compute total area of the 3D mesh
//------------------------------------------------------------------------------
double totalArea( Polyhedron & poly )
{
    double sum = 0.0;
    for ( Facet_iterator fi = poly.facets_begin(); fi != poly.facets_end(); ++fi ) {
	Facet_handle fh = fi;
	Halfedge_handle hhR = fh->halfedge();
	Halfedge_handle hhL = hhR->prev()->opposite();
	Vector3 vecR = hhR->vertex()->point() - hhR->opposite()->vertex()->point();
	Vector3 vecL = hhL->vertex()->point() - hhL->opposite()->vertex()->point();
	Vector3 oprd = cross_product( vecL, vecR );
	sum += 0.5 * sqrt( oprd.squared_length() );
	// cerr << " fh->id = " << fh->id() << " sum = " << sum << endl;
    }

    cerr << " total area = " << sum << endl;
    return sum;
}

//------------------------------------------------------------------------------
//	To verify the compactness of features (peak and saddle) over its' incidence edges
//  The higher the density indicate more saddle points or peak points surrounding it
//------------------------------------------------------------------------------
double densityAtFeature( Vertex_handle & vh, int feature )
{
    assert( vh->label() == feature );
    Halfedge_vertex_circulator hvc;
    // First scan to count the number of cut edges
    hvc = vh->vertex_begin();
    int nFeatures = 0, nIncidents = 0;
    do {
	Vertex_handle vhC = hvc->opposite()->vertex();
	nIncidents++; // edges that incident to the "feature" vertex
	if ( vhC->label() == feature ) {
	    nFeatures++; // check for the neighboring "feature" vertex
	}
    } while ( ++hvc != vh->vertex_begin() );

//	cout << "feature: " << feature << endl;
//	cout << "nFeatures: " << nFeatures << endl;
//	cout << "nIncidents: " << nIncidents << endl;
	
    return ( double )nFeatures/( double )nIncidents;
}


//------------------------------------------------------------------------------
//	To check whether the given dual edge can be perfect matching or vice versa
//------------------------------------------------------------------------------
bool fixDual( Graph & dual, Halfedge_handle & hh )
{
    const int nD = 4;
    Halfedge_handle curH = hh;
    // int curS = curH->label();

    // diamond[0]:DNext, diamond[1]:DPrev, diamond[2]:LNext, diamond[3]:LPrev
    Halfedge_handle diaH[ 4 ];
    int diaS[ 4 ], srcVID[ 4 ], tarVID[ 4 ], srcDeg[ 4 ], tarDeg[ 4 ];
    // Check if the edge can be chosen as a cut edge
    diaH[ 0 ] = curH->opposite()->prev();
    diaH[ 1 ] = curH->opposite()->next();
    diaH[ 2 ] = curH->next();
    diaH[ 3 ] = curH->prev();
    for ( int i = 0; i < nD; ++i ) {
#ifdef DEBUG
	if ( diaH[ i ]->label() != diaH[ i ]->opposite()->label() ) {
	    cerr << " halfedge IDA = " << diaH[ i ]->id() 
		 << " halfedge IDB = " << diaH[ i ]->opposite()->id() << endl;
	    cerr << " halfedge labelA = " << diaH[ i ]->label() 
		 << " halfedge labelB = " << diaH[ i ]->opposite()->label() << endl;
	}
#endif	// DEBUG
	diaS[ i ]  = diaH[ i ]->label();
	srcVID[ i ] = diaH[ i ]->facet()->id();
	tarVID[ i ] = diaH[ i ]->opposite()->facet()->id();
    }
    if ( ( curH->label() != FIXED_EDGE ) && ( curH->label() != INVALID_EDGE ) &&
	 ( diaS[ 0 ] != FIXED_EDGE ) && ( diaS[ 1 ] != FIXED_EDGE ) &&
	 ( diaS[ 2 ] != FIXED_EDGE ) && ( diaS[ 3 ] != FIXED_EDGE ) ) {
	for ( int i = 0; i < nD; ++i )
	    if ( diaS[ i ] != INVALID_EDGE ) {
		remove_edge( srcVID[ i ], tarVID[ i ], dual );
#ifdef MY_DEBUG
		cerr << " remove edge : " << srcVID[ i ] << " == " << tarVID[ i ] << endl;
#endif	// MY_DEBUG
	    }
	
	bool success = true;
	
	for ( int i = 0; i < nD; ++i ) {
	    srcDeg[ i ] = out_degree( srcVID[ i ], dual );
	    tarDeg[ i ] = out_degree( tarVID[ i ], dual );
	    if ( ( srcDeg[ i ] <= 0 ) || ( tarDeg[ i ] <= 0 ) ) success = false;
	}
	    
	if ( success ) success = isPerfect( dual );
	
	if ( success ) {
	    curH->label() = curH->opposite()->label() = FIXED_EDGE;
	    for ( int i = 0; i < nD; ++i ) {
#ifdef DEBUG
		cerr << " Dual edge : " << diaH[ i ]->facet()->id() << " == "
		     << diaH[ i ]->opposite()->facet()->id() << endl;
#endif	// DEBUG
		if ( diaS[ i ] != INVALID_EDGE ) {
		    diaH[ i ]->label() = diaH[ i ]->opposite()->label() = INVALID_EDGE;
#ifdef DEBUG
		    cerr << " Invalidate the edge " << endl;
#endif	// DEBUG
		}
	    }
	    return true;
	}
	else {
	    for ( int i = 0; i < nD; ++i ) {
		if ( diaS[ i ] != INVALID_EDGE ) {
		    add_edge( srcVID[ i ], tarVID[ i ], diaH[ i ]->id(), dual );
		}
	    }
	    return false;
	}
    }
    else return false;
}


//------------------------------------------------------------------------------
//	To check whether the given dual edge can be perfect matching or vice versa
//------------------------------------------------------------------------------
bool fixWeight( Graph & dual, Halfedge_handle & hh, double value )
{
    hh->weight() = hh->opposite()->weight()	= value;
    hh->fixed()	= hh->opposite()->fixed()	= true;

    EdgeIDMap weightMap = get( edge_weight, dual );
    boost::graph_traits< Graph >::edge_iterator ei, ei_end;
    int count = 0;
    for( tie( ei, ei_end ) = edges( dual ); ei != ei_end; ++ei ) {
	if ( count == hh->id() ) {
	    cerr << " Now weight value for the edge ID " << count << " is set from "
		 << weightMap[ *ei ] << " to " << hh->weight() << endl;
	    weightMap[ *ei ] = hh->weight();
	}
	count++;
    }

    return true;
}

//------------------------------------------------------------------------------
//	Validate the type of vertices: saddle vertices, peak vertices and default vertices
//------------------------------------------------------------------------------
void labelByCurvature( Polyhedron & poly )
{
    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	vi->label() = DEFAULT_LABEL;
	if ( isSaddle( vi ) ) {
#ifdef DEBUG
	    cerr << " This is a saddle point : " << vi->id() << endl;
#endif	// DEBUG
	    vi->label() = SADDLE_VERTEX;
	}
	else if ( isPeak( vi ) ) {
	    vi->label() = PEAK_VERTEX;
	}
	else {
	    vi->label() = DEFAULT_LABEL;
	}
    }
}


#ifdef NO_NEED
int randomSeek( Graph & dual, stack< Vertex_handle > & saddle,
		int nEdges, 
		const Halfedge_vertex_circulator & cur, 
		const Halfedge_vertex_circulator & last )
{
    const int nD = 4;

    static int count = 0;
    // const int stopCount = 1;
    // const int stopCount = 7;
    // const int stopCount = 8;
    // const int stopCount = 78;
    // const int stopCount = 129; // complete
    // const int stopCount = 10000000; // complete

    Vertex_handle vhC;
    Halfedge_vertex_circulator hvcC, hvcS, hvcE;
    int nCuts;

    if ( nEdges == 0 ) {
	if ( saddle.empty() ) {
	    cerr << "Constrained perfect matching succeeded at the iteration count " 
		 << count << endl;
	    return TRUE;
	}
	else {
	    vhC = saddle.top();
	    int check = checkSaddle( vhC );
	    if ( check == -1 ) return FALSE;
	    else {
		nCuts = 0;
		hvcC = vhC->vertex_begin();
		do {
		    if ( hvcC->label() == FIXED_EDGE ) nCuts++;
		} while ( ++hvcC != vhC->vertex_begin() );
		hvcS = vhC->vertex_begin();
		hvcE = vhC->vertex_begin();
	    }
	}
    }
    // Do not attach "else" here
    if ( nEdges == 1 ) {
	nCuts = 1;
	vhC = saddle.top();
	hvcS = cur;
	hvcE = last;
    }
    // Do not attach "else" here
    if ( nEdges >= 2 ) {
	Vertex_handle store = saddle.top();

	int isCompatible = isSurrounds( store );
	if ( !isCompatible ) return FALSE;

	saddle.pop();
	// We do not care about hvcS & hvcE nwhen nEdges == 0
	int check = randomSeek( dual, saddle, 0, hvcS, hvcE );
	if ( check ) return TRUE;
	else {
	    saddle.push( store );
	    return FALSE;
	}

    }

//------------------------------------------------------------------------------
//	
//------------------------------------------------------------------------------
    // Pick up one of the remaing edges for a cut edge
    hvcC = hvcS;
    int nInc = 0;
    do {
#ifdef MY_DEBUG
	cerr << " Edge ID [ " << nInc << " ] => " << hvcC->id() << endl;
#endif	// MY_DEBUG
	Halfedge_handle curH = hvcC;
	int curS = curH->label();

#ifdef MYDEBUG
	cerr << " At Vertex No. " << vh->id() << " Handle Edge No. " << curH->id() << ":  nCuts = " << nCuts << endl;
#endif	// MYDEBUG

	// diamond[0]:DNext, diamond[1]:OPrev, diamond[2]:LNext, diamond[3]:LPrev
	Halfedge_handle diaH[ 4 ];
	int diaS[ 4 ], srcVID[ 4 ], tarVID[ 4 ], srcDeg[ 4 ], tarDeg[ 4 ];
	// Check if the edge can be chosen as a cut edge
	diaH[ 0 ] = curH->opposite()->prev();
	diaH[ 1 ] = curH->opposite()->next();
	diaH[ 2 ] = curH->next();
	diaH[ 3 ] = curH->prev();
	for ( int i = 0; i < nD; ++i ) {
	    diaS[ i ]  = diaH[ i ]->label();
	    srcVID[ i ] = diaH[ i ]->facet()->id();
	    tarVID[ i ] = diaH[ i ]->opposite()->facet()->id();
	}
	if ( ( curH->label() != FIXED_EDGE ) && ( curH->label() != INVALID_EDGE ) &&
	     ( diaS[ 0 ] != FIXED_EDGE ) && ( diaS[ 1 ] != FIXED_EDGE ) &&
	     ( diaS[ 2 ] != FIXED_EDGE ) && ( diaS[ 3 ] != FIXED_EDGE ) ) {
	    for ( int i = 0; i < nD; ++i )
		if ( diaS[ i ] != INVALID_EDGE ) {
		    remove_edge( srcVID[ i ], tarVID[ i ], dual );
#ifdef MY_DEBUG
		    cerr << " remove edge : " << srcVID[ i ] << " == " << tarVID[ i ] << endl;
#endif	// MY_DEBUG
		}

	    bool success = true;

	    for ( int i = 0; i < nD; ++i ) {
		srcDeg[ i ] = out_degree( srcVID[ i ], dual );
		tarDeg[ i ] = out_degree( tarVID[ i ], dual );
#ifdef MYDEBUG
		cerr << "src[ " << srcVID[ i ] << " ] : "  << srcDeg[i] << endl;
		cerr << "tar[ " << tarVID[ i ] << " ] : "  << tarDeg[i] << endl;
#endif	// MYDEBUG
		if ( ( srcDeg[ i ] <= 0 ) || ( tarDeg[ i ] <= 0 ) ) success = false;
	    }
	    
	    if ( success ) success = isPerfect( dual );

	    if ( success ) {
#ifdef MYDEBUG
		cerr << "split along the edge with two saddles : " << curH->id() << endl;
#endif	// MYDEBUG
		curH->label() = curH->opposite()->label() = FIXED_EDGE;
#ifdef MYDEBUG
		cerr << " Edge ID : " << curH->id() << " becomes FIXED_EDGE" << endl;
#endif	// MYDEBUG
		for ( int i = 0; i < nD; ++i )
		    if ( diaS[ i ] != INVALID_EDGE )
			diaH[ i ]->label() = diaH[ i ]->opposite()->label() = INVALID_EDGE;

		// Post-processing
		nCuts++;

		Vertex_handle store;
		store = saddle.top();
#ifdef DEBUG
		store->debug() = nCuts;
#endif	// DEBUG

		int check = randomSeek( dual, saddle, nCuts, hvcC, hvcE );
		if ( check ) return TRUE;
		else {
		    // here
		    if ( ++count % 10000 == 0 ) cerr << " count = " << count << endl;
		    // if ( ( count++ == stopCount ) ) return TRUE;
#ifdef MYDEBUG
		    cerr << "########## failed in splitting along the edge : ";
		    cerr << curH->vertex()->id() << " == " << curH->opposite()->vertex()->id()
			 << " at Vertex No. " << vh->id() << endl;
#endif	// MYDEBUG
		    // if ( saddle.size() > 17 ) cerr << "SIZE = " << saddle.size() << endl;

#ifdef NO_NEED
		    if ( nCuts == 2 ) saddle.push( store );
#endif	// NO_NEED
		    nCuts--;
		    curH->label() = curH->opposite()->label() = curS;
		    for ( int i = 0; i < nD; ++i ) {
			if ( diaS[ i ] != INVALID_EDGE ) {
			    diaH[ i ]->label() = diaH[ i ]->opposite()->label() = diaS[ i ];
			    add_edge( srcVID[ i ], tarVID[ i ], diaH[ i ]->id(), dual );
#ifdef MY_DEBUG
			    cerr << " add edge : " << srcVID[ i ]
				 << " == " << tarVID[ i ] << endl;
#endif	// MY_DEBUG
			}
		    }
		}
	    }
	    else {
		// here
		if ( ++count % 10000 == 0 ) cerr << " count = " << count << endl;
		// if ( ( count++ == stopCount ) ) return TRUE;
		// curH->label() = curH->opposite()->label() = curS;
		for ( int i = 0; i < nD; ++i ) {
		    if ( diaS[ i ] != INVALID_EDGE ) {
			// diaH[ i ]->label() = diaH[ i ]->opposite()->label() = diaS[ i ];
			add_edge( srcVID[ i ], tarVID[ i ], diaH[ i ]->id(), dual );
#ifdef MY_DEBUG
			cerr << " add edge : " << srcVID[ i ]
			     << " == " << srcVID[ i ] << endl;
#endif	// MY_DEBUG
		    }
		}
	    }
	}
	nInc++;
    } while ( ++hvcC != hvcE );

    return FALSE;
}
#endif	// NO_NEED
    

void sequentialSeek( Graph & dual, stack< Vertex_handle > & set, unsigned char feature, int limit )
{
    const int nD = 4;

    // static int count = 0;
    // const int stopCount = 1;
    // const int stopCount = 7;
    // const int stopCount = 8;
    // const int stopCount = 78;
    // const int stopCount = 129; // complete
    // const int stopCount = 10000000; // complete

    Vertex_handle vhC;
    Halfedge_vertex_circulator hvcC, hvcS, hvcE;
    int nCuts; // number of cutting lines

    while ( ! set.empty() ) {

	vhC = set.top();
	nCuts = 0;
	hvcC = vhC->vertex_begin();
	do {
	    if ( hvcC->label() == FIXED_EDGE ) nCuts++; // perfect matching edges
	} while ( ++hvcC != vhC->vertex_begin() );
	hvcS = vhC->vertex_begin();
	hvcE = vhC->vertex_begin();
	
	// Do not attach "else" here
	if ( nCuts >= limit ) {
	    cerr << " nCuts = " << nCuts << " The saddle already has more than 2 edge cuts" << endl;
	    set.pop();
	    continue;
	}
	
//------------------------------------------------------------------------------
//	Seek in sequence for appropriate cutting edges that incident to saddle feature points
//------------------------------------------------------------------------------
	// Pick up one of the remaining edges for a cut edge
	hvcC = vhC->vertex_begin();
	do {
	    Halfedge_handle curH = hvcC;
	    int curS = curH->label();
	    
#ifdef MYDEBUG
	    cerr << " At Vertex No. " << vh->id() << " Handle Edge No. " << curH->id() << ":  nCuts = " << nCuts << endl;
#endif	// MYDEBUG
	    
	// diamond[0]:DNext, diamond[1]:OPrev, diamond[2]:LNext, diamond[3]:LPrev
	    Halfedge_handle diaH[ 4 ];
	    int diaS[ 4 ], srcVID[ 4 ], tarVID[ 4 ];
	    // int srcDeg[ 4 ], tarDeg[ 4 ];
	    // Check if the edge can be chosen as a cut edge
	    diaH[ 0 ] = curH->opposite()->prev();
	    diaH[ 1 ] = curH->opposite()->next();
	    diaH[ 2 ] = curH->next();
	    diaH[ 3 ] = curH->prev();
	    for ( int i = 0; i < nD; ++i ) {
		diaS[ i ]  = diaH[ i ]->label();
		srcVID[ i ] = diaH[ i ]->facet()->id();
		tarVID[ i ] = diaH[ i ]->opposite()->facet()->id();
	    }
	    if ( ( curH->label() != FIXED_EDGE ) && ( curH->label() != INVALID_EDGE ) &&
		 ( diaS[ 0 ] != FIXED_EDGE ) && ( diaS[ 1 ] != FIXED_EDGE ) &&
		 ( diaS[ 2 ] != FIXED_EDGE ) && ( diaS[ 3 ] != FIXED_EDGE ) ) {
		for ( int i = 0; i < nD; ++i )
		    if ( diaS[ i ] != INVALID_EDGE ) {
			remove_edge( srcVID[ i ], tarVID[ i ], dual );
#ifdef MY_DEBUG
			cerr << " remove edge : " << srcVID[ i ] << " == " << tarVID[ i ] << endl;
#endif	// MY_DEBUG
		    }

		
		bool success = isPerfect( dual );
		
		if ( success ) {
#ifdef MYDEBUG
		    cerr << "split along the edge around the feature point : " << curH->vertex()->id() << endl;
#endif	// MYDEBUG
		    curH->label() = curH->opposite()->label() = FIXED_EDGE;
#ifdef MYDEBUG
		    cerr << " Edge ID : " << curH->id() << " becomes FIXED_EDGE" << endl;
#endif	// MYDEBUG
		    for ( int i = 0; i < nD; ++i )
			if ( diaS[ i ] != INVALID_EDGE )
			    diaH[ i ]->label() = diaH[ i ]->opposite()->label() = INVALID_EDGE;
		    
		    // Post-processing
		    nCuts++;
		}
		else {
		    curH->label() = curH->opposite()->label() = curS;
		    for ( int i = 0; i < nD; ++i ) {
			if ( diaS[ i ] != INVALID_EDGE ) {
			    diaH[ i ]->label() = diaH[ i ]->opposite()->label() = diaS[ i ];
			    add_edge( srcVID[ i ], tarVID[ i ], diaH[ i ]->id(), dual );
#ifdef MY_DEBUG
			    cerr << " add edge : " << srcVID[ i ]
				 << " == " << tarVID[ i ] << endl;
#endif	// MY_DEBUG
			}
		    }
		}
	    }
	} while ( ( ++hvcC != hvcE ) && ( nCuts < limit ) );
	
	cerr << " nSaddles = " << set.size() << " The feature point gets more than " << limit << " edge cuts" << endl;
	set.pop();
    }
}
    

//------------------------------------------------------------------------------
//	Cutting path features
//------------------------------------------------------------------------------
void splitFeatures( Polyhedron & poly, Graph & dual, unsigned char feature, int limit )
{
    vector< bool > visited;

    multimap< double, Vertex_handle > ratio;
    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	visited.push_back( false );
	if ( vi->label() == feature ) {
	    Halfedge_vertex_circulator hvc = vi->vertex_begin();
	    double cur = densityAtFeature( vi, feature ); // the compactness 
//	    ratio.insert( make_pair< double, Vertex_handle >( cur, vi ) );
        ratio.insert( make_pair( cur, vi ) );
	    // cerr << " Vertex ID : " << vi->id() << " Density = " << curRatio << endl;
		// cout << "densityAtFeature: " << cur << endl;
	}
    }
	// ratio: keep track of the features' info (including the size)

    // vector< Vertex_handle > order;
    vector< vector< Vertex_handle > > component;
    multimap< double, Vertex_handle >::iterator it = ratio.end();
    if ( it == ratio.begin() ) return;
    else it--;
    while ( true ) {
	queue < Vertex_handle > candidate;
	vector< Vertex_handle > order;
	while ( ( visited[ it->second->id() ] ) && ( it != ratio.begin() ) ) it--;
	Vertex_handle start;
	if ( it == ratio.begin() ) {
	    if ( visited[ it->second->id() ] == false ) start = it->second;
	    else break;
	}
	else {
	    start = it->second;
	}
	
	visited[ start->id() ] = true;
	
	order.push_back( start );
	candidate.push( start );

	while ( ! candidate.empty() ) {
	    Vertex_handle vh = candidate.front();
	    candidate.pop();
	    Halfedge_vertex_circulator hvc = vh->vertex_begin();
#ifdef MYDEBUG
	    cerr << " Handling Vertex No. " << vh->id() << endl;
#endif	// MYDEBUG
	    do {
		if ( ( hvc->opposite()->vertex()->label() == feature ) &&
		     ( ! visited[ hvc->opposite()->vertex()->id() ] ) ) {
		    visited[ hvc->opposite()->vertex()->id() ] = true;
		    order.push_back( hvc->opposite()->vertex() );
		    candidate.push( hvc->opposite()->vertex() );
		}
	    } while ( ++hvc != vh->vertex_begin() );
#ifdef MYDEBUG
	    cerr << " No. of candidates = " << candidate.size() << endl;
#endif	// MYDEBUG
	}
	    
	component.push_back( order );
    }
	
    cerr << " Number of connected region components = " << component.size() << endl;

    int check = true;
    for ( unsigned int k = 0; k < component.size(); ++k ) {
	cerr << " Component No. " << k 
	     << " number of hyperbolics = " << component[ k ].size() << endl;
	stack< Vertex_handle > set;
	for ( int i = ( int )component[ k ].size() - 1; i >= 0; --i ) {
	    // cerr << " i = " << i << " order = " << order[ i ]->id();
	    set.push( component[ k ][ i ] );
	    // cerr << " siddle size = " << saddle.size() << endl;
	}
	Halfedge_vertex_circulator dummyC, dummyL;
	sequentialSeek( dual, set, feature, limit );
    }
    
    if ( check ) {
	cerr << "Constrained perfect matching succeeded!!" << endl;
	return;
    }
    else assert( FALSE );
}


//------------------------------------------------------------------------------
//	Lay out the cutting paths over the peak points then saddle points
//------------------------------------------------------------------------------
void applyConstraints( Polyhedron & poly, Graph & dual )
{
    splitFeatures( poly, dual, PEAK_VERTEX, PEAK_CUTS );
    splitFeatures( poly, dual, SADDLE_VERTEX, SADDLE_CUTS );
}
