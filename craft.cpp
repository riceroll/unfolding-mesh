//------------------------------------------------------------------------------
//
//	craft.cpp 
//		fundamental operations for papercraft
//
//------------------------------------------------------------------------------
#include <cmath>
#include <cassert>
#include <iostream>
#include <iomanip>
#include "common.h"


using namespace std;

//------------------------------------------------------------------------------
//
//	Fundamental geometric operations
//
//------------------------------------------------------------------------------
Point2 centerTriangle( const Triangle2 & t )
{
    Vector2 a = t[ 0 ] - CGAL::ORIGIN;
    Vector2 b = t[ 1 ] - CGAL::ORIGIN;
    Vector2 c = t[ 2 ] - CGAL::ORIGIN;
    Point2  r = CGAL::ORIGIN + (a+b+c)/3.0;
    return r;
}

bool isIntersected( const Triangle2 & triA, const Triangle2 & triB )
{
    const double shrink = SHRINKAGE_RATIO;

    // 1st stage: rough check with bounding boxes
    Bbox2 boxA = triA.bbox();
    Bbox2 boxB = triB.bbox();
    if ( ! do_overlap( boxA, boxB ) ) return false;

    // 2nd stage: rigorous check
    Vector2 cenA = centerTriangle( triA ) - CGAL::ORIGIN;
    Vector2 cenB = centerTriangle( triB ) - CGAL::ORIGIN;
    for ( unsigned int i = 0; i < NUM_SIDES; ++i ) {
	Point2 origA = CGAL::ORIGIN + shrink * ( ( triA[ i ] - CGAL::ORIGIN ) - cenA ) + cenA;
	Point2 destA = CGAL::ORIGIN + shrink * ( ( triA[i+1] - CGAL::ORIGIN ) - cenA ) + cenA;
	Segment2 segA( origA, destA );
	for ( unsigned int j = 0; j < NUM_SIDES; ++j ) {
	    Point2 origB = CGAL::ORIGIN + shrink * ( ( triB[ j ] - CGAL::ORIGIN ) - cenB ) + cenB;
	    Point2 destB = CGAL::ORIGIN + shrink * ( ( triB[j+1] - CGAL::ORIGIN ) - cenB ) + cenB;
	    Segment2 segB( origB, destB );
	    if ( do_intersect( segA, segB ) ) {
#ifdef MYDEBUG
		cerr << " Intersection with Triangle No. " << k << " was detected." << endl;
		cerr << " triA = " << triA << endl;
		cerr << " triB = " << triB << endl;
		cerr << " cenA = " << cenA << " cenB = " << cenB << endl;
		cerr << " segA " << segA << endl;
		cerr << " segB " << segB << endl;
#endif	// MYDEBUG
		return true;
	    }
	}
    }
    return false;
}


bool isIntersected( const Triangle2 & tri, const vector< Triangle2 > & strip )
{
    const double shrink = SHRINKAGE_RATIO;

    // 1st stage: rough check with bounding boxes
    Bbox2 boxA = tri.bbox();
    Bbox2 boxO;
    for ( unsigned int k = 0; k < strip.size(); ++k )
	boxO = boxO + strip[ k ].bbox();
    if ( ! do_overlap( boxA, boxO ) ) {
	return false;
    }

    for ( unsigned int k = 0; k < strip.size(); ++k ) {
	// Vector2 cenI( 0.0f, 0.0f ), cenJ( 0.0f, 0.0f );
	// for ( unsigned int i = 0; i < NUM_SIDES; ++i ) cenI = cenI + ( strip[ k ][ i ] - CGAL::ORIGIN );
	// cenI = cenI / ( double )NUM_SIDES;
	// for ( unsigned int j = 0; j < NUM_SIDES; ++j ) cenJ = cenJ + ( tri[ j ] - CGAL::ORIGIN );
	// cenJ = cenJ / ( double )NUM_SIDES;
	Vector2 cenI = centerTriangle( strip[ k ] ) - CGAL::ORIGIN;
	Vector2 cenJ = centerTriangle( tri ) - CGAL::ORIGIN;
	for ( unsigned int i = 0; i < NUM_SIDES; ++i ) {

	    // 2nd stage: smaller bounding box 
	    Bbox2 boxB = strip[ k ].bbox();
	    if ( do_overlap( boxA, boxB ) ) {

		// 3rd stage: rigorous check
		Point2 origI = CGAL::ORIGIN + shrink * ( ( strip[ k ][ i ] - CGAL::ORIGIN ) - cenI ) + cenI;
		Point2 destI = CGAL::ORIGIN + shrink * ( ( strip[ k ][i+1] - CGAL::ORIGIN ) - cenI ) + cenI;
		Segment2 segI( origI, destI );
		for ( unsigned int j = 0; j < NUM_SIDES; ++j ) {
		    Point2 origJ = CGAL::ORIGIN + shrink * ( ( tri[ j ] - CGAL::ORIGIN ) - cenJ ) + cenJ;
		    Point2 destJ = CGAL::ORIGIN + shrink * ( ( tri[j+1] - CGAL::ORIGIN ) - cenJ ) + cenJ;
		    Segment2 segJ( origJ, destJ );
		    if ( do_intersect( segI, segJ ) ) {
#ifdef MYDEBUG
			cerr << " Intersection with Triangle No. " << k << " was detected." << endl;
			cerr << " triI = " << strip[k] << endl;
			cerr << " triJ = " << tri << endl;
			cerr << " cenI = " << cenI << " cenJ = " << cenJ << endl;
			cerr << " segI " << segI << endl;
			cerr << " segJ " << segJ << endl;
#endif	// MYDEBUG
			return true;
		    }
		}
	    }
	}
    }
    return false;
}


//------------------------------------------------------------------------------
//
//	Fundamental topological operations
//
//------------------------------------------------------------------------------
Halfedge_handle sharedHalfedge( const Facet_handle & fhI, const Facet_handle & fhO )
{
    if ( fhI->halfedge()->prev()->opposite()->facet()->id() == fhO->id() ) 
	return fhI->halfedge()->prev();
    else if ( fhI->halfedge()->opposite()->facet()->id() == fhO->id() ) 
	return fhI->halfedge();
    else if ( fhI->halfedge()->next()->opposite()->facet()->id() == fhO->id() ) 
	return fhI->halfedge()->next();
    else {
	cerr << " fhO->id() = " << fhO->id() << endl;
	cerr << " fhI->halfedge()->prev()->opposite()->facet()->id() = " 
	     << fhI->halfedge()->prev()->opposite()->facet()->id() << endl;
	cerr << " fhI->halfedge()->opposite()->facet()->id() = " 
	     << fhI->halfedge()->opposite()->facet()->id() << endl;
	cerr << " fhI->halfedge()->next()->opposite()->facet()->id() = " 
	     << fhI->halfedge()->next()->opposite()->facet()->id() << endl;
	cerr << " fhI = " << fhI->halfedge()->prev()->vertex()->id()
	     << " -- " << fhI->halfedge()->vertex()->id()
	     << " -- " << fhI->halfedge()->next()->vertex()->id() << endl;
	cerr << " fhO = " << fhO->halfedge()->prev()->vertex()->id()
	     << " -- " << fhO->halfedge()->vertex()->id()
	     << " -- " << fhO->halfedge()->next()->vertex()->id() << endl;
	assert( false );
    }
}


// Explicity count the number of cut edges at every vertex
void countNCuts( Polyhedron & poly )
{
    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	int numberOfCuts = 0;
	Halfedge_vertex_circulator hvc = vi->vertex_begin();
	do {
	    if ( ! hvc->connect() ) numberOfCuts++;
	} while ( ++hvc != vi->vertex_begin() );
	vi->nCuts()	= numberOfCuts;
	// cerr << "Vertex ID : " << vi->id() << " Number of cut edges = " << vi->nCuts() << endl;
    }
}


void initEdgeAttr( Polyhedron & poly )
{
    // mark the current boundary edges as visited
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	if ( hi->connect() ) hi->visit() = false;
	else hi->visit() = true;
    }
}



Halfedge_handle findPathEnd( Polyhedron & poly )
{
    // Find one endpoint of the vertex spanning tree
    Halfedge_handle initH;
    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	Halfedge_vertex_circulator hvc = vi->vertex_begin();
	Halfedge_handle hh = hvc->opposite();
	int numOfCuts = 0;
	do {
	    if ( ! hh->connect() ) {
		if ( numOfCuts == 0 ) initH = hh;

		numOfCuts++;
	    }
	} while ( ++hvc != vi->vertex_begin() );

	if ( ( ( numOfCuts == 1 ) || ( numOfCuts >= 3 ) ) &&
	     ( initH->path() == NO_INDEX ) ) {
#ifdef MYDEBUG
	    cerr << " Starting vertex : " << vi->id() << endl;
	    cerr << " Starting half edge : " << hh->opposite()->vertex()->id()
		 << " == " << hh->vertex()->id() << endl;
	    getchar();
#endif	// MYDEBUG
	    return initH;
	}
    }
    return NULL;
}


Halfedge_handle findVisited( Polyhedron & poly )
{
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	Halfedge_handle hh = hi;
	if ( ( hh->path() == NO_INDEX ) && ( hh->opposite()->path() != NO_INDEX ) ) {
#ifdef MYDEBUG
	    cerr << " Starting half edge : " << hh->vertex()->id() << " == " << hh->id() << endl;
	    getchar();
#endif	// MYDEBUG
	    return hh;
	}
    }
    return NULL;
}					  


Halfedge_handle nextBoundary( Halfedge_handle & prevH )
{
#ifdef MYDEBUG
    cerr << " previous half edge : " << prevH->opposite()->vertex()->id()
	 << " -- " << prevH->vertex()->id() << endl;
#endif	// MYDEBUG

    Halfedge_handle nextH;
    Vertex_handle vh = prevH->vertex();
    Halfedge_vertex_circulator hvc = vh->vertex_begin();
    do {
#ifdef MYDEBUG
	cerr << " current half edge : " << hvc->opposite()->vertex()->id()
	     << " -- " << hvc->vertex()->id() << endl;
#endif	// MYDEBUG
	if ( ( hvc->vertex() == prevH->vertex() ) &&
	     ( hvc->opposite()->vertex() == prevH->opposite()->vertex() ) )
	    break;
    } while ( ++hvc != vh->vertex_begin() );
    Halfedge_vertex_circulator stop = hvc;
    do { 
	hvc++;
#ifdef MYDEBUG
	cerr << " visited half edge : " << hvc->opposite()->vertex()->id()
	     << " -- " << hvc->vertex()->id() << endl;
#endif	// MYDEBUG
	if ( ! hvc->connect() ) {
#ifdef MYDEBUG
	    cerr << " split " << endl;
#endif	// MYDEBUG
	    return hvc->opposite();
	}
    } while ( hvc != stop );
    return NULL;
}


// Retrieve boundary of the unfolded patch
// The interior faces are assumed at the left side of the half edges in the
// array.
// Returns the number of boundary runs
int extractBoundary	( Polyhedron & poly,  vector< vector< Halfedge_handle > > & bound )
{
    // Initalize boundary path IDs
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	hi->path() = NO_INDEX;
    }
    
    bound.clear();

    bool doStart = true;
    bool isInitial = true;
    int pathID = -1;
    while ( true ) { 
	Halfedge_handle initH;
	if ( doStart ) {
	    initH = findPathEnd( poly );
	    doStart = false;
	}
	else {
	    initH = findVisited( poly );
	}

	if ( initH == NULL ) break;

	vector< Halfedge_handle > curBound;

	curBound.push_back( initH );

	Halfedge_handle curH = initH, prevH;
	// assign path IDs when traversing the boundary edges
	if ( curH->opposite()->path() != NO_INDEX ) {
	    curH->path() = curH->opposite()->path();
	    // curH->orient() = false;
	    curH->orient() = ! curH->opposite()->orient();
	}
	else {
	    pathID++;
	    curH->path() = pathID;
	    if ( isInitial ) {
		curH->orient() = false;
		isInitial = false;
	    }
	    else {
		curH->orient() = true;
	    }
	}

	prevH   = curH;
	curH	= nextBoundary( curH );
	while ( curH != initH ) {

	    curBound.push_back( curH );

	    int nCuts = curH->opposite()->vertex()->nCuts();
	    assert( nCuts != 0 );
	    if ( curH->opposite()->path() != NO_INDEX ) {
		curH->path() = curH->opposite()->path();
		// curH->orient() = false;
		curH->orient() = ! curH->opposite()->orient();
	    }
	    else if ( nCuts == 2 ) {
		curH->path() = pathID;
		// curH->orient() = true;
		curH->orient() = prevH->orient();
	    }
	    else {
		pathID++;
		curH->path() = pathID;
		curH->orient() = true;
	    }

	    prevH	= curH;
	    curH	= nextBoundary( curH );
	} while ( curH != initH );

	bound.push_back( curBound );
    }    
    pathID++;

    return pathID;
}


void updateAttr( Polyhedron & poly, 
		 vector< vector< Halfedge_handle > > & bound,
		 int & nRuns,
		 Attribute & attr )
{
    attr.nRuns()	= nRuns;
    attr.nColors()	= nRuns * 4 / 3;

    int nStars = 0;
    int nHyperbolics = 0;
    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	Vertex_handle vh = vi;
	if ( vh->nCuts() >= 3 ) {
	    nStars++;
	    if ( vh->label() == SADDLE_VERTEX ) {
		// cerr << "saddle has multiple cut edges" << endl;
		nHyperbolics++;
	    }
	}
    }
    attr.nStars()		= nStars;
    attr.nHyperbolics()	= nHyperbolics;
    
    double sum = 0.0;
    // Initalize boundary path IDs
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	Halfedge_handle hh = hi;
	if ( ( ! hh->connect() ) &&
	     ( hh->opposite()->vertex()->id() < hh->vertex()->id() ) ) {
	    sum += hh->weight();
	}
    }
    attr.sumWeights() = sum;
    
//------------------------------------------------------------------------------
//	Tracking the length of each cut path
//------------------------------------------------------------------------------
    vector< int > pathLength( nRuns );
    for ( int m = 0; m < nRuns; ++m ) pathLength[ m ] = 0;
    
    for ( unsigned int i = 0; i < bound.size(); ++i ) {
	unsigned int baseIndex = 0;
	int prevID = bound[ i ][ 0 ]->path();
	for ( unsigned int j = 1; j < bound[ i ].size(); ++j ) {
#ifdef MYDEBUG
	    cerr << "[P] i = " << i << " j = " << j << " path = " << bound[ i ][ j ]->path() << endl;
#endif	// MYDEBUG
	    if ( ( bound[ i ][ j ]->vertex()->nCuts() == 1 ) ||
		 ( bound[ i ][ j ]->vertex()->nCuts() >= 3 ) ) {
		baseIndex = j;
		break;
	    }
	}
	
	unsigned int k = baseIndex;
	prevID = bound[ i ][ k ]->path();
	int count = 1;
#ifdef MYDEBUG
	cerr << "[A] i = " << i << " k = " << k << " path = " << prevID << " count = " << count << endl;
#endif	// MYDEBUG
	while ( true ) {
	    k = ( k + 1 ) % bound[ i ].size();
	    
	    if ( ( bound[ i ][ k ]->vertex()->nCuts() == 1 ) ||
		 ( bound[ i ][ k ]->vertex()->nCuts() >= 3 ) ) {
		if ( count > pathLength[ prevID ] ) pathLength[ prevID ] = count;
		count = 1;
	    }
	    else {
		count++;
	    }
#ifdef MYDEBUG
	    cerr << "[B] i = " << i << " k = " << k << " path = " << prevID << " count = " << count << endl;
#endif	// MYDEBUG
	    
	    if ( k == baseIndex ) break;
	}
    }
    
    int maxLength = 0;
    for ( int m = 0; m < nRuns; ++m ) {
#ifdef MYDEBUG
	cerr << " pathID = " << setw( 3 ) << m << " Number of edges = " << setw( 4 ) << pathLength[ m ] << endl;
#endif	// MYDEBUG
	if ( pathLength[ m ] > maxLength ) {
	    maxLength = pathLength[ m ];
	}
    }
    attr.maxLength() = maxLength;
#ifdef MYDEBUG
    cerr << " maxLength = " << maxLength << endl;#
#endif	// MYDEBUG
						     
						     
//------------------------------------------------------------------------------
//	Computing the dual edge length between a pair of cut edges
//------------------------------------------------------------------------------

    if ( bound.size() <= DETAILED_CHECK_LIMIT ) {
	double aveDist = 0.0;
	double maxDist = 0.0;
	int nPairs = 0;
	for ( unsigned int i = 0; i < bound.size(); ++i ) {
	    for ( unsigned int j = 0; j < bound[ i ].size(); ++j ) {
		Halfedge_handle origH = bound[ i ][ j ];
		Halfedge_handle destH = bound[ i ][ j ]->opposite();
		Facet_handle origF = origH->facet();
		Facet_handle destF = destH->facet();
		assert( origF->id() != destF->id() );
		if ( origF->id() < destF->id() ) {
		    double curDist = ( double )dualLength( origF, destF, poly );
		    if ( curDist > maxDist ) maxDist = curDist;
		    aveDist += curDist;
		    nPairs++;
#ifdef DEBUG
		    cerr << " Dual distance between " << origF->id() << " and " << destF->id()
			 << " = " << curDist << " average = " << aveDist << endl;
#endif	// DEBUG
		}
	    }
	}
	
	attr.maxDist() = maxDist/(double)poly.size_of_facets();
	attr.aveDist() = aveDist/(double)(poly.size_of_facets()*nPairs);
#ifdef DEBUG
	cerr << " maximum distance of dual edges = " << attr.maxDist() << endl;
	cerr << " avelrage distance of dual edges = " << attr.aveDist() << endl;
#endif	// DEBUG
    }
    else {
	cerr << " SKIPPED " << endl;
	attr.maxDist() = 1.0;
	attr.aveDist() = 1.0;
    }
	
#ifdef MYDEBUG
    cerr << " cutPath : nRun = " << setw( 5 ) << attr.nRuns()
	 << " sumWeights = " << attr.sumWeights() << endl;
    for ( unsigned int k = 0; k < bound.size(); ++k ) {
	cerr << " ****** Boundary No. " << setw( 3 ) << k << endl;
	for ( unsigned int m = 0; m < bound[ k ].size(); ++m ) {
	    cerr << " boundary edge [ " << setw( 3 ) << m << " ] : "
		 << setw( 3 ) << bound[ k ][ m ]->opposite()->vertex()->id()
		 << " -- " << setw( 3 ) << bound[ k ][ m ]->vertex()->id()
		 << " ID : " << setw( 3 ) << bound[ k ][ m ]->path()
		 << endl;
	}
    }
#endif	// MYDEBUG
}


//------------------------------------------------------------------------------
//	Boundary Edge Reordering
//------------------------------------------------------------------------------
void evalBoundarySharpness( vector< vector< Halfedge_handle > > & vvhh, 
			    vector< vector< double > > & vvd )
{
    vvd.clear();

    for ( unsigned int k = 0; k < vvhh.size(); ++k ) {
	vector< double > seq;
	for ( unsigned int m = 0; m < vvhh[ k ].size(); ++m ) {
	    seq.push_back( fabs( evalConvexity( vvhh[ k ][ m ] ) ) );
	}
	vvd.push_back( seq );
    }
}


void propagateBoundarySharpness( vector< Halfedge_handle > & vhh, vector< double > & vd )
{
    assert( vhh.size() == vd.size() );

    unsigned int N = vhh.size();
    for ( unsigned int m = 0; m < N; ++m ) {
	Halfedge_handle curHH = vhh[ m ];
	unsigned int nFore = 0, nBack = 0;
	for ( unsigned int j = 1; j < N; j++ ) {
	    if ( vhh[ ( m + j ) % N ]->opposite() == curHH ) break;
	    else nFore++;
	}
	for ( unsigned int j = 1; j < N; j++ ) {
	    if ( vhh[ ( m + N - j ) % N ]->opposite() == curHH ) break;
	    else nBack++;
	}
	// cerr << " m = " << m << " nFore = " << nFore << " nBack = " << nBack << endl;
	if ( nFore <= nBack ) {
	    for ( unsigned int j = 1; j <= nFore; j++ ) {
		if ( vd[ ( m + j ) % N ] < vd[ m ] ) {
		    vd[ ( m + j ) % N ] = vd[ m ];
		}
	    }
	}
	else {
	    for ( unsigned int j = 1; j <= nBack; j++ ) {
		if ( vd[ ( m + N - j ) % N ] < vd[ m ] ) {
		    vd[ ( m + N - j ) % N ] = vd[ m ];
		}
	    }
	}
    }
}


void reorderBoundaryEdges( vector< Halfedge_handle > & vhh, vector< double > & vd )
{
    assert( vhh.size() == vd.size() );

    vector< int > match, invmatch;
    vector< Halfedge_handle > vhhT = vhh;
    vector< double > vdT = vd;

#ifdef DEBUG
//------------------------------------------------------------------------------
//	Show the topological surgery encoding
//------------------------------------------------------------------------------
    cerr << "BEFORE: " << endl;
    for ( unsigned int m = 0; m < vhh.size(); ++m ) {
	cerr << "[" << vhh[ m ]->path() << "]";
	if ( vhh[ m ]->orient() ) cerr << "^(+1)";
	else cerr << "^(-1)";
	cerr << "(" << vd[ m ] << ")";
	cerr << "=";
    }
#endif	// DEBUG

//------------------------------------------------------------------------------
    while ( vhhT.size() > 0 ) {
	int maxPath = NO_INDEX;
	double maxSharpness = -100.0;
	for ( unsigned int m = 0; m < vhhT.size(); ++m ) {
	    if ( ( vhhT[ m ]->path() == vhhT[ (m+1)%vhhT.size() ]->path() ) &&
		 ( ( ( vhhT[ m ]->orient() ) && ( ! vhhT[ (m+1)%vhhT.size() ]->orient() ) ) || 
		   ( ( ! vhhT[ m ]->orient() ) && ( vhhT[ (m+1)%vhhT.size() ]->orient() ) ) ) ) {
		// cerr << " m = " << m << " passes" << endl;
		if ( vdT[ m ] > maxSharpness ) {
		    maxSharpness = vdT[ m ];
		    maxPath = vhhT[ m ]->path();
		}
	    }
	}

	// Sphere case
	if ( maxPath != NO_INDEX ) {
	    match.push_back( maxPath );
	    for ( int m = (int)vhhT.size() - 1; m >= 0; --m ) {
		// cerr << " m = " << m << endl;
		if ( vhhT[ m ]->path() == maxPath ) {
		    // cerr << " Deleting m = " << m << endl;
		    vhhT.erase( vhhT.begin() + m );
		    vdT.erase( vdT.begin() + m );
		}
	    }
	    // cerr << "The end" << endl;
	}
	// Torus case
	else {
	    for ( unsigned int m = 0; m < vhhT.size(); ++m ) {
		bool isExist = false;
		for ( unsigned int k = 0; k < match.size(); ++k ) {
		    if ( match[ k ] == vhhT[ m ]->path() ) {
			isExist = true;
			break;
		    }
		}
		if ( ! isExist ) match.push_back( vhhT[ m ]->path() );
	    }
	    break;
	}
    }

    invmatch.resize( match.size() );
    for ( unsigned int k = 0; k < match.size(); ++k ) {
	invmatch[ match[ k ] ] = k;
    }

    // Reordering
    for ( unsigned int m = 0; m < vhh.size(); ++m ) {
	// cerr << " Change path ID " << vhh[ m ]->path() << " into " << invmatch[ vhh[ m ]->path() ] << endl;
	vhh[ m ]->path() = invmatch[ vhh[ m ]->path() ];
    }

#ifdef SKIP
//------------------------------------------------------------------------------
//	Show the topological surgery encoding
//------------------------------------------------------------------------------
    cerr << "AFTER: " << endl;
    for ( unsigned int m = 0; m < vhh.size(); ++m ) {
	cerr << "[" << vhh[ m ]->path() << "]";
	if ( vhh[ m ]->orient() ) cerr << "^(+1)";
	else cerr << "^(-1)";
	cerr << "(" << vd[ m ] << ")";
	cerr << "=";
    }
#endif	// SKIP
}



// Outline the boundary edges of each unfoled patch as unvisited
void traverseBoundary( Polyhedron & poly, Attribute & attr, vector< vector< Halfedge_handle > > & bound )
{
    int nRuns;

    nRuns = extractBoundary	( poly, bound ); 
    updateAttr	( poly, bound, nRuns, attr );
    vector< vector< double > > sharpness;
    evalBoundarySharpness( bound, sharpness );
    if ( bound.size() == 1 ) {
	propagateBoundarySharpness( bound[ 0 ], sharpness[ 0 ] );
	reorderBoundaryEdges( bound[ 0 ], sharpness[ 0 ] );
    }
}

// Outline the boundary edges of each unfoled patch as unvisited
void traverseBoundary( Polyhedron & poly, Attribute & attr )
{
    // int nRuns;
    vector< vector< Halfedge_handle > > bound;

    // nRuns = extractBoundary	( poly, bound ); 
    // updateAttr	( poly, bound, nRuns, attr );
    traverseBoundary( poly, attr, bound );
}


// Outline the boundary edges of each unfoled patch as unvisited
void labelBoundary( Polyhedron & poly, Attribute & attr )
{
    int nRuns;
    vector< vector< Halfedge_handle > > bound;

    initEdgeAttr	( poly );
    nRuns = extractBoundary	( poly, bound ); 
    updateAttr	( poly, bound, nRuns, attr );
    vector< vector< double > > sharpness;
    evalBoundarySharpness( bound, sharpness );
    if ( bound.size() == 1 ) {
	propagateBoundarySharpness( bound[ 0 ], sharpness[ 0 ] );
	reorderBoundaryEdges( bound[ 0 ], sharpness[ 0 ] );
    }
}



//------------------------------------------------------------------------------
//
//	Functions for traversing dual meshes
//
//------------------------------------------------------------------------------
// Find the mutual intersection between a pair of unfolded triangular faces
void findIntersections( vector< Facet_handle > & vfh,
			vector< pair< Facet_handle, Facet_handle > > & endfaces ) 
{
    // initialize the list
    // This does not applied. 
    // endfaces.clear();

    // for each face
    for ( unsigned int i = 0; i < vfh.size(); ++i ) {
	Facet_handle fi = vfh[ i ];
	for ( unsigned int j = i+1; j < vfh.size(); ++j ) {
	    Facet_handle fj = vfh[ j ];
	    if ( isIntersected( fi->triangle(), fj->triangle() ) ) {
//		endfaces.push_back( make_pair< Facet_handle, Facet_handle >( fi, fj ) ); error happend
			endfaces.push_back( make_pair( fi, fj ) );
	    }
	}
    }

#ifdef DEBUG
    cerr << "Total number of endfaces pair = " << endfaces.size() << endl;
    for ( unsigned int k = 0; k < endfaces.size(); ++k ) {
	cerr << "Intersected face pair [ " << setw( 3 ) << k << " ]: "
	     << setw( 4 ) << endfaces[ k ].first->id() << " -- " 
	     << setw( 4 ) << endfaces[ k ].second->id()
	     << endl;
    }
#endif	// DEBUG
}    


// Find the mutual intersection between a pair of unfolded triangular faces
void findIntersections( vector< Facet_handle > & setA, vector< Facet_handle > & setB,
			vector< pair< Facet_handle, Facet_handle > > & endfaces ) 
{
    // initialize the list
    // This does not applied. 
    // endfaces.clear();

    // for each face
    for ( unsigned int i = 0; i < setA.size(); ++i ) {
	Facet_handle fi = setA[ i ];
	for ( unsigned int j = 0; j < setB.size(); ++j ) {
	    Facet_handle fj = setB[ j ];
	    if ( isIntersected( fi->triangle(), fj->triangle() ) ) {
//		endfaces.push_back( make_pair< Facet_handle, Facet_handle >( fi, fj ) );
			endfaces.push_back( make_pair( fi, fj ) );
	    }
	}
    }
}    


// Find the path between two face that have intersection with each other.
bool findDualPath( const Facet_handle & fhO, // origin face of the path
		   const Facet_handle & fhD, // destination face of the path
		   Polyhedron & poly,
		   vector< Facet_handle > & path )
{    
    // initialize the path
    path.clear();
    // list of pointers to the previous face
    vector< Facet_handle >		back;
    // list of current wavefront faces
    vector< Halfedge_handle >		stepP, stepN;
    
#ifdef MYDEBUG
    cerr << " Origin : " << fhO->id() << " --> Destination : " << fhD->id() << endl;
#endif	// MYDEBUG

    // initialize the capacity of the list of previous faces
    back.resize( poly.size_of_facets() );
    for ( unsigned int k = 0; k < back.size(); ++k ) back[ k ] = NULL;

    // initialize the starting point
    Halfedge_handle hhC = fhD->halfedge();
    Halfedge_handle hhL = hhC->prev();
    Halfedge_handle hhR = hhC->next();

    if ( hhC->connect() ) {
	back[ hhC->opposite()->facet()->id() ] = hhC->facet();
	stepN.push_back( hhC->opposite() );
    }
    if ( hhL->connect() ) {
	back[ hhL->opposite()->facet()->id() ] = hhL->facet();
	stepN.push_back( hhL->opposite() );
    }
    if ( hhR->connect() ) {
	back[ hhR->opposite()->facet()->id() ] = hhR->facet();
	stepN.push_back( hhR->opposite() );
    }

    // step forward up to the destination face
    while ( stepN.size() != 0 ) {

	stepP = stepN;
	stepN.clear();
	
	for ( unsigned int k = 0; k < stepP.size(); ++k ) {
	    if ( stepP[ k ]->facet()->id() == fhO->id() ) {
		// the search reaches the other end face
		int curFID = fhO->id();
#ifdef MYDEBUG
		cerr << "*** Path (Face IDs) : " << ends;
#endif	// MYDEBUG
		path.push_back( fhO );
#ifdef MYDEBUG
		cerr << setw( 4 ) << fhO->id();
#endif	// MYDEBUG
		do {
		    if ( back[ curFID ] == NULL ) {
			cerr << " Invalid path at " << curFID << endl;
		    }
		    path.push_back( back[ curFID ] );
#ifdef MYDEBUG
		    cerr << setw( 4 ) << back[ curFID ]->id();
#endif	// MYDEBUG
		    curFID = back[ curFID ]->id();
		} while ( curFID != fhD->id() );
#ifdef MYDEBUG
		cerr << endl;
#endif	// MYDEBUG
		return true;
	    }
	}

#ifdef MYDEBUG
	for ( unsigned int k = 0; k < stepP.size(); ++k ) {
	    cerr << "Last Step[ " << setw( 3 ) << k << " ]: ";
	    if ( stepP[ k ] == NULL ) cerr << "NULL";
	    else cerr << setw( 4 ) << stepP[ k ]->facet()->id();
	    cerr << endl;
	}
	getchar();
#endif	// MYDEBUG
	
	// This is necessary because step.size() will be incremented in this loop.
	// unsigned int limit = stepP.size();
	for ( unsigned int k = 0; k < stepP.size(); ++k ) {
	    Halfedge_handle curL = stepP[ k ]->prev();
	    Halfedge_handle curR = stepP[ k ]->next();
	    
#ifdef MYDEBUG
	    cerr << " curL : " << curL->opposite()->vertex()->id()
		 << " == " << curL->vertex()->id() << endl;
	    cerr << " curR : " << curR->opposite()->vertex()->id()
		 << " == " << curR->vertex()->id() << endl;
#endif	// MYDEBUG
	    
	    bool isBoth = false;
	    if ( curL->connect() && curR->connect() )
		isBoth = true;
	    
	    if ( isBoth ) {
		back[ curR->opposite()->facet()->id() ] = curR->facet();
		stepN.push_back( curR->opposite() );
		back[ curL->opposite()->facet()->id() ] = curL->facet();
		stepN.push_back( curL->opposite() );
	    }
	    else if ( curL->connect() ) {
		back[ curL->opposite()->facet()->id() ] = curL->facet();
		stepN.push_back( curL->opposite() );
	    }
	    else if ( curR->connect() ) {
		back[ curR->opposite()->facet()->id() ] = curR->facet();
		stepN.push_back( curR->opposite() );
	    }
	}
    }
    return false;
}


// Find the path between two face that have intersection with each other.
unsigned int dualLength( const Facet_handle & fhO, // origin face of the path
			 const Facet_handle & fhD, // destination face of the path
			 Polyhedron & poly )
{    
    vector< Facet_handle > path;
    if ( findDualPath( fhO, fhD, poly, path ) ) {
	return path.size();
    }
    else {
	return poly.size_of_facets();
    }
}



// Find the path between two face that have intersection with each other.
void fillPatch( const Facet_handle & fh, // starting face of the path
		vector< Facet_handle > & patch )
{    
    // initialize the patch
    patch.clear();
    // list of current wavefront faces
    vector< Halfedge_handle >		stepP, stepN;
    
    patch.push_back( fh );

    // initialize the starting point
    Halfedge_handle hhC = fh->halfedge();
    Halfedge_handle hhL = hhC->prev();
    Halfedge_handle hhR = hhC->next();

    if ( hhC->connect() ) {
	stepN.push_back( hhC->opposite() );
	patch.push_back( hhC->opposite()->facet() );
    }
    if ( hhL->connect() ) {
	stepN.push_back( hhL->opposite() );
	patch.push_back( hhL->opposite()->facet() );
    }
    if ( hhR->connect() ) {
	stepN.push_back( hhR->opposite() );
	patch.push_back( hhR->opposite()->facet() );
    }

    // step forward up to the destination face
    while ( stepN.size() != 0 ) {

	stepP = stepN;
	stepN.clear();

	for ( unsigned int k = 0; k < stepP.size(); ++k ) {
	    Halfedge_handle curL = stepP[ k ]->prev();
	    Halfedge_handle curR = stepP[ k ]->next();
	    
#ifdef MYDEBUG
	    cerr << " curL : " << curL->opposite()->vertex()->id()
		 << " == " << curL->vertex()->id() << endl;
	    cerr << " curR : " << curR->opposite()->vertex()->id()
		 << " == " << curR->vertex()->id() << endl;
#endif	// MYDEBUG
	    
	    if ( curL->connect() && curR->connect() ) {
		stepN.push_back( curR->opposite() );
		patch.push_back( curR->opposite()->facet() );
		stepN.push_back( curL->opposite() );
		patch.push_back( curL->opposite()->facet() );
	    }
	    else if ( curL->connect() ) {
		stepN.push_back( curL->opposite() );
		patch.push_back( curL->opposite()->facet() );
	    }
	    else if ( curR->connect() ) {
		stepN.push_back( curR->opposite() );
		patch.push_back( curR->opposite()->facet() );
	    }
	}
    }
    return;
}


// Minimum cover set
void minimumCover( Polyhedron & poly,
		   vector< pair< Facet_handle, Facet_handle > > & endfaces,
		   vector< vector< Facet_handle > > & pathlist,
		   vector< Halfedge_handle > & cut )
{
    const int unitDegree = 1;
    
    // initialize the list of cut edges
    cut.clear();

    // prepare the list of boolean values for representing whether path is cut
    // or not.
    vector< bool > isCut( endfaces.size() );
    for ( unsigned int i = 0; i < isCut.size(); ++i ) isCut[ i ] = false;
    unsigned int nSplits = 0;

    cerr << " Size of isCut = " << isCut.size() << endl;

    while ( nSplits < endfaces.size() ) {

	// prepare the list of edge weights that count the number of paths over them
	vector< int > degree( poly.size_of_halfedges()/2 );
	for ( unsigned int k = 0; k < degree.size(); ++k ) degree[ k ] = 0;
	
	// count the number of paths over each edge
	for ( unsigned int i = 0; i < isCut.size(); ++i ) {
	    if ( ! isCut[ i ] ) {
		for ( unsigned int j = 0; j < pathlist[ i ].size() - 1; ++j ) {
		    Halfedge_handle hh = sharedHalfedge( pathlist[ i ][ j ], pathlist[ i ][ j + 1 ] );
		    degree[ hh->id() ] += unitDegree;
#ifdef MYDEBUG
		    cerr << "Edge : " << hh->opposite()->vertex()->id() 
			 << " == " << hh->vertex()->id() 
			 << " ::: degree[ " << hh->id() << " ] = " << degree[ hh->id() ] << endl;
#endif	// MYDEBUG
		}
	    }
	}

	// find the edge having the largest weight
	int maxDegree = 1;
	// double maxWeight = 0.0;
	double maxWeight = -1.0e+8;
	Halfedge_handle maxHH = NULL;
	for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
#ifdef MYDEBUG
	    cerr << " Edge : " << hi->opposite()->vertex()->id() 
		 << " == " << hi->vertex()->id() << " now checking." << endl;
#endif	// MYDEBUG
	    if ( maxDegree < degree[ hi->id() ] ) {
		maxDegree = degree[ hi->id() ];
		maxWeight = hi->weight();
		maxHH = hi;
#ifdef MYDEBUG
		cerr << "[A] degree = " << maxDegree << " weight = " << maxWeight
		     << " Edge " << maxHH->opposite()->vertex()->id() << " == " 
		     << maxHH->vertex()->id() << endl;
#endif	// MYDEBUG
	    }
	    else if ( maxDegree == degree[ hi->id() ] ) {
		if ( maxWeight < hi->weight() ) {
		    maxWeight = hi->weight();
		    maxHH = hi;
#ifdef MYDEBUG
		    cerr << "[B] degree = " << maxDegree << " weight = " << maxWeight
			 << " Edge " << maxHH->opposite()->vertex()->id() << " == " 
			 << maxHH->vertex()->id() << endl;
#endif	// MYDEBUG
		}
	    }
	}

	assert( maxHH != NULL );
#ifdef MYDEBUG
	cerr << " Select the dual edge : " 
	     << " number of overlaps = " << maxWeight << " :: "
	     << setw( 3 ) << maxHH->face()->id()
	     << " == " << setw( 3 ) << maxHH->opposite()->face()->id() << endl;
#endif	// MYDEBUG

	cut.push_back( maxHH );

	// remove the path that contains the cut edge
	for ( unsigned int i = 0; i < isCut.size(); ++i ) {
	    if ( ! isCut[ i ] ) {
		for ( unsigned int j = 0; j < pathlist[ i ].size() - 1; ++j ) {
		    Halfedge_handle hh = sharedHalfedge( pathlist[ i ][ j ], pathlist[ i ][ j + 1 ] );
		    if ( hh->id() == maxHH->id() ) {
			isCut[ i ] = true;
			nSplits++;
		    }
		    if ( isCut[ i ] ) break;
		}
	    }
	}
#ifdef MYDEBUG
	cerr << " nSplits = " << nSplits << endl;
	for ( unsigned int i = 0; i < isCut.size(); ++i ) {
	    if ( ! isCut[ i ] ) {
		cerr << "+++ Path[ " << setw( 3 ) << i << " ] : ";
		for ( unsigned int k = 0; k < pathlist[ i ].size(); ++k ) 
		    cerr << setw( 4 ) << pathlist[ i ][ k ]->id();
		cerr << endl;
	    }
	}
#endif	// MYDEBUG
    }

#ifdef MYDEBUG
    cerr << "List of minimum cover set : size = " << cut.size() << endl;
    for ( unsigned int m = 0; m < cut.size(); ++m ) {
	cerr << "[ " << setw( 3 ) << m << " ] : " << setw( 4 ) 
	     << cut[ m ]->facet()->id() << " == " << cut[ m ]->opposite()->facet()->id()
	     << "  (VID : " << setw( 3 ) << cut[ m ]->opposite()->vertex()->id()
	     << " -- " << setw( 3 ) << cut[ m ]->vertex()->id() << ")" << endl;
    }
#endif	// MYDEBUG

    return;
}


// Minimum cover set
bool constrainedMinimumCover( Polyhedron & poly,
			      vector< pair< Facet_handle, Facet_handle > > & endfaces,
			      vector< vector< Facet_handle > > & pathlist,
			      vector< Halfedge_handle > & cut )
{
    const int unitDegree = 1;
    
    // initialize the list of cut edges
    cut.clear();

    // prepare the list of boolean values for representing whether path is cut
    // or not.
    vector< bool > isCut( endfaces.size() );
    for ( unsigned int i = 0; i < isCut.size(); ++i ) isCut[ i ] = false;
    unsigned int nSplits = 0;

#ifdef MYDEBUG
    cerr << " Size of isCut = " << isCut.size() << endl;
#endif	// MYDEBUG

    while ( nSplits < endfaces.size() ) {

	// prepare the list of edge weights that count the number of paths over them
	vector< int > degree( poly.size_of_halfedges()/2 );
	for ( unsigned int k = 0; k < degree.size(); ++k ) degree[ k ] = 0;
	
	// count the number of paths over each edge
	for ( unsigned int i = 0; i < isCut.size(); ++i ) {
	    if ( ! isCut[ i ] ) {
		for ( unsigned int j = 0; j < pathlist[ i ].size() - 1; ++j ) {
		    Halfedge_handle hh = sharedHalfedge( pathlist[ i ][ j ], pathlist[ i ][ j + 1 ] );
		    degree[ hh->id() ] += unitDegree;
		}
	    }
	}

	// find the edge having the largest weight
	int maxDegree = 1;
	double maxWeight = 0.0;
	Halfedge_handle maxHH = NULL;
	for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	    if ( hi->visit() ) {
#ifdef MYDEBUG
		cerr << " Edge : " << hi->opposite()->vertex()->id() 
		     << " == " << hi->vertex()->id() << " has been rejected" << endl;
#endif	// MYDEBUG
		; // do nothing
	    }
	    else if ( maxDegree < degree[ hi->id() ] ) {
		maxDegree = degree[ hi->id() ];
		maxWeight = hi->weight();
		maxHH = hi;
#ifdef MYDEBUG
		cerr << "[A] degree = " << maxDegree << " weight = " << maxWeight
		     << " Edge " << maxHH->opposite()->vertex()->id() << " == " 
		     << maxHH->vertex()->id() << endl;
#endif	// MYDEBUG
	    }
	    else if ( maxDegree == degree[ hi->id() ] ) {
		if ( maxWeight < hi->weight() ) {
		    maxWeight = hi->weight();
		    maxHH = hi;
#ifdef MYDEBUG
		    cerr << "[B] degree = " << maxDegree << " weight = " << maxWeight
			 << " Edge " << maxHH->opposite()->vertex()->id() << " == " 
			 << maxHH->vertex()->id() << endl;
#endif	// MYDEBUG
		}
	    }
	}

	if ( maxHH == NULL ) return false;

	cut.push_back( maxHH );

	// remove the path that contains the cut edge
	for ( unsigned int i = 0; i < isCut.size(); ++i ) {
	    if ( ! isCut[ i ] ) {
		for ( unsigned int j = 0; j < pathlist[ i ].size() - 1; ++j ) {
		    Halfedge_handle hh = sharedHalfedge( pathlist[ i ][ j ], pathlist[ i ][ j + 1 ] );
		    if ( hh->id() == maxHH->id() ) {
			isCut[ i ] = true;
			nSplits++;
		    }
		    if ( isCut[ i ] ) break;
		}
	    }
	}
    }

#ifdef MYDEBUG
    cerr << "List of constrained minimum cover set : size = " << cut.size() << endl;
    for ( unsigned int m = 0; m < cut.size(); ++m ) {
	cerr << "[ " << setw( 3 ) << m << " ] : " << setw( 4 ) 
	     << cut[ m ]->facet()->id() << " == " << cut[ m ]->opposite()->facet()->id()
	     << "  (VID : " << setw( 3 ) << cut[ m ]->opposite()->vertex()->id()
	     << " -- " << setw( 3 ) << cut[ m ]->vertex()->id() << ")" << endl;
    }
#endif	// MYDEBUG

    return true;
}


// intentionally cut the strips at the specified edge
void cutStrips( vector< vector< Facet_handle > > & patch,
		Halfedge_handle & cut )
{
#ifdef MYDEBUG
    cerr << "Input cut edge : "
	 << cut->facet()->id()
	 << " == " << cut->opposite()->facet()->id() << " (VID : "
	 << cut->opposite()->vertex()->id()
	 << " -- " << cut->vertex()->id() << ") " << endl;
#endif	// MYDEBUG

    //------------------------------------------------------------------------------
    //	If the cut edge is in the middle of a triangular strip
    //------------------------------------------------------------------------------
    for ( unsigned int i = 0; i < patch.size(); ++i ) {
	for ( unsigned int j = 0; j < patch[ i ].size(); ++j ) {
	    
	    Facet_handle fh = patch[ i ][ j ];
	    // check all the three halfedges
	    Halfedge_handle hhC = fh->halfedge();
	    Halfedge_handle hhL = hhC->prev();
	    Halfedge_handle hhR = hhC->next();
	    
	    vector< Halfedge_handle > cand;
	    
	    if ( hhC->connect() ) cand.push_back( hhC );
	    if ( hhL->connect() ) cand.push_back( hhL );
	    if ( hhR->connect() ) cand.push_back( hhR );
	    
	    for ( unsigned k = 0; k < cand.size(); ++k ) {
		if ( ( cand[ k ] == cut ) || ( cand[ k ]->opposite() == cut ) ) {

		    //------------------------------------------------------------------------------
		    cut->connect() = cut->opposite()->connect() = false;
		    cut->opposite()->vertex()->nCuts()++;
		    cut->vertex()->nCuts()++;
		    //------------------------------------------------------------------------------

		    vector< Facet_handle >		patchF, patchL;
		    fillPatch( cand[ k ]->facet(), patchF );
		    fillPatch( cand[ k ]->opposite()->facet(), patchL );

#ifdef MYDEBUG
		    cerr << "[patchF] : size = " << patchF.size() << endl;
		    for ( unsigned int m = 0; m < patchF.size(); ++m ) {
			cerr << " m = " << setw( 3 ) << m << " FaceID = " << setw( 3 ) << patchF[ m ]->id() << endl;
		    }

		    cerr << "[patchL] : size = " << patchL.size() << endl;
		    for ( unsigned int m = 0; m < patchL.size(); ++m ) {
			cerr << " m = " << setw( 3 ) << m << " FaceID = " << setw( 3 ) << patchL[ m ]->id() << endl;
		    }
#endif	// MYDEBUG

		    patch[ i ] = patchF;
		    unsigned int curID = patch.size();
		    for ( unsigned int m = 0; m < patchL.size(); ++m ) {
			if ( patchL[ m ]->halfedge()->connect() ) {
			    patchL[ m ]->halfedge()->cycle() = curID;
			    patchL[ m ]->halfedge()->opposite()->cycle() = curID;
			}
			if ( patchL[ m ]->halfedge()->prev()->connect() ) {
			    patchL[ m ]->halfedge()->prev()->cycle() = curID;
			    patchL[ m ]->halfedge()->prev()->opposite()->cycle() = curID;
			}
			if ( patchL[ m ]->halfedge()->next()->connect() ) {
			    patchL[ m ]->halfedge()->next()->cycle() = curID;
			    patchL[ m ]->halfedge()->next()->opposite()->cycle() = curID;
			}
			patchL[ m ]->piece() = curID;
		    }
		    patch.push_back( patchL );
		    
		    
#ifdef MYDEBUG
		    for ( unsigned int k = 0; k < formerHH.size(); ++k )
			cerr << "First  : [" << k << "] : "
			     << formerHH[ k ]->opposite()->vertex()->id()
			     << " -- " << formerHH[ k ]->vertex()->id() 
			     << " edgeID = " 
			     << formerHH[ k ]->cycle() 
			     << " faceID = " << formerHH[ k ]->face()->piece() << endl;
		    
		    for ( unsigned int k = 0; k < latterHH.size(); ++k )
			cerr << "Second : [" << k << "] : "
			     << latterHH[ k ]->opposite()->vertex()->id()
			     << " -- " << latterHH[ k ]->vertex()->id() 
			     << " edgeID = " 
			     << latterHH[ k ]->cycle() 
			     << " faceID = " << latterHH[ k ]->face()->piece() << endl;
#endif	// MYDEBUG

#ifdef MYDEBUG
		    for ( unsigned int m = 0; m < spine.size(); ++m ) {
			cerr << "Strip  : [" << m << "] : " 
			     << " nEdges = " << spine[ m ].size() << endl;
			for ( unsigned n = 0; n < spine[ m ].size(); ++n ) { 
			    cerr << " Edge : " << spine[ m ][ n ]->opposite()->vertex()->id()
				 << " -- " << spine[ m ][ n ]->vertex()->id() 
				 << " edgeID = " 
				 << spine[ m ][ n ]->cycle() 
				 << " faceID = " << spine[ m ][ n ]->facet()->piece() << endl;
			}
		    }
#endif	// MYDEBUG
		    
#ifdef MYDEBUG
		    cerr << "[OO] Inserted the intentional Cut at : " << setw( 3 )
			 << cut->opposite()->vertex()->id() << " -- " << setw( 3 )
			 << cut->vertex()->id() << endl;
#endif	// MYDEBUG
		    return;
		}
	    }
	}
    }
#ifdef MYDEBUG
    cerr << "[XX] No need to have the intentional Cut at : " << setw( 3 )
	 << cut->opposite()->vertex()->id() << " -- " << setw( 3 )
	 << cut->vertex()->id() << endl;
#endif	// MYDEBUG
}


//	resolve overlaps on the triangular strip
void resolveOverlaps( Polyhedron & poly,
		      vector< vector< Facet_handle > > & patch )
{
    vector< pair< Facet_handle, Facet_handle > > endfaces;
    vector< Facet_handle > path;
    vector< vector< Facet_handle > > pathlist;
    vector< Halfedge_handle > cutList;
    vector< Halfedge_handle > joint;

    // This should be called here.
    endfaces.clear();
    for ( unsigned int i = 0; i < patch.size(); ++i ) 
	findIntersections( patch[ i ], endfaces );

    pathlist.clear();
    for ( unsigned int k = 0; k < endfaces.size(); ++k ) {
	if ( findDualPath( endfaces[ k ].first, endfaces[ k ].second,
			  poly, path ) ) {
	    pathlist.push_back( path );
	}
    }

    minimumCover( poly, endfaces, pathlist, cutList );

#ifdef MYDEBUG
    cerr << "Number of cuts = " << cutList.size() << endl;
#endif	// MYDEBUG

    // The followings must be checked 
    // MSTtoStrips( poly, spine, joint );
    // double limitCos = INNER_PROD_ALL;
    // MSTMerge( poly, spine, joint, patch, sheet, bound );
    for ( unsigned int m = 0; m < cutList.size(); ++m ) {
	// cutStrips( spine, patch, joint, cutList[ m ] );
	cutStrips( patch, cutList[ m ] );
    }

    return;
}		      


// resolve overlaps on the triangular strip
bool resplitOverlaps( Polyhedron & poly,
		      vector< vector< Facet_handle > > & patch )
{
    vector< pair< Facet_handle, Facet_handle > > endfaces;
    vector< Facet_handle > path;
    vector< vector< Facet_handle > > pathlist;
    vector< Halfedge_handle > cutList;
    vector< Halfedge_handle > joint;

    // This should be called here.
    endfaces.clear();
    for ( unsigned int i = 0; i < patch.size(); ++i ) 
	findIntersections( patch[ i ], endfaces );

    pathlist.clear();
    for ( unsigned int k = 0; k < endfaces.size(); ++k ) {
	if ( findDualPath( endfaces[ k ].first, endfaces[ k ].second,
			  poly, path ) ) {
	    pathlist.push_back( path );
	}
    }

    bool canFind = constrainedMinimumCover( poly, endfaces, pathlist, cutList );

    //#ifdef MYDEBUG
    cerr << " Number of cuts = " << cutList.size() << endl;
    //#endif	// MYDEBUG

    if ( ! canFind ) return false;

    // The followings must be checked 
    // MSTtoStrips( poly, spine, joint );
    // double limitCos = INNER_PROD_ALL;
    // MSTMerge( poly, spine, joint, patch, sheet, bound );
    for ( unsigned int m = 0; m < cutList.size(); ++m ) {
	// cutStrips( spine, patch, joint, cutList[ m ] );
	cutStrips( patch, cutList[ m ] );
	cutList[ m ]->visit() = cutList[ m ]->opposite()->visit() = true;
    }

#ifdef NO_NEED
    embedCycles( patch, sheet, bound );
#endif	// NO_NEED

    return true;
}		      


//	resolve overlaps on the triangular strip
bool resplitOverlaps( Polyhedron & poly,
		      vector< vector< Facet_handle > > & patch,
		      vector< Facet_handle > & setA, vector< Facet_handle > & setB )
{
    vector< pair< Facet_handle, Facet_handle > > endfaces;
    vector< Facet_handle > path;
    vector< vector< Facet_handle > > pathlist;
    vector< Halfedge_handle > cutList;
    vector< Halfedge_handle > joint;

    // This should be called here.
    endfaces.clear();
    findIntersections( setA, setB, endfaces );

    pathlist.clear();
    for ( unsigned int k = 0; k < endfaces.size(); ++k ) {
	if ( findDualPath( endfaces[ k ].first, endfaces[ k ].second,
			  poly, path ) ) {
	    pathlist.push_back( path );
	}
    }

    bool canFind = constrainedMinimumCover( poly, endfaces, pathlist, cutList );

    //#ifdef MYDEBUG
    cerr << " Number of cuts = " << cutList.size() << endl;
    //#endif	// MYDEBUG

    if ( ! canFind ) return false;

    // The followings must be checked 
    // MSTtoStrips( poly, spine, joint );
    // double limitCos = INNER_PROD_ALL;
    // MSTMerge( poly, spine, joint, patch, sheet, bound );
    for ( unsigned int m = 0; m < cutList.size(); ++m ) {
	// cutStrips( spine, patch, joint, cutList[ m ] );
	cutStrips( patch, cutList[ m ] );
	cutList[ m ]->visit() = cutList[ m ]->opposite()->visit() = true;
    }

#ifdef NO_NEED
    embedCycles( patch, sheet, bound );
#endif	// NO_NEED

    return true;
}		      
