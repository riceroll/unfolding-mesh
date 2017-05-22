//------------------------------------------------------------------------------
//
//	project.cpp
//
//------------------------------------------------------------------------------
#include "common.h"
#include <CGAL/bounding_box.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

#include <iostream>
#include <iomanip>
#include <list>
#include <stack>

using namespace std;

typedef Kernel::Iso_rectangle_2 Rectangle2;


#define CUT_WITH_SELF_INTERSECTION


#ifdef OLD
Point2 coordTriangle( const Halfedge_handle & hh, const Facet_handle & fh )
{
    if ( fh->halfedge()->prev() == hh ) 
	return fh->triangle()[ 0 ];
    else if ( fh->halfedge() == hh ) 
	return fh->triangle()[ 1 ];
    else if ( fh->halfedge()->next() == hh ) 
	return fh->triangle()[ 2 ];
    else
	assert( false );
}
#endif	// OLD

//------------------------------------------------------------------------------
//	To update the coordinate for the new split triangle 
//------------------------------------------------------------------------------
void rewindTriangle( Halfedge_handle & oldbase, const Halfedge_handle & newbase )
{
    const int nD = NUM_SIDES;

    if ( oldbase->facet() != newbase->facet() ) {
	cerr << " Old face ID " << oldbase->facet()->id() << endl;
	cerr << " New face ID " << newbase->facet()->id() << endl;
	cerr << " Old face piece ID " << oldbase->facet()->piece() << endl;
	cerr << " New face piece ID " << newbase->facet()->piece() << endl;
	cerr << " Old base edge " << oldbase->vertex()->id() << " == "
	     << oldbase->opposite()->vertex()->id() << endl;
	cerr << " New base edge " << newbase->vertex()->id() << " == "
	     << newbase->opposite()->vertex()->id() << endl;
	assert( oldbase->facet() == newbase->facet() );
    }

    Point2 newcoord[ NUM_SIDES ]; // NUM_SIDES=3=nD

    if ( newbase == oldbase->prev() ) {
	for ( int k = 0; k < nD; ++k ) 
	    newcoord[ k ] = oldbase->facet()->triangle()[ k+2 ];
    }
    else if ( newbase == oldbase->next() ) {
	for ( int k = 0; k < nD; ++k ) 
	    newcoord[ k ] = oldbase->facet()->triangle()[ k+1 ];
    }
    else {
	for ( int k = 0; k < nD; ++k ) 
	    newcoord[ k ] = oldbase->facet()->triangle()[ k ];
    }

    oldbase->facet()->triangle() = Triangle2( newcoord[ 0 ], newcoord[ 1 ], newcoord[ 2 ] );
}


//------------------------------------------------------------------------------
//	seek a sequence of boundary of a patch
//------------------------------------------------------------------------------
void boundaryHalfedges( vector< Halfedge_handle > & spine, vector< Halfedge_handle > & boundary )
{
    Halfedge_handle hh = spine[ 0 ];
    boundary.clear();
    boundary.push_back( hh );
    // correct right boundary halfedges
    for ( unsigned int i = 1; i < spine.size(); ++i ) {
	if ( hh->prev()->opposite() == spine[ i ] ) {
	    boundary.push_back( hh->next() );
	}
	hh = spine[ i ];
    }
    boundary.push_back( hh->next() );
    boundary.push_back( hh->prev() );
    for ( int i = spine.size() - 2; i >= 0; --i ) {
	if ( hh->opposite()->prev() == spine[ i ] ) {
	    boundary.push_back( hh->opposite()->next() );
	}
	hh = spine[ i ];
    }

#ifdef MYDEBUG
    if ( boundary.size() < 10 ) {
	cerr << " Boundary : ";
	for ( unsigned int k = 0; k < boundary.size(); ++k )
	    cerr << " " << boundary[ k ]->vertex()->id() << ends;
	cerr << endl;
    }
#endif	// MYDEBUG
}


int countStitches( vector< Halfedge_handle > & bound )
{
    int nStitches = 0;

    for ( unsigned int i = 0; i < bound.size(); ++i ) {
	bound[ i ]->facet()->piece() = -100;
    }
    for ( unsigned int i = 0; i < bound.size(); ++i ) {
#ifdef MYDEBUG
	cerr << " i = " << i << " nCuts0 = " << bound[ i ]->vertex()->nCuts()
	     << " nCuts1 = " << bound[ i ]->opposite()->vertex()->nCuts()
	     << " piece1 = " << bound[ i ]->facet()->piece()
	     << " piece2 = " << bound[ i ]->opposite()->facet()->piece()
	     << endl;
#endif	// MYDEBUG
	if ( ( bound[ i ]->vertex()->nCuts() >= 3 ) && 
	     ( bound[ i ]->opposite()->vertex()->nCuts() >= 3 ) &&
	     ( bound[ i ]->facet()->piece() != bound[ i ]->opposite()->facet()->piece() ) ) {
	    nStitches++;
	}
    }
    for ( unsigned int i = 0; i < bound.size(); ++i ) {
	bound[ i ]->facet()->piece() = NO_INDEX;
    }

    return nStitches;
}


//------------------------------------------------------------------------------
//	To unfold one sequence of triangle strips
//------------------------------------------------------------------------------
void projectEachCycle( vector< Halfedge_handle > & cycleS,
		       vector< vector< Halfedge_handle > > & cycleM,
		       vector< vector< Triangle2 > > & stripM,
		       double limitCos )
{
    int nextDir = 0;		// A(side[0]): -1, B(side[1]): 1
    
    cycleM.clear();
    stripM.clear();
    vector< Halfedge_handle > cycleC;
    vector< Triangle2 > stripC;
    int m = 0;

    Point2  coordP, coordN;
    Vector2 baseV,  curV;

    // cycleS: the cycle index
    for ( unsigned int i = 0; i < cycleS.size(); ++i ) {
	
	Halfedge_handle hh = cycleS[ i ];

	Point3 coord   [ 3 ];
	double side    [ 3 ];
	double sinAngle, cosAngle;
	Point2 corner  [ 3 ];

	// The 3 vertices for a triangle
	coord[ 0 ] = hh->opposite()->vertex()->point();
	coord[ 1 ] = hh->vertex()->point();
	coord[ 2 ] = hh->next()->vertex()->point();

	// Length of each side: to calculate the length of each side
	for ( int k = 0; k < 3; ++k ) 
		// side[0]=a, side[1]=b, side[2]=c
	    side[ k ] = sqrt( ( coord[ (k+2)%3 ] - coord[ (k+1)%3 ] ).squared_length() );

	// Law of cosine
	cosAngle = ( SQUARE( side[ 1 ] ) + SQUARE( side[ 2 ] ) - SQUARE( side[ 0 ] ) ) / ( 2.0*side[ 1 ]*side[ 2 ] );
	cosAngle = MIN2( 1.0f, MAX2( -1.0f, cosAngle ) );
	sinAngle = sqrt( 1.0 - SQUARE( cosAngle ) );

	// Initial coordinate for the unfolded triangle strip
	corner[ 0 ] = CGAL::ORIGIN;
	corner[ 1 ] = corner[ 0 ] + Vector2( side[ 2 ], 0.0 );
	corner[ 2 ] = corner[ 0 ] + Vector2( side[ 1 ]*cosAngle, side[ 1 ]*sinAngle );

	// A triangle strip
	Triangle2 t( corner[ 0 ], corner[ 1 ], corner[ 2 ] );

	// cerr << "[Before] m = " << m << " triangle = " << t << endl;
	// Next new vertices in 2D unfolding
	if ( m != 0 ) {
	    Vector2 coordTrans;
	    Vector2 baseVec;
	    double baseLength;
	    double sinRot, cosRot;
	    if ( nextDir == -1 ) {
		// cerr << " Direction B " << endl;
		coordTrans = stripC[ m-1 ][ 2 ] - CGAL::ORIGIN;
		baseVec = stripC[ m-1 ][ 1 ] - stripC[ m-1 ][ 2 ];
	    }
	    else if ( nextDir == 1 ) {
		// cerr << " Direction C " << endl;
		coordTrans = stripC[ m-1 ][ 0 ] - CGAL::ORIGIN;
		baseVec = stripC[ m-1 ][ 2 ] - stripC[ m-1 ][ 0 ];
	    }
	    else assert( FALSE );

	    // cerr << " baseVec = " << baseVec << endl;
	    baseLength = sqrt( baseVec.squared_length() );
	    sinRot = baseVec[ 1 ] / baseLength;
	    cosRot = baseVec[ 0 ] / baseLength;
	    // cerr << " sinRot = " << sinRot << " cosRot = " << cosRot << endl;

	    Transformation2 rotate( CGAL::ROTATION, sinRot, cosRot );
	    Transformation2 translate( CGAL::TRANSLATION, coordTrans );
	    Transformation2 composite = translate * rotate;
	    t = t.transform( composite );
	    // cerr << " coordTrans = " << coordTrans << endl;
	    // cerr << " sinRot = " << sinRot << " cosRot = " << cosRot << endl;
	    // cerr << "[B] size = " << cycleC.size() << " stripC.size() = " << stripC.size() << endl;
	}
	// cerr << "[After ] m = " << m << " triangle = " << t << endl;
	// getchar();

	// To compute the inner product of first dual line with the current dual line
	if ( m == 0 ) {
	    ; 			// do nothing
	}
	else if ( m == 1 ) {
	    coordP = centerTriangle( stripC[ 0 ] );
	    coordN = centerTriangle( t );
	    // cerr << " m = " << m << " coordP = " << coordP << endl;
	    // cerr << " m = " << m << " coordN = " << coordN << endl;
	    baseV  = coordN - coordP;
	    // cerr << " m = " << m << " baseV = "  << baseV << endl;
	}
	else {
	    coordP = coordN;
	    coordN = centerTriangle( t );
	    curV   = coordN - coordP;

	    // cerr << " m = " << m << " coordP = " << coordP << endl;
	    // cerr << " m = " << m << " coordN = " << coordN << endl;
	    // cerr << " m = " << m << " baseV = " << baseV << endl;
	    // cerr << " m = " << m << " curV = " << curV << endl;

	    // cerr << " i = " << i << " m = " << m << " baseV * curV " << baseV * curV << endl;

	    // If the computed inner product is more than defined limitCos (60 or 90 degree), cut and split it into other cycle
	    if ( baseV * curV < limitCos ) {
		// #ifdef CUT_WITH_SELF_INTERSECTION
		cycleM.push_back( cycleC );
		stripM.push_back( stripC );
		cycleC.clear();
		stripC.clear();
		// cerr << " Cut out the strip " << endl;
		// getchar();
		m = 0;
		// #endif	// CUT_WITH_SELF_INTERSECTION
	    }
	}

	cycleC.push_back( hh );
	stripC.push_back( t );
	m++;

	// If we reach the end half edge of the array, leave the loop.
	if ( i == cycleS.size() - 1 ) break;

#ifdef MYDEBUG
	cerr << "cycle[ " << setw( 3 ) << i << " ] = " << cycleS[ i ]->id() 
	     << " : size = " << ( int )cycleS.size() << endl;
	cerr << "cycle[ " << setw( 3 ) << (i+1)%cycleS.size() << " ] = "
	     << cycleS[ (i+1)%cycleS.size() ]->id() 
	     << " : size = " << ( int )cycleS.size() << endl;
	cerr << " hh->id() = " << hh->id() << endl;
	cerr << " hh->opposite()->id() = " << hh->opposite()->id() << endl;
	cerr << " hh->prev()->id() = " << hh->prev()->id() << endl;
	cerr << " hh->prev()->opposite()->id() = " << hh->prev()->opposite()->id() << endl;
	cerr << " hh->next()->id() = " << hh->next()->id() << endl;
	cerr << " hh->next()->opposite()->id() = " << hh->next()->opposite()->id() << endl;
#endif	// MYDEBUG

	if ( cycleS[ (i+1)%cycleS.size() ] == hh->next()->opposite() ) {
	    nextDir = -1;
	}
	else if ( cycleS[ (i+1)%cycleS.size() ] == hh->prev()->opposite() ) {
	    nextDir = 1;
	}
	else assert( FALSE ); 

    }
    cycleM.push_back( cycleC );
    stripM.push_back( stripC );
    assert( ( int )cycleM.size() <= ( int )cycleS.size() );
}


//------------------------------------------------------------------------------
//	Translation for each cycle
//------------------------------------------------------------------------------
void translateEachCycle( vector< Facet_handle > & subpatch )
{
    int nPoints = 0;
    Vector2 sum( 0.0, 0.0 );
    Vector2 ave;

    for ( unsigned int i = 0; i < subpatch.size(); ++i ) {
	Halfedge_handle hv[ NUM_SIDES ];
	hv[ 0 ] = subpatch[ i ]->halfedge()->prev();
	hv[ 1 ] = subpatch[ i ]->halfedge();
	hv[ 2 ] = subpatch[ i ]->halfedge()->next();
	for ( unsigned int j = 0; j < NUM_SIDES; ++j ) {
	    if ( ! hv[ j ]->connect() ) {
#ifdef MYDEBUG
		cerr << " HalfEdge: " 
		     << hv[ j ]->opposite()->vertex()->id() << " == " 
		     << hv[ j ]->vertex()->id() << " is not connected " << endl;
		cerr << " Face ID : " << hv[ j ]->facet()->id() << endl;
#endif	// MYDEBUG
		sum = sum + ( subpatch[ i ]->triangle()[ j ] - CGAL::ORIGIN );
		nPoints++;
	    }
#ifdef MYDEBUG
	    else {
		cerr << " HalfEdge: " 
		     << hv[ j ]->opposite()->vertex()->id() << " == " 
		     << hv[ j ]->vertex()->id() << " is connected " << endl;
		cerr << " Face ID : " << hv[ j ]->facet()->id() << endl;
	    }
#endif	// MYDEBUG
	}
    }
    ave = sum / ( double )nPoints;
    if ( nPoints != ( int )( subpatch.size() + 2 ) ) {
	cerr << " nFaces = " << subpatch.size() << " nPoints = " << nPoints << endl;
	assert( nPoints == ( int )( subpatch.size() + 2 ) );
    }

    Transformation2 translate( CGAL::TRANSLATION, -ave );
    for ( unsigned int i = 0; i < subpatch.size(); ++i ) {
	subpatch[ i ]->triangle() = subpatch[ i ]->triangle().transform( translate );
    }
}

//------------------------------------------------------------------------------
//	Rotation for each cycle using covariance calculation
//------------------------------------------------------------------------------
void rotateEachCycle( vector< Facet_handle > & subpatch )
{
    /* Data has been centered already */
    /*
      X =       ( x_{11}, x_{12} )
                ( x_{21}, x_{22} )
                ...
                ( x_{n1}, x_{n2} )
                d = dim, n = num
    */
    unsigned int num	= subpatch.size() + 2;
    unsigned int dim	= 2;
    double * data	= new double [ num*dim ];

    int nPoints = 0;
    for ( unsigned int i = 0; i < subpatch.size(); ++i ) {
	Halfedge_handle hv[ NUM_SIDES ];
	hv[ 0 ] = subpatch[ i ]->halfedge()->prev();
	hv[ 1 ] = subpatch[ i ]->halfedge();
	hv[ 2 ] = subpatch[ i ]->halfedge()->next();
	for ( unsigned int j = 0; j < NUM_SIDES; ++j ) {
	    if ( ! hv[ j ]->connect() ) {
		data[ nPoints*dim + 0 ] = subpatch[ i ]->triangle()[ j ][ 0 ];
		data[ nPoints*dim + 1 ] = subpatch[ i ]->triangle()[ j ][ 1 ];
		nPoints++;
	    }
	}
    }

    assert( nPoints == ( int )num );

    /*****************************************
      analyze the engenstructure of X^T X
     *****************************************/

    /* define the matrix X */
    gsl_matrix_view X = gsl_matrix_view_array( data, num, dim );

    /* memory reallocation */
    // proj = ( double * )realloc( proj, sizeof( double ) * PDIM * num );

    /* calculate the covariance matrix B */
    gsl_matrix * B              = gsl_matrix_alloc( dim, dim );
    gsl_blas_dgemm( CblasTrans, CblasNoTrans, 1.0,
                    &(X.matrix),
                    &(X.matrix), 0.0,
                    B );
       
    /* divided by the number of samples */
    gsl_matrix_scale( B, 1.0/(double)num );

    gsl_vector * eigenVal       = gsl_vector_alloc( dim );
    gsl_matrix * eigenVec       = gsl_matrix_alloc( dim, dim );
    gsl_eigen_symmv_workspace * w
                                = gsl_eigen_symmv_alloc( dim );

    // eigenanalysis of the matrix B
    gsl_eigen_symmv( B, eigenVal, eigenVec, w );
    // release the memory of w
    gsl_eigen_symmv_free( w );
    // sort eigenvalues in a descending order
    gsl_eigen_symmv_sort( eigenVal, eigenVec, GSL_EIGEN_SORT_VAL_DESC );

    gsl_vector_view eachVec = gsl_matrix_column( eigenVec, 0 );
    double cosRot = gsl_matrix_get( eigenVec, 0, 1 );
    double sinRot = gsl_matrix_get( eigenVec, 1, 1 );
#ifdef DEBUG
    cerr << " 2nd axis : " << cosRot << " , " << sinRot << endl;
#endif	// DEBUG

    Transformation2 rotate( CGAL::ROTATION, -sinRot, cosRot );
    for ( unsigned int i = 0; i < subpatch.size(); ++i ) {
	subpatch[ i ]->triangle() = subpatch[ i ]->triangle().transform( rotate );
    }
}

//------------------------------------------------------------------------------
//	Transformation for each patch
//------------------------------------------------------------------------------
void transformEachCycle( vector< Facet_handle > & subpatch )
{
    translateEachCycle	( subpatch );
    rotateEachCycle	( subpatch );
}

//------------------------------------------------------------------------------
//	Bouding box drawing
//------------------------------------------------------------------------------
void boundEachCycle( vector< Facet_handle > & subpatch, Bbox2 & box )
{
    list< Point2 > corner;
    for ( unsigned int i = 0; i < subpatch.size(); ++i )
	for ( unsigned int j = 0; j < NUM_SIDES; ++j ) {
	    corner.push_back( subpatch[ i ]->triangle()[ j ] );
	    // cerr << " i = " << i << " j = " << j << " corner = " << subpatch[ i ]->triangle()[ j ] << endl;
	}
    Rectangle2 rect = CGAL::bounding_box( corner.begin(), corner.end() );

    box = rect.bbox();
    // cerr << " Bounding box = " << box << endl;
}


void resetCycleIDs( vector< vector< Facet_handle > > & patch )
{
    // Assign the strip ID to the sequence of halfedges and faces
    for ( unsigned int i = 0; i < patch.size(); ++i ) {
	for ( unsigned int j = 0; j < patch[ i ].size(); ++j ) {
	    patch[ i ][ j ]->piece() = i;
	    Halfedge_handle hv[ 3 ];
	    hv[ 0 ] = patch[ i ][ j ]->halfedge()->prev();
	    hv[ 1 ] = patch[ i ][ j ]->halfedge();
	    hv[ 2 ] = patch[ i ][ j ]->halfedge()->next();
	    for ( unsigned int k = 0; k < NUM_SIDES; ++k ) {
		Facet_handle fhI = hv[ k ]->facet();
		assert( fhI == patch[ i ][ j ] );
		Facet_handle fhO = hv[ k ]->opposite()->facet();
		if ( ( fhI->piece() == fhO->piece() ) && hv[ k ]->connect() ) {
		    hv[ k ]->cycle() = i;
		    hv[ k ]->opposite()->cycle() = i;
		}
	    }
	}
    }
}


//------------------------------------------------------------------------------
//	Adjust the unfolded cycles properly
//------------------------------------------------------------------------------
void boundCycles( vector< int > & cycleID,
		  vector< vector< Facet_handle > > & patch,
		  vector< Bbox2 > & bound )
{
    multimap< double, int > sortMap;

    cycleID.clear();
    bound.clear();

    // Transform each strip for better 2D arrangement of the unfolded patterns
    for ( unsigned int i = 0; i < patch.size(); ++i ) {
#ifdef MYDEBUG
	cerr << " i = " << i << " patch[ i ].size() = " << patch[ i ].size() << endl;
#endif	// MYDEBUG
	transformEachCycle( patch[ i ] ); // transformation for each patch
	Bbox2 localbox;
	boundEachCycle( patch[ i ], localbox ); // bounding box drawing
	bound.push_back( localbox );
	sortMap.insert( make_pair< double, int >( localbox.ymax() - localbox.ymin(), i ) );
    }
    
    multimap< double, int >::iterator it = sortMap.end();
    do {
	it--;
	cycleID.push_back( it->second );
    } while ( it != sortMap.begin() );
    sortMap.clear();
}

//------------------------------------------------------------------------------
//	Layout the unfolded patches (with bounding box) properly on the diplay window
//------------------------------------------------------------------------------
void layoutSingle( vector< int > & cycleID,
		   vector< vector< Facet_handle > > & patch,
		   Bbox2 & sheet, vector< Bbox2 > & bound )
{
    int isFit;
    const double ratio = MARGIN_RATIO;
    stack< Point2 > cornerP, cornerN;
    double smallgap = ratio * ( sheet.xmax() - sheet.xmin() );

#ifdef MYDEBUG
    for ( unsigned int k = 0; k < bound.size(); ++k )
	cerr << " bound[ " << k << " ] = " << bound[ k ] << endl;
#endif	// MYDEBUG

    do {
	isFit = true;
	int dir = 1;		// 1 for left, -1 for right

	while ( ! cornerN.empty() ) cornerN.pop();
	// initialization
	Point2 topRight( sheet.xmax(), sheet.ymax() );
	Point2 topLeft ( sheet.xmin(), sheet.ymax() );
	cornerN.push( topRight );
	cornerN.push( topLeft  );

	cornerP = cornerN;

	double xLimit = cornerP.top().x() + smallgap;
	double yLimit = cornerP.top().y() - smallgap;
	double yLimitNext = yLimit;
	cornerP.pop();
	
	for ( unsigned int k = 0; k < cycleID.size(); ++k ) {
	    int id = cycleID[ k ];

	    double width	= bound[ id ].xmax() - bound[ id ].xmin();
	    double height	= bound[ id ].ymax() - bound[ id ].ymin();
	    Vector2 shift;
	    Bbox2 curBound;

#ifdef MYDEBUG
	    cerr << "************************************************************" << endl;
	    cerr << " Handling Cycle ID: " << id << endl;
	    cerr << "************************************************************" << endl;
	    cerr << " xLimit = " << xLimit << " yLimit = " << yLimit << endl;
	    cerr << " width = " << width << " height = " << height << endl;
	    cerr << " cornerP.top() = " << cornerP.top() << endl;
#endif	// MYDEBUG

	    while ( true ) {
		// Check the bottom boundary:
		// If the unfoled patch cannot fit into the given sheet size, go
		// out of this loop with the flag "failure"
		if ( ( yLimit - height < sheet.ymin() ) ||
		    // This is very important!!
		    // added on 2010/09/17
		     ( width > sheet.xmax() - sheet.xmin() - 2.0*smallgap ) ||
		     ( height > sheet.ymax() - sheet.ymin() - 2.0*smallgap ) ) {
#ifdef MYDEBUG
		    cerr << " Case A : " << endl;
		    cerr << " Enlarge the sheet and try again" << endl;
#endif	// MYDEBUG
		    isFit = false;
		    while ( ! cornerP.empty() ) cornerP.pop();
		    while ( ! cornerN.empty() ) cornerN.pop();
		    break;
		}

#ifdef MYDEBUG
		cerr << " xLimit = " << xLimit << " yLimit = " << yLimit << endl;
		cerr << " width = " << width << " height = " << height << endl;
		cerr << " cornerP.top() = " << cornerP.top() << endl;
#endif	// MYDEBUG
		    
		// Check the right boundary:
		// If the unfolded patch exceeds the right boundary of the
		// sheet, then change the orientation of the patch layout.
		if ( ( dir == 1 ) && ( xLimit + width > sheet.xmax() ) ) {
#ifdef MYDEBUG
		    cerr << " Case B : " << endl;
#endif	// MYDEBUG
		    dir = -1;
		    xLimit = sheet.xmax() - smallgap;
		    cornerP = cornerN;
		    while ( ! cornerN.empty() ) cornerN.pop();
		    yLimit = yLimitNext = cornerP.top().y() - smallgap;

		    Point2 mark( sheet.xmax(), cornerP.top().y() );
		    cornerN.push( mark );
#ifdef MYDEBUG
		    cerr << " mark = " << mark << " is pushed. " << endl;
#endif	// MYDEBUG

		    cornerP.pop();
#ifdef MYDEBUG
		    cerr << " Change direction into Left" << endl;
		    cerr << " xLimit = " << xLimit << " yLimit = " << yLimit << endl;
#endif	// MYDEBUG
		    continue;
		}
		
		// Check the left boundary
		// If the unfolded patch exceeds the left boundary of the
		// sheet, then change the orientation of the patch layout.
		else if ( ( dir == -1 ) && ( xLimit - width < sheet.xmin() ) ) {
#ifdef MYDEBUG
		    cerr << " Case C : " << endl;
#endif	// MYDEBUG
		    dir = 1;
		    xLimit = sheet.xmin() + smallgap;
		    cornerP = cornerN;
		    while ( ! cornerN.empty() ) cornerN.pop();
		    yLimit = yLimitNext = cornerP.top().y() - smallgap;

		    Point2 mark( sheet.xmin(), cornerP.top().y() );
		    cornerN.push( mark );
#ifdef MYDEBUG
		    cerr << " mark = " << mark << " is pushed. " << endl;
#endif	// MYDEBUG

		    cornerP.pop();
#ifdef MYDEBUG
		    cerr << " Change direction into Right" << endl;
		    cerr << " xLimit = " << xLimit << " yLimit = " << yLimit << endl;
#endif	// MYDEBUG
		    continue;
		}
		
		else if ( ( dir == 1 ) && ( xLimit + width > cornerP.top().x() ) ) {
#ifdef MYDEBUG
		    cerr << " Case D : " << endl;
		    cerr << " xLimit + width = " << xLimit + width << endl;
		    cerr << " cornerP.top() = " << cornerP.top() << endl;
#endif	// MYDEBUG
		    yLimitNext = cornerP.top().y() - smallgap;
		    yLimit = MIN2( yLimit, yLimitNext );
		    cornerP.pop();
		    continue;
		}

		else if ( ( dir == -1 ) && ( xLimit - width < cornerP.top().x() ) ) {
#ifdef MYDEBUG
		    cerr << " Case E : " << endl;
		    cerr << " xLimit - width = " << xLimit - width << endl;
		    cerr << " cornerP.top() = " << cornerP.top() << endl;
#endif	// MYDEBUG
		    yLimitNext = cornerP.top().y() - smallgap;
		    yLimit = MIN2( yLimit, yLimitNext );
		    cornerP.pop();
		    continue;
		}

#ifdef MYDEBUG
		cerr << " Main Case : " << endl;
#endif	// MYDEBUG

		Point2 center( 0.5*(bound[id].xmin()+bound[id].xmax()), 0.5*(bound[id].ymin()+bound[id].ymax()) );
		Point2 target( xLimit + 0.5*dir*width, yLimit - 0.5*height );
		shift = target - center;
		
		// Check the possible overlaps
		// The overlaps here are expected to be already eliminated.
		bool isOverlapped = false;
		curBound = Bbox2( bound[ id ].xmin() + shift.x(), bound[ id ].ymin() + shift.y(),
				  bound[ id ].xmax() + shift.x(), bound[ id ].ymax() + shift.y() );
		// If the layout process goes in the right direction
		if ( dir == 1 ) {
		    Point2 mark( bound[ id ].xmax() + shift.x(), bound[ id ].ymin() + shift.y() );
#ifdef MYDEBUG
		    cerr << " mark = " << mark << " is pushed. " << endl;
#endif	// MYDEBUG
		    cornerN.push( mark );
		}
		else if ( dir == -1 ) {
		    Point2 mark( bound[ id ].xmin() + shift.x(), bound[ id ].ymin() + shift.y() );
#ifdef MYDEBUG
		    cerr << " mark = " << mark << " is pushed. " << endl;
#endif	// MYDEBUG
		    cornerN.push( mark );
		}

		yLimit = yLimitNext;

		for ( unsigned int m = 0; m < k; ++m ) {
		    if ( do_overlap( bound[ cycleID[ m ] ], curBound ) ) {
			isOverlapped = true;
			cerr << " bound[ " << cycleID[ m ] << " ] : " << bound[ cycleID[ m ] ] << endl;
			cerr << " curbound : " << curBound << endl;
			assert( false );
		    }
		}
		break;
	    }

	    // If the unfolded patch cannot fit into he given sheet size,
	    // enlarge the sheet and try the layout process from the
	    // begining again.
	    if ( !isFit ) {
		double curSize = sheet.xmax() - sheet.xmin();
//------------------------------------------------------------------------------
//		This part is specific for a signle patch
//------------------------------------------------------------------------------
		double nextMargin = ratio;
		double nextSize = smallgap;
		double nextWidth = 0.0;
		unsigned int m = 0;
		do {
		    double eachWidth = bound[ cycleID[ m ] ].xmax() - bound[ cycleID[ m ] ].xmin();
		    nextSize += eachWidth;
		    nextWidth += eachWidth;
		    nextSize += smallgap;
		    nextMargin += ratio;
		    m++;
		} while ( ( nextSize < curSize + 1.0e-8 ) && ( m < bound.size() ) );
		double paperWidth = nextWidth / ( 1.0 - nextMargin ) + 1.0e-8; 
		double paperHeight = SQRT2 * paperWidth;
		double initHeight = bound[ cycleID[ 0 ] ].ymax() - bound[ cycleID[ 0 ] ].ymin();
		initHeight /= ( 1.0 - 2 * ratio );
#ifdef MYDEBUG
		cerr << "[UNFIT]" << endl;
		cerr << " paperWidth = " << paperWidth << " paperHeight = " << paperHeight;
		cerr << " initHeight = " << initHeight << endl;
#endif	// MYDEBUG
		if ( initHeight > paperHeight ) {
		    sheet = Bbox2( -( 0.5 * initHeight / SQRT2 ),
				   -  0.5 * initHeight,
				    ( 0.5 * initHeight / SQRT2 ),
				      0.5 * initHeight );
		}
		else {
		    sheet = Bbox2( -( 0.5 * paperWidth ),
				   -( 0.5 * paperWidth ) * SQRT2,
				    ( 0.5 * paperWidth ),
				    ( 0.5 * paperWidth ) * SQRT2 );
		}
		smallgap = ratio * ( sheet.xmax() - sheet.xmin() );
		break;
	    }
	    else {
#ifdef MYDEBUG
		cerr << "[FIT!!]" << endl;
		cerr << " Shift = " << shift << " id = " << id << " patch.size() = " << patch[ id ].size() << endl;
#endif	// MYDEBUG
		Transformation2 translate( CGAL::TRANSLATION, shift );
		for ( unsigned int j = 0; j < patch[ id ].size(); ++j )
		    patch[ id ][ j ]->triangle() = 
			patch[ id ][ j ]->triangle().transform( translate );
		bound[ id ] = curBound;
		xLimit = xLimit + dir*(width+smallgap);
		// if ( yLimitN > yLimitP - height ) yLimitN = yLimitP - height;
		// yLimit = yLimit - height;
	    }
	} 
    } while ( !isFit );


#ifdef DEBUG
    cerr << "##### Final bounding box ##### : "
	 << " xmin = " << sheet.xmin()
	 << " ymin = " << sheet.ymin()
	 << " xmax = " << sheet.xmax()
	 << " ymax = " << sheet.ymax() << endl;
#endif	// DEBUG
}


void layoutMultiple( vector< int > & cycleID,
		     vector< vector< Facet_handle > > & patch,
		     Bbox2 & sheet, vector< Bbox2 > & bound )
{

    // cerr << "########## Now in layoutMultiple ##########" << endl;

    int isFit;
    const double ratio = MARGIN_RATIO;
    stack< Point2 > cornerP, cornerN;
    double smallgap = ratio * ( sheet.xmax() - sheet.xmin() );

#ifdef MYDEBUG
    for ( unsigned int k = 0; k < bound.size(); ++k )
	cerr << " bound[ " << k << " ] = " << bound[ k ] << endl;
#endif	// MYDEBUG

    do {
	isFit = true;
	int dir = 1;		// 1 for left, -1 for right

	while ( ! cornerN.empty() ) cornerN.pop();
	// initialization
	Point2 topRight( sheet.xmax(), sheet.ymax() );
	Point2 topLeft ( sheet.xmin(), sheet.ymax() );
	cornerN.push( topRight );
	cornerN.push( topLeft  );

	cornerP = cornerN;

	double xLimit = cornerP.top().x() + smallgap;
	double yLimit = cornerP.top().y() - smallgap;
	double yLimitNext = yLimit;
	cornerP.pop();
	
	for ( unsigned int k = 0; k < cycleID.size(); ++k ) {
	    int id = cycleID[ k ];

	    double width	= bound[ id ].xmax() - bound[ id ].xmin();
	    double height	= bound[ id ].ymax() - bound[ id ].ymin();
	    Vector2 shift;
	    Bbox2 curBound;

#ifdef MYDEBUG
	    cerr << "************************************************************" << endl;
	    cerr << " Handling Cycle ID: " << id << endl;
	    cerr << "************************************************************" << endl;
	    cerr << " xLimit = " << xLimit << " yLimit = " << yLimit << endl;
	    cerr << " width = " << width << " height = " << height << endl;
	    cerr << " cornerP.top() = " << cornerP.top() << endl;
#endif	// MYDEBUG

	    while ( true ) {
		// Check the bottom boundary:
		// If the unfoled patch cannot fit into the given sheet size, go
		// out of this loop with the flag "failure"
		if ( ( yLimit - height < sheet.ymin() ) ||
		    // This is very important!!
		    // added on 2010/09/17
		     ( width > sheet.xmax() - sheet.xmin() - 2.0*smallgap ) ||
		     ( height > sheet.ymax() - sheet.ymin() - 2.0*smallgap ) ) {
#ifdef MYDEBUG
		    cerr << " Case A : " << endl;
		    cerr << " Enlarge the sheet and try again" << endl;
#endif	// MYDEBUG
		    isFit = false;
		    while ( ! cornerP.empty() ) cornerP.pop();
		    while ( ! cornerN.empty() ) cornerN.pop();
		    break;
		}

#ifdef MYDEBUG
		cerr << " xLimit = " << xLimit << " yLimit = " << yLimit << endl;
		cerr << " width = " << width << " height = " << height << endl;
		cerr << " cornerP.top() = " << cornerP.top() << endl;
#endif	// MYDEBUG
		    
		// Check the right boundary:
		// If the unfolded patch exceeds the right boundary of the
		// sheet, then change the orientation of the patch layout.
		if ( ( dir == 1 ) && ( xLimit + width > sheet.xmax() ) ) {
#ifdef MYDEBUG
		    cerr << " Case B : " << endl;
#endif	// MYDEBUG
		    dir = -1;
		    xLimit = sheet.xmax() - smallgap;
		    cornerP = cornerN;
		    while ( ! cornerN.empty() ) cornerN.pop();
		    yLimit = yLimitNext = cornerP.top().y() - smallgap;

		    Point2 mark( sheet.xmax(), cornerP.top().y() );
		    cornerN.push( mark );
#ifdef MYDEBUG
		    cerr << " mark = " << mark << " is pushed. " << endl;
#endif	// MYDEBUG

		    cornerP.pop();
#ifdef MYDEBUG
		    cerr << " Change direction into Left" << endl;
		    cerr << " xLimit = " << xLimit << " yLimit = " << yLimit << endl;
#endif	// MYDEBUG
		    continue;
		}
		
		// Check the left boundary
		// If the unfolded patch exceeds the left boundary of the
		// sheet, then change the orientation of the patch layout.
		else if ( ( dir == -1 ) && ( xLimit - width < sheet.xmin() ) ) {
#ifdef MYDEBUG
		    cerr << " Case C : " << endl;
#endif	// MYDEBUG
		    dir = 1;
		    xLimit = sheet.xmin() + smallgap;
		    cornerP = cornerN;
		    while ( ! cornerN.empty() ) cornerN.pop();
		    yLimit = yLimitNext = cornerP.top().y() - smallgap;

		    Point2 mark( sheet.xmin(), cornerP.top().y() );
		    cornerN.push( mark );
#ifdef MYDEBUG
		    cerr << " mark = " << mark << " is pushed. " << endl;
#endif	// MYDEBUG

		    cornerP.pop();
#ifdef MYDEBUG
		    cerr << " Change direction into Right" << endl;
		    cerr << " xLimit = " << xLimit << " yLimit = " << yLimit << endl;
#endif	// MYDEBUG
		    continue;
		}
		
		else if ( ( dir == 1 ) && ( xLimit + width > cornerP.top().x() ) ) {
#ifdef MYDEBUG
		    cerr << " Case D : " << endl;
		    cerr << " xLimit + width = " << xLimit + width << endl;
		    cerr << " cornerP.top() = " << cornerP.top() << endl;
#endif	// MYDEBUG
		    yLimitNext = cornerP.top().y() - smallgap;
		    yLimit = MIN2( yLimit, yLimitNext );
		    cornerP.pop();
		    continue;
		}

		else if ( ( dir == -1 ) && ( xLimit - width < cornerP.top().x() ) ) {
#ifdef MYDEBUG
		    cerr << " Case E : " << endl;
		    cerr << " xLimit - width = " << xLimit - width << endl;
		    cerr << " cornerP.top() = " << cornerP.top() << endl;
#endif	// MYDEBUG
		    yLimitNext = cornerP.top().y() - smallgap;
		    yLimit = MIN2( yLimit, yLimitNext );
		    cornerP.pop();
		    continue;
		}

#ifdef MYDEBUG
		cerr << " Main Case : " << endl;
#endif	// MYDEBUG

		Point2 center( 0.5*(bound[id].xmin()+bound[id].xmax()), 0.5*(bound[id].ymin()+bound[id].ymax()) );
		Point2 target( xLimit + 0.5*dir*width, yLimit - 0.5*height );
		shift = target - center;
		
		// Check the possible overlaps
		// The overlaps here are expected to be already eliminated.
		bool isOverlapped = false;
		curBound = Bbox2( bound[ id ].xmin() + shift.x(), bound[ id ].ymin() + shift.y(),
				  bound[ id ].xmax() + shift.x(), bound[ id ].ymax() + shift.y() );
		// If the layout process goes in the right direction
		if ( dir == 1 ) {
		    Point2 mark( bound[ id ].xmax() + shift.x(), bound[ id ].ymin() + shift.y() );
#ifdef MYDEBUG
		    cerr << " mark = " << mark << " is pushed. " << endl;
#endif	// MYDEBUG
		    cornerN.push( mark );
		}
		else if ( dir == -1 ) {
		    Point2 mark( bound[ id ].xmin() + shift.x(), bound[ id ].ymin() + shift.y() );
#ifdef MYDEBUG
		    cerr << " mark = " << mark << " is pushed. " << endl;
#endif	// MYDEBUG
		    cornerN.push( mark );
		}

		yLimit = yLimitNext;

		for ( unsigned int m = 0; m < k; ++m ) {
		    if ( do_overlap( bound[ cycleID[ m ] ], curBound ) ) {
			isOverlapped = true;
			cerr << " bound[ " << cycleID[ m ] << " ] : " << bound[ cycleID[ m ] ] << endl;
			cerr << " curbound : " << curBound << endl;
			assert( false );
		    }
		}
		break;
	    }

	    // If the unfolded patch cannot fit into he given sheet size,
	    // enlarge the sheet and try the layout process from the
	    // begining again.
	    if ( !isFit ) {
		double curSize = sheet.xmax() - sheet.xmin();
//------------------------------------------------------------------------------
//		This part is specific for multiple patches
//------------------------------------------------------------------------------
		double nextMargin = ratio;
		double nextSize = curSize + EXTEND_PAPER_WIDTH;
		double nextWidth = nextSize;
//------------------------------------------------------------------------------
		double paperWidth = nextWidth / ( 1.0 - nextMargin ) + 1.0e-8; 
		double paperHeight = SQRT2 * paperWidth;
		double initHeight = bound[ cycleID[ 0 ] ].ymax() - bound[ cycleID[ 0 ] ].ymin();
		initHeight /= ( 1.0 - 2 * ratio );
#ifdef MYDEBUG
		cerr << " paperWidth = " << paperWidth << " paperHeight = " << paperHeight;
		cerr << " initHeight = " << initHeight << endl;
#endif	// MYDEBUG
		if ( initHeight > paperHeight ) {
		    sheet = Bbox2( -( 0.5 * initHeight / SQRT2 ),
				   -  0.5 * initHeight,
				    ( 0.5 * initHeight / SQRT2 ),
				      0.5 * initHeight );
		}
		else {
		    sheet = Bbox2( -( 0.5 * paperWidth ),
				   -( 0.5 * paperWidth ) * SQRT2,
				    ( 0.5 * paperWidth ),
				    ( 0.5 * paperWidth ) * SQRT2 );
		}
		smallgap = ratio * ( sheet.xmax() - sheet.xmin() );
		break;
	    }
	    else {
#ifdef MYDEBUG
		cerr << " Shift = " << shift << " id = " << id << " patch.size() = " << patch[ id ].size() << endl;
#endif	// MYDEBUG
		Transformation2 translate( CGAL::TRANSLATION, shift );
		for ( unsigned int j = 0; j < patch[ id ].size(); ++j )
		    patch[ id ][ j ]->triangle() = 
			patch[ id ][ j ]->triangle().transform( translate );
		bound[ id ] = curBound;
		xLimit = xLimit + dir*(width+smallgap);
		// if ( yLimitN > yLimitP - height ) yLimitN = yLimitP - height;
		// yLimit = yLimit - height;
	    }
	} 
    } while ( !isFit );


#ifdef DEBUG
    cerr << "##### Final bounding box ##### : "
	 << " xmin = " << sheet.xmin()
	 << " ymin = " << sheet.ymin()
	 << " xmax = " << sheet.xmax()
	 << " ymax = " << sheet.ymax() << endl;
#endif	// DEBUG
}


//------------------------------------------------------------------------------
//	arrange unfolded patches from their 2D projections
//------------------------------------------------------------------------------
void arrangeSingle( vector< vector< Facet_handle > > & patch,
		    Bbox2 & sheet, vector< Bbox2 > & bound )
{
    vector< int > cycleID;
    boundCycles( cycleID, patch, bound );
    sheet = Bbox2( -DEFAULT_PAPER_WIDTH, -DEFAULT_PAPER_HEIGHT, 
		   DEFAULT_PAPER_WIDTH,  DEFAULT_PAPER_HEIGHT );
    layoutSingle( cycleID, patch, sheet, bound );
}


void arrangeMultiple( vector< vector< Facet_handle > > & patch,
		      Bbox2 & sheet, vector< Bbox2 > & bound )
{
    vector< int > cycleID;
    boundCycles( cycleID, patch, bound );
    sheet = Bbox2( -DEFAULT_PAPER_WIDTH, -DEFAULT_PAPER_HEIGHT, 
		   DEFAULT_PAPER_WIDTH,  DEFAULT_PAPER_HEIGHT );
    layoutMultiple( cycleID, patch, sheet, bound );
}


//------------------------------------------------------------------------------
//	The 2D unfoldings for all the sequence of triangle strips
//------------------------------------------------------------------------------
void projectCycles( vector< vector< Halfedge_handle > > & spine,
		    vector< vector< Facet_handle > > & patch,
		    Bbox2 & sheet, vector< Bbox2 > & bound,
		    double limitCos )
{
    vector< vector< Halfedge_handle > > temp = spine;
    spine.clear();
    patch.clear();

    // temp.size(): the number of disjoint cycles
    for ( unsigned int i = 0; i < temp.size(); ++i ) {
#ifdef MYDEBUG
	cerr << " Handling the original cycle No. " << i << " Number of triangles = " << temp[ i ].size() << endl;
#endif	// MYDEBUG
	vector< vector< Halfedge_handle > > splitCycle; // the split sub-cycles
	vector< vector< Triangle2 > > splitStrip; // the split sub-strips
	// projectEachCycle( temp[ i ], splitCycle, splitStrip );
	projectEachCycle( temp[ i ], splitCycle, splitStrip, limitCos ); // to unfold one sequence of triangle strips
#ifdef MYDEBUG
	cerr << "******************************" << endl;
	cerr << " Cycle No. " << i << " has " << temp[ i ].size() << " triangles in total." << endl;
	cerr << " Cycle No. " << i << " is split into " << splitCycle.size() << " subcycles." << endl;
#endif	// MYDEBUG
	assert( splitCycle.size() == splitStrip.size() );
	for ( unsigned k = 0; k < splitCycle.size(); ++k ) {
	    vector< Facet_handle > splitPatch;
#ifdef MYDEBUG
	    cerr << "==> Subcycle No. " << k << " has " << splitStrip[ k ].size() << " triangles." << endl;
#endif	// MYDEBUG
	    // Increment the number of edge cuts
	    splitCycle[ k ][ 0 ]->vertex()->nCuts()++;
	    splitCycle[ k ][ 0 ]->opposite()->vertex()->nCuts()++;
	    for ( unsigned int m = 0; m < splitStrip[ k ].size(); ++m ) {
		splitCycle[ k ][ m ]->facet()->triangle() = splitStrip[ k ][ m ];
		rewindTriangle( splitCycle[ k ][ m ], splitCycle[ k ][ m ]->facet()->halfedge() );
		splitPatch.push_back( splitCycle[ k ][ m ]->facet() );
	    }
	    spine.push_back( splitCycle[ k ] );
	    patch.push_back( splitPatch );

#ifdef MYDEBUG
	    cerr << " Removing the edge : " 
		 << splitCycle[ k ][ 0 ]->facet()->id() << " == "
		 << splitCycle[ k ][ 0 ]->opposite()->facet()->id() << endl;
#endif	// MYDEBUG

	    splitCycle[ k ][ 0 ]->connect() = splitCycle[ k ][ 0 ]->opposite()->connect() = false;
	}
    }

    temp.clear();

    resetCycleIDs ( patch );
    cerr << "Number of unfolded patterns = " << spine.size() << endl;
    arrangeSingle( patch, sheet, bound );
}

