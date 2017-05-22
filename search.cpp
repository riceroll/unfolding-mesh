//------------------------------------------------------------------------------
//
//	search.cpp
//
//------------------------------------------------------------------------------
#include "common.h"

// #define USE_QUEUE
#define USE_STACK

#include <iostream>
#include <iomanip>

#ifdef USE_QUEUE
#include <queue>
#endif	// USE_QUEUE

#ifdef USE_STACK
#include <stack>
#endif	// USE_STACK



Vertex_handle sharpestVertex( Polyhedron & poly, vector< bool > & visit )
{
    double minAngle = 8.0 * M_PI;
    Vertex_handle ret;
    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	if ( ! visit[ vi->id() ] ) {
	    double curAngle = centralAngle( vi );
	    if ( curAngle < minAngle ) {
		minAngle = curAngle;
		ret = vi;
	    }
	}
    }

    return ret;
}


Halfedge_handle nextStart( Polyhedron & poly, vector< bool > & visit )
{
    Vertex_handle vhS = sharpestVertex( poly, visit );
    cerr << " Sharpest vertex ID : " << vhS->id() << endl;
    Halfedge_handle ret;
    Halfedge_vertex_circulator hvc = vhS->vertex_begin();
    do {
	if ( ! hvc->connect() ) {
	    ret = hvc->opposite();
	    return ret;
	}
    } while ( ++hvc != vhS->vertex_begin() );
    assert( false );
}


void nextRoute( const Halfedge_handle & base, vector< Halfedge_handle > & route )
{
    route.clear();
    Halfedge_handle cur = base;
    do {
	cur = cur->opposite()->next();
	// cerr << " Cheking : "; printHalfedge( cur );
	// if ( cur->connect() ) cerr << " Connect " << endl;
	// if ( cur->visit() ) cerr << " Visit " << endl;
	if ( ( ! cur->connect() ) && ( ! cur->visit() ) ) {
#ifdef MYDEBUG
	    cerr << " registered seam : ";
	    printHalfedge( cur );
#endif	// MYDEBUG
	    route.push_back( cur );
	}
    } while ( cur != base );
}


void depthFirst( const Halfedge_handle & hhS, 
		 const Halfedge_handle & hhC,
		 vector< Halfedge_handle > & seam,
		 vector< bool > & markV,
		 Bbox2 & sheet,
		 vector< vector< Facet_handle > > & patch,
		 vector< Bbox2 > & bound,
		 Polyhedron & poly, 
		 void (*redraw)( void ) )
{
    // static int count = 0;
    // const int limit = 0;
    // unsigned int k;

    if ( hhC == hhS ) return;

    markV[ hhC->vertex()->id() ] = true;
    hhC->visit() = true;
    seam.push_back( hhC );

    vector< Halfedge_handle > cand;
    nextRoute( hhC->opposite(), cand );

    if ( cand.size() == 0 ) {
	depthFirst( hhS, hhC->opposite(), seam, markV, sheet, patch, bound, poly, redraw );
    }
    else {
	unsigned int k = 0;
	while ( k < cand.size() - 1 ) {
	    vector< Triangle2 > small;
#ifdef MYDEBUG
	    cerr << " try seam : ";
	    printHalfedge( cand[ k ] );
#endif	// MYDEBUG
	    if ( canEmbed( cand[ k ], patch, small ) ) {
		stitch( cand[ k ], patch, small );
		
		resetCycleIDs( patch );
		arrangeSingle( patch, sheet, bound );
#ifdef MYDEBUG
		cerr << " Stitched via seam : ";
		printHalfedge( cand[ k ] );
#endif	// MYDEBUG
		
#ifdef DEBUG
		cerr << " Seam size = " << seam.size() << endl;
		for ( unsigned int i = 0; i < seam.size(); ++i ) {
		    cerr << " Seam[ " << setw( 3 ) << i << " ] : ";
		    printHalfedge( seam[ i ] );
		}
#endif	// DEBUG
		
		++k;
	    }
	    else {
#ifdef MYDEBUG
		cerr << "\a\a\a";
		cerr << "%%%%%%%%%%%%%%%%%%%% Cannot Stitched seam : ";
		printHalfedge( cand[ k ] );
#endif	// MYDEBUG
		
		break;
	    }
	}

#ifdef MYDEBUG
	cerr << " Pushed seam : ";
	printHalfedge( cand[ k ] );
#endif	// MYDEBUG
	depthFirst( hhS, cand[ k ], seam, markV, sheet, patch, bound, poly, redraw );
    }
}


void traverse( Outline & cur, 
#ifdef USE_QUEUE
	       queue< Outline > & next,
#endif	// USE_QUEUE
#ifdef USE_STACK
	       stack< Outline > & next,
#endif	// USE_STACK
	       vector< bool > & mark,
	       vector< vector< Facet_handle > > & patch,
	       Polyhedron & poly )
{
    // static int count = 0;
    // const int limit = 0;
    Halfedge_handle hhC = cur.start();
    vector< Triangle2 > small;
    
    // preprocessing: marking visited vertex 
    for ( unsigned int i = 0; i < cur.seam().size(); ++i ) {
        mark[ cur.seam()[ i ]->opposite()->vertex()->id() ] = true;
	mark[ cur.seam()[ i ]->vertex()->id() ] = true;
    }
    mark[ hhC->opposite()->vertex()->id() ] = true;
    // preprocessing: stitching edges
    for ( unsigned int i = 0; i < cur.tape().size(); ++i ) {
#ifdef MYDEBUG
	cerr << " Connecting the edge : ";
	printHalfedge( cur.tape()[ i ] );
#endif	// MYDEBUG
	if ( canEmbed( cur.tape()[ i ], patch, small ) ) {
	    stitch( cur.tape()[ i ], patch, small );
	    resetCycleIDs( patch );
	}
	else assert( false );
    }
    // preprocessing: visited edges
    for ( unsigned int i = 0; i < cur.seam().size(); ++i ) {
	cur.seam()[ i ]->visit() = true;
    }

    while ( ! hhC->visit() ) {
    // while ( hhC != cur.target() ) {
    // do {
#ifdef MYDEBUG
	cerr << " Current seam : "; printHalfedge( hhC );
	cerr << " Final seam : "; printHalfedge( cur.target() );
	// getchar();
#endif	// MYDEBUG
	
	mark[ hhC->vertex()->id() ] = true;
	hhC->visit() = true;
	cur.seam().push_back( hhC );

	vector< Halfedge_handle > cand;
	nextRoute( hhC->opposite(), cand );

	if ( cand.size() == 0 ) {
	    hhC = hhC->opposite();
	}
	else {
	    unsigned int k = 0;
	    while ( k < cand.size() - 1 ) {
#ifdef MYDEBUG
		cerr << " try seam : ";
		printHalfedge( cand[ k ] );
#endif	// MYDEBUG
		if ( canEmbed( cand[ k ], patch, small ) ) {

		    Outline remainder;
		    remainder.start() = cand[ k ];
		    remainder.seam() = cur.seam();
		    remainder.tape() = cur.tape();
		    if ( next.size() < MAX_RECORDS ) next.push( remainder );

		    cur.tape().push_back( cand[ k ] );
		    stitch( cand[ k ], patch, small );
		    resetCycleIDs( patch );
		
#ifdef MYDEBUG
		    cerr << " Stitched via seam : ";
		    printHalfedge( cand[ k ] );
#endif	// MYDEBUG
		
		    ++k;
		}
		else {
#ifdef MYDEBUG
		    cerr << "\a\a\a";
		    cerr << "%%%%%%%%%%%%%%%%%%%% Cannot Stitched seam : ";
		    printHalfedge( cand[ k ] );
		
#endif	// MYDEBUG
		    break;
		}
	    }

#ifdef MYDEBUG
	    cerr << " Pushed seam : ";
	    printHalfedge( cand[ k ] );
#endif	// MYDEBUG
	    
	    hhC = cand[ k ];
	}
    }
	// } while ( hhC != cur.target() );
}


void divide( Polyhedron & poly,
	     Bbox2 & sheet,
	     vector< vector< Facet_handle > > & patch,
	     vector< Bbox2 > & bound )
{
    vector< bool > mark( poly.size_of_vertices() );
    vector< Halfedge_handle > seam;
    Outline start;
#ifdef USE_QUEUE
    queue< Outline > next;
#endif	// USE_QUEUE
#ifdef USE_STACK
    stack< Outline > next;
#endif	// USE_STACK
    
    vector< vector< Facet_handle > > patch2, patchF;
    vector< bool > mark2;
    vector< int > nCuts( poly.size_of_vertices() ), nCutsF( poly.size_of_vertices() );
    vector< bool > connect( poly.size_of_halfedges() ), connectF( poly.size_of_halfedges() );
    vector< int > piece( poly.size_of_facets() ), pieceF( poly.size_of_facets() );
    vector< Halfedge_handle > tapeF;

    for ( unsigned int i = 0; i < mark.size(); ++i ) mark[ i ] = false;

    // while ( true ) {
    Halfedge_handle hhP = nextStart( poly, mark );
    
    cerr << " Picked edge : "; printHalfedge( hhP );
    
    vector< Halfedge_handle > route;
    nextRoute( hhP, route );
    for ( unsigned int k = 0; k < route.size(); ++k ) {
	cerr << " Route [ " << k << " ] = ";
	printHalfedge( route[ k ] );
    }
    
    start.start() = route[ 0 ];
    start.seam().clear();
    start.tape().clear();
    mark[ start.start()->opposite()->vertex()->id() ] = true;
    
    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi )
	nCuts[ vi->id() ] = vi->nCuts();
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi )
	if ( hi->opposite()->vertex()->id() < hi->vertex()->id() )
	    connect[ hi->id() ] = hi->connect();
    for ( Facet_iterator fi = poly.facets_begin(); fi != poly.facets_end(); ++fi )
	piece[ fi->id() ] = fi->piece();

    unsigned int minCount = 0;
    unsigned int nComponents = patch.size();

    next.push( start );

    int curCount = 0; 
    do { 
#ifdef USE_QUEUE
	start = next.front();
#endif	// USE_QUEUE
#ifdef USE_STACK
	start = next.top();
#endif	// USE_STACK
	next.pop();

	// preprocessing
	patch2 = patch;
	mark2 = mark;
	for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi )
	    hi->visit() = false;

	// main
#ifdef MYDEBUG
	cerr << "#################### start edge [ " << setw( 3 ) << curCount << " ] : ";
	printHalfedge( start.source() );
	cerr << "#################### final edge [ " << setw( 3 ) << curCount << " ] : ";
	printHalfedge( start.target() );
	cerr << "#################### seam edges [ " << setw( 3 ) << curCount << " ] : ";
	for ( unsigned int k = 0; k < start.seam().size(); ++k ) {
	    cerr << " Seam[ " << setw( 3 ) << k << " ] : ";
	    printHalfedge( start.seam()[ k ] );
	}
#endif	// MYDEBUG

	resetCycleIDs( patch2 );
	traverse( start, next, mark2, patch2, poly );

	if ( patch2.size() < nComponents ) {
	    patchF = patch2;
	    nComponents = patchF.size();
	    cerr << "************* Found small number of unfoleded patterns N = " << nComponents
		 << " at count " << curCount << endl;
	    minCount = curCount;
	    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi )
		nCutsF[ vi->id() ] = vi->nCuts();
	    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi )
		if ( hi->opposite()->vertex()->id() < hi->vertex()->id() )
		    connectF[ hi->id() ] = hi->connect();
	    for ( Facet_iterator fi = poly.facets_begin(); fi != poly.facets_end(); ++fi )
		pieceF[ fi->id() ] = fi->piece();
	    tapeF = start.tape();
	}

	if ( curCount % 10 == 0 ) 
	    cerr << " Current iteration count = " << setw( 4 ) << curCount
		 << " nComponents = " << setw( 3 ) << nComponents
		 << " next size = " << setw( 4 ) << next.size() << endl;
	curCount++;

	// postprocessing
	for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi )
	    vi->nCuts() = nCuts[ vi->id() ];
	for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi )
	    hi->connect() = connect[ hi->id() ];
	for ( Facet_iterator fi = poly.facets_begin(); fi != poly.facets_end(); ++fi )
	    fi->piece() = piece[ fi->id() ];

	if ( nComponents == 1 ) break;
	if ( curCount >= MAX_SEAMS ) break;
    } while( ! next.empty() );

    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi )
	hi->visit() = false;
    // This must be avoided.
    // patch = patchF;
    resetCycleIDs( patch );
#ifdef MYDEBUG
    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi )
	vi->nCuts() = nCutsF[ vi->id() ];
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi )
	hi->connect() = connectF[ hi->id() ];
    for ( Facet_iterator fi = poly.facets_begin(); fi != poly.facets_end(); ++fi )
	fi->piece() = pieceF[ fi->id() ];
#endif	// MYDEBUG
    // preprocessing: stitching edges
    vector< Triangle2 > small;
    for ( unsigned int i = 0; i < tapeF.size(); ++i ) {
#ifdef MYDEBUG
	cerr << " Connecting the edge : ";
	printHalfedge( tapeF[ i ] );
#endif	// MYDEBUG
	if ( canEmbed( tapeF[ i ], patch, small ) ) {
	    stitch( tapeF[ i ], patch, small );
	    resetCycleIDs( patch );
	}
	else assert( false );
    }

    cerr << " Number of next candidate outlines = " << next.size() << endl;
#ifdef MYDEBUG
    count = 0;
    while ( ! next.empty() ) {
	cerr << " Start[ " << setw( 3 ) << count++ << " ] = ";
	printHalfedge( next.front().source() );
	for ( unsigned int i = 0; i < next.front().seam().size(); ++i ) {
	    cerr << " Seam[ " << setw( 3 ) << i << " ] : ";
	    printHalfedge( next.front().seam()[ i ] );
	}
	for ( unsigned int i = 0; i < next.front().tape().size(); ++i ) {
	    cerr << " Tape[ " << setw( 3 ) << i << " ] : ";
	    printHalfedge( next.front().tape()[ i ] );
	}
	next.pop();
    }
#endif	// MYDEBUG

    arrangeSingle( patch, sheet, bound );

    cerr << "\a\a\a";
}
