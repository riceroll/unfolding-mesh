//------------------------------------------------------------------------------
//
//	stitch.cpp
//
//------------------------------------------------------------------------------
#include "common.h"
#include <vector>

//------------------------------------------------------------------------------
//	To validate the stitchability for each edge
//------------------------------------------------------------------------------

// Check if the edge is locally stitchable or not
bool isLocally( const Halfedge_handle & hh )
{
    Vertex_handle vhO = hh->vertex();
    Vertex_handle vhD = hh->opposite()->vertex();

#ifdef DEBUG
    cerr << " hh->facet()->piece() = " << hh->facet()->piece();
    cerr << " hh->opposite()->facet()->piece() = " << hh->opposite()->facet()->piece() << endl;
    if ( vhO->label() == SADDLE_VERTEX ) cerr << " vhO is saddle";
    else cerr << " vhO is not saddle";
    cerr << " vhO->nCuts = " << vhO->nCuts() << endl;
    if ( vhD->label() == SADDLE_VERTEX ) cerr << " vhD is saddle";
    else cerr << " vhD is not saddle";
    cerr << " vhD->nCuts = " << vhD->nCuts() << endl;
#endif	// DEBUG

    if ( 
	( ( hh->facet()->piece() != hh->opposite()->facet()->piece() ) &&
	  ( hh->facet()->piece() != NO_INDEX ) )
	&&
	( ( vhO->label() != SADDLE_VERTEX ) ||
	  ( ( vhO->label() == SADDLE_VERTEX ) && ( vhO->nCuts() >= 3 ) ) )
	&&
	( ( vhD->label() != SADDLE_VERTEX ) ||
	  ( ( vhD->label() == SADDLE_VERTEX ) && ( vhD->nCuts() >= 3 ) ) )
	 ) {
	return true;
    }
    else {
	return false;
    }
}


Segment2 sideTriangle( Halfedge_handle & key, Facet_handle & fh )
{
    Point2 coordS, coordT;

    if ( fh->halfedge() == key ) {
	coordS = key->facet()->triangle()[ 0 ];
	coordT = key->facet()->triangle()[ 1 ];
    }
    else if ( fh->halfedge()->next() == key ) {
	// right shoulder
	coordS = key->facet()->triangle()[ 1 ];
	coordT = key->facet()->triangle()[ 2 ];
    }
    else if ( fh->halfedge()->prev() == key ) {
	// left shoulder
	coordS = key->facet()->triangle()[ 2 ];
	coordT = key->facet()->triangle()[ 0 ];
    }
    else assert( false );

    return Segment2( coordS, coordT );
}


//------------------------------------------------------------------------------
// Seeks the projected 2D coordinates of the given half edge
//------------------------------------------------------------------------------
Segment2 seekStrip( Halfedge_handle & key, vector< Facet_handle > & subpatch )
{
#ifdef MYDEBUG
    cerr << " key = ";
    printHalfedge( key );
    cerr << " key->facet()->id() " << key->facet()->id() << endl;
    cerr << " subpatch.size() = " << subpatch.size() << endl;
    cerr << " key->facet()->piece() = " << key->facet()->piece() << endl;
    cerr << " subpatch[ 0 ]->piece() = " << subpatch[ 0 ]->piece() << endl;
#endif	// MYDEBUG

    if ( key->facet()->piece() != subpatch[ 0 ]->piece() ) {
	cerr << " key->facet()->piece() = " << key->facet()->piece() << endl;
	cerr << " subpatch[ 0 ]->piece() = " << subpatch[ 0 ]->piece() << endl;
	assert( key->facet()->piece() == subpatch[ 0 ]->piece() );
    }

    int facetID = key->facet()->id();

    Facet_handle fh;
    for ( unsigned int k = 0; k < subpatch.size(); ++k ) {
	if ( subpatch[ k ]->id() == facetID ) {
	    fh = subpatch[ k ];
	    break;
	}
    }

    return sideTriangle( key, fh );
}


//------------------------------------------------------------------------------
// Seeks the projected 2D coordinates of the given half edge
//------------------------------------------------------------------------------
Segment2 splitStrip( Halfedge_handle & key, vector< Facet_handle > & subpatch, vector< Facet_handle > & newly )
{
    assert( key->facet()->piece() == subpatch[ 0 ]->piece() );

    int facetID = key->facet()->id();
    vector< Facet_handle > overall = subpatch;
    subpatch.clear();
    newly.clear();

    Facet_handle fh;
    for ( unsigned int k = 0; k < overall.size(); ++k ) {
	if ( overall[ k ]->id() == facetID ) {
	    fh = overall[ k ];
	    subpatch.push_back( key->facet() );
	    cerr << " insert edge : " << key->vertex()->id() << " == " << key->opposite()->vertex()->id() << endl;

	    for ( unsigned int m = k+1; m < overall.size(); ++m ) {
		subpatch.push_back( overall[ m ] );
	    }
	    for ( int m = k-1; m >= 0; --m ) {
		newly.push_back( overall[ m ] );
	    }
	    break;
	}
    }

    return sideTriangle( key, fh );
}


//------------------------------------------------------------------------------
//	Global self-intersection validation
//------------------------------------------------------------------------------
void transformStrip( Halfedge_handle & bridge,
		     vector< vector< Facet_handle > > & patch,
		     vector< Triangle2 > & transformed )
{
    Facet_handle	fhI	= bridge->facet();
    Facet_handle	fhO	= bridge->opposite()->facet();
    // Vertex_handle	vhD	= bridge->vertex();
    // Vertex_handle	vhS	= bridge->opposite()->vertex();
    
    // A is the large unfolded pattern
    Halfedge_handle	hhA	= bridge->opposite();
    Halfedge_handle	hhB	= bridge;
    int			cycleIDA = fhO->piece();
    int			cycleIDB = fhI->piece();
    Segment2		edgeA	= seekStrip( hhA, patch[ cycleIDA ] );
    Segment2		edgeB	= seekStrip( hhB, patch[ cycleIDB ] );

#ifdef MYDEBUG
    cerr << " edgeA = " << edgeA << endl;
    cerr << " edgeB = " << edgeB << endl;
#endif	// MYDEBUG

    Vector2		vecA	= edgeA.source() - edgeA.target();
    Vector2		vecB	= edgeB.target() - edgeB.source();
    double		lengthA = sqrt( edgeA.squared_length() );
    double		lengthB = sqrt( edgeB.squared_length() );

    Vector2		position = edgeB.source() - CGAL::ORIGIN;
    Vector2		displacement = edgeA.target() - edgeB.source();
    double		cosRot	= MIN2( 1.0f, MAX2( -1.0f, (vecA*vecB)/(lengthA*lengthB) ) );
    double		sinRot	= sqrt( 1.0 - SQUARE( cosRot ) );

#ifdef MYDEBUG
    cerr << " position = " << position << endl;
    cerr << " displacement = " << displacement << endl;
    cerr << " cosRot = " << cosRot << endl;
    cerr << " sinRot = " << sinRot << endl;
#endif	// MYDEBUG

	// Returns the vector perpendicular to vecA in clockwise or counterclockwise orientation
    if ( vecA.perpendicular( CGAL::COUNTERCLOCKWISE )*vecB > 0 ) sinRot *= (-1.0);

    transformed.clear();
    for ( unsigned int k = 0; k < patch[ cycleIDB ].size(); ++k )
	transformed.push_back( patch[ cycleIDB ][ k ]->triangle() );

    Transformation2 translate0( CGAL::TRANSLATION, -position );
    Transformation2 rotate( CGAL::ROTATION, sinRot, cosRot );
    Transformation2 translate1( CGAL::TRANSLATION, position );
    Transformation2 translate2( CGAL::TRANSLATION, displacement );
    Transformation2 composite = translate2*translate1*rotate*translate0;

    Segment2 edgeN = edgeB.transform( composite );
#ifdef MYDEBUG
    cerr << " edgeN = " << edgeN << endl;
#endif	// MYDEBUG

    for ( unsigned int j = 0; j < transformed.size(); ++j )
	transformed[ j ] = transformed[ j ].transform( composite );
}


bool isGlobally( vector< Facet_handle > & fixed,
		 vector< Triangle2 > & transformed )
{
    for ( unsigned int i = 0; i < fixed.size(); ++i ) {
	if ( isIntersected( fixed[ i ]->triangle(), transformed ) ) {
#ifdef MYDEBUG
	    cerr << " Edge : " << bridge->opposite()->vertex()->id()
		 << " == " << bridge->vertex()->id()
		 << ", Triangle : "
		 << patch[ cycleID ][ i ]->halfedge()->prev()->vertex()->id() << " == "
		 << patch[ cycleID ][ i ]->halfedge()->vertex()->id() << " == "
		 << patch[ cycleID ][ i ]->halfedge()->next()->vertex()->id()
		 << " has self-intersection." << endl;
#endif	// MYDEBUG
	    return false;
	}
    }

    return true;
}


bool doEmbed( Halfedge_handle & bridge,
	      vector< vector< Facet_handle > > & patch,
	      vector< Triangle2 > & transformed )
{
    int thisID = bridge->opposite()->facet()->piece();
    int thatID = bridge->facet()->piece();

    // Unfoled patch with a smaller number of faces should be transformed
    if ( patch[ thisID ].size() < patch[ thatID ].size() ) {
	bridge = bridge->opposite();
    }

    transformStrip( bridge, patch, transformed );

    // ID of the cycle to be fixed
    int	cycleID = bridge->opposite()->facet()->piece();

#ifdef DEBUG
    cerr << " patch[ cycleID ].size = " << patch[ cycleID ].size()
	 << " transformed.size = " << transformed.size() << endl;
#endif	// DEBUG
    // assert( patch[ cycleID ].size() >= transformed.size() );  
    // if ( ( int )transformed.size() > 100 ) cerr << "L";

    if ( ! isGlobally( patch[ cycleID ], transformed ) ) {
	return false;
    }

    return true;
}    


bool canEmbed( Halfedge_handle & bridge,
	       vector< vector< Facet_handle > > & patch,
	       vector< Triangle2 > & transformed )
{
    // If the edge has been manually specified as a cut edge
    if ( bridge->fixed() ) {
	return false;
    }

    if ( ! isLocally( bridge ) ) {
	// cerr << " Numbers of cuts around saddles are inappropriate" << endl;
	return false;
    }

    int thisID = bridge->opposite()->facet()->piece();
    int thatID = bridge->facet()->piece();

    // Unfoled patch with a smaller number of faces should be transformed
    if ( patch[ thisID ].size() < patch[ thatID ].size() ) {
	bridge = bridge->opposite();
    }

    transformStrip( bridge, patch, transformed );

    // ID of the cycle to be fixed
    int	cycleID = bridge->opposite()->facet()->piece();

#ifdef DEBUG
    cerr << " patch[ cycleID ].size = " << patch[ cycleID ].size()
	 << " transformed.size = " << transformed.size() << endl;
#endif	// DEBUG
    // assert( patch[ cycleID ].size() >= transformed.size() );  
    // if ( ( int )transformed.size() > 100 ) cerr << "L";

    if ( ! isGlobally( patch[ cycleID ], transformed ) ) {
	return false;
    }

    return true;
}    


//------------------------------------------------------------------------------
//	Stitching operation
//------------------------------------------------------------------------------
void stitch( Halfedge_handle & bridge,
	     vector< vector< Facet_handle > > & patch,
	     vector< Triangle2 > & transformed )
{
    // comment out for the use of MST-based unfolding
    // assert( isLocally( bridge ) );

    Facet_handle	fhI	= bridge->facet();
    Facet_handle	fhO	= bridge->opposite()->facet();
    Vertex_handle	vhD	= bridge->vertex();
    Vertex_handle	vhS	= bridge->opposite()->vertex();
    int			cycleIDA = fhO->piece();
    int			cycleIDB = fhI->piece();

#ifdef MYDEBUG
    cerr << " Bridge : " << bridge->opposite()->vertex()->id()
	 << " -- " << bridge->vertex()->id() << endl;
#endif	// MYDEBUG

    if ( cycleIDA == cycleIDB ) {
	cerr << " The same cycle IDs : " << cycleIDA << " vs." << cycleIDB << endl;
	assert( cycleIDA != cycleIDB );
    }

    // parameter adjustments
    bridge->connect() = bridge->opposite()->connect() = true;
    vhS->nCuts()--;
    vhD->nCuts()--;
    
    assert( patch[ cycleIDB ].size() == transformed.size() );
    for ( unsigned int k = 0; k < patch[ cycleIDB ].size(); ++k )
	patch[ cycleIDB ][ k ]->triangle() = transformed[ k ];

    for ( unsigned int k = 0; k < patch[ cycleIDB ].size(); ++k ) {
	patch[ cycleIDA ].push_back( patch[ cycleIDB ][ k ] );
    }
    patch[ cycleIDB ].clear();
    patch.erase( patch.begin() + cycleIDB );
}    


//------------------------------------------------------------------------------
//	
bool mergeBySize( Polyhedron & poly, 
		  Bbox2 & sheet,
		  vector< vector< Facet_handle > > & patch,
		  vector< Bbox2 > & bound )
{
    if ( patch.size() <= 1 ) {
	cerr << "We cannot merge unfolded patterns any more." << endl;
	return false;
    }
    
    multimap< int, Halfedge_handle > mmap;
    const int D = 3;
    Halfedge_handle hh[ 3 ];
    vector< int > nStitch;
    
    for ( unsigned int i = 0; i < patch.size(); ++i ) {
	int count = 0;
	for ( unsigned int j = 0; j < patch[ i ].size(); ++j ) {
	    hh[ 0 ] = patch[ i ][ j ]->halfedge()->prev();
	    hh[ 1 ] = patch[ i ][ j ]->halfedge();
	    hh[ 2 ] = patch[ i ][ j ]->halfedge()->next();
	    for ( int k = 0; k < D; ++k ) {
		if ( isLocally( hh[ k ] ) ) count++; // Stitchable validation for each edge
	    }
	}
	nStitch.push_back( count );
    }
    
#ifdef MYDEBUG
    for ( unsigned int i = 0; i < nStitch.size(); ++i )
	cerr << " nStitch[ " << i << "] = " << nStitch[ i ] << endl;
#endif	// MYDEBUG
    
    for ( unsigned int i = 0; i < patch.size(); ++i ) {
	int idI = i;
	for ( unsigned int j = 0; j < patch[ i ].size(); ++j ) {
	    hh[ 0 ] = patch[ i ][ j ]->halfedge()->prev();
	    hh[ 1 ] = patch[ i ][ j ]->halfedge();
	    hh[ 2 ] = patch[ i ][ j ]->halfedge()->next();
	    for ( int k = 0; k < D; ++k ) {
		if ( isLocally( hh[ k ] ) ) {
		    int idO = hh[ k ]->opposite()->facet()->piece();
		    if ( nStitch[ idI ] <= nStitch[ idO ] ) {
//			mmap.insert( make_pair< int, Halfedge_handle >( nStitch[ idI ], hh[ k ] ) );	error
				mmap.insert( make_pair( nStitch[ idI ], hh[ k ] ) );
		    }
		}
	    }
	}
    }
    
    multimap< int, Halfedge_handle >::iterator it;

#ifdef MYDEBUG
    it = mmap.begin();
    while ( it != mmap.end() ) {
	cerr << " minTris = " << it->first << " Edge : " 
	     << it->second->opposite()->vertex()->id() << " == " << it->second->vertex()->id() << endl;
	it++;
    }
#endif	// MYDEBUG

    // int count = 0;
    vector< Triangle2 > small;
    it = mmap.begin();
    // Pick the first bridge edge
    while ( it != mmap.end() ) {
	if ( canEmbed( it->second, patch, small ) ) {
	    stitch( it->second, patch, small );
#ifdef MYDEBUG
	    cerr << "[Stitched] minTris = " << it->first << " Edge : " 
		 << it->second->opposite()->vertex()->id() << " == " << it->second->vertex()->id() << endl;
#endif	// MYDEBUG
	    cerr << "*";
	    resetCycleIDs( patch );
#ifdef MYDEBUG
	    cerr << "Number of unfolded patterns = " << patch.size() << endl;
#endif	// MYDEBUG
	    return true;
	}
	else cerr << ".";
	it++;
    }
    cerr << endl;
#ifdef MYDEBUG
    cerr << " Cannot find edges that can be stitched" << endl;
    it = mmap.begin();
    while ( it != mmap.end() ) {
	cerr << " minTris = " << it->first << " Edge : " 
	     << it->second->opposite()->vertex()->id() << " == " << it->second->vertex()->id() << endl;
	it++;
    }
#endif	// MYDEBUG
    return false;
}


//------------------------------------------------------------------------------
//	
bool mergeByWeight( Polyhedron & poly, 
		    Bbox2 & sheet,
		    vector< vector< Facet_handle > > & patch,
		    vector< Bbox2 > & bound )
{
    if ( patch.size() <= 1 ) {
	cerr << "We cannot merge unfolded patterns any more." << endl;
	return false;
    }
    
    multimap< double, Halfedge_handle > mmap;
    const int D = 3;
    Halfedge_handle hh[ 3 ];
    vector< int > nStitch;
    
    for ( unsigned int i = 0; i < patch.size(); ++i ) {
	int count = 0;
	for ( unsigned int j = 0; j < patch[ i ].size(); ++j ) {
	    hh[ 0 ] = patch[ i ][ j ]->halfedge()->prev();
	    hh[ 1 ] = patch[ i ][ j ]->halfedge();
	    hh[ 2 ] = patch[ i ][ j ]->halfedge()->next();
	    for ( int k = 0; k < D; ++k ) {
		if ( isLocally( hh[ k ] ) ) count++; // Stitchable validation for each edge
	    }
	}
	nStitch.push_back( count );
    }
    
#ifdef MYDEBUG
    for ( unsigned int i = 0; i < nStitch.size(); ++i )
	cerr << " nStitch[ " << i << "] = " << nStitch[ i ] << endl;
#endif	// MYDEBUG
    
    for ( unsigned int i = 0; i < patch.size(); ++i ) {
	int idI = i;
	for ( unsigned int j = 0; j < patch[ i ].size(); ++j ) {
	    hh[ 0 ] = patch[ i ][ j ]->halfedge()->prev();
	    hh[ 1 ] = patch[ i ][ j ]->halfedge();
	    hh[ 2 ] = patch[ i ][ j ]->halfedge()->next();
	    for ( int k = 0; k < D; ++k ) {
		if ( isLocally( hh[ k ] ) ) {
		    int idO = hh[ k ]->opposite()->facet()->piece();
		    if ( nStitch[ idI ] <= nStitch[ idO ] ) {
//			mmap.insert( make_pair< int, Halfedge_handle >( hh[ k ]->weight(), hh[ k ] ) );	error
				mmap.insert( make_pair( hh[ k ]->weight(), hh[ k ] ) );
		    }
		}
	    }
	}
    }

    multimap< double, Halfedge_handle >::iterator it;

#ifdef MYDEBUG
    it = mmap.begin();
    while ( it != mmap.end() ) {
	cerr << " minTris = " << it->first << " Edge : " 
	     << it->second->opposite()->vertex()->id() << " == " << it->second->vertex()->id() << endl;
	it++;
    }
#endif	// MYDEBUG

    vector< Triangle2 > small;
    it = mmap.begin();
    while ( it != mmap.end() ) {
	cerr << "*";
	if ( canEmbed( it->second, patch, small ) ) {
	    stitch( it->second, patch, small );
#ifdef MYDEBUG
	    cerr << "[Stitched] minTris = " << it->first << " Edge : " 
		 << it->second->opposite()->vertex()->id() << " == " << it->second->vertex()->id() << endl;
#endif	// MYDEBUG
	    resetCycleIDs( patch );
#ifdef MYDEBUG
	    cerr << "Number of unfolded patterns = " << patch.size() << endl;
#endif	// MYDEBUG
	    return true;
	}
	it++;
    }
    cerr << endl;
#ifdef MYDEBUG
    cerr << " Cannot find edges that can be stitched" << endl;
    it = mmap.begin();
    while ( it != mmap.end() ) {
	cerr << " minTris = " << it->first << " Edge : " 
	     << it->second->opposite()->vertex()->id() << " == " << it->second->vertex()->id() << endl;
	it++;
    }
#endif	// MYDEBUG
    return false;
}


//------------------------------------------------------------------------------
//	
//------------------------------------------------------------------------------
bool mergeByOne( Polyhedron & poly, 
		 Bbox2 & sheet,
		 vector< vector< Facet_handle > > & patch,
		 vector< Bbox2 > & bound )
{
    return mergeBySize( poly, sheet, patch, bound );
    // return mergeByWeight( poly, sheet, patch, bound );
}

//------------------------------------------------------------------------------
//	
//------------------------------------------------------------------------------
void greedy( Polyhedron & poly, 
	     Bbox2 & sheet,
	     vector< vector< Facet_handle > > & patch,
	     vector< Bbox2 > & bound )
{
    while ( true ) {
	if ( ! mergeByOne( poly, sheet, patch, bound ) ) break;
    }

    arrangeSingle( patch, sheet, bound );
}


//------------------------------------------------------------------------------
//	Re-Split and Re-Stitch
//------------------------------------------------------------------------------
bool mergeWithCut( Polyhedron & poly, 
		   Bbox2 & sheet,
		   vector< vector< Facet_handle > > & patch,
		   vector< Bbox2 > & bound )
{
    if ( patch.size() <= 1 ) {
	cerr << "We cannot merge unfolded patterns any more." << endl;
	return false;
    }
    
    multimap< int, Halfedge_handle > mmap;
    const int D = 3;
    Halfedge_handle hh[ 3 ];
    vector< int > nStitch;
    
    for ( unsigned int i = 0; i < patch.size(); ++i ) {
	int count = 0;
	for ( unsigned int j = 0; j < patch[ i ].size(); ++j ) {
	    hh[ 0 ] = patch[ i ][ j ]->halfedge()->prev();
	    hh[ 1 ] = patch[ i ][ j ]->halfedge();
	    hh[ 2 ] = patch[ i ][ j ]->halfedge()->next();
	    for ( int k = 0; k < D; ++k ) {
		// if ( isLocally( hh[ k ] ) ) count++; // Stitchable validation for each edge
		if ( ! hh[ k ]->connect() ) count++; // Stitchable validation for each edge
	    }
	}
	if ( count == 0 ) {
	    cerr << " cycle ID = " << i << endl;
	    assert( count != 0 );
	}
	nStitch.push_back( count );
    }
    
#ifdef MYDEBUG
    for ( unsigned int i = 0; i < nStitch.size(); ++i )
	cerr << " nStitch[ " << i << "] = " << nStitch[ i ] << endl;
#endif	// MYDEBUG
    
    for ( unsigned int i = 0; i < patch.size(); ++i ) {
	int idI = i;
	for ( unsigned int j = 0; j < patch[ i ].size(); ++j ) {
	    hh[ 0 ] = patch[ i ][ j ]->halfedge()->prev();
	    hh[ 1 ] = patch[ i ][ j ]->halfedge();
	    hh[ 2 ] = patch[ i ][ j ]->halfedge()->next();
	    for ( int k = 0; k < D; ++k ) {
		// if ( isLocally( hh[ k ] ) ) {
		int idO = hh[ k ]->opposite()->facet()->piece();
		if ( ( ! hh[ k ]->connect() ) && ( idI != idO ) ) {
		    if ( nStitch[ idI ] <= nStitch[ idO ] ) {
			mmap.insert( make_pair( nStitch[ idI ], hh[ k ] ) );
		    }
		}
	    }
	}
    }

    if ( mmap.size() == 0 ) {
	cerr << "\a\a\a" << ends;
	cerr << "There are no stitchable edges" << endl;
	return false;
    }
    
    multimap< int, Halfedge_handle >::iterator it;
#ifdef MYDEBUG
    it = mmap.begin();
    while ( it != mmap.end() ) {
	cerr << " minTris = " << it->first << " Edge : " 
	     << it->second->opposite()->vertex()->id() << " == " << it->second->vertex()->id() << endl;
	it++;
    }
#endif	// MYDEBUG

    vector< Triangle2 > small;
    it = mmap.begin();
    bool isConnected = false;
    while ( it != mmap.end() ) {
	// Here we do not care the patch has global intersection or not
	cerr << "Stitch the edge : " << it->second->opposite()->vertex()->id()
	     << " == " << it->second->vertex()->id() << endl;
	Halfedge_handle bridge = it->second;
	doEmbed( bridge, patch, small );
	vector< Facet_handle > setA = patch[ bridge->facet()->piece() ];
	vector< Facet_handle > setB = patch[ bridge->opposite()->facet()->piece() ];
	stitch( bridge, patch, small );

	resetCycleIDs( patch );
	bool success = resplitOverlaps( poly, patch, setA, setB );
	if ( success ) {
	    cerr << "Successful!!" << endl;
	    isConnected = true;
	    break;
	}
	else {
	    cerr << "Failed!!" << endl;
	    cutStrips( patch, bridge );
	    bridge->visit() = bridge->opposite()->visit() = true;
	}
	it++;
    }

    return isConnected;
}


void reconnectByOne( Polyhedron & poly, 
		     Bbox2 & sheet,
		     vector< vector< Facet_handle > > & patch,
		     vector< Bbox2 > & bound )
{
    cerr << "Reconnecting the unfolded patches" << endl;
    bool success = mergeWithCut	( poly, sheet, patch, bound );
    arrangeSingle		( patch, sheet, bound );

    if ( success ) cerr << "reconnection successful" << endl;
    else cerr << "reconnection failed!!" << "\a\a\a" << endl;
}    


void reconnect( Polyhedron & poly, 
		Bbox2 & sheet,
		vector< vector< Facet_handle > > & patch,
		vector< Bbox2 > & bound,
		Attribute & attr,
		void (*preview)( void ) )
{
    vector< vector < Facet_handle > >	sPatch,	  fPatch;
    vector< vector < Triangle2 > >	sTriplet, fTriplet;
    vector< bool >			sConnect, fConnect;
    vector< int >			sNCuts,   fNCuts;
    Bbox2				sSheet,   fSheet;


    // store the current status
    countNCuts( poly );
//------------------------------------------------------------------------------
//	Store the current status
//------------------------------------------------------------------------------
    sPatch		= patch;
    sSheet		= sheet;

    sTriplet.resize( patch.size() );
    for ( unsigned int i = 0; i < patch.size(); ++i ) {
	sTriplet[ i ].resize( patch[ i ].size() );
	for ( unsigned int j = 0; j < patch[ i ].size(); ++j ) {
	    sTriplet[ i ][ j ] = patch[ i ][ j ]->triangle();
	}
    }

    sConnect.resize( poly.size_of_halfedges() );
    fConnect.resize( poly.size_of_halfedges() );
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi )
	sConnect[ hi->id() ]	= hi->connect();
    sNCuts.resize( poly.size_of_vertices() );
    fNCuts.resize( poly.size_of_vertices() );
    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi )
	sNCuts[ vi->id() ]	= vi->nCuts();

    int minPatches = ( int )patch.size();
    int maxFaces = 0;
    for ( unsigned int i = 0; i < patch.size(); ++i )
	if ( ( int )patch[ i ].size() > maxFaces )
	    maxFaces = ( int )patch[ i ].size();

    // call this first
    labelBoundary( poly, attr );

    int bestN = 0;
    const int n = MAX_RETRIALS;
    // const int n = 4;
    // greedy( poly, sheet, patch, bound );
    for ( int k = 1; k <= n; ++k ) {
	for ( int m = 0; m < k; ++m ) {
	    if ( ! mergeWithCut( poly, sheet, patch, bound ) ) break;
	    else
		;
		// arrangeCycles( patch, sheet, bound );
	    greedy( poly, sheet, patch, bound );
	}
	// greedy( poly, sheet, patch, bound );

	int curPatches = ( int )patch.size();
	int curFaces = 0;
	for ( unsigned int i = 0; i < patch.size(); ++i )
	    if ( ( int )patch[ i ].size() > curFaces )
		curFaces = ( int )patch[ i ].size();

	cerr << " k = " << k << " #patches = " << curPatches << " #faces = " << curFaces << endl;
	if ( curPatches < minPatches ) {
	    minPatches = curPatches;
	    maxFaces = curFaces;
	    bestN = k;

	    bound.clear();		// This is very important
	    arrangeSingle( patch, sheet, bound );
	    labelBoundary( poly, attr );
	    (*preview)();
	    
	    fPatch		= patch;
	    fSheet		= sheet;
	    fTriplet.resize( patch.size() );
	    for ( unsigned int i = 0; i < patch.size(); ++i ) {
		fTriplet[ i ].resize( patch[ i ].size() );
		for ( unsigned int j = 0; j < patch[ i ].size(); ++j )
		    fTriplet[ i ][ j ] = patch[ i ][ j ]->triangle();
	    }
	    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi )
		fConnect[ hi->id() ]	= hi->connect();
	    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi )
		fNCuts[ vi->id() ]	= vi->nCuts();
	    cerr << "Go out of the loop" << endl;
	    break;
	}
	else if ( ( curPatches == minPatches ) && ( curFaces > maxFaces ) ) {
	    minPatches = curPatches;
	    maxFaces = curFaces;
	    bestN = k;

	    bound.clear();		// This is very important
	    arrangeSingle( patch, sheet, bound );
	    (*preview)();

	    fPatch		= patch;
	    fSheet		= sheet;
	    fTriplet.resize( patch.size() );
	    for ( unsigned int i = 0; i < patch.size(); ++i ) {
		fTriplet[ i ].resize( patch[ i ].size() );
		for ( unsigned int j = 0; j < patch[ i ].size(); ++j )
		    fTriplet[ i ][ j ] = patch[ i ][ j ]->triangle();
	    }
	    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi )
		fConnect[ hi->id() ]	= hi->connect();
	    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi )
		fNCuts[ vi->id() ]	= vi->nCuts();
	}
	
	// restore the original mesh status
	patch = sPatch;
	for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	    hi->connect() = sConnect[ hi->id() ];
	}
	for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	    vi->nCuts() = sNCuts[ vi->id() ];
	}
	for ( unsigned int i = 0; i < patch.size(); ++i ) {
	    for ( unsigned int j = 0; j < patch[ i ].size(); ++j ) {
		patch[ i ][ j ]->triangle() = sTriplet[ i ][ j ];
	    }
	}
	resetCycleIDs( patch );
	countNCuts( poly );
	sheet = sSheet;
	// call this first
	labelBoundary( poly, attr );
    }	

    if ( bestN != 0 ) {
	// restore the final mesh status
	patch = fPatch;
	for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	    hi->connect() = fConnect[ hi->id() ];
	}
	for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	    vi->nCuts() = fNCuts[ vi->id() ];
	}
	for ( unsigned int i = 0; i < patch.size(); ++i ) {
	    for ( unsigned int j = 0; j < patch[ i ].size(); ++j ) {
		patch[ i ][ j ]->triangle() = fTriplet[ i ][ j ];
	    }
	}
    }
    else {
	// restore the original mesh status
	patch = sPatch;
	for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	    hi->connect() = sConnect[ hi->id() ];
	}
	for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	    vi->nCuts() = sNCuts[ vi->id() ];
	}
	for ( unsigned int i = 0; i < patch.size(); ++i ) {
	    for ( unsigned int j = 0; j < patch[ i ].size(); ++j ) {
		patch[ i ][ j ]->triangle() = sTriplet[ i ][ j ];
	    }
	}
    }

    resetCycleIDs( patch );
    countNCuts( poly );
    
    bound.clear();		// This is very important
    arrangeSingle( patch, sheet, bound );
    // call this at last
    labelBoundary( poly, attr );
    cerr << "\a\a\a";
    if ( bestN == 0 ) cerr << "Rearrangement failed" << endl;
    else cerr << "Rearrangement successful" << endl;
    cerr << "[Final: ] bestN = " << bestN << " #patches = " << minPatches << " #faces = " << maxFaces << endl;
    cerr << "%%%%% paper size = "
	 << ( sheet.xmax() - sheet.xmin() ) * ( sheet.ymax() - sheet.ymin() )
	 << "%%%%%" << endl;
}
