//------------------------------------------------------------------------------
//
//	mst.cpp (minimum spanning tree)
//
//------------------------------------------------------------------------------
#include "common.h"
#include "ui.h"

#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <CGAL/Polyhedron_3.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

#include <iostream>
#include <iomanip>
#include <cassert>
using namespace std;


#define WEIGHT_BASED_SPLIT

//------------------------------------------------------------------------------
//	MST computation
//------------------------------------------------------------------------------
void calcMST( Graph & dual, Polyhedron & poly ) 
{	     
    // Initialization of edge connections
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	hi->connect() = false;
    }

    EdgeIDMap weight = get( edge_weight, dual );
    std::vector < EdgeDescriptor > spanning_tree;
    
    kruskal_minimum_spanning_tree( dual, std::back_inserter( spanning_tree ) );
    std::cout << "Number of edges for spanning tree: " << spanning_tree.size() << std::endl;
	
#ifdef MYDEBUG
    std::cout << "Print the edges in the MST: " << std::endl;
#endif	// MYDEBUG
    for ( std::vector < EdgeDescriptor >::iterator spanning_ei = spanning_tree.begin(); 
	  spanning_ei!=spanning_tree.end(); ++spanning_ei ) {		
#ifdef MYDEBUG
	// std::cout << *spanning_ei << std::endl;
	std::cout << source( *spanning_ei, dual ) << " <--> " << target( *spanning_ei, dual ) 
		  << " with weight of " << weight[ *spanning_ei ] << std::endl;
#endif	// MYDEBUG
	
	for ( Halfedge_iterator hi  = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	    unsigned int origID = hi->facet()->id();
	    unsigned int destID = hi->opposite()->facet()->id();
	    
	    if( origID < destID ) {
		if ( ( source( *spanning_ei, dual ) == origID ) &&
		     ( target( *spanning_ei, dual ) == destID ) ) {
#ifdef MYDEBUG
		    cerr << "hi: " << hi->vertex()->id() << ", hi->oppposite: " 
			 << hi->opposite()->vertex()->id() << endl;
#endif	// MYDEBUG
		    hi->match() = hi->opposite()->match() = KRUSKAL_EDGE;
		    hi->connect() = hi->opposite()->connect() = true;
		}
	    }		
	}			
	//cerr << temp.at(0)->vertex()->id() << "," << temp.at(0)->next()->vertex()->id() << "," << temp.at(0)->prev()->vertex()->id() << endl;
	//cerr << temp.at(1)->vertex()->id() << "," << temp.at(1)->next()->vertex()->id() << "," << temp.at(1)->prev()->vertex()->id() << endl;
	
    }
}

//------------------------------------------------------------------------------
//	Extract the list of triangular strips from the MST
//------------------------------------------------------------------------------
void incTreeEdges( const Halfedge_handle & start, vector< Halfedge_handle > & inc )
{
    inc.clear();
    Halfedge_handle cur = start;
    do {
	cur = cur->next();
	if ( cur->connect() ) {
	    inc.push_back( cur );
	}
    } while ( cur != start );

#ifdef MYDEBUG
    for ( unsigned int i = 0; i < inc.size(); ++i ) {
	cerr << "[ " << setw( 3 ) << i << " ]: cur->id() = " << cur->id()
	     << " start->id() = " << start->id() << endl;
    }
#endif	// MYDEBUG
}


void MSTtoStrips( Polyhedron & poly,
		  vector< vector< Halfedge_handle > > & spine,
		  vector< Halfedge_handle > & stitch )
{
    // GraphTraits::edge_iterator ei, ei_end;
    // VertexIDMap VIDmap = get( vertex_index, tree );
    int stripID = 0;
    // Initialize the input array of triangular strips
    spine.clear();
    stitch.clear();

    // Find one endpoint of the MST as the starting point for tracing

    // flag for the existence of a triangular strip in MST.
    bool isExist = false;
    Halfedge_handle startHH;
    
    // for each dual edge in the MST
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	if ( hi->connect() ) {
	    // Extracting the deadend halfedge of the minimum spanning tree
	    // (There must be at least 2.)
	    vector< Halfedge_handle > srcVHH, tarVHH; 
	    incTreeEdges( hi, srcVHH );
	    incTreeEdges( hi->opposite(), tarVHH );
	    if ( srcVHH.size() == 1 ) {
		// we can start a new trianuglar strip from that dual vertex
		isExist = true;
		startHH = hi;
		break;
	    }
	    else if ( tarVHH.size() == 1 ) {
		isExist = true;
		startHH = hi->opposite();
		break;
	    }
	}
    }

    //------------------------------------------------------------------------------
    //	If we cannot find the endpoints of the MST...
    //------------------------------------------------------------------------------
    if ( ! isExist ) {
	cerr << "We cannot find the endpoints of the MST!!" << endl;
	assert( isExist );
    }
    
    //------------------------------------------------------------------------------
    //	If we can successfully find one endpoint of the MST...
    //------------------------------------------------------------------------------

#ifdef MYDEBUG
    cerr << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    cerr << " starting MST edge : " << startHH->opposite()->vertex()->id() << " == " 
	 << startHH->vertex()->id() << endl;
    cerr << " starting MST dual edge : " << startHH->facet()->id() << " == " 
	 << startHH->opposite()->facet()->id() << endl;
#endif	// MYDEBUG

    //------------------------------------------------------------------------------
    //	Start to track a triangular strip on the remaing part of the MST
    //------------------------------------------------------------------------------

    // Now we are ready to split the MST into a set of line components.

    // Put the initial starting halfedge to the queue
    queue< Halfedge_handle > queueHH;
    // Insert a dummy half edge, either previous and next to the current edge is OK.
    queueHH.push( startHH->prev() );

    // Line components to be traced
    vector< Halfedge_handle > stripHH;
    // Array of branching halfedges
    vector< Halfedge_handle > branchHH;

    while ( ! queueHH.empty() ) {
	Halfedge_handle curHH = queueHH.front();
	queueHH.pop();
	stripHH.clear();
	stripHH.push_back( curHH );
	
#ifdef MYDEBUG
	cerr << " base MST edge : " 
	     << curHH->opposite()->vertex()->id() << " == " 
	     << curHH->vertex()->id() << endl;
	cerr << " base MST dual edge : "
	     << curHH->facet()->id() << " == " 
	     << curHH->opposite()->facet()->id() << endl;
#endif	// MYDEBUG
	    
	while ( true ) {

#ifdef MYDEBUG
	    cerr << "------------------------------" << endl;
	    cerr << " current MST edge : "
		 << curHH->opposite()->vertex()->id() << " == "
		 << curHH->vertex()->id() << endl;
	    cerr << " current MST dual edge : "
		 << curHH->facet()->id() << " == "
		 << curHH->opposite()->facet()->id() << endl;
#endif	// MYDEBUG

	    // If a next connected half edge exists,
	    incTreeEdges( curHH, branchHH );

#ifdef MYDEBUG
	    cerr << " number of next possible MST dual edges = " << branchHH.size() << endl;
#endif	// MYDEBUG

	    // In case of a straight dual line
	    if ( ( int )branchHH.size() == 1 ) {
		// Update the current half edge   
		curHH = branchHH[ 0 ]->opposite();
		stripHH.push_back( curHH );
		// Disable the connection flag of the current half edge for the search in the MST
		curHH->opposite()->connect() = curHH->connect() = false;
#ifdef MYDEBUG
		cerr << " 1st next MST dual edge (now disconnected): "
		     << branchHH[ 0 ]->opposite()->facet()->id() << " == "
		     << branchHH[ 0 ]->facet()->id() << endl;
#endif	// MYDEBUG
	    }
	    // In case of branching dual lines 
	    else if ( ( int )branchHH.size() == 2 ) {
		Halfedge_handle storeHH;
#ifdef WEIGHT_BASED_SPLIT
		if ( branchHH[ 0 ]->weight() < branchHH[ 1 ]->weight() ) {
		    // Update the current half edge   
		    curHH = branchHH[ 0 ]->opposite();
		    storeHH = branchHH[ 1 ]->opposite();
		}
		else {
		    // Update the current half edge   
		    curHH = branchHH[ 1 ]->opposite();
		    storeHH = branchHH[ 0 ]->opposite();
		}
#else	// WEIGHT_BASED_SPLIT
		double length[ 2 ];
		for ( unsigned int k = 0; k < 2; ++k ) {
		    Vector3 vec = branchHH[k]->vertex()->point() - branchHH[k]->opposite()->vertex()->point();
		    length[ k ] = sqrt( vec.squared_length() );
		}
		if ( length[ 0 ] > length[ 1 ] ) {
		    // Update the current half edge   
		    curHH = branchHH[ 0 ]->opposite();
		    storeHH = branchHH[ 1 ]->opposite();
		}
		else {
		    // Update the current half edge   
		    curHH = branchHH[ 1 ]->opposite();
		    storeHH = branchHH[ 0 ]->opposite();
		}
#endif	// WEIGHT_BASED_SPLIT
		stripHH.push_back( curHH );
		// Disable the connection flag of the current half edge for the search in the MST
		curHH->opposite()->connect() = curHH->connect() = false;

		queueHH.push( storeHH );
		stitch.push_back( storeHH );
		// Disable the connection flag of the branching half edge for the search in the MST
		storeHH->connect() = storeHH->opposite()->connect() = false;
	    }
	    else {
#ifdef MYDEBUG
		cerr << "##############################" << endl;
		cerr << "Here is a deadend of the MST." << endl;
		cerr << "##############################" << endl;
#endif	// MYDEBUG
		// If the MST ends at this face, go out of the loop
		break;
	    }
	}

#ifdef MYDEBUG
	cerr << "##############################" << endl;
	for ( unsigned int i = 0; i < stripVHH.size(); ++i ) {
	    cerr << "[ " << setw( 3 ) << i << " ] = halfedge ID: " << setw( 5 ) << stripVHH[ i ]->id() << endl;
	}
#endif	// MYDEBUG

	// Numbering the current sequence of half edges
	for ( unsigned int i = 0; i < stripHH.size(); i++ ) {
	    stripHH[ i ]->cycle() = stripHH[ i ]->opposite()->cycle() = stripID;
	}

#ifdef MYDEBUG
	for ( unsigned int i = 0; i < stripVHH.size(); i++ ) {
	    cerr << " i = " << i;
	    if ( (stripVHH[ i ]->connect()) && (stripVHH[ i ]->opposite()->connect()) )
		cerr << " both connected" << endl;
	    else if ( (!stripVHH[ i ]->connect()) && (!stripVHH[ i ]->opposite()->connect()) ) 
		cerr << " both disconnected" << endl;
	    else
		cerr << " illegal status" << endl;
	}
#endif	// MYDEBUG

	spine.push_back( stripHH );
	stripID++;
    }

    //------------------------------------------------------------------------------
    //	Recover the connection flags of half edges for projecting cycles
    //------------------------------------------------------------------------------
    for ( unsigned int k = 0; k < spine.size(); ++k ) {
	// starting from 1 to size
	for ( unsigned int m = 1; m < spine[ k ].size(); ++m ) {
	    spine[ k ][ m ]->connect() = spine[ k ][ m ]->opposite()->connect() = true;
	}
    }

    cerr << "*** Number of edges to be stitched = " << stitch.size() << endl;
}

void MSTtoBranches( Polyhedron & poly,
		    vector< vector< Halfedge_handle > > & spine,
		    vector< Halfedge_handle > & stitch )
{
    int stripID = 0;
    // Initialize the input array of triangular strips
    spine.clear();
    stitch.clear();

    // Find one endpoint of the MST as the starting point for tracing

    // flag for the existence of a triangular strip in MST.
    bool isExist = false;
    Halfedge_handle startHH;
    
    // for each dual edge in the MST
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	if ( hi->connect() ) {
	    // Extracting the deadend halfedge of the minimum spanning tree
	    // (There must be at least 2.)
	    vector< Halfedge_handle > srcVHH, tarVHH; 
	    incTreeEdges( hi, srcVHH );
	    incTreeEdges( hi->opposite(), tarVHH );
	    if ( srcVHH.size() == 1 ) {
		// we can start a new trianuglar strip from that dual vertex
		isExist = true;
		startHH = hi;
		break;
	    }
	    else if ( tarVHH.size() == 1 ) {
		isExist = true;
		startHH = hi->opposite();
		break;
	    }
	}
    }

    //------------------------------------------------------------------------------
    //	If we cannot find the endpoints of the MST...
    //------------------------------------------------------------------------------
    if ( ! isExist ) {
	cerr << "We cannot find the endpoints of the MST!!" << endl;
	assert( isExist );
    }
    
    //------------------------------------------------------------------------------
    //	If we can successfully find one endpoint of the MST...
    //------------------------------------------------------------------------------

#ifdef MYDEBUG
    cerr << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    cerr << " starting MST edge : " << startHH->opposite()->vertex()->id() << " == " 
	 << startHH->vertex()->id() << endl;
    cerr << " starting MST dual edge : " << startHH->facet()->id() << " == " 
	 << startHH->opposite()->facet()->id() << endl;
#endif	// MYDEBUG

    //------------------------------------------------------------------------------
    //	Start to track a triangular strip on the remaing part of the MST
    //------------------------------------------------------------------------------
    // Now we are ready to split the MST into a set of line components.

    // Put the initial starting halfedge to the queue
    queue< Halfedge_handle > queueHH;
    // Insert a dummy half edge, either previous and next to the current edge is OK.
    queueHH.push( startHH->prev() );

    // Line components to be traced
    vector< Halfedge_handle > stripHH;
    // Array of branching halfedges
    vector< Halfedge_handle > branchHH;

    while ( ! queueHH.empty() ) {
	Halfedge_handle curHH = queueHH.front();
	queueHH.pop();
	stripHH.clear();
	stripHH.push_back( curHH );
	
#ifdef MYDEBUG
	cerr << " base MST edge : " 
	     << curHH->opposite()->vertex()->id() << " == " 
	     << curHH->vertex()->id() << endl;
	cerr << " base MST dual edge : "
	     << curHH->facet()->id() << " == " 
	     << curHH->opposite()->facet()->id() << endl;
#endif	// MYDEBUG

	while ( true ) {

#ifdef MYDEBUG
	    cerr << "------------------------------" << endl;
	    cerr << " current MST edge : "
		 << curHH->opposite()->vertex()->id() << " == "
		 << curHH->vertex()->id() << endl;
	    cerr << " current MST dual edge : "
		 << curHH->facet()->id() << " == "
		 << curHH->opposite()->facet()->id() << endl;
#endif	// MYDEBUG

	    // If a next connected half edge exists,
	    incTreeEdges( curHH, branchHH );
	    
#ifdef MYDEBUG
	    cerr << " number of next possible MST dual edges = " << branchHH.size() << endl;
#endif	// MYDEBUG

	    // In case of a straight dual line
	    if ( ( int )branchHH.size() == 1 ) {
		// Update the current half edge   
		curHH = branchHH[ 0 ]->opposite();
		stripHH.push_back( curHH );
		// Disable the connection flag of the current half edge for the search in the MST
		curHH->opposite()->connect() = curHH->connect() = false;
#ifdef MYDEBUG
		cerr << " next MST dual edge (now disconnected): "
		     << branchHH[ 0 ]->opposite()->facet()->id() << " == "
		     << branchHH[ 0 ]->facet()->id() << endl;
#endif	// MYDEBUG
	    }
	    // In case of branching dual lines 
	    else if ( ( int )branchHH.size() == 2 ) {
		queueHH.push( branchHH[ 0 ]->opposite() );
		stitch.push_back( branchHH[ 0 ]->opposite() );
		branchHH[ 0 ]->opposite()->connect() = branchHH[ 0 ]->connect() = false;
		queueHH.push( branchHH[ 1 ]->opposite() );
		stitch.push_back( branchHH[ 1 ]->opposite() );
		branchHH[ 1 ]->opposite()->connect() = branchHH[ 1 ]->connect() = false;
#ifdef MYDEBUG
		cerr << " 1st next MST dual edge (now disconnected): "
		     << branchHH[ 0 ]->opposite()->facet()->id() << " == "
		     << branchHH[ 0 ]->facet()->id() << endl;
		cerr << " 2nd next MST dual edge (now disconnected): "
		     << branchHH[ 1 ]->opposite()->facet()->id() << " == "
		     << branchHH[ 1 ]->facet()->id() << endl;
#endif	// MYDEBUG

		curHH = NULL;
		// If the MST stops at this face, go out of the loop
		break;
	    }
	    else {
		curHH = NULL;
		// cerr << " the other endpoint" << endl;
		// If the MST ends at this face, go out of the loop
		break;
	    }
	}

	// Numbering the current sequence of half edges
	for ( unsigned int i = 0; i < stripHH.size(); i++ ) {
	    stripHH[ i ]->cycle() = stripHH[ i ]->opposite()->cycle() = stripID;
	}
	spine.push_back( stripHH );
	stripID++;
    }

    //------------------------------------------------------------------------------
    //	Recover the connection flags of half edges for projecting cycles
    //------------------------------------------------------------------------------
    for ( unsigned int k = 0; k < spine.size(); ++k ) {
	// starting from 1 to size
	for ( unsigned int m = 1; m < spine[ k ].size(); ++m ) {
	    spine[ k ][ m ]->connect() = spine[ k ][ m ]->opposite()->connect() = true;
	}
    }

    cerr << "*** Number of edges to be stitched = " << stitch.size() << endl;
}


void MSTtoST( Polyhedron & poly )
{
    double threshold = 0.0;
    multimap< double, Halfedge_handle > mmap;

#ifdef DEBUG
    // Initialization of edge connections
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	Halfedge_handle hh = hi;
	if ( ( hh->connect() ) && 
	     ( hh->opposite()->vertex()->id() < hh->vertex()->id() ) ) {
	    mmap.insert( make_pair< double, Halfedge_handle >( -hh->weight(), hh ) );
	}
    }

    double ave, stdev;
    weightStatistics( poly, ave, stdev );
    
    multimap< double, Halfedge_handle >::iterator it = mmap.begin();
    unsigned int counter = 0;
    while ( it != mmap.end() ) {
	// cerr << " counter = " << counter << " weight = " << it->first << endl;
	if ( it->first > threshold ) {
	    cerr << " it->first = " << it->first << " threshold = " << threshold << endl;
	    break;
	}
	Halfedge_handle hh = it->second;
	hh->connect() = hh->opposite()->connect() = false;
	counter++;
	it++;
    }
    // cerr << " counter = " << counter << endl;
#endif	// DEBUG

    // Initialization of edge connections
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	Halfedge_handle hh = hi;
	if ( ( hh->connect() ) && 
	     ( hh->opposite()->vertex()->id() < hh->vertex()->id() ) ) {
//	    mmap.insert( make_pair< double, Halfedge_handle >( hh->weight(), hh ) );

		mmap.insert( make_pair( hh->weight(), hh ) );
	}
    }

    double ave, stdev;
    weightStatistics( poly, ave, stdev );
    threshold = ave - mstDeviation * stdev;
    
    multimap< double, Halfedge_handle >::iterator it = mmap.end();
    unsigned int counter = 0;
    do {
	it--;
	if ( it->first < threshold ) {
	    cerr << " it->first = " << it->first << " threshold = " << threshold << endl;
	    break;
	}
	Halfedge_handle hh = it->second;
	hh->connect() = hh->opposite()->connect() = false;
	counter++;
    } while ( it != mmap.begin() );
}					      


void STtoStrips( Polyhedron & poly,
		 vector< vector< Halfedge_handle > > & spine,
		 vector< Halfedge_handle > & stitch )
{
    int patchID = 0;
    // Initialize the input array of triangular strips
    spine.clear();
    vector< bool > facetFlag( poly.size_of_facets() );

    for ( unsigned int k = 0; k < facetFlag.size(); ++k ) facetFlag[ k ] = false;

    while ( true ) {
	// Find one endpoint of the ST as the starting point for tracing

	// flag for the existence of a triangular strip in MST.
	bool isExist = false;
	Halfedge_handle startHH;
    
	// for each dual edge in the MST
	for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	    if ( hi->connect() ) {
		// Extracting the deadend halfedge of the minimum spanning tree
		// (There must be at least 2.)
		vector< Halfedge_handle > srcVHH, tarVHH; 
		incTreeEdges( hi, srcVHH );
		incTreeEdges( hi->opposite(), tarVHH );
		if ( srcVHH.size() == 1 ) {
		    // we can start a new trianuglar strip from that dual vertex
		    isExist = true;
		    startHH = hi;
		}
		else if ( tarVHH.size() == 1 ) {
		    isExist = true;
		    startHH = hi->opposite();
		}
	    }

	    if ( isExist ) break;
	}

	//------------------------------------------------------------------------------
	//	If we cannot find the endpoints of the ST,
	//	get out of the loop
	//------------------------------------------------------------------------------
	if ( ! isExist ) break;
    
	//------------------------------------------------------------------------------
	//	If we can successfully find one endpoint of the ST,
	//	start to track the corresponding spanning tree
	//------------------------------------------------------------------------------
	
	// Now we are ready to track the ST

	// Put the initial starting halfedge to the queue
	queue< Halfedge_handle > queueHH;
	// Insert a dummy half edge, either previous and next to the current edge is OK.
	queueHH.push( startHH->prev() );

	// Line components to be traced
	vector< Halfedge_handle > stripHH;
	// Array of branching halfedges
	vector< Halfedge_handle > branchHH;

	while ( ! queueHH.empty() ) {
	    Halfedge_handle curHH = queueHH.front();
	    queueHH.pop();
	    stripHH.clear();
	    stripHH.push_back( curHH );
	
	    while ( true ) {
		// If a next connected half edge exists,
		incTreeEdges( curHH, branchHH );

		// In case of a straight dual line
		if ( ( int )branchHH.size() == 1 ) {
		    // Update the current half edge   
		    curHH = branchHH[ 0 ]->opposite();
		    stripHH.push_back( curHH );
		    // Disable the connection flag of the current half edge for the search in the MST
		    curHH->opposite()->connect() = curHH->connect() = false;
		}
		// In case of branching dual lines 
		else if ( ( int )branchHH.size() == 2 ) {
		    Halfedge_handle storeHH;
#ifdef WEIGHT_BASED_SPLIT
		    if ( branchHH[ 0 ]->weight() < branchHH[ 1 ]->weight() ) {
			// Update the current half edge   
			curHH = branchHH[ 0 ]->opposite();
			storeHH = branchHH[ 1 ]->opposite();
		    }
		    else {
			// Update the current half edge   
			curHH = branchHH[ 1 ]->opposite();
			storeHH = branchHH[ 0 ]->opposite();
		}
#else	// WEIGHT_BASED_SPLIT
		    double length[ 2 ];
		    for ( unsigned int k = 0; k < 2; ++k ) {
			Vector3 vec = branchHH[k]->vertex()->point() - branchHH[k]->opposite()->vertex()->point();
			length[ k ] = sqrt( vec.squared_length() );
		    }
		    if ( length[ 0 ] > length[ 1 ] ) {
			// Update the current half edge   
			curHH = branchHH[ 0 ]->opposite();
			storeHH = branchHH[ 1 ]->opposite();
		    }
		    else {
			// Update the current half edge   
			curHH = branchHH[ 1 ]->opposite();
			storeHH = branchHH[ 0 ]->opposite();
		    }
#endif	// WEIGHT_BASED_SPLIT
		    stripHH.push_back( curHH );
		    // Disable the connection flag of the current half edge for the search in the MST
		    curHH->opposite()->connect() = curHH->connect() = false;
		    
		    queueHH.push( storeHH );
		    stitch.push_back( storeHH );
		    // Disable the connection flag of the branching half edge for the search in the MST
		    storeHH->connect() = storeHH->opposite()->connect() = false;
		}
		else {
		    // If the MST ends at this face, go out of the loop
		    break;
		}
	    }

	    // Numbering the current sequence of half edges
	    for ( unsigned int i = 0; i < stripHH.size(); i++ ) {
		stripHH[ i ]->cycle() = stripHH[ i ]->opposite()->cycle() = patchID;
		// Caution !!: Visited faces are marked here
		facetFlag[ stripHH[ i ]->facet()->id() ] = true;
	    }
	    
	    spine.push_back( stripHH );
	    patchID++;
	    
	    // cerr << " Still " << setw( 3 ) << queueHH.size() << " remaining starting points" << endl;
	}
    }
	      
    //------------------------------------------------------------------------------
    //	Recover the connection flags of half edges for projecting cycles
    //------------------------------------------------------------------------------
    for ( unsigned int k = 0; k < spine.size(); ++k ) {
	// starting from 1 to size
	for ( unsigned int m = 1; m < spine[ k ].size(); ++m ) {
	    spine[ k ][ m ]->connect() = spine[ k ][ m ]->opposite()->connect() = true;
	}
    }


    //------------------------------------------------------------------------------
    //	Caution!!: Isolated faces are contained in this case.
    //------------------------------------------------------------------------------
    for ( Facet_iterator fi = poly.facets_begin(); fi != poly.facets_end(); ++fi ) {
	Facet_handle fh = fi;
	// If the face is not visited yet,
	if ( ! facetFlag[ fh->id() ] ) {


	    // bool isIsolated = true;
	    // Halfedge_facet_circulator hfc = fi->facet_begin();
	    // do {
	    // Halfedge_handle hh = hfc;
	    // if ( hh->cycle() != NO_INDEX ) isIsolated = false;
	    // if ( ! isIsolated ) break;
	    // } while ( ++hfc != fi->facet_begin() );

	    // if ( isIsolated ) {

	    vector< Halfedge_handle > single;
	    // cerr << " Facet ID " << setw( 5 ) << fh->id() << " is isolated." << endl;
	    // Any halfedge among three is OK.
	    single.push_back( fh->halfedge() );
	    spine.push_back( single );
	}
    }
}



//------------------------------------------------------------------------------
//	Unfolding the 3D mesh into a set of MST strips
//------------------------------------------------------------------------------
void MSTSplit( Graph & tree, Polyhedron & poly,
	       vector< vector< Halfedge_handle > > & spine,
	       vector< Halfedge_handle > & joint )
{
    calcMST( tree, poly );
    MSTtoStrips( poly, spine, joint );
}


//------------------------------------------------------------------------------
//	Construct the overall MST unfolded patterns
//------------------------------------------------------------------------------
void MSTMerge( Polyhedron & poly,
	       vector< Halfedge_handle > & joint,
	       vector< vector< Facet_handle > > & patch,
	       Bbox2 & sheet, vector< Bbox2 > & bound )
{
    vector< Triangle2 > small;
    for ( unsigned int k = 0; k < joint.size(); ++k ) {
#ifdef MYDEBUG
	cerr << " Joint at k = " << k << endl;
#endif	// MYDEBUG
	small.clear();
	transformStrip( joint[ k ], patch, small );
	stitch( joint[ k ], patch, small );
	resetCycleIDs( patch );
    }
    arrangeSingle( patch, sheet, bound );
}



//------------------------------------------------------------------------------
//	Resolve intersections on the MST-based unfold patch
//------------------------------------------------------------------------------
void MSTDissolve( Polyhedron & poly,
		  vector< vector< Facet_handle > > & patch,
		  Bbox2 & sheet, vector< Bbox2 > & bound )
{
    vector< pair< Facet_handle, Facet_handle > > endfaces;
    vector< Facet_handle > path;
    vector< vector< Facet_handle > > pathlist;
    vector< Halfedge_handle > cutList;
    // vector< Halfedge_handle > joint;

    vector< Facet_handle > allfaces;
    for ( unsigned int i = 0; i < patch.size(); ++i ) 
	for ( unsigned int j = 0; j < patch[ i ].size(); ++j )
	    allfaces.push_back( patch[ i ][ j ] );
    endfaces.clear();
    findIntersections( allfaces, endfaces );
    allfaces.clear();

    pathlist.clear();
    for ( unsigned int k = 0; k < endfaces.size(); ++k ) {
	if ( findDualPath( endfaces[ k ].first, endfaces[ k ].second,
			   poly, path ) ) {
	    pathlist.push_back( path );
	}
    }

    minimumCover( poly, endfaces, pathlist, cutList );

    for ( unsigned int m = 0; m < cutList.size(); ++m ) {
	cutStrips( patch, cutList[ m ] );
    }
    arrangeSingle( patch, sheet, bound );

#ifdef MYDEBUG
    cerr << "%%%%% paper size = "
	 << ( sheet.xmax() - sheet.xmin() ) * ( sheet.ymax() - sheet.ymin() )
	 << "%%%%%" << endl;
#else	// MYDEBUG
    cerr << "%%%%% paper width = "
	 << ( sheet.xmax() - sheet.xmin() ) << "%%%%%" << endl;
#endif	// MYDEBUG
}


void MSTStrips( Graph & tree, Polyhedron & poly,
		 vector< vector< Halfedge_handle > > & spine,
		 vector< vector< Facet_handle > > & patch,
		 double limitCos,
		 Bbox2 & sheet, vector< Bbox2 > & bound )
{
    vector< pair< Facet_handle, Facet_handle > > endfaces;
    vector< Facet_handle > path;
    vector< vector< Facet_handle > > pathlist;
    vector< Halfedge_handle > cutList;
    vector< Halfedge_handle > joint;

    calcMST( tree, poly );
    MSTtoStrips( poly, spine, joint );
    //MSTtoBranches( poly, spine, joint );
    projectCycles( spine, patch, sheet, bound, limitCos );

    resolveOverlaps( poly, patch );

    arrangeSingle( patch, sheet, bound );
}


void MSTPieces( Graph & tree, Polyhedron & poly,
		 vector< vector< Halfedge_handle > > & spine,
		 vector< vector< Facet_handle > > & patch,
		 double limitCos,
		 Bbox2 & sheet, vector< Bbox2 > & bound )
{
    vector< pair< Facet_handle, Facet_handle > > endfaces;
    vector< Facet_handle > path;
    vector< vector< Facet_handle > > pathlist;
    vector< Halfedge_handle > cutList;
    vector< Halfedge_handle > joint;

    calcMST( tree, poly );
    // MSTtoStrips( poly, spine, joint );
    MSTtoBranches( poly, spine, joint );
    projectCycles( spine, patch, sheet, bound, limitCos );

    resolveOverlaps( poly, patch );

    arrangeSingle( patch, sheet, bound );
}


//------------------------------------------------------------------------------
//	Minimum spanning tree
//------------------------------------------------------------------------------
void minimumSpanningTree( Graph & tree, Polyhedron & poly,
			  vector< vector< Facet_handle > > & patch,
			  Bbox2 & sheet, vector< Bbox2 > & bound )
{
    vector< Halfedge_handle > joint;
    vector< vector< Halfedge_handle > >	skeleton;

    MSTSplit		( partial, mesh, skeleton, joint );
    double limitCos	= INNER_PROD_ALL;
    projectCycles	( skeleton, pattern, paper, bound, limitCos );
    MSTMerge		( mesh, joint, pattern, paper, bound );
    MSTDissolve		( mesh, pattern, paper, bound );
    arrangeMultiple	( pattern, paper, bound );
}


//------------------------------------------------------------------------------
//	Retrieve the spanning trees using the input threshold for the edge
//	weights
//------------------------------------------------------------------------------
void spanningTrees( Graph & tree, Polyhedron & poly,
		    vector< vector< Facet_handle > > & patch,
		    Bbox2 & sheet, vector< Bbox2 > & bound )
{
    vector< pair< Facet_handle, Facet_handle > > endfaces;
    vector< Facet_handle > path;
    vector< vector< Facet_handle > > pathlist;
    vector< Halfedge_handle > cutList;
    vector< Halfedge_handle > joint;
    vector< vector< Halfedge_handle > >	skeleton;

    calcMST( tree, poly );
    MSTtoST( poly );

    STtoStrips( poly, skeleton, joint );

#ifdef MYDEBUG
    for ( unsigned int i = 0; i < spine.size(); ++i ) {
	cerr << "[ " << setw( 3 ) << i << " ] Num = " << setw( 3 ) << spine[ i ].size() << endl;
	for ( unsigned int j = 0; j < spine[ i ].size(); ++j ) {
	    cerr << "[ " << setw( 3 ) << i << " ][ = " << setw( 3 ) << j << " ] : " << ends;
	    cerr << " Face ID = " << setw( 3 ) << spine[ i ][ j ]->facet()->id() << " | ";
	    printHalfedge( spine[ i ][ j ] );
	}
    }
#endif	// MYDEBUG

    double limitCos = INNER_PROD_ALL;
    projectCycles	( skeleton, patch, sheet, bound, limitCos );

    MSTMerge		( poly, joint, patch, sheet, bound );

    resolveOverlaps	( poly, patch );

    arrangeSingle	( patch, sheet, bound );
}



