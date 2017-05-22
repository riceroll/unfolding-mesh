//------------------------------------------------------------------------------
//
//	matching.cpp
//
//------------------------------------------------------------------------------
#include "common.h"
#include <boost/graph/max_cardinality_matching.hpp>

//------------------------------------------------------------------------------
//	Decomposition to single faces
//------------------------------------------------------------------------------
void singleFaces( Polyhedron & poly,
		  vector< vector< Halfedge_handle > > & spine )
{
    spine.clear();

    // initialization
    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	vi->nCuts()	= 0;
    }

    int count = 0;
    for ( Facet_iterator fi = poly.facets_begin(); fi != poly.facets_end(); ++fi ) {
	fi->piece()	= count;
	count++;
	vector< Halfedge_handle > each;
	each.push_back( fi->halfedge() );
	spine.push_back( each );
    }

    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	// unsigned int origID = hi->facet()->id();
	// unsigned int destID = hi->opposite()->facet()->id();
	hi->cycle() = NO_INDEX;
	hi->connect() = false;
    }

    countNCuts( poly );
}

//------------------------------------------------------------------------------
//	Perfect matching validation
//------------------------------------------------------------------------------
bool isPerfect( Graph & graph )
{
    int nVertices = num_vertices( graph );
    vector< VertexDescriptor > bridge( nVertices );
    bool success = checked_edmonds_maximum_cardinality_matching( graph, &bridge[0] );
    int nMatchings = matching_size( graph, &bridge[0] );
    return ( success && ( nMatchings == nVertices/2 ) );
}

//------------------------------------------------------------------------------
//	Perfect matching operations
//------------------------------------------------------------------------------
void perfectMatching( Graph & dual, Polyhedron & poly,
		      vector< vector< Halfedge_handle > > & spine )
{
    int nVertices = num_vertices( dual ); // boost:number of vertices in the dual graph
    vector< VertexDescriptor > bridge( nVertices );

    bool success = checked_edmonds_maximum_cardinality_matching( dual, &bridge[0] );
    assert( success );

    cerr << " Found a perfect matching of the size " << matching_size( dual, &bridge[0] ) << endl;
    cerr << " Number of dual vertices " << num_vertices( dual ) << endl;

    // initialization
    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	vi->nCuts()	= 0;
    }
    for ( Facet_iterator fi = poly.facets_begin(); fi != poly.facets_end(); ++fi ) {
	fi->piece()	= NO_INDEX;
    }

    // First identify the edges involved in perfect matching which is in "bridge" variable
	// Perfect matching edges equals to the cutting edges
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	unsigned int origID = hi->facet()->id();
	unsigned int destID = hi->opposite()->facet()->id();
	hi->cycle() = NO_INDEX;
	if ( origID < destID ) {
		// a perfect matching
	    if ( bridge[ origID ] == destID ) {
		hi->match() = hi->opposite()->match() = MATCHED_EDGE;
		// increment the number of cut edges
		hi->vertex()->nCuts()++;
		hi->opposite()->vertex()->nCuts()++;
		hi->connect() = hi->opposite()->connect() = false;
	    }
		// a non-perfect matching
	    else {
		hi->connect() = hi->opposite()->connect() = true;
		hi->match() = hi->opposite()->match() = UNMATCHED_EDGE;
	    }
	}
    }


    // Next trace the DUAL CYCLES on the mesh: the connected edges of non-perfect matching (disjoint cycles)
    int nCycles = 0;
    // group.clear();
    spine.clear(); 
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	if ( ( hi->match() == UNMATCHED_EDGE ) && ( hi->cycle() == NO_INDEX ) ) {
	    vector< Halfedge_handle > temp, circ;
	    vector< Facet_handle > temporaryF, circularF;
	    
	    Halfedge_handle curH = hi->opposite(); //current halfedge
	    Halfedge_handle startH = curH; // start halfedge
	    Facet_handle startF = curH->facet(); // start facet
	    
	    // while ( curH->facet() != startF ) {
		// while ( curH != startH ) {
	    do {
		temp.push_back( curH );
		curH->cycle() = curH->opposite()->cycle() = nCycles;
		temporaryF.push_back( curH->facet() );
		Halfedge_handle halfL = curH->prev();
		Halfedge_handle halfR = curH->next();
		// search only the UNMATCH EDGES
		if ( halfL->match() == UNMATCHED_EDGE ) {
		    curH = halfL->opposite();
		}
		else if ( halfR->match() == UNMATCHED_EDGE ) {
		    curH = halfR->opposite();
		}
		else assert( FALSE );
	    } while ( curH != startH );

	    // search the most sharp junction along the dual cycle
	    double maxAngle = 0.0;
	    unsigned int base = 0;
	    //for ( unsigned int i = 0; i < temporaryF.size(); ++i ) {
	    for ( unsigned int i = 0; i < temp.size(); ++i ) {
		// Point pointP = temporaryF[ i ]->center();
		// Point pointC  = temporaryF[ (i+1)%temporaryF.size() ]->center();
		// Point pointN = temporaryF[ (i+2)%temporaryF.size() ]->center();
		Point3 point0 = temp[         i         ]->facet()->center();
		Point3 point1 = temp[ (i+1)%temp.size() ]->facet()->center();
		Point3 point2 = temp[ (i+2)%temp.size() ]->facet()->center();
		Point3 point3 = temp[ (i+3)%temp.size() ]->facet()->center();
		Vector3 wingP = point0 - point1;
		Vector3 wingN = point3 - point2;
		double lengthP = sqrt( wingP.squared_length() );
		double lengthN = sqrt( wingN.squared_length() );
		double innerProd = ( wingP * wingN ) / ( lengthP*lengthN ); // the angle
		double curAngle = fabs( M_PI - acos( innerProd ) );
		if ( curAngle > maxAngle ) {
#ifdef DEBUG
		    cerr << " maxAngle = " << maxAngle << " curAngle = " << curAngle
			 << " base = " << (i+1)%temporaryF.size() << endl;
#endif	// DEBUG
		    maxAngle = curAngle;
		    base = (i+2)%temporaryF.size();
		}
	    }

		// circ: connected cycle for the unmatch edges
		// spine: consists of all the circs
	    for ( unsigned int i = 0; i < temporaryF.size(); ++i ) {
		circ.push_back( temp[ (i+base)%temp.size() ] );
		circularF.push_back( temporaryF[ (i + base)%temporaryF.size() ] );
	    }
	    spine.push_back( circ );
	    // group.push_back( circularF );
	    nCycles++;
	}
    }
    cerr << " Total number of dual cycles = " << nCycles << endl;
    // cerr << " Total number of dual cycles = " << group.size() << endl;
}



void weightedMatching( Graph & dual, Polyhedron & poly,
		       vector< vector< Halfedge_handle > > & spine )
{
    vector< bool > isVisit( poly.size_of_halfedges() );
    multimap< double, Halfedge_handle > pqueue;

    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	if ( hi->opposite()->vertex()->id() < hi->vertex()->id() ) {
//	    pqueue.insert( make_pair< double, Halfedge_handle >( hi->weight(), hi ) );	error
		pqueue.insert( make_pair( hi->weight(), hi ) );
	}
    }

    multimap< double, Halfedge_handle >::iterator it = pqueue.end();
    cerr << endl;
    while ( it != pqueue.begin() ) {
	it--;
	fixDual( dual, it->second );
#ifdef DEBUG
	if ( fixDual( dual, it->second ) )
	    cerr << "o" << it->first;
	else
	    cerr << "x" << it->first;
#endif	// DEBUG
    }
    cerr << endl;
    perfectMatching( dual, poly, spine );
}
    
