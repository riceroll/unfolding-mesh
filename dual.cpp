//------------------------------------------------------------------------------
//
//	dual.cpp
//
//------------------------------------------------------------------------------
#include "common.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

#include <iostream>
#include <iomanip>

#ifndef	SMOOTHING_LAMBDA
#define	SMOOTHING_LAMBDA	(0.6307)
// #define	_SMOOTHING_LAMBDA	(0.6000)
// #define	_SMOOTHING_LAMBDA	(0.5000)
// #define	_SMOOTHING_LAMBDA	(0.3000)
// #define	_SMOOTHING_LAMBDA	(0.8000)
#endif // _SMOOTHING_LAMBDA

#ifndef	SMOOTHING_MU
#define	SMOOTHING_MU		(-0.673155945)
// #define	_SMOOTHING_MU		(_SMOOTHING_LAMBDA/(0.1*_SMOOTHING_LAMBDA-1.0))
// #define	_SMOOTHING_MU		(0.0)
#endif // SMOOTHING_MU


void printHalfedge( const Halfedge_handle & hh )
{
    cerr << setw( 3 ) << hh->opposite()->vertex()->id() 
	 << " == " << setw( 3 ) << hh->vertex()->id() << endl;
}

//------------------------------------------------------------------------------
//	Normalize the 3D triangulated mesh with the display window
//------------------------------------------------------------------------------
void normalizeMesh( Polyhedron & poly, vector< Segment3 > & bone )
{
    Vector3 sum, ave;

    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	sum = sum + ( vi->point() - CGAL::ORIGIN );
    }
    ave = sum / ( double )poly.size_of_vertices();

    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	vi->point() = vi->point() - ave;
    }
    
    cerr << " ave = " << ave << endl;
    Transformation3 translate( CGAL::TRANSLATION, -ave );

    // fabs:absolute values-no negative values
    double sideMax = 0.0;
    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	if ( fabs( vi->point().x() ) > sideMax ) sideMax = fabs( vi->point().x() );
	if ( fabs( vi->point().y() ) > sideMax ) sideMax = fabs( vi->point().y() );
	if ( fabs( vi->point().z() ) > sideMax ) sideMax = fabs( vi->point().z() );
    }
    
	// sideMax: the largest number
    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	vi->point() = CGAL::ORIGIN + ( vi->point() - CGAL::ORIGIN ) / sideMax;
    }

    Transformation3 scale( CGAL::SCALING, 1.0/sideMax );

    Transformation3 composite = scale * translate;
	    
    for ( unsigned int k = 0; k < bone.size(); ++k ) {
	bone[ k ] = bone[ k ].transform( composite );
    }
}

//------------------------------------------------------------------------------
//	Initialize the dual center 
//------------------------------------------------------------------------------
void initAttrs( Polyhedron & poly )
{
    // assign IDs to vertices
    int vID = 0;
    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	vi->id()	= vID++;
	vi->label()	= DEFAULT_LABEL;
	vi->nCuts()	= 0;
    }
    cerr << "Total number of vertices = " << vID << endl;

    // assign IDs to faces: the dual vertex and the dual coordinates-the center coordinates
    int fID = 0;
    for ( Facet_iterator fi = poly.facets_begin(); fi != poly.facets_end(); ++fi ) {
	fi->id()	= fID++;
	// fi->label()	= DEFAULT_LABEL;
	fi->piece()	= NO_INDEX;
	Vector3 sum( 0.0, 0.0, 0.0 );
	Halfedge_facet_circulator hfc = fi->facet_begin();
	do {
	    sum = sum + ( hfc->vertex()->point() - CGAL::ORIGIN );
	} while ( ++hfc != fi->facet_begin() );
	sum = sum / ( double )CGAL::circulator_size( hfc );
	fi->center() = CGAL::ORIGIN + sum;
	// cerr << " Barycenter No. " << fid-1 << " : " << bary[ fid-1 ] << endl;
    }
    cerr << "Total number of facets = " << fID << endl;

    // initialize the halfedge IDs
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	hi->id() = NO_INDEX;
    }

    // assign the same IDs to the identical/ opposite halfedges
    int eID = 0;
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	if ( hi->id() == NO_INDEX ) {
	    assert( hi->opposite()->id() == NO_INDEX );
	    hi->id() = eID;
	    hi->opposite()->id() = eID;
	    eID++;
	    hi->label()		= hi->opposite()->label()	= DEFAULT_LABEL;
	    hi->match()		= hi->opposite()->match()	= DEFAULT_LABEL;
	    hi->cycle()		= hi->opposite()->cycle()	= NO_INDEX;
	    hi->weight()	= hi->opposite()->weight()	= 0.0;
	    hi->connect()	= hi->opposite()->connect()	= true;
	    hi->visit()		= hi->opposite()->visit()	= false;
	    // We individually handle hi->fixed
	}
    }
    cerr << "Total number of edges = " << eID << endl;
}


void renewAttrs( Polyhedron & poly )
{
    // assign IDs to vertices
    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	vi->nCuts()	= 0;
    }

    // assign IDs to faces: the dual vertex and the dual coordinates-the center coordinates
    for ( Facet_iterator fi = poly.facets_begin(); fi != poly.facets_end(); ++fi ) {
	fi->piece()	= NO_INDEX;
	Vector3 sum( 0.0, 0.0, 0.0 );
	Halfedge_facet_circulator hfc = fi->facet_begin();
	do {
	    sum = sum + ( hfc->vertex()->point() - CGAL::ORIGIN );
	} while ( ++hfc != fi->facet_begin() );
	sum = sum / ( double )CGAL::circulator_size( hfc );
	fi->center() = CGAL::ORIGIN + sum;
    }

    // assign the same IDs to the identical/ opposite halfedges
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	hi->label()		= DEFAULT_LABEL;
	hi->match()		= DEFAULT_LABEL;
	hi->cycle()		= NO_INDEX;
	hi->connect()		= true;
	hi->visit()		= false;
	hi->path()		= NO_INDEX;
	hi->orient()		= true;
	// We individually handle hi->fixed
    }
}


void initWeights( Polyhedron & poly )
{
    // assign the same IDs to the identical/ opposite halfedges
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	hi->fixed()		= false;
    }
}


//------------------------------------------------------------------------------
//	Compute and store the results of dual graphs into appropriate variables
//
//  computeDual( mesh, whole, bary, mid ); 
//------------------------------------------------------------------------------
void computeDual( Polyhedron & poly, Graph & dual,
		  vector< Point3 > & bary, vector< Point3 > & mid )
{
    // initialize the array of face barycenters: refer to the point at which the gravitational forces exerted by 2 objects are equal

    bary.clear();
    bary.resize( poly.size_of_facets() );

    // initialize the dual graph with dual vertices
    dual.clear();
    int fid = 0;
    for ( Facet_iterator fi = poly.facets_begin(); fi != poly.facets_end(); ++fi ) {
	// add a vertex: each face equals to each dual vertex
	add_vertex( dual );
	bary[ fid++ ] = fi->center(); 
	// cerr << " Barycenter No. " << fid-1 << " : " << bary[ fid-1 ] << endl;
    }

    int size_of_edges = poly.size_of_halfedges() / 2; //redundant

    // initialize the array of midpoints
    mid.clear();
    mid.resize( size_of_edges );

    // initialize the array of lengths
    // length.clear();
    // length.resize( size_of_edges );

    // construct the connecitivity of the dual graph
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	int origID = hi->facet()->id(); // the face
	int destID = hi->opposite()->facet()->id(); // the neighbor face
	Point3 origCoord = hi->vertex()->point();
	Point3 destCoord = hi->opposite()->vertex()->point();
	Vector3 dispVec = origCoord - destCoord;
	Point3 midCoord = destCoord + 0.5 * dispVec;
	// double curLength = sqrt( dispVec.squared_length() ); // the length between two connected vertices
#ifdef DEBUG
	cerr << origID << " -- " << destID << endl;
#endif	// DEBUG	
	
	if ( origID < destID )
		add_edge( origID, destID, hi->id(), dual );  
		
	mid[ hi->id() ] = midCoord;
	hi->mid() = hi->opposite()->mid() = midCoord;
	hi->weight() = hi->opposite()->weight() = 1.0;
    }
}

//------------------------------------------------------------------------------
//	Align the mesh
//------------------------------------------------------------------------------
Vector3 principalAxis( Polyhedron & poly )
{
    int num = poly.size_of_vertices();

    // initialization here is very important!!
    Point3 ave( 0.0, 0.0, 0.0 );
    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	ave = ave + ( vi->point() - CGAL::ORIGIN );
    }
    ave = CGAL::ORIGIN + ( ave - CGAL::ORIGIN )/( double )num;

    unsigned int dim	= 3;
    double * data	= new double [ num*dim ];

    int nPoints = 0;
    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	data[ nPoints * dim + 0 ] = vi->point().x() - ave.x();
	data[ nPoints * dim + 1 ] = vi->point().y() - ave.y();
	data[ nPoints * dim + 2 ] = vi->point().z() - ave.z();
	nPoints++;
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

    gsl_vector * eVal       = gsl_vector_alloc( dim );
    gsl_matrix * eVec       = gsl_matrix_alloc( dim, dim );
    gsl_eigen_symmv_workspace * w
                                = gsl_eigen_symmv_alloc( dim );

    // eigenanalysis of the matrix B
    gsl_eigen_symmv( B, eVal, eVec, w );
    // release the memory of w
    gsl_eigen_symmv_free( w );
    // sort eigenvalues in a descending order
    gsl_eigen_symmv_sort( eVal, eVec, GSL_EIGEN_SORT_VAL_DESC );

#ifdef MYDEBUG
    for ( unsigned int i = 0; i < dim; ++i ) {
	cerr << "Eigenvalue No. " << i << " = " << gsl_vector_get( eVal, i ) << endl;
	cerr << "Eigenvector No. " << i << endl;
	for ( unsigned int j = 0; j < dim; ++j ) {
	    length += gsl_matrix_get( eVec, i, j )*gsl_matrix_get( eVec, i, j );
	}
	cerr << " length = " << length << endl;
    }
#endif	// MYDEBUG

    Vector3 ref( gsl_matrix_get( eVec, 0, 0 ),
		 gsl_matrix_get( eVec, 0, 1 ),
		 gsl_matrix_get( eVec, 0, 2 ) );
    return ref;

#ifdef DEBUG
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
#endif	// DEBUG
}


void transformMesh( Polyhedron & poly, Transformation3 & map )
{
    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	Point3 store = vi->point().transform( map );
	vi->point() = store;
    }
}


void rotateMesh( Polyhedron & poly, char axis, double radius )
{
    switch ( axis ) {
      case 'x':
	  {
	      Transformation3 map( 1.0,		0.0,		0.0,
				   0.0,		cos(radius),	-sin(radius),
				   0.0,		sin(radius),	cos(radius) );
	      transformMesh( poly, map );
	  }
	  break;
      case 'y':
	  {
	      Transformation3 map( cos(radius),	0.0,		sin(radius),
				   0.0,		1.0,		0.0,
				   -sin(radius),	0.0,		cos(radius) );
	      transformMesh( poly, map );
	  }
	  break;
      case 'z':
	  {
	      Transformation3 map( cos(radius),	-sin(radius),	0.0,
				   sin(radius),	cos(radius),	0.0,
				   0.0,		0.0,		1.0 );
	      transformMesh( poly, map );
	  }
	  break;
      default:
	  cerr << "Illegal axis : " << axis << endl;
	  break;
    }
}


void alignMesh( Polyhedron & poly )
{
    int num = poly.size_of_vertices();

    // initialization here is very important!!
    Point3 ave( 0.0, 0.0, 0.0 );
    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	ave = ave + ( vi->point() - CGAL::ORIGIN );
    }
    ave = CGAL::ORIGIN + ( ave - CGAL::ORIGIN )/( double )num;

    unsigned int dim	= 3;
    double * data	= new double [ num*dim ];

    int nPoints = 0;
    for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	data[ nPoints * dim + 0 ] = vi->point().x() - ave.x();
	data[ nPoints * dim + 1 ] = vi->point().y() - ave.y();
	data[ nPoints * dim + 2 ] = vi->point().z() - ave.z();
	nPoints++;
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

    gsl_vector * eVal       = gsl_vector_alloc( dim );
    gsl_matrix * eVec       = gsl_matrix_alloc( dim, dim );
    gsl_eigen_symmv_workspace * w
                                = gsl_eigen_symmv_alloc( dim );

    // eigenanalysis of the matrix B
    gsl_eigen_symmv( B, eVal, eVec, w );
    // release the memory of w
    gsl_eigen_symmv_free( w );
    // sort eigenvalues in a descending order
    gsl_eigen_symmv_sort( eVal, eVec, GSL_EIGEN_SORT_VAL_DESC );

    // #ifdef MYDEBUG
    for ( unsigned int i = 0; i < dim; ++i ) {
	cerr << "Eigenvalue No. " << i << " = " << gsl_vector_get( eVal, i ) << endl;
	cerr << "Eigenvector No. " << i << endl;
	double length = 0.0;
	for ( unsigned int j = 0; j < dim; ++j ) {
	    cerr << gsl_matrix_get( eVec, i, j ) << " ";
	    length += gsl_matrix_get( eVec, i, j )*gsl_matrix_get( eVec, i, j );
	}
	cerr << " length = " << length << endl;
    }
    // #endif	// MYDEBUG

    // 0 1 2, 0 2 1, 
    

    Transformation3 map;
    map = 
	Transformation3( gsl_matrix_get(eVec,1,0), gsl_matrix_get(eVec,0,0), gsl_matrix_get(eVec,2,0),
			 gsl_matrix_get(eVec,1,1), gsl_matrix_get(eVec,0,1), gsl_matrix_get(eVec,2,1),
			 gsl_matrix_get(eVec,1,2), gsl_matrix_get(eVec,0,2), gsl_matrix_get(eVec,2,2) );
    if ( map.is_odd() ) {
	cerr << " Transformation matrix reflected" << endl;
	map = 
	    Transformation3( gsl_matrix_get(eVec,1,0), gsl_matrix_get(eVec,0,0), -gsl_matrix_get(eVec,2,0),
			     gsl_matrix_get(eVec,1,1), gsl_matrix_get(eVec,0,1), -gsl_matrix_get(eVec,2,1),
			     gsl_matrix_get(eVec,1,2), gsl_matrix_get(eVec,0,2), -gsl_matrix_get(eVec,2,2) );
    }

    for ( unsigned int i = 0; i < dim; ++i ) {
	cerr << "| ";
	for ( unsigned int j = 0; j < dim; ++j ) {
	    cerr << map.cartesian( i, j ) << " ";
	}
	cerr << "|" << endl;
    }

    transformMesh( poly, map );

    return;
}

void moveVertex( Vertex_handle & vh, Point3 & coord, double coef )
{
    // Initialization is very important.
    Vector3 average( 0.0, 0.0, 0.0 );
    int degree = 0;
    Halfedge_vertex_circulator hvc = vh->vertex_begin();
    do {
	average = average + ( hvc->opposite()->vertex()->point() - CGAL::ORIGIN );
	degree++;
#ifdef MYDEBUG
	cerr << "point = " << hvc->opposite()->vertex()->point() << endl;
	cerr << "average = " << average << endl;
#endif	// MYDEBUG
    } while ( --hvc != vh->vertex_begin() );

    average = average / ( double )degree;
#ifdef MYDEBUG
    cerr << "Final average = " << average << endl;
#endif	// MYDEBUG

    Vector3 offset = average - ( vh->point() - CGAL::ORIGIN );
#ifdef MYDEBUG
    cerr << "Final offset = " << offset << endl;
#endif	// MYDEBUG

    coord = vh->point() + coef * offset;
}


void smoothMesh( Polyhedron & poly, unsigned int nTimes )
{
    int nV = poly.size_of_vertices();
    const double lambda = SMOOTHING_LAMBDA;
    const double mu	= SMOOTHING_MU;

    vector< Point3 > shrink  ( nV );
    vector< Point3 > expand  ( nV );

    for ( unsigned int k = 0; k < nTimes; ++k ) {
	// copy the vertex coordinates
	for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	    shrink  [ vi->id() ] = vi->point();
	}
	// shrinking stage
	for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	    moveVertex( vi, shrink[ vi->id() ], lambda );
	}
	// copy back the vertex coordinates
	for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	    vi->point()		= shrink[ vi->id() ];
	    expand[ vi->id() ]	= shrink[ vi->id() ];
	}
	// expanding stage
	for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	    moveVertex( vi, expand[ vi->id() ], mu );
	}
	// copy back the vertex coordinates
	for ( Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi ) {
	    vi->point()		= expand[ vi->id() ];
	}
    }
}


