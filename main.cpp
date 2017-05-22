//------------------------------------------------------------------------------
//
//	main.cpp
//
//------------------------------------------------------------------------------
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <cstdio>
#include <cmath>
#include <cassert>

#define __MAC__

#ifdef __MAC__
// OpenCV
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/opencv.hpp>

#endif	// __MAC__

#include "common.h"
#include "gene.h"
#include "ui.h"
#include "timer.h"

using namespace std;

#ifdef __LINUX__
// OpenCV
#include <opencv/cv.h>
#include <opencv/highgui.h>
#endif	// __LINUX__

#ifndef M_PI
#define M_PI	(3.14159265358979323846)	/* pi */
#endif	// M_PI

//#define TINY_WINDOWS
//#define SMALL_WINDOWS
//#define MEDIUM_WINDOWS
#define LAPTOP_WINDOWS
// #define LARGE_WINDOWS
//#define SCREENCAST_WINDOWS

#define ASCII_PPM
// #define BINARY_PPM

// #define DISPLAY_OPTION

// This is necessary for the Linux environment
#define AVOID_DOUBLE_FREE

#ifdef TINY_WINDOWS
#define DEFAULT_MESH_WIDTH	 (256)
#define DEFAULT_MESH_HEIGHT	 (256)
#define DEFAULT_SHEET_WIDTH	 (362)
#define DEFAULT_SHEET_HEIGHT	 (256)
#endif	// TINY_WINDOWS
#ifdef SMALL_WINDOWS
#define DEFAULT_MESH_WIDTH	 (512)
#define DEFAULT_MESH_HEIGHT	 (512)
//#define DEFAULT_SHEET_WIDTH	 (362)
#define DEFAULT_SHEET_WIDTH	 (724)
#define DEFAULT_SHEET_HEIGHT	 (512)
#endif	// SMALL_WINDOWS
#ifdef MEDIUM_WINDOWS
#define DEFAULT_MESH_WIDTH	 (768)
#define DEFAULT_MESH_HEIGHT	 (768)
// #define DEFAULT_SHEET_WIDTH	 (534)
#define DEFAULT_SHEET_WIDTH	(1086)
#define DEFAULT_SHEET_HEIGHT	 (768)
#endif	// MEDIUM_WINDOWS
#ifdef LAPTOP_WINDOWS
// 600 848
#define DEFAULT_MESH_WIDTH	 (600)
#define DEFAULT_MESH_HEIGHT	 (600)
#define DEFAULT_SHEET_WIDTH	 (848)
#define DEFAULT_SHEET_HEIGHT	 (600)
#endif	// LAPTOP_WINDOWS
#ifdef LARGE_WINDOWS
#define DEFAULT_MESH_WIDTH	(1024)
#define DEFAULT_MESH_HEIGHT	(1024)
// #define DEFAULT_SHEET_WIDTH	 (724)
#define DEFAULT_SHEET_WIDTH	(1448)
#define DEFAULT_SHEET_HEIGHT	(1024)
#endif	// LARGE_WINDOWS
#ifdef SCREENCAST_WINDOWS
#define DEFAULT_MESH_WIDTH	(512)
#define DEFAULT_MESH_HEIGHT	(512)
#define DEFAULT_SHEET_WIDTH	(1448)
#define DEFAULT_SHEET_HEIGHT	(1024)
#endif	// SCREENCAST_WINDOWS

#define KEY_ESC		(27)

#define MESH_ONLY			(200)
#define WITH_DUAL			(210)
#define WITH_MATCHING			(220)
#define PIECE_LAYOUT			(230)
#define WITH_MST			(240)
#define MEDIAL_AXIS			(250)
#define CRAFT_RENDERING			(260)
#define PATTERN_RENDERING		(290)

#define SADDLE_ON			(310)
#define SADDLE_OFF			(320)

#define LABEL_ON			(330)
#define LABEL_OFF			(340)

#define ROTATE_XAXIS			(351)
#define ROTATE_YAXIS			(352)
#define ROTATE_ZAXIS			(353)
#define ALIGN_MESH			(355)
#define SMOOTH_MESH			(361)

#define UNIFORM_ASSIGNMENT		(411)
#define RANDOM_ASSIGNMENT		(412)
#define MINPERIMETER_ASSIGNMENT		(421)
#define MAXPERIMETER_ASSIGNMENT		(422)
#define FLATSPANNING_ASSIGNMENT		(431)
#define COILSPANNING_ASSIGNMENT		(432)
#define MINDIHEDRAL_ASSIGNMENT		(441)
#define MAXDIHEDRAL_ASSIGNMENT		(442)
#define MINBLENDING_ASSIGNMENT		(451)
#define MAXBLENDING_ASSIGNMENT		(452)
#define MINCURVATURE_ASSIGNMENT		(461)
#define MAXCURVATURE_ASSIGNMENT		(462)
#define MINCONCAVITY_ASSIGNMENT		(471)
#define MAXCONCAVITY_ASSIGNMENT		(472)
#define HYPERBOLICAVE_ASSIGNMENT	(481)
#define HYPERBOLICMAX_ASSIGNMENT	(482)
#define HYPERBOLICMIN_ASSIGNMENT	(483)
#define SKELFLAT_ASSIGNMENT		(491)
#define SKELCOIL_ASSIGNMENT		(492)

#define FLIP_EDGES			(510)

#define SPANNING_TREES			(520)
#define NORMAL_SMALL			(521)
#define NORMAL_MEDIUM			(522)
#define NORMAL_LARGE			(523)
#define BUMPY_SMALL			(531)
#define BUMPY_MEDIUM			(532)
#define BUMPY_LARGE			(533)
#define STRIP_SMALL			(541)
#define STRIP_MEDIUM			(542)
#define STRIP_LARGE			(543)
#define PIECE_SMALL			(551)
#define PIECE_MEDIUM			(552)
#define PIECE_LARGE			(553)
#define MATCHING_SMALL			(561)
#define MATCHING_MEDIUM			(562)
#define MATCHING_LARGE			(563)
#define SINGLE_FACE			(581)

#define	NUMBER_ONLY			(701)
#define	NUMBER_AREA			(702)
#define	NUMBER_RUN			(703)
#define	NUMBER_STAR			(704)
#define	NUMBER_WEIGHT			(705)
#define	NUMBER_AREA_RUN			(706)
#define	NUMBER_AREA_STAR		(707)
#define	NUMBER_AREA_BIAS		(708)
#define	NUMBER_AREA_LENGTH		(709)
#define	NUMBER_MAX			(710)
#define	NUMBER_AREA_MAX			(711)
#define	NUMBER_AVE			(712)
#define	NUMBER_AREA_AVE			(713)
#define	TEST_MEASURE			(714)

#define	GREEDY_STITCH			(751)
#define	EXHAUSTIVE_STITCH		(752)
#define GA_STITCH			(753)
#define ONESTEP_STITCH			(754)
#define ONESTEP_REMERGE			(755)
#define GREEDY_REMERGE			(756)

#define DEVIATION_050			(788)
#define DEVIATION_055			(789)
#define DEVIATION_060			(790)
#define DEVIATION_065			(791)
#define DEVIATION_070			(792)
#define DEVIATION_075			(793)
#define DEVIATION_080			(794)
#define DEVIATION_085			(795)

#define MINIMUM_SPANNINGTREE		(800)

#define LOAD_MESH			(810)
#define SAVE_MESH			(820)

#define CHANGE_VIEWPOINT		(850)

#define CAPTURE_MESH			(900)
#define CAPTURE_SHEET			(910)

#define QUIT				(999)

using namespace std;
using namespace boost;

//------------------------------------------------------------------------------
//	Variables
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//	Global variables
//------------------------------------------------------------------------------
// ID
// original triangular mesh
Polyhedron				mesh;
// dual graph of the original mesh
Graph					whole;
// temporary dual graph
Graph					partial;
// array of dual cycles
vector< vector< Facet_handle > >	pattern;
// array of 2D bounding boxes
vector< Bbox2 >				bound;
// paper size for mesh undolding	
Bbox2					paper( -DEFAULT_PAPER_WIDTH, -DEFAULT_PAPER_HEIGHT ,
					        DEFAULT_PAPER_WIDTH,  DEFAULT_PAPER_HEIGHT );
// weight assignment type for the dual edges
Weight_type				weight_type;

// Boundary cut path of the unfoled patches
Attribute				attr;

// curve skeleton
vector< Segment3 >			bone;
vector< Point3 >			refer;

//------------------------------------------------------------------------------
//	Tentative use
//------------------------------------------------------------------------------
// face barycenters [FID]
vector< Point3 >			bary;
// midpoints of the original edges [EID]
vector< Point3 >			mid;

//------------------------------------------------------------------------------
//	Mesh window
//------------------------------------------------------------------------------
// window IDs
int					win_mesh;

// window sizes
int					width_mesh	= DEFAULT_MESH_WIDTH;
int					height_mesh	= DEFAULT_MESH_HEIGHT;

// projection parameters
double					fovy		= 12.0;
double					aspect		= 1.0;
double					near		= 0.1;
double					far		= 100.0;

double					view_distance	= 10.0;
double					view_azimuth	= 0.0;	/*179.429688;*///45.0;   /*angle of object*/
//double					view_incidence	= 0.0; /*90.484375*/ //60.0;
//double				view_incidence	= 45.0; /*90.484375*/ //60.0;
double					view_incidence	= 90.0; /*90.484375*/ //60.0;
double					view_twist	= 0.0;

double					translate_x = 0.0, translate_y = 0.0, translate_z = 0.0;

// mouse positions and its button status
int					left_button = 0, middle_button = 0, right_button = 0;
int					pointer_x, pointer_y;

// drawing mode
//static int drawing_mode		= 1;
int					drawing_mode	= 4;
bool					saddle_flag	= false;
bool					label_flag	= false;
bool					number_flag	= false;
bool					reorder_flag	= false;
//static bool saddle_flag		= 0;

//------------------------------------------------------------------------------
//	Sheet window
//------------------------------------------------------------------------------
// window IDs
int				win_sheet;

// window sizes
int				width_sheet		= DEFAULT_SHEET_WIDTH;
int				height_sheet		= DEFAULT_SHEET_HEIGHT;

// strip splitting limit
double				limitCos		= INNER_PROD_DEGREE90;

double				mstDeviation		= DEFAULT_DEVIATION;


void redisplayAll	( void );
void displayMesh	( void );
void displaySheet	( void );
void displayBoth	( void );
void captureMesh	( char * name );
void captureSheet	( char * name );
void loadMesh		( Polyhedron & poly, char * filename );
void saveMesh		( Polyhedron & poly, char * filename );

//------------------------------------------------------------------------------
//	for debugging
//------------------------------------------------------------------------------
void printGraph( Graph & g )
{
    typedef property_map< Graph, vertex_index_t >::type IndexMap;
    IndexMap index = get( vertex_index, g );

//------------------------------------------------------------------------------
//	Accessing the Vertex Set
//------------------------------------------------------------------------------
    cerr << "vertices(g) = ";
    typedef graph_traits< Graph >::vertex_iterator vertex_iter;
    std::pair< vertex_iter, vertex_iter > vp;
    for ( vp = vertices( g ); vp.first != vp.second; ++vp.first )
	cerr << index[ *vp.first ] << " ";
    cerr << endl;
    
//------------------------------------------------------------------------------
//	Accessing the Edge Set
//------------------------------------------------------------------------------
    cerr << "edges(g) = ";
    graph_traits< Graph >::edge_iterator ei, ei_end;
    for ( tie( ei, ei_end ) = edges(g); ei != ei_end; ++ei ) 
	// std::cout << "(" << index[source(*ei. g)] << "," << index[target(*ei, g)] << ")";
	cerr << "(" << index[source(*ei,g)] << "," << index[target(*ei,g)] << ")";
    cerr << endl;
}


//------------------------------------------------------------------------------
//
//	For Mesh Window 
//
//------------------------------------------------------------------------------
void displayMesh( void )
{
    glutSetWindow( win_mesh );

    glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    gluPerspective( fovy, aspect, near, far );

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    glPushMatrix();
    glTranslated( 0, 0, -view_distance );
    glRotated( -view_twist, 0.0, 0.0, 1.0 );
    glRotated( -view_incidence, 1.0, 0.0, 0.0 );
    glRotated( -view_azimuth, 0.0, 0.0, 1.0 );
    glTranslated( translate_x, translate_y, translate_z );

    // cerr << " drawing_mode = " << drawing_mode << endl;
    switch ( drawing_mode ) {
      case 1:			// mesh only
	  hidden( GL_RENDER );
	  break;
      case 2:			// mesh with dual graph
	  duality( GL_RENDER );
	  hidden( GL_RENDER );
	  break;
      case 3:
	  // hidden( GL_RENDER );
	  layout( GL_RENDER );
	  break;
      case 4:	  
	  mst( GL_RENDER );
	  duality( GL_RENDER );
	  hidden( GL_RENDER );
	  break;
      case 5:	  
	  wireframe( GL_RENDER );
	  medial( GL_RENDER );
	  break;
      case 6:	  
	  cutedges( GL_RENDER );
	  break;
      case 7:
      case 8:
      case 9:
	  hidden( GL_RENDER );
	  break;
      default:
	  cerr << "Illegal drawing mode" << endl;
	  break;
    }

    if ( drawing_mode != 5 ) {
	if ( saddle_flag ) saddle( GL_RENDER );
	if ( label_flag ) annotate3D( GL_RENDER );
    }

    glPopMatrix();

    glFlush();
    glutSwapBuffers();
}


/* reshape the display size */
void reshapeMesh( int w, int h )
{
    width_mesh = w;
    height_mesh = h;
    glViewport ( 0, 0, width_mesh, height_mesh ); 
    glMatrixMode ( GL_PROJECTION );
    glLoadIdentity();
    gluPerspective( fovy, (double)width_mesh/(double)height_mesh, near, far );
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
}


// dump the window image in BMP format
void captureMesh( char * name )
{
    // int				i, j;		/* loop counters */
#ifdef AVOID_DOUBLE_FREE
    static unsigned char	 *buf_mesh = NULL;
    static IplImage*		image_mesh = NULL;	// Mesh image
#else	// AVOID_DOUBLE_FREE
    unsigned char		 *buf_mesh = NULL;
    IplImage*			image_mesh = NULL;	// Mesh image
#endif	// AVOID_DOUBLE_FREE
    int				h = ( int )height_mesh;
    int				w = ( int )width_mesh;

#ifdef OLD
    /* open the file */
    ofstream ofs( name, ios::binary );
    if ( ! ofs ) {
	cerr << "cannot open the file: " << name << endl;
	return;
    }
#endif	// OLD
    
    glutSetWindow( win_mesh );

    if ( buf_mesh == NULL ) buf_mesh = new unsigned char [ w * h * 3 ];
    glReadPixels( 0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, buf_mesh );
#ifdef OLD
#ifdef BINARY_PPM
    ofs << "P6" << endl;
#endif	// BINARY_PPM
#ifdef ASCII_PPM
    ofs << "P3" << endl;
#endif	// ASCII_PPM
    ofs << w << " " << h << endl;
    ofs << "255" << endl;
    for ( i = 0; i < h; i++ ) {
	for( j = 0; j < w; j++ ) {
#ifdef BINARY_PPM
	    ofs.write( (char*)&(buf_mesh[ ( (h-1-i) * w + j ) * 3 + 0 ]), sizeof( unsigned char ) );
	    ofs.write( (char*)&(buf_mesh[ ( (h-1-i) * w + j ) * 3 + 1 ]), sizeof( unsigned char ) );
	    ofs.write( (char*)&(buf_mesh[ ( (h-1-i) * w + j ) * 3 + 2 ]), sizeof( unsigned char ) );
#endif	// BINARY_PPM
#ifdef ASCII_PPM
	    ofs << setw(4) << (int)buf_mesh[ ( (h-1-i) * w + j ) * 3 + 0 ];
	    ofs << setw(4) << (int)buf_mesh[ ( (h-1-i) * w + j ) * 3 + 1 ];
	    ofs << setw(4) << (int)buf_mesh[ ( (h-1-i) * w + j ) * 3 + 2 ];
#endif	// ASCII_PPM
	}
#ifdef ASCII_PPM
	ofs << endl;
#endif	// ASCII_PPM
    }
    ofs.close();
#endif	// OLD

    if ( image_mesh != NULL ) {
	cvReleaseImage( &image_mesh );
	image_mesh = NULL;
    }
    image_mesh = cvCreateImage( cvSize( w, h ), IPL_DEPTH_8U, 3 );
    memcpy( image_mesh->imageData, buf_mesh, image_mesh->imageSize );

//    cvCvtColor( image_mesh, image_mesh, CV_BGR2RGB ); error
    cvCvtColor( image_mesh, image_mesh, CV_BGR2RGB );
    cvFlip( image_mesh, NULL, 0 );
    cvSaveImage( name, image_mesh );

    cerr << "Capturing the mesh window .. done" << endl;
#ifndef AVOID_DOUBLE_FREE
    delete [] buf_mesh;
    buf_mesh = NULL;
    cvReleaseImage( &image_mesh );
    image_mesh = NULL;
#endif	// AVOID_DOUBLE_FREE
}


void processHits( int nHits, unsigned int * buffer, int button )
{
    unsigned int * ptr = NULL; //, names;
    float minDepth = 1000.0;
    int hitID = NO_INDEX;

    cerr << "**** Here is processHits" << endl;

    ptr = buffer;

    for ( int i = 0; i < nHits; ++i ) { // for each bit
	if ( ptr[ 0 ] != 1 ) {
	    cerr << " Number of names for hit = " << ( int )ptr[ 0 ] << endl;
	    assert( ptr[ 0 ] == 1 );
	}
	float curDepth = (float)ptr[ 1 ]/0xffffffff;
	int curID = ( int )ptr[ 3 ];
	// #ifdef DEBUG
	cerr << " i = " << i 
	     << " [0]: " << ptr[ 0 ]
	     << " [1]: " << ptr[ 1 ]
	     << " [2]: " << ptr[ 2 ]
	     << " [3]: " << ptr[ 3 ] << endl;
	//#endif	// DEBUG
	if ( ( curDepth < minDepth ) && ( curID != NO_INDEX ) ) {
	    minDepth = curDepth;
	    hitID = ptr[ 3 ];
	}
	ptr += 4;
    }

    cerr << " hitID = " << hitID << " depth = " << minDepth << endl;

    if ( hitID != NO_INDEX ) {
	// initialize the halfedge IDs
	for ( Halfedge_iterator hi = mesh.halfedges_begin(); hi != mesh.halfedges_end(); ++hi ) {
	    if ( hi->id() == hitID ) {
		bool success;
		switch ( button ) {
		  case GLUT_MIDDLE_BUTTON:
		      success = fixWeight( partial, hi, 1.0 );
		      break;
		  case GLUT_RIGHT_BUTTON:
		      success = fixWeight( partial, hi, 0.0 );
		      break;
		}
		if ( success ) {
		    cerr << "Successfully assigned the specified weight value to the edge!!" << endl;
		    // cerr << "########### fitness function value = " << fitness( mesh ) << endl;
		}
		else {
		    cerr << "Could not handle the edge!!" << endl;
		}
	    }
	}
    }
}

void pickEdge( int x, int y, int button )
{
    unsigned int selectBuf[ BUFFER_SIZE ];
    int nHits;
    int viewport[ 4 ];

    glGetIntegerv( GL_VIEWPORT, viewport );

    //  Picking -- begin
    glSelectBuffer( BUFFER_SIZE, selectBuf );
    glRenderMode( GL_SELECT );

    glInitNames();
    // glLoadName( -1 );
    
    glMatrixMode( GL_PROJECTION );
    glPushMatrix(); // <====
    glLoadIdentity();
    // create small picking region near cursor location
    gluPickMatrix( (double)x, (double)(viewport[3]-y), 5.0, 5.0, viewport );
    // gluPickMatrix( (double)x, (double)(viewport[3]-y), 1.0, 1.0, viewport );
    gluPerspective( fovy, aspect, near, far );

    glMatrixMode( GL_MODELVIEW );
    glPushMatrix();
    glLoadIdentity();

    glPushMatrix();
    glTranslated( 0, 0, -view_distance );
    glRotated( -view_twist, 0.0, 0.0, 1.0 );
    glRotated( -view_incidence, 1.0, 0.0, 0.0 );
    glRotated( -view_azimuth, 0.0, 0.0, 1.0 );

    // duality( GL_SELECT );
    hidden( GL_SELECT );
    
    glPopMatrix();

    glMatrixMode( GL_PROJECTION );
    glPopMatrix();
    glMatrixMode( GL_MODELVIEW );
    glPopMatrix();

    glFlush();

    nHits = glRenderMode( GL_RENDER );
    processHits( nHits, selectBuf, button );

    redisplayAll();
}


/* handle the mouse motion */
void motionMesh( int x, int y )
{
    const double scale = 0.1;
    if ( left_button ) {
	/* The left button is used for the menu selection. */
	//  glReadPixels(x,width-y,1,1,GL_DEPTH_COMPONENT,GL_FLOAT,&z);
    }
    else if ( middle_button ) {
	view_distance -= scale * ( double )( y - pointer_y );
	pointer_x = x;
	pointer_y = y;
	// glutPostRedisplay();
	redisplayAll();
    }
    else if ( right_button ) {
	const int ratio = 2;
	view_azimuth -=
	    ( ( double )( x - pointer_x ) / ( double )( ratio * width_mesh ) ) * 180.0;
	view_incidence -=
	    ( ( double )( y - pointer_y ) / ( double )( ratio * height_mesh ) ) * 90.0;
	pointer_x = x;
	pointer_y = y;
	// glutPostRedisplay();
	redisplayAll();
    }
}

/* handle the mouse clicking */
void mouseMesh( int press, int state, int x, int y )
{
    if ( press == GLUT_LEFT_BUTTON ) {
	if ( state == GLUT_DOWN ) {
	    pointer_x = x;
	    pointer_y = y;
	    left_button = 1;
	}
	else {
	    left_button = 0;  
	}
    }
    else if ( press == GLUT_MIDDLE_BUTTON ) {
	if ( state == GLUT_DOWN ) {
	    int m = glutGetModifiers();
	    if ( m == GLUT_ACTIVE_SHIFT ) pickEdge( x, y, GLUT_MIDDLE_BUTTON );
	    else {
		pointer_x = x;
		pointer_y = y;
		middle_button = 1;
	    }
	}
	else {
	    middle_button = 0;
	}
    }
    else if ( press == GLUT_RIGHT_BUTTON ) {
	if ( state == GLUT_DOWN ) {
	    int m = glutGetModifiers();
	    if ( m == GLUT_ACTIVE_SHIFT ) pickEdge( x, y, GLUT_RIGHT_BUTTON );
	    else {
		pointer_x = x;
		pointer_y = y;
		right_button = 1;
	    }
	}
	else {
	    right_button = 0;
	}
    }
    // glutPostRedisplay();
    redisplayAll();
}


void keyboardMesh( unsigned char key, int x, int y )
{
    char buf[ 256 ], *ptr;
    char filename[ 256 ];
    // istringstream istr;
    static vector< Halfedge_handle > joint;
    vector< vector< Halfedge_handle > > skeleton;

    switch ( key ) {
      case 'R':
	  reorder_flag = !reorder_flag;
	  break;
      case 'N':
	  number_flag = !number_flag;
	  break;
      case 'v': // greedy search for unfolding
	  mergeByOne		( mesh, paper, pattern, bound );
	  arrangeMultiple	( pattern, paper, bound );
	  break;
      case 't': // exhaustive search for unfolding
	  greedy( mesh, paper, pattern, bound );
	  arrangeMultiple	( pattern, paper, bound );
	  labelBoundary		( mesh, attr );
	  break;
      case 'o': // optimized search for unfolding
	  optimize		( mesh, paper, pattern, bound, attr, displayBoth );
	  reportTime();
	  labelBoundary		( mesh, attr );
	  break;
      case 'e': // exhaustive search for unfolding
	  exhaustive		( mesh, paper, pattern, bound, attr, displayBoth );
	  labelBoundary		( mesh, attr );
	  break;
      case '%': // aligh the mesh
	  smoothMesh		( mesh );
	  renewAttrs		( mesh );
	  computeDual		( mesh, whole, bary, mid );
	  labelByCurvature	( mesh );
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  break;
      case 'd': // aligh the mesh
	  alignMesh		( mesh );
	  renewAttrs		( mesh );
	  computeDual		( mesh, whole, bary, mid );
	  labelByCurvature	( mesh );
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  break;
      case 'X': // rotate along the x-axis
	  rotateMesh		( mesh, 'x', 0.5*M_PI );
	  renewAttrs		( mesh );
	  computeDual		( mesh, whole, bary, mid );
	  labelByCurvature	( mesh );
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  break;
      case 'Y': // rotate along the y-axis
	  rotateMesh		( mesh, 'y', 0.5*M_PI );
	  renewAttrs		( mesh );
	  computeDual		( mesh, whole, bary, mid );
	  labelByCurvature	( mesh );
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  break;
      case 'Z': // rotate along the z-axis
	  rotateMesh		( mesh, 'z', 0.5*M_PI );
	  renewAttrs		( mesh );
	  computeDual		( mesh, whole, bary, mid );
	  labelByCurvature	( mesh );
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  break;
      case 'w':
	  fprintf		( stderr, " Displacement ( %6.3f, %6.3f, %6.3f ) : ",
				  translate_x, translate_y, translate_z );
	  ptr = fgets		( buf, sizeof( buf ), stdin );
	  sscanf		( buf, "%lf %lf %lf", 
				  &translate_x, &translate_y, &translate_z );
	  break;
      case 'S':
	  fprintf( stderr, " Input file name : " );
	  ptr = fgets( buf, sizeof( buf ), stdin );
	  sscanf( buf, "%s", filename );
	  saveMesh		( mesh, filename );
	  break;
      case 'L':
	  fprintf( stderr, " Input file name : " );
	  ptr = fgets( buf, sizeof( buf ), stdin );
	  sscanf( buf, "%s", filename );
	  loadMesh		( mesh, filename );
	  // normalize the mesh size
	  normalizeMesh		( mesh, bone );
	  // initialize the dual graph from the input 3D mesh
	  initAttrs		( mesh );
	  // calcluate total sum of the face aeras
	  attr.totalArea()	= totalArea( mesh );
	  // construct the connectivity of the dual graph
	  computeDual		( mesh, whole, bary, mid );
	  labelByCurvature	( mesh );
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  break;
      case 'y': // reconnect
	  reconnectByOne( mesh, paper, pattern, bound );
	  break;
      case 'g': // reconnect
	  reconnect( mesh, paper, pattern, bound, attr, displayBoth );
	  reportTime();
	  cerr << " occupation ratio = " << attr.totalArea()/( (paper.xmax()-paper.xmin())*
							       (paper.ymax()-paper.ymin()) ) << endl;
	  break;
      case '1':
	  drawing_mode = 1;
	  break;
      case '2':
	  drawing_mode = 2;
	  break;
      case '3':
	  drawing_mode = 3;
	  break;
      case '4':
	  drawing_mode = 4;
	  break;
      case '5':
	  drawing_mode = 5;
	  break;
      case '6':
	  drawing_mode = 6;
	  break;
      case '7':
	  drawing_mode = 7;
	  break;
      case '8':
	  drawing_mode = 8;
	  break;
      case '9':
	  drawing_mode = 9;
	  break;
      case 'a':
	  renewAttrs( mesh );
	  partial = whole;
	  perfectMatching( partial, mesh, skeleton );
	  projectCycles( skeleton, pattern, paper, bound, limitCos );
	  break;
      case 'm':			// MST
	  renewAttrs		( mesh );
	  computeDual		( mesh, whole, bary, mid ); 
	  labelByCurvature	( mesh );
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  minimumSpanningTree	( partial, mesh, pattern, paper, bound );
	  arrangeMultiple	( pattern, paper, bound );
	  countNCuts		( mesh );
	  labelBoundary		( mesh, attr );
	  cerr << "MST: Number of components = " << pattern.size() << endl;
	  cerr << " occupation ratio = " << attr.totalArea()/( (paper.xmax()-paper.xmin())*
							       (paper.ymax()-paper.ymin()) ) << endl;
	  drawing_mode = 3;
	  break;
      case 'x':			// MST Strip (Large)
	  renewAttrs( mesh );
	  computeDual( mesh, whole, bary, mid );
	  labelByCurvature( mesh );
	  partial = whole;
	  assignWeights( mesh, partial, weight_type );
	  MSTStrips( partial, mesh, skeleton, pattern, INNER_PROD_ALL, paper, bound );
	  countNCuts( mesh );
	  drawing_mode = 3;
	  break;
      case 'z':			// MST Piece (Large)
	  renewAttrs( mesh );
	  computeDual( mesh, whole, bary, mid );
	  labelByCurvature( mesh );
	  partial = whole;
	  assignWeights( mesh, partial, weight_type );
	  MSTPieces( partial, mesh, skeleton, pattern, INNER_PROD_ALL, paper, bound );
	  countNCuts( mesh );
	  drawing_mode = 3;
	  break;
      case 's':			// Decomposition into spanning trees by cutting
				// ridge edges
	  renewAttrs		( mesh );
	  computeDual		( mesh, whole, bary, mid );
	  labelByCurvature	( mesh );
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  spanningTrees		( partial, mesh, pattern, paper, bound );
	  countNCuts		( mesh );
	  labelBoundary		( mesh, attr );
	  drawing_mode = 3;
	  break;
      case 'r':
	  // renewAttrs		( mesh );
	  // computeDual		( mesh, whole, bary, mid );
	  // labelByCurvature	( mesh );
	  // partial = whole;
	  // assignWeights		( mesh, partial, weight_type );
	  spanningTrees		( partial, mesh, pattern, paper, bound );
	  countNCuts		( mesh );
	  labelBoundary		( mesh, attr );
	  drawing_mode = 3;
	  break;
      case 'c':
	  fprintf( stderr, " Input file name : " );
	  ptr = fgets( buf, sizeof( buf ), stdin );
	  sscanf( buf, "%s", filename );
	  captureMesh( filename );
	  break;
      case 'u':
	  fprintf( stderr, " Input file name : " );
	  ptr = fgets( buf, sizeof( buf ), stdin );
	  sscanf( buf, "%s", filename );
	  captureSheet( filename );
	  break;
      case 'i':
	  checkInETime();
	  checkInCPUTime();
	  break;
      case 'q':
      case 27:
	  exit( 0 );
	  break;
    }
    // glutPostRedisplay();
    redisplayAll();
}


//------------------------------------------------------------------------------
//
//	For Sheet Window 
//
//------------------------------------------------------------------------------
void displaySheet( void )
{
    glutSetWindow( win_sheet );

    glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    glRotatef( 90.0, 0.0, 0.0, 1.0 );
    gluOrtho2D( paper.xmin(), paper.xmax(), paper.ymin(), paper.ymax() );

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    glPushMatrix();

    switch ( drawing_mode ) {
      case 1:
      case 2:
      case 3:
      case 4:
      case 5:
	  strip( GL_RENDER );
	  break;
      case 6:
	  craft( GL_RENDER );
	  break;
      case 7:
	  outline( GL_RENDER );
	  break;
      case 8:
	  folding( GL_RENDER );
	  break;
      case 9:
	  outline( GL_RENDER );
	  folding( GL_RENDER );
	  break;
      default:
	  break;
    }
    
    glPopMatrix();

    glutSwapBuffers();
    glFlush();
}


/* reshape the display size */
void reshapeSheet( int w, int h )
{
    width_sheet = w;
    height_sheet = h;
    glViewport ( 0, 0, width_sheet, height_sheet ); 
    glMatrixMode ( GL_PROJECTION );
    glLoadIdentity();
    gluOrtho2D( paper.xmin(), paper.xmax(), paper.ymin(), paper.ymax() );

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
}


/* dump the window image */
void captureSheet( char * name )
{
    // int				i, j;		/* loop counters */
#ifdef AVOID_DOUBLE_FREE
    static unsigned char	*buf_sheet = NULL;
    static IplImage*		image_sheet = NULL;	// Sheet image
#else	// AVOID_DOUBLE_FREE
    unsigned char	*buf_sheet = NULL;
    IplImage*		image_sheet = NULL;	// Sheet image
#endif	// AVOID_DOUBLE_FREE
    int				h = ( int )height_sheet;
    int				w = ( int )width_sheet;

#ifdef OLD
    /* open the file */
    ofstream ofs( name, ios::binary );
    if ( ! ofs ) {
	cerr << "cannot open the file: " << name << endl;
	return;
    }
#endif	// OLD
    
    glutSetWindow( win_sheet );

    if ( buf_sheet == NULL ) buf_sheet = new unsigned char [ w * h * 3 ];
    glReadPixels( 0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, buf_sheet );
#ifdef OLD
#ifdef BINARY_PPM
    ofs << "P6" << endl;
#endif	// BINARY_PPM
#ifdef ASCII_PPM
    ofs << "P3" << endl;
#endif	// ASCII_PPM
    ofs << w << "\t" << h << endl;
    ofs << "255" << endl;
    for ( i = 0; i < h; i++ ) {
	for( j = 0; j < w; j++ ) {
#ifdef BINARY_PPM
	    ofs.write( (char*)&(buf_sheet[ ( (h-1-i) * w + j ) * 3 + 0 ]), sizeof( unsigned char ) );
	    ofs.write( (char*)&(buf_sheet[ ( (h-1-i) * w + j ) * 3 + 1 ]), sizeof( unsigned char ) );
	    ofs.write( (char*)&(buf_sheet[ ( (h-1-i) * w + j ) * 3 + 2 ]), sizeof( unsigned char ) );
#endif	// BINARY_PPM
#ifdef ASCII_PPM
	    ofs << setw(4) << (int)buf_sheet[ ( (h-1-i) * w + j ) * 3 + 0 ];
	    ofs << setw(4) << (int)buf_sheet[ ( (h-1-i) * w + j ) * 3 + 1 ];
	    ofs << setw(4) << (int)buf_sheet[ ( (h-1-i) * w + j ) * 3 + 2 ];
#endif	// ASCII_PPM
	}
#ifdef ASCII_PPM
	ofs << endl;
#endif	// ASCII_PPM
    }
    ofs.close();
#endif	// OLD


    if ( image_sheet != NULL ) {
	cvReleaseImage( &image_sheet );
	image_sheet = NULL;
    }
    image_sheet = cvCreateImage( cvSize( w, h ), IPL_DEPTH_8U, 3 );
    memcpy( image_sheet->imageData, buf_sheet, image_sheet->imageSize );

    cvCvtColor( image_sheet, image_sheet, CV_BGR2RGB );
    cvFlip( image_sheet, NULL, 0 );
    cvSaveImage( name, image_sheet );

    cerr << "Capturing the sheet window .. done" << endl;
#ifndef AVOID_DOUBLE_FREE
    delete [] buf_sheet;
    buf_sheet = NULL;
    cvReleaseImage( &image_sheet );
    image_sheet = NULL;
#endif	// AVOID_DOUBLE_FREE
}    


/* handle the mouse motion */
void motionSheet( int x, int y )
{
#ifdef DEBUG
    const double scale = 0.1;
    if ( left_button ) {
	/* The left button is used for the menu selection. */
	//  glReadPixels(x,width-y,1,1,GL_DEPTH_COMPONENT,GL_FLOAT,&z);
    }
    else if ( middle_button ) {
	view_distance -= scale * ( double )( y - pointer_y );
	pointer_x = x;
	pointer_y = y;
    }
    else if ( right_button ) {
	const int ratio = 2;
	view_azimuth -=
	    ( ( double )( x - pointer_x ) / ( double )( ratio * width ) ) * 180.0;
	view_incidence -=
	    ( ( double )( y - pointer_y ) / ( double )( ratio * height ) ) * 90.0;
	pointer_x = x;
	pointer_y = y;
    }
#endif	// DEBUG
    // glutPostRedisplay();
    redisplayAll();
}

/* handle the mouse clicking */
void mouseSheet( int press, int state, int x, int y )
{
#ifdef DEBUG
    if ( press == GLUT_LEFT_BUTTON ) {
	if ( state == GLUT_DOWN ) {
	    pointer_x = x;
	    pointer_y = y;
	    left_button = 1;
	}
	else {
	    left_button = 0;  
	}
    }
    else if ( press == GLUT_MIDDLE_BUTTON ) {
	if ( state == GLUT_DOWN ) {
	    pointer_x = x;
	    pointer_y = y;
	    middle_button = 1;
	}
	else {
	    middle_button = 0;
	}
    }
    else if ( press == GLUT_RIGHT_BUTTON ) {
	if ( state == GLUT_DOWN ) {
	    pointer_x = x;
	    pointer_y = y;
	    right_button = 1;
	}
	else {
	    right_button = 0;
	}
    }
#endif	// DEBUG
    // glutPostRedisplay();
    redisplayAll();
}


void keyboardSheet( unsigned char key, int x, int y )
{
    redisplayAll();
}


//------------------------------------------------------------------------------
//	Menu initialization
//------------------------------------------------------------------------------
/* handling the menu selection */
void selectItem( int mode )
{
    char buf[ 256 ], *ptr;
    char filename[ 256 ];
    // istringstream istr;
    vector< Halfedge_handle > joint;
    vector< vector< Halfedge_handle > >	skeleton;

    /* int id,i; */
    switch( mode ) {
#ifdef DEBUG
      case LOAD_VOLUME:
	  fprintf( stderr, " Input file name : " );
	  fgets( buf, sizeof( buf ), stdin );
	  sscanf( buf, "%s", filename );
	  load( filename );
	  dirty = TRUE;
	  break;
      case REDRAW_CASE:
	  dirty = TRUE;
	  break;
#endif	// DEBUG
	  // for changing drawing modes
      case MESH_ONLY:
	  drawing_mode = 1;
	  break;
      case WITH_DUAL:
	  drawing_mode = 2;
	  break;
      case PIECE_LAYOUT:
	  drawing_mode = 3;
	  break;
      case WITH_MST:
	  drawing_mode = 4;
	  break;
      case MEDIAL_AXIS:
	  drawing_mode = 5;
	  break;
      case CRAFT_RENDERING:
	  drawing_mode = 6;
	  break;
      case PATTERN_RENDERING:
	  drawing_mode = 9;
	  break;
//------------------------------------------------------------------------------
//	for changing the flag for saddles
//------------------------------------------------------------------------------
      case SADDLE_ON:
	  saddle_flag = 1;
	  break;
      case SADDLE_OFF:
	  saddle_flag = 0;
	  break;
//------------------------------------------------------------------------------
//	for displaying the labels
//------------------------------------------------------------------------------
      case LABEL_ON:
	  label_flag = 1;
	  break;
      case LABEL_OFF:
	  label_flag = 0;
	  break;
//------------------------------------------------------------------------------
//	for transforming the mesh geometry
//------------------------------------------------------------------------------
      case ROTATE_XAXIS: // rotate along the x-axis
	  rotateMesh		( mesh, 'x', 0.5*M_PI );
	  renewAttrs		( mesh );
	  computeDual		( mesh, whole, bary, mid );
	  labelByCurvature	( mesh );
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  break;
      case ROTATE_YAXIS: // rotate along the y-axis
	  rotateMesh		( mesh, 'y', 0.5*M_PI );
	  renewAttrs		( mesh );
	  computeDual		( mesh, whole, bary, mid );
	  labelByCurvature	( mesh );
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  break;
      case ROTATE_ZAXIS: // rotate along the z-axis
	  rotateMesh		( mesh, 'z', 0.5*M_PI );
	  renewAttrs		( mesh );
	  computeDual		( mesh, whole, bary, mid );
	  labelByCurvature	( mesh );
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  break;
      case ALIGN_MESH: // aligh the mesh
	  alignMesh		( mesh );
	  renewAttrs		( mesh );
	  computeDual		( mesh, whole, bary, mid );
	  labelByCurvature	( mesh );
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  break;
      case SMOOTH_MESH: // smooth the mesh
	  smoothMesh		( mesh );
	  renewAttrs		( mesh );
	  computeDual		( mesh, whole, bary, mid );
	  labelByCurvature	( mesh );
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  break;
//------------------------------------------------------------------------------
//	for assigning specif weight values to dual edges
//------------------------------------------------------------------------------
      case UNIFORM_ASSIGNMENT:
	  weight_type	= UNIFORM_WEIGHT;
	  initWeights		( mesh );
	  cerr << "Uniform weight values assigned to dual edges." << endl;
	  break;
      case RANDOM_ASSIGNMENT:
	  weight_type	= RANDOM_WEIGHT;
	  initWeights		( mesh );
	  cerr << "Random weight values assigned to dual edges." << endl;
	  break;
      case MINPERIMETER_ASSIGNMENT:
	  weight_type	= MINIMUM_PERIMETER;
	  initWeights		( mesh );
	  cerr << "Minimum perimeter weight values assigned to dual edges." << endl;
	  break;
      case MAXPERIMETER_ASSIGNMENT:
	  weight_type	= MAXIMUM_PERIMETER;
	  initWeights		( mesh );
	  cerr << "Maximum perimeter weight values assigned to dual edges." << endl;
	  break;
      case FLATSPANNING_ASSIGNMENT:
	  weight_type	= FLAT_SPANNING;
	  initWeights		( mesh );
	  cerr << "Flat spanning weight values assigned to dual edges." << endl;
	  break;
      case COILSPANNING_ASSIGNMENT:
	  weight_type	= COIL_SPANNING;
	  initWeights		( mesh );
	  cerr << "Coil spanning weight values assigned to dual edges." << endl;
	  break;
      case MINDIHEDRAL_ASSIGNMENT:
	  weight_type	= MINIMUM_DIHEDRAL;
	  initWeights		( mesh );
	  cerr << "Mininum dihedral angle weight values assigned to dual edges." << endl;
	  cerr << "(Flat region is more likely to be split.)" << endl;
	  break;
      case MAXDIHEDRAL_ASSIGNMENT:
	  weight_type	= MAXIMUM_DIHEDRAL;
	  initWeights		( mesh );
	  cerr << "Maximum dihedral angle weight values assigned to dual edges." << endl;
	  cerr << "(Sharp region is more likely to be split.)" << endl;
	  break;
      case MINBLENDING_ASSIGNMENT:
	  weight_type	= MINIMUM_BLENDING;
	  initWeights		( mesh );
	  cerr << "Minimum perimeter and flat spanning weights blended and assigned to dual edges." << endl;
	  break;
      case MAXBLENDING_ASSIGNMENT:
	  weight_type	= MAXIMUM_BLENDING;
	  initWeights		( mesh );
	  cerr << "Maximum perimeter and flat spanning weights bleended and assigned to dual edges." << endl;
	  break;
      case MINCURVATURE_ASSIGNMENT:
	  weight_type	= MINIMUM_CURVATURE;
	  initWeights		( mesh );
	  cerr << "Minimum curvature weight values assigned to dual edges." << endl;
	  break;
      case MAXCURVATURE_ASSIGNMENT:
	  weight_type	= MAXIMUM_CURVATURE;
	  initWeights		( mesh );
	  cerr << "Maximum curvature weight values assigned to dual edges." << endl;
	  break;
      case MINCONCAVITY_ASSIGNMENT:
	  weight_type	= MINIMUM_CONCAVITY;
	  initWeights		( mesh );
	  cerr << "Minimum hollowness weight values assigned to dual edges." << endl;
	  break;
      case MAXCONCAVITY_ASSIGNMENT:
	  weight_type	= MAXIMUM_CONCAVITY;
	  initWeights		( mesh );
	  cerr << "Maximum hollowness weight values assigned to dual edges." << endl;
	  break;
      case HYPERBOLICAVE_ASSIGNMENT:
	  weight_type	= MINIMUM_HYPERBOLIC_AVE;
	  initWeights		( mesh );
	  cerr << "Mininum average hyperbolic weight values assigned to dual edges." << endl;
	  break;
      case HYPERBOLICMAX_ASSIGNMENT:
	  weight_type	= MINIMUM_HYPERBOLIC_MAX;
	  initWeights		( mesh );
	  cerr << "Minimum maximum hyperbolic weight values assigned to dual edges." << endl;
	  break;
      case HYPERBOLICMIN_ASSIGNMENT:
	  weight_type	= MINIMUM_HYPERBOLIC_MIN;
	  initWeights		( mesh );
	  cerr << "Mininum minimum hyperbolic weight values assigned to dual edges." << endl;
	  break;
      case SKELFLAT_ASSIGNMENT:
	  weight_type	= SKELETON_FLAT;
	  initWeights		( mesh );
	  cerr << "Skeleton flat assigned to dual edges." << endl;
	  break;
      case SKELCOIL_ASSIGNMENT:
	  weight_type	= SKELETON_FLAT;
	  initWeights		( mesh );
	  cerr << "Skeleton coil assigned to dual edges." << endl;
	  break;
//------------------------------------------------------------------------------
//	
//------------------------------------------------------------------------------
      case SINGLE_FACE:
	  renewAttrs		( mesh );
	  computeDual		( mesh, whole, bary, mid );
	  labelByCurvature	( mesh );
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  singleFaces		( mesh, skeleton );
	  projectCycles		( skeleton, pattern, paper, bound, INNER_PROD_DEGREE60 );
	  arrangeMultiple	( pattern, paper, bound );
	  break;
      case MATCHING_SMALL: 
	  renewAttrs		( mesh );
	  computeDual		( mesh, whole, bary, mid );
	  labelByCurvature	( mesh );
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  weightedMatching	( partial, mesh, skeleton );
	  projectCycles		( skeleton, pattern, paper, bound, INNER_PROD_DEGREE60 );
	  resolveOverlaps	( mesh, pattern );
	  arrangeMultiple	( pattern, paper, bound );
	  countNCuts		( mesh );
	  break;
      case MATCHING_MEDIUM: 
	  renewAttrs		( mesh ); 
	  computeDual		( mesh, whole, bary, mid ); 
	  labelByCurvature	( mesh ); 
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  weightedMatching	( partial, mesh, skeleton );
	  projectCycles		( skeleton, pattern, paper, bound, INNER_PROD_DEGREE90 );
	  resolveOverlaps	( mesh, pattern );
	  arrangeMultiple	( pattern, paper, bound );
	  countNCuts		( mesh );
	  break;
      case MATCHING_LARGE:
	  renewAttrs		( mesh );
	  computeDual		( mesh, whole, bary, mid );
	  labelByCurvature	( mesh );
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  weightedMatching	( partial, mesh, skeleton );
	  projectCycles		( skeleton, pattern, paper, bound, INNER_PROD_ALL );
	  resolveOverlaps	( mesh, pattern );
	  arrangeMultiple	( pattern, paper, bound );
	  countNCuts		( mesh );
	  break;
      case STRIP_SMALL:		// MST-based layout (large patches) 
	  renewAttrs		( mesh );
	  computeDual		( mesh, whole, bary, mid );
	  labelByCurvature	( mesh );
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  MSTStrips		( partial, mesh, skeleton, pattern, INNER_PROD_DEGREE60, paper, bound );
	  arrangeMultiple	( pattern, paper, bound );
	  countNCuts		( mesh );
	  drawing_mode = 3;
	  break;
      case STRIP_MEDIUM:		// MST-based layout (large patches) 
	  renewAttrs		( mesh );
	  computeDual		( mesh, whole, bary, mid );
	  labelByCurvature	( mesh );
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  MSTStrips		( partial, mesh, skeleton, pattern, INNER_PROD_DEGREE90, paper, bound );
	  arrangeMultiple	( pattern, paper, bound );
	  countNCuts		( mesh );
	  drawing_mode = 3;
	  break;
      case STRIP_LARGE:		// MST-based layout (large patches) 
	  renewAttrs		( mesh );
	  computeDual		( mesh, whole, bary, mid );
	  labelByCurvature	( mesh );
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  MSTStrips		( partial, mesh, skeleton, pattern, INNER_PROD_ALL, paper, bound );
	  countNCuts		( mesh );
	  drawing_mode = 3;
	  break;
      case PIECE_SMALL:		// MST-based layout (large patches) 
	  renewAttrs		( mesh );
	  computeDual		( mesh, whole, bary, mid );
	  labelByCurvature	( mesh );
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  MSTPieces		( partial, mesh, skeleton, pattern, INNER_PROD_DEGREE60, paper, bound );
	  arrangeMultiple	( pattern, paper, bound );
	  countNCuts		( mesh );
	  drawing_mode = 3;
	  break;
      case PIECE_MEDIUM:	// MST-based layout (large patches) 
	  renewAttrs		( mesh );
	  computeDual		( mesh, whole, bary, mid );
	  labelByCurvature	( mesh );
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  MSTPieces		( partial, mesh, skeleton, pattern, INNER_PROD_DEGREE90, paper, bound );
	  arrangeMultiple	( pattern, paper, bound );
	  countNCuts		( mesh );
	  drawing_mode = 3;
	  break;
      case PIECE_LARGE:		// MST-based layout (large patches) 
	  renewAttrs		( mesh );
	  computeDual		( mesh, whole, bary, mid );
	  labelByCurvature	( mesh );
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  MSTPieces		( partial, mesh, skeleton, pattern, INNER_PROD_ALL, paper, bound );
	  arrangeMultiple	( pattern, paper, bound );
	  countNCuts		( mesh );
	  drawing_mode = 3;
	  break;
      case NORMAL_SMALL: 
	  renewAttrs		( mesh );
	  computeDual		( mesh, whole, bary, mid );
	  labelByCurvature	( mesh );
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  perfectMatching	( partial, mesh, skeleton );
	  projectCycles		( skeleton, pattern, paper, bound, INNER_PROD_DEGREE60 );
	  arrangeMultiple	( pattern, paper, bound );
	  break;
      case NORMAL_MEDIUM: 
	  renewAttrs		( mesh ); 
	  computeDual		( mesh, whole, bary, mid ); 
	  labelByCurvature	( mesh ); 
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  perfectMatching	( partial, mesh, skeleton );
	  projectCycles		( skeleton, pattern, paper, bound, INNER_PROD_DEGREE90 );
	  arrangeMultiple	( pattern, paper, bound );
	  break;
      case NORMAL_LARGE:
	  renewAttrs		( mesh );
	  computeDual		( mesh, whole, bary, mid );
	  labelByCurvature	( mesh );
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  perfectMatching	( partial, mesh, skeleton );
	  projectCycles		( skeleton, pattern, paper, bound, INNER_PROD_ALL );
	  resolveOverlaps	( mesh, pattern );
	  arrangeMultiple	( pattern, paper, bound );
	  break;
	  // Constraints are considered in BUMPY_SMALL, BUMPY_MEDIUM, BUMPY_LARGE
      case BUMPY_SMALL:	// lay out the cut path over the saddle points
	  renewAttrs		( mesh );
	  computeDual		( mesh, whole, bary, mid );
	  labelByCurvature	( mesh );
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  applyConstraints	( mesh, partial );
	  perfectMatching	( partial, mesh, skeleton );
	  projectCycles		( skeleton, pattern, paper, bound, INNER_PROD_DEGREE60 );
	  arrangeMultiple	( pattern, paper, bound );
	  break;
      case BUMPY_MEDIUM: // Lay out the cut path over the saddle points
	  // Initialize the dual center
	  renewAttrs		( mesh );
	  // Compute and store the results of dual graphs into appropriate variables
	  computeDual		( mesh, whole, bary, mid );
	  // Validate the type of vertices: saddle vertices, peak vertices and default vertices
	  labelByCurvature	( mesh );
	  // For perfect matching purposes
	  partial = whole; 
	  // Weight assignment
	  assignWeights		( mesh, partial, weight_type );
	  // Constraint based on saddle points
	  applyConstraints	( mesh, partial ); 
	  // Perfect matching operations
	  perfectMatching	( partial, mesh, skeleton ); 
	  // The dual angle condition for splitting: 90 degree to avoid local self-intersection
	  // The 2D unfoldings according to the limitCos
	  projectCycles		( skeleton, pattern, paper, bound, INNER_PROD_DEGREE90 ); 
	  break;
      case BUMPY_LARGE:	// lay out the cut path over the saddle points
	  renewAttrs		( mesh );
	  computeDual		( mesh, whole, bary, mid );
	  labelByCurvature	( mesh );
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  applyConstraints	( mesh, partial );
	  perfectMatching	( partial, mesh, skeleton );
	  projectCycles		( skeleton, pattern, paper, bound, INNER_PROD_ALL );
	  resolveOverlaps	( mesh, pattern );
	  arrangeMultiple	( pattern, paper, bound );
	  countNCuts		( mesh );
	  break;
      case SPANNING_TREES:	// Decomposition into spanning trees by cutting
				// ridge edges
	  renewAttrs		( mesh );
	  computeDual		( mesh, whole, bary, mid );
	  labelByCurvature	( mesh );
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  spanningTrees		( partial, mesh, pattern, paper, bound );
	  arrangeMultiple	( pattern, paper, bound );
	  countNCuts		( mesh );
	  labelBoundary		( mesh, attr );
	  drawing_mode = 3;
	  break;
//------------------------------------------------------------------------------
//	Stitching operations
//------------------------------------------------------------------------------
      case GREEDY_STITCH:	// for debug
	  greedy		( mesh, paper, pattern, bound );
	  arrangeMultiple	( pattern, paper, bound );
	  break;
      case EXHAUSTIVE_STITCH:
	  divide		( mesh, paper, pattern, bound );
	  arrangeMultiple	( pattern, paper, bound );
	  break;
	  // Prautzsch's paper-minimum spanning tree
      case GA_STITCH:
	  optimize		( mesh, paper, pattern, bound, attr, displayBoth );
	  arrangeMultiple	( pattern, paper, bound );
	  reportTime();
	  labelBoundary		( mesh, attr );
	  drawing_mode = 3;
	  break;
      case ONESTEP_STITCH:
	  mergeByOne		( mesh, paper, pattern, bound );
	  arrangeMultiple	( pattern, paper, bound );
	  break;
      case ONESTEP_REMERGE:
	  reconnectByOne	( mesh, paper, pattern, bound );
	  arrangeMultiple	( pattern, paper, bound );
	  break;
      case GREEDY_REMERGE:
	  reconnect		( mesh, paper, pattern, bound, attr, displayBoth );
	  arrangeMultiple	( pattern, paper, bound );
	  reportTime();
	  cerr << " occupation ratio = " << attr.totalArea()/( (paper.xmax()-paper.xmin())*
							       (paper.ymax()-paper.ymin()) ) << endl;
	  cerr << "Time passed = " << checkOutETime() << endl;
	  cerr << "CPU Time    = " << checkOutCPUTime() << endl;
	  break;
//------------------------------------------------------------------------------
//	Objective function selection
//------------------------------------------------------------------------------
      case NUMBER_ONLY:
	  Genome::setFitness( number_only );
	  cerr << "Number only" << endl;
	  break;
      case NUMBER_AREA:
	  Genome::setFitness( number_area );
	  cerr << "Number + Area" << endl;
	  break;
      case NUMBER_RUN:
	  Genome::setFitness( number_run );
	  cerr << "Number + Run" << endl;
	  break;
      case NUMBER_STAR:
	  Genome::setFitness( number_star );
	  cerr << "Number + Star" << endl;
	  break;
      case NUMBER_WEIGHT:
	  Genome::setFitness( number_weight );
	  cerr << "Number + Weight" << endl;
	  break;
      case NUMBER_AREA_RUN:
	  Genome::setFitness( number_area_run );
	  cerr << "Number + Area + Run" << endl;
	  break;
      case NUMBER_AREA_STAR:
	  Genome::setFitness( number_area_star );
	  cerr << "Number + Area + Star" << endl;
	  break;
      case NUMBER_AREA_BIAS:
	  Genome::setFitness( number_area_bias );
	  cerr << "Number + Area + Biased Star" << endl;
	  break;
      case NUMBER_AREA_LENGTH:
	  Genome::setFitness( number_area_length );
	  cerr << "Number + Area + Length" << endl;
	  break;
//------------------------------------------------------------------------------
      case NUMBER_MAX:
	  Genome::setFitness( number_max );
	  cerr << "Number + DualMax" << endl;
	  break;
      case NUMBER_AREA_MAX:
	  Genome::setFitness( number_area_max );
	  cerr << "Number + Area + DualMax" << endl;
	  break;
      case NUMBER_AVE:
	  Genome::setFitness( number_ave );
	  cerr << "Number + DualAverage" << endl;
	  break;
      case NUMBER_AREA_AVE:
	  Genome::setFitness( number_area_ave );
	  cerr << "Number + Area + DualAverage" << endl;
	  break;
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
      case TEST_MEASURE:
	  Genome::setFitness( test_measure );
	  cerr << "Test Measure" << endl;
	  break;
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
      case DEVIATION_050:
	  mstDeviation = 0.00;
	  break;
      case DEVIATION_055:
	  mstDeviation = 0.12;
	  break;
      case DEVIATION_060:
	  mstDeviation = 0.23;
	  break;
      case DEVIATION_065:
	  mstDeviation = 0.38;
	  break;
      case DEVIATION_070:
	  mstDeviation = 0.52;
	  break;
      case DEVIATION_075:
	  mstDeviation = 0.67;
	  break;
      case DEVIATION_080:
	  mstDeviation = 0.84;
	  break;
      case DEVIATION_085:
	  mstDeviation = 1.03;
	  break;
//------------------------------------------------------------------------------
//	
//------------------------------------------------------------------------------
      case MINIMUM_SPANNINGTREE:
	  renewAttrs		( mesh );
	  computeDual		( mesh, whole, bary, mid ); 
	  labelByCurvature	( mesh );
	  partial = whole;
	  assignWeights		( mesh, partial, weight_type );
	  minimumSpanningTree	( partial, mesh, pattern, paper, bound );
	  arrangeMultiple	( pattern, paper, bound );
	  countNCuts		( mesh );
	  labelBoundary		( mesh, attr );
	  cerr << "MST: Number of components = " << pattern.size() << endl;
	  cerr << " occupation ratio = " << attr.totalArea()/( (paper.xmax()-paper.xmin())*
							       (paper.ymax()-paper.ymin()) ) << endl;
	  drawing_mode = 3;
	  break;
//------------------------------------------------------------------------------
//	
//------------------------------------------------------------------------------
      case CHANGE_VIEWPOINT:
	  fprintf( stderr,"Input view parameters: (distance = %6.3f):", view_distance );
	  ptr = fgets( buf, sizeof(buf), stdin );
	  sscanf( buf, "%lf", &view_distance );
	  fprintf( stderr,"Input view parameters: (incidence = %6.3f):", view_incidence );
	  ptr = fgets( buf, sizeof(buf), stdin );
	  sscanf( buf, "%lf", &view_incidence );
	  fprintf( stderr,"Input view parameters: (azimuth = %6.3f):", view_azimuth );
	  ptr = fgets( buf, sizeof(buf), stdin );
	  sscanf( buf, "%lf", &view_azimuth );
	  break;
//------------------------------------------------------------------------------
//	
//------------------------------------------------------------------------------
      case LOAD_MESH:
	  fprintf( stderr, " Input file name : " );
	  ptr = fgets( buf, sizeof( buf ), stdin );
	  sscanf( buf, "%s", filename );
	  loadMesh( mesh, filename );
	  // normalize the mesh size
	  normalizeMesh	( mesh, bone );
	  // initialize the dual graph from the input 3D mesh
	  initAttrs		( mesh );
	  // calculate the total area
	  attr.totalArea() = totalArea( mesh );
	  // construct the connectivity of the dual graph
	  computeDual		( mesh, whole, bary, mid );
	  labelByCurvature	( mesh );
	  partial = whole;
	  assignWeights	( mesh, partial, weight_type );
	  break;
      case SAVE_MESH:
	  fprintf( stderr, " Input file name : " );
	  ptr = fgets( buf, sizeof( buf ), stdin );
	  sscanf( buf, "%s", filename );
	  saveMesh( mesh, filename );
	  break;
//------------------------------------------------------------------------------
//	
//------------------------------------------------------------------------------
      case CAPTURE_MESH:
	  fprintf( stderr, "Capturing the Mesh Window... Input file name : " );
	  ptr = fgets( buf, sizeof( buf ), stdin );
	  sscanf( buf, "%s", filename );
	  captureMesh( filename );
	  return;
      case CAPTURE_SHEET:
	  fprintf( stderr, "Capturing the Sheet Window.. Input file name : " );
	  ptr = fgets( buf, sizeof( buf ), stdin );
	  sscanf( buf, "%s", filename );
	  captureSheet( filename );
	  return;
//------------------------------------------------------------------------------
//	
//------------------------------------------------------------------------------
      case QUIT:
	  exit(0);
      default:
	  break;
    }
    // glutPostRedisplay();
    redisplayAll();
}


//------------------------------------------------------------------------------
//	Common functions
//------------------------------------------------------------------------------
void redisplayAll( void )
{
    glutSetWindow( win_mesh );
    glutPostRedisplay();
    glutSetWindow( win_sheet );
    glutPostRedisplay();
}

void displayBoth( void )
{
    // glutSetWindow( win_mesh );
    displayMesh();
    // glutSetWindow( win_sheet );
    displaySheet();
}

//------------------------------------------------------------------------------
//	Special function
//------------------------------------------------------------------------------
void loadMesh( Polyhedron & poly, char * filename )
{
    ifstream ifs( filename );

    if ( !ifs) {
//    if ( ifs == NULL ) { error
	cerr << "Cannot open the file : " << filename << endl;
	exit( 0 );
    }

    mesh.clear();
    ifs >> mesh;
    mesh.normalize_border();
    ifs.close();

    cerr << "Loading done." << endl;
}


void saveMesh( Polyhedron & poly, char * filename )
{
    ofstream ofs( filename );

//    if ( ofs == NULL ) { error
    if ( !ofs  ) {
        cerr << "Cannot open the file : " << filename << endl;
	exit( 0 );
    }

    ofs << mesh;
    ofs.close();

    cerr << "Saving done." << endl;
}


void loadSkeleton( char * name )
{
    ifstream ifs( name );
    char buf[ 256 ];
    istringstream istr;
    vector< Point3 > plot;
    double eps = 1.0e-4;

//    if ( ifs == NULL ) { error
    if ( !ifs ) {
	cerr << "No skeleton file : " << name << endl;
	return;
    }

    ifs.getline( buf, sizeof(buf) );
    istr.clear();
    istr.str( buf );
    assert( ( buf[ 0 ] == 'O' ) && ( buf[ 1 ] == 'F' ) && ( buf[ 2 ] == 'F' ) );

    ifs.getline( buf, sizeof(buf) );
    istr.clear();
    istr.str( buf );
    int nV, nE, nZ;
    istr >> nV >> nE >> nZ;
    cerr << "#vertices = " << nV << "#edges = " << nE << endl;
    
    plot.clear();
    for ( int k = 0; k < nV; ++k ) {
	ifs.getline( buf, sizeof(buf) );
	istr.clear();
	istr.str( buf );
	Point3 coord;
	istr >> coord;
	plot.push_back( coord );
	// cerr << " k = " << k << " coord = " << coord << endl;
    }

    bone.clear();
    for ( int k = 0; k < nE; ++k ) {
	ifs.getline( buf, sizeof(buf) );
	istr.clear();
	istr.str( buf );
	int n, idA, idB;
	istr >> n >> idA >> idB;
	assert( n == 2 );
	// cerr << " k = " << k << " pair = " << idA << " -- " << idB << endl;
	Segment3 line = Segment3( plot[ idA ], plot[ idB ] );
	if ( line.squared_length() > eps ) {
	    bone.push_back( line );
	}
	else {
	    cerr << "Too small bone segment : " << endl;
	    cerr << "Source: " << line.source() << endl;
	    cerr << "Target: " << line.target() << endl;
	    cerr << endl;
	}
    }

    plot.clear();
    ifs.close();
}    


//------------------------------------------------------------------------------
//	Main function
//------------------------------------------------------------------------------
int main( int argc, char *argv[] ) 
{
    if ( ( argc <= 1 ) || ( argc >= 4 ) ) {
	cerr << "Usage : " << argv[ 0 ] << " <off file name> (<skeleton file name>)" << endl;
	exit( 0 );
    }

    loadMesh( mesh, argv[ 1 ] );
    if ( argc == 3 ) loadSkeleton( argv[ 2 ] );

	
//------------------------------------------------------------------------------
//	Several procedures
//------------------------------------------------------------------------------
    // normalize the mesh size
    normalizeMesh	( mesh, bone );
    // initialize the dual graph from the input 3D mesh
    initAttrs		( mesh );
    // calculate the total area
    attr.totalArea() = totalArea( mesh );
    // construct the connectivity of the dual graph
    computeDual		( mesh, whole, bary, mid );
    labelByCurvature	( mesh );
    partial = whole;
    assignWeights	( mesh, partial, weight_type );

//------------------------------------------------------------------------------
//	OpenGL setup
//------------------------------------------------------------------------------

    if ( mesh.size_of_border_edges() != 0 ) {
	cerr << "The input object has border edges. Cannot unfold." << endl;
	exit( 0 );
    }

    glutInit( &argc, argv );                /* Initialize GLUT */

//------------------------------------------------------------------------------
//	Genetic algorithm setup
//------------------------------------------------------------------------------
    // Genome::setFitness( number_area );
    // Genome::setFitness( number_area_length );
    // Genome::setFitness( number_area_run );
    Genome::setFitness( number_area_ave );

//------------------------------------------------------------------------------
//	Mesh window
//------------------------------------------------------------------------------
    glutInitDisplayMode ( GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );
    glutInitWindowSize ( width_mesh, height_mesh ); 
    glutInitWindowPosition( 0, 0 );
    win_mesh = glutCreateWindow ( "Mesh" );

    initMesh();
    // Default edge weight type
    weight_type	= MINIMUM_PERIMETER;
    // weight_type	= MINIMUM_CONCAVITY;

    glutDisplayFunc( displayMesh ); 
    glutReshapeFunc( reshapeMesh );
    glutKeyboardFunc( keyboardMesh );
    glutMouseFunc( mouseMesh );
    glutMotionFunc( motionMesh );

    int menu_render	= glutCreateMenu( selectItem );
    glutAddMenuEntry    ( "[1] Mesh Only",		MESH_ONLY );
    glutAddMenuEntry    ( "[2] With Dual",		WITH_DUAL );
    glutAddMenuEntry    ( "[3] Piece Layout",		PIECE_LAYOUT );
    glutAddMenuEntry	( "[4] With MST",		WITH_MST );
    // glutAddMenuEntry	( "[5] Curve skeleton",		MEDIAL_AXIS );
    glutAddMenuEntry	( "[6] Papercraft",		CRAFT_RENDERING );
    glutAddMenuEntry	( "[9] Pattern",		PATTERN_RENDERING );

#ifdef DISPLAY_OPTION
    int menu_saddle	= glutCreateMenu( selectItem );
    glutAddMenuEntry    ( "On",				SADDLE_ON );
    glutAddMenuEntry    ( "Off",			SADDLE_OFF );
#endif	// DISPLAY_OPTION

#ifdef DISPLAY_OPTION
    int menu_label	= glutCreateMenu( selectItem );
    glutAddMenuEntry    ( "On",				LABEL_ON );
    glutAddMenuEntry    ( "Off",			LABEL_OFF );
#endif	// DISPLAY_OPTION

#ifdef DISPLAY_OPTION
    int menu_transform	= glutCreateMenu( selectItem );
    glutAddMenuEntry    ( "Rotate X-axis",		ROTATE_XAXIS );
    glutAddMenuEntry    ( "Rotate Y-axis",		ROTATE_YAXIS );
    glutAddMenuEntry    ( "Rotate Z-axis",		ROTATE_ZAXIS );
    glutAddMenuEntry    ( "Align Mesh",			ALIGN_MESH );
    glutAddMenuEntry    ( "Smooth Mesh",		SMOOTH_MESH );
#endif	// DISPLAY_OPTION

    int menu_split	= glutCreateMenu( selectItem );
    glutAddMenuEntry	( "[s] Spanning Trees",		SPANNING_TREES );
#ifdef DISPLAY_OPTION
    glutAddMenuEntry	( "MST Strip (Small)",		STRIP_SMALL );
    glutAddMenuEntry	( "MST Strip (Medium)",		STRIP_MEDIUM );
#endif	// DISPLAY_OPTION
    glutAddMenuEntry	( "[x] MST Strip (Large)",	STRIP_LARGE );
#ifdef DISPLAY_OPTION
    glutAddMenuEntry	( "MST Piece (Small)",		PIECE_SMALL );
    glutAddMenuEntry	( "MST Piece (Medium)",		PIECE_MEDIUM );
#endif	// DISPLAY_OPTION
    glutAddMenuEntry	( "[z] MST Piece (Large)",	PIECE_LARGE );
#ifdef DISPLAY_OPTION
    glutAddMenuEntry	( "PM-Based (Small)",		MATCHING_SMALL );
    glutAddMenuEntry	( "PM-Based (Medium)",		MATCHING_MEDIUM );
    glutAddMenuEntry	( "PM-Based (Large)",		MATCHING_LARGE );
#endif	// DISPLAY_OPTION
    glutAddMenuEntry	( "Single Faces",		SINGLE_FACE );
#ifdef OLD
    glutAddMenuEntry	( "Normal (Small)",		NORMAL_SMALL );
    glutAddMenuEntry	( "Normal (Medium)",		NORMAL_MEDIUM );
    glutAddMenuEntry	( "Normal (Large)",		NORMAL_LARGE );
    glutAddMenuEntry	( "Bump-Based (Small)",		BUMPY_SMALL );
    glutAddMenuEntry	( "Bump-Based (Medium)",	BUMPY_MEDIUM );
    glutAddMenuEntry	( "Bump-Based (Large)",		BUMPY_LARGE );
#endif	// OLD

    int menu_objective	= glutCreateMenu( selectItem );
    glutAddMenuEntry	( "Number Only",		NUMBER_ONLY );
    glutAddMenuEntry	( "Number+Area",		NUMBER_AREA );
    glutAddMenuEntry	( "Number+Runs",		NUMBER_RUN );
    // glutAddMenuEntry	( "Number+Stars",		NUMBER_STAR );
    // glutAddMenuEntry	( "Number+Weights",		NUMBER_WEIGHT );
    // glutAddMenuEntry	( "Number+Area+Runs",		NUMBER_AREA_RUN );
    // glutAddMenuEntry	( "Number+Area+Stars",		NUMBER_AREA_STAR );
    // glutAddMenuEntry	( "Number+Area+Biased Stars",	NUMBER_AREA_BIAS );
    // glutAddMenuEntry	( "Number+Area+Length",		NUMBER_AREA_LENGTH );
    // glutAddMenuEntry	( "Number+DualMax",		NUMBER_MAX );
    // glutAddMenuEntry	( "Number+Area+DualMax",	NUMBER_AREA_MAX );
    glutAddMenuEntry	( "Number+DualAverage",		NUMBER_AVE );
    glutAddMenuEntry	( "Number+Area+DualAverage",	NUMBER_AREA_AVE );
    // glutAddMenuEntry	( "Test Measure",		TEST_MEASURE );

    int menu_deviation	= glutCreateMenu( selectItem );
    glutAddMenuEntry	( "50%",			DEVIATION_050 );
    glutAddMenuEntry	( "55%",			DEVIATION_055 );
    glutAddMenuEntry	( "60%",			DEVIATION_060 );
    glutAddMenuEntry	( "65%",			DEVIATION_065 );
    glutAddMenuEntry	( "70%",			DEVIATION_070 );
    glutAddMenuEntry	( "75%(Default)",		DEVIATION_075 );
    glutAddMenuEntry	( "80%",			DEVIATION_080 );
    glutAddMenuEntry	( "85%",			DEVIATION_085 );

    int menu_stitch	= glutCreateMenu( selectItem );
    glutAddMenuEntry	( "[o] GA-based Stitch",	GA_STITCH );
    // glutAddMenuEntry	( "[t] Greedy Stitch",		GREEDY_STITCH );
    // glutAddMenuEntry	( "[v] One-Step Stitch",	ONESTEP_STITCH );
    // glutAddMenuEntry	( "[y] One-Step Remerge",	ONESTEP_REMERGE );
    glutAddMenuEntry	( "[g] Greedy Remerge",		GREEDY_REMERGE );
	
    int menu_weight	= glutCreateMenu( selectItem );
    glutAddMenuEntry    ( "Uniform",		UNIFORM_ASSIGNMENT );
    // glutAddMenuEntry    ( "Random",		RANDOM_ASSIGNMENT );
    glutAddMenuEntry    ( "Minimum Perimeter",	MINPERIMETER_ASSIGNMENT );
    // glutAddMenuEntry    ( "Maximum Perimeter",	MAXPERIMETER_ASSIGNMENT );
    glutAddMenuEntry    ( "Flat Spanning",	FLATSPANNING_ASSIGNMENT );
#ifdef DISPLAY_OPTION
    // glutAddMenuEntry    ( "Coil Spanning",	COILSPANNING_ASSIGNMENT );
    // glutAddMenuEntry    ( "Flat Region Split",	MINDIHEDRAL_ASSIGNMENT );
    // glutAddMenuEntry    ( "Sharp Region Split",	MAXDIHEDRAL_ASSIGNMENT );
    glutAddMenuEntry    ( "Minimum Blending",	MINBLENDING_ASSIGNMENT );
    // glutAddMenuEntry    ( "Maximum Blending",	MAXBLENDING_ASSIGNMENT );
    glutAddMenuEntry    ( "Minimum Curvature",	MINCURVATURE_ASSIGNMENT );
    // glutAddMenuEntry    ( "Maximum Curvatrre",	MAXCURVATURE_ASSIGNMENT );
    glutAddMenuEntry    ( "Minimum Concavity",	MINCONCAVITY_ASSIGNMENT );
    // glutAddMenuEntry    ( "Maximum Concavity",	MAXCONCAVITY_ASSIGNMENT );
    glutAddMenuEntry    ( "Hyperbolic Average",	HYPERBOLICAVE_ASSIGNMENT );
    glutAddMenuEntry    ( "Hyperbolic Maximum",	HYPERBOLICMAX_ASSIGNMENT );
    glutAddMenuEntry    ( "Hyperbolic Minimum",	HYPERBOLICMIN_ASSIGNMENT );
    glutAddMenuEntry    ( "Skeleton Flat",	SKELFLAT_ASSIGNMENT );
    // glutAddMenuEntry    ( "Skeleton Coil",	SKELCOIL_ASSIGNMENT );
#endif	// DISPLAY_OPTION

    int menu_capture	= glutCreateMenu( selectItem );
    glutAddMenuEntry	( "[c] Mesh Window",	CAPTURE_MESH );
    glutAddMenuEntry	( "[u] Sheet Window",	CAPTURE_SHEET );

    int menu_fileIO	= glutCreateMenu( selectItem );
    glutAddMenuEntry	( "[L] Load Mesh",	LOAD_MESH );
    glutAddMenuEntry	( "[S] Save Mesh",	SAVE_MESH );

//------------------------------------------------------------------------------
//	Main menu definition
//------------------------------------------------------------------------------
    // int menu_top	=
    glutCreateMenu	( selectItem );
    glutAddSubMenu	( "Render",		menu_render );
#ifdef DISPLAY_OPTION
    glutAddSubMenu	( "Saddle On/Off",	menu_saddle );
    glutAddSubMenu	( "Label On/Off",	menu_label );
    glutAddSubMenu	( "Transform",		menu_transform );
#endif	// DISPLAY_OPTION
    glutAddSubMenu	( "Split",		menu_split );
    glutAddSubMenu	( "Criteria for GA",	menu_objective );
    glutAddSubMenu	( "Stitch",		menu_stitch );
    glutAddSubMenu	( "Dual edge weights",	menu_weight );
    glutAddSubMenu	( "ST Deviation",	menu_deviation );
	
    glutAddMenuEntry	( "[m] MST",		MINIMUM_SPANNINGTREE);
    glutAddMenuEntry	( "Set viewpoint",	CHANGE_VIEWPOINT );

    glutAddSubMenu	( "File I/O",		menu_fileIO );
    glutAddSubMenu	( "Capture",		menu_capture );

    glutAddMenuEntry	( "Quit",		QUIT );
    glutAttachMenu( GLUT_LEFT_BUTTON );


//------------------------------------------------------------------------------
//	Sheet window
//------------------------------------------------------------------------------
    glutInitDisplayMode ( GLUT_DOUBLE | GLUT_RGB );
    glutInitWindowSize ( width_sheet, height_sheet ); 
    glutInitWindowPosition( 100, 50 );

#ifdef LARGE_WINDOWS
    glutInitWindowPosition( 450,  0 );
#endif	// LARGE_WINDOWS
#ifdef SCREENCAST_WINDOWS
    glutInitWindowPosition( 472,  0 );
#endif	// SCREENCAST_WINDOWS


    // glutInitWindowPosition( width_mesh+10-500, 0 );
    // glutInitWindowPosition( 376, 0 );
    win_sheet = glutCreateWindow ( "Sheet" );

    initSheet();

    glutDisplayFunc( displaySheet ); 
    glutReshapeFunc( reshapeSheet );
    // glutKeyboardFunc( keyboardSheet );
    glutKeyboardFunc( keyboardMesh );
    glutMouseFunc( mouseSheet );
    glutMotionFunc( motionSheet );


//------------------------------------------------------------------------------
//	Main Loop
//------------------------------------------------------------------------------
    checkInETime();
    checkInCPUTime();

    glutMainLoop();                       /* Start GLUT event-processing loop */
    
    return 0;
}
