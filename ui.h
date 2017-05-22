//#include <GL/glut.h> mac virsion is as follow
#include <GLUT/glut.h>

//------------------------------------------------------------------------------
//	Geometric Variables
//------------------------------------------------------------------------------
// original triangular mesh
extern Polyhedron				mesh;
// dual graph of the original mesh
extern Graph					whole;
// partial dual graph for computing the perfect matching
extern Graph					partial;
// connectivity graph
// face barycenters [FID]
extern vector< Point3 >				bary;
//  midpoints of the original edges [EID]
extern vector< Point3 >				mid;
// lengths of the original edges [EID]
//extern vector< double >				length; // <- to be removed
// array of dual cycles
extern vector< vector< Halfedge_handle > >	skeleton;
extern vector< vector< Facet_handle > >		pattern;
// array of 2D bounding boxes
extern vector< Bbox2 >				bound;
// paper size for mesh undolding	
extern Bbox2					paper;
// Boundary cut path of the unfoled patches
extern Attribute				attr;

// curve skeleton
extern vector< Segment3 >			bone;
extern vector< Point3 >				refer;


//------------------------------------------------------------------------------
//	UI Variables
//------------------------------------------------------------------------------

// Mesh window
extern int					win_mesh;

// window sizes
extern int					width_mesh;
extern int					height_mesh;

// projection parameters
extern double					fovy;
extern double					aspect;
extern double					near;
extern double					far;

extern double					view_distance;
extern double					view_azimuth;
extern double					view_incidence;
extern double					view_twist;

extern double					translate_x, translate_y, translate_z;

// mouse positions and its button status
extern int					left_button, middle_button, right_button;
extern int					pointer_x, pointer_y;

// drawing mode
extern int					drawing_mode;
extern bool					saddle_flag;
extern bool					label_flag;
extern bool					number_flag;
extern bool					reorder_flag;


// Sheet window

// window IDs
extern int					win_sheet;

// window sizes
extern int					width_sheet;
extern int					height_sheet;


extern double					mstDeviation;


//------------------------------------------------------------------------------
//	Functions in meshui.cpp
//------------------------------------------------------------------------------
extern void	initMesh	( void );
extern void	wireframe	( GLenum mode );
extern void	hidden		( GLenum mode );
extern void	duality		( GLenum mode );
extern void	perfect		( GLenum mode );
extern void	layout		( GLenum mode );
extern void	domain		( GLenum mode );
extern void	connectivity	( GLenum mode );
extern void	saddle		( GLenum mode );
extern void	annotate3D	( GLenum mode );

extern void	mst		( GLenum mode );

extern void	medial		( GLenum mode );

extern void	cutedges	( GLenum mode );

//------------------------------------------------------------------------------
//	Functions in sheetui.cpp
//------------------------------------------------------------------------------
extern void	initSheet	( void );
extern void	strip		( GLenum mode ); 
extern void	craft		( GLenum mode );
extern void	outline		( GLenum mode );
extern void	folding		( GLenum mode );

//------------------------------------------------------------------------------
//	Functions in render.cpp
//------------------------------------------------------------------------------
extern void	initLights	( void );
extern void	redMaterial	( void );
extern void	orangeMaterial	( void );
extern void	greenMaterial	( void );
extern void	cyanMaterial	( void );
extern void	blueMaterial	( void );
extern void	magentaMaterial	( void );
extern void	grayMaterial	( void );
extern void	changeNumberToColor
				( int step, int number, double rgb[ 3 ] );
