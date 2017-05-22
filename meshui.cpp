//------------------------------------------------------------------------------
//
//	For Mesh Window 
//
//------------------------------------------------------------------------------
#include <sstream>

#include "common.h"
#include "ui.h"

//------------------------------------------------------------------------------
//	Open GL initialization
//------------------------------------------------------------------------------
void initMesh( void ) 
{
    glClearColor ( 0.8, 0.8, 0.8, 0.0 );
    /* glShadeModel ( GL_FLAT ); */
    /* glEnable( GL_LIGHTING ); */
    /* glEnable( GL_AUTO_NORMAL ); */
    /* glEnable( GL_NORMALIZE ); */

    /* Enable Z-buffer */
    glEnable( GL_AUTO_NORMAL );
    glEnable( GL_NORMALIZE );
    glEnable( GL_DEPTH_TEST );
    glEnable( GL_NORMALIZE );
    glDepthFunc( GL_LESS );

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    gluPerspective( fovy, aspect, near, far );

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

    initLights();
}


//------------------------------------------------------------------------------
//	Drawing functions
//------------------------------------------------------------------------------
void string3D( const char *str, double x, double y, double z )
{
    const double l = 1.20;
    const double s = 1.15;
    // glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
    // glColor4f( 1.0f, 1.0f, 0.2f, 1.0f );
    // glColor4f( 1.0f, 0.45f, 0.0f, 1.0f );
    
    // glRasterPos3f( 1.10*x, 1.10*y, 1.10*z );
    glRasterPos3f( l*x, l*y, l*z );

    // glDisable( GL_TEXTURE_2D );
    // glDisable( GL_DEPTH_TEST );
    // glDepthFunc( GL_ALWAYS );
    for (; *str != 0; str++) {
	// glutBitmapCharacter( FUTL_FONT_TYPE, *str );
	// GLUT_BITMAP_8_BY_13
	// GLUT_BITMAP_9_BY_15
	// GLUT_BITMAP_TIMES_ROMAN_10
	// GLUT_BITMAP_TIMES_ROMAN_24
	// GLUT_BITMAP_HELVETICA_10
	// GLUT_BITMAP_HELVETICA_12
	// GLUT_BITMAP_HELVETICA_18
	// glutBitmapCharacter( GLUT_BITMAP_HELVETICA_10, *str );
	glutBitmapCharacter( GLUT_BITMAP_HELVETICA_18, *str );
	// glutBitmapCharacter( GLUT_BITMAP_TIMES_ROMAN_24, *str );
    }
    // glEnable( GL_DEPTH_TEST );
    // glDepthFunc( GL_LESS );

    glLineWidth( 1.0 );
    glBegin( GL_LINES );
    glVertex3d( s*x, s*y, s*z );
    glVertex3d(   x,   y,   z );
    glEnd();
}


void label3D( int label, Point3 coord )
{
    ostringstream ostr;
    ostr << label << ends;
    string3D( ostr.str().c_str(), coord.x(), coord.y(), coord.z() );
}


void wireframe( GLenum mode )
{
    // for enabling antialiasing
    glEnable( GL_LINE_SMOOTH );
    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
    glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );
    glEnable( GL_DEPTH_TEST );

    glDisable( GL_LIGHTING );
    glLineWidth( 2.0 );
    glPointSize( 1.0 );

    glColor3d( 0.4, 0.4, 0.4 );

    for ( Facet_iterator fi = mesh.facets_begin(); fi != mesh.facets_end(); ++fi ) {
        Halfedge_facet_circulator hfc = fi->facet_begin();
        // Facets in polyhedral surfaces are at least triangles.
        CGAL_assertion( CGAL::circulator_size( hfc ) >= 3 );
	glBegin( GL_LINE_LOOP );
	do {
            glVertex3d( hfc->vertex()->point().x(),
			hfc->vertex()->point().y(),
			hfc->vertex()->point().z() );
        } while ( ++hfc != fi->facet_begin() );
	glEnd();
    }

    glDisable( GL_LINE_SMOOTH );
    glDisable( GL_BLEND );

    return;
}    


void hidden( GLenum mode )
{
    // for enabling antialiasing
    glEnable( GL_LINE_SMOOTH );
    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
    glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );

    glDisable( GL_LIGHTING );
    // glEnable( GL_DEPTH_TEST );
    glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
    // draw boundary lines of faces in black
    glPolygonOffset( 0.1, 0.1 );
    glColor3d( 0.0, 0.0, 0.0 );
    glLineWidth( 2.0 );
    glPointSize( 1.0 );

    if ( mode == GL_SELECT ) {
	glInitNames();
	glPushName( -1 );
    }

    for ( Facet_iterator fi = mesh.facets_begin(); fi != mesh.facets_end(); ++fi ) {
        Halfedge_facet_circulator hfc = fi->facet_begin();
        // Facets in polyhedral surfaces are at least triangles.
        CGAL_assertion( CGAL::circulator_size( hfc ) >= 3 );
	glBegin( GL_POLYGON );
	do {
            glVertex3d( hfc->vertex()->point().x(),
			hfc->vertex()->point().y(),
			hfc->vertex()->point().z() );
        } while ( ++hfc != fi->facet_begin() );
	glEnd();
    }

    // for disabling antialiasing
    glDisable( GL_LINE_SMOOTH );
    glDisable( GL_BLEND );
    glEnable( GL_DEPTH_TEST );

    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
    glEnable( GL_POLYGON_OFFSET_FILL );
    glPolygonOffset( 1.0, 1.0 );
    // draw interior of faces in gray
    glColor3d( 0.8, 0.8, 0.8 );

    for ( Facet_iterator fi = mesh.facets_begin(); fi != mesh.facets_end(); ++fi ) {
        Halfedge_facet_circulator hfc = fi->facet_begin();
        // Facets in polyhedral surfaces are at least triangles.
        CGAL_assertion( CGAL::circulator_size( hfc ) >= 3 );
	glBegin( GL_POLYGON );
	do {
            glVertex3d( hfc->vertex()->point().x(),
			hfc->vertex()->point().y(),
			hfc->vertex()->point().z() );
        } while ( ++hfc != fi->facet_begin() );
	glEnd();
    }

    glDisable( GL_POLYGON_OFFSET_FILL );

    if ( mode == GL_SELECT ) {
	glLineWidth( 1.0 );
	for ( Halfedge_iterator hi = mesh.halfedges_begin(); hi != mesh.halfedges_end(); ++hi ) {
	    glLoadName( hi->id() );
	    glBegin( GL_LINES );
	    glVertex3d( hi->vertex()->point().x(),
			hi->vertex()->point().y(),
			hi->vertex()->point().z() );
	    glVertex3d( hi->opposite()->vertex()->point().x(),
			hi->opposite()->vertex()->point().y(),
			hi->opposite()->vertex()->point().z() );
	    glEnd();
	    glLoadName( NO_INDEX );
	}
	glPopName();
	glInitNames();
    }


    return;
}    


#ifdef NO_NEED
void domain( GLenum mode )
{
    Graph * ptrG = &partial;

    // for enabling antialiasing
    glEnable( GL_LINE_SMOOTH );
    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
    glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );

    glDisable( GL_LIGHTING );
    // glEnable( GL_DEPTH_TEST );
    glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
    // draw boundary lines of faces in black
    glColor3d( 0.0, 0.6, 0.0 );
    glLineWidth( 2.0 );
    glPointSize( 1.0 );

    VertexIDMap VIDMap = get( vertex_index, (*ptrG) );
    EdgeIDMap EIDMap = get( edge_weight, (*ptrG) );
    graph_traits< Graph >::edge_iterator iterE, lastE;
    for ( tie( iterE, lastE ) = edges( (*ptrG) ); iterE != lastE; ++iterE ) {

#ifdef DEBUG
	glColor3d( 0.0, 0.6, 0.0 );
#endif	// DEBUG
        int srcIDF = VIDMap[ source( *iterE, (*ptrG) ) ];
        int tarIDF = VIDMap[ target( *iterE, (*ptrG) ) ];
        int midIDH = EIDMap[ *iterE ];
        glBegin( GL_LINE_STRIP );
        glVertex3d( bary[ srcIDF ].x(), bary[ srcIDF ].y(), bary[ srcIDF ].z() );
        glVertex3d( mid [ midIDH ].x(), mid [ midIDH ].y(), mid [ midIDH ].z() );
        glVertex3d( bary[ tarIDF ].x(), bary[ tarIDF ].y(), bary[ tarIDF ].z() );
        glEnd();

#ifdef DEBUG
	glColor4f( 0.0f, 0.0f, 0.8f, 1.0f );
	label3D( strIDF, bary[ srcIDF ] );
	label3D( tarIDF, bary[ tarIDF ] );
#endif	// DEBUG
    }


    // for disabling antialiasing
    glDisable( GL_LINE_SMOOTH );
    glDisable( GL_BLEND );
    glEnable( GL_DEPTH_TEST );

    return;
}    
#endif	// NO_NEED


void duality( GLenum mode )
{
    // for enabling antialiasing
    glEnable( GL_LINE_SMOOTH );
    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
    glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );

    glDisable( GL_LIGHTING );
    // glEnable( GL_DEPTH_TEST );
    glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
    // draw boundary lines of faces in black
    glColor3d( 0.2, 0.6, 0.2 );
    glLineWidth( 2.0 );
    glPointSize( 1.0 );

#ifdef SKIP
    if ( mode == GL_SELECT ) {
	glInitNames();
	// glLoadName( -1 );
	glPushName( -1 );
    }
#endif	// SKIP

    // initialize the halfedge IDs
    for ( Halfedge_iterator hi = mesh.halfedges_begin(); hi != mesh.halfedges_end(); ++hi ) {
	if ( hi->facet()->id() < hi->opposite()->facet()->id() ) {
	    if ( hi->label() == FIXED_EDGE ) {
		glColor3d( 1.0, 0.0, 0.0 );
		// cerr << " Drawing FIXED_EDGE No. " << hi->id() << endl;
	    }
	    else if ( hi->label() == INVALID_EDGE ) {
		glColor3d( 0.0, 0.0, 1.0 );
		// cerr << " Drawing INVALID_EDGE No. " << hi->id() << endl;
	    }
	    else {
		glColor3d( 0.3, 0.8, 0.3 );
	    }
	    Point3 * mid = &( hi->mid() );
	    Point3 * src = &( hi->facet()->center() );
	    Point3 * tar = &( hi->opposite()->facet()->center() );
#ifdef SKIP
	    if ( mode == GL_SELECT ) glLoadName( hi->id() );
#endif	// SKIP
	    glBegin( GL_LINE_STRIP );
	    glVertex3d( src->x(), src->y(), src->z() );
	    glVertex3d( mid->x(), mid->y(), mid->z() );
	    glVertex3d( tar->x(), tar->y(), tar->z() );
	    glEnd();
#ifdef SKIP
	    if ( mode == GL_SELECT ) glLoadName( NO_INDEX );
#endif	// SKIP
	}
    }

#ifdef SKIP
    if ( mode == GL_SELECT ) {
	glPopName();
	glInitNames();
    }
#endif	// SKIP

    // for disabling antialiasing
    glDisable( GL_LINE_SMOOTH );
    glDisable( GL_BLEND );
    glEnable( GL_DEPTH_TEST );

    return;
}    


void perfect( GLenum mode )
{
    // for enabling antialiasing
    glEnable( GL_LINE_SMOOTH );
    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
    glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );

    glDisable( GL_LIGHTING );
    // glEnable( GL_DEPTH_TEST );
    glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
    // draw boundary lines of faces in black
    glColor3d( 1.0, 0.8, 0.0 );
    glLineWidth( 2.0 );
    glPointSize( 1.0 );

    // initialize the halfedge IDs
    for ( Halfedge_iterator hi = mesh.halfedges_begin(); hi != mesh.halfedges_end(); ++hi ) {
	if ( hi->facet()->id() < hi->opposite()->facet()->id() ) {
	    // Matched edges (cut lines) are drawn in gray
	    if ( hi->match() == MATCHED_EDGE ) {
		glLineWidth( 2.0 );
		glColor3d( 0.5, 0.5, 0.5 );
	    }
	    // Unmatched edge (strip skeleton lines) are drawn in red
	    else {
		glLineWidth( 4.0 );
		//double rgb[ 3 ];
		//changeNumberToColor( 8, hi->cycle(), rgb );
		//glColor3dv( rgb );
		glColor3d( 0.8, 0.0, 0.0 );
	    }
	    Point3 * mid = &( hi->mid() );
	    Point3 * src = &( hi->facet()->center() );
	    Point3 * tar = &( hi->opposite()->facet()->center() );
	    glBegin( GL_LINE_STRIP );
	    glVertex3d( src->x(), src->y(), src->z() );
	    glVertex3d( mid->x(), mid->y(), mid->z() );
	    glVertex3d( tar->x(), tar->y(), tar->z() );
	    glEnd();
	}
    }

    // for disabling antialiasing
    glDisable( GL_LINE_SMOOTH );
    glDisable( GL_BLEND );
    glEnable( GL_DEPTH_TEST );

    return;
}    

// To display the KRUSKAL MINIMUM SPANNING TREE
void mst( GLenum mode )
{
    // for enabling antialiasing
    glEnable( GL_LINE_SMOOTH );
    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
    glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );

    glDisable( GL_LIGHTING );
    // glEnable( GL_DEPTH_TEST );
    glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
    // draw boundary lines of faces in black
    glColor3d( 1.0, 0.8, 0.0 );
    glLineWidth( 2.0 );
    glPointSize( 1.0 );

    // initialize the halfedge IDs
    for ( Halfedge_iterator hi = mesh.halfedges_begin(); hi != mesh.halfedges_end(); ++hi ) {
	if ( hi->facet()->id() < hi->opposite()->facet()->id() ) {
	    if ( hi->match() == KRUSKAL_EDGE ) {
		glLineWidth( 4.0 );
		glColor3d( 1.0, 0.0, 0.5 );
	    }
	    else {
		glLineWidth( 2.0 );
		glColor3d( 0.5, 0.5, 0.5 );
	    }
	    Point3 * mid = &( hi->mid() );
	    Point3 * src = &( hi->facet()->center() );
	    Point3 * tar = &( hi->opposite()->facet()->center() );
	    glBegin( GL_LINE_STRIP );
	    glVertex3d( src->x(), src->y(), src->z() );
	    glVertex3d( mid->x(), mid->y(), mid->z() );
	    glVertex3d( tar->x(), tar->y(), tar->z() );
	    glEnd();
	}
    }

    // for disabling antialiasing
    glDisable( GL_LINE_SMOOTH );
    glDisable( GL_BLEND );
    glEnable( GL_DEPTH_TEST );

    return;
}

void layout( GLenum mode )
{
    const double boundaryRatio = 0.3;

    // for enabling antialiasing
    glEnable( GL_LINE_SMOOTH );
    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
    glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );

    glDisable( GL_LIGHTING );
    // glEnable( GL_DEPTH_TEST );
    glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
    // draw boundary lines of faces in black
    glPolygonOffset( 0.1, 0.1 );
    glLineWidth( 2.0 );
    glPointSize( 1.0 );

//------------------------------------------------------------------------------
//	Draw each face boundary with thin black lines
//------------------------------------------------------------------------------
    glColor3d( 0.0, 0.0, 0.0 );
    for ( Facet_iterator fi = mesh.facets_begin(); fi != mesh.facets_end(); ++fi ) {
        Halfedge_facet_circulator hfc = fi->facet_begin();
        // Facets in polyhedral surfaces are at least triangles.
        // CGAL_assertion( CGAL::circulator_size( hfc ) >= NUM_SIKDS );
	glBegin( GL_POLYGON );
	do {
            glVertex3d( hfc->vertex()->point().x(),
			hfc->vertex()->point().y(),
			hfc->vertex()->point().z() );
        } while ( ++hfc != fi->facet_begin() );
	glEnd();
    }

//------------------------------------------------------------------------------
//	Drawing each face boundary with thick color lines
//------------------------------------------------------------------------------
    // for disabling antialiasing
    glDisable( GL_LINE_SMOOTH );
    glDisable( GL_BLEND );
    glEnable( GL_DEPTH_TEST );

    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
    glEnable( GL_POLYGON_OFFSET_FILL );
    glPolygonOffset( 1.0, 1.0 );

    // glDisable( GL_LIGHTING );
    // glShadeModel( GL_FLAT );
    for ( Facet_iterator fi = mesh.facets_begin(); fi != mesh.facets_end(); ++fi ) {

	Halfedge_facet_circulator hfc = fi->facet_begin();
	// Facets in polyhedral surfaces are at least triangles.
	// CGAL_assertion( CGAL::circulator_size( hfc ) >= NUM_SIDES );
	Point3	prevB, nextB, prevT, nextT, centerF;
	Vector3	prevV, nextV;
	double rgb[ 3 ];
	do {
	    prevB = hfc->vertex()->point();
	    nextB = hfc->next()->vertex()->point();
	    centerF = fi->center();
	    nextV = centerF - hfc->next()->vertex()->point();
	    prevV = centerF - hfc->vertex()->point();
	    nextT = hfc->next()->vertex()->point() + boundaryRatio * nextV;
	    prevT = hfc->vertex()->point() + boundaryRatio * prevV;
	    
	    rgb[ 0 ] = 0.8; rgb[ 1 ] = 0.8; rgb[ 2 ] = 0.8;
	    if ( fi->piece() != NO_INDEX ) {
		changeNumberToColor( 8, fi->piece(), rgb );
	    }
	    glColor3dv( rgb );
	    glBegin( GL_POLYGON );
	    glVertex3d( prevB.x(), prevB.y(), prevB.z() );
	    glVertex3d( nextB.x(), nextB.y(), nextB.z() );
	    glVertex3d( nextT.x(), nextT.y(), nextT.z() );
	    glVertex3d( prevT.x(), prevT.y(), prevT.z() );
	    glEnd();
	    
	    glColor3d( 0.8, 0.8, 0.8 );
	    glBegin( GL_POLYGON );
	    glVertex3d( prevT.x(), prevT.y(), prevT.z() );
	    glVertex3d( nextT.x(), nextT.y(), nextT.z() );
	    glVertex3d( centerF.x(), centerF.y(), centerF.z() );
	    glEnd();
	    
	    glFlush();
	} while ( ++hfc != fi->facet_begin() );
    }

    glPolygonOffset( 0.1, 0.1 );
    glDisable( GL_POLYGON_OFFSET_FILL );

//------------------------------------------------------------------------------
//	Draw each dual edge according to its connection status
//------------------------------------------------------------------------------
    EdgeDescriptor ed;
    // initialize the halfedge IDs
    for ( Halfedge_iterator hi = mesh.halfedges_begin(); hi != mesh.halfedges_end(); ++hi ) {

	// if the destination vertex has a smaller ID than the source vertex
	if ( hi->facet()->id() < hi->opposite()->facet()->id() ) {

	    // if 
	    // the end faces has the same piece ID &&
	    // the dual edge has the same cycle ID with the face piece ID &&
	    // the dual edge actually connects the two end faces &&
	    // the end faces has the valid piece ID
	    if ( ( hi->facet()->piece() == hi->opposite()->facet()->piece() ) &&
		 ( hi->cycle() == hi->facet()->piece() ) &&
		 ( hi->connect() ) &&
		 ( hi->facet()->piece() != NO_INDEX ) ) {

		glLineWidth( 4.0 );
		double rgb[ 3 ];
		changeNumberToColor( 8, hi->cycle(), rgb );
		glColor3dv( rgb );
		// glColor3d( 0.4, 0.4, 0.4 );
		Point3 * mid = &( hi->mid() );
		Point3 * src = &( hi->facet()->center() );
		Point3 * tar = &( hi->opposite()->facet()->center() );
		glBegin( GL_LINE_STRIP );
		glVertex3d( src->x(), src->y(), src->z() );
		glVertex3d( mid->x(), mid->y(), mid->z() );
		glVertex3d( tar->x(), tar->y(), tar->z() );
		glEnd();
	    }
	    else if ( ( ! hi->connect() ) &&
		      ( ( hi->vertex()->label() != SADDLE_VERTEX ) ||
			( ( hi->vertex()->label() == SADDLE_VERTEX ) &&
			  ( hi->vertex()->nCuts() >= 3 ) ) )
		      &&
		      ( ( hi->opposite()->vertex()->label() != SADDLE_VERTEX ) ||
			( ( hi->opposite()->vertex()->label() == SADDLE_VERTEX ) &&
			  ( hi->opposite()->vertex()->nCuts() >= 3 ) ) ) 
		      &&
		      ( ( hi->facet()->piece() != hi->opposite()->facet()->piece() ) &&
			( hi->facet()->piece() != NO_INDEX ) )
		      ) {
		glColor3d( 0.4, 0.4, 0.4 );
		Point3 * mid = &( hi->mid() );
		Point3 * src = &( hi->facet()->center() );
		Point3 * tar = &( hi->opposite()->facet()->center() );
		glBegin( GL_LINE_STRIP );
		glVertex3d( src->x(), src->y(), src->z() );
		glVertex3d( mid->x(), mid->y(), mid->z() );
		glVertex3d( tar->x(), tar->y(), tar->z() );
		glEnd();
	    }
	}
   }

    return;
}    


void saddle( GLenum mode )
{
    glShadeModel( GL_SMOOTH );
    glEnable( GL_LIGHTING );
    for ( Vertex_iterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi ) {
	Point3 * p = &( vi->point() );

	if ( vi->label() == SADDLE_VERTEX ) {
	    switch ( vi->nCuts() ) {
	      case 0:
		  grayMaterial(); 
		  break;
	      case 1:
		  blueMaterial();
		  break;
	      case 2:
		  cyanMaterial();
		  break;
	      case 3:
		  greenMaterial();
		  break;
	      case 4:
		  orangeMaterial();
		  break;
	      case 5:
		  redMaterial();
		  break;
	      case 6:
		  magentaMaterial();
		  break;
	      default:
		  grayMaterial();
		  // cerr << "=====> This saddle has more than 7 edge cuts!!" << endl;
		  // assert( false );
		  break;
	    }
	    glPushMatrix();
	    glTranslated( p->x(), p->y(), p->z() );
	    glutSolidSphere( 0.03, 12, 12 );
	    // glVertex3d( p->x(), p->y(), p->z() );
	    glPopMatrix();


	}

    }
    glDisable( GL_LIGHTING );
    // glEnd();

}


void annotate3D( GLenum mode )
{
    glColor4f( 0.0f, 0.5f, 1.0f, 0.0f );
    for ( Vertex_iterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi ) {
	label3D( vi->id(), vi->point() );
    }

    glColor4f( 1.0f, 0.0f, 0.0f, 0.0f );
    for ( Facet_iterator fi = mesh.facets_begin(); fi != mesh.facets_end(); ++fi) {
	label3D( fi->id(), fi->center() );	
    }
}



// To display the curve skeleton
void medial( GLenum mode )
{
    // for enabling antialiasing
    glEnable( GL_LINE_SMOOTH );
    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
    glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );
    glEnable( GL_DEPTH_TEST );

    glDisable( GL_LIGHTING );
    // glEnable( GL_DEPTH_TEST );
    // draw boundary lines of faces in black
    glColor3d( 1.0, 0.0, 0.0 );
    glLineWidth( 8.0 );

    glBegin( GL_LINES );
    for ( unsigned int k = 0; k < bone.size(); ++k ) {
	glVertex3d( bone[ k ].source().x(), 
		    bone[ k ].source().y(),
		    bone[ k ].source().z() );
	glVertex3d( bone[ k ].target().x(), 
		    bone[ k ].target().y(),
		    bone[ k ].target().z() );
    }
    glEnd();

    // for disabling antialiasing
    glDisable( GL_LINE_SMOOTH );
    glDisable( GL_BLEND );

    return;
}


#ifdef YET
void drawArrow3D( const Point3 & orig, const Point3 & dest, const double & scale )
{
    Vector3	dir = dest - orig;
    Vector3	unit = dir / sqrt( dir.squared_length() );
    double	theta, phi;
    double	arrow_length, arrow_x, arrow_y;
 
    if ( fabs( dir.z() ) < EPS ) {
        if ( dir.y() > 0.0 ) theta = 90.0;
	else theta = -90.0;
    }
    else {
        theta = 180.0/M_PI * atan2( dir.y(), dir.z() );
    }

    // Setting the arrow head
    arrow_length = scale * sqrt( dir.squared_length() );
    arrow_x = 0.1 * scale;
    arrow_y = arrow_x * tan( 0.349065850 );
 
    glPushMatrix();
 
    // Translate the source of the arrow to the coordinate oigin
    glTranslatef( orig.x(), orig.y(), orig.z() );	
    // Rotate the arrow along the z-axis
    // glRotatef( theta, 1.0, 0.0, 0.0 );
 
//	glColor3f(0.0, 0.8, 0.2);
    glBegin( GL_LINES );
    glVertex3d( 0.0, 0.0, 0.0 );
    glVertex3d( 0.0, 0.0, arrow_length );
    glEnd();

    // glutSolidCone( arrow_x, arrow_y, 32, 16 );
 
    glPopMatrix();
}
#endif	// YET



void cutedges( GLenum mode )
{
    int nColors;
    if ( reorder_flag ) nColors = attr.nColors();
    else nColors = DEFAULT_COLORS;

    // for enabling antialiasing
    glEnable( GL_LINE_SMOOTH );
    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
    glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );

    glDisable( GL_LIGHTING );
    // glEnable( GL_DEPTH_TEST );
    glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
    // draw boundary lines of faces in black
    glPolygonOffset( 0.1, 0.1 );
    glLineWidth( 2.0 );
    glPointSize( 1.0 );

#ifdef SKIP
//------------------------------------------------------------------------------
//	Draw each face boundary with thin black lines
//------------------------------------------------------------------------------
    glColor3d( 0.0, 0.0, 0.0 );
    for ( Facet_iterator fi = mesh.facets_begin(); fi != mesh.facets_end(); ++fi ) {
        Halfedge_facet_circulator hfc = fi->facet_begin();
        // Facets in polyhedral surfaces are at least triangles.
        // CGAL_assertion( CGAL::circulator_size( hfc ) >= 3 );
	glBegin( GL_POLYGON );
	do {
            glVertex3d( hfc->vertex()->point().x(),
			hfc->vertex()->point().y(),
			hfc->vertex()->point().z() );
        } while ( ++hfc != fi->facet_begin() );
	glEnd();
    }
#endif	// SKIP

//------------------------------------------------------------------------------
//	Draw cut paths
//------------------------------------------------------------------------------
    glEnable( GL_LINE_STIPPLE );

    EdgeDescriptor ed;
    // initialize the halfedge IDs
    for ( Halfedge_iterator hi = mesh.halfedges_begin(); hi != mesh.halfedges_end(); ++hi ) {

	// if the destination vertex has a smaller ID than the source vertex
	if ( hi->vertex()->id() < hi->opposite()->vertex()->id() ) {

	    if ( isConvex( hi ) ) {
		// glDisable( GL_LINE_STIPPLE );
		glLineStipple( 1, 0xFFFF );
		// glColor3d( 1.0, 0.0, 0.0 );
	    }
	    else {
		// glEnable( GL_LINE_STIPPLE );
		glLineStipple( 2, 0xAAAA );
		// glColor3d( 0.0, 1.0, 0.0 );
	    }
	    
	    if ( ! hi->connect() ) {
		glLineWidth( 8.0 );
		double rgb[ 3 ];
		changeNumberToColor( nColors, hi->path(), rgb );
		glColor3dv( rgb );
		// glColor3d( 0.0, 1.0, 0.0 );
		Point3 * orig = &( hi->vertex()->point() );
		Point3 * dest = &( hi->opposite()->vertex()->point() );
		glBegin( GL_LINES );
		glVertex3d( orig->x(), orig->y(), orig->z() );
		glVertex3d( dest->x(), dest->y(), dest->z() );
		glEnd();
		// drawArrow3D( *orig, *dest, 1.0 );
	    }
	    else {
		glLineWidth( 2.0 );
		// double rgb[ 3 ];
		// changeNumberToColor( 8, hi->cycle(), rgb );
		// glColor3dv( rgb );
		glColor3d( 0.0, 0.0, 0.0 );
		Point3 * orig = &( hi->vertex()->point() );
		Point3 * dest = &( hi->opposite()->vertex()->point() );
		glBegin( GL_LINES );
		glVertex3d( orig->x(), orig->y(), orig->z() );
		glVertex3d( dest->x(), dest->y(), dest->z() );
		glEnd();
	    }
	}
    }
    glDisable( GL_LINE_STIPPLE );


//------------------------------------------------------------------------------
//	Fill the interior of the faces
//------------------------------------------------------------------------------
    // for disabling antialiasing
    glDisable( GL_LINE_SMOOTH );
    glDisable( GL_BLEND );
    glEnable( GL_DEPTH_TEST );

    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
    glEnable( GL_POLYGON_OFFSET_FILL );
    glPolygonOffset( 1.0, 1.0 );
    // draw interior of faces in gray
    glColor3d( 0.8, 0.8, 0.8 );

    for ( Facet_iterator fi = mesh.facets_begin(); fi != mesh.facets_end(); ++fi ) {
        Halfedge_facet_circulator hfc = fi->facet_begin();
        // Facets in polyhedral surfaces are at least triangles.
        CGAL_assertion( CGAL::circulator_size( hfc ) >= 3 );
	glBegin( GL_POLYGON );
	do {
            glVertex3d( hfc->vertex()->point().x(),
			hfc->vertex()->point().y(),
			hfc->vertex()->point().z() );
        } while ( ++hfc != fi->facet_begin() );
	glEnd();
    }

    glDisable( GL_POLYGON_OFFSET_FILL );

    return;
}    
