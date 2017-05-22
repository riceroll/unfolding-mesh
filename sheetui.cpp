//------------------------------------------------------------------------------
//
//	For Sheet Window 
//
//------------------------------------------------------------------------------
#include <iomanip>
#include <sstream>
#include <cmath>

#include "common.h"
#include "ui.h"

#ifdef LABEL_INSIDE_STRIP
#define	DISPLACEMENT_RATIO	(0.1)
#define	ANNOTATION_RATIO	(0.4)
#else	// LABEL_INSIDE_STRIP
#define	DISPLACEMENT_RATIO	(0.0)
#define	ANNOTATION_RATIO	(1.0)
#endif	// LABEL_INSIDE_STRIP

//------------------------------------------------------------------------------
//	Open GL initialization
//------------------------------------------------------------------------------
void initSheet( void ) 
{
    glClearColor ( 1.0, 1.0, 1.0, 0.0 );
    /* glShadeModel ( GL_FLAT ); */
    /* glEnable( GL_LIGHTING ); */
    /* glEnable( GL_AUTO_NORMAL ); */
    /* glEnable( GL_NORMALIZE ); */

    /* Enable Z-buffer */
    glEnable( GL_AUTO_NORMAL );
    glEnable( GL_NORMALIZE );
    // glEnable( GL_DEPTH_TEST );
    // glEnable( GL_NORMALIZE );
    // glDepthFunc( GL_LESS );

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    // gluPerspective( fovy, aspect, near, far );
    gluOrtho2D( paper.xmin(), paper.xmax(), paper.ymin(), paper.ymax() );
    
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
}


//------------------------------------------------------------------------------
//	Drawing functions
//------------------------------------------------------------------------------
void string2D( const char *str, double x, double y )
{
    // glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
    // glColor4f( 1.0f, 1.0f, 0.2f, 1.0f );
    // glColor4f( 1.0f, 0.45f, 0.0f, 1.0f );
    
    glRasterPos2f( x, y );
    // glRasterPos2f( x-0.05, y-0.02 );

    for (; *str != 0; str++) {
	// glutBitmapCharacter( FUTL_FONT_TYPE, *str );
	// GLUT_BITMAP_8_BY_13
	// GLUT_BITMAP_9_BY_15
	// GLUT_BITMAP_TIMES_ROMAN_10
	// GLUT_BITMAP_TIMES_ROMAN_24
	// GLUT_BITMAP_HELVETICA_10
	// GLUT_BITMAP_HELVETICA_12
	// GLUT_BITMAP_HELVETICA_18
	// glutBitmapCharacter( GLUT_BITMAP_8_BY_13, *str );
	glutBitmapCharacter( GLUT_BITMAP_HELVETICA_10, *str );
	// glutBitmapCharacter( GLUT_BITMAP_HELVETICA_18, *str );
	// glutBitmapCharacter( GLUT_BITMAP_TIMES_ROMAN_24, *str );
    }
}


void label2D( int label, const Point2 coord )
{
    ostringstream ostr;
    ostr << " " << label << ends;
    // ostr << " " << setw( 3 ) << label << ends;
    string2D( ostr.str().c_str(), coord.x(), coord.y() );
}


Point2 midpointTriangle( const Halfedge_handle & hh )
{
    Point2 coordS, coordT, midpoint;

    if ( hh == hh->facet()->halfedge() ) {
	coordS = hh->facet()->triangle()[ 0 ];
	coordT = hh->facet()->triangle()[ 1 ];
    }
    else if ( hh == hh->facet()->halfedge()->next() ) {
	// right shoulder
	coordS = hh->facet()->triangle()[ 1 ];
	coordT = hh->facet()->triangle()[ 2 ];
    }
    else if ( hh == hh->facet()->halfedge()->prev() ) {
	// left shoulder
	coordS = hh->facet()->triangle()[ 2 ];
	coordT = hh->facet()->triangle()[ 0 ];
    }
    else assert( false );

    midpoint = CGAL::ORIGIN + 0.5 * ( ( coordS - CGAL::ORIGIN ) + ( coordT - CGAL::ORIGIN ) );

    return midpoint;
}


void colorByCuts( int nCuts )
{
    switch ( nCuts ) {
      case 0:			// gray
	  glColor3d( 0.6, 0.6, 0.6 );
	  break;
      case 1:			// blue
	  glColor3d( 0.0, 0.0, 1.0 );
	  break;
      case 2:			// cyan
	  glColor3d( 0.0, 0.6, 1.0 );
	  break;
      case 3:			// green
	  glColor3d( 0.0, 1.0, 0.0 );
	  break;
      case 4:			// orange
	  glColor3d( 1.0, 0.6, 0.0 );
	  break;
      case 5:			// red
	  glColor3d( 1.0, 0.0, 0.0 );
	  break;
      case 6:			// magenta
	  glColor3d( 1.0, 0.0, 1.0 );
	  break;
      default:
	  glColor3d( 0.6, 0.6, 0.6 );
	  // cerr << "=====> This saddle has more than 7 edge cuts!!" << endl;
	  // assert( false );
    }
}


void possibleStitch( GLenum mode, const Halfedge_handle & hh )
{
    if ( isLocally( hh ) ) {
	Point2 coord = centerTriangle( hh->facet()->triangle() );
	Point2 midpoint = midpointTriangle( hh );
	glBegin( GL_LINES );
	glVertex2d( coord.x(), coord.y() );
	glVertex2d( midpoint.x(), midpoint.y() );
	glEnd();
    }
}


void strip( GLenum mode )
{
    glEnable( GL_LINE_SMOOTH );
    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
    glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );

    // cerr << " Number of unfolded patterns = " << pattern.size() << endl;
    for ( unsigned int i = 0; i < pattern.size(); ++i ) {

//------------------------------------------------------------------------------
//	Draw triangle sides
//------------------------------------------------------------------------------
	glLineWidth( 2.0 );
	double rgb[ 3 ];
	changeNumberToColor( 8, i, rgb );
	glColor3dv( rgb );
	for ( unsigned int j = 0; j < pattern[ i ].size(); ++j ) {
	    glBegin( GL_LINE_LOOP );
	    for ( unsigned int k = 0; k < NUM_SIDES; ++k )
		glVertex2d( pattern[ i ][ j ]->triangle()[ k ].x(), 
			    pattern[ i ][ j ]->triangle()[ k ].y() );
	    glEnd();
	}

//------------------------------------------------------------------------------
//	Draw the skeleton
//------------------------------------------------------------------------------
	glLineWidth( 4.0 );
	glColor3dv( rgb );
	Point2 coordI, coordO;
	Halfedge_handle hv[ NUM_SIDES ];
	for ( unsigned int j = 0; j < pattern[ i ].size(); ++j ) {
	    // calculate the barycenter of the triangle pattern[ i ][ j ];
	    coordI = centerTriangle( pattern[ i ][ j ]->triangle() );
	    // pointer to the base edge of the triangle
	    Halfedge_handle baseH = pattern[ i ][ j ]->halfedge();
	    // pointers to the three side edges
	    hv[ 0 ] = baseH->prev();
	    hv[ 1 ] = baseH;
	    hv[ 2 ] = baseH->next();
	    for ( unsigned int k = 0; k < NUM_SIDES; ++k ) {
		// If this unfolded pattern includes this edge and connects the
		// pair of adjacent triangles along that edge
		if ( ( hv[ k ]->cycle() == ( int )i ) &&
		     ( hv[ k ]->connect() ) &&
		     // if the ID of this face is smaller then the other one
		     ( pattern[ i ][ j ]->id() < hv[ k ]->opposite()->facet()->id() ) ) {
		    // calculate the barycenter of the other opponent triangle
		    coordO = centerTriangle( hv[ k ]->opposite()->facet()->triangle() );
		    // Draw the line connecting two barycenters
		    glBegin( GL_LINES );
		    glVertex2d( coordI.x(), coordI.y() );
		    glVertex2d( coordO.x(), coordO.y() );
		    glEnd();
		}
	    }
	}

//------------------------------------------------------------------------------
//	Draw the possible connection
//------------------------------------------------------------------------------
	glLineWidth( 2.0 );
	glColor3d( 0.6, 0.6, 0.6 );
	int N = ( int )pattern[ i ].size();
	for ( int j = 0; j < N; ++j ) {
	    // pointer to the base edge of the triangle
	    Halfedge_handle baseH = pattern[ i ][ j ]->halfedge();
	    // draw a line if the patterns can be stitched along that edge
	    possibleStitch( mode, baseH->prev() );
	    possibleStitch( mode, baseH );
	    possibleStitch( mode, baseH->next() );
	}

//------------------------------------------------------------------------------
//	Mark saddle verteices
//------------------------------------------------------------------------------
	if ( saddle_flag ) {
	    glPointSize( 8.0 );
	    Halfedge_handle hv[ NUM_SIDES ];
	    for ( unsigned int j = 0; j < pattern[ i ].size(); ++j ) {
		// pointer to the base edge of the triangle
		Halfedge_handle baseH = pattern[ i ][ j ]->halfedge();
		// pointers to the three side edges
		hv[ 0 ] = baseH->prev();
		hv[ 1 ] = baseH;
		hv[ 2 ] = baseH->next();
		for ( unsigned int k = 0; k < NUM_SIDES; ++k ) {
		    if ( hv[ k ]->vertex()->label() == SADDLE_VERTEX ) {
			colorByCuts( hv[ k ]->vertex()->nCuts() );
			glBegin( GL_POINTS );
			glVertex2d( pattern[ i ][ j ]->triangle()[ k ].x(), 
				    pattern[ i ][ j ]->triangle()[ k ].y() ); 
			glEnd();
		    }
		}
	    }
	}

//------------------------------------------------------------------------------
//	Draw vertex ID labels
//------------------------------------------------------------------------------
	if ( label_flag ) {
	    Halfedge_handle hv[ NUM_SIDES ];
	    for ( unsigned int j = 0; j < pattern[ i ].size(); ++j ) {
		// pointer to the base edge of the triangle
		Halfedge_handle baseH = pattern[ i ][ j ]->halfedge();
		// pointers to the three side edges
		hv[ 0 ] = baseH->prev();
		hv[ 1 ] = baseH;
		hv[ 2 ] = baseH->next();

		glColor3d( 1.0, 0.0, 0.0 );
		Point2 center = centerTriangle( pattern[ i ][ j ]->triangle() );
		Vector2	move = 0.02 * Vector2( -2.0, -0.5 );
		Point2 plot = center + move;
		label2D( hv[ 0 ]->facet()->id(), plot );

		glColor3d( 0.0, 0.0, 0.0 );
		for ( unsigned int k = 0; k < NUM_SIDES; ++k ) {
		    Vector2	span = pattern[ i ][ j ]->triangle()[ k ] - center;
		    Vector2	shift = DISPLACEMENT_RATIO * sqrt( span.squared_length() ) * Vector2( -2.0, -0.5 );
		    Point2	coord = center + ANNOTATION_RATIO * span + shift;
		    label2D( hv[ k ]->vertex()->id(), coord );
		    // label2D( hv[ k ]->vertex()->id(), pattern[ i ][ j ]->triangle()[ k ] );
		}

	    }
	}
//------------------------------------------------------------------------------
//	Draw the bounding box
//------------------------------------------------------------------------------
	glColor3d( 0.7, 0.7, 0.7 );
	glLineWidth( 1.0 );
	glBegin( GL_LINE_LOOP );
	glVertex2d( bound[ i ].xmin(), bound[ i ].ymin() );
	glVertex2d( bound[ i ].xmax(), bound[ i ].ymin() );
	glVertex2d( bound[ i ].xmax(), bound[ i ].ymax() );
	glVertex2d( bound[ i ].xmin(), bound[ i ].ymax() );
	glEnd();
    }

    glDisable( GL_LINE_SMOOTH );
    glDisable( GL_BLEND );

}


void drawArrow2D( const Point2 & orig, const Point2 & dest, const double & scale )
{
    Vector2	dir = dest - orig;
    double	theta;
    double	arrow_length, arrow_x, arrow_y;
 
    if ( fabs( dir.x() ) < EPS ) {
        if ( dir.y() > 0.0 ) theta = 90.0;
	else theta = -90.0;
    }
    else {
        theta = 180.0/M_PI * atan2( dir.y(), dir.x() );
    }
    
    // Setting the arrow head
    arrow_length = scale * sqrt( dir.squared_length() );
    // If the size of the arrow head is too small,
    // fix its size by referring to the length of the corresponding arrow
    // if ( arrow_length < scale * 0.2 ) {
    // arrow_x = 0.2 * arrow_length;
    //}
    // Use the fixed-sized arrow head if the corresponding arrow is enough long
    // else {                        
    // arrow_x = 0.2 * scale;
    arrow_x = 0.1 * scale; // <-
    // arrow_x = 0.3 * scale;
    // }
    // arrow_y = arrow_x * tan( 0.174532925 );
    arrow_y = arrow_x * tan( 0.349065850 );
 
    glPushMatrix();
 
    // Translate the source of the arrow to the coordinate oigin
    glTranslatef( orig.x(), orig.y(), 0.0 );	
    // Rotate the arrow along the z-axis
    glRotatef( theta, 0.0, 0.0, 1.0 );
 
//	glColor3f(0.0, 0.8, 0.2);
    glBegin( GL_LINES );
    glVertex2d( 0.0, 0.0 );
    glVertex2d( arrow_length, 0.0 );
    glVertex2d( arrow_length, 0.0 );
    glVertex2d( arrow_length - arrow_x,  arrow_y );
    glVertex2d( arrow_length, 0.0 );
    glVertex2d( arrow_length - arrow_x, -arrow_y );
    glEnd();
 
    glPopMatrix();
}



void craft( GLenum mode )
{
    // 100 faces
    // const double diffx = -0.015, diffy = 0.030, diffr = 0.035;
    // 300 faces
    const double diffx = -0.015, diffy = 0.020, diffr = 0.030;
    Halfedge_handle hv[ NUM_SIDES ];
    int nColors;
    if ( reorder_flag ) nColors = attr.nColors();
    else nColors = DEFAULT_COLORS;

    glEnable( GL_LINE_SMOOTH );
    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
    glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );

    // cerr << " Number of unfolded patterns = " << pattern.size() << endl;
    glEnable( GL_LINE_STIPPLE );
    // glDisable( GL_LINE_STIPPLE );
    for ( unsigned int i = 0; i < pattern.size(); ++i ) {
//------------------------------------------------------------------------------
//	Draw triangle sides
//------------------------------------------------------------------------------
	for ( unsigned int j = 0; j < pattern[ i ].size(); ++j ) {
	    // Adjust halfedge
	    // f->halfedge() corresponds to triangle()[ 0 ] - triangle()[ 1 ]
	    // f->halfedge()->vertex() corresponds to triangle()[ 1 ];
	    Halfedge_handle baseH = pattern[ i ][ j ]->halfedge();
	    hv[ 0 ] = baseH;
	    hv[ 1 ] = baseH->next();
	    hv[ 2 ] = baseH->prev();

	    for ( unsigned int k = 0; k < NUM_SIDES; ++k ) {

		// drawType	0 : do nothing
		//		1 : internal edge
		//		2 : boundary edge
		int drawType;
		// If the edge is on the boundary
		if ( hv[ k ]->path() != NO_INDEX ) {
		    double rgb[ 3 ];
		    // changeNumberToColor( 8, hv[ k ]->path(), rgb );
		    changeNumberToColor( nColors, hv[ k ]->path(), rgb );
		    glColor3dv( rgb );
		    drawType = 2;
		}
		// If the edge is internal
		else {
		    // draw only when the source vertex ID is less than that of
		    // the target vertex
		    if ( hv[ k ]->opposite()->vertex()->id() < hv[ k ]->vertex()->id() ) {
			drawType = 1;
			glColor3d( 0.0, 0.0, 0.0 );
		    }
		    else drawType = 0;
		}

		// If the edge should be drawn, ...
		if ( drawType > 0 ) {
		    if ( isConvex( hv[ k ] ) ) {
			// glDisable( GL_LINE_STIPPLE );
			glLineStipple( 1, 0xFFFF );
			// glColor3d( 1.0, 0.0, 0.0 );
		    }
		    else {
			// glEnable( GL_LINE_STIPPLE );
			glLineStipple( 2, 0xAAAA );
			// glColor3d( 0.0, 1.0, 0.0 );
		    }
		
		    switch ( drawType ) {
		      case 1:	// internal edge
			  glLineWidth( 2.0 );
			  glBegin( GL_LINES );
			  glVertex2d( pattern[ i ][ j ]->triangle()[ k%NUM_SIDES ].x(), 
				      pattern[ i ][ j ]->triangle()[ k%NUM_SIDES ].y() );
			  glVertex2d( pattern[ i ][ j ]->triangle()[ (k+1)%NUM_SIDES ].x(), 
				      pattern[ i ][ j ]->triangle()[ (k+1)%NUM_SIDES ].y() );
			  glEnd();
			  break;
		      case 2:	// boundary edge
			  glLineWidth( 4.0 );
			  if ( ( hv[ k ]->orient() ) && 
			       ( hv[ k ]->vertex()->nCuts() != 2 ) ) {
			      drawArrow2D( pattern[ i ][ j ]->triangle()[ k%NUM_SIDES ],
					   pattern[ i ][ j ]->triangle()[ (k+1)%NUM_SIDES ],
					   1.00 );
			  }
			  else if ( ( ! hv[ k ]->orient() ) &&
				    ( hv[ k ]->opposite()->vertex()->nCuts() != 2 ) ) {
			      drawArrow2D( pattern[ i ][ j ]->triangle()[ (k+1)%NUM_SIDES ],
					   pattern[ i ][ j ]->triangle()[ k%NUM_SIDES ],
					   1.00 );
			  }
			  else {
			      glBegin( GL_LINES );
			      glVertex2d( pattern[ i ][ j ]->triangle()[ k%NUM_SIDES ].x(), 
					  pattern[ i ][ j ]->triangle()[ k%NUM_SIDES ].y() );
			      glVertex2d( pattern[ i ][ j ]->triangle()[ (k+1)%NUM_SIDES ].x(), 
					  pattern[ i ][ j ]->triangle()[ (k+1)%NUM_SIDES ].y() );
			      glEnd();
			  }
			  if ( number_flag ) {
			      Point2 mid = 
				  CGAL::ORIGIN +
				  0.5 * ( pattern[ i ][ j ]->triangle()[ k%NUM_SIDES ] - CGAL::ORIGIN ) +
				  0.5 * ( pattern[ i ][ j ]->triangle()[ (k+1)%NUM_SIDES ] - CGAL::ORIGIN );
			      Point2 cen = centerTriangle( pattern[ i ][ j ]->triangle() );
			      Vector2 dir = mid - cen;
			      Vector2 unit = dir / sqrt( dir.squared_length() );
			      Point2 plot = mid - diffr * unit;
			      ostringstream ostr;
			      ostr << setw( 2 ) << hv[ k ]->path() << ends;
			      // int nC = strlen( ostr.str().c_str() );
			      // glPushMatrix();
			      // glRotatef( 60.0, 0.0, 0.0, 1.0 );
			      // cerr << " nC = " << nC << endl;
			      glColor3d( 0.0, 0.0, 0.0 );
			      string2D( ostr.str().c_str(), plot.x()+diffx, plot.y()+diffy );
			      // string2D( ostr.str().c_str(), plot.x(), plot.y() );
			      // string2D( ostr.str().c_str(), mid.x()+diffx,  mid.y()+diffy );
			      // glPopMatrix();
			      glLineWidth( 1.0 );
			      glLineStipple( 1, 0xFFFF );
			      glBegin( GL_LINES );
			      glVertex2d( plot.x()+1.5*diffx, plot.y()+diffy );
			      glVertex2d( plot.x()+1.5*diffx, plot.y()-diffy );
			      glEnd();
			  }
			  break;
		    }
		}
	    }
	}

//------------------------------------------------------------------------------
//	DRAW vertex ID labels
//------------------------------------------------------------------------------
	if ( label_flag ) {
	    Halfedge_handle hv[ NUM_SIDES ];
	    for ( unsigned int j = 0; j < pattern[ i ].size(); ++j ) {
		// pointer to the base edge of the triangle
		Halfedge_handle baseH = pattern[ i ][ j ]->halfedge();
		// pointers to the three side edges
		hv[ 0 ] = baseH;
		hv[ 1 ] = baseH->next();
		hv[ 2 ] = baseH->prev();

		glColor3d( 1.0, 0.0, 0.0 );
		Point2 center = centerTriangle( pattern[ i ][ j ]->triangle() );
		Vector2	move = 0.02 * Vector2( -2.0, -0.5 );
		Point2 plot = center + move;
		label2D( hv[ 0 ]->facet()->id(), plot );

		glColor3d( 0.0, 0.0, 1.0 );
		for ( unsigned int k = 0; k < NUM_SIDES; ++k ) {
		    Vector2	span = pattern[ i ][ j ]->triangle()[ k ] - center;
		    Vector2	shift = DISPLACEMENT_RATIO * sqrt( span.squared_length() ) * Vector2( -2.0, -0.5 );
		    Point2	coord = center + ANNOTATION_RATIO * span + shift;
		    label2D( hv[ k ]->opposite()->vertex()->id(), coord );
		    // label2D( hv[ k ]->vertex()->id(), pattern[ i ][ j ]->triangle()[ k ] );
		}

	    }
	}

//------------------------------------------------------------------------------
//	Draw the bounding box
//------------------------------------------------------------------------------
//#ifdef SKIP
	glColor3d( 0.7, 0.7, 0.7 );
	glLineWidth( 1.0 );
	glBegin( GL_LINE_LOOP );
	glVertex2d( bound[ i ].xmin(), bound[ i ].ymin() );
	glVertex2d( bound[ i ].xmax(), bound[ i ].ymin() );
	glVertex2d( bound[ i ].xmax(), bound[ i ].ymax() );
	glVertex2d( bound[ i ].xmin(), bound[ i ].ymax() );
	glEnd();
	//#endif	// SKIP
    }
    glDisable( GL_LINE_STIPPLE );

//------------------------------------------------------------------------------
//	Mark the dead ends of cut paths
//------------------------------------------------------------------------------
    for ( unsigned int i = 0; i < pattern.size(); ++i ) {
	glColor3d( 0.0, 0.0, 0.0 );
	glPointSize( 10.0 );
	glBegin( GL_POINTS );
	for ( unsigned int j = 0; j < pattern[ i ].size(); ++j ) {
	    Halfedge_handle baseH = pattern[ i ][ j ]->halfedge();
	    hv[ 0 ] = baseH;
	    hv[ 1 ] = baseH->next();
	    hv[ 2 ] = baseH->prev();

	    for ( unsigned int k = 0; k < NUM_SIDES; ++k ) {
		// If the edge is on the boundary
		if ( hv[ k ]->path() != NO_INDEX ) {
		    // Check the source vertex of the boundary edge for marking
		    // dead ends of the cut paths
		    int nCuts = hv[ k ]->opposite()->vertex()->nCuts();
		    if ( ( nCuts == 1 ) || ( nCuts >= 3 ) ) {
			glVertex2d( pattern[ i ][ j ]->triangle()[ k%NUM_SIDES ].x(), 
				    pattern[ i ][ j ]->triangle()[ k%NUM_SIDES ].y() );
		    }
		}
	    }
	}
	glEnd();
    }

    glDisable( GL_LINE_SMOOTH );
    glDisable( GL_BLEND );
}


void outline( GLenum mode )
{
    Halfedge_handle hv[ NUM_SIDES ];
    int nColors;
    if ( reorder_flag ) nColors = attr.nColors();
    else nColors = DEFAULT_COLORS;

    glEnable( GL_LINE_SMOOTH );
    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
    glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );

    glEnable( GL_LINE_STIPPLE );

    for ( unsigned int i = 0; i < pattern.size(); ++i ) {
//------------------------------------------------------------------------------
//	Draw triangle sides
//------------------------------------------------------------------------------
	for ( unsigned int j = 0; j < pattern[ i ].size(); ++j ) {
	    // Adjust halfedge
	    // f->halfedge() corresponds to triangle()[ 0 ] - triangle()[ 1 ]
	    // f->halfedge()->vertex() corresponds to triangle()[ 1 ];
	    Halfedge_handle baseH = pattern[ i ][ j ]->halfedge();
	    hv[ 0 ] = baseH;
	    hv[ 1 ] = baseH->next();
	    hv[ 2 ] = baseH->prev();

	    for ( unsigned int k = 0; k < NUM_SIDES; ++k ) {

		// drawType	0 : do nothing
		//		1 : internal edge
		//		2 : boundary edge
		int drawType;
		// If the edge is on the boundary
		if ( hv[ k ]->path() != NO_INDEX ) {
		    drawType = 2;
		}
		// If the edge is internal
		else {
		    // draw only when the source vertex ID is less than that of
		    // the target vertex
		    if ( hv[ k ]->opposite()->vertex()->id() < hv[ k ]->vertex()->id() ) {
			drawType = 1;
		    }
		    else drawType = 0;
		}

		if ( isConvex( hv[ k ] ) ) {
		    // glDisable( GL_LINE_STIPPLE );
		    glLineStipple( 1, 0xFFFF );
		    // glColor3d( 1.0, 0.0, 0.0 );
		}
		else {
		    // glEnable( GL_LINE_STIPPLE );
		    glLineStipple( 2, 0xAAAA );
		    // glLineStipple( 2, 0xF0F0 );
		    // glColor3d( 0.0, 1.0, 0.0 );
		}

		// If the edge should be drawn, ...
		switch ( drawType ) {
		  case 0:	// do nothing
		  case 1:	// internal edge
		      break;
		  case 2:	// boundary edge
		      glLineWidth( 4.0 );
		      glColor3d( 0.0, 0.0, 0.0 );
		      glBegin( GL_LINES );
		      glVertex2d( pattern[ i ][ j ]->triangle()[ k%NUM_SIDES ].x(), 
				  pattern[ i ][ j ]->triangle()[ k%NUM_SIDES ].y() );
		      glVertex2d( pattern[ i ][ j ]->triangle()[ (k+1)%NUM_SIDES ].x(), 
				  pattern[ i ][ j ]->triangle()[ (k+1)%NUM_SIDES ].y() );
		      glEnd();
		      break;
		}
	    }
	}

//------------------------------------------------------------------------------
//	Draw the bounding box
//------------------------------------------------------------------------------
	glColor3d( 0.7, 0.7, 0.7 );
	glLineWidth( 1.0 );
	glBegin( GL_LINE_LOOP );
	glVertex2d( bound[ i ].xmin(), bound[ i ].ymin() );
	glVertex2d( bound[ i ].xmax(), bound[ i ].ymin() );
	glVertex2d( bound[ i ].xmax(), bound[ i ].ymax() );
	glVertex2d( bound[ i ].xmin(), bound[ i ].ymax() );
	glEnd();
    }

    glDisable( GL_LINE_SMOOTH );
    glDisable( GL_BLEND );
}


void folding( GLenum mode )
{
    Halfedge_handle hv[ NUM_SIDES ];
    int nColors;
    if ( reorder_flag ) nColors = attr.nColors();
    else nColors = DEFAULT_COLORS;

    glEnable( GL_LINE_SMOOTH );
    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
    glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );

    glEnable( GL_LINE_STIPPLE );

    for ( unsigned int i = 0; i < pattern.size(); ++i ) {
//------------------------------------------------------------------------------
//	Draw triangle sides
//------------------------------------------------------------------------------
	for ( unsigned int j = 0; j < pattern[ i ].size(); ++j ) {
	    // Adjust halfedge
	    // f->halfedge() corresponds to triangle()[ 0 ] - triangle()[ 1 ]
	    // f->halfedge()->vertex() corresponds to triangle()[ 1 ];
	    Halfedge_handle baseH = pattern[ i ][ j ]->halfedge();
	    hv[ 0 ] = baseH;
	    hv[ 1 ] = baseH->next();
	    hv[ 2 ] = baseH->prev();

	    for ( unsigned int k = 0; k < NUM_SIDES; ++k ) {

		// drawType	0 : do nothing
		//		1 : internal edge
		//		2 : boundary edge
		int drawType;
		// If the edge is on the boundary
		if ( hv[ k ]->path() != NO_INDEX ) {
		    drawType = 2;
		}
		// If the edge is internal
		else {
		    // draw only when the source vertex ID is less than that of
		    // the target vertex
		    if ( hv[ k ]->opposite()->vertex()->id() < hv[ k ]->vertex()->id() ) {
			drawType = 1;
		    }
		    else drawType = 0;
		}

		if ( isConvex( hv[ k ] ) ) {
		    // glDisable( GL_LINE_STIPPLE );
		    glLineStipple( 1, 0xFFFF );
		    // glColor3d( 1.0, 0.0, 0.0 );
		}
		else {
		    // glEnable( GL_LINE_STIPPLE );
		    glLineStipple( 2, 0xAAAA );
		    // glLineStipple( 2, 0xF0F0 );
		    // glColor3d( 0.0, 1.0, 0.0 );
		}

		// If the edge should be drawn, ...
		switch ( drawType ) {
		  case 0:	// do nothing
		  case 2:	// boundary edge
		      break;
		  case 1:	// internal edge
		      glLineWidth( 2.0 );
		      glColor3d( 0.0, 0.0, 0.0 );
		      glBegin( GL_LINES );
		      glVertex2d( pattern[ i ][ j ]->triangle()[ k%NUM_SIDES ].x(), 
				  pattern[ i ][ j ]->triangle()[ k%NUM_SIDES ].y() );
		      glVertex2d( pattern[ i ][ j ]->triangle()[ (k+1)%NUM_SIDES ].x(), 
				  pattern[ i ][ j ]->triangle()[ (k+1)%NUM_SIDES ].y() );
		      glEnd();
		      break;
		}
	    }
	}

//------------------------------------------------------------------------------
//	Draw the bounding box
//------------------------------------------------------------------------------
	glColor3d( 0.7, 0.7, 0.7 );
	glLineWidth( 1.0 );
	glBegin( GL_LINE_LOOP );
	glVertex2d( bound[ i ].xmin(), bound[ i ].ymin() );
	glVertex2d( bound[ i ].xmax(), bound[ i ].ymin() );
	glVertex2d( bound[ i ].xmax(), bound[ i ].ymax() );
	glVertex2d( bound[ i ].xmin(), bound[ i ].ymax() );
	glEnd();
    }

    glDisable( GL_LINE_STIPPLE );
    glDisable( GL_LINE_SMOOTH );
    glDisable( GL_BLEND );
}

