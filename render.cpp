//------------------------------------------------------------------------------
//
//	render.cpp
//
//------------------------------------------------------------------------------
#include "common.h"
#include "ui.h"
#include <iostream>
#include <GLUT/glut.h>

using namespace std;

/* normal light model */
float normalLM_ambient[] =      {  0.2,  0.2,  0.2,  1.0 };

/* normal material and light */
float red_material_ambient[] =     {  0.3,  0.0,  0.0,  1.0 };
float red_material_diffuse[] =     {  0.5,  0.0,  0.0,  1.0 };
float red_material_emission[] =    {  0.2,  0.0,  0.0,  1.0 };
float red_material_specular[] =    {  0.9,  0.0,  0.0,  1.0 };
float red_material_shininess[] =   {  5.0 };

float orange_material_ambient[] =     {  0.3,  0.2,  0.0,  1.0 };
float orange_material_diffuse[] =     {  0.5,  0.3,  0.0,  1.0 };
float orange_material_emission[] =    {  0.2,  0.1,  0.0,  1.0 };
float orange_material_specular[] =    {  0.9,  0.6,  0.0,  1.0 };
float orange_material_shininess[] =   {  5.0 };

float green_material_ambient[] =     {  0.0,  0.3,  0.0,  1.0 };
float green_material_diffuse[] =     {  0.0,  0.5,  0.0,  1.0 };
float green_material_emission[] =    {  0.0,  0.2,  0.0,  1.0 };
float green_material_specular[] =    {  0.0,  0.9,  0.0,  1.0 };
float green_material_shininess[] =   {  5.0 };

float cyan_material_ambient[] =     {  0.0,  0.2,  0.3,  1.0 };
float cyan_material_diffuse[] =     {  0.0,  0.3,  0.5,  1.0 };
float cyan_material_emission[] =    {  0.0,  0.1,  0.2,  1.0 };
float cyan_material_specular[] =    {  0.0,  0.6,  0.9,  1.0 };
float cyan_material_shininess[] =   {  5.0 };

float blue_material_ambient[] =     {  0.0,  0.0,  0.3,  1.0 };
float blue_material_diffuse[] =     {  0.0,  0.0,  0.5,  1.0 };
float blue_material_emission[] =    {  0.0,  0.0,  0.2,  1.0 };
float blue_material_specular[] =    {  0.0,  0.0,  0.9,  1.0 };
float blue_material_shininess[] =   {  5.0 };

float magenta_material_ambient[] =     {  0.3,  0.0,  0.3,  1.0 };
float magenta_material_diffuse[] =     {  0.5,  0.0,  0.5,  1.0 };
float magenta_material_emission[] =    {  0.2,  0.0,  0.2,  1.0 };
float magenta_material_specular[] =    {  0.9,  0.0,  0.9,  1.0 };
float magenta_material_shininess[] =   {  5.0 };

float gray_material_ambient[] =     {  0.3,  0.3,  0.3,  1.0 };
float gray_material_diffuse[] =     {  0.5,  0.5,  0.5,  1.0 };
float gray_material_emission[] =    {  0.2,  0.2,  0.2,  1.0 };
float gray_material_specular[] =    {  0.9,  0.9,  0.9,  1.0 };
float gray_material_shininess[] =   {  5.0 };

float lightL_ambient[] =        {  0.2,  0.2,  0.2,  1.0 };
float lightL_diffuse[] =        {  0.6,  0.6,  0.6,  1.0 };
float lightL_specular[] =       {  0.6,  0.6,  0.6,  1.0 };
//float lightL_position[] =     { -1.0,  0.0,  0.05,  0.0 };
float lightL_position[] =       { -3.0, -2.0, -1.0,  0.0 };

float lightR_ambient[] =        {  0.2,  0.2,  0.2,  1.0 };
float lightR_diffuse[] =        {  0.6,  0.6,  0.6,  1.0 };
float lightR_specular[] =       {  0.6,  0.6,  0.6,  1.0 };
//float lightR_position[] =     {  0.0,  0.0,  0.05,  0.0 };
float lightR_position[] =       {  3.0,  2.0,  1.0,  0.0 };

void initLights( void )
{
    glLightfv( GL_LIGHT0, GL_AMBIENT,               lightL_ambient );
    glLightfv( GL_LIGHT0, GL_DIFFUSE,               lightL_diffuse );
    glLightfv( GL_LIGHT0, GL_SPECULAR,              lightL_specular );
    glLightfv( GL_LIGHT0, GL_POSITION,              lightL_position );
    
    glLightfv( GL_LIGHT1, GL_AMBIENT,               lightR_ambient );
    glLightfv( GL_LIGHT1, GL_DIFFUSE,               lightR_diffuse );
    glLightfv( GL_LIGHT1, GL_SPECULAR,              lightR_specular );
    glLightfv( GL_LIGHT1, GL_POSITION,              lightR_position );
    
    glEnable( GL_LIGHTING );
    glEnable ( GL_LIGHT0 );
    glEnable ( GL_LIGHT1 );
    glDisable( GL_LIGHT2 );
    glDisable( GL_LIGHT3 );
    
    glLightModelfv( GL_LIGHT_MODEL_AMBIENT, normalLM_ambient );
}


void redMaterial( void )
{
    glMaterialfv( GL_FRONT, GL_AMBIENT,             red_material_ambient );
    glMaterialfv( GL_FRONT, GL_DIFFUSE,             red_material_diffuse );
    glMaterialfv( GL_FRONT, GL_SPECULAR,            red_material_specular );
    glMaterialfv( GL_FRONT, GL_SHININESS,           red_material_shininess );
    glMaterialfv( GL_FRONT, GL_EMISSION,            red_material_emission );
}


void orangeMaterial( void )
{
    glMaterialfv( GL_FRONT, GL_AMBIENT,             orange_material_ambient );
    glMaterialfv( GL_FRONT, GL_DIFFUSE,             orange_material_diffuse );
    glMaterialfv( GL_FRONT, GL_SPECULAR,            orange_material_specular );
    glMaterialfv( GL_FRONT, GL_SHININESS,           orange_material_shininess );
    glMaterialfv( GL_FRONT, GL_EMISSION,            orange_material_emission );
}


void greenMaterial( void )
{
    glMaterialfv( GL_FRONT, GL_AMBIENT,             green_material_ambient );
    glMaterialfv( GL_FRONT, GL_DIFFUSE,             green_material_diffuse );
    glMaterialfv( GL_FRONT, GL_SPECULAR,            green_material_specular );
    glMaterialfv( GL_FRONT, GL_SHININESS,           green_material_shininess );
    glMaterialfv( GL_FRONT, GL_EMISSION,            green_material_emission );
}


void cyanMaterial( void )
{
    glMaterialfv( GL_FRONT, GL_AMBIENT,             cyan_material_ambient );
    glMaterialfv( GL_FRONT, GL_DIFFUSE,             cyan_material_diffuse );
    glMaterialfv( GL_FRONT, GL_SPECULAR,            cyan_material_specular );
    glMaterialfv( GL_FRONT, GL_SHININESS,           cyan_material_shininess );
    glMaterialfv( GL_FRONT, GL_EMISSION,            cyan_material_emission );
}


void blueMaterial( void )
{
    glMaterialfv( GL_FRONT, GL_AMBIENT,             blue_material_ambient );
    glMaterialfv( GL_FRONT, GL_DIFFUSE,             blue_material_diffuse );
    glMaterialfv( GL_FRONT, GL_SPECULAR,            blue_material_specular );
    glMaterialfv( GL_FRONT, GL_SHININESS,           blue_material_shininess );
    glMaterialfv( GL_FRONT, GL_EMISSION,            blue_material_emission );
}


void magentaMaterial( void )
{
    glMaterialfv( GL_FRONT, GL_AMBIENT,             magenta_material_ambient );
    glMaterialfv( GL_FRONT, GL_DIFFUSE,             magenta_material_diffuse );
    glMaterialfv( GL_FRONT, GL_SPECULAR,            magenta_material_specular );
    glMaterialfv( GL_FRONT, GL_SHININESS,           magenta_material_shininess );
    glMaterialfv( GL_FRONT, GL_EMISSION,            magenta_material_emission );
}


void grayMaterial( void )
{
    glMaterialfv( GL_FRONT, GL_AMBIENT,             gray_material_ambient );
    glMaterialfv( GL_FRONT, GL_DIFFUSE,             gray_material_diffuse );
    glMaterialfv( GL_FRONT, GL_SPECULAR,            gray_material_specular );
    glMaterialfv( GL_FRONT, GL_SHININESS,           gray_material_shininess );
    glMaterialfv( GL_FRONT, GL_EMISSION,            gray_material_emission );
}


static void changeHSVToColor( double hue, double saturation, double value, double rgb[ 3 ] )
{
    int r,g,b;
    int region;
    float fraction;
    int min, max, up, down;

    // check and revise the input data
    while ( ( hue > 360.0f ) || ( hue < 0.0f ) ) {
	if ( hue >= 360.0f ) hue -= 360.0f;
	else if ( hue < 0.0f ) hue += 360.0f;
    }
    if ( saturation>1.0f ) saturation = 1.0f;
    else if ( saturation<0.0f ) saturation=0.0f;
    if ( value>1.0f ) value = 1.0f;
    else if ( value < 0.0f ) value = 0.0f;

    max = (int)( value*255 );

    if ( saturation==0.0f ) {
	r = max;
	g = max;
	b = max;
    }
    else{
	region		= (int)(hue/60.0f);
	fraction	= hue/60.0f-region;
	min		= (int)( max * (1.0f-saturation) );
	up		= min + (int)( fraction*max*saturation );
	down		= max - (int)( fraction*max*saturation );

	switch (region) {
	  case 0:
	      r= max; g= up; b= min; break; // red -> yellow
	  case 1:
	      r=down; g= max; b= min; break; // yellow -> green
	  case 2:
	      r= min; g= max; b= up; break; // green -> cyan
	  case 3:
	      r= min; g=down; b= max; break; // cyan -> blue
	  case 4:
	      r= up; g= min; b= max; break; // blue -> magenta
	  case 5:
	      r= max; g= min; b=down; break; // magenta -> red
	}
    }

    rgb[ 0 ] = ( double )r/255.0;
    rgb[ 1 ] = ( double )g/255.0;
    rgb[ 2 ] = ( double )b/255.0;

    return;
}


void changeNumberToColor( int step, int number, double rgb[ 3 ] )
{
    if ( number < 0 ) {
	cerr << "number = " << number << endl;
	assert( number >= 0 );
    }
    // if ( number < 0 ) number = 0;
    const int nOffsets = 7;
    double offset[ nOffsets ] = {
	0.50, 0.25, 0.75, 0.125, 0.375, 0.625, 0.875
	// 0.50, 0.75, 0.375, 0.125, 0.625, 0.250, 0.875
    };
    int quotient = number / step;
    int residue = number % step;
    int index = quotient % nOffsets;
    double hue = ( 360.0f / step ) * ( ( double )residue + offset[ index ] );
#ifdef DEBUG
    cerr << " quotient = " << quotient;
    cerr << " residue = " << residue;
    cerr << " offset = " << offset;
    cerr << " hue = " << hue;
    cerr << endl;
#endif	// DEBUG
    changeHSVToColor ( hue, 1.0f, 1.0f, rgb );
}
