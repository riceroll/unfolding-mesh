#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <algorithm>

#define NEW_GA

using namespace std;

#include "common.h"
#include "random.h"
#include "gene.h"

static vector< Halfedge_handle >		pEdges;
static Polyhedron				* pMesh;
static Bbox2					* pSheet;
static vector< vector < Facet_handle > >	* pPatch;
static vector< Bbox2 >				* pBound;
static Attribute				* pAttr;
static vector< vector < Facet_handle > >	sPatch,   fPatch;
static vector< vector < Triangle2 > >		sTriplet, fTriplet;
static vector< bool >				sConnect, fConnect;
static vector< int >				sNCuts,   fNCuts;
static Bbox2					sSheet,   fSheet;
static double					bestScore;

static void (*redraw)( void ) = NULL;

double	(*Genome::_fitness)( Genome & );

//------------------------------------------------------------------------------
//	
//------------------------------------------------------------------------------
// random number generator
// This is likely to genreate smaller numbers
int triRand( int max )
{
    int i, j, k, l, r_max = 0;

    for ( i = 0; i < max; ++i )
	r_max += i + 1;
    // i = lrand48() % r_max;
	i = rand() % r_max;
    for( j = max, k = max - 1, l = 0; j <= r_max; j += k, k--, l++ )
	if ( i < j ) return l;

    return 0;
}

//------------------------------------------------------------------------------
//	
//------------------------------------------------------------------------------
Genome & Genome::copy( const Genome & obj )
{
    if ( this != &obj ) {
	_dirty		= obj._dirty;
	_score		= obj._score;
	_chromosome	= obj._chromosome;
	_separator	= obj._separator;
    }
    return *this;
}

void Genome::init( void )
{
    int n = ( int )_chromosome.size();
    for ( int k = 0; k < n; ++k )
	_chromosome[ k ] = k;
    for ( int k = 0; k < n; ++k ) {
	// int p = triRand( n );
	// int q = triRand( n );
	int p = RandomInt( 0, n-1 );
	int q = RandomInt( 0, n-1 );
	SWAPINT( _chromosome[p], _chromosome[q] );
    }
    _separator = 0;
}

void Genome::eval( void )
{
    if ( _dirty ) {
	_score = (*_fitness)( *this );
	clearDirty();
    }
}


//------------------------------------------------------------------------------
//	
//------------------------------------------------------------------------------
int commonHead( Genome & p1, Genome & p2, 
		vector< unsigned int > & corres1, vector< unsigned int > & corres2 )
{
    
#ifdef MYDEBUG
    cerr << "[Before] p1 = " << p1;
    cerr << "[Before] p2 = " << p2;
#endif	// MYDEBUG

    // create: store
    vector< unsigned int > store( p1.chromosome().size() );
    
    // create: seq1, seq2, sorted1, sorted2
    vector< unsigned int > seq1, seq2, sorted1, sorted2;
    for ( unsigned int k = 0; k < p1.separator(); ++k ) {
	seq1.push_back( p1.chromosome()[ k ] );
    }
    sorted1 = seq1;
    sort( sorted1.begin(), sorted1.end() );
    
#ifdef MYDEBUG
    cerr << " sorted1 = ";
    for ( unsigned int k = 0; k < sorted1.size(); ++k )
	cerr << setw( 3 ) << sorted1[ k ];
    cerr << endl;
#endif	// MYDEBUG
	
    for ( unsigned int k = 0; k < p2.separator(); ++k ) {
	seq2.push_back( p2.chromosome()[ k ] );
    }
    sorted2 = seq2;
    sort( sorted2.begin(), sorted2.end() );
    
#ifdef MYDEBUG
    cerr << " sorted2 = ";
    for ( unsigned int k = 0; k < sorted2.size(); ++k )
	cerr << setw( 3 ) << sorted2[ k ];
    cerr << endl;
#endif	// MYDEBUG

    // create: pres_end
    vector< unsigned int >::iterator pres_end =
	set_intersection( sorted1.begin(), sorted1.end(),
			  sorted2.begin(), sorted2.end(),
			  store.begin() );
    // release: sorted1, sorted2
    sorted1.clear();
    sorted2.clear();

    // create: sect
    vector< unsigned int > sect;
    for ( vector< unsigned int >::iterator iter = store.begin(); iter != pres_end; ++iter ) {
	sect.push_back( *iter );
    }
    // release: store
    store.clear();

#ifdef MYDEBUG
    cerr << " sect = ";
    for ( unsigned int k = 0; k < sect.size(); ++k )
	cerr << setw( 3 ) << sect[ k ];
    cerr << endl;
#endif	// MYDEBUG
    unsigned int limit = sect.size();

    // create: idF1, idL1, posF1, posL1, reorder1
    vector< unsigned int > idF1, idL1, posF1, posL1, reorder1;
    for ( unsigned int k = 0; k < seq1.size(); ++k ) {
	vector< unsigned int >::iterator iter = find( sect.begin(), sect.end(), seq1[ k ] );
	if ( iter != sect.end() ) {
	    idF1.push_back ( seq1[ k ] );
	    posF1.push_back( k );
	}
	else {
	    idL1.push_back ( seq1[ k ] );
	    posL1.push_back( k );
	}
    }
    // initialize: corres1
    corres1.clear();
    for ( unsigned int k = 0; k < idF1.size(); ++k ) {
	reorder1.push_back( idF1 [ k ] );
	corres1.push_back ( posF1[ k ] );
    }
    for ( unsigned int k = 0; k < idL1.size(); ++k ) {
	reorder1.push_back( idL1 [ k ] );
	corres1.push_back ( posL1[ k ] );
    }
    // release: idF1, idL1, posF1, posL1
    idF1.clear();
    idL1.clear();
    posF1.clear();
    posL1.clear();

#ifdef MYDEBUG
    cerr << " reorder1 = ";
    for ( unsigned int k = 0; k < reorder1.size(); ++k )
	cerr << setw( 3 ) << reorder1[ k ];
    cerr << endl;
    cerr << " corres1  = ";
    for ( unsigned int k = 0; k < corres1.size(); ++k )
	cerr << setw( 3 ) << corres1[ k ];
    cerr << endl;
#endif	// MYDEBUG

    // create: idF2, idL2, posF2, posL2, reorder2
    vector< unsigned int > idF2, idL2, posF2, posL2, reorder2;
    for ( unsigned int k = 0; k < seq2.size(); ++k ) {
	vector< unsigned int >::iterator iter = find( sect.begin(), sect.end(), seq2[ k ] );
	if ( iter != sect.end() ) {
	    idF2.push_back ( seq2[ k ] );
	    posF2.push_back( k );
	}
	else {
	    idL2.push_back ( seq2[ k ] );
	    posL2.push_back( k );
	}
    }
    // initialize: corres2
    corres2.clear();
    for ( unsigned int k = 0; k < idF2.size(); ++k ) {
	reorder2.push_back( idF2 [ k ] );
	corres2.push_back ( posF2[ k ] );
    }
    for ( unsigned int k = 0; k < idL2.size(); ++k ) {
	reorder2.push_back( idL2 [ k ] );
	corres2.push_back ( posL2[ k ] );
    }
    // release: idF2, idL2, posF2, posL2
    idF2.clear();
    idL2.clear();
    posF2.clear();
    posL2.clear();

#ifdef MYDEBUG
    cerr << " reorder2 = ";
    for ( unsigned int k = 0; k < reorder2.size(); ++k )
	cerr << setw( 3 ) << reorder2[ k ];
    cerr << endl;
    cerr << " corres2  = ";
    for ( unsigned int k = 0; k < corres2.size(); ++k )
	cerr << setw( 3 ) << corres2[ k ];
    cerr << endl;
#endif	// MYDEBUG

    // release: sect
    sect.clear();

    // #ifdef MYDEBUG
    assert( reorder1.size() == p1.separator() );
    assert( reorder2.size() == p2.separator() );
    // #endif	// MYDEBUG

    for ( unsigned int k = 0; k < p1.separator(); ++k )
	p1.chromosome()[ k ] = reorder1[ k ];
    for ( unsigned int k = 0; k < p2.separator(); ++k )
	p2.chromosome()[ k ] = reorder2[ k ];

    // release: reorder1, reorder2
    reorder1.clear();
    reorder2.clear();
	
#ifdef MYDEBUG
    cerr << "[After]  p1 = " << p1;
    cerr << "[After]  p2 = " << p2;
    cerr << " limit = " << limit << endl;
#endif	// MYDEBUG

    return limit;
}


int commonTail( Genome & p1, Genome & p2, 
		vector< unsigned int > & corres1, vector< unsigned int > & corres2 )
{
    
#ifdef MYDEBUG
    cerr << "[Before] p1 = " << p1;
    cerr << "[Before] p2 = " << p2;
#endif	// MYDEBUG

    vector< unsigned int > store( p1.chromosome().size() );
    
    vector< unsigned int > seq1, seq2, sorted1, sorted2;
    for ( unsigned int k = p1.separator(); k < p1.chromosome().size(); ++k ) {
	seq1.push_back( p1.chromosome()[ k ] );
    }
    sorted1 = seq1;
    sort( sorted1.begin(), sorted1.end() );
    
#ifdef MYDEBUG
    cerr << " sorted1 = ";
    for ( unsigned int k = 0; k < sorted1.size(); ++k )
	cerr << setw( 3 ) << sorted1[ k ];
    cerr << endl;
#endif	// MYDEBUG
	
    for ( unsigned int k = p2.separator(); k < p2.chromosome().size(); ++k ) {
	seq2.push_back( p2.chromosome()[ k ] );
    }
    sorted2 = seq2;
    sort( sorted2.begin(), sorted2.end() );
    
#ifdef MYDEBUG
    cerr << " sorted2 = ";
    for ( unsigned int k = 0; k < sorted2.size(); ++k )
	cerr << setw( 3 ) << sorted2[ k ];
    cerr << endl;
#endif	// MYDEBUG

    vector< unsigned int >::iterator pres_end =
	set_intersection( sorted1.begin(), sorted1.end(),
			  sorted2.begin(), sorted2.end(),
			  store.begin() );
    sorted1.clear();
    sorted2.clear();

    vector< unsigned int > sect;
    for ( vector< unsigned int >::iterator iter = store.begin(); iter != pres_end; ++iter ) {
	sect.push_back( *iter );
    }
    store.clear();

#ifdef MYDEBUG
    cerr << " sect = ";
    for ( unsigned int k = 0; k < sect.size(); ++k )
	cerr << setw( 3 ) << sect[ k ];
    cerr << endl;
#endif	// MYDEBUG
    unsigned int limit = p1.chromosome().size() - sect.size();

    vector< unsigned int > idF1, idL1, posF1, posL1, reorder1;
    for ( unsigned int k = 0; k < seq1.size(); ++k ) {
	vector< unsigned int >::iterator iter = find( sect.begin(), sect.end(), seq1[ k ] );
	if ( iter != sect.end() ) {
	    idL1.push_back( seq1[ k ] );
	    posL1.push_back( k );
	}
	else {
	    idF1.push_back( seq1[ k ] );
	    posF1.push_back( k );
	}
    }
    corres1.clear();
    for ( unsigned int k = 0; k < idF1.size(); ++k ) {
	reorder1.push_back( idF1 [ k ] );
	corres1.push_back ( p1.separator()+posF1[ k ] );
    }
    for ( unsigned int k = 0; k < idL1.size(); ++k ) {
	reorder1.push_back( idL1 [ k ] );
	corres1.push_back ( p1.separator()+posL1[ k ] );
    }
    idF1.clear();
    idL1.clear();
    posF1.clear();
    posL1.clear();

#ifdef MYDEBUG
    cerr << " reorder1 = ";
    for ( unsigned int k = 0; k < reorder1.size(); ++k )
	cerr << setw( 3 ) << reorder1[ k ];
    cerr << endl;
    cerr << " corres1  = ";
    for ( unsigned int k = 0; k < corres1.size(); ++k )
	cerr << setw( 3 ) << corres1[ k ];
    cerr << endl;
#endif	// MYDEBUG

    vector< unsigned int > idF2, idL2, posF2, posL2, reorder2;
    for ( unsigned int k = 0; k < seq2.size(); ++k ) {
	vector< unsigned int >::iterator iter = find( sect.begin(), sect.end(), seq2[ k ] );
	if ( iter != sect.end() ) {
	    idL2.push_back( seq2[ k ] );
	    posL2.push_back( k );
	}
	else {
	    idF2.push_back( seq2[ k ] );
	    posF2.push_back( k );
	}
    }
    corres2.clear();
    for ( unsigned int k = 0; k < idF2.size(); ++k ) {
	reorder2.push_back( idF2 [ k ] );
	corres2.push_back ( p2.separator()+posF2[ k ] );
    }
    for ( unsigned int k = 0; k < idL2.size(); ++k ) {
	reorder2.push_back( idL2 [ k ] );
	corres2.push_back ( p2.separator()+posL2[ k ] );
    }
    idF2.clear();
    idL2.clear();
    posF2.clear();
    posL2.clear();

#ifdef MYDEBUG
    cerr << " reorder2 = ";
    for ( unsigned int k = 0; k < reorder2.size(); ++k )
	cerr << setw( 3 ) << reorder2[ k ];
    cerr << endl;
    cerr << " corres2  = ";
    for ( unsigned int k = 0; k < corres2.size(); ++k )
	cerr << setw( 3 ) << corres2[ k ];
    cerr << endl;
#endif	// MYDEBUG

    sect.clear();

    for ( unsigned int k = p1.separator(); k < p1.chromosome().size(); ++k )
	p1.chromosome()[ k ] = reorder1[ k-p1.separator() ];
    for ( unsigned int k = p2.separator(); k < p2.chromosome().size(); ++k )
	p2.chromosome()[ k ] = reorder2[ k-p2.separator() ];
	
    reorder1.clear();
    reorder2.clear();

#ifdef MYDEBUG
    cerr << "[After]  p1 = " << p1;
    cerr << "[After]  p2 = " << p2;
    cerr << " limit = " << limit << endl;
#endif	// MYDEBUG

    return limit;
}


//------------------------------------------------------------------------------
//	
//------------------------------------------------------------------------------
// Parital Match Crossover
void crossover( const Genome & p1, const Genome & p2, Genome * c1, Genome * c2 )
{
    const Genome &mom = DYN_CAST(const Genome &, p1);
    const Genome &dad = DYN_CAST(const Genome &, p2);

#ifdef NO_NEED
    if ( mom.chromosome().size() != dad.chromosome().size() ) {
	cerr << "Not equivalent lengths in crossover" << endl;
	assert( mom.chromosome().size() == dad.chromosome().size() );
    }
#endif	// NO_NEED

    int a = RandomInt( 0, mom.chromosome().size() );
    int b = RandomInt( 0, dad.chromosome().size() );
    if ( b < a ) SWAPINT( a, b );
    const unsigned int * index;
    int i, j;

    if ( c1 != NULL ) {
	Genome & sis = DYN_CAST( Genome &, *c1 );
	// sis.GAList<T>::copy(mom);
	sis.copy( mom );
	// GAListIter<T> diter(dad);
	// index = diter.warp(a);
	for ( i = a;  i < b; ++i ) {
	    index = &( dad.chromosome()[ i ] );
	    for ( j = 0; ( j < (int)sis.chromosome().size() ) && ( sis.chromosome()[ j ] != *index ); ++j );
	    unsigned int temp = sis.chromosome()[ i ];
	    sis.chromosome()[ i ] = sis.chromosome()[ j ];
	    sis.chromosome()[ j ] = temp;
	}
    }

    if ( c2 != NULL ) {
	Genome & bro = DYN_CAST( Genome &, *c2 );
	// bro.GAList<T>::copy(dad);
	bro.copy( dad );
	// GAListIter<T> miter(mom);
	// index = miter.warp(a);
	for ( i = a; i < b; i++ ) {
	    index = &( mom.chromosome()[ i ] );
	    for ( j = 0; ( j < (int)bro.chromosome().size() ) && ( bro.chromosome()[ j ] != *index ); ++j );
	    unsigned int temp = bro.chromosome()[ i ];
	    bro.chromosome()[ i ] = bro.chromosome()[ j ];
	    bro.chromosome()[ j ] = temp;
	}
    }

#ifdef DEBUG
    cerr << " a = " << a << " b = " << b << endl;
    cerr << " mom = " << mom;
    cerr << " dad = " << dad;
    cerr << " c1 = " << *c1;
    cerr << " c2 = " << *c2;
#endif	// DEBUG
}


//------------------------------------------------------------------------------
//	
//------------------------------------------------------------------------------
// Parital Match Crossover
void newcrossover( const Genome & p1, const Genome & p2, Genome * c1, Genome * c2 )
{
    const Genome &mom = DYN_CAST(const Genome &, p1);
    const Genome &dad = DYN_CAST(const Genome &, p2);

#ifdef NO_NEED
    if ( mom.chromosome().size() != dad.chromosome().size() ) {
	cerr << "Not equivalent lengths in crossover" << endl;
	assert( mom.chromosome().size() == dad.chromosome().size() );
    }
#endif	// NO_NEED

//------------------------------------------------------------------------------
//	compute the head commor part (where the order of IDs must be preserved)
//------------------------------------------------------------------------------
    Genome newMom, newDad;
    newMom.copy( mom );
    newDad.copy( dad );
    // #ifdef SKIP
    vector< unsigned int > corresMF, corresDF;
    int limitF = commonHead( newMom, newDad, corresMF, corresDF );
    assert( newMom.separator() == corresMF.size() );
    assert( newDad.separator() == corresDF.size() );
    // #endif	// SKIP

    // #ifdef SKIP
    vector< unsigned int > corresML, corresDL;
    int limitL = commonTail( newMom, newDad, corresML, corresDL );
    assert( newMom.chromosome().size() - newMom.separator() == corresML.size() );
    assert( newDad.chromosome().size() - newDad.separator() == corresDL.size() );
    // #endif	// SKIP

    int a = RandomInt( limitF, limitL );
    int b = RandomInt( limitF, limitL );
    if ( b < a ) SWAPINT( a, b );
    const unsigned int * index;
    int i, j;

    Genome sis, bro;
    if ( c1 != NULL ) {
	// Genome & sis = DYN_CAST( Genome &, *c1 );
	sis.copy( newMom );
	for ( i = a;  i < b; ++i ) {
	    index = &( newDad.chromosome()[ i ] );
	    for ( j = 0; ( j < (int)sis.chromosome().size() ) && ( sis.chromosome()[ j ] != *index ); ++j );
	    unsigned int temp = sis.chromosome()[ i ];
	    sis.chromosome()[ i ] = sis.chromosome()[ j ];
	    sis.chromosome()[ j ] = temp;
	}
    }

    if ( c2 != NULL ) {
	// Genome & bro = DYN_CAST( Genome &, *c2 );
	bro.copy( newDad );
	for ( i = a; i < b; i++ ) {
	    index = &( newMom.chromosome()[ i ] );
	    for ( j = 0; ( j < (int)bro.chromosome().size() ) && ( bro.chromosome()[ j ] != *index ); ++j );
	    unsigned int temp = bro.chromosome()[ i ];
	    bro.chromosome()[ i ] = bro.chromosome()[ j ];
	    bro.chromosome()[ j ] = temp;
	}
    }

#ifdef MYDEBUG
    cerr << " newMom = " << newMom;
    cerr << " newDad = " << newDad;

    cerr << " a = " << a << " b = " << b << endl;

    cerr << " sis = " << sis;
    cerr << " bro = " << bro;
#endif	// MYDEBUG

    Genome & newSis = DYN_CAST( Genome &, *c1 );
    Genome & newBro = DYN_CAST( Genome &, *c2 );

    newSis.copy( sis );
    for ( unsigned int k = 0; k < corresMF.size(); ++k )
	newSis.chromosome()[ corresMF[ k ] ] = sis.chromosome()[ k ];
    for ( unsigned int k = 0; k < corresML.size(); ++k )
	newSis.chromosome()[ corresML[ k ] ] = sis.chromosome()[ sis.separator() + k ];
    newBro.copy( bro );
    for ( unsigned int k = 0; k < corresDF.size(); ++k )
	newBro.chromosome()[ corresDF[ k ] ] = bro.chromosome()[ k ];
    for ( unsigned int k = 0; k < corresDL.size(); ++k )
	newBro.chromosome()[ corresDL[ k ] ] = bro.chromosome()[ bro.separator() + k ];

#ifdef MYDEBUG
    cerr << " newSis = " << newSis;
    cerr << " newBro = " << newBro;
#endif	// MYDEBUG

}

//------------------------------------------------------------------------------
//	
//------------------------------------------------------------------------------
bool Genome::mutate( float pmut )
{
    register int n, i;

    // cerr << "Original = " << *this;

    if ( pmut <= 0.0 ) return false;

    n = ( int )_chromosome.size();
    float nMut = pmut * STA_CAST(float, n);
    nMut *= 0.5;		// swapping one node swaps another as well
    if ( nMut < 1.0 ) {		// we have to do a flip test for each node
	nMut = 0;
	for ( i = 0; i < n; ++i ) {
	    if ( FlipCoin( pmut ) ) {
		int j = RandomInt( 0, n - 1 );
		SWAPINT( _chromosome[ i ], _chromosome[ j ] );
		nMut += 1.0;
	    }
	}
    }
    else {				// only nuke the number of nodes we need to
	for ( i = 0; i < nMut; i++ ) {
	    int p = RandomInt( 0, n - 1 );
	    int q = RandomInt( 0, n - 1 );
	    SWAPINT( _chromosome[ p ], _chromosome[ q ] );
	}
    }

    // cerr << "Mutated = " << *this;
    if ( nMut > 1.0e-10 ) return true;
    else return false;
}


void Genome::newmutate( void )
{
    register int n, i, j;

    n = ( int )_chromosome.size();
    i = ( int )RandomInt( 0, _separator );
    if ( i < 0 ) i = 0;
    if ( i >= ( int )_separator ) i = ( int )_separator - 1;
    j = ( int )RandomInt( _separator, n );
    if ( j < ( int )_separator ) j = ( int )_separator;
    if ( j >= n ) j = n - 1;

    if ( ! ( ( 0 <= i ) && ( i < n ) ) ) {
	cerr << " i = " << i << " n = " << n << endl;
	assert( ( 0 <= i ) && ( i < n ) );
    }
    if ( ! ( ( 0 <= j ) && ( j < n ) ) ) {
	cerr << " j = " << j << " n = " << n << endl;
	assert( ( 0 <= j ) && ( j < n ) );
    }
    SWAPINT( _chromosome[ i ], _chromosome[ j ] );
}


ostream & operator << ( ostream & stream, const Genome & obj )
{
    stream << " Score = " << obj._score << " : ";
    
    if ( obj._dirty ) stream << "[ ON]";
    else stream << "[OFF]";
    
    stream << " Chromosome : ";
    for ( unsigned int k = 0; k < obj._chromosome.size(); ++k ) {
	stream << setw( 4 ) << obj._chromosome[ k ];
    }
    stream << endl;
    stream << " Separator = " << obj._separator << endl;
    
    return stream;
}


//------------------------------------------------------------------------------
//	Popoluation
//------------------------------------------------------------------------------
void Population::_update( void )
{
    double maxScore = _pool[ _pool.size()-1 ]->score();
    double minScore = _pool[ 0 ]->score();

    _psum[ 0 ] = -_pool[ 0 ]->score() + maxScore + minScore;

    for ( unsigned int k = 1; k < _pool.size(); ++k )
	_psum[ k ] = -_pool[ k ]->score() + maxScore + minScore + _psum[ k-1 ];
    for ( unsigned k = 0; k < _pool.size(); ++k )
	_psum[ k ] /= _psum[ _pool.size() - 1 ];
}


void Population::_prep( void )
{
    if ( _ready ) return;

    _update();
    _ready = true;
}

Genome & Population::_rouletteWheel( void )
{
    float cutoff;
    int i;
#ifdef NEW_GA
    static int oldindex = -1;
#endif	// NEW_GA

    cutoff = RandomFloat();
    int n = ( int )_pool.size();
    int lower = 0; 
    int upper = n - 1;
    while ( upper >= lower ) {
	i = lower + ( upper - lower )/2;
	assert( i >= 0 && i < n );
	if ( _psum[ i ] > cutoff )
	    upper = i - 1;
	else
	    lower = i + 1;
    }
    lower = MIN2( n - 1, lower );
    lower = MAX2( 0, lower );
    
#ifdef DEBUG
    for ( unsigned int k = 0; k < _pool.size(); ++k ) {
	cerr << " k = " << k << " Score = " << _pool[ k ]->score() << " psum = " << _psum[ k ] << endl;
    }
#endif	// DEBUG

#ifdef NEW_GA
    if ( lower == oldindex ) {
	// cerr << "The same index as the previous" << endl;
	lower = ( lower + 1 ) % n;
    }
    oldindex = lower;
#endif	// NEW_GA

#ifdef MYDEBUG
    cerr << " cutOff = " << cutoff;
    cerr << " index = " << lower << endl;
#endif	// MYDEBUG

    return (*_pool[ lower ]);
}


void Population::insert( Genome & genome )
{
    // cerr << "Evaluating the new genome" << endl;
    genome.eval();

    Genome * ptrGenome = new Genome();
    // ptrGenome->resize( _nChromosomes );
    ptrGenome->copy( genome );
    vector< Genome * >::iterator it = _pool.begin();
    // int count = 0;
    while ( it != _pool.end() ) {
	// cerr << " count = " << count++ << endl;
	if ( (*it)->score() > ptrGenome->score() ) {
	    _pool.insert( it, ptrGenome );
	    return;
	}
	it++;
    }
    _pool.insert( _pool.end(), ptrGenome );
}

void Population::sort( void )
{
    unsigned int size = _pool.size();
    Genome* temp;

    // sort the chromosomes according to their costs
    // Wow!! : This is not the best sorting algorithm, should be replace with
    // better one later
    for ( unsigned int i = 0; i < size - 1; ++i )
	for ( unsigned int j = i + 1; j < size; ++j )
	    if ( _pool[ i ]->score() > _pool[ j ]->score() ) {
		temp		= _pool[ i ];
		_pool[ i ]	= _pool[ j ];
		_pool[ j ]	= temp;
	    }
}


void Population::step( void )
{
    unsigned int k;
    bool mut;
    // int c1, c2;
    Genome *mom, *dad;          // tmp holders for selected genomes
    
    // Generate the individuals in the temporary population from individuals in 
    // the main population.
    
    int nChanges = ( int )floor( POPULATION_COUNT * pREPLACEMENT );
    vector< Genome * > tmpGenome( nChanges );
    for ( unsigned int i = 0; i < tmpGenome.size(); ++i ) {
	tmpGenome[ i ] = new Genome;
	tmpGenome[ i ]->resize( _nChromosomes );
    }
    // takes care of odd population
    for ( k = 0; k < tmpGenome.size() - 1; k += 2 ) {	
        mom = &(select());  
	dad = &(select());

	// c1 = c2 = 0;

	if ( FlipCoin( pCROSSOVER ) ) {
	    crossover( *mom, *dad, tmpGenome[ k ], tmpGenome[ k+1 ] );
	    // c1 = c2 = 1;
	}
	else {
	    tmpGenome[  k  ]->copy( *mom );
	    tmpGenome[ k+1 ]->copy( *dad );
	}
	mut = tmpGenome[  k  ]->mutate( pMUTATION );
        // if ( mut > 0 ) c1 = 1;
	mut = tmpGenome[ k+1 ]->mutate( pMUTATION );
	// if ( mut > 0 ) c2 = 1;
    }

    if ( tmpGenome.size() % 2 != 0 ) {	// do the remaining population member
	mom = &(select());  
	dad = &(select());

	// c1 = 0;
	if ( FlipCoin( pCROSSOVER ) ) {
	    crossover( *mom, *dad, tmpGenome[ k ], (Genome*)NULL );
	    // c1 = 1;
	}
	else {
	    if ( RandomBit() )
		tmpGenome[ k ]->copy( *mom );
	    else
		tmpGenome[ k ]->copy( *dad );
	}
	mut = tmpGenome[  k  ]->mutate( pMUTATION );
	// if ( mut > 0 ) c1 = 1;
    }

    // Replace the worst genomes in the main population with all of the individuals
    // we just created.  Notice that we invoke the population's add member with a
    // genome pointer rather than reference.  This way we don't force a clone of
    // the genome - we just let the population take over.  Then we take it back by
    // doing a remove then a replace in the tmp population.
    
    for ( unsigned int i = 0; i < tmpGenome.size(); i++ ) {
        tmpGenome[ i ]->setDirty();
	insert( (*tmpGenome[ i ]) );
    }    
	      
    // the individuals in tmpPop are all owned by pop, but tmpPop does not know 
    // that.  so we use replace to take the individuals from the pop and stick 
    // them back into tmpPop
    for ( unsigned int i = 0; i < tmpGenome.size(); i++ ) {
	delete tmpGenome[ i ];
	tmpGenome[ i ] = NULL;
    }
    tmpGenome.resize( 0 );

    unsigned int limit = _pool.size();
    // cerr << " limit = " << limit << endl;
    for ( unsigned int i = limit - 1; i >= POPULATION_COUNT; i-- ) {
	delete _pool[ i ];
	_pool[ i ] = NULL;
	_pool.pop_back();
    }

    // cerr << "Now the population = " << (*this);

    assert( _pool.size() == POPULATION_COUNT );
}


void Population::newstep( void )
{
    unsigned int k;
    // int c1, c2;
    Genome *mom, *dad;          // tmp holders for selected genomes
    
    // Generate the individuals in the temporary population from individuals in 
    // the main population.
    
    int nChanges = ( int )floor( POPULATION_COUNT * pREPLACEMENT );
    vector< Genome * > tmpGenome( nChanges );
    for ( unsigned int i = 0; i < tmpGenome.size(); ++i ) {
	tmpGenome[ i ] = new Genome;
	tmpGenome[ i ]->resize( _nChromosomes );
    }
    // takes care of odd population
    for ( k = 0; k < tmpGenome.size() - 1; k += 2 ) {	
        mom = &(select());  
	dad = &(select());

#ifdef MYDEBUG
	cerr << " Mom : " << *mom;
	cerr << " Dad : " << *dad;
#endif	// MYDEBUG
	// c1 = c2 = 0;

	if ( FlipCoin( pCROSSOVER ) ) {
	    newcrossover( *mom, *dad, tmpGenome[ k ], tmpGenome[ k+1 ] );
	    // c1 = c2 = 1;
	}
	else {
	    tmpGenome[  k  ]->copy( *mom );
	    tmpGenome[ k+1 ]->copy( *dad );
	}

	if ( FlipCoin( pMUTATION ) ) {
	    tmpGenome[  k  ]->newmutate();
	}
	if ( FlipCoin( pMUTATION ) ) {
	    tmpGenome[ k+1 ]->newmutate();
	}
    }

    if ( tmpGenome.size() % 2 != 0 ) {	// do the remaining population member
	mom = &(select());  
	dad = &(select());

	// c1 = 0;
	if ( FlipCoin( pCROSSOVER ) ) {
	    newcrossover( *mom, *dad, tmpGenome[ k ], (Genome*)NULL );
	    // c1 = 1;
	}
	else {
	    if ( RandomBit() )
		tmpGenome[ k ]->copy( *mom );
	    else
		tmpGenome[ k ]->copy( *dad );
	}

	if ( FlipCoin( pMUTATION ) ) {
	    tmpGenome[ k ]->newmutate();
	}
	// if ( mut > 0 ) c1 = 1;
    }

    // Replace the worst genomes in the main population with all of the individuals
    // we just created.  Notice that we invoke the population's add member with a
    // genome pointer rather than reference.  This way we don't force a clone of
    // the genome - we just let the population take over.  Then we take it back by
    // doing a remove then a replace in the tmp population.
    
    for ( unsigned int i = 0; i < tmpGenome.size(); i++ ) {
        tmpGenome[ i ]->setDirty();
	insert( (*tmpGenome[ i ]) );
    }    
	      
    // the individuals in tmpPop are all owned by pop, but tmpPop does not know 
    // that.  so we use replace to take the individuals from the pop and stick 
    // them back into tmpPop
    for ( unsigned int i = 0; i < tmpGenome.size(); i++ ) {
	delete tmpGenome[ i ];
	tmpGenome[ i ] = NULL;
    }
    tmpGenome.clear();

    unsigned int limit = _pool.size();
    // cerr << " limit = " << limit << endl;
    for ( unsigned int i = limit - 1; i >= POPULATION_COUNT; i-- ) {
	delete _pool[ i ];
	_pool[ i ] = NULL;
	_pool.pop_back();
    }

    // cerr << "Now the population = " << (*this);

    assert( _pool.size() == POPULATION_COUNT );
}

ostream & operator << ( ostream & stream, const Population & obj )
{
    for ( unsigned int k = 0; k < obj._pool.size(); ++k ) {
	stream << "[" << setw( 4 ) << k << "] : " << *(obj._pool[ k ]);
    }
    cerr << endl;
    return stream;
}





//------------------------------------------------------------------------------
//	
//	
//	
//	
//------------------------------------------------------------------------------
#ifdef OLD
// Calculate the cost of the input chromosome
double orderChanges( Genome & x )
{
    double cost = 0.0;
    for ( unsigned int i = 0; i < x.chromosome().size(); ++i ) {
	for ( unsigned int j = i+1; j < x.chromosome().size(); ++j ) {
	    if ( (x.chromosome()[ i ]) > (x.chromosome()[ j ]) ) cost += 1.0;
	}
    }
    return cost;
}
#endif	// OLD


double averageWeight( Polyhedron & poly )
{
    double sum = 0.0;
    int num = 0;
    const double penalty = 1.0;
    for ( Halfedge_iterator hi = poly.halfedges_begin(); hi != poly.halfedges_end(); ++hi ) {
	if ( hi->connect() ) sum += hi->weight();
	else sum += penalty;
	num++;
    }
    return sum/(double)num;
}


void applyGenome( Genome & genome )
{
    vector< Triangle2 > small;
    vector< int > accept, reject;
    vector< unsigned int > replacement;
    
    for ( unsigned int i = 0; i < (unsigned int)genome.chromosome().size(); ++i ) {
	int id = genome.chromosome()[ i ];
	Halfedge_handle hh = pEdges[ id ];
	if ( canEmbed( hh, *(pPatch), small ) ) {
	    cerr << ".";
	    stitch( hh, *(pPatch), small );
	    resetCycleIDs( *(pPatch) );
	    accept.push_back( id );
	}
	else {
	    reject.push_back( id );
	}

	if ( pPatch->size() <= 1 ) {

	    //------------------------------------------------------------------------------
	    //	Actually calculate the layout over the sheet
	    //------------------------------------------------------------------------------
	    // pBound->clear();	// This is very important
	    // This is necessary for evaluating the paper size.
	    arrangeSingle( *(pPatch), *(pSheet), *(pBound) );
	    traverseBoundary( *(pMesh), *(pAttr) );

	    for ( unsigned int k = 0; k < accept.size(); ++k )
		replacement.push_back( accept[ k ] );
	    for ( unsigned int k = 0; k < reject.size(); ++k )
		replacement.push_back( reject[ k ] );
	    for ( unsigned int k = i+1; k < ( unsigned int )genome.chromosome().size(); ++k )
		replacement.push_back( genome.chromosome()[ k ] );
	    assert( genome.chromosome().size() == replacement.size() );
	    genome.chromosome() = replacement;
	    genome.setSeparator( accept.size() );

	    return;
	}
    }

    //------------------------------------------------------------------------------
    //	Actually calculate the layout over the sheet
    //------------------------------------------------------------------------------
    // pBound->clear();	// This is very important
    // This is necessary for evaluating the paper size.
    arrangeSingle( *(pPatch), *(pSheet), *(pBound) );
    traverseBoundary( *(pMesh), *(pAttr) );

    for ( unsigned int k = 0; k < accept.size(); ++k )
	replacement.push_back( accept[ k ] );
    for ( unsigned int k = 0; k < reject.size(); ++k )
	replacement.push_back( reject[ k ] );
    assert( genome.chromosome().size() == replacement.size() );
    genome.chromosome() = replacement;
    genome.setSeparator( accept.size() );

    return;
}


//------------------------------------------------------------------------------
//	Objective functions
//------------------------------------------------------------------------------
// Number of patches and its associated infomation
double cluster( void )
{
    //	Number of components
    double nPatches	= ( double )pPatch->size();
    // unsigned int totalFaces = pMesh->size_of_facets();
    unsigned int maxFaces = 0;
    for ( unsigned int k = 0; k < pPatch->size(); ++k ) {
	if ( maxFaces < (*pPatch)[ k ].size() )
	    maxFaces = (*pPatch)[ k ].size();
    }
    double config = 1.0 + nPatches;

    if ( config < 0.0 ) {
	cerr << " nPatches = " << nPatches << " config = " << config << endl;
    }

    return config;
}


// Number of patches and its associated infomation
double configuration( void )
{
    //	Number of components
    double nPatches	= ( double )pPatch->size();
    unsigned int totalFaces = pMesh->size_of_facets();
    unsigned int maxFaces = 0;
    for ( unsigned int k = 0; k < pPatch->size(); ++k ) {
	if ( maxFaces < (*pPatch)[ k ].size() )
	    maxFaces = (*pPatch)[ k ].size();
    }
    double ratio = ( double )maxFaces/( double )totalFaces;
    double config = 1.0 + nPatches - ratio;

    if ( config < 0.0 ) {
	cerr << " nPatches = " << nPatches << " ratio = " << ratio << " config = " << config << endl;
    }

    return config;
}

// Number of patches and its associated infomation
double normalizedArea( void )
{
    double area		=  
	( double )( (pSheet->xmax() - pSheet->xmin())*
		    (pSheet->ymax() - pSheet->ymin()) );
    double fragment = 1.0 - pAttr->totalArea()/area;

    if ( ( fragment < 0.0 ) || ( fragment > 1.0 ) ) {
	cerr << " area = " << area << " totalArea = " << pAttr->totalArea();
	cerr << " fragment = " << fragment << endl;
	assert( false );
    }

    return fragment;
}






//------------------------------------------------------------------------------
//	Post process
//------------------------------------------------------------------------------
void postprocess( Genome & genome, double score )
{
    // cerr << " postProcess : score " << score << " bestScore = " << bestScore << endl;
    if ( bestScore > score ) {
	bestScore = score;

	cerr << endl;
	cerr << " Genome sequence : ";
	for ( unsigned int m = 0; m < genome.chromosome().size(); ++ m )
	    cerr << setw( 4 ) << genome.chromosome()[ m ];
	cerr << endl;
	cerr << " objection function score = " << score << endl;

	// pBound->clear();	// This is very important.
	arrangeSingle( *(pPatch), *(pSheet), *(pBound) );
	traverseBoundary( *(pMesh), *(pAttr) );
	cerr << " nRuns = " << pAttr->nRuns()
	     << " nStars = " << pAttr->nStars()
	     << " nHyperbolics = " << pAttr->nHyperbolics()
	     << " maxLength = " << pAttr->maxLength()
	     << " sumWeights = " << pAttr->sumWeights() << endl;

	(*redraw)();

	fPatch		= *pPatch;
	fSheet		= *pSheet;
	fTriplet.resize( pPatch->size() );
	for ( unsigned int i = 0; i < pPatch->size(); ++i ) {
	    fTriplet[ i ].resize( (*pPatch)[ i ].size() );
	    for ( unsigned int j = 0; j < (*pPatch)[ i ].size(); ++j )
		fTriplet[ i ][ j ] = (*pPatch)[ i ][ j ]->triangle();
	}
	for ( Halfedge_iterator hi = pMesh->halfedges_begin(); hi != pMesh->halfedges_end(); ++hi )
	    fConnect[ hi->id() ]	= hi->connect();
	for ( Vertex_iterator vi = pMesh->vertices_begin(); vi != pMesh->vertices_end(); ++vi )
	    fNCuts[ vi->id() ]		= vi->nCuts();
    }

    // restore the original mesh status
    *(pPatch) = sPatch;
    for ( Halfedge_iterator hi = pMesh->halfedges_begin(); hi != pMesh->halfedges_end(); ++hi ) {
	hi->connect() = sConnect[ hi->id() ];
    }
    for ( Vertex_iterator vi = pMesh->vertices_begin(); vi != pMesh->vertices_end(); ++vi ) {
	vi->nCuts() = sNCuts[ vi->id() ];
    }
    resetCycleIDs( *pPatch );
    countNCuts( *pMesh );
    *(pSheet) = sSheet;
    // pBound->clear();	// This is very important
    arrangeSingle( *(pPatch), *(pSheet), *(pBound) );
}

double test_measure( Genome & genome )
{
//------------------------------------------------------------------------------
//	Preprocessing
//------------------------------------------------------------------------------
    cerr << endl << flush;
    cerr << "*";
    applyGenome( genome );

//------------------------------------------------------------------------------
//	Score estimation
//------------------------------------------------------------------------------
    //double config	= configuration();
    double config	= cluster();
    //double area		= normalizedArea();
    // Total dual distance
    //assert( pAttr->aveDist() <= 1.0 );
    //double aveDist = pAttr->aveDist();

//------------------------------------------------------------------------------
//	Total score
//------------------------------------------------------------------------------
    double score =
	WEIGHT_NUMBER * config;
    //+
    //WEIGHT_AREA * area +
    //WEIGHT_DUAL * aveDist;


//------------------------------------------------------------------------------
//	Postprocessing
//------------------------------------------------------------------------------
    postprocess( genome, score );

    return score;
}

double number_only( Genome & genome )
{
//------------------------------------------------------------------------------
//	Preprocessing
//------------------------------------------------------------------------------
    cerr << endl << flush;
    cerr << "*";
    applyGenome( genome );

//------------------------------------------------------------------------------
//	Score estimation
//------------------------------------------------------------------------------
    //	Number of components
    double config	= configuration();

//------------------------------------------------------------------------------
//	Total score
//------------------------------------------------------------------------------
    double score =
	WEIGHT_NUMBER * config;

//------------------------------------------------------------------------------
//	Postprocessing
//------------------------------------------------------------------------------
    postprocess( genome, score );
    
    return score;
}


double number_area( Genome & genome )
{
//------------------------------------------------------------------------------
//	Preprocessing
//------------------------------------------------------------------------------
    cerr << endl << flush;
    cerr << "*";
    applyGenome( genome );

//------------------------------------------------------------------------------
//	Score estimation
//------------------------------------------------------------------------------
    double config	= configuration();
    double area		= normalizedArea();

//------------------------------------------------------------------------------
//	Total score
//------------------------------------------------------------------------------
    double score =
	WEIGHT_NUMBER * config +
	WEIGHT_AREA * area;

//------------------------------------------------------------------------------
//	Postprocessing
//------------------------------------------------------------------------------
    postprocess( genome, score );

    return score;
}


double number_run( Genome & genome )
{
//------------------------------------------------------------------------------
//	Preprocessing
//------------------------------------------------------------------------------
    cerr << endl << flush;
    cerr << "*";
    applyGenome( genome );

//------------------------------------------------------------------------------
//	Score estimation
//------------------------------------------------------------------------------
    double config	= configuration();
    // Number of runs
    int nRuns		= pAttr->nRuns();

//------------------------------------------------------------------------------
//	Total score
//------------------------------------------------------------------------------
    double score =
	WEIGHT_NUMBER * config +
	WEIGHT_RUN * ( double )nRuns/( double )pMesh->size_of_halfedges();

//------------------------------------------------------------------------------
//	Postprocessing
//------------------------------------------------------------------------------
    postprocess( genome, score );

    return score;
}


double number_star( Genome & genome )
{
//------------------------------------------------------------------------------
//	Preprocessing
//------------------------------------------------------------------------------
    cerr << endl << flush;
    cerr << "*";
    applyGenome( genome );

//------------------------------------------------------------------------------
//	Score estimation
//------------------------------------------------------------------------------
    double config	= configuration();
    // Number of runs
    double rStars	= ( double )pAttr->nStars()/( double )pMesh->size_of_vertices();

//------------------------------------------------------------------------------
//	Total score
//------------------------------------------------------------------------------
    double score =
	WEIGHT_NUMBER * config +
	WEIGHT_STAR * ( double )rStars;

//------------------------------------------------------------------------------
//	Postprocessing
//------------------------------------------------------------------------------
    postprocess( genome, score );

    return score;
}


double number_weight( Genome & genome )
{
//------------------------------------------------------------------------------
//	Preprocessing
//------------------------------------------------------------------------------
    cerr << endl << flush;
    cerr << "*";
    applyGenome( genome );

//------------------------------------------------------------------------------
//	Score estimation
//------------------------------------------------------------------------------
    double config	= configuration();
    // Sum of weights
    int sumWeights	= pAttr->sumWeights();

//------------------------------------------------------------------------------
//	Total score
//------------------------------------------------------------------------------
    double score = 
	WEIGHT_NUMBER * config +
	(-0.05) * sumWeights;

//------------------------------------------------------------------------------
//	Postprocessing
//------------------------------------------------------------------------------
    postprocess( genome, score );

    return score;
}


double number_area_run( Genome & genome )
{
//------------------------------------------------------------------------------
//	Preprocessing
//------------------------------------------------------------------------------
    cerr << endl << flush;
    cerr << "*";
    applyGenome( genome );

//------------------------------------------------------------------------------
//	Score estimation
//------------------------------------------------------------------------------
    double config	= configuration();
    double area		= normalizedArea();
    // Number of runs
    // int nRuns		= pAttr->nRuns();
    // maximum length of boundary runs
    int maxLength	= pAttr->maxLength();

//------------------------------------------------------------------------------
//	Total score
//------------------------------------------------------------------------------
    double score =
	WEIGHT_NUMBER * config +
	WEIGHT_AREA * area +
	WEIGHT_RUN * ( double )maxLength/( double )pMesh->size_of_halfedges();
#ifdef SKIP
    double score =
	WEIGHT_AREA * config +
	10.0 * WEIGHT_NUMBER * area +
	WEIGHT_RUN * ( double )maxLength/( double )pMesh->size_of_halfedges();
#endif	// SKIP


//------------------------------------------------------------------------------
//	Postprocessing
//------------------------------------------------------------------------------
    postprocess( genome, score );

    return score;
}


double number_area_star( Genome & genome )
{
//------------------------------------------------------------------------------
//	Preprocessing
//------------------------------------------------------------------------------
    cerr << endl << flush;
    cerr << "*";
    applyGenome( genome );

//------------------------------------------------------------------------------
//	Score estimation
//------------------------------------------------------------------------------
    double config	= configuration();
    double area		= normalizedArea();
    // Number of stars
    double rStars	= ( double )pAttr->nStars()/( double )pMesh->size_of_vertices();

//------------------------------------------------------------------------------
//	Total score
//------------------------------------------------------------------------------
    cerr << " config = " << config << " area = " << area << " rStars = " << rStars << endl;
    double score =
	WEIGHT_NUMBER * config +
	WEIGHT_AREA * area +
	WEIGHT_STAR * rStars;

//------------------------------------------------------------------------------
//	Postprocessing
//------------------------------------------------------------------------------
    postprocess( genome, score );

    return score;
}


double number_area_bias( Genome & genome )
{
//------------------------------------------------------------------------------
//	Preprocessing
//------------------------------------------------------------------------------
    cerr << endl << flush;
    cerr << "*";
    applyGenome( genome );

//------------------------------------------------------------------------------
//	Score estimation
//------------------------------------------------------------------------------
    double config	= configuration();
    double area		= normalizedArea();
    // Number of stars
    double rHyperbolics	= ( double )pAttr->nHyperbolics()/( double )pMesh->size_of_vertices();
    double rElliptics	= ( double )(pAttr->nStars() - pAttr->nHyperbolics())/( double )pMesh->size_of_vertices();

//------------------------------------------------------------------------------
//	Total score
//------------------------------------------------------------------------------
#ifdef MYDEBUG
    cerr << " config = " << config << " area = " << area << endl;
    cerr << " rElliptics = " << rElliptics << " rHyperbolics = " << rHyperbolics << endl;
#endif	// MYDEBUG
    double score =
	WEIGHT_NUMBER * config +
	WEIGHT_AREA * area +
	WEIGHT_STAR * ( double )( rElliptics + HYPERBOLIC_PENALTY * rHyperbolics );


//------------------------------------------------------------------------------
//	Postprocessing
//------------------------------------------------------------------------------
    postprocess( genome, score );

    return score;
}


double number_area_length( Genome & genome )
{
//------------------------------------------------------------------------------
//	Preprocessing
//------------------------------------------------------------------------------
    cerr << endl << flush;
    cerr << "*";
    applyGenome( genome );

//------------------------------------------------------------------------------
//	Score estimation
//------------------------------------------------------------------------------
    double config	= configuration();
    double area		= normalizedArea();
    // Number of stars
    double rHyperbolics	= ( double )pAttr->nHyperbolics()/( double )pMesh->size_of_vertices();
    double rElliptics	= ( double )(pAttr->nStars() - pAttr->nHyperbolics())/( double )pMesh->size_of_vertices();
    // maximum length of boundary runs
    int maxLength	= pAttr->maxLength();

//------------------------------------------------------------------------------
//	Total score
//------------------------------------------------------------------------------
#ifdef MYDEBUG
    cerr << " config = " << config << " area = " << area << endl;
    cerr << " rElliptics = " << rElliptics << " rHyperbolics = " << rHyperbolics << " maxLength = " << maxLength << endl;
#endif	// MYDEBUG
    double score =
	WEIGHT_NUMBER * config +
	WEIGHT_AREA * area +
	WEIGHT_STAR * ( double )( rElliptics + HYPERBOLIC_PENALTY * rHyperbolics ) +
	WEIGHT_RUN * ( double )maxLength/( double )pMesh->size_of_halfedges();

//------------------------------------------------------------------------------
//	Postprocessing
//------------------------------------------------------------------------------
    postprocess( genome, score );

    return score;
}


double number_max( Genome & genome )
{
//------------------------------------------------------------------------------
//	Preprocessing
//------------------------------------------------------------------------------
    cerr << endl << flush;
    cerr << "*";
    applyGenome( genome );

//------------------------------------------------------------------------------
//	Score estimation
//------------------------------------------------------------------------------
    double config	= configuration();
    // Max dual distance
    double maxDist	= pAttr->maxDist();

//------------------------------------------------------------------------------
//	Total score
//------------------------------------------------------------------------------
    double score =
	WEIGHT_NUMBER * config +
	WEIGHT_DUAL * maxDist;

//------------------------------------------------------------------------------
//	Postprocessing
//------------------------------------------------------------------------------
    postprocess( genome, score );

    return score;
}


double number_area_max( Genome & genome )
{
//------------------------------------------------------------------------------
//	Preprocessing
//------------------------------------------------------------------------------
    cerr << endl;
    cerr << "*";
    applyGenome( genome );

//------------------------------------------------------------------------------
//	Score estimation
//------------------------------------------------------------------------------
    double config	= configuration();
    double area		= normalizedArea();
    // Max dual distance
    double maxDist	= pAttr->maxDist();

//------------------------------------------------------------------------------
//	Total score
//------------------------------------------------------------------------------
    double score =
	WEIGHT_NUMBER * config +
	WEIGHT_AREA * area +
	WEIGHT_DUAL * maxDist;


//------------------------------------------------------------------------------
//	Postprocessing
//------------------------------------------------------------------------------
    postprocess( genome, score );

    return score;
}



double number_ave( Genome & genome )
{
//------------------------------------------------------------------------------
//	Preprocessing
//------------------------------------------------------------------------------
    cerr << endl << flush;
    cerr << "*";
    applyGenome( genome );

//------------------------------------------------------------------------------
//	Score estimation
//------------------------------------------------------------------------------
    double config	= configuration();
    // Total dual distance
    assert( pAttr->aveDist() <= 1.0 );
    double aveDist = pAttr->aveDist();

//------------------------------------------------------------------------------
//	Total score
//------------------------------------------------------------------------------
    double score =
	WEIGHT_NUMBER * config +
	WEIGHT_DUAL * aveDist;

//------------------------------------------------------------------------------
//	Postprocessing
//------------------------------------------------------------------------------
    postprocess( genome, score );

    return score;
}


double number_area_ave( Genome & genome )
{
//------------------------------------------------------------------------------
//	Preprocessing
//------------------------------------------------------------------------------
    cerr << endl << flush;
    cerr << "*";
    applyGenome( genome );

//------------------------------------------------------------------------------
//	Score estimation
//------------------------------------------------------------------------------
    double config	= configuration();
    double area		= normalizedArea();
    // Total dual distance
    assert( pAttr->aveDist() <= 1.0 );
    double aveDist = pAttr->aveDist();

//------------------------------------------------------------------------------
//	Total score
//------------------------------------------------------------------------------
    double score =
	WEIGHT_NUMBER * config +
	WEIGHT_AREA * area +
	WEIGHT_DUAL * aveDist;


//------------------------------------------------------------------------------
//	Postprocessing
//------------------------------------------------------------------------------
    postprocess( genome, score );

    return score;
}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//	Optimization function
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int optimize( Polyhedron & poly, 
	      Bbox2 & sheet,
	      vector< vector< Facet_handle > > & patch,
	      vector< Bbox2 > & bound,
	      Attribute & attr,
	      void (*preview)(void) )
{
    if ( patch.size() <= 1 ) {
	cerr << "We cannot merge unfolded patterns any more." << endl;
	return true;
    }

    redraw = preview;

    countNCuts( poly );
//------------------------------------------------------------------------------
//	Store the current status
//------------------------------------------------------------------------------
    pMesh		= &poly;
    pSheet		= &sheet;
    pPatch		= &patch;
    pBound		= &bound;
    pAttr		= &attr;

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
    bestScore = WEIGHT_NUMBER * ( double )poly.size_of_facets();


//------------------------------------------------------------------------------
//	Collect a set of stitchable edges
//------------------------------------------------------------------------------
    const int D = 3;
    map< int, Halfedge_handle > emap;
    map< int, Halfedge_handle >::iterator it;
    Halfedge_handle hh[ 3 ];

    for ( unsigned int i = 0; i < patch.size(); ++i ) {
	for ( unsigned int j = 0; j < patch[ i ].size(); ++j ) {
	    hh[ 0 ] = patch[ i ][ j ]->halfedge()->prev();
	    hh[ 1 ] = patch[ i ][ j ]->halfedge();
	    hh[ 2 ] = patch[ i ][ j ]->halfedge()->next();
	    for ( int k = 0; k < D; ++k ) {
		if ( isLocally( hh[ k ] ) ) {
		    int idCI = hh[ k ]->facet()->piece();
		    int idCY = hh[ k ]->opposite()->facet()->piece();
		    assert( idCI == ( int )i );
		    if ( patch[ i ].size() <= patch[ idCY ].size() ) {
			it = emap.find( hh[ k ]->id() );
			if ( it == emap.end() )
			    emap.insert( pair< int, Halfedge_handle >( hh[ k ]->id(), hh[ k ] ) );
		    }
		}
	    }
	}
    }
    
    unsigned int nStitches = 0;

    pEdges.clear();
    cerr << "Number of stichable edges = " << emap.size() << endl;
    it = emap.begin();
    while ( it != emap.end() ) {
#ifdef DEBUG
	cerr << "Edge ID [ " << setw( 3 ) << it->first 
	     << " ] : " << setw( 3 ) << it->second->opposite()->vertex()->id()
	     << " -- " << setw( 3 ) << it->second->vertex()->id() << endl;
#endif	// DEBUG
	pEdges.push_back( it->second );
	it++;
    }
    nStitches = pEdges.size();
    cerr << "Size of stichable edge list = " << pEdges.size() << endl;

//------------------------------------------------------------------------------
//	Begin the optimization using the Genetic Algorithm
//------------------------------------------------------------------------------
    Population ga;

    // cerr << "nStitches = " << nStitches << endl;
    // Genome::setFitness( objective );
    ga.init( nStitches );
    // cerr << "Population : " << endl << ga << endl;

    int nRepeats = 0;
    float curScore = ( float )poly.size_of_facets(), prevScore;
    int count = 0;
    while ( true ) {
	// while ( true ) {
	prevScore = curScore;
	// cerr << " generation = " << ga.generation() << endl;
	// cout << "Population = " << endl;
	// cout << ga.population() << endl;
#ifdef NEW_GA
	ga.newstep();
#else	// NEW_GA
	ga.step();
#endif	// NEW_GA
	
	curScore = ga.bestIndividual().score();
	cerr << endl;
	cerr << "####################" << endl;
	cerr << setw( 5 ) << ++count << " " << curScore << endl;
	cerr.flush();

	// if ( fabs( prevScore - curScore ) < 0.005 ) nRepeats++;
	if ( fabs( prevScore - curScore ) < 0.0001 ) nRepeats++;
	else nRepeats = 0;
	cerr << " nRepeats = " << nRepeats << endl;

	// getchar();
	if ( nRepeats >= 5 ) break;
	// if ( nRepeats >= 2 ) break;
    }

    Genome bestGenome = ga.bestIndividual();
    cerr << "the minimal number of unfolded patches is " << bestGenome.score() << "\n";
    cerr << "this can be obtained from the sequence of edge stitching\n";
    cerr << bestGenome << "\n\n";

#ifdef NO_NEED
    // cerr << "Before Genome : " << genome << endl;
    applyGenome( bestGenome );
#endif	// NO_NEED

    cerr << "[A] pPatch->size = " << pPatch->size() << endl;

    // restore the final mesh status
    *pPatch = fPatch;

    cerr << "[B] pPatch->size = " << pPatch->size() << endl;

    for ( Halfedge_iterator hi = pMesh->halfedges_begin(); hi != pMesh->halfedges_end(); ++hi ) {
	hi->connect() = fConnect[ hi->id() ];
    }
    for ( Vertex_iterator vi = pMesh->vertices_begin(); vi != pMesh->vertices_end(); ++vi ) {
	vi->nCuts() = fNCuts[ vi->id() ];
    }
    for ( unsigned int i = 0; i < pPatch->size(); ++i ) {
	for ( unsigned int j = 0; j < (*pPatch)[ i ].size(); ++j ) {
	    (*pPatch)[ i ][ j ]->triangle() = fTriplet[ i ][ j ];
	}
    }
    resetCycleIDs( *pPatch );
    countNCuts( *pMesh );

    bound.clear();		// This is very important
    arrangeSingle( *(pPatch), *(pSheet), *(pBound) );
    vector< vector< Halfedge_handle > > encode;
    traverseBoundary( *(pMesh), *(pAttr), encode );

    ga.clear();

//------------------------------------------------------------------------------
//	Processing boundary edge orders
//------------------------------------------------------------------------------
    for ( unsigned int k = 0; k < encode.size(); ++k ) {
	cerr << " Encode No. " << k << endl;
	for ( unsigned int m = 0; m < encode[ k ].size(); ++m ) {
	    cerr << "[" << encode[ k ][ m ]->path() << "]";
	    if ( encode[ k ][ m ]->orient() ) cerr << "^(+1)";
	    else cerr << "^(-1)";
	    cerr << "=";
	}
    }
    cerr << "##### Final bounding box ##### : "
	 << " xmin = " << pSheet->xmin()
	 << " ymin = " << pSheet->ymin()
	 << " xmax = " << pSheet->xmax()
	 << " ymax = " << pSheet->ymax() << endl;
#ifdef DEBUG
    cerr << "%%%%% paper size = "
	 << ( pSheet->xmax() - pSheet->xmin() ) * ( pSheet->ymax() - pSheet->ymin() )
	 << "%%%%%" << endl;
#endif	// DEBUG
    cerr << "%%%%% paper width = " << ( pSheet->xmax() - pSheet->xmin() ) << "%%%%%" << endl;
    cerr << " nRuns = " << pAttr->nRuns()
	 << " nStars = " << pAttr->nStars()
	 << " nHyperbolics = " << pAttr->nHyperbolics()
	 << " maxLength = " << pAttr->maxLength()
	 << " sumWeights = " << pAttr->sumWeights() << endl;
    cerr << " occupation ratio = " << pAttr->totalArea()/( (pSheet->xmax()-pSheet->xmin())*
							   (pSheet->ymax()-pSheet->ymin()) );
    cerr << " maxDist = " << pAttr->maxDist() 
	 << " aveDist = " << pAttr->aveDist() << endl;

    cerr << "\a\a\a";
    return true;
}


int exhaustive( Polyhedron & poly, 
		Bbox2 & sheet,
		vector< vector< Facet_handle > > & patch,
		vector< Bbox2 > & bound,
		Attribute & attr,
		void (*preview)(void) )
{
    if ( patch.size() <= 1 ) {
	cerr << "We cannot merge unfolded patterns any more." << endl;
	return true;
    }

    cerr << "Exhaustive search" << endl;

    redraw = preview;

    countNCuts( poly );
//------------------------------------------------------------------------------
//	Store the current status
//------------------------------------------------------------------------------
    pMesh		= &poly;
    pSheet		= &sheet;
    pPatch		= &patch;
    pBound		= &bound;
    pAttr		= &attr;

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
    bestScore = ( double )poly.size_of_facets();


//------------------------------------------------------------------------------
//	Collect a set of stitchable edges
//------------------------------------------------------------------------------
    const int D = 3;
    map< int, Halfedge_handle > emap;
    map< int, Halfedge_handle >::iterator it;
    Halfedge_handle hh[ 3 ];

    for ( unsigned int i = 0; i < patch.size(); ++i ) {
	for ( unsigned int j = 0; j < patch[ i ].size(); ++j ) {
	    hh[ 0 ] = patch[ i ][ j ]->halfedge()->prev();
	    hh[ 1 ] = patch[ i ][ j ]->halfedge();
	    hh[ 2 ] = patch[ i ][ j ]->halfedge()->next();
	    for ( int k = 0; k < D; ++k ) {
		if ( isLocally( hh[ k ] ) ) {
		    int idCI = hh[ k ]->facet()->piece();
		    int idCY = hh[ k ]->opposite()->facet()->piece();
		    assert( idCI == ( int )i );
		    if ( patch[ i ].size() <= patch[ idCY ].size() ) {
			it = emap.find( hh[ k ]->id() );
			if ( it == emap.end() )
			    emap.insert( pair< int, Halfedge_handle >( hh[ k ]->id(), hh[ k ] ) );
		    }
		}
	    }
	}
    }
    
    unsigned int nStitches = 0;

    pEdges.clear();
    cerr << "Number of stichable edges = " << emap.size() << endl;
    it = emap.begin();
    while ( it != emap.end() ) {
#ifdef DEBUG
	cerr << "Edge ID [ " << setw( 3 ) << it->first 
	     << " ] : " << setw( 3 ) << it->second->opposite()->vertex()->id()
	     << " -- " << setw( 3 ) << it->second->vertex()->id() << endl;
#endif	// DEBUG
	pEdges.push_back( it->second );
	it++;
    }
    nStitches = pEdges.size();
    cerr << "Size of stichable edge list = " << pEdges.size() << endl;

//------------------------------------------------------------------------------
//	Begin the optimization using the Genetic Algorithm
//------------------------------------------------------------------------------
    Genome genome;
    genome.resize( nStitches );
    // Genome::setFitness( objective );
    // reset
    vector< unsigned int > seq( nStitches );
    for ( unsigned int k = 0; k < nStitches; ++k )
	seq[ k ] = k;

    int count = 0;
    do {
	for ( unsigned int k = 0; k < nStitches; ++k )
	    genome.chromosome()[ nStitches - k - 1 ] = seq[ k ];

	cerr << endl << " count = " << count;
	cerr << " Genome sequence : ";
	for ( unsigned int m = 0; m < genome.chromosome().size(); ++m )
	    cerr << setw( 4 ) << genome.chromosome()[ m ];
	cerr << endl;

	genome.setDirty();
	genome.eval();
	count++;
	
    } while ( next_permutation( seq.begin(), seq.end() ) );

    // restore the final mesh status
    *pPatch = fPatch;
    for ( Halfedge_iterator hi = pMesh->halfedges_begin(); hi != pMesh->halfedges_end(); ++hi ) {
	hi->connect() = fConnect[ hi->id() ];
    }
    for ( Vertex_iterator vi = pMesh->vertices_begin(); vi != pMesh->vertices_end(); ++vi ) {
	vi->nCuts() = fNCuts[ vi->id() ];
    }
    for ( unsigned int i = 0; i < pPatch->size(); ++i ) {
	for ( unsigned int j = 0; j < (*pPatch)[ i ].size(); ++j ) {
	    (*pPatch)[ i ][ j ]->triangle() = fTriplet[ i ][ j ];
	}
    }
    resetCycleIDs( *pPatch );
    countNCuts( *pMesh );

    // pBound->clear();		// This is very important
    arrangeSingle( *(pPatch), *(pSheet), *(pBound) );

#ifdef DEBUG
    cerr << "##### Final bounding box ##### : "
	 << " xmin = " << pSheet->xmin()
	 << " ymin = " << pSheet->ymin()
	 << " xmax = " << pSheet->xmax()
	 << " ymax = " << pSheet->ymax() << endl;
    cerr << "%%%%% paper size = "
	 << ( pSheet->xmax() - pSheet->xmin() ) * ( pSheet->ymax() - pSheet->ymin() )
	 << "%%%%%" << endl;
    cerr << "%%%%% paper width = " << ( pSheet->xmax() - pSheet->xmin() ) << "%%%%%" << endl;
#endif	// DEBUG

    return true;
}



