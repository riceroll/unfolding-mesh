//------------------------------------------------------------------------------
//	
//	gene.h
//	
//------------------------------------------------------------------------------
#ifndef __GENE_H_
#define __GENE_H_

#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
using namespace std;

#define HYPERBOLIC_PENALTY	(10.0)

#define WEIGHT_NUMBER		(10.0)

#define WEIGHT_AREA		(1.0)
//#define WEIGHT_AREA		(0.1)

//#define WEIGHT_RUN		(0.5)
#define WEIGHT_RUN		(1.0)

//#define WEIGHT_STAR		(1.0)
#define WEIGHT_STAR		(0.1)

//#define WEIGHT_DUAL		(1.0)
#define WEIGHT_DUAL		(0.1)

// Number of individuals
//#define POPULATION_COUNT 5000
// #define POPULATION_COUNT	(16)
// #define POPULATION_COUNT	(32)
#define POPULATION_COUNT	(64) //<-
// #define POPULATION_COUNT	(128)
// #define POPULATION_COUNT	(256)
// Probability for mutation
// #define pMUTATION		(0.05)
#define pMUTATION		(0.1) //<-
// #define pMUTATION		(0.2)
// #define pMUTATION		(0.5)
// #define pMUTATION		(0.8)
// Probability for crossover
// #define pCROSSOVER		(0.95)
#define pCROSSOVER		(0.9) //<-
// #define pCROSSOVER		(0.8)
// #define pCROSSOVER		(0.3)
// #define pCROSSOVER		(0.0)
// Probability of replacement crossover
#define pREPLACEMENT		(0.5)
// Penalty value
#define PENALTY		(1.0e+10)


#define SWAPINT(a,b)	{int tmp=a; a=b; b=tmp;}
#define SWAPUINT(a,b)	{unsigned int tmp=a; a=b; b=tmp;}


class Genome {

private:
    bool			_dirty;
    double			_score;		// fitting function value
    vector< unsigned int >	_chromosome;	// flag for cutting constraints
    unsigned int		_separator;
    static double		(*_fitness)( Genome & );

protected:

public:
    Genome() {
	_dirty		= true;
	_score		= 1.0e+10;
	_chromosome.clear();
	_separator	= 0;
    }
    Genome( const Genome & obj ) {
	_dirty		= obj._dirty;
	_score		= obj._score;
	_chromosome	= obj._chromosome;
	_separator	= obj._separator;
    }
    ~Genome() {
	_chromosome.clear();
    }

    static void setFitness( double (*__fitness)( Genome & ) ) { Genome::_fitness = __fitness; }

    Genome & resize( unsigned int count ) {
	_chromosome.resize( count );
	_dirty		= true;
	_score		= 0.0;
	_separator	= 0;
	return *this;
    }
    void clear( void ) {
	_chromosome.clear();
	_dirty		= true;
	_score		= 0.0;
	_separator	= 0;
	return;
    }

    void setDirty( void )		{ _dirty = true; }
    void clearDirty( void )		{ _dirty = false; }
    const bool & dirty( void ) const	{ return _dirty; }

    const vector< unsigned int > &	chromosome( void ) const {
	return _chromosome;
    }
    vector< unsigned int > &		chromosome( void ) {
	return _chromosome;
    }

    const double & score( void ) const {
	return _score;
    }

    const unsigned int separator( void ) const {
	return _separator;
    }
    void setSeparator( unsigned int __separator ) {
	_separator = __separator;
    }

    Genome & copy( const Genome & obj );
    
    void init( void );
    void eval( void );

    friend int commonHead( Genome & p1, Genome & p2, 
			   vector< unsigned int > & corres1, vector< unsigned int > & corres2 );
    friend int commonTail( Genome & p1, Genome & p2, 
			   vector< unsigned int > & corres1, vector< unsigned int > & corres2 );


    friend void crossover( const Genome & p1, const Genome & p2, Genome * c1, Genome * c2 );
    friend void newcrossover( const Genome & p1, const Genome & p2, Genome * c1, Genome * c2 );

    bool mutate( float ratio );
    void newmutate( void );

    Genome & operator = ( const Genome & obj ) {
	if ( this != &obj ) {
	    _dirty	= obj._dirty;
	    _score	= obj._score;
	    _chromosome	= obj._chromosome;
	    _separator	= obj._separator;
	}
	return *this;
    }

    friend ostream & operator << ( ostream & stream, const Genome & obj );
};


class Population {

private:
    unsigned int	_nChromosomes;
    vector< Genome * >	_pool;
    bool		_ready;
    vector< float >	_psum;

protected:
    void	_update( void );
    void	_prep( void );
    Genome &	_rouletteWheel( void );

public:
    Population() {
	_pool.resize( POPULATION_COUNT );
	_psum.resize( POPULATION_COUNT );
	_ready = false;
    }

    Population( const Population & obj ) {
	cerr << "Illegal call" << endl;
	assert( false );
	unsigned int n = MIN2( _pool.size(), obj._pool.size() );
	for ( unsigned int k = 0; k < n; ++k )
	    _pool[ k ] = obj._pool[ k ];
    }

    Population & init( int __nChromosomes ) {
	_nChromosomes = __nChromosomes;
	for ( unsigned int k = 0; k < _pool.size(); ++k ) {
	    _pool[ k ] = new Genome;
	    _pool[ k ]->resize( _nChromosomes );
	    _pool[ k ]->init();
	    _pool[ k ]->eval();
	}
	sort();
	return *this;
    }

    Population & clear( void ) {
	_nChromosomes = 0;
	for ( unsigned int k = 0; k < _pool.size(); ++k ) {
	    _pool[ k ]->clear();
	    delete _pool[ k ];
	}
	return *this;
    }

    const unsigned int & nChromosomes( void ) const {
	return _nChromosomes;
    }

    const Genome & bestIndividual( void ) const {
	return (*_pool[ 0 ]);
    }

    void insert( Genome & obj );

    Genome & select( void ) {
	if ( ! _ready ) _prep();
	return _rouletteWheel();
    }

    void sort( void );

    void step( void );
    void newstep( void );

    friend ostream & operator << ( ostream & stream, const Population & obj );
};



extern double	number_only		( Genome & genome );
extern double	number_area		( Genome & genome );
extern double	number_run		( Genome & genome );
extern double	number_star		( Genome & genome );
extern double	number_weight		( Genome & genome );
extern double	number_area_run		( Genome & genome );
extern double	number_area_star	( Genome & genome );
extern double	number_area_bias	( Genome & genome );
extern double	number_area_length	( Genome & genome );

extern double	number_max		( Genome & genome );
extern double	number_area_max		( Genome & genome );
extern double	number_ave		( Genome & genome );
extern double	number_area_ave		( Genome & genome );

extern double	test_measure		( Genome & genome );

#endif	// __GENE_H_
