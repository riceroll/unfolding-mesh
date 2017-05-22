#include "random.h"

static void bitseed(unsigned int seed=1);

// Return a string indicating which random number generator was compiled-in to
// the library.
const char* GetRNG( void ) {
  return "RAN2";
}

// Seed the random number generator with an appropriate value.  We seed both 
// the random number generator and the random bit generator.  Set the seed only
// if a seed is not specified.  If a seed is specified, then set the seed to 
// the specified value and use it.  We remember the seed so that multiple calls
// to this function with the same seed do not reset the generator.  Subsequent
// calls to this function with a different seed will initialize the generator
// to the new seed.  Multiple calls with a value of 0 do nothing (we do *not*
// re-seed the generator because 0 is the default value and we don't want 
// people to re-seed the generator inadvertantly).
//   Some systems return a long as the return value for time, so we need to be
// sure to get whatever variation from it that we can since our seed is only an
// unsigned int.
static unsigned int seed=0;

unsigned int GetRandomSeed( void ) { return seed; }

#ifdef SKIP
void GARandomSeed(unsigned int s) {
  if(s == 0 && seed == 0) {
    unsigned long int tmp;
    while(seed == 0) {
      tmp = time(NULL) _GA_PID;
      for(unsigned int i=0; i<GALIB_BITS_IN_WORD*sizeof(unsigned int); i++)
	seed += (tmp & (1 << i));
    }
    _GA_RND_SEED (seed); 
    bitseed(seed);
  }
  else if(s != 0 && seed != s) {
    seed = s;
    _GA_RND_SEED (seed); 
    bitseed(seed);
  }
}
#endif	// SKIP

// Similar to setting the random seed, but this one sets it as long as the
// specified seed is non-zero.
void
ResetRNG(unsigned int s) {
  if(s != 0) {
    seed = s;
    _GA_RND_SEED (seed); 
    bitseed(seed);
  }
}

#ifdef SKIP
// Return a number from a unit Gaussian distribution.  The mean is 0 and the
// standard deviation is 1.0.
//   First we generate two uniformly random variables inside the complex unit 
// circle.  Then we transform these into Gaussians using the Box-Muller 
// transformation.  This method is described in Numerical Recipes in C 
// ISBN 0-521-43108-5 at http://world.std.com/~nr
//   When we find a number, we also find its twin, so we cache that here so 
// that every other call is a lookup rather than a calculation.  (I think GNU 
// does this in their implementations as well, but I don't remember for 
// certain.)
double GAUnitGaussian( void ) {
    static GABoolean cached=gaFalse;
  static double cachevalue;
  if(cached == gaTrue){
    cached = gaFalse;
    return cachevalue;
  }

  double rsquare, factor, var1, var2;
  do{
    var1 = 2.0 * GARandomDouble() - 1.0;
    var2 = 2.0 * GARandomDouble() - 1.0;
    rsquare = var1*var1 + var2*var2;
  } while(rsquare >= 1.0 || rsquare == 0.0);

  double val = -2.0 * log(rsquare) / rsquare;
  if(val > 0.0) factor = sqrt(val);
  else           factor = 0.0;	// should not happen, but might due to roundoff

  cachevalue = var1 * factor;
  cached = gaTrue;

  return (var2 * factor);
}
#endif	// SKIP



// This is the random bit generator Method II from numerical recipes in C.  The
// seed determines where in the cycle of numbers the generator will start, so
// we don't need full 'long' precision in the argument to the seed function.

#define IB1 1
#define IB2 2
#define IB5 16
#define IB18 131072L
#define MASK (IB1+IB2+IB5)

static unsigned long iseed;

void bitseed( unsigned int seed ) {
    iseed = seed;
}

int RandomBit()
{
    if (iseed & IB18) {
	iseed=((iseed ^ MASK) << 1) | IB1;
	return 1;
    } else {
	iseed <<= 1;
	return 0;
    }
}

#undef MASK
#undef IB18
#undef IB5
#undef IB2
#undef IB1



#define IM1 2147483563L
#define IM2 2147483399L
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014L
#define IA2 40692L
#define IQ1 53668L
#define IQ2 52774L
#define IR1 12211L
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

static long idum2=123456789;
static long iy=0;
static long iv[NTAB];
static long idum=0;

void 
sran2(unsigned int seed) {
  int j;
  long k;

  idum = STA_CAST(long,seed);
  if (idum == 0) idum=1;
  if (idum < 0) idum = -idum;
  idum2=(idum);
  for (j=NTAB+7;j>=0;j--) {
    k=(idum)/IQ1;
    idum=IA1*(idum-k*IQ1)-k*IR1;
    if (idum < 0) idum += IM1;
    if (j < NTAB) iv[j] = idum;
  }
  iy=iv[0];
}

float
ran2() {
  int j;
  long k;
  float temp;

  k=(idum)/IQ1;
  idum=IA1*(idum-k*IQ1)-k*IR1;
  if (idum < 0) idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = idum;
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
