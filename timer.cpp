#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>

static double startETime, startCPUTime;

double getETime( void )
{
    struct timeval tv;
    gettimeofday( &tv, NULL );
    return tv.tv_sec + ( double )tv.tv_usec*1e-6;
}

void checkInETime( void )
{
    startETime = getETime();
}


double checkOutETime( void )
{
    double curETime = getETime();
    return ( curETime - startETime );
}    


double getCPUTime( void )
{
    struct rusage RU;
    getrusage( RUSAGE_SELF, &RU );
    return RU.ru_utime.tv_sec + ( double )RU.ru_utime.tv_usec*1e-6;
}


void checkInCPUTime( void )
{
    startCPUTime = getCPUTime();
}


double checkOutCPUTime( void )
{
    double curCPUTime = getCPUTime();
    return ( curCPUTime - startCPUTime );
}    

void reportTime( void )
{
    fprintf( stderr, "Time passsed = %8.3f, CPU time = %8.3f\n", checkOutETime(),checkOutCPUTime() );
}
