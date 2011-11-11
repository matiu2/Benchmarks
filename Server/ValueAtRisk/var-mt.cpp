#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>

#define NS 10
#define CORES 8

inline uint32_t initialSeed( size_t index )
{
  uint32_t seed = 1;
  uint32_t mult = 300773;
  size_t mask = 1;
  while ( mask != 0 )
  {
    if ( index & mask )
      seed = (seed * mult) % 1073741824;
    mult = (mult * mult) % 1073741824;
    mask <<= 1;
  }
  return seed;
}

inline double nextUniform01( uint32_t &seed )
{
  seed = (seed * 300773) % 1073741824;
  //printf( "%u\n", (unsigned)seed );
  return double(seed) / double(1073741824.0);
}

inline double randomNormal( uint32_t &seed )
{
  double x1, x2, w;
  do
  {
    x1 = 2.0 * nextUniform01( seed ) - 1.0;
    x2 = 2.0 * nextUniform01( seed ) - 1.0;
    w = x1 * x1 + x2 * x2;
  } while ( w >= 1.0 );
  w = sqrt( (-2.0 * log(w)) / w );
  return x1 * w;
}

inline void randomNormalVec( double vec[NS], uint32_t &seed )
{
  for ( size_t i=0; i<NS; ++i )
    vec[i] = randomNormal( seed );
}

inline void multMatVec( double const mat[NS][NS], double const vec[NS], double res[NS] )
{
  for ( size_t i=0; i<NS; ++i )
  {
    res[i] = 0.0;
    for ( size_t j=0; j<NS; ++j )
      res[i] += mat[i][j] * vec[j];
  }
}

inline double runTrial(
  size_t index,
  size_t numTradingDays,
  double dt,
  double sqrtDT,
  double const choleskyTrans[NS][NS],
  double const drifts[NS]
  )
{
  uint32_t seed = initialSeed( 4096 * (1+index) );

  double amounts[NS];
  for ( size_t i=0; i<NS; ++i )
    amounts[i] = 100.0;

  for ( size_t day=0; day<numTradingDays; ++day )
  {
    double Z[NS];
    randomNormalVec( Z, seed );
    double X[NS];
    multMatVec( choleskyTrans, Z, X );
    for ( size_t i=0; i<NS; ++i )
      amounts[i] *= exp(drifts[i]*dt + X[i]*sqrtDT);
  }

  double value = 0.0;
  for ( size_t i=0; i<NS; ++i )
    value += amounts[i];
  return value;
}

struct FixedArgs
{
  size_t numTradingDays;
  double dt;
  double sqrtDT;
  double choleskyTrans[NS][NS];
  double drifts[NS];
};

struct Args
{
  FixedArgs *fixedArgs;
  size_t startIndex;
  size_t endIndex;
  double *trialResults;
};

void *threadEntry( void *_args )
{
  Args const *args = (Args const *)_args;
  for ( size_t index = args->startIndex;
    index != args->endIndex; ++index )
  {
    args->trialResults[index] = runTrial(
      index,
      args->fixedArgs->numTradingDays,
      args->fixedArgs->dt,
      args->fixedArgs->sqrtDT,
      args->fixedArgs->choleskyTrans,
      args->fixedArgs->drifts
    );
  }
  return 0;
}

inline void trans( double A[NS][NS], double B[NS][NS] )
{
  for ( size_t i=0; i<NS; ++i )
  {
    for ( size_t j=0; j<NS; ++j )
    {
      B[i][j] = A[j][i];
    }
  }
}

inline void multMatMat( double A[NS][NS], double B[NS][NS], double R[NS][NS] )
{
  for ( size_t i=0; i<NS; ++i )
  {
    for ( size_t j=0; j<NS; ++j )
    {
      R[i][j] = 0.0;
      for ( size_t k=0; k<NS; ++k )
        R[i][j] += A[i][k] * B[k][j];
    }
  }
}

void randomCorrelation( double R[NS][NS], uint32_t &seed )
{
  double T[NS][NS];
  for ( size_t i=0; i<NS; ++i )
  {
    for ( size_t j=0; j<NS; ++j )
    {
      T[i][j] = randomNormal( seed );
    }
  }

  for ( size_t j=0; j<NS; ++j )
  {
    double sqSum = 0.0;
    for ( size_t i=0; i<NS; ++i )
    {
      sqSum += T[i][j] * T[i][j];
    }
    double norm = sqrt( sqSum );
    for ( size_t i=0; i<NS; ++i )
      T[i][j] /= norm;
  }

  double TTrans[NS][NS];
  trans( T, TTrans );

  multMatMat( TTrans, T, R );
}

void computeCholeskyTrans( double A[NS][NS], double B[NS][NS] )
{
  for ( size_t i=0; i<NS; ++i )
    for ( size_t j=0; j<NS; ++j )
      B[i][j] = 0.0;

  for ( size_t i=0; i<NS; ++i )
  {
    for ( size_t j=0; j<i+1; ++j )
    {
      double s = 0.0;
      for ( size_t k=0; k<j; ++k )
        s += B[i][k] * B[j][k];
      if ( i == j )
        B[i][i] = sqrt( A[i][i] - s );
      else
        B[i][j] = 1.0 / B[j][j] * (A[i][j] - s);
    }
  }
}

int doubleCompare( void const *_lhs, void const *_rhs )
{
  double lhs = *(double const *)_lhs;
  double rhs = *(double const *)_rhs;
  if ( lhs < rhs )
    return -1;
  else if ( lhs > rhs )
    return 1;
  else return 0;
}

int main( int argc, char **argv )
{
  size_t const numTrials = 1048576;
  double *trialResults = new double[numTrials];

  FixedArgs fixedArgs;
  fixedArgs.numTradingDays = 252;
  fixedArgs.dt = 1.0 / fixedArgs.numTradingDays;
  fixedArgs.sqrtDT = sqrt( fixedArgs.dt );

  double priceMeans[NS];
  for ( size_t i=0; i<NS; ++i )
    priceMeans[i] = 25.0/fixedArgs.numTradingDays;

  double priceDevs[NS];
  for ( size_t i=0; i<NS; ++i )
    priceDevs[i] = 25.0/fixedArgs.numTradingDays;

  uint32_t seed = initialSeed(0);

  double priceCorrelations[NS][NS];
  randomCorrelation( priceCorrelations, seed );

  double priceCovariance[NS][NS];
  for ( size_t i=0; i<NS; ++i )
  {
    for ( size_t j=0; j<NS; ++j )
    {
      priceCovariance[i][j] = priceDevs[i] * priceDevs[j] * priceCorrelations[i][j];
    }
  }

  computeCholeskyTrans( priceCovariance, fixedArgs.choleskyTrans );

  for ( size_t i=0; i<NS; ++i )
    fixedArgs.drifts[i] = priceMeans[i] - priceCovariance[i][i]/2.0;

  Args args[CORES];
  pthread_t threads[CORES-1];
  for ( size_t core=0; core<CORES; ++core )
  {
    if ( core == 0 )
      args[core].startIndex = 0;
    else
      args[core].startIndex = args[core-1].endIndex;
    if ( core+1 == CORES )
      args[core].endIndex = numTrials;
    else
      args[core].endIndex = args[core].startIndex + numTrials/CORES;
    args[core].fixedArgs = &fixedArgs;
    args[core].trialResults = trialResults;
    if ( core+1 == CORES )
      threadEntry( &args[core] );
    else
      pthread_create( &threads[core], 0, &threadEntry, &args[core] );
  }

  for ( size_t core=0; core<CORES-1; ++core )
    pthread_join( threads[core], 0 );

  qsort( trialResults, numTrials, sizeof(double), doubleCompare );

  printf( "VaR = %.16f\n", 100.0*NS - trialResults[(size_t)floor( 0.05 * numTrials )] );

  delete [] trialResults;

  return 0;
}
