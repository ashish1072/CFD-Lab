#include "uvp.h"
#include<math.h>
#include<stdio.h>

void calculate_dt(
  double Re,
  double tau,
  double *dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V
)
{

    int i;
    int j;
    double maxU = 0.0;
    double maxV = 0.0;

    for(i = 1 ; i <= imax ; i++)  //find umax and vmax
        for(j = 1 ; j <= jmax ; j++)
        {
            if (U[i][j] > maxU)
              maxU = U[i][j];

            if (V[i][j] > maxU)
              maxV = V[i][j];
        }
    if ( tau < 0)
      *dt = 0.05;
    else
      *dt = tau * ( fmin ( fmin ( (Re/2) * ( 1 / ( (1/pow(dx,2) ) + (1/pow(dy,2) ) ) ) , ( dx / maxU) ) , ( dy / maxV) ) );

}


void calculate_fg(
  double Re,
  double GX,
  double GY,
  double alpha,
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G
)
{

  int i;
  int j;

// F Calculation

for (i=1; i<= imax-1; i++)
  for (j=1; j<=jmax; j++)
      F[i][j]= U[i][j] + dt*
             (  1/Re*((1.0/(dx*dx)*( U[i+1][j] - 2*U[i][j] + U[i-1][j] )) + (1.0/(dy*dy)*( U[i][j+1] - 2*U[i][j] + U[i][j-1] )))

         - (1.0/dx*0.25*( (U[i][j] + U[i+1][j])*(U[i][j] + U[i+1][j]) - (U[i-1][j] + U[i][j])*(U[i-1][j] + U[i][j]))
             + alpha*1.0/dx*0.25*( fabs(U[i][j] + U[i+1][j])*(U[i][j] - U[i+1][j]) - fabs(U[i-1][j] + U[i][j])*(U[i-1][j] - U[i][j])) )

             - (1.0/dy*0.25*( (V[i][j] + V[i+1][j])*(U[i][j] + U[i][j+1]) - (V[i][j-1] + V[i+1][j-1])*(U[i][j-1] + U[i][j]))
             + alpha*1.0/dy*0.25*( fabs(V[i][j] + V[i+1][j])*(U[i][j] - U[i][j+1]) - fabs(V[i][j-1] + V[i+1][j-1])*(U[i][j-1] - U[i][j])) )

         + GX );

/* G  Calculation */
for (i=1; i<= imax; i++)
  for (j=1; j<=jmax-1; j++)
     G[i][j]= V[i][j] + dt*
              (  1/Re*((1.0/(dx*dx)*( V[i+1][j] - 2*V[i][j] + V[i-1][j] )) + (1.0/(dy*dy)*( V[i][j+1] - 2*V[i][j] + V[i][j-1] )))

          - (1.0/dy*0.25*( (V[i][j] + V[i][j+1])*(V[i][j] + V[i][j+1]) - (V[i][j-1] + V[i][j])*(V[i][j-1] + V[i][j]))
              + alpha*1.0/dy*0.25*( fabs(V[i][j] + V[i][j+1])*(V[i][j] - V[i][j+1]) - fabs(V[i][j-1] + V[i][j])*(V[i][j-1] - V[i][j])) )

              - (1.0/dx*0.25*( (U[i][j] + U[i][j+1])*(V[i][j] + V[i+1][j]) - (U[i-1][j] + U[i-1][j+1])*(V[i-1][j] + V[i][j]))
              + alpha*1.0/dx*0.25*( fabs(U[i][j] + U[i][j+1])*(V[i][j] - V[i+1][j]) - fabs(U[i-1][j] + U[i-1][j+1])*(V[i-1][j] - V[i][j])) )

            + GY );

/* Boundry Values for F */

  for( j = 1 ; j <= jmax ; j++ )
  {
    F[0][j] = U[0][j];
    F[imax][j] = U[imax][j];
  }

/* Boundry Values for G */

  for( i = 1 ; i <= imax ; i++ )
  {
    G[i][0] = V[i][0];
    G[i][jmax] = V[i][jmax];
  }
}

void calculate_rs(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **F,
  double **G,
  double **RS
)
{
	int i;
  int j;

  for( i = 1 ; i <= imax ; i++ )
    for( j = 1 ; j <= jmax ; j++ )
	   RS[i][j] = (1.0/dt)*( (1.0/dx)*( F[i][j] - F[i-1][j])  +  (1.0/dy)*( G[i][j] - G[i][j-1]) );

}

void calculate_uv(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G,
  double **P
)
{
	int i;
  int j;

  for( i = 1 ; i <= imax-1 ; i++ )
	 for( j = 1 ; j <= jmax ; j++ )
	  U[i][j] = F[i][j] - (dt/dx)*( P[i+1][j] - P[i][j]);

	for( i = 1 ; i <= imax ; i++ )
	 for( j = 1 ; j <= jmax-1 ; j++ )
	  V[i][j] = G[i][j] - (dt/dy)*( P[i][j+1] - P[i][j]);
}
