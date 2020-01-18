#include "sor.h"
#include <math.h>

void sor(
  double omg,
  double dx,
  double dy,
  int    imax,
  int    jmax,
  double **P,
  double **RS,
  double *res,
  int **Flag)
  {
  int i,j;
  double rloc;
  double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));
  int num_of_fluid_cells = 0;


  // Pressure values at Boundary/Obstacle cells

for(i = 0 ; i <= imax + 1 ; i++ )
		for(j = 0 ; j <= jmax + 1 ; j++ )
			if(Flag[i][j] > 16)
      {
				switch(Flag[i][j]&480)
				{
					case 32:                                  // For Northern neighbor
                  P[i][j] = P[i][j+1];
                  break;
					case 64:                                  //For Southern Neighbor
                  P[i][j] = P[i][j-1];
                  break;
					case 128:                                 // For Western Neighbor
                  P[i][j] = P[i-1][j];
                  break;
					case 256:                                 //For Eastern Neighbor
                  P[i][j] = P[i+1][j];
                  break;
					case 288:                                 //For NorthEastern Neighbor
                  P[i][j] = (P[i][j+1] + P[i+1][j])/2;
                  break;
					case 192:                                  //For SouthWestern Neighbor
                  P[i][j] = (P[i-1][j] + P[i][j-1])/2;
                  break;
					case 320:                                  //For SouthEastern Neighbor
                  P[i][j] = (P[i+1][j] + P[i][j-1])/2;
                  break;
					case 160:                                  //For NorthWestern Neighbor
                  P[i][j] = (P[i][j+1] + P[i-1][j])/2;
                  break;
         /* default:
                  P[i][j] = 0;
                  break; */
				}
      }


  /* SOR iteration */
  for(i = 1; i <= imax; i++) {
    for(j = 1; j<=jmax ; j++) {
      if(Flag[i][j]&1)
        P[i][j] = (1.0-omg)*P[i][j]
              + coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);
    }
  }

//Calculate number of fluid cells

  /* for( i = 1 ; i <= imax ; i++ )         
    for( j = 1 ; j <= jmax ; j++ )
      if(Flag[i][j]&1)
        num_of_fluid_cells++; */


  /* compute the residual */
  rloc = 0;
  for(i = 1; i <= imax; i++)
  {
    for(j = 1; j <= jmax ; j++)
    {
      if(Flag[i][j]&1)
       {
        num_of_fluid_cells++;                //***** can simply put num of cells in below residual calc loop
        rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
                ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);
       }
      }
    }

  rloc = rloc/(num_of_fluid_cells);
  rloc = sqrt(rloc);
  /* set residual */
  *res = rloc;



}

