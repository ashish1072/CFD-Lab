#include "sor.h"
#include <math.h>
#include "parallel.h"

void sor(
  double omg,
  double dx,
  double dy,
  int imax,
  int jmax,
  int    il,
  int    ir,
  int    jb,
  int    jt,
  int rank_l,
  int rank_r, 
  int rank_b, 
  int rank_t,
  double **P,
  double **RS,
  double *res,
  double *bufSend,
  double *bufRecv,
  MPI_Status *status,
  int chunk
) {
  int i,j;
  double p_rloc;
  double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));
  double rloc;
  /* Set Pressure boundary values */

// Left Boundary
  if (rank_l == MPI_PROC_NULL)
  {  for(j = 1; j <= jt-jb +1; j++)
     { P[0][j] = P[1][j];}
  }

// Right Boundary
  if (rank_r == MPI_PROC_NULL)
  {  for(j = 1; j <= jt-jb +1; j++)
     { P[ir-il +2][j] = P[ir-il +1][j];}
  }

// Top Boundary
  if (rank_t == MPI_PROC_NULL)
  {  for(i = 1; i <= ir-il +1; i++)
     { P[i][jt-jb +2] = P[i][jt-jb +1];}
  }

// Bottom Boundary
  if (rank_b == MPI_PROC_NULL)
  {  for(i = 1; i <= ir-il +1; i++)
     { P[i][0] = P[i][1];}
  }
  /* SOR iteration */
  for(i = 1; i <= ir-il +1; i++) {
    for(j = 1; j<= jt-jb +1; j++) {
      P[i][j] = (1.0-omg)*P[i][j]
              + coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);
    }
  }
  // exchanging pressure values between subdomains
  pressure_comm( P, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, bufSend, bufRecv, status, chunk);


  /* compute the residual */
  for(i = 1; i <= ir-il +1; i++) {
    for(j = 1; j <= jt-jb +1; j++) {
      p_rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
              ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);
    }
  }
  // sum of individual process residuals and broadcasting it to all the processes  
  MPI_Allreduce( &p_rloc, &rloc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  rloc = rloc/(imax*jmax);                                   
  rloc = sqrt(rloc);
  //printf("rloc = %f\n",rloc);
  
  
  /* set residual */
  *res = rloc;


}


