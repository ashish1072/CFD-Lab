#include "boundary_val.h"
#include"mpi.h"

void boundaryvalues(int il,int ir,int jb,int jt, int rank_l, int rank_r, int rank_b, int rank_t, double **U, double **V)
{
  int i;
  int j;

  // Boundary Conditions on Left Boundary
  if(rank_l == MPI_PROC_NULL)
	{		for ( j = 1 ; j <= jt-jb +1 ; j++ )
				{ U[1][j] = 0.0;
		      V[0][j] = -V[1][j];
			  }
	}

  // Boundary Conditions on Right Boundary
	if(rank_r == MPI_PROC_NULL)
	{   for ( j = 1 ; j <= jt-jb +1 ; j++ )
			  { U[ir-il +2][j] = 0.0;
				  V[ir-il +2][j] = -V[ir-il +1][j];
				}
	}

  // Boundary Conditions on Top Boundary
	if(rank_t == MPI_PROC_NULL)
	{   for ( i = 0 ; i <= ir-il +2 ; i++ )
			  { U[i][jt-jb +2] = 2 - U[i][jt-jb +1];
				  V[i][jt-jb +2] = 0.0;
				}
	}
		
	// Boundary Conditions on Bottom Boundary
	if(rank_b == MPI_PROC_NULL)
	{   for ( i = 1 ; i <= ir-il +1 ; i++ )
			  { U[i][0] = -U[i][1];
				  V[i][1] = 0.0;
				}
	}


}
