#include "uvp.h"
#include<math.h>
#include<stdio.h>
#include<mpi.h>
#include "parallel.h"

void calculate_dt(        
  double Re,
  double tau,
  double *dt,
  double dx,
  double dy,
  int il,                             
  int ir,                                 
  int jb,
  int jt,
  double **U,
  double **V
)
{

    int i;
    int j;
    double proc_maxU = 0.0;
    double proc_maxV = 0.0;
    int myrank;

    for(i = 0 ; i <= (ir-il +2) ; i++)                 
    {   for(j = 0 ; j <= (jt-jb +2) ; j++)
        {
            if (U[i][j] > proc_maxU)
              proc_maxU = U[i][j];

            if (V[i][j] > proc_maxU)
              proc_maxV = V[i][j];
        }
    }
    double maxU = 0.0;
    double maxV = 0.0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    //printf("PROC Umax = %f PROC Vmax= %f   PID = %d\n\n ", proc_maxU,proc_maxV, myrank);
    MPI_Allreduce( &proc_maxU, &maxU, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); // maximum U sent to all the processes.
    MPI_Allreduce( &proc_maxV, &maxV, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); // maximum V sent to all the processes.
    //printf("PROC Umax = %f PROC Vmax= %f\n\n ", proc_maxU,proc_maxV);
    if ( tau < 0)
      *dt = 0.05;
    else
      *dt = tau * ( fmin ( fmin ( (Re/2) * ( 1 / ( (1/pow(dx,2) ) + (1/pow(dy,2) ) ) ) , ( dx / maxU) ) , ( dy / maxV) ) );
    //printf("Umax = %f, Vmax= %f, dt = %f, PID = %d\n\n", maxU, maxV, *dt,myrank);
}


void calculate_fg(
  double Re,
  double GX,
  double GY,
  double alpha,
  double dt,
  double dx,
  double dy,
  int il,
  int ir,
  int jb,
  int jt,
  int rank_l,
  int rank_r,
  int rank_b,
  int rank_t,
  double **U,
  double **V,
  double **F,
  double **G
)
{

  int i;
  int j;
  //int fx_start=0, fx_end=0, fy_start=0,fy_end=0,gx_start=0, gx_end=0, gy_start=0,gy_end=0;
  
// F Calculation
for (i = 1; i < ir+1-(il-2) ; i++){
  for (j = 1; j < jt+1-(jb-1) ; j++){
  if (rank_l==MPI_PROC_NULL && rank_r==MPI_PROC_NULL && (i==1 || i==ir-(il-2))){
        F[1][j] = U[1][j];
        F[ir-(il-2)][j] = U[ir-(il-2)][j];
  }
  else if (rank_l==MPI_PROC_NULL && i==1){
    F[1][j] = U[1][j];
  }    
  else if(rank_r==MPI_PROC_NULL && i==ir-(il-2)){
    F[ir-(il-2)][j] = U[ir-(il-2)][j];
  }
  else{    
      F[i][j]= U[i][j] + dt*
             (  1/Re*((1.0/(dx*dx)*( U[i+1][j] - 2*U[i][j] + U[i-1][j] )) + (1.0/(dy*dy)*( U[i][j+1] - 2*U[i][j] + U[i][j-1] )))

         - (1.0/dx*0.25*( (U[i][j] + U[i+1][j])*(U[i][j] + U[i+1][j]) - (U[i-1][j] + U[i][j])*(U[i-1][j] + U[i][j]))
             + alpha*1.0/dx*0.25*( fabs(U[i][j] + U[i+1][j])*(U[i][j] - U[i+1][j]) - fabs(U[i-1][j] + U[i][j])*(U[i-1][j] - U[i][j])) )

             - (1.0/dy*0.25*( (V[i-1][j+1] + V[i][j+1])*(U[i][j] + U[i][j+1]) - (V[i-1][j] + V[i][j])*(U[i][j-1] + U[i][j]))
             + alpha*1.0/dy*0.25*( fabs(V[i-1][j+1] + V[i][j+1])*(U[i][j] - U[i][j+1]) - fabs(V[i-1][j] + V[i][j])*(U[i][j-1] - U[i][j])) )

         + GX );
         }
        }
      }
      
    

/* G  Calculation */
for (i = 1 ; i < ir+1-(il-1); i++){
  for (j = 1 ; j < jt+1-(jb-2); j++){
  if (rank_t==MPI_PROC_NULL && rank_b==MPI_PROC_NULL && (j==1 || j==jt-(jb-2))){
        G[i][1] = V[i][1];
        G[i][jt-(jb-2)] = V[i][jt-(jb-2)];
  }
  else if (rank_t==MPI_PROC_NULL && j==1){
    G[i][1] = V[i][1];
  }    
  else if(rank_r==MPI_PROC_NULL && i==ir-(il-2)){
    G[i][jt-(jb-2)] = V[i][jt-(jb-2)];
  }
  else{    
     
     
     G[i][j]= V[i][j] + dt*
              (  1/Re*((1.0/(dx*dx)*( V[i+1][j] - 2*V[i][j] + V[i-1][j] )) + (1.0/(dy*dy)*( V[i][j+1] - 2*V[i][j] + V[i][j-1] )))

          - (1.0/dy*0.25*( (V[i][j] + V[i][j+1])*(V[i][j] + V[i][j+1]) - (V[i][j-1] + V[i][j])*(V[i][j-1] + V[i][j]))
              + alpha*1.0/dy*0.25*( fabs(V[i][j] + V[i][j+1])*(V[i][j] - V[i][j+1]) - fabs(V[i][j-1] + V[i][j])*(V[i][j-1] - V[i][j])) )

              - (1.0/dx*0.25*( (U[i+1][j-1] + U[i+1][j])*(V[i][j] + V[i+1][j]) - (U[i][j-1] + U[i][j])*(V[i-1][j] + V[i][j]))
              + alpha*1.0/dx*0.25*( fabs(U[i+1][j-1] + U[i+1][j])*(V[i][j] - V[i+1][j]) - fabs(U[i][j-1] + U[i][j])*(V[i-1][j] - V[i][j])) )

            + GY );
  }
  }
  }

}

void calculate_rs(
  double dt,
  double dx,
  double dy,
  int il,
  int ir,
  int jb,
  int jt,
  double **F,
  double **G,
  double **RS
)
{
	int i;
  int j;

  for( i = 1 ; i <= ir-il+1 ; i++ )
    for( j = 1 ; j <= jt-jb+1 ; j++ )
	   RS[i][j] = (1.0/dt)*( (1.0/dx)*( F[i+1][j] - F[i][j])  +  (1.0/dy)*( G[i][j+1] - G[i][j]) );

}

void calculate_uv(
  double dt,
  double dx,
  double dy,
  int il,
  int ir,
  int jb,
  int jt,
  int rank_l,
  int rank_r,
  int rank_b,
  int rank_t,
  MPI_Status *status,
  int chunk,
  double **U,
  double **V,
  double **F,
  double **G,
  double **P,
  double *bufSend,
  double *bufRecv
)
{
	int i;
  int j;
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  //printf("I am calculating UV, my rank is %d\n",myrank);
  for( i = 1 ; i <= ir-il+2 ; i++ ){
	 for( j = 1 ; j <= jt-jb+1 ; j++ ){
	  U[i][j] = F[i][j] - (dt/dx)*( P[i][j] - P[i-1][j]);
    }
    }

	for( i = 1 ; i <= ir-il+1 ; i++ )
	 for( j = 1 ; j <= jt-jb+2 ; j++ )
	  V[i][j] = G[i][j] - (dt/dy)*( P[i][j] - P[i][j-1]);
  uv_comm(U, V, il, ir, jb,jt, rank_l,rank_r,rank_b,rank_t,bufSend,bufRecv,status, chunk);
 // printf("UV calculated for rank %d\n",myrank);
}
