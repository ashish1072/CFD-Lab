#include "helper.h"
#include "visual.h"
#include "boundary_val.h"
#include "init.h"
#include "sor.h"
#include "uvp.h"
#include <stdio.h>
#include<mpi.h>
#include"parallel.h"
#include<sys/stat.h>

void print(double ** mat,int il, int ir, int jb, int jt)
{
  for(int i = 0 ; i < ir-il+1; ++i )
  {
    printf("\n");
    for( int j = 0 ; j < jt-jb+1; ++j )
      printf("%f ",mat[i][j]);
  }
}
/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed using the matrix() command
 * - create the initial setup init_uvp(), init_flag(), output_uvp()
 * - perform the main loop
 * - trailer: destroy memory allocated and do some statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop the following big steps are done (for some of the
 * operations a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */
int main(int argn, char** args){

    const char* file = "cavity100.dat"; 
    char problem[100];
    double Re;                /* reynolds number   */
    double UI;                /* velocity x-direction */
    double VI;                /* velocity y-direction */
    double PI;                /* pressure */
    double GX;                /* gravitation x-direction */
    double GY;                /* gravitation y-direction */
    double t_end;             /* end time */
    double xlength;           /* length of the domain x-dir.*/
    double ylength;           /* length of the domain y-dir.*/
    double dt;                /* time step */
    double dx;                /* length of a cell x-dir. */
    double dy;                /* length of a cell y-dir. */
    int  imax;                /* number of cells x-direction*/
    int  jmax;                /* number of cells y-direction*/
    int iproc;                /* Number of processes in i dir*/
    int jproc;                /* Number of processes in j dir */
    double alpha;             /* uppwind differencing factor*/
    double omg;               /* relaxation factor */
    double tau;               /* safety factor for time step*/
    int  itermax;             /* max. number of iterations  */
    double eps;               /* accuracy bound for pressure*/
    double dt_value;          /* time for output */
    int myrank;               /* Process Id*/
    int rank_l;               /* Left neighbor */
    int rank_r;               /* right neighbor */
    int rank_t;               /* Top neighbor */
    int rank_b;               /* Bottom neighbor */
    int omg_i;                /* ith index of process */
    int omg_j;                /* jth index of process */
    int il;                   
    int ir;                   
    int jb;
    int jt;
    double **U = NULL;
    double **V = NULL;
    double **P = NULL;
    double **F = NULL;
    double **G = NULL;
    double **RS = NULL;
    int num_proc;
    int it;
    double t = 0.0;
    int n = 0;
    double res;
    int sor_break = 0;
    int visual_iter = 0;
    struct stat st = {0};
    char sol_folder[80];
    char sol_directory[80];    
    MPI_Status status;
    int chunk = 0;

    MPI_Init(&argn, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    /* Reading parameters from .dat file */
    read_parameters(file, problem, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau, &itermax, &eps, &dt_value, &iproc, &jproc);
    
    /* Initialization of parallel processes with x,y dimensions and neighbors */
    init_parallel(iproc, jproc, imax, jmax, &myrank, &il, &ir, &jb, &jt, &rank_l, &rank_r, &rank_b, &rank_t ,&omg_i ,&omg_j , &dx, &dy, xlength, ylength);

    /* Memory allocation of Matrices for U,V,P Calulation*/
    P = matrix(0, ir-il+2, 0, jt-jb+2);
    U = matrix(0, ir-il+3 , 0 , jt-jb+2);
    V = matrix(0, ir-il+2, 0, jt-jb+3);
    F = matrix(0, ir-il+3 , 0 , jt-jb+2);
    G = matrix(0, ir-il+2, 0, jt-jb+3);
    RS = matrix(0, ir-il+2 , 0, jt-jb+2);
	
    double *bufSend = (double*)malloc((ir-il+5)*sizeof(double)); 
    double *bufRecv = (double*)malloc((ir-il+5)*sizeof(double));
    
    //Make solution folder
    sprintf( sol_folder,"Solution_%s",problem);
    if (stat(sol_folder, &st) == -1) 
    	mkdir(sol_folder, 0700);
	
    sprintf( sol_directory,"Solution_%s/sol", problem);

    /* Initialization of U,V,P Matrices with UI, VI, PI values*/
    init_uvp(UI, VI, PI, il, ir, jb, jt, U, V, P);

    while( t < t_end)
    { /* setting up of velocity boundary values before each iteration step*/
      boundaryvalues(il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, U, V);
      calculate_fg(Re, GX, GY, omg, dt, dx, dy, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, U, V, F, G);
      calculate_rs(dt, dx, dy, il, ir, jb, jt, F, G, RS);

      res = 50;  // random value of residual that is greater than eps
      it = 0;
      printf("time = %f and dt= %f\n",t,dt);
      
      /* Calculation of pressure using SOR */
      while( it < itermax && res > eps)
      {
        sor(omg, dx, dy, imax, jmax, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, P, RS, &res, bufSend, bufRecv, &status, chunk);
        it++;
      }

      if(it == itermax)
	  sor_break++;

      /* U and V calculation using discritized N-S Equation*/
      calculate_uv(dt, dx, dy, il, ir, jb, jt,rank_l,rank_r,rank_b,rank_t,&status,chunk, U, V, F, G, P, bufSend,bufRecv);

      if(t >= visual_iter*dt_value)
     {/* generation of vtk files for visualization */
       output_uvp(U,V,P,il,ir,jb, jt,omg_i,omg_j,dx, dy, sol_directory, n);
       visual_iter++;
  	 }
  	  /* calculation of time-step size*/
      calculate_dt(Re, tau, &dt, dx, dy, il, ir, jb, jt, U, V);
      
      t += dt;
      n++;
    }
    
    /* Allocated memory freed*/
    free_matrix(P, 0, ir-il+2, 0, jt-jb+2);
    free_matrix(U, 0, ir-il+3, 0, jt-jb+2);
    free_matrix(V, 0, ir-il+2, 0, jt-jb+3);
    free_matrix(F, 0, ir-il+3, 0, jt-jb+2);
    free_matrix(G, 0, ir-il+2, 0, jt-jb+3);
    free_matrix(RS, 0, ir-il+2 , 0, jt-jb+2);
    free(bufSend);
    free(bufRecv);
   
  Programm_Stop("End of the Program"); 
  // End MPI Session  
  MPI_Finalize();     
  
  return -1;
}
