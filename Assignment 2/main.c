#include "helper.h"
#include "visual.h"
#include "boundary_val.h"
#include "init.h"
#include "sor.h"
#include "uvp.h"
#include <stdio.h>

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

    const char *problem = "karmanvortex.dat";     /* The scenario to be computed */
    const char *geometry = "geometry.pgm";        /* The geometry file */
    double Re;				                            /* reynolds number   */
    double UI;             	                      /* velocity x-direction */
    double VI;                                    /* velocity y-direction */
    double PI;                                    /* pressure */
    double GX;    			           	              /* gravitation x-direction */
    double GY;               				              /* gravitation y-direction */
    double t_end;            				              /* end time */
    double xlength;          				              /* length of the domain x-dir.*/
    double ylength;         				              /* length of the domain y-dir.*/
    double dt;                                    /* time step */
    double dx;               				              /* length of a cell x-dir. */
    double dy;              				              /* length of a cell y-dir. */
    int  imax;              				              /* number of cells x-direction*/
    int  jmax;              				              /* number of cells y-direction*/
    double alpha;            			                /* uppwind differencing factor*/
    double omg;                                   /* relaxation factor */
    double tau;                                   /* safety factor for time step*/
    int itermax;                                  /* max. number of iterations  */
    double eps;                                   /* accuracy bound for pressure*/
    double dt_value;
    double **U = NULL;
    double **V = NULL;
    double **P = NULL;
    double **F = NULL;
    double **G = NULL;
    double **RS = NULL;
    int **Flag = NULL;
    int it;
    double t = 0.0;
    int n = 0;
    double res;
   
    read_parameters(problem, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau, &itermax, &eps, &dt_value );
    
    U = matrix(0, imax , 0, jmax + 1);
    V = matrix(0, imax + 1, 0, jmax );
    P = matrix(0, imax + 1, 0, jmax + 1);
    F = matrix(0, imax +1, 0, jmax + 1);
    G = matrix(0, imax + 1, 0, jmax +1);
    RS = matrix(0, imax + 1, 0, jmax + 1);
    Flag = imatrix(0, imax + 1, 0, jmax + 1);
    int visual_iter = 0;
    int num_file;
    int i,j;

    init_uvp(UI, VI, PI, imax, jmax, U, V, P, Flag);

    printf("U, V, P initialised... \n");

    init_flag(geometry, imax, jmax, Flag);

    printf(" Flag initialised... \n");

    for(i = 0 ; i <= imax + 1 ; i++)
    {
      printf("\n");
      for(j = 0 ; j <= jmax + 1 ; j++)
        printf("%d ",Flag[i][j]);
    }

    while( t < t_end)
    {
      calculate_dt(Re, tau, &dt, dx, dy, imax, jmax, U, V);
      printf("t = %f ,dt = %f, ",t,dt);
      //printf("dt done \n");
      
      boundaryvalues(imax, jmax, U, V, Flag);
     // printf("Boundary done \n");

      spec_boundary_val(imax, jmax ,U ,V ,Flag );
     // printf("Spec boundary done \n");

      calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, U, V, F, G, Flag);
      //printf("fg done \n");

      calculate_rs(dt, dx, dy, imax, jmax, F, G, RS,Flag);
      //printf("rs done \n");


      
      res = 50;  // random value of residual that is greater than eps
      it = 0;

      
      while( it < itermax && res > eps)
      {
        sor(omg, dx, dy, imax, jmax, P, RS, &res, Flag);
        it++;
      }

      printf("SOR itertions = %d, residual = %lf \n", it-1, res);

      calculate_uv(dt, dx, dy, imax, jmax, U, V, F, G, P,Flag);
      
      num_file = visual_iter*(dt_value);
      if(t >= num_file)
      {
      	write_vtkFile("karmanvortex", n, xlength, ylength, imax, jmax, dx, dy, U, V, P);
        visual_iter++;
      }

      t += (dt);
      n++;
    }
    
    free_matrix( U, 0, imax + 1, 0, jmax + 1 );
    free_matrix( V, 0, imax + 1, 0, jmax + 1 );
    free_matrix( P, 0, imax + 1, 0, jmax + 1 );
    free_matrix( F, 0, imax + 1, 0, jmax + 1 );
    free_matrix( G, 0, imax + 1, 0, jmax + 1 );
    free_matrix( RS, 0, imax + 1, 0, jmax + 1 );
    free_imatrix(Flag,0, imax + 1,0 , jmax + 1);
    
  return -1;
}
