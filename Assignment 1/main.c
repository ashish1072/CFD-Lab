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

    const char *file = "cavity100.dat";
    double *Re = malloc(sizeof(double*));                /* reynolds number   */
    double *UI = malloc(sizeof(double*));             /* velocity x-direction */
    double *VI = malloc(sizeof(double*));                /* velocity y-direction */
    double *PI = malloc(sizeof(double*));                /* pressure */
    double *GX = malloc(sizeof(double*));                /* gravitation x-direction */
    double *GY = malloc(sizeof(double*));                /* gravitation y-direction */
    double *t_end = malloc(sizeof(double*));             /* end time */
    double *xlength = malloc(sizeof(double*));           /* length of the domain x-dir.*/
    double *ylength = malloc(sizeof(double*));           /* length of the domain y-dir.*/
    double *dt = malloc(sizeof(double*));                /* time step */
    double *dx = malloc(sizeof(double*));                /* length of a cell x-dir. */
    double *dy = malloc(sizeof(double*));               /* length of a cell y-dir. */
    int  *imax = malloc(sizeof(int*));                /* number of cells x-direction*/
    int  *jmax = malloc(sizeof(int*));                /* number of cells y-direction*/
    double *alpha = malloc(sizeof(double*));             /* uppwind differencing factor*/
    double *omg = malloc(sizeof(double*));               /* relaxation factor */
    double *tau = malloc(sizeof(double*));               /* safety factor for time step*/
    int  *itermax = malloc(sizeof(int*));             /* max. number of iterations  */
    double *eps = malloc(sizeof(double*));               /* accuracy bound for pressure*/
    double *dt_value = malloc(sizeof(double*));
    double **U = NULL;
    double **V = NULL;
    double **P = NULL;
    double **F = NULL;
    double **G = NULL;
    double **RS = NULL;
    int it;
    double t = 0.0;
    int n = 0;
    double *res = malloc(sizeof(double*));

    read_parameters(file, Re, UI, VI, PI, GX, GY, t_end, xlength, ylength, dt, dx, dy, imax, jmax, alpha, omg, tau, itermax, eps, dt_value );
    
    U = matrix(0, (*imax) + 2, 0, (*jmax) + 2);
    V = matrix(0, (*imax) + 2, 0, (*jmax) + 2);
    P = matrix(0, (*imax) + 2, 0, (*jmax) + 2);
    F = matrix(0, (*imax) + 2, 0, (*jmax) + 2);
    G = matrix(0, (*imax) + 2, 0, (*jmax) + 2);
    RS = matrix(0, (*imax) + 2, 0, (*jmax) + 2);

    init_uvp(*UI, *VI, *PI, *imax, *jmax, U, V, P);

    while( t < *t_end)
    {

      calculate_dt(*Re, *tau, dt, *dx, *dy, *imax, *jmax, U, V);
      boundaryvalues(*imax, *jmax, U, V);
      calculate_fg(*Re, *GX, *GY, *alpha, *dt, *dx, *dy, *imax, *jmax, U, V, F, G);
      calculate_rs(*dt, *dx, *dy, *imax, *jmax, F, G, RS);

      *res = 50;  // random value of residual that is greater than eps
      it = 0;
      while( it < *itermax && *res > *eps)
      {
        sor(*omg, *dx, *dy, *imax, *jmax, P, RS, res);
        it++;
      }

      calculate_uv(*dt, *dx, *dy, *imax, *jmax, U, V, F, G, P);

      write_vtkFile("szProblem", n, *xlength, *ylength, *imax, *jmax, *dx, *dy, U, V, P);

      t += (*dt);
      n++;
    }
    
    free_matrix( U, 0, (*imax) + 2, 0, (*jmax) + 2 );
    free_matrix( V, 0, (*imax) + 2, 0, (*jmax) + 2 );
    free_matrix( P, 0, (*imax) + 2, 0, (*jmax) + 2 );
    free_matrix( F, 0, (*imax) + 2, 0, (*jmax) + 2 );
    free_matrix( G, 0, (*imax) + 2, 0, (*jmax) + 2 );
    free_matrix( RS, 0, (*imax) + 2, 0, (*jmax) + 2 );

  return -1;
}
