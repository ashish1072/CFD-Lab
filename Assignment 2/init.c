#include "helper.h"
#include "init.h"
#include<string.h>

int read_parameters( const char *problem,       /* name of the file */
                    double *Re,                /* reynolds number   */
                    double *UI,                /* velocity x-direction */
                    double *VI,                /* velocity y-direction */
                    double *PI,                /* pressure */
                    double *GX,                /* gravitation x-direction */
                    double *GY,                /* gravitation y-direction */
                    double *t_end,             /* end time */
                    double *xlength,           /* length of the domain x-dir.*/
                    double *ylength,           /* length of the domain y-dir.*/
                    double *dt,                /* time step */
                    double *dx,                /* length of a cell x-dir. */
                    double *dy,                /* length of a cell y-dir. */
                    int  *imax,                /* number of cells x-direction*/
                    int  *jmax,                /* number of cells y-direction*/
                    double *alpha,             /* uppwind differencing factor*/
                    double *omg,               /* relaxation factor */
                    double *tau,               /* safety factor for time step*/
                    int  *itermax,             /* max. number of iterations  */
		                                       /* for pressure per time step */
                    double *eps,               /* accuracy bound for pressure*/
		    double *dt_value)                  /* time for output */
{



   READ_DOUBLE(problem, *xlength );
   READ_DOUBLE(problem, *ylength );
   READ_DOUBLE(problem, *dt    );
   READ_DOUBLE(problem, *t_end );
   READ_DOUBLE(problem, *tau   );
   READ_DOUBLE(problem, *dt_value );
   READ_INT   (problem, *itermax );
   READ_DOUBLE(problem, *eps   );
   READ_DOUBLE(problem, *omg   );
   READ_DOUBLE(problem, *alpha );
   READ_INT   (problem, *imax );
   READ_INT   (problem, *jmax );
   READ_DOUBLE(problem, *Re    );
   READ_DOUBLE(problem, *GX );
   READ_DOUBLE(problem, *GY );
   READ_DOUBLE(problem, *PI );
   READ_DOUBLE(problem, *UI );
   READ_DOUBLE(problem, *VI );
   //READ_STRING(problem,*geometry);


   *dx = (*xlength / (double)(*imax) );
   *dy = *ylength / (double)(*jmax);

   return 1;
}

void init_uvp(double UI,double VI,double PI,int imax,int jmax,double **U,double **V,double **P, int **Flag)
{
    int i;
    int j;

    for(i = 0 ; i <= imax ; i++)                    // Initialises values of fluid cells only
        for(j = 0 ; j <= jmax +1 ; j++)
            U[i][j] = UI;

    for(i = 0 ; i <= imax + 1 ; i++)                 // Initialises values of fluid cells only
        for(j = 0 ; j <= jmax ; j++)     
            V[i][j] = VI;

    for(i = 0 ; i <= imax + 1 ; i++)                 // Initialises values of fluid cells only
        for(j = 0 ; j <= jmax + 1 ; j++) 
            P[i][j] = PI;
}


void init_flag(const char* geometry, int imax, int jmax, int **Flag)
{

        int i;
        int j;
        int nofn = 0;      // no. of neighbors
        int **image_values = imatrix(0, imax + 1, 0 ,jmax + 1);
        image_values = read_pgm(geometry);

        for(i = 0 ; i <= imax + 1; i++)                        //Loop to set the Fluid or Boundary field
            for(j = 0 ; j <= jmax + 1; j++)
                switch(image_values[i][j])
                {
                    case 0:
                            Flag[i][j] = Flag[i][j] | 1;               //Fluid Cell
                            break;
                    case 1:
                            Flag[i][j] = Flag[i][j] | 2;               //No-slip 
                            break;
                    case 2:
                            Flag[i][j] = Flag[i][j] | 4;               //Free-slip
                            break;
                    case 3:
                            Flag[i][j] = Flag[i][j] | 8;               //Out-flow 
                            break;
                    case 4:
                            Flag[i][j] = Flag[i][j] | 16;              //In-flow
                            break;
                }

        for(i = 0 ; i <= imax + 1 ; i++)                      //Loop to set the Fluid Neighbors
            for(j = 0 ; j <= jmax + 1 ; j++)   
                if(!(Flag[i][j]&1))                           // Not a fluid cell 
                {
                    if((j!=jmax+1) && (Flag[i][j+1]&1) )               
                        { Flag[i][j] = Flag[i][j] | 32; nofn++;}          //North Neighbor
                                            
                    if(j!=0 && (Flag[i][j-1]&1) )
                        { Flag[i][j] = Flag[i][j] | 64; nofn++;}          //South Neighbor        
                   
                    if(i!=0 && (Flag[i-1][j]&1) )
                        { Flag[i][j] = Flag[i][j] | 128; nofn++;}         //West Neighbor
                  
                    if((i!=imax+1) && (Flag[i+1][j]&1) )
                        { Flag[i][j] = Flag[i][j] | 256; nofn++;}         //East Neighbor
                    
                    assert(nofn < 3);              // If number of fluid neighbors is greater than 2, sends error
                    nofn = 0;
                }
            
}
