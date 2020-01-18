#include "helper.h"
#include "init.h"
#include <stdio.h>
#include<assert.h>

void read_parameters( const char *szFileName,       /* name of the file */
			char *problem,
			char *geometry,
		    int  *imax,                /* number of cells x-direction*/
            int  *jmax,                /* number of cells y-direction*/ 
		    double *xlength,           /* length of the domain x-dir.*/
            double *ylength,           /* length of the domain y-dir.*/
		    double *dt,                /* time step */
		    double *t_end,             /* end time */
		    double *tau,               /* safety factor for time step*/
		    double *dt_value,	       /* time for output */
		    double *eps,               /* accuracy bound for pressure*/
		    double *omg,               /* relaxation factor */
		    double *alpha,             /* uppwind differencing factor*/
            int  *itermax,             /* max. number of iterations  */
		    double *GX,                /* gravitation x-direction */
            double *GY,                /* gravitation y-direction */
		    double *Re,                /* reynolds number   */
			double *Pr,
			double *U_INFLOW,
			double *V_INFLOW,
			double *T_INFLOW,
			double *Th, 
			double *Tc,
			double *ks, 
			double *L,
		    double *UI,                /* velocity x-direction */
            double *VI,                /* velocity y-direction */
            double *PI,                /* pressure */
       		double *TI,
		    double *T_h,
		    double *T_c,
		    double *beta,
		    double *dx,                /* length of a cell x-dir. */
            double *dy                 /* length of a cell y-dir. */,
		    double *x_origin,
		    double *y_origin,
		    char* precice_config,
            char* participant_name,
		    char* mesh_name,
		    char* read_data_name,
		    char* write_data_name)           
{
   
   //READ_STRING( szFileName, *problem );
   //READ_STRING( szFileName, geometry );

   READ_INT   ( szFileName, *imax );
   READ_INT   ( szFileName, *jmax );

   READ_DOUBLE( szFileName, *xlength );
   READ_DOUBLE( szFileName, *ylength );

   READ_DOUBLE( szFileName, *dt    );
   READ_DOUBLE( szFileName, *t_end );
   READ_DOUBLE( szFileName, *tau   );
   READ_DOUBLE( szFileName, *dt_value );
	READ_DOUBLE( szFileName, *U_INFLOW   );
	READ_DOUBLE( szFileName, *V_INFLOW   );
   READ_DOUBLE( szFileName, *eps   );
   READ_DOUBLE( szFileName, *omg   );
   READ_DOUBLE( szFileName, *alpha );
   READ_DOUBLE( szFileName, *T_INFLOW );
   READ_INT   ( szFileName, *itermax );
	READ_DOUBLE( szFileName, *Th );
	READ_DOUBLE( szFileName, *Tc );
	READ_DOUBLE( szFileName, *ks );
	READ_DOUBLE( szFileName, *L );
   READ_DOUBLE( szFileName, *GX );
   READ_DOUBLE( szFileName, *GY );
   READ_DOUBLE( szFileName, *Re );
   READ_DOUBLE( szFileName, *Pr );

   READ_DOUBLE( szFileName, *UI );
   READ_DOUBLE( szFileName, *VI );
   READ_DOUBLE( szFileName, *PI );
   READ_DOUBLE( szFileName, *TI );
   READ_DOUBLE( szFileName, *T_h );
   READ_DOUBLE( szFileName, *T_c );
   READ_DOUBLE( szFileName, *beta );
	
   READ_STRING( szFileName, problem);
   READ_STRING( szFileName, geometry);

   READ_DOUBLE( szFileName, *x_origin );
   READ_DOUBLE( szFileName, *y_origin );
   READ_STRING( szFileName, precice_config);
   READ_STRING( szFileName, participant_name);
   READ_STRING( szFileName, mesh_name);
   READ_STRING( szFileName, read_data_name);
   READ_STRING( szFileName, write_data_name);

   *dx = *xlength / (double)(*imax);
   *dy = *ylength / (double)(*jmax);

   

}


void init_uvp(double UI, double VI, double PI, int imax, int jmax,
		 double** U, double** V, double** P, int** flag)
{
	printf("Status: initialization of U,V,P \n");
	
	for(int i=0; i<imax; i++){
		for(int j=0; j<jmax; j++){
			if(flag[i][j]&(1<<0)){

				U[i][j] = UI;
				V[i][j] = VI;
				P[i][j] = PI;
			}
		}
	}
	printf("Status: U,V,P matrices initialized... \n \n");
}

void init_uvpt(double UI, double VI, double PI, double TI, int imax, int jmax,
		 double** U, double** V, double** P, double** T, int** flag)
{
	printf("Status: Starting initialization of U,V,P,T ... \n");
	
	for(int i=0; i<imax; i++){
		for(int j=0; j<jmax; j++){
			if(flag[i][j]&(1<<0)){

				U[i][j] = UI;
				V[i][j] = VI;
				P[i][j] = PI;
				T[i][j] = TI;
			}
		}
	}
	printf("Status: U,V,P,T matrices initialized... \n \n");
}

int  isfluid(int pic){
	if((pic == 2)||(pic == 3)||(pic == 6)) {return 1;}
		else {return 0;}
}


void init_flag(char* problem, char* geometry, int imax, int jmax, int **flag, int* num_of_fluid_cells, int* num_of_coupling_cells)
{
	int fluid = 0;
        int coupling = 0;
	int neighbors;
	int i;
	int j;
	int **pic = imatrix(0,imax-1,0,jmax-1);
	pic = read_pgm(geometry);
	
	for (int i=0; i<imax; i++)
	{
		for (int j=0; j<jmax; j++)
		{		flag[i][j] = 0;
		switch(pic[i][j])
		{
			case 0: //noslip
			flag[i][j] = 1<<1;
			break;

			case 1: //freeslip
			flag[i][j] = 1<<2;
			break;

			case 2: //outflow
			flag[i][j] = 1<<3;
			break;

			case 3: //inflow
			flag[i][j] = 1<<4;
			break;

			case 4: //coupling
			flag[i][j] = 1<<9;
			coupling++;
			break;
			
			case 6://fluid
			flag[i][j] = 1<<0;
			fluid++;
			break;
			
		}
		neighbors = 0;
			if(!isfluid(pic[i][j])) //set boundaries if not 
			{
		//		printf("setting neighbours\n");
				if(i<imax-1 && pic[i+1][j]==6)
				{
				flag[i][j] |= 1<<8;
				neighbors++;
				}
				if( i>0 && pic[i-1][j]==6)
				{
				flag[i][j] |= 1<<7;
				neighbors++;
				}
				if(j<jmax-1 && pic[i][j+1]==6)
				{
				flag[i][j] |= 1<<5;
				neighbors++;
				}
				if(j>0 && pic[i][j-1]==6)
				{
				flag[i][j] |= 1<<6;
				neighbors++;
				}
			}
		}

	}

	*num_of_fluid_cells = fluid;
	*num_of_coupling_cells = coupling;
	free_imatrix(pic, 0,imax-1,0,jmax-1);
	


}
