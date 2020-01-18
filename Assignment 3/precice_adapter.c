#include "precice_adapter.h"
#include "precice/SolverInterfaceC.h"
#include <stdio.h>
#include <stdlib.h>

int *precice_set_interface_vertices(int imax, int jmax, double dx, double dy, double x_origin, double y_origin,
                                    int num_coupling_cells, int meshID, int **Flag)
{
	int dimension = precicec_getDimensions();
	double *vertices = (double*) malloc(num_coupling_cells*512*sizeof(double));
	int i;
	int j;
	int index = 0;
  int *ids=(int*) malloc(num_coupling_cells*sizeof(int));
      //bottom wall
    for(int i = 0; i<imax; i++)
        if (Flag[i][0]&(1<<9)) {
            vertices[index]= x_origin + ((double)i-0.5)*dx;
            index++;
            vertices[index] = y_origin;
            index++;
            vertices[index] = 0;
            index;
        }
       
    //top wall
    for(int i = 0; i<imax; i++)
        if (Flag[i][jmax-1]&(1<<9)) {
            vertices[index]= x_origin + ((double)i-0.5)*dx;
            index++;
            vertices[index] = y_origin + ((double)jmax-1)*dy;
            index++;
            vertices[index] = 0;
            index++;
        }
        
    //left wall
    for(int j = 1 ; j< jmax-1; j++)
        if (Flag[0][j]&(1<<9)) {
            vertices[index]= x_origin;
            index++;
            vertices[index] = y_origin + ((double)j-0.5)*dy;
            index++;
            vertices[index] = 0;
            index++;
        }
        
    //right wall
    for(int j = 1 ; j< jmax-1; j++)
        if (Flag[imax-1][j]&(1<<9)) {
            vertices[index]= x_origin + ((double)imax-1)*dx;
            index++;
            vertices[index] = y_origin + ((double)j-0.5)*dy;
            index++;
            vertices[index] = 0;
            index++;
}

	for(i = 1 ; i < imax-1 ; i++){
		for(j = 1 ; j < jmax-1 ; j++){
        if(Flag[i][j]&(1<<9) && (Flag[i][j+1]&1<<0)){ // north fluid cell
          
            vertices[index] = x_origin + ((double)i-0.5)*dx;
   //         printf("vertices[%d]=%f\n",index,vertices[index]);
            index++;
            vertices[index] = y_origin + ((double)j-1.0)*dy;
     //       printf("vertices[%d]=%f\n",index,vertices[index]);
            index++;
            vertices[index] = 0.0;
      //      printf("vertices[%d]=%f\n",index,vertices[index]);
            index++;
        }
          if(Flag[i][j]&(1<<9) && (Flag[i][j-1]&1<<0)){  // South Fluid Cell
            vertices[index] = x_origin + ((double)i -0.5)*dx;
            index++;
            vertices[index] = y_origin + ((double)j+1.0)*dy;
            index++;
            vertices[index] = 0.0;
            index++; 
        }



		}
    }
		precicec_setMeshVertices(meshID,num_coupling_cells,vertices,ids);
    free(vertices);
    return ids;			
}





void precice_write_temperature(int imax, int jmax, double L, double ks, double Th, double Tc,int num_coupling_cells, double *temperature, int *vertexIDs,
                               int temperatureID, double **TEMP, int **Flag)
{
    int i,j;
    int n = 0;
    double Tinf = Tc;    
    //moving over the walls
    //bottom wall
    for( i = 0; i<imax; i++)
        if (Flag[i][0]&(1<<9)) {
            temperature[n] = TEMP[i][1]*(Th-Tc) + Tinf;
  //          printf("Writing Temp Flag[i][j]>>5 = %d, temp[%d]=%f\n", Flag[i][0]>>5,n,temperature[n]);
            n++;
            
        }
    //top wall
    for(i = 0; i<imax; i++)
        if (Flag[i][jmax-1]&(1<<9)) {
            temperature[n] = TEMP[i][jmax-2]*(Th-Tc) + Tinf;
  //          printf("Writing Temp Flag[i][j]>>5 = %d, temp[%d]=%f\n", Flag[i][jmax-1]>>5,n,temperature[n]);
            n++;
        }
    //left wall
    for(int j = 1 ; j< jmax-1; j++)
        if (Flag[0][j]&(1<<9)) {
            temperature[n] = TEMP[1][j]*(Th-Tc) + Tinf;
//  printf("Writing Temp Flag[i][j]>>5 = %d, temp[%d]=%f\n", Flag[0][j]>>5,n,temperature[n]);
            n++;
//            printf("Writing Temp Flag[i][j]>>5 = %d, temp[%d]=%f\n", Flag[i][j]>>5,n,temperature[n]);
        }
    //right wall
    for(int j = 1 ; j< jmax-1; j++)
        if (Flag[imax-1][j]&(1<<9)) {
            temperature[n] = TEMP[imax-2][j]*(Th-Tc) + Tinf;
 //  printf("Writing Temp Flag[i][j]>>5 = %d, temp[%d]=%f\n", Flag[imax-1][j],n,temperature[n]);
            n++;
 //           printf("Writing Temp Flag[i][j]>>5 = %d, temp[%d]=%f\n", Flag[i][j]>>5,n,temperature[n]);
        }

    //inside the domain
    for(int j = 1 ; j< jmax-1; j++)
        for(int i = 1; i<imax-1; i++)
        {
            //if bottom is fluid
            if ( (Flag[i][j]&(1<<9))&&(Flag[i][j-1]&1<<0) ) {
                temperature[n] = TEMP[i][j-1]*(Th-Tc) + Tinf;
  //  printf("Writing Temp Flag[i][j] = %d, temp[%d]=%f\n", Flag[i][j],n,temperature[n]);
                n++;
  //              printf("Writing Temp Flag[i][j]>>5 = %d, temp[%d]=%f\n", Flag[i][j]>>5,n,temperature[n]);
            }
            //if top is fluid
            if ( (Flag[i][j]&(1<<9))&&(Flag[i][j+1]&1<<0) ) {
                temperature[n] = TEMP[i][j+1]*(Th-Tc) + Tinf;
   //printf("Writing Temp Flag[i][j]>>5 = %d, temp[%d]=%f\n", Flag[i][j]>>5,n,temperature[n]);
                n++;
 //               printf("Writing Temp Flag[i][j]>>5 = %d, temp[%d]=%f\n", Flag[i][j]>>5,n,temperature[n]);
            }
        }
precicec_writeBlockScalarData(temperatureID, num_coupling_cells, vertexIDs, temperature);
    
    
    
    
    /*
    for(i = 0 ; i < imax ; i++){
		for(j = 0 ; j < jmax ; j++){
			if(Flag[i][j]&(1<<9))
			{
        if (choice == 1 || choice == 2){
        switch(Flag[i][j]>>5){
          case 24: // east fluid cell and self coupling
            temperature[n] = TEMP[i+1][j]*(Th-Tc) + Tinf;
            n++;
            break;
          
          case 20:  /// West Fluid cell
            temperature[n] = TEMP[i-1][j]*(Th-Tc) + Tinf;
            n++;
            break;
          
          case 17:  // north fluid cell
            temperature[n] = TEMP[i][j+1]*(Th-Tc) + Tinf;
            printf("Writing Temp Flag[i][j]>>5 = %d, temp[%d]=%f\n", Flag[i][j]>>5,n,temperature[n]);
            n++;
            break;
            
          case 18:  // South Fluid Cell
            temperature[n] = TEMP[i][j-1]*(Th-Tc) + Tinf;
            n++;
            break;

        }
        }
        else{

        }
        }
        } */
      }









 // Left Boundary
/*      for( j= 0; j<= jmax-1; j++)
     {
         if(FLAG[0][j] & 512)                // If cell is a coupling cell
           { temperature[n] = TEMP[1][j];
             n++;}
     }

// Right Boundary
     for( j= 0; j<=jmax-1; j++)
     {
         if(FLAG[imax-1][j] & 512)               
           { temperature[n] = TEMP[imax-2][j];
             n++;}
     }

// Top boundary
     for( i= 1; i<= imax-2; i++)
     {
         if(FLAG[i][jmax-1] & 512)               
           { temperature[n] = TEMP[i][jmax-2];
             n++;}
     }
    
// Bottom boundary
    for( i= 1; i<= imax-2; i++)
     {
         if(FLAG[i][0] & 512)               
           { temperature[n] = TEMP[i][1];
             n++;}
     }
*/

void set_coupling_boundary(int imax, int jmax,double L, double ks, double Th, double Tc, double dx, double dy, double *heatflux, double **TEMP, int **Flag)
{
  
    int i,j;
    double dT=Th-Tc;
    int n = 0;
    //moving over the walls
    //bottom wall
    for(int i = 0; i<imax; i++)
        if (Flag[i][0]&(1<<9)) {
            TEMP[i][0] = TEMP[i][1] + dy*L*heatflux[n]/(ks*dT);
            n++;
        }
    //top wall
    for(int i = 0; i<imax; i++)
        if (Flag[i][jmax-1]&(1<<9)) {
            TEMP[i][jmax-1] = TEMP[i][jmax-2] + dy*heatflux[n]*L/(ks*dT);
            n++;
        }
    //left wall
    for(int j = 1 ; j< jmax-1; j++)
        if (Flag[0][j]&(1<<9)) {
            TEMP[0][j] = TEMP[1][j] + dx*heatflux[n]*L/(ks*dT);
            n++;
        }
    //right wall
    for(int j = 1 ; j< jmax-1; j++)
        if (Flag[imax-1][j]&(1<<9)) {
            TEMP[imax-1][j] = TEMP[imax-2][j] + dx*heatflux[n]*L/(ks*dT);
            n++;
        }

    //inside the domain
    for(int j = 1 ; j< jmax-1; j++)
        for(int i = 1; i<imax-1; i++)
        {
            //if bottom is fluid
            if ( (Flag[i][j]&(1<<9))&&(Flag[i][j-1]&1<<0) ) {
                TEMP[i][j] = TEMP[i][j-1] + dy*heatflux[n]*L/(ks*dT);
                n++;
            }
            //if top is fluid
            if ( (Flag[i][j]&(1<<9))&&(Flag[i][j+1]&1<<0) ) {
                TEMP[i][j] = TEMP[i][j+1] + dy*heatflux[n]*L/(ks*dT);
                n++;
            }
}

  
  
  
  
  
  
  
  
  /* 
  for(i = 0 ; i < imax ; i++){
		for(j = 0 ; j < jmax ; j++){
			if(Flag[i][j]&(1<<9))
			{
   //     printf("Reading Flux Flag[i][j]>>5 = %d\n", Flag[i][j]>>5);
     //   printf("Writing Temp Flag[i][j]>>5 = %d", Flag[i][j]>>5);
        switch(Flag[i][j]>>5){
          case 24: // east fluid cell and self coupling
            TEMP[i][j] = TEMP[i+1][j]  + dx*L*heatflux[n]/(ks*dT);
             n++;
            break;
          
          case 20:  /// West Fluid cell
            TEMP[i][j] = TEMP[i-1][j]  + dx*L*heatflux[n]/(ks*dT);
             n++;
             break;
          
          case 17:  // north fluid cell
            TEMP[i][j] = TEMP[i][j+1]  + dx*L*heatflux[n]/(ks*dT);
            printf("Writing Heatflux Flag[i][j]>>5 = %d, hf[%d]=%f, TEMP[%d][%d]=%f\n", Flag[i][j]>>5,n,heatflux[n],i,j,TEMP[i][j]);
            n++;
            break;
            
          case 18:  // South Fluid Cell
            TEMP[i][j] = TEMP[i][j-1]  + dx*L*heatflux[n]/(ks*dT);
             n++;
            break;

        }
        }
        }
      }

  */

 // Left Boundary
/*     for( j= 0; j<= jmax-1; j++)
     {
         if(FLAG[0][j] & 512)                       // If cell is a coupling cell
           { TEMP[0][j] = TEMP[1][j]  + dx*heatflux[n];
             n++;}
     }

// Right Boundary
     for( j= 0; j<= jmax-1; j++)
     {
         if(FLAG[imax-1][j] & 512)               
           { TEMP[imax-1][j] = TEMP[imax-2][j]  + dx*heatflux[n];
             n++;}
     }

// Top boundary
     for( i= 1; i<= imax-2; i++)
     {
         if(FLAG[i][jmax-1] & 512)               
           { TEMP[i][jmax-1] = TEMP[i][jmax-2] + dy*heatflux[n];
             n++;}
     }
    
// Bottom boundary
    for( i= 1; i<= imax-2; i++)
     {
         if(FLAG[i][0] & 512)               
           { TEMP[i][0] = TEMP[i][1] + dy*heatflux[n];
             n++;}
     }
 */
}

void write_checkpoint(double time, 
                      double** U, 
                      double **V, 
                      double** T, 
                      double* time_cp, 
                      double** U_cp, 
                      double** V_cp, 
                      double **T_cp, 
                      int imax, 
                      int jmax
                      ){
                      for (int i=0;i<imax;++i){
                        for (int j=0;j<jmax;++j){
                          U_cp[i][j]=U[i][j];
                          V_cp[i][j] = V[i][j];
                          T_cp[i][j] = T[i][j];
                          *time_cp = time;
                        }
                      }
}

void restore_checkpoint(double *time, 
                        double **U, 
                        double **V, 
                        double **TEMP, 
                        double time_cp, 
                        double **U_cp, 
                        double **V_cp, 
                        double **TEMP_cp, 
                        int imax, 
                        int jmax)
{
int i;
int j;

for(i = 0 ; i < imax ; i++ ){
  for(j = 0 ; j < jmax ; j++)
  {
    U[i][j] = U_cp[i][j];
    V[i][j] = V_cp[i][j];
    TEMP[i][j] = TEMP_cp[i][j];
    *time = time_cp;
    }
  }
}

