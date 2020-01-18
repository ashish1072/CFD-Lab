#include "helper.h"
#include "visual.h"
#include "init.h"
#include"uvp.h"
#include"boundary_val.h"
#include"sor.h"
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "precice/SolverInterfaceC.h"
#include "precice_adapter.h"



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
int main(int argn, char** args)
{

	
					
					

    //define parameter variables
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
    double alpha;             /* uppwind differencing factor*/
    double omg;               /* relaxation factor */
    double tau;               /* safety factor for time step*/
    int  itermax;             /* max. number of iterations  */
    				          /* for pressure per time step */
    double eps;               /* accuracy bound for pressure*/
    double dt_value;          /* time for output */
    double Pr;                /* Prandtl number */
    double TI;                /*Initial value of Temp*/
    double T_h;
    double T_c;
    double beta;
    int choice;
		char* geometry = (char*)(malloc(sizeof(char)*300));
	char* problem = (char*)(malloc(sizeof(char)*200));

	 char*  precice_config= (char*)(malloc(sizeof(char)*400));
    char* participant_name = (char*)(malloc(sizeof(char)*400));
    char* mesh_name = (char*)(malloc(sizeof(char)*400));
    char* read_data_name = (char*)(malloc(sizeof(char)*400));
    char* write_data_name = (char*)(malloc(sizeof(char)*400));
	double x_origin;
	double y_origin;
 
    double **P = NULL;
    double **U = NULL;
    double **V = NULL;
    double **F = NULL;
    double **G = NULL;
    double **RS = NULL;
    int **flag = NULL;
    double **T = NULL;
    double **T1 = NULL;
    double t = 0;             // iteration variables
    int n = 0;
    int visual_iter = 0;
    int it;                   // Sor iterable
	double res;               // residual
    struct stat st = {0};     // Folders for vtk files
	char sol_folder[80];
	char sol_directory[80];
    int sor_break = 0;        // number of times sor does not converge
     int num_of_coupling_cells;
    int num_of_fluid_cells;
	double U_INFLOW;
	double V_INFLOW;
	double T_INFLOW;
	int i;
    int j;
	double Th,Tc,ks,L;

    printf("Select the problem from the list below by typing 1-6 \n");      // Menu
	printf("1. Heated Plate \n");
	printf("2. Convection \n");
	printf("3. Heat Exchanger F1 \n");
	printf("4. Heat Exchanger F2\n" );
	scanf("%d",&choice);
	printf("setting the choice\n");
	//int dummy = 0;
	int xhstart;
	int xhend ;
	int yhstart ;
	int yhend ;
	int xcstart ;
	int xcend;
	int ycstart;
	int ycend;
	double time_cp;
	const char* filename = "0";
	switch(choice)
	{
		case 1:
		        filename = "configs/heated-plate.dat";
				break;
		case 2:
				filename = "configs/convection.dat";
				break;
		case 3:
			    filename = "configs/F1-heat-exchange.dat";
				break;
		case 4:
			    filename = "configs/F2-heat-exchange.dat";
				break;
    }
    //Read and assign the parameter values from data file
      read_parameters(filename, problem, geometry, &imax, &jmax, &xlength, &ylength,
			&dt, &t_end, &tau, &dt_value, &eps, &omg, &alpha, &itermax,
			&GX, &GY, &Re, &Pr,&U_INFLOW,&V_INFLOW, &T_INFLOW, &Th, &Tc, &ks, &L, &UI, &VI, &PI, &TI, &T_h, &T_c, &beta, &dx, &dy, &x_origin, &y_origin, precice_config, participant_name, mesh_name, read_data_name, write_data_name);
	    int include_temp = 1;
		
switch(choice)
	{
		case 3:

				xhstart = 0;
				xhend = 0;
				yhstart = 0;
				yhend = jmax-1;
				xcstart= imax-1;
				xcend = imax-1;
				ycstart = 0;
				ycend = jmax-1;
				break;
		case 4:
				xhstart = 0;
				xhend = 0;
				yhstart = 0;
				yhend = jmax-1;
				xcstart= imax-1;
				xcend = imax-1;
				ycstart = 0;
				ycend = jmax-1;
				break;
        case 5:
				
				xhstart = 0;
				xhend = 0;
				yhstart = 0;
				yhend = jmax-1;
				xcstart= imax-1;
				xcend = imax-1;
				ycstart = 0;
				ycend = jmax-1;
				break;
		case 6:

				xhstart = 0;
				xhend = imax-1;
				yhstart = 0;
				yhend = 0;
				xcstart = 0;
				xcend = imax-1;
				ycstart= jmax-1;
				ycend = jmax-1;
				break;
    }

    //Allocate the matrices for P(pressure), U(velocity_x), V(velocity_y), F, and G on heap
    printf("Status: Starting matrix allocation... \n");
    P = matrix(0, imax-1, 0, jmax-1);
    U = matrix(0, imax-1, 0, jmax-1);
    V = matrix(0, imax-1, 0, jmax-1);
    F = matrix(0, imax-1, 0, jmax-1);
    G = matrix(0, imax-1, 0, jmax-1);
    RS = matrix(0, imax-1, 0, jmax-1);
    flag = imatrix(0, imax-1, 0, jmax-1);
	double **T_set_h = matrix(0, imax-1, 0, jmax-1);
	double **T_set_c = matrix(0, imax-1, 0, jmax-1);
    double **T_cp = matrix(0, imax-1, 0, jmax-1);
	double **U_cp = matrix(0, imax-1, 0, jmax-1);
	double **V_cp = matrix(0, imax-1, 0, jmax-1);
	if(include_temp)
	{
		T = matrix(0, imax-1, 0, jmax-1);
		T1= matrix(0, imax-1, 0, jmax-1);
	}
    const char* coric = precicec_actionReadIterationCheckpoint();
	const char* cowic = precicec_actionWriteIterationCheckpoint();
    //Initilize flags
    init_flag(problem,geometry, imax, jmax, flag, &num_of_fluid_cells, &num_of_coupling_cells);
    //Initialize the U, V, P, in case T
    if(include_temp){
		init_uvpt(UI, VI, PI, TI, imax, jmax, U, V, P, T, flag);
		}
	else{
		init_uvp(UI, VI, PI, imax, jmax, U, V, P, flag);
	}
	//Make solution folder
	sprintf( sol_folder,"Solution_%s",problem);
	if (stat(sol_folder, &st) == -1) 
    		mkdir(sol_folder, 0700);
	
	sprintf( sol_directory,"Solution_%s/sol", problem);

    printf("Status: Starting the simulation...\n");
	precicec_createSolverInterface(participant_name,precice_config,0,1);
	int dim = precicec_getDimensions();
	int meshID = precicec_getMeshID(mesh_name);
	int* vertexIDs = precice_set_interface_vertices(imax,jmax,dx,dy,x_origin,y_origin,num_of_coupling_cells,meshID,flag);
	int temperatureID = precicec_getDataID(write_data_name,meshID);
	int heatFluxID = precicec_getDataID(read_data_name,meshID);
	double* heatfluxCoupled = (double*)malloc(num_of_coupling_cells*sizeof(double));
	double precice_dt = precicec_initialize();
	double* temperature = (double*)malloc(num_of_coupling_cells*sizeof(double));
	precice_write_temperature(imax,jmax,L, ks, Th, Tc,num_of_coupling_cells,temperature,vertexIDs,temperatureID,T,flag);
	precicec_initialize_data();
	precicec_readBlockScalarData(heatFluxID,num_of_coupling_cells,vertexIDs,heatfluxCoupled);
	while (precicec_isCouplingOngoing())  
	{
       	if(precicec_isActionRequired(cowic)){
			write_checkpoint(t, U, V, T, &time_cp, U_cp, V_cp, T_cp, imax, jmax);
			precicec_fulfilledAction(cowic);
}
		calculate_dt(Re,tau,&dt,dx,dy,imax,jmax, U, V, Pr, include_temp);
   		printf("t = %f ,dt = %f, \n",t,dt);
		dt = fmin(dt,precice_dt);
		boundaryvalues(imax, jmax, U, V, flag);
		if(include_temp){
		//	if (choice == 3 || choice ==4 || choice ==5 || choice ==5){ /// WS 2 cases turn off id not required
		//	set_dirichlet(xhstart,xhend,yhstart,yhend,T_set_h,T,flag);
	//		set_dirichlet(xcstart,xcend,ycstart,ycend,T_set_c,T,flag);
		//	}
			set_temp_boundary(imax, jmax,  T_INFLOW, T, flag);
			set_coupling_boundary(imax,jmax,L, ks, Th, Tc,dx,dy,heatfluxCoupled,T,flag);
			calculate_temp(T, U, V, flag, Re, Pr, dt, dx, dy, alpha, imax, jmax);
			}
    	spec_boundary_val(imax, jmax, U_INFLOW,V_INFLOW,U, V, flag);

    	calculate_fg(Re,GX,GY,alpha,dt,dx,dy,imax,jmax,U,V,F,G,flag, beta, T, include_temp);

    	calculate_rs(dt,dx,dy,imax,jmax,F,G,RS,flag);
		it = 0;
		res = 10.0;

    	while( it < itermax && res > eps)
    	{
    		sor(omg,dx,dy,imax,jmax,P,RS,&res,flag);
			++it;
    	} 

		printf("SOR itertions = %d \n", it-1);

		if(it == itermax)
			sor_break++;
		calculate_uv(dt,dx,dy,imax,jmax,U,V,F,G,P,flag);
		precice_write_temperature(imax,jmax,L, ks, Th, Tc,num_of_coupling_cells,temperature,vertexIDs,temperatureID,T,flag);
		precice_dt = precicec_advance(dt);
		
		if(precicec_isActionRequired(coric)){ // timestep not converged
			restore_checkpoint(&t, U, V, T, time_cp, U_cp, V_cp, T_cp, imax, jmax);
			precicec_fulfilledAction(coric);
		}
		else{ // timestep converged
			t =t+ dt;
			n = n+ 1;
}		
		precicec_readBlockScalarData(heatFluxID,num_of_coupling_cells,vertexIDs,heatfluxCoupled);

/*		if(!include_temp)
			noslip(U, V, P, flag, imax, jmax);
		else
			noslip2(U, V, P, T, flag, imax, jmax);
 */
		if ((t >= visual_iter*dt_value) && (t!=0.0))
  		{
   			write_vtkFile(sol_directory ,n ,xlength ,ylength ,imax-2 ,jmax-2 ,
							dx ,dy ,U ,V ,P,T,include_temp,x_origin,y_origin);

			printf("Result at %f seconds \n",visual_iter*dt_value);
    		visual_iter++;
    		continue;
  		}
    	
		
    }
	precicec_finalize();
	
    

    //Free memory
    free_matrix( P, 0, imax-1, 0, jmax-1);
    free_matrix( U, 0, imax-1, 0, jmax-1);
    free_matrix( V, 0, imax-1, 0, jmax-1);
    free_matrix( F, 0, imax-1, 0, jmax-1);
    free_matrix( G, 0, imax-1, 0, jmax-1);
    free_matrix(RS, 0, imax-1, 0, jmax-1);
    free_imatrix(flag, 0, imax-1, 0, jmax-1);
	free(heatfluxCoupled);
	free(precice_config);
    free(mesh_name);
    free(read_data_name);
    free(write_data_name);
    free(participant_name);
	free(temperature);
	if(include_temp) { free_matrix(T, 0, imax-1, 0, jmax-1);
			   free_matrix(T1, 0, imax-1, 0, jmax-1); }
	free(geometry);
	free(problem);
    printf("Status: allocated memory released...\n \n");
    printf("WARNING: SOR did not converge %d times out of %d iterations. \n",sor_break,n-1);
	printf("Status: Program end\n");
  return -1;

}
