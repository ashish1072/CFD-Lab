#include "parallel.h"


void Program_Message(char *txt)
/* produces a stderr text output  */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
}


void Programm_Sync(char *txt)
/* produces a stderr textoutput and synchronize all processes */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                             /* synchronize output */  
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
}


void Programm_Stop(char *txt)
/* all processes will produce a text output, be synchronized and finished */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                           /* synchronize output */
   fprintf(stderr,"-STOP- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   exit(1);
}

/* Initialization of all subdomains */
void init_parallel(int iproc, int jproc, int imax, int jmax, int *myrank,
 int *il, int *ir, int *jb, int *jt, int *rank_l, int *rank_r, int *rank_b, int *rank_t,
  int *omg_i, int *omg_j,double *dx, double *dy, int xlength, int ylength)
{
   *omg_i = ((*myrank)%iproc) + 1;
   *omg_j = ((*myrank)/iproc) + 1;

   *il = (*omg_i - 1)*(imax/iproc) + 1;
   *ir = (*omg_i)*(imax/iproc);

   *jb = (*omg_j - 1)*(jmax/jproc) + 1;
   *jt = (*omg_j)*(jmax/jproc);

   if(*il == 1)
      *rank_l = MPI_PROC_NULL;
   else
      *rank_l = *myrank - 1;

   if(*ir == imax)
      *rank_r = MPI_PROC_NULL;
   else
      *rank_r = *myrank + 1;
   
   if(*jb == 1)
      *rank_b = MPI_PROC_NULL;
   else
      *rank_b = *myrank - 2;
   
   if(*jt == jmax)
      *rank_t = MPI_PROC_NULL;
   else
      *rank_t = *myrank + 2;
   
  // *dx = xlength/ (double)(*ir-*il+1);
  // *dy = ylength/ (double)(*jt-*jb+1);

   MPI_Barrier(MPI_COMM_WORLD);

   printf("Thread_id: %d omg_ij: %d%d \nil: %d, ir: %d, jb: %d, jt: %d \n",*myrank,*omg_i,*omg_j, *il,*ir,*jb,*jt);

   printf("my rank %d, l_rank: %d, r_rank: %d, b_rank: %d, t_rank: %d \n \n", *myrank,*rank_l,*rank_r,*rank_b,*rank_t);
}



/* Exchange of pressure values at subdomain boundaries */
void pressure_comm(double **P,int il,int ir,int jb,int jt,
                   int rank_l, int rank_r, int rank_b, int rank_t, 
                   double *bufSend, double *bufRecv, MPI_Status *status, int chunk)
{ 

   // Sending pressure values to the left subdomain
         if (rank_l != MPI_PROC_NULL)
         { for (int j = 0; j< (jt-jb +1); j++ )
               {
                  bufSend[j] = P[1][j+1];
               }
                  MPI_Send( bufSend, (jt-jb +5), MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD );
         }

         if (rank_r != MPI_PROC_NULL)
         {  MPI_Recv( bufRecv, (jt-jb +5), MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD, status);    //****MPI_Status *status
            
            for (int j = 0; j< (jt-jb +1); j++ )
               {
                  P[ir-il+2][j+1] = bufRecv[j];
               }
         }

   
   //Sending pressure values to the right subdomain
         if (rank_r != MPI_PROC_NULL)
         { for (int j = 0; j< (jt-jb +1); j++ )
               {
                  bufSend[j] = P[ir-il+1][j+1];
               }
                  MPI_Send( bufSend, (jt-jb +5), MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD );
         }

         if (rank_l != MPI_PROC_NULL)
         {  MPI_Recv( bufRecv, (jt-jb +5), MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD, status);
            
            for (int j = 0; j< (jt-jb +1); j++ )
               {
                  P[0][j+1] = bufRecv[j];
               }
         }

   // Sending pressure values to the top subdomain
         if (rank_t != MPI_PROC_NULL)
         { for (int i = 0; i< (ir-il +1); i++ )
               {
                  bufSend[i] = P[i+1][jt-jb +1];
               }
                  MPI_Send( bufSend, (ir-il +5), MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD );
         }

         if (rank_b != MPI_PROC_NULL)
         {  MPI_Recv( bufRecv, (ir-il +5), MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD, status);
            
            for (int i = 0; i< (ir-il +1); i++ )
               {
                  P[i+1][0] = bufRecv[i];
               }
         }

    // Sending pressure values to the bottom subdomain
         if (rank_b != MPI_PROC_NULL)
         { for (int i = 0; i< (ir-il +1); i++ )
               {
                  bufSend[i] = P[i+1][1];
               }
                  MPI_Send( bufSend, (ir-il +5), MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD );
         }

         if (rank_t != MPI_PROC_NULL)
         {  MPI_Recv( bufRecv, (ir-il +5), MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD, status);
            
            for (int i = 0; i< (ir-il +1); i++ )
               {
                  P[i+1][jt-jb +2] = bufRecv[i];
               }
         }

}



void uv_comm(double** U, 
             double** V, 
             int il, 
             int ir, 
             int jb, 
             int jt, 
             int rank_l,
             int rank_r,
             int rank_b,
             int rank_t,
             double* bufSend,
             double* bufRecv,
             MPI_Status *status,
             int chunk
             ){
                // right-left-top-bottom
                int myrank;
                MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
                //printf("Communicating UV\n");
                if (rank_r != MPI_PROC_NULL){  //right
                   for (int i=1;i<=jt-jb+1;i++){  //for U                   
                     bufSend[i-1] = U[ir-1-(il-2)][i];   ////// this need rechecking
                   }
                   MPI_Send(bufSend, jt-jb+2, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD);
                }
                  if (rank_l!=MPI_PROC_NULL){ // if left is sending
                     MPI_Recv(bufRecv, jt-jb+2, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD, status);
                     for (int i=1;i<=jt-jb+1;i++){
                        U[0][i] = bufRecv[i-1];
                     }
                  }
            
                  // for V
                  if (rank_r != MPI_PROC_NULL){ 
                  for (int i=1;i<=jt-(jb-2);i++){  
                     bufSend[i-1] = V[ir-(il-1)][i];
                   }
                   MPI_Send(bufSend, jt-jb+2, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD);
                   }
                  if (rank_l!=MPI_PROC_NULL){ // if left is sending
                     MPI_Recv(bufRecv, jt-jb+2, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD, status);
                     for (int i=1;i<=jt-jb+2;i++){
                        V[0][i] = bufRecv[i-1];
                     }
                  }
             
               // left 
               if (rank_l != MPI_PROC_NULL){
                   for (int i=1;i<=jt-jb+1;i++){
                     bufSend[i-1] = U[2][i];
                   }
                   MPI_Send(bufSend, jt-jb+2, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD);
               }
                  if (rank_r!=MPI_PROC_NULL){
                     MPI_Recv(bufRecv, jt-jb+2, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD, status);
                     for (int i=1;i<=jt-jb+1;i++){
                        U[ir+1-(il-2)][i] = bufRecv[i-1];
                     }
                  }
          
                  // for V
                  if (rank_l != MPI_PROC_NULL){
                  for (int i=1;i<=jt-(jb-1)+1;i++){  
                     bufSend[i-1] = V[1][i];
                   }
                   MPI_Send(bufSend, jt-jb+2, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD);
                   }
                  if (rank_r!=MPI_PROC_NULL){ // if left is sending
                     MPI_Recv(bufRecv, jt-jb+2, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD, status);
                     for (int i=1;i<=jt-(jb-1)+1;i++){
                        V[ir+1-(il-1)][i] = bufRecv[i-1];
                     }
                  }
                

                // top 
               if (rank_t != MPI_PROC_NULL){
                   for (int i=1;i<=ir-(il-1)+1;i++){  //preparing to send
                     bufSend[i-1] = U[i][jt-(jb-1)];
                   }
                  
                   MPI_Send(bufSend, ir-il+2, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD);  //sending
               }
                  if (rank_b != MPI_PROC_NULL){  // if bottom sending receive
                  
                     MPI_Recv(bufRecv, ir-il+2, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD, status);
                     for (int i=1;i<=ir-(il-1)+1;i++){
                        U[i][0] = bufRecv[i-1];
                     }
                  }
            // for V
          if (rank_t != MPI_PROC_NULL){
                  for (int i=1;i<=ir-il+1;i++){ 
                     bufSend[i-1] = V[i][jt-1-(jb-2)];
                   }
                   
                   MPI_Send(bufSend, ir-il+2, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD);
                   }
                  if (rank_b!=MPI_PROC_NULL){ // if bottom is sending
                     
                     MPI_Recv(bufRecv, ir-il+2, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD, status);
                     for (int i=1;i<=ir-il+1;i++){
                        V[i][0] = bufRecv[i-1];
                     }
                  }
               



                // bottom
               if (rank_b != MPI_PROC_NULL){
                   for (int i=1;i<=ir-(il-1)+1;i++){  //preparing to send
                     bufSend[i-1] = U[i][1];
                   }
                   
                   MPI_Send(bufSend, ir-il+2, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD);  //sending
               }
                  if (rank_t != MPI_PROC_NULL){  // if top sending receive
                  
                     MPI_Recv(bufRecv, ir-il+2, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD, status);
                     for (int i=1;i<=ir-(il-1)+1;i++){
                        U[i][jt+1-(jb-1)] = bufRecv[i-1];
                     }
                  }
          // for V
             if (rank_b != MPI_PROC_NULL){
                  for (int i=1;i<=ir-il+1;i++){ 
                     bufSend[i-1] = V[i][2];
                   }
                   
                   MPI_Send(bufSend, ir-il+2, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD);
             }
                  if (rank_t!=MPI_PROC_NULL){ // if bottom is sending
                 
                     MPI_Recv(bufRecv, ir-il+2, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD, status);
                     for (int i=1;i<=ir-il+1;i++){
                        V[i][jt+1-(jb-2)] = bufRecv[i-1];
                     }
                  }
                

             }
