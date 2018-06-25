#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sys/time.h>
#include <mpi.h>

#define epsilon 1.e-8

using namespace std;

int main (int argc, char* argv[]){

  int M,N;

  string T,P,Db;
  M = atoi(argv[1]);
  N = atoi(argv[2]);

  if(argc > 3){

    T = argv[3];
    if(argc > 4){
      P = argv[4];
      if(argc > 5){
        Db = argv[5];
      }
    }
  }
 // cout<<T<<P<<endl;
  


  //Read from file matrix, if not available, app quit
  //Already transposed
  //

    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double *U_t;
    double *local_Alphas,*local_Gammas;
    double *Alphas, *Gammas;

    const int ROW_BLOCK = (M+size-1) / size;

    U_t = new double[N * N];
    local_Alphas = new double[ROW_BLOCK];
    local_Gammas = new double[ROW_BLOCK * N];

    if(rank == 0)
    {
        Alphas = new double[ROW_BLOCK * size];
        Gammas = new double[(ROW_BLOCK * size) * N];
    }



    if(rank == 0)
    {
      ifstream matrixfile("matrix");
      if(!(matrixfile.is_open())){
        cout<<"Error: file not found"<<endl;
        return 0;
      }

      for(int i = 0; i < M; i++){
        for(int j =0; j < N; j++){
          matrixfile >> U_t[i*N+j];
        }
      }
      matrixfile.close();
    }

	
	struct timeval ts, tt;

	gettimeofday(&ts, NULL);
    //broadcast to all process
    MPI_Bcast(U_t, N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* Reductions */


    int start = rank * ROW_BLOCK;
    int end = (rank < size-1) ? (start + ROW_BLOCK) : M;
    double tmp;
    for(int i = start; i < end; i++)
    {
        tmp = 0.0;
        for(int j = 0; j < N; j ++)
            tmp += U_t[i*N+j] * U_t[i*N+j];
        local_Alphas[i-start] = tmp;
        
        for(int j = 0; j < M; j++)
        {
            tmp = 0.0;
            for(int k = 0; k < N; k ++)
                tmp += U_t[i*N+k] * U_t[j*N+k];
            local_Gammas[(i-start)*N + j] = tmp;
        }
    }

    //gather Alphas;
    MPI_Gather(local_Alphas, ROW_BLOCK, MPI_DOUBLE, Alphas, ROW_BLOCK, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //gather Gammas;
    MPI_Gather(local_Gammas, ROW_BLOCK*N, MPI_DOUBLE, Gammas, ROW_BLOCK*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	gettimeofday(&tt, NULL);

  if( rank == 0 && (T=="-t" || P =="-t")){
    double elapsedTime = (tt.tv_sec - ts.tv_sec) * 1000.0;
    elapsedTime += (tt.tv_usec - ts.tv_usec) / 1000.0;
    cout<<"Time: "<<elapsedTime<<" ms."<<endl<<endl;
	}


// fix final result


  //Output time and iterations

  // Output the matrixes for debug
  
    if(rank == 0)
    {

    if(T== "-p" || P == "-p"){
    cout<<"Alphas"<<endl<<endl;
    for(int i =0; i<M; i++){

      for(int j =0; j<N;j++){
    		    
      	cout<<Alphas[i]<<"  ";
      }
      cout<<endl;
    }

    cout<<endl<<"Betas"<<endl<<endl;
    for(int i =0; i<M; i++){

     for(int j=0; j<N;j++){	  
        cout<<Alphas[j]<<"  ";
     }
     cout<<endl;
    }

    cout<<endl<<"Gammas"<<endl<<endl;
    for(int i =0; i<M; i++){
      for(int j =0; j<N; j++){

         cout<<Gammas[i*N+j]<<"  ";
      
       }
      cout<<endl;
    }

    }

    //Generate files for debug purpouse
     if(Db == "-d" || T == "-d" || P == "-d"){


      ofstream Af;
      //file for Matrix A
      Af.open("AlphasMPI.mat"); 
/*      Af<<"# Created from debug\n# name: A\n# type: matrix\n# rows: "<<M<<"\n# columns: "<<N<<"\n";*/

      Af<<M<<"  "<<N;
      for(int i = 0; i<M;i++){
        for(int j =0; j<N;j++){
          Af<<" "<<Alphas[i];
        }
        Af<<"\n";
      }
      
      Af.close();

      ofstream Uf;

      //File for Matrix U
      Uf.open("BetasMPI.mat");
/*      Uf<<"# Created from debug\n# name: Ugpu\n# type: matrix\n# rows: "<<M<<"\n# columns: "<<N<<"\n";*/
      
      for(int i = 0; i<M;i++){
        for(int j =0; j<N;j++){
          Uf<<" "<< Alphas[j];
        }
        Uf<<"\n";
      }
      Uf.close();

      ofstream Vf;
      //File for Matrix V
      Vf.open("GammasMPI.mat");
/*      Vf<<"# Created from debug\n# name: Vgpu\n# type: matrix\n# rows: "<<M<<"\n# columns: "<<N<<"\n";*/

      for(int i = 0; i<M;i++){
        for(int j =0; j<N;j++){
          Vf<<" "<<Gammas[i*N+j];
        }
        Vf<<"\n";
      }
      

      Vf.close();

      ofstream Sf;


   }  
    }

    
   
   delete [] local_Alphas;
   delete [] local_Gammas;

   delete [] U_t;
   if(rank == 0)
   {
        delete [] Alphas;
        delete [] Gammas;
   }
   MPI_Finalize();

  return 0;
}
