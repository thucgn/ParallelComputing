/*************************************************************************
	> File Name: mpi_parallelio.cpp
	> Author: 
	> Mail: 
	> Created Time: 四  6/21 20:08:21 2018
 ************************************************************************/

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include "mpi.h"
#include "util.h"

using namespace std;

int read_file(std::string, int* buf)
{
    int tmp;
    FILE* f = fopen("pio", "rb");
    while(fread(&tmp, sizeof(int), 1, f) > 0)
        cout << tmp << " ";
    cout << endl;
    fclose(f);
    
}

int main(int argc, char** argv)
{

    int rank, size;

    //array
    const int N = 128;
    int array[N] = {0};
    int recv_array[N] = {0};
    for(int i = 0; i < N; i++)
        array[i] = i;

    MPI_File fh;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const int BLOCK_SIZE = N / size;

    MPI_Datatype etype, filetype;
    MPI_Offset disp;

    etype = MPI_INT;
    MPI_Type_contiguous(BLOCK_SIZE, MPI_INT, &filetype);
    MPI_Type_commit(&filetype);

    MPI_File_open(MPI_COMM_WORLD, "pio", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    disp = rank * BLOCK_SIZE * sizeof(int);

    MPI_File_set_view(fh, disp, etype, filetype, "native", MPI_INFO_NULL);
    
    MPI_Status status;
    MPI_File_write_all(fh, &array[rank * BLOCK_SIZE], BLOCK_SIZE, etype, &status);

    MPI_File_close(&fh);

    //MPI_Gather(&array[rank*2], 2, contig, &ret[rank], 2, contig, 0, MPI_COMM_WORLD); 
    MPI_Gather(&array[rank * BLOCK_SIZE], 1, filetype, recv_array, 1, filetype, 0, MPI_COMM_WORLD);

    if(rank == 0)
    {
        for(int i = 0;i < N; i ++)
            cout << recv_array[i] << " ";
        cout << endl;
    }

    int *data = new int[N+1];
    int cnt = 0;

    if(rank == 0)
    {
        int tmp;
        FILE* f = fopen("pio", "rb");

        while(fread(&tmp, sizeof(int), 1, f) > 0)
        {
            cout << tmp << " ";
            data[cnt++] = tmp;
        }
        cout << endl;
        fclose(f);

        bool isequal = true;
        if(cnt != N)
            isequal = false;
        for(int i = 0;i < N; i ++)
            if(recv_array[i] != data[i])
                isequal = false;
        if(isequal)
            cout << "Right! : The gathered data is equal to data writed in file" << endl;
        else 
            cout << "ERROR! : The gathered data is not equal to data writed in file" << endl;
            
    }

    MPI_Finalize();

    return 0;
}

