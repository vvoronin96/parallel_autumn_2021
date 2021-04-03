/**********************************************************************
 * MPI-based matrix multiplication AxB=C
 *********************************************************************/


#include<iostream>
#include<stdlib.h>
#include<time.h>
#include <chrono>
#include <stdio.h>
#include "mpi.h"
using namespace std;
const int N = 900;
const int root = 0;
MPI_Status status;

double firstMatrix[N][N], secondMatrix[N][N], resultMatrix[N][N];

int main(int argc, char** argv)
{
    double totalTime = 0.0;

    int countThreads, rank, countWorker, source, targetWorker, rows, offset, i, j, k;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &countThreads);

    countWorker = countThreads - 1;

    /*---------------------------- master ----------------------------*/
    if (rank == root) {
        srand(time(0));
        
        //root thread generate two matrixes
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                firstMatrix[i][j] = rand() % 10;
                secondMatrix[i][j] = rand() % 10;
            }
        }

        auto t1 = std::chrono::high_resolution_clock::now();

        /* send matrix data to the worker tasks */
        rows = N / countWorker; // how many rows will be processed by the one worker
        offset = 0;

        for (targetWorker = 1; targetWorker <= countWorker; targetWorker++)
        {
            MPI_Send(&offset, 1, MPI_INT, targetWorker, 1, MPI_COMM_WORLD);
            MPI_Send(&rows, 1, MPI_INT, targetWorker, 1, MPI_COMM_WORLD);
            MPI_Send(&firstMatrix[offset][0], rows * N, MPI_DOUBLE, targetWorker, 1, MPI_COMM_WORLD);
            MPI_Send(&secondMatrix, N * N, MPI_DOUBLE, targetWorker, 1, MPI_COMM_WORLD);
            offset = offset + rows;
        }

        /* wait for results from all worker tasks */
        for (source = 1; source <= countWorker; source++)
        {
            MPI_Recv(&offset, 1, MPI_INT, source, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(&rows, 1, MPI_INT, source, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(&resultMatrix[offset][0], rows * N, MPI_DOUBLE, source, 2, MPI_COMM_WORLD, &status);
        }

        auto t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::micro> fp_ms = t2 - t1;
        totalTime = fp_ms.count();

        /*printf("Here is the result matrix:\n");
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++)
                printf("%6.2f   ", resultMatrix[i][j]);
            printf("\n");
        }*/
        cout <<"Time spent for order "<<N<<"|\t"<< totalTime <<" ms";
    }

    /*---------------------------- worker----------------------------*/
    if (rank > 0) {
        MPI_Recv(&offset, 1, MPI_INT, root, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&rows, 1, MPI_INT, root, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&firstMatrix, rows * N, MPI_DOUBLE, root, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&secondMatrix, N * N, MPI_DOUBLE, root, 1, MPI_COMM_WORLD, &status);

        /* Matrix multiplication */
        for (k = 0; k < N; k++)
            for (i = 0; i < rows; i++) {
                resultMatrix[i][k] = 0.0;
                for (j = 0; j < N; j++)
                    resultMatrix[i][k] = resultMatrix[i][k] + firstMatrix[i][j] * secondMatrix[j][k];
            }


        MPI_Send(&offset, 1, MPI_INT, root, 2, MPI_COMM_WORLD);
        MPI_Send(&rows, 1, MPI_INT, root, 2, MPI_COMM_WORLD);
        MPI_Send(&resultMatrix, rows * N, MPI_DOUBLE, root, 2, MPI_COMM_WORLD);
    }

    MPI_Finalize();
}