#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <assert.h>
#include <time.h>
#include "mpi.h"


#define GRID_WIDTH 1048576
#define GRID_SIZ 1024
#define SRAND_VALUE 1985

int main(int argc, char **argv)
{
    struct timeval inicio,inicio2,final,final2;
    gettimeofday(&inicio, NULL);
    int i = 0;
    int ResultGrid[GRID_WIDTH];
    srand(SRAND_VALUE);

    for(i = 0; i<GRID_SIZ; i++) {
        for(int j = 0; j<GRID_SIZ; j++) {
            ResultGrid[(i*GRID_SIZ)+j]= rand() % 2;
        }
    }
    int size;
    int Rank, j;
    int it = 0;
    int Gen = 200;
    int result = 0;

    gettimeofday(&inicio2, NULL);
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);

    assert(GRID_SIZ % size == 0);
    int *v = (int *)malloc(GRID_SIZ * ((GRID_SIZ / size) + 2) * sizeof(int));
    for (it = 0; it < Gen; it++)
    {
        j = GRID_SIZ;
        for (i = Rank * (GRID_WIDTH / size); i < (Rank + 1) * (GRID_WIDTH / size); i++)
        {
            v[j] = ResultGrid[i];
            j++;
        }

        if (size != 1)
        {
            int rcv1[GRID_SIZ],rcv2[GRID_SIZ],sd1[GRID_SIZ],sd2[GRID_SIZ];
            if (Rank % 2 == 0)
            {
                for (i = 0; i < GRID_SIZ; i++)
                {
                    sd1[i] = v[i + GRID_SIZ];
                }
                result = (Rank - 1)%size;
                if(result< 0 )
                    result = result + size;           
                MPI_Ssend(&sd1, GRID_SIZ, MPI_INTEGER,result , 1, MPI_COMM_WORLD);
                for (i = 0; i < GRID_SIZ; i++)
                {
                    sd2[i] = v[(GRID_SIZ * (GRID_SIZ / size)) + i];
                }
                result = (Rank + 1)%size;
                if(result< 0 )
                    result = result + size;                
                MPI_Ssend(&sd2, GRID_SIZ, MPI_INTEGER, result, 1, MPI_COMM_WORLD);
            }
            else
            {
                result = (Rank + 1)%size;
                if(result< 0 )
                    result = result + size; 
                MPI_Recv(&rcv2, GRID_SIZ, MPI_INTEGER, result, 1, MPI_COMM_WORLD, &status);
                result = (Rank - 1)%size;
                if(result< 0 )
                    result = result + size; 
                MPI_Recv(&rcv1, GRID_SIZ, MPI_INTEGER, result, 1, MPI_COMM_WORLD, &status);
            }
            if (Rank % 2 == 0)
            {
                result = (Rank + 1)%size;
                if(result< 0 )
                    result = result + size; 
                MPI_Recv(&rcv2, GRID_SIZ, MPI_INTEGER, result, 1, MPI_COMM_WORLD, &status);
                result = (Rank - 1)%size;
                if(result< 0 )
                    result = result + size; 
                MPI_Recv(&rcv1, GRID_SIZ, MPI_INTEGER, result, 1, MPI_COMM_WORLD, &status);
            }
            else
            {
                for (i = 0; i < GRID_SIZ; i++)
                {
                    sd1[i] = v[i + GRID_SIZ];
                }
                result = (Rank - 1)%size;
                if(result< 0 )
                    result = result + size; 
                MPI_Ssend(&sd1, GRID_SIZ, MPI_INTEGER, result, 1, MPI_COMM_WORLD);
                for (i = 0; i < GRID_SIZ; i++)
                {
                    sd2[i] = v[(GRID_SIZ * (GRID_SIZ / size)) + i];
                }
                result = (Rank + 1)%size;
                if(result< 0 )
                    result = result + size; 
                MPI_Ssend(&sd2, GRID_SIZ, MPI_INTEGER, result, 1, MPI_COMM_WORLD);
            }
            for (i = 0; i < GRID_SIZ; i++)
            {
                v[i] = rcv1[i];
                v[(GRID_SIZ * ((GRID_SIZ / size) + 1)) + i] = rcv2[i];
            }
        }
        else
        {
            for (i = 0; i < GRID_SIZ; i++)
            {
                v[i + GRID_WIDTH + GRID_SIZ] = ResultGrid[i];
            }
            for (i = GRID_WIDTH; i < GRID_WIDTH + GRID_SIZ; i++)
            {
                v[i - GRID_WIDTH] = ResultGrid[i - GRID_SIZ];
            }
        }
        int * new_g = (int *)malloc(GRID_SIZ * ((GRID_SIZ / size)) * sizeof(int));

        for (i = GRID_SIZ; i < GRID_SIZ * ((GRID_SIZ / size) + 1); i++)
        {
            int tot = GRID_SIZ * (GRID_SIZ / size) + 2;
            int rst = i / GRID_SIZ;
            int rest = i % GRID_SIZ;
                result = (rst - 1)%tot;
                if(result< 0 )
                    result = result + tot; 
            int result_ant = result;
                    result = (rst + 1)%tot;
                if(result< 0 )
                    result = result + tot;
            int prox_result = result;
                    result = (rest + 1)%GRID_SIZ;
                if(result< 0 )
                    result = result + GRID_SIZ;
            int prox_rest = result;
                    result = (rest - 1)%GRID_SIZ;
                if(result< 0 )
                    result = result + GRID_SIZ;            
            int rest_ant = result;
   

            int Cells_Alive= v[result_ant * GRID_SIZ + rest_ant] + v[result_ant * GRID_SIZ + rest] + v[result_ant * GRID_SIZ + prox_rest] + v[rst * GRID_SIZ + rest_ant] + v[rst * GRID_SIZ + prox_rest] + v[prox_result * GRID_SIZ + rest_ant] + v[prox_result * GRID_SIZ + rest] + v[prox_result * GRID_SIZ + prox_rest];
                if(v[i] ==1){
                if (Cells_Alive< 2 ||Cells_Alive>3 )
                    new_g[i - GRID_SIZ] = 0;
                else
                {
                   new_g[i - GRID_SIZ] = 1;
                }
                }
                else{    
                if (Cells_Alive== 3)
                    new_g[i - GRID_SIZ] = 1;
                else 
                    new_g[i - GRID_SIZ] = 0;
                }
        }

        j = 0;
        for (i = Rank * (GRID_WIDTH / size); i < (Rank + 1) * (GRID_WIDTH / size); i++)
        {
            ResultGrid[i] = new_g[j];
            j++;
        }
        MPI_Gather(new_g, GRID_SIZ * (GRID_SIZ / size), MPI_INTEGER, &ResultGrid, GRID_SIZ * (GRID_SIZ / size), MPI_INTEGER, 0, MPI_COMM_WORLD);
    }
    free(v);
    MPI_Finalize(); 
    if(Rank == 0){
        gettimeofday(&final2, NULL);
        int time = (int) (1000 * (final2.tv_sec - inicio2.tv_sec) + (final2.tv_usec - inicio2.tv_usec) / 1000);

        printf(" Tempo Decorrido seção paralela: %d ms\n",time );

        int Cells = 0;
        for(int a =0; a <GRID_SIZ/2; a++){
            for(int b = 0; b < GRID_SIZ; b++){
                if (ResultGrid[(a*GRID_SIZ)+b])
                    Cells = Cells + 1;
                if (ResultGrid[((GRID_SIZ-1-a)*GRID_SIZ)+b])
                    Cells = Cells + 1;
            }
        }

        printf(" Celulas vivas: %d \n",Cells);
        gettimeofday(&final, NULL);
        int time2 = (int) (1000 * (final.tv_sec - inicio.tv_sec) + (final.tv_usec - inicio.tv_usec) / 1000);
        printf(" Tempo total Decorrido: %d ms\n",time2 );
    }
}