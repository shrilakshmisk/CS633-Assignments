#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include <math.h>
#include "mpi.h"

// #pragma prutor-mpi-args: -np 12 -ppn 4
// #pragma prutor-mpi-sysargs: 4 16 1 42 9
double	**data = NULL, **data2 = NULL;

// function to swap matrices without copying data
void swap_matrices (double ***data, double ***data2){
	double **temp = *data;
	*data = *data2;
	*data2 = temp;
}

// returns data[i][j] for valid i, j. else returns 0
double mat (int i, int j, int side){
	if(i >= 0 && i < side && j >= 0 && j < side){
        return data[i][j];
    }
    else{
        return 0;
    }
}

// prints the matrix for debugging
void print(int side){
    for (int i = 0; i < side; i++){
        for (int j = 0; j < side; j++){
            printf("%lf ", data[i][j]);
        }
        printf("\n");
    }
}

int main(int argc, char *argv[]) 
{
    // default values
	int N = 512*512, P = 12, steps = 10, seed = 42, stencil = 5, P_x = 4, P_y, side, myrank, num;	
    bool left_n, right_n, up_n, down_n;
    double t_start, t_end, t_diff;
    double *left_recv = NULL, *right_recv = NULL, *up_recv = NULL, *down_recv = NULL, *left_send = NULL, *right_send = NULL, *up_send = NULL, *down_send = NULL, *int_left_recv = NULL, *int_right_recv = NULL, *int_up_recv = NULL, *int_down_recv = NULL, *down_recv_leader = NULL, *down_send_leader = NULL, *up_recv_leader = NULL, *up_send_leader = NULL;

    // taking arguments from command line
	if(argc == 5){
		P_x = atoi(argv[1]);
		N = atoi(argv[2]);
		steps = atoi(argv[3]);
		seed = atoi(argv[4]);
// 		stencil = atoi(argv[5]);
	}
    else{
        printf("Please provide 5 arguments\n");
        // printf("%d \n", argc);
        return 0;
    }
    
    stencil = 9;
    
    char hostname[MPI_MAX_PROCESSOR_NAME];
	MPI_Init (&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Status status;
    int len;
    MPI_Get_processor_name (hostname, &len);
    // int coreID = sched_getcpu();
    
    // if(myrank == 0 || myrank == 4 || myrank == 8)printf("rank %d running on %s\n", myrank, hostname);
    

    P_y = P/P_x;
	side = (int)sqrt(N);
	if(stencil == 5){
        num = 1; // num is number of neighbours to include in stencil computation in each direction
    }
    else if(stencil == 9){
        num = 2;
    }

    int leader = 0;
    
    // allocating memory to matrices
	data = (double **)malloc(side * sizeof(double*));
	data2 = (double **)malloc(side * sizeof(double*));
	for(int i = 0; i < side; i++){
		data[i] = (double *)malloc(side * sizeof(double));
		data2[i] = (double *)malloc(side * sizeof(double));
	}

    // initialising the matrix
	for(int i = 0; i < side; i++){
		for(int j = 0; j < side; j++){
			srand(seed*(myrank + 10));
			data[i][j] = abs(rand() + (i*rand() + j*myrank))/100;
            // data[i][j] = 10*i + j;
		}
	}
	
// 	MPI_Barrier(MPI_COMM_WORLD);
    // print(side);
	// preprocessing the existence of neighbours for a particular process
    if(myrank % P_x == 0){
        left_n = false;
    }
    else{
        left_n = true;
    }
    if(myrank % P_x == P_x - 1){
        right_n = false;
    }
    else{
        right_n = true;
    }
    if(myrank / P_x == 0){
        up_n = false;
    }
    else{
        up_n = true;
    }
    if(myrank / P_x == P_y - 1){
        down_n = false;
    } 
    else{
        down_n = true;
    }
    
    
    // allocating memory based on neighbours and stencil to send and receive buffers
	if(left_n){
		left_recv = (double *)malloc(side * num * sizeof(double));
// 		int_left_recv = (double *)malloc(side * num * sizeof(double));
		left_send = (double *)malloc(side * num * sizeof(double));
	}
	if(right_n){
		right_recv = (double *)malloc(side * num * sizeof(double));
// 		int_right_recv = (double *)malloc(side * num * sizeof(double));
		right_send = (double *)malloc(side * num * sizeof(double));
	}
	if(up_n){
		up_recv = (double *)malloc(side * num * sizeof(double));
// 		int_up_recv = (double *)malloc(side * num * sizeof(double));
		up_send	= (double *)malloc(side * num * sizeof(double));
	}
	if(down_n){
		down_recv = (double *)malloc(side * num * sizeof(double));
// 		int_down_recv = (double *)malloc(side * num * sizeof(double));
		down_send = (double *)malloc(side * num * sizeof(double));
	}
	
	// Calling MPI_Barrier so that previous initializations and preprocessing do not contribute in blocking of earlier processes in further Send and Recv calls
    MPI_Barrier(MPI_COMM_WORLD);
    t_start = MPI_Wtime(); 
    
    // the communication + computation loop begins
    for(int step = 0; step < steps; step++){
        
        if(left_n){
            int l = 0;
            for(int j = 0; j < num; j++){
                for(int i = 0; i < side; i++){
                    // we pack all the data required in a single array, even for stencil = 9
                    MPI_Pack(&data[i][j], 1, MPI_DOUBLE, left_send, side * num * sizeof(double), &l, MPI_COMM_WORLD);
                }
            }
        }  
        
        if(right_n){
            int r = 0;
            for(int j = side - num; j < side; j++){
                for(int i = 0; i < side; i++){
                    MPI_Pack(&data[i][j], 1, MPI_DOUBLE, right_send, side * num * sizeof(double), &r, MPI_COMM_WORLD);
                }
            }
        } 

        if(up_n){
            int u = 0;
            for(int i = 0; i < num; i++){
                for(int j = 0; j < side; j++){
                    MPI_Pack(&data[i][j], 1, MPI_DOUBLE, up_send, side * num * sizeof(double), &u, MPI_COMM_WORLD); 
                }
            }
        }

        if(down_n){
            int d = 0;
            for(int i = side - num; i < side; i++){
                for(int j = 0; j < side; j++){
                    MPI_Pack(&data[i][j], 1, MPI_DOUBLE, down_send, side * num * sizeof(double), &d, MPI_COMM_WORLD); 
                }
            }
        }
        
        // The following commented code is Odd - Even algorithm for Nearest Neighbour communication
        
        // if(((myrank % P_x) % 2) == 0 && right_n){ 
        //     MPI_Send(right_send, side * num * sizeof(double), MPI_PACKED, myrank + 1, myrank + 1, MPI_COMM_WORLD);
        //     MPI_Recv(right_recv, side * num * sizeof(double), MPI_PACKED, myrank + 1, myrank, MPI_COMM_WORLD, &status);
        // }
        // else if (((myrank % P_x) % 2) == 1 && left_n){
        //     MPI_Recv(left_recv, side * num * sizeof(double), MPI_PACKED, myrank - 1, myrank, MPI_COMM_WORLD, &status);
        //     MPI_Send(left_send, side * num * sizeof(double), MPI_PACKED, myrank - 1, myrank-1, MPI_COMM_WORLD);
        // }
        // if(((myrank % P_x) % 2) == 0 && left_n){
        //     MPI_Send(left_send, side * num * sizeof(double), MPI_PACKED, myrank - 1, myrank - 1, MPI_COMM_WORLD);
        //     MPI_Recv(left_recv, side * num * sizeof(double), MPI_PACKED, myrank - 1, myrank, MPI_COMM_WORLD, &status);
        // }
        // else if (((myrank % P_x) % 2) == 1 && right_n){
        //     MPI_Recv(right_recv, side * num * sizeof(double), MPI_PACKED, myrank + 1, myrank, MPI_COMM_WORLD, &status);
        //     MPI_Send(right_send, side * num * sizeof(double), MPI_PACKED, myrank + 1, myrank + 1, MPI_COMM_WORLD);
        // }
        // if(((myrank / P_x) % 2 )== 0 && down_n){
        //     MPI_Send(down_send, side * num * sizeof(double), MPI_PACKED, myrank + P_x, myrank + P_x, MPI_COMM_WORLD);
        //     MPI_Recv(down_recv, side * num * sizeof(double), MPI_PACKED, myrank + P_x, myrank, MPI_COMM_WORLD, &status);
        // }
        // else if (((myrank / P_x) % 2) == 1 && up_n){
        //     MPI_Recv(up_recv, side * num * sizeof(double), MPI_PACKED, myrank - P_x, myrank, MPI_COMM_WORLD, &status);
        //     MPI_Send(up_send, side * num * sizeof(double), MPI_PACKED, myrank - P_x, myrank - P_x, MPI_COMM_WORLD);
        // }
        // if(((myrank / P_x) % 2) == 0 && up_n){
        //     MPI_Send(up_send, side * num * sizeof(double), MPI_PACKED, myrank - P_x, myrank - P_x, MPI_COMM_WORLD);
        //     MPI_Recv(up_recv, side * num * sizeof(double), MPI_PACKED, myrank - P_x, myrank, MPI_COMM_WORLD, &status);
        // }
        // else if (((myrank / P_x) % 2) == 1 && down_n){
        //     MPI_Recv(down_recv, side * num * sizeof(double), MPI_PACKED, myrank + P_x, myrank, MPI_COMM_WORLD, &status);
        //     MPI_Send(down_send, side * num * sizeof(double), MPI_PACKED, myrank + P_x, myrank + P_x, MPI_COMM_WORLD);
        // }
        
        // Right Sends followed by Left Sends algorithm for Nearest Neighbours communication 
        
        // Left - Right communication with appropriate tags
        if(right_n){
            MPI_Send(right_send, side * num * sizeof(double), MPI_PACKED, myrank + 1, myrank + 1, MPI_COMM_WORLD);
            MPI_Recv(right_recv, side * num * sizeof(double), MPI_PACKED, myrank + 1, myrank, MPI_COMM_WORLD, &status);
        }

        if(left_n){
            MPI_Recv(left_recv, side * num * sizeof(double), MPI_PACKED, myrank - 1, myrank, MPI_COMM_WORLD, &status);
            MPI_Send(left_send, side * num * sizeof(double), MPI_PACKED, myrank - 1, myrank - 1, MPI_COMM_WORLD);
        }
        // Up - Down communication with appropriate tags
        if(down_n){
            MPI_Send(down_send, side * num * sizeof(double), MPI_PACKED, myrank + P_x, myrank + P_x, MPI_COMM_WORLD);
            MPI_Recv(down_recv, side * num * sizeof(double), MPI_PACKED, myrank + P_x, myrank, MPI_COMM_WORLD, &status);
        }

        if(up_n){
            MPI_Recv(up_recv, side * num * sizeof(double), MPI_PACKED, myrank - P_x, myrank, MPI_COMM_WORLD, &status);
            MPI_Send(up_send, side * num * sizeof(double), MPI_PACKED, myrank - P_x, myrank - P_x, MPI_COMM_WORLD);
        }
        
        // stencil computation
		for(int i = 0; i < side; i++){
			for(int j = 0; j < side; j++){
				double sum = data[i][j];
                double num_n = (double)stencil; 
                if(stencil == 5){
                    
                    // Halo cells
                    if(j == 0){
                        if(left_n){
                            sum += left_recv[i];
                        }
                        else{
                            num_n-=1;
                        }
                    }
                    
                    if(j == side - 1){
                        if(right_n){
                            sum += right_recv[i];
                        }                         
                        else{
                            num_n-=1;                            
                        }
                    }
                        
                    if(i == 0){
                        if(up_n){
                            sum += up_recv[j];
                        }
                        else{
                            num_n-=1;
                        }
                    }

                    if(i == side - 1){
                        if(down_n){
                            sum += down_recv[j];
                        }
                        else{
                            num_n-=1;
                        }
                    }

                    sum = sum + mat(i, j - 1, side) + mat(i, j + 1, side) + mat(i - 1, j, side) + mat(i + 1, j, side); // Internal cells
                    
                    data2[i][j] = sum/num_n;
                }    
                if(stencil == 9){
                    // Halo cells 
                    if(j == 0){
                        if(left_n){
                            sum += left_recv[i];
                            sum += left_recv[i + side];
                        }
                        else{
                            num_n-=2;
                        }
                    }
                    
                    if(j == 1){
                        if(left_n){
                            sum += left_recv[i + side];
                        }
                        else{
                            num_n-=1;
                        }
                    }

                    if(j == side - 1){
                        if(right_n){
                            sum += right_recv[i];
                            sum += right_recv[i + side];
                        }
                        else{
                            num_n-=2;
                        }
                    }

                    if(j == side - 2){
                        if(right_n){
                            sum += right_recv[i];
                        }
                        else{
                            num_n-=1;
                        }
                    }

                    if(i == 0){
                        if(up_n){
                            sum += up_recv[j];
                            sum += up_recv[j + side];
                        }
                        else{
                            num_n-=2;
                        }
                    }
                    
                    if(i == 1){
                        if(up_n){
                            sum += up_recv[j + side];
                        }
                        else{
                            num_n-=1;
                        }
                    }

                    if(i == side - 1){
                        if(down_n){
                            sum += down_recv[j];
                            sum += down_recv[j + side];
                        }
                        else{
                            num_n-=2;
                        }
                    }

                    if(i == side - 2){
                        if(down_n){
                            sum += down_recv[j];
                        }
                        else{
                            num_n-=1;
                        }
                    }
                    sum = sum + mat(i, j - 1, side) + mat(i, j - 2, side) + mat(i, j + 1, side) + mat(i, j + 2, side) + mat(i - 1, j, side) + mat(i - 2, j, side) + mat(i + 1, j, side) + mat(i + 2, j, side); // Internal cells 

                    data2[i][j] = sum/num_n;
                }
			}
		}
		swap_matrices(&data, &data2); // getting new values to data matrix
	}
    t_end = MPI_Wtime() - t_start; // ending time 
    
    MPI_Reduce(&t_end, &t_diff, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); // reporting maximum time taken among all processes
    if(myrank == 0) printf("Time without leader = %lf\n", t_diff);
    // if(myrank == 0) printf("Data without leader = %lf\n", data[0][0]);
    // MPI_Barrier(MPI_COMM_WORLD);
    // print(side);

    
    // Freeing all memory allocations
    
    if (left_n) {
        free(left_recv);
        left_recv = NULL;
        free(left_send);
        left_send = NULL;
    }
    if (right_n) {
        free(right_recv);
        right_recv = NULL;
        free(right_send);
        right_send = NULL;
    }
    if (up_n) {
        free(up_recv);
        up_recv = NULL;
        free(up_send);
        up_send = NULL;
    }
    if (down_n) {
        free(down_recv);
        down_recv = NULL;
        free(down_send);
        down_send = NULL;
    }

    for (int i = 0; i < side; i++) {
        free(data[i]);
        free(data2[i]);
    }
    free(data);
    data = NULL;
    free(data2);
    data2 = NULL;





    ///////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////



    leader = 1;
    MPI_Barrier(MPI_COMM_WORLD);

    data = (double **)malloc(side * sizeof(double*));
	data2 = (double **)malloc(side * sizeof(double*));
	for(int i = 0; i < side; i++){
		data[i] = (double *)malloc(side * sizeof(double));
		data2[i] = (double *)malloc(side * sizeof(double));
	}

    // initialising the matrix
	for(int i = 0; i < side; i++){
		for(int j = 0; j < side; j++){
			srand(seed*(myrank + 10));
			data[i][j] = abs(rand() + (i*rand() + j*myrank))/100;
            // data[i][j] = 10*i + j;
		}
	}
// 	print(side);

	
	// preprocessing the existence of neighbours for a particular process
    if(myrank % P_x == 0){
        left_n = false;
    }
    else{
        left_n = true;
    }
    if(myrank % P_x == P_x - 1){
        right_n = false;
    }
    else{
        right_n = true;
    }
    if(myrank / P_x == 0){
        up_n = false;
    }
    else{
        up_n = true;
    }
    if(myrank / P_x == P_y - 1){
        down_n = false;
    } 
    else{
        down_n = true;
    }
    
    
    // allocating memory based on neighbours and stencil to send and receive buffers
	if(left_n){
		left_recv = (double *)malloc(side * num * sizeof(double));
// 		int_left_recv = (double *)malloc(side * num * sizeof(double));
		left_send = (double *)malloc(side * num * sizeof(double));
	}
	if(right_n){
		right_recv = (double *)malloc(side * num * sizeof(double));
// 		int_right_recv = (double *)malloc(side * num * sizeof(double));
		right_send = (double *)malloc(side * num * sizeof(double));
	}
	if(up_n){  
		up_recv = (double *)malloc(side * num * sizeof(double));
// 		int_up_recv = (double *)malloc(side * num * sizeof(double));
		up_send	= (double *)malloc(side * num * sizeof(double));
	}
    if(up_n && myrank % P_x == 0){ // leader rank
        up_recv_leader = (double *)malloc(P_x * side * num * sizeof(double));
// 		int_up_recv = (double *)malloc(side * num * sizeof(double));
		up_send_leader	= (double *)malloc(P_x * side * num * sizeof(double));
    }
	if(down_n){ 
		down_recv = (double *)malloc(side * num * sizeof(double));
// 		int_down_recv = (double *)malloc(side * num * sizeof(double));
		down_send = (double *)malloc(side * num * sizeof(double));
	}
    if(down_n && myrank % P_x == 0){ // leader rank
        down_recv_leader = (double *)malloc(P_x * side * num * sizeof(double));
// 		int_down_recv = (double *)malloc(side * num * sizeof(double));
		down_send_leader = (double *)malloc(P_x * side * num * sizeof(double));
    }

    // have to make subcommunicators

    int ranks[P_x];
    int start_rank = myrank - myrank % P_x;
    for(int i = 0; i < P_x; i++){
        ranks[i] = start_rank + i;
    }
    MPI_Group g_group;
    MPI_Comm_group (MPI_COMM_WORLD, &g_group);

    MPI_Group row_group;
    MPI_Group_incl (g_group, P_x, ranks, &row_group);

    MPI_Comm row_comm;
    MPI_Comm_create_group (MPI_COMM_WORLD, row_group, myrank / P_x, &row_comm);

    // int leader_ranks[P_y];
    // for(int i = 0; i < P_y; i ++){
    //     leader_ranks[i] = i * P_x;
    // }
    // MPI_Group leader_group;
    // if(myrank % P_x == 0){
    //     MPI_Group_incl (g_group, P_y, leader_ranks, &leader_group)
    // }
    // MPI_Comm leader_comm;
    // MPI_Comm_create_group (MPI_COMM_WORLD, leader_group, &leader_comm);
	
	// Calling MPI_Barrier so that previous initializations and preprocessing do not contribute in blocking of earlier processes in further Send and Recv calls
    MPI_Barrier(MPI_COMM_WORLD);
    t_start = MPI_Wtime(); 
    
    // the communication + computation loop begins
    for(int step = 0; step < steps; step++){
        
        if(left_n){
            int l = 0;
            for(int j = 0; j < num; j++){
                for(int i = 0; i < side; i++){
                    // we pack all the data required in a single array, even for stencil = 9
                    MPI_Pack(&data[i][j], 1, MPI_DOUBLE, left_send, side * num * sizeof(double), &l, MPI_COMM_WORLD);
                }
            }
        }  
        
        if(right_n){
            int r = 0;
            for(int j = side - num; j < side; j++){
                for(int i = 0; i < side; i++){
                    MPI_Pack(&data[i][j], 1, MPI_DOUBLE, right_send, side * num * sizeof(double), &r, MPI_COMM_WORLD);
                }
            }
        } 

        if(up_n){ 
            int u = 0;
            for(int i = 0; i < num; i++){
                for(int j = 0; j < side; j++){
                    MPI_Pack(&data[i][j], 1, MPI_DOUBLE, up_send, side * num * sizeof(double), &u, MPI_COMM_WORLD); 
                }
            }
        }

        if(down_n){ 
            int d = 0;
            for(int i = side - num; i < side; i++){
                for(int j = 0; j < side; j++){
                    MPI_Pack(&data[i][j], 1, MPI_DOUBLE, down_send, side * num * sizeof(double), &d, MPI_COMM_WORLD); 
                }
            }
        }
        
        // The following commented code is Odd - Even algorithm for Nearest Neighbour communication
        
        // if(((myrank % P_x) % 2) == 0 && right_n){ 
        //     MPI_Send(right_send, side * num * sizeof(double), MPI_PACKED, myrank + 1, myrank + 1, MPI_COMM_WORLD);
        //     MPI_Recv(right_recv, side * num * sizeof(double), MPI_PACKED, myrank + 1, myrank, MPI_COMM_WORLD, &status);
        // }
        // else if (((myrank % P_x) % 2) == 1 && left_n){
        //     MPI_Recv(left_recv, side * num * sizeof(double), MPI_PACKED, myrank - 1, myrank, MPI_COMM_WORLD, &status);
        //     MPI_Send(left_send, side * num * sizeof(double), MPI_PACKED, myrank - 1, myrank-1, MPI_COMM_WORLD);
        // }
        // if(((myrank % P_x) % 2) == 0 && left_n){
        //     MPI_Send(left_send, side * num * sizeof(double), MPI_PACKED, myrank - 1, myrank - 1, MPI_COMM_WORLD);
        //     MPI_Recv(left_recv, side * num * sizeof(double), MPI_PACKED, myrank - 1, myrank, MPI_COMM_WORLD, &status);
        // }
        // else if (((myrank % P_x) % 2) == 1 && right_n){
        //     MPI_Recv(right_recv, side * num * sizeof(double), MPI_PACKED, myrank + 1, myrank, MPI_COMM_WORLD, &status);
        //     MPI_Send(right_send, side * num * sizeof(double), MPI_PACKED, myrank + 1, myrank + 1, MPI_COMM_WORLD);
        // }
        // if(((myrank / P_x) % 2 )== 0 && down_n){
        //     MPI_Send(down_send, side * num * sizeof(double), MPI_PACKED, myrank + P_x, myrank + P_x, MPI_COMM_WORLD);
        //     MPI_Recv(down_recv, side * num * sizeof(double), MPI_PACKED, myrank + P_x, myrank, MPI_COMM_WORLD, &status);
        // }
        // else if (((myrank / P_x) % 2) == 1 && up_n){
        //     MPI_Recv(up_recv, side * num * sizeof(double), MPI_PACKED, myrank - P_x, myrank, MPI_COMM_WORLD, &status);
        //     MPI_Send(up_send, side * num * sizeof(double), MPI_PACKED, myrank - P_x, myrank - P_x, MPI_COMM_WORLD);
        // }
        // if(((myrank / P_x) % 2) == 0 && up_n){
        //     MPI_Send(up_send, side * num * sizeof(double), MPI_PACKED, myrank - P_x, myrank - P_x, MPI_COMM_WORLD);
        //     MPI_Recv(up_recv, side * num * sizeof(double), MPI_PACKED, myrank - P_x, myrank, MPI_COMM_WORLD, &status);
        // }
        // else if (((myrank / P_x) % 2) == 1 && down_n){
        //     MPI_Recv(down_recv, side * num * sizeof(double), MPI_PACKED, myrank + P_x, myrank, MPI_COMM_WORLD, &status);
        //     MPI_Send(down_send, side * num * sizeof(double), MPI_PACKED, myrank + P_x, myrank + P_x, MPI_COMM_WORLD);
        // }
        
        // Right Sends followed by Left Sends algorithm for Nearest Neighbours communication 
        
        // Left - Right communication with appropriate tags
        if(right_n){
            MPI_Send(right_send, side * num * sizeof(double), MPI_PACKED, myrank + 1, myrank + 1, MPI_COMM_WORLD);
            MPI_Recv(right_recv, side * num * sizeof(double), MPI_PACKED, myrank + 1, myrank, MPI_COMM_WORLD, &status);
        }

        if(left_n){
            MPI_Recv(left_recv, side * num * sizeof(double), MPI_PACKED, myrank - 1, myrank, MPI_COMM_WORLD, &status);
            MPI_Send(left_send, side * num * sizeof(double), MPI_PACKED, myrank - 1, myrank - 1, MPI_COMM_WORLD);
        }
        // Up - Down communication with appropriate tags
        if(down_n){ // gather to leader, then send to leader + 1, then scatter from leader + 1 
            MPI_Gather(down_send, side * num * sizeof(double), MPI_PACKED, down_send_leader, side * num * sizeof(double), MPI_PACKED, 0, row_comm);
            if(myrank % P_x == 0){
                MPI_Send(down_send_leader, P_x * side * num * sizeof(double), MPI_PACKED, myrank + P_x, myrank + P_x, MPI_COMM_WORLD);
                MPI_Recv(down_recv_leader, P_x * side * num * sizeof(double), MPI_PACKED, myrank + P_x, myrank, MPI_COMM_WORLD, &status);
            }
            MPI_Scatter(down_recv_leader, side * num * sizeof(double), MPI_PACKED, down_recv, side * num * sizeof(double), MPI_PACKED, 0, row_comm);
        }

        if(up_n){ // gather to leader, then send to leader - 1, then scatter from leader - 1 
            MPI_Gather(up_send, side * num * sizeof(double), MPI_PACKED, up_send_leader, side * num * sizeof(double), MPI_PACKED, 0, row_comm);
            if(myrank % P_x == 0){
                MPI_Recv(up_recv_leader, P_x * side * num * sizeof(double), MPI_PACKED, myrank - P_x, myrank, MPI_COMM_WORLD, &status);
                MPI_Send(up_send_leader, P_x * side * num * sizeof(double), MPI_PACKED, myrank - P_x, myrank - P_x, MPI_COMM_WORLD);
            }
            MPI_Scatter(up_recv_leader, side * num * sizeof(double), MPI_PACKED, up_recv, side * num * sizeof(double), MPI_PACKED, 0, row_comm);      
        }
        
        // stencil computation
		for(int i = 0; i < side; i++){
			for(int j = 0; j < side; j++){
				double sum = data[i][j];
                double num_n = (double)stencil; 
                if(stencil == 5){
                    
                    // Halo cells
                    if(j == 0){
                        if(left_n){
                            sum += left_recv[i];
                        }
                        else{
                            num_n-=1;
                        }
                    }
                    
                    if(j == side - 1){
                        if(right_n){
                            sum += right_recv[i];
                        }                         
                        else{
                            num_n-=1;                            
                        }
                    }
                        
                    if(i == 0){
                        if(up_n){
                            sum += up_recv[j]; 
                        }
                        else{
                            num_n-=1;
                        }
                    }

                    if(i == side - 1){
                        if(down_n){
                            sum += down_recv[j]; 
                        }
                        else{
                            num_n-=1;
                        }
                    }

                    sum = sum + mat(i, j - 1, side) + mat(i, j + 1, side) + mat(i - 1, j, side) + mat(i + 1, j, side); // Internal cells
                    
                    data2[i][j] = sum/num_n;
                }    
                if(stencil == 9){
                    // Halo cells 
                    if(j == 0){
                        if(left_n){
                            sum += left_recv[i];
                            sum += left_recv[i + side];
                        }
                        else{
                            num_n-=2;
                        }
                    }
                    
                    if(j == 1){
                        if(left_n){
                            sum += left_recv[i + side];
                        }
                        else{
                            num_n-=1;
                        }
                    }

                    if(j == side - 1){
                        if(right_n){
                            sum += right_recv[i];
                            sum += right_recv[i + side];
                        }
                        else{
                            num_n-=2;
                        }
                    }

                    if(j == side - 2){
                        if(right_n){
                            sum += right_recv[i];
                        }
                        else{
                            num_n-=1;
                        }
                    }

                    if(i == 0){
                        if(up_n){
                            sum += up_recv[j]; 
                            sum += up_recv[j + side]; 
                        }
                        else{
                            num_n-=2;
                        }
                    }
                    
                    if(i == 1){
                        if(up_n){
                            sum += up_recv[j + side]; 
                        }
                        else{
                            num_n-=1;
                        }
                    }

                    if(i == side - 1){
                        if(down_n){
                            sum += down_recv[j]; 
                            sum += down_recv[j + side]; 
                        }
                        else{
                            num_n-=2;
                        }
                    }

                    if(i == side - 2){
                        if(down_n){
                            sum += down_recv[j]; 
                        }
                        else{
                            num_n-=1;
                        }
                    }
                    sum = sum + mat(i, j - 1, side) + mat(i, j - 2, side) + mat(i, j + 1, side) + mat(i, j + 2, side) + mat(i - 1, j, side) + mat(i - 2, j, side) + mat(i + 1, j, side) + mat(i + 2, j, side); // Internal cells 

                    data2[i][j] = sum/num_n;
                }
			}
		}
		swap_matrices(&data, &data2); // getting new values to data matrix
	}
    t_end = MPI_Wtime() - t_start; // ending time 
    
    MPI_Reduce(&t_end, &t_diff, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); // reporting maximum time taken among all processes
    if(myrank == 0) printf("Time with leader = %lf\n", t_diff);
    if(myrank == 0) printf("Data = %lf\n", data[0][0]);
    // MPI_Barrier(MPI_COMM_WORLD);
// 	print(side);

    
    // Freeing all memory allocations
    // have to change this based on what we allocate
    
    if (left_n) {
        free(left_recv);
        left_recv = NULL;
        free(left_send);
        left_send = NULL;
    }
    if (right_n) {
        free(right_recv);
        right_recv = NULL;
        free(right_send);
        right_send = NULL;
    }
    if (up_n) {
        free(up_recv);
        up_recv = NULL;
        free(up_send);
        up_send = NULL;
    }
    if (down_n) {
        free(down_recv);
        down_recv = NULL;
        free(down_send);
        down_send = NULL;
    }
    if(down_n && myrank % P_x == 0){ 
        free(down_recv_leader);
        down_recv_leader = NULL;
        free(down_send_leader);
        down_send_leader = NULL;
    }
    if(up_n && myrank % P_x == 0){ 
        free(up_recv_leader);
        up_recv_leader = NULL;
        free(up_send_leader);
        up_send_leader = NULL;
    }

    for (int i = 0; i < side; i++) {
        free(data[i]);
        free(data2[i]);
    }
    free(data);
    data = NULL;
    free(data2);
    data2 = NULL;






  	MPI_Finalize();
	return 0;
}