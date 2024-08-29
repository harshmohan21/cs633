#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

int main(int argc, char *argv[])
{
    int myrank, size;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int Px = atoi(argv[1]);
    int N = atoi(argv[2]);
    int num_time_steps = atoi(argv[3]);
    int x = num_time_steps;
    int seed = atoi(argv[4]);
    int stencil = atoi(argv[5]);
    int n = sqrt(N);
    int Py = size / Px;

    double data_set[n][n] ;

    srand(seed * (myrank + 10));
    // printf("ORIGINIAL DATA of rank %d \n",myrank);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            data_set[i][j] = abs(rand() + (i * rand() + j * myrank)) / 100;
            // data_set[i][j] = myrank;
            // printf("%f ",data_set[i][j]);
        }
        // printf("\n");
    }
    double sTime, eTime, time, maxTime;
    sTime = MPI_Wtime();
    while(num_time_steps--)
    {
        // 
        int pos1 = 0;
        int pos2 = 0;
        int pos3 = 0;
        int pos4 = 0;
        double temp[n][n];
        if (stencil == 5) {

            double send_buff_top[n], send_buff_bot[n];
            double send_buff_left[n], send_buff_right[n];
            double recv_buff_top[n], recv_buff_bot[n];
            double recv_buff_left[n], recv_buff_right[n];

            // Pack and send top and bottom borders
            for (int i = 0; i < n; i++) {
                MPI_Pack(&data_set[0][i], 1, MPI_DOUBLE, send_buff_top, n * sizeof(double), &pos1, MPI_COMM_WORLD);
                MPI_Pack(&data_set[n - 1][i], 1, MPI_DOUBLE, send_buff_bot, n * sizeof(double), &pos2, MPI_COMM_WORLD);
            }

            if ((myrank/Py) % 2 == 0 && (myrank/Py) < Px-1) 
            {
                // Send/recv bottom neighbor from even ranks
                MPI_Send(send_buff_bot, n * sizeof(double), MPI_PACKED, myrank + Py, myrank+Py, MPI_COMM_WORLD);
                MPI_Recv(recv_buff_bot, n * sizeof(double), MPI_DOUBLE, myrank + Py, myrank, MPI_COMM_WORLD, &status);
            }

            else if ((myrank/Py) % 2 != 0 && (myrank/Py) > 0) // 0 is anyway excluded here
            {
                // Send/recv left neighbor
                MPI_Recv(recv_buff_top, n * sizeof(double), MPI_DOUBLE, myrank - Py, myrank, MPI_COMM_WORLD, &status);
                MPI_Send(send_buff_top, n * sizeof(double), MPI_PACKED, myrank - Py, myrank-Py, MPI_COMM_WORLD);
            }

            if ((myrank/Py) % 2 != 0 && (myrank/Py) < Px-1) 
            {
                MPI_Send(send_buff_bot, n * sizeof(double), MPI_PACKED, myrank +Py, myrank+Py, MPI_COMM_WORLD);
                MPI_Recv(recv_buff_bot, n * sizeof(double), MPI_DOUBLE, myrank + Py, myrank, MPI_COMM_WORLD, &status);
            }

            else if ((myrank/Py) % 2 == 0 && (myrank/Py) > 0) // 0 is anyway excluded here
            {
                MPI_Recv(recv_buff_top, n * sizeof(double), MPI_DOUBLE, myrank - Py, myrank, MPI_COMM_WORLD, &status);
                MPI_Send(send_buff_top, n * sizeof(double), MPI_PACKED, myrank - Py, myrank-Py, MPI_COMM_WORLD);
            }
            // Pack and send left and right borders
            for (int i = 0; i < n; i++) {
                MPI_Pack(&data_set[i][0], 1, MPI_DOUBLE, send_buff_left, n * sizeof(double), &pos3, MPI_COMM_WORLD);
                MPI_Pack(&data_set[i][n - 1], 1, MPI_DOUBLE, send_buff_right, n * sizeof(double), &pos4, MPI_COMM_WORLD);
            }

            if ((myrank%Py) % 2 == 0 && (myrank%Py) < Py-1) // Py-1 is anyway odd (assume)
            {
                // Send/recv right neighbor from even ranks
                MPI_Send(send_buff_right, n * sizeof(double), MPI_PACKED, myrank + 1, myrank+1, MPI_COMM_WORLD);
                MPI_Recv(recv_buff_right, n * sizeof(double), MPI_DOUBLE, myrank + 1, myrank, MPI_COMM_WORLD, &status);
            }

            else if ((myrank%Py) % 2 != 0 && (myrank%Py) > 0) // 0 is anyway excluded here
            {
                // Send/recv left neighbor
                MPI_Recv(recv_buff_left, n * sizeof(double), MPI_DOUBLE, myrank - 1, myrank, MPI_COMM_WORLD, &status);
                MPI_Send(send_buff_left, n * sizeof(double), MPI_PACKED, myrank - 1, myrank-1, MPI_COMM_WORLD);
            }

            if ((myrank%Py) % 2 != 0 && (myrank%Py) < Py-1) // Py-1 is anyway odd (assume)
            {
                MPI_Send(send_buff_right, n * sizeof(double), MPI_PACKED, myrank +1, myrank+1, MPI_COMM_WORLD);
                MPI_Recv(recv_buff_right, n * sizeof(double), MPI_DOUBLE, myrank + 1, myrank, MPI_COMM_WORLD, &status);
            }

            else if ((myrank%Py) % 2 == 0 && (myrank%Py) > 0) // 0 is anyway excluded here
            {
                MPI_Recv(recv_buff_left, n * sizeof(double), MPI_DOUBLE, myrank - 1, myrank, MPI_COMM_WORLD, &status);
                MPI_Send(send_buff_left, n * sizeof(double), MPI_PACKED, myrank - 1, myrank-1, MPI_COMM_WORLD);
            }
        
            // Compute average and update temp
            // printf("Result data_set of rank %d\n",myrank);
            for (int row = 0; row < n; row++) {
                for (int col = 0; col < n; col++) {
                    int sum_sum_count = 1;
                    temp[row][col] = data_set[row][col];

                    if (col > 0)
                    {
                        temp[row][col] += data_set[row][col - 1] ;
                        sum_sum_count++;
                    }
                    else if (col == 0 && myrank % Py != 0){
                        temp[row][col] += recv_buff_left[row];
                        sum_sum_count++;
                    }

                    if (col < n - 1){
                        temp[row][col] += data_set[row][col + 1];
                        sum_sum_count++;
                    }
                    else if (col == n - 1 && myrank % Py != Py - 1){
                        temp[row][col] += recv_buff_right[row];
                        sum_sum_count++;
                    }
                    if (row > 0){
                        temp[row][col] += data_set[row - 1][col];
                        sum_sum_count++;
                    }
                    else if (row == 0 && myrank >= Py){
                        temp[row][col] += recv_buff_top[col];
                        sum_sum_count++;
                    }
                    if (row < n - 1){
                        temp[row][col] += data_set[row + 1][col];
                        sum_sum_count++;
                    }
                    else if (row == n - 1 && myrank < size - Py){
                        temp[row][col] += recv_buff_bot[col];
                        sum_sum_count++;
                    }

                    temp[row][col] /= sum_sum_count;
                    // printf("%f ",temp[row][col]);
                }
    // printf("\n");
            }

        }

        else if (stencil == 9) {
            double send_buff_top[2*n], send_buff_bot[2*n];
            double send_buff_left[2*n], send_buff_right[2*n];
            double recv_buff_top[2*n], recv_buff_bot[2*n];
            double recv_buff_left[2*n], recv_buff_right[2*n];

            // Pack and send top and bottom borders
            for (int j = 0; j < 2; j++) {
                for (int i = 0; i < n; i++) {
                    MPI_Pack(&data_set[j][i], 1, MPI_DOUBLE, send_buff_top, n * sizeof(double), &pos1, MPI_COMM_WORLD);
                    MPI_Pack(&data_set[n - j - 1][i], 1, MPI_DOUBLE, send_buff_bot, n * sizeof(double), &pos2, MPI_COMM_WORLD);
                }
            }

            if ((myrank/Py) % 2 == 0 && (myrank/Py) < Px-1) 
            {
                // Send/recv bottom neighbor from even ranks
                MPI_Send(send_buff_bot, n * sizeof(double), MPI_PACKED, myrank + Py, myrank+Py, MPI_COMM_WORLD);
                MPI_Recv(recv_buff_bot, n * sizeof(double), MPI_DOUBLE, myrank + Py, myrank, MPI_COMM_WORLD, &status);
            }

            else if ((myrank/Py) % 2 != 0 && (myrank/Py) > 0) // 0 is anyway excluded here
            {
                // Send/recv top neighbor
                MPI_Recv(recv_buff_top, n * sizeof(double), MPI_DOUBLE, myrank - Py, myrank, MPI_COMM_WORLD, &status);
                MPI_Send(send_buff_top, n * sizeof(double), MPI_PACKED, myrank - Py, myrank-Py, MPI_COMM_WORLD);
            }

            if ((myrank/Py) % 2 != 0 && (myrank/Py) < Px-1) 
            {
                MPI_Send(send_buff_bot, n * sizeof(double), MPI_PACKED, myrank +Py, myrank+Py, MPI_COMM_WORLD);
                MPI_Recv(recv_buff_bot, n * sizeof(double), MPI_DOUBLE, myrank + Py, myrank, MPI_COMM_WORLD, &status);
            }

            else if ((myrank/Py) % 2 == 0 && (myrank/Py) > 0) // 0 is anyway excluded here
            {
                MPI_Recv(recv_buff_top, n * sizeof(double), MPI_DOUBLE, myrank - Py, myrank, MPI_COMM_WORLD, &status);
                MPI_Send(send_buff_top, n * sizeof(double), MPI_PACKED, myrank - Py, myrank-Py, MPI_COMM_WORLD);
            }
            // Pack and send left and right borders
            for (int j = 0; j < 2; j++) {
                for (int i = 0; i < n; i++) {
                    MPI_Pack(&data_set[i][j], 1, MPI_DOUBLE, send_buff_left, n * sizeof(double), &pos3, MPI_COMM_WORLD);
                    MPI_Pack(&data_set[i][n - j - 1], 1, MPI_DOUBLE, send_buff_right, n * sizeof(double), &pos4, MPI_COMM_WORLD);
                }
            }

            if ((myrank%Py) % 2 == 0 && (myrank%Py) < Py-1) // Py-1 is anyway odd (assume)
            {
                // Send/recv right neighbor from even ranks
                MPI_Send(send_buff_right, n * sizeof(double), MPI_PACKED, myrank + 1, myrank+1, MPI_COMM_WORLD);
                MPI_Recv(recv_buff_right, n * sizeof(double), MPI_DOUBLE, myrank + 1, myrank, MPI_COMM_WORLD, &status);
            }

            else if ((myrank%Py) % 2 != 0 && (myrank%Py) > 0) // 0 is anyway excluded here
            {
                // Send/recv left neighbor
                MPI_Recv(recv_buff_left, n * sizeof(double), MPI_DOUBLE, myrank - 1, myrank, MPI_COMM_WORLD, &status);
                MPI_Send(send_buff_left, n * sizeof(double), MPI_PACKED, myrank - 1, myrank-1, MPI_COMM_WORLD);
            }

            if ((myrank%Py) % 2 != 0 && (myrank%Py) < Py-1) // Py-1 is anyway odd (assume)
            {
                MPI_Send(send_buff_right, n * sizeof(double), MPI_PACKED, myrank +1, myrank+1, MPI_COMM_WORLD);
                MPI_Recv(recv_buff_right, n * sizeof(double), MPI_DOUBLE, myrank + 1, myrank, MPI_COMM_WORLD, &status);
            }

            else if ((myrank%Py) % 2 == 0 && (myrank%Py) > 0) // 0 is anyway excluded here
            {
                MPI_Recv(recv_buff_left, n * sizeof(double), MPI_DOUBLE, myrank - 1, myrank, MPI_COMM_WORLD, &status);
                MPI_Send(send_buff_left, n * sizeof(double), MPI_PACKED, myrank - 1, myrank-1, MPI_COMM_WORLD);
            }

            // Compute average and update temp
            for (int row = 0; row < n; row++) {
                for (int col = 0; col < n; col++) {
                    int sum_count = 1;
                    temp[row][col] = data_set[row][col];

                    if (col > 0){
                        temp[row][col] += data_set[row][col - 1];
                        sum_count++;
                    }
                    else if (col == 0 && myrank % Py != 0){
                        temp[row][col] += recv_buff_left[row];
                        sum_count++;
                    }
                    if (col > 1)
                    {
                        temp[row][col] += data_set[row][col - 2];
                        sum_count++;
                    }
                    else if (col == 1 && myrank % Py != 0){
                        temp[row][col] += recv_buff_left[row];
                        sum_count++;
                    }
                    else if (col == 0 && myrank % Py != 0){
                        temp[row][col] += recv_buff_left[row + n];
                        sum_count++;
                    }

                    if (col < n - 2){
                        temp[row][col] += data_set[row][col + 2];
                        sum_count++;
                    }
                    else if (col == n - 2 && myrank % Py != Py - 1){
                        temp[row][col] += recv_buff_right[row];
                        sum_count++;
                    }
                    else if (col == n - 1 && myrank % Py != Py - 1){
                        temp[row][col] += recv_buff_right[row + n];
                        sum_count++;
                    }

                    if (col < n - 1){
                        temp[row][col] += data_set[row][col + 1];
                        sum_count++;
                    }
                    else if (col == n - 1 && myrank % Py != Py - 1){
                        temp[row][col] += recv_buff_right[row];
                        sum_count++;
                    }

                    if (row > 0){
                        temp[row][col] += data_set[row - 1][col];
                        sum_count++;
                    }
                    else if (row == 0 && myrank >= Py){
                        temp[row][col] += recv_buff_top[col];
                        sum_count++;
                    }
                    if (row > 1){
                        temp[row][col] += data_set[row - 2][col];
                        sum_count++;
                    }
                    else if (row == 1 && myrank >= Py){
                        temp[row][col] += recv_buff_top[col];
                        sum_count++;
                    }
                    else if (row == 0 && myrank >= Py){
                        temp[row][col] += recv_buff_top[col + n];
                        sum_count++;
                    }
 
                    if (row < n - 1){
                        temp[row][col] += data_set[row + 1][col];
                        sum_count++;
                    }
                    else if (row == n - 1 && myrank < size - Py){
                        temp[row][col] += recv_buff_bot[col];
                        sum_count++;
                    }
                   if (row < n - 2){
                        temp[row][col] += data_set[row + 2][col];
                        sum_count++;
                    }
                    else if (row == n - 2 && myrank < size - Py){
                        temp[row][col] += recv_buff_bot[col];
                        sum_count++;
                    }
                    else if (row == n - 1 && myrank < size - Py){
                        temp[row][col] += recv_buff_bot[col + n];
                        sum_count++;
                    }
                    temp[row][col] /= sum_count;
                }
            }
        }
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < size; j++)
            {
                data_set[i][j] = temp[i][j];
            }
        }
    }   
    eTime = MPI_Wtime();
    time = eTime - sTime; 
    // The maximum execution time across all processes is determined using MPI_Reduce
    MPI_Reduce(&time, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if(myrank == 0)
        printf("Time : %f  \n",time);
    MPI_Finalize();
}