#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define WIDTH 1200
#define HEIGHT 1200
#define MAX_ITER 1000

int mandelbrot(double x, double y, double cx, double cy) {
    int iter = 0;
    double zr = x, zi = y, zr_new, zi_new;

    while (zr * zr + zi * zi < 4.0 && iter < MAX_ITER) {
        zr_new = zr * zr - zi * zi + cx;
        zi_new = 2.0 * zr * zi + cy;
        zr = zr_new;
        zi = zi_new;
        iter++;
    }

    return iter;
}

int main(int argc, char **argv) {
    double commtime;
    double start_time = MPI_Wtime();
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double xmin = -2.0, xmax = 2.0, ymin = -2.0, ymax = 2.0;
    double dx = (xmax - xmin) / WIDTH;
    double dy = (ymax - ymin) / HEIGHT;

    double cx_min = -1.5, cx_max = 0.5, cy_min = -1.0, cy_max = 1.0;
    double dcx = (cx_max - cx_min) / WIDTH;
    double dcy = (cy_max - cy_min) / HEIGHT;

    int rows_per_task = HEIGHT / size;
    int start_row = rank * rows_per_task;
    int end_row = (rank == size - 1) ? HEIGHT : start_row + rows_per_task;

    int *buffer = (int *) malloc(WIDTH * rows_per_task * sizeof(int));

    for (int row = start_row; row < end_row; row++) {
        for (int col = 0; col < WIDTH; col++) {
            double cx = cx_min + col * dcx;
            double cy = cy_min + row * dcy;
            buffer[(row - start_row) * WIDTH + col] = mandelbrot(0.0, 0.0, cx, cy);
        }
    }

    MPI_Request request;
    if (rank == 0) {
        int *output = (int *) malloc(WIDTH * HEIGHT * sizeof(int));
        for (int i = 0; i < rows_per_task * WIDTH; i++) {
            output[i] = buffer[i];
        }
        double starttime = MPI_Wtime();
        for (int i = 1; i < size; i++) {
  
            MPI_Irecv(output + i * rows_per_task * WIDTH, rows_per_task * WIDTH, MPI_INT, i, 0, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, MPI_STATUS_IGNORE);
        }
        double endtime = MPI_Wtime();
        commtime = endtime-starttime;

     
        FILE *fp = fopen("mandelbrot3.pgm", "wb");
        fprintf(fp, "P2\n%d %d\n%d\n", WIDTH, HEIGHT, MAX_ITER);
        for (int row = 0; row < HEIGHT; row++) {
            for (int col = 0; col < WIDTH; col++) {
                fprintf(fp, "%d ", output[row * WIDTH + col]);
            }
            fprintf(fp, "\n"); }
        fclose(fp);

        free(output);
    } else {
    	double starttime2 = MPI_Wtime();
   
        MPI_Isend(buffer, rows_per_task * WIDTH, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
	    double endtime2 = MPI_Wtime();
	    commtime += endtime2-starttime2;
    }
    double end_time = MPI_Wtime();
    double elapsed_time = end_time-start_time;
    printf("Elapsed time: %f seconds\n", elapsed_time);
    free(buffer);
    
    printf("commtime: %f seconds\n", commtime);
    MPI_Finalize();
    return 0;
}

 