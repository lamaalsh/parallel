


#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define WIDTH 800
#define HEIGHT 800

int main(int argc, char** argv) {
    double commtime;
    double start_time = MPI_Wtime();
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double real_min = -1.5;
    double real_max = 0.5;
    double imag_min = -1.0;
    double imag_max = 1.0;
    int max_iter = 500;

    int rows_per_process = HEIGHT / size;
    int start_row = rows_per_process * rank;
    int end_row = start_row + rows_per_process;

    if (rank == size - 1) {
        end_row = HEIGHT;
    }

    int* row = (int*)malloc(sizeof(int) * WIDTH);
    int* data = (int*)malloc(sizeof(int) * WIDTH * rows_per_process);

    for (int y = start_row; y < end_row; y++) {
        for (int x = 0; x < WIDTH; x++) {
            double c_real = real_min + (real_max - real_min) * x / WIDTH;
            double c_imag = imag_min + (imag_max - imag_min) * y / HEIGHT;

            double z_real = 0.0;
            double z_imag = 0.0;

            int iter = 0;
            while (z_real * z_real + z_imag * z_imag < 4.0 && iter < max_iter) {
                double next_z_real = z_real * z_real - z_imag * z_imag + c_real;
                double next_z_imag = 2.0 * z_real * z_imag + c_imag;
                z_real = next_z_real;
                z_imag = next_z_imag;
                iter++;
            }

            if (iter == max_iter) {
                row[x] = 0;
            } else {
                row[x] = iter % 256;
            }
        }
        int row_index = (y - start_row) * WIDTH;
        for (int x = 0; x < WIDTH; x++) {
            data[row_index + x] = row[x];
        }
    }

    free(row);

    int* final_data = NULL;
    if (rank == 0) {
        final_data = (int*)malloc(sizeof(int) * WIDTH * HEIGHT);
    }
    
    double starttime = MPI_Wtime();
    MPI_Gather(data, WIDTH * rows_per_process, MPI_INT, final_data, WIDTH * rows_per_process, MPI_INT, 0, MPI_COMM_WORLD);
    double endtime = MPI_Wtime();
    commtime = endtime-starttime;

    free(data);

    if (rank == 0) {
        FILE* fp = fopen("namnom.pgm", "wb");
        fprintf(fp, "P5\n%d %d\n255\n", WIDTH, HEIGHT);
        for (int i = 0; i < WIDTH * HEIGHT; i++) {
            fputc(final_data[i], fp);
        }
        fclose(fp);
        free(final_data);
    }
    double end_time = MPI_Wtime();
    double elapsed_time = end_time-start_time;
    printf("Elapsed time: %f seconds\n", elapsed_time);

    printf("commtime: %f seconds\n", commtime);
    MPI_Finalize();
    return 0;
}
