# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <complex.h>
# include <mpi.h>

# include "hhg_mpi.h"

int main(){
    int my_rank, comm_sz, i, n_t = (7337);

    double local_alpha;
    double _Complex * local_cur_den;
    double _Complex * total_cur_den;
    int source;

    total_cur_den = ( double _Complex * ) calloc ( (7337)  , sizeof ( double _Complex ));

    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    local_alpha = 1.0+(my_rank*0.0);

    local_cur_den = hhg(local_alpha);

    if (my_rank != 0) {
        MPI_Send(local_cur_den, n_t, MPI_C_DOUBLE_COMPLEX, 0, 0, MPI_COMM_WORLD);
    }

    else {

        for ( i=0 ; i<n_t ; i++ ) {
            total_cur_den[i]=total_cur_den[i]+local_cur_den[i];
        }

        for (source = 1; source < comm_sz; source++) {
            MPI_Recv(local_cur_den, n_t, MPI_C_DOUBLE_COMPLEX, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for (i=0;i<n_t;i++) {
                total_cur_den[i]=total_cur_den[i]+local_cur_den[i];
            }

        }
    }

    if (my_rank == 0){
        for (i=0;i<n_t;i++){
            printf("%G\n",creal(total_cur_den[i]));
        }
    }

    free(local_cur_den);
    free(total_cur_den);
    MPI_Finalize();
    return 0;

}
