#include <stdlib.h> 
#include <mpi.h> 
#include <random>

using namespace std;

int main(int argc, char* argv[])
{
	int n = 144000;
	int my_rank, size;
	int snd_buf, rcv_buf;
	int right, left;
	int sum, i;
	int npes;
	int myrank;
	int nlocal;
	int* elmnts;
	int* relmnts;
	int oddrank;
	int evenrank;
	int* wspace;
	MPI_Comm new_comm;
	int dims[max_dims],
		periods[max_dims],
		reorder;

	MPI_Aint    snd_displs[2], rcv_displs[2];
	int         snd_counts[2], rcv_counts[2];
	MPI_Datatype snd_types[2], rcv_types[2];

	MPI_Status  status;
	MPI_Request request;


	MPI_Init(&argc, &argv);

	/* Get process info. */
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	/* Set cartesian topology. */
	dims[0] = size;
	periods[0] = 1;
	reorder = 1;

	MPI_Cart_create(MPI_COMM_WORLD, max_dims, dims, periods,
		reorder, &new_comm);

	/* Get coords */
	MPI_Comm_rank(new_comm, &my_rank);
	/* MPI_Cart_coords(new_comm, my_rank, max_dims, my_coords); */

	/* Get nearest neighbour rank. */
	MPI_Cart_shift(new_comm, 0, 1, &left, &right);

	for (i = 0; i < nlocal; i++)
		elmnts[i] = rand() % (1000 + 1);

	double t1, t2, dt;
	t1 = MPI_Wtime();

	qsort(elmnts, nlocal, sizeof(int), IncOrder);

	if (myrank % 2 == 0) {
		oddrank = myrank - 1;
		evenrank = myrank + 1;
	}
	else {
		oddrank = myrank + 1;
		evenrank = myrank - 1;
	}
	if (oddrank == -1 || oddrank == npes)
		oddrank = MPI_PROC_NULL;
	if (evenrank == -1 || evenrank == npes)
		evenrank = MPI_PROC_NULL;

	for (i = 0; i < npes - 1; i++) {
		if (i % 2 == 1)
			MPI_Sendrecv(elmnts, nlocal, MPI_INT, oddrank, 1, relmnts,
				nlocal, MPI_INT, oddrank, 1, MPI_COMM_WORLD, &status);
		else
			MPI_Sendrecv(elmnts, nlocal, MPI_INT, evenrank, 1, relmnts,
				nlocal, MPI_INT, evenrank, 1, MPI_COMM_WORLD, &status);

		MPI_Neighbor_alltoallw(MPI_BOTTOM, snd_counts, snd_displs, snd_types,
			MPI_BOTTOM, rcv_counts, rcv_displs, rcv_types, new_comm);
		CompareSplit(nlocal, elmnts, relmnts, wspace,
			myrank < status.MPI_SOURCE);

	}
	t2 = MPI_Wtime();
	count << "The sort time: " << t2 - t1 << endl;
	MPI_Finalize();
}

void CompareSplit(int nlocal, int* elmnts, int* relmnts, int* wspace,
	int keepsmall)
{
	int i, j, k;

	for (i = 0; i < nlocal; i++)
		wspace[i] = elmnts[i];

	if (keepsmall) {
		for (i = j = k = 0; k < nlocal; k++) {
			if (j == nlocal || (i < nlocal && wspace[i] < relmnts[j]))
				elmnts[k] = wspace[i++];
			else
				elmnts[k] = relmnts[j++];

		}

	}
	else {
		for (i = k = nlocal - 1, j = nlocal - 1; k >= 0; k--) {
			if (j == 0 || (i >= 0 && wspace[i] >= relmnts[j]))
				elmnts[k] = wspace[i--];
			else
				elmnts[k] = relmnts[j--];

		}

	}
}

int IncOrder(const void* e1, const void* e2)
{
	return (*((int*)e1) - *((int*)e2));
}

