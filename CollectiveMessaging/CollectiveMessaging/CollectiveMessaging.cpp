#include <mpi.h>
#include <iostream>

using namespace std;

int main(int* argc, char** argv) {
	int ProcNum, ProcRank, RecvRank, M = 4, currentLoop = M + 1, message, totalSum = 0;
	MPI_Status Status;


	MPI_Init(argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	while (M != 0) {
		if (ProcRank == 0) {

			cout << "Loop: " << currentLoop - M << endl;
			message = ProcRank;
			MPI_Send(&message, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
			MPI_Recv(&RecvRank, 1, MPI_INT, ProcNum - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
			cout << "Senior process: " << ProcRank << " recieve message from: " << RecvRank << endl;

			M--;
		}
		else {
			message = ProcRank;
			MPI_Recv(&RecvRank, 1, MPI_INT, ProcRank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
			cout << "Junior process: " << ProcRank << " recieve message from: " << RecvRank << endl;

			if (ProcRank + 1 != ProcNum)
				MPI_Send(&message, 1, MPI_INT, ProcRank + 1, 0, MPI_COMM_WORLD);
			else
				MPI_Send(&message, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
		MPI_Reduce(&message, &totalSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	}


	MPI_Barrier(MPI_COMM_WORLD);

	if (ProcRank == 0) {
		cout << "Total sum " << M << " loop of process: " << totalSum << endl;
	}

	MPI_Finalize();

	return 0;
}