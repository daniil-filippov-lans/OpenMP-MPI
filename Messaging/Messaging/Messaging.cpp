#include <mpi.h>
#include <iostream>

using namespace std;

int main(int* argc, char** argv) {
	//cmd: mpiexec -n 4 Messaging.exe
	//cd: C:\Users\danim\Desktop\remote\6\MPI\OpenMP-MPI\Messaging\Debug

	int ProcNum, ProcRank, RecvRank, M = 4, currentLoop = M + 1;
	MPI_Status Status;


	MPI_Init(argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	while (M != 0) {
		if (ProcRank == 0) {

			cout << "Loop: " << currentLoop - M << endl;
			MPI_Send(&ProcRank, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
			MPI_Recv(&RecvRank, 1, MPI_INT, ProcNum - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
			cout << "Senior process: " << ProcRank << " recieve message from: " << RecvRank << endl;

			M--;
		}
		else {
			MPI_Recv(&RecvRank, 1, MPI_INT, ProcRank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
			cout << "Junior process: " << ProcRank << " recieve message from: " << RecvRank << endl;

			if (ProcRank + 1 != ProcNum)
				MPI_Send(&ProcRank, 1, MPI_INT, ProcRank + 1, 0, MPI_COMM_WORLD);
			else
				MPI_Send(&ProcRank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
	}

	MPI_Finalize();

	return 0;
}
