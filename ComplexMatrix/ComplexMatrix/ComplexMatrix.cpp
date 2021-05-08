#include <math.h>;
#include <stdio.h>;
#include <stdlib.h>;
#include <complex>;
#include "mpi.h";

using namespace std;

int main(int argc, char* argv[]) {

	double x[100], TotalSum, ProcSum = 0.0;
	int ProcRank, ProcNum, N = 100;
	MPI_Status Status;
	complex<double> z(1.0, 2.0);

	// initialization 
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	// preparation
	//if (ProcRank == 0) DataInitialization(x, N);

	// send data to all proc
	MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// calculating the partial sum of each process
	// on each proc sum vec from i1 to i2

	int k = N / ProcNum;
	int i1 = k * ProcRank;
	int i2 = k * (ProcRank + 1);
	if (ProcRank == ProcNum - 1) i2 = N;
	for (int i = i1; i < i2; i++)
		ProcSum = ProcSum + x[i];
	// recive part sum in proc rank 0
	if (ProcRank == 0) {
		TotalSum = ProcSum;
		for (int i = 1; i < ProcNum; i++) {
			MPI_Recv(&ProcSum, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,
				&Status);
			TotalSum = TotalSum + ProcSum;
		}
	}
	else // all proc send sum to 0 
		MPI_Send(&ProcSum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	// result
	if (ProcRank == 0)
		printf("\nTotal Sum = %10.2f", TotalSum);
	MPI_Finalize();

	return 0;
}