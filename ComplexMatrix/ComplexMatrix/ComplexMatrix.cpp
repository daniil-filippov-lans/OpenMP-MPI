#define MSMPI_NO_DEPRECATE_20
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mpi.h"
#include <clocale>
#include <iostream>

using namespace std;

int randInt() {
	return 1 + rand() % 5;
}
struct Complex {
	int Re;
	int Im;
};
Complex MultiplicationComplex(Complex a, Complex b) {
	Complex c;
	c.Re = a.Re * b.Re - a.Im * b.Im;
	c.Im = a.Re * b.Im + b.Re * a.Im;
	return c;
}

Complex* initMatrix(int N) {
	Complex* M = new Complex [N * N];
	return M;
}
void randMatrix(Complex* A, int N) {
	for (int i = 0; i < N * N; i++)	{
		A[i].Re = randInt();
		A[i].Im = randInt();
	}
}
void printMatrix(Complex* A, int N) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			cout << "(" << A[i * N + j].Re << ", " << A[i * N + j].Im << ") ";
		}		
	}
}
void MultiplicationMatrix(Complex* A, Complex* B, Complex* C, int N) {
	Complex Sum, Mult;
	Sum.Re = 0;
	Sum.Im = 0;

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++)	{
			Sum.Re = 0;
			Sum.Im = 0;
			for (int k = 0; k < N; k++) {
				Mult = MultiplicationComplex(A[i * N + k], B[k * N + j]);
				Sum.Re += Mult.Re;
				Sum.Im += Mult.Im;
			}
			C[i * N + j].Re = Sum.Re;
			C[i * N + j].Im = Sum.Im;
		}
	}	
}

int main(int argc, char* argv[]) {
	int A = atoi(argv[1]); int N = atoi(argv[2]); // À - matrix amount; N - matrix size 

	int pair = A / 2 + 1;
	bool even = 1;
	if (A % 2 == 1)
		even = 0;

	Complex* matrixA = initMatrix(N); Complex* matrixB = initMatrix(N);	Complex* matrixC = initMatrix(N);
	Complex* arrayMatrix = new Complex[N * N * pair];

	int ProcNum, ProcRank;
	MPI_Status Status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (ProcNum < pair) {
		MPI_Finalize();
		delete[] matrixA; delete[] matrixB; delete[] matrixC;
		delete[] arrayMatrix;
		return 0;
	}

	int const count = 2;

	int blocklens[count]; blocklens[0] = 1;	blocklens[1] = 1;

	MPI_Aint indices[count];
	Complex ComplexNum;
	MPI_Address(&ComplexNum.Re, &indices[0]);
	MPI_Address(&ComplexNum.Im, &indices[1]);
	indices[1] = indices[1] - indices[0];
	indices[0] = 0;

	MPI_Datatype old_types[count]; old_types[0] = MPI_INT; old_types[1] = MPI_INT;
	MPI_Datatype MPI_INT_COMPLEX_TYPE;

	MPI_Type_struct(count, blocklens, indices, old_types, &MPI_INT_COMPLEX_TYPE);
	MPI_Type_commit(&MPI_INT_COMPLEX_TYPE);

	MPI_Group group;
	MPI_Comm_group(MPI_COMM_WORLD, &group);

	int* ranks = new int[pair];
	for (int i = 0; i < pair; i++) {
		ranks[i] = i;
	}
	
	MPI_Group complex_group;
	MPI_Group_incl(group, pair, ranks, &complex_group);

	MPI_Comm COMM;
	MPI_Comm_create(MPI_COMM_WORLD, complex_group, &COMM);

	if (ProcRank == 0) {
		cout << "Matrix size: " << N << "x" << N << " amount multiplication: " << A << endl;
		randMatrix(matrixA, N);
		printMatrix(matrixA, N);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Bcast(matrixA, N * N, MPI_INT_COMPLEX_TYPE, 0, MPI_COMM_WORLD);

	if (ProcRank >= 1 && ProcRank < pair) {
		MultiplicationMatrix(matrixA, matrixA, matrixB, N);
		cout << "  Rank " << ProcRank << " - multi done"<< endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if (ProcRank >= 0 && ProcRank < pair) {
		MPI_Gather(matrixB, N * N, MPI_INT_COMPLEX_TYPE, arrayMatrix, N * N, MPI_INT_COMPLEX_TYPE, 0, COMM);
	}

	if (ProcRank == 0) {
		if (even) {
			for (int i = 0; i < N * N; i++)	{
				matrixA[i].Re = arrayMatrix[N * N + i].Re;
				matrixA[i].Im = arrayMatrix[N * N + i].Im;
			}
		} 

		for (int i = even + 1; i < pair; i++) {
			for (int j = 0; j < N * N; j++) {
				matrixB[j].Re = arrayMatrix[i * N * N + j].Re;
				matrixB[j].Im = arrayMatrix[i * N * N + j].Im;
			}
			MultiplicationMatrix(matrixA, matrixB, matrixC, N);
			for (int j = 0; j < N * N; j++) {
				matrixA[j].Re = matrixC[j].Re;
				matrixA[j].Im = matrixC[j].Im;
			}
		}

		cout << "Result Matrix:" << endl;
		printMatrix(matrixA, N);

		MPI_Type_free(&MPI_INT_COMPLEX_TYPE);
		MPI_Group_free(&group);
		MPI_Group_free(&complex_group);
		MPI_Comm_free(&COMM);
	}

	MPI_Finalize();
	delete[] matrixA; delete[] matrixB; delete[] matrixC;
	delete[] arrayMatrix;
	delete[] ranks;
	return 0;
}