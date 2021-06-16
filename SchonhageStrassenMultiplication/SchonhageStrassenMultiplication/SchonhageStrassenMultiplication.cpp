#include <iostream>
#include <vector>
#include <complex>
#include <mpi.h>
#include <windows.h>

using namespace std;

const double PI = acos(-1);

// fast fourier transform
void FFT(vector<complex<double>>& A, bool back) {
	int N = (int)A.size();

	if (N == 1)  return;

	vector<complex<double>> a0(N / 2), a1(N / 2);
	for (int i = 0, j = 0; i < N; i += 2, j++) {
		a0[j] = A[i];
		a1[j] = A[i + 1];
	}

	FFT(a0, back);
	FFT(a1, back);

	double ang;
	if (back) {
		ang = 2 * PI / N * (-1);
	}
	else {
		ang = 2 * PI / N * 1;
	}

	complex<double> w(1);
	complex<double> wn(cos(ang), sin(ang));
	for (int i = 0; i < N / 2; i++) {
		A[i] = a0[i] + w * a1[i];
		A[i + N / 2] = a0[i] - w * a1[i];
		if (back)
			A[i] /= 2, A[i + N / 2] /= 2;
		w *= wn;
	}
}

void multiply(vector<int> A, vector<int> B, vector<int>& result) {
	vector<complex<double>> FA = vector<complex<double>>(A.begin(), A.end());
	vector<complex<double>> FB = vector<complex<double>>(B.begin(), B.end());
	int n = 1;

	while (n < max(A.size(), B.size())) {
		n *= 2;
	}
	n *= 2;

	while (FA.size() != n) {
		FA.emplace(FA.begin(), complex<double>(0, 0));
	}
	while (FB.size() != n) {
		FB.emplace(FB.begin(), complex<double>(0, 0));
	}

	FFT(FA, false);
	FFT(FB, false);

	for (int i = 0; i < n; i++)
	{
		FA[i] *= FB[i];
	}

	FFT(FA, true);

	result.resize(n, 0);

	for (int i = 0; i < n - 1; i++)
	{
		result[i + 1] = int(FA[i].real() + 0.5);
	}

	for (int i = result.size() - 1; i > 0; i--) {
		if (result[i] >= 10) {
			result[i - 1] += result[i] / 10;
			result[i] = result[i] % 10;
		}
	}

}

int main(int* argc, char** argv) {
	MPI_Status status;
	MPI_Init(argc, &argv);
	int ProcNum, ProcRank;

	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	int INDEX[] = { 2,5,8,9,10,11,12 };
	int EDGES[] = { 1,2,3,4,0,5,6,0,1,1,2,2 };

	MPI_Comm F_COMM;
	MPI_Graph_create(MPI_COMM_WORLD, 7, INDEX, EDGES, 1, &F_COMM);

	int neighborsCount;
	MPI_Graph_neighbors_count(F_COMM, ProcRank, &neighborsCount);

	int* neighbors = new int[neighborsCount];
	MPI_Graph_neighbors(F_COMM, ProcRank, neighborsCount, neighbors);

	int long_num_size = 2;
	MPI_Datatype MPI_LONG_NUM;
	
	MPI_Type_contiguous(long_num_size * ProcNum, MPI_INT, &MPI_LONG_NUM);
	MPI_Type_commit(&MPI_LONG_NUM);

	if (neighborsCount == 1)
	{
		vector<int> A_vector = { 1, ProcRank };
		vector<int> B_vector = { ProcRank, 1 };
		vector<int> result;

		multiply(A_vector, B_vector, result);

		for (vector<int>::iterator i = result.begin(); i != result.end(); ++i) {
			printf("%d ", *i);
		}

		MPI_Send(result.data(), 1, MPI_LONG_NUM, neighbors[0], ProcRank, F_COMM);
	}

	if (neighborsCount == 3)
	{
		long_num_size *= 2;
		int* A = new int[long_num_size];
		int* B = new int[long_num_size];

		if (ProcRank == 1) {
			MPI_Recv(A, 1, MPI_LONG_NUM, neighbors[0], 3, F_COMM, &status);
			MPI_Recv(B, 1, MPI_LONG_NUM, neighbors[1], 4, F_COMM, &status);
		}
		if (ProcRank == 2) {
			MPI_Recv(A, 1, MPI_LONG_NUM, neighbors[0], 5, F_COMM, &status);
			MPI_Recv(B, 1, MPI_LONG_NUM, neighbors[1], 6, F_COMM, &status);
		}

		vector<int> A_vector(A, A + long_num_size);
		vector<int> B_vector(B, B + long_num_size);
		vector<int> result(long_num_size * 2, 0);

		multiply(A_vector, B_vector, result);

		for (vector<int>::iterator i = result.begin(); i != result.end(); ++i) {
			printf("%d ", *i);
		}

		MPI_Send(result.data(), 1, MPI_LONG_NUM, neighbors[2], ProcRank, F_COMM);
		MPI_Send(result.data(), 1, MPI_LONG_NUM, neighbors[2], ProcRank, F_COMM);
	}

	if (neighborsCount == 2)
	{
		long_num_size *= 4;
		int* A = new int[long_num_size];
		int* B = new int[long_num_size];

		MPI_Recv(A, 1, MPI_LONG_NUM, neighbors[0], 1, F_COMM, &status);
		MPI_Recv(B, 1, MPI_LONG_NUM, neighbors[1], 2, F_COMM, &status);

		vector<int> A_vector(A, A + long_num_size);
		vector<int> B_vector(B, B + long_num_size);
		vector<int> result(long_num_size * 2, 0);

		multiply(A_vector, B_vector, result);

		for (vector<int>::iterator i = result.begin(); i != result.end(); ++i) {
			printf("%d ", *i);
		}
	}

	MPI_Type_free(&MPI_LONG_NUM);
	MPI_Comm_free(&F_COMM);

	MPI_Finalize();
}