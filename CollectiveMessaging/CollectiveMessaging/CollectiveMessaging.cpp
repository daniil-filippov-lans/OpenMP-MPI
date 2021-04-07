#include <mpi.h>
#include <iostream>

using namespace std;

// lab 1
void Messaging(int M)
{
	int ProcNum, ProcRank, RecvRank, currentLoop = M + 1;
	MPI_Status Status;

	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	while (M > 0) 
	{
		if (ProcRank == 0)
		{
			cout << "Loop: " << currentLoop - M << endl;
			MPI_Send(&ProcRank, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
			MPI_Recv(&RecvRank, 1, MPI_INT, ProcNum - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
			cout << "#1 process: " << ProcRank << " from: " << RecvRank << endl;
		}
		else 
		{
			MPI_Recv(&RecvRank, 1, MPI_INT, ProcRank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
			cout << "Process: " << ProcRank << " from: " << RecvRank << endl;

			if (ProcRank + 1 != ProcNum)
				MPI_Send(&ProcRank, 1, MPI_INT, ProcRank + 1, 0, MPI_COMM_WORLD);
			else
				MPI_Send(&ProcRank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
		M--;
	}
}
// lab 2
void CollectiveMessaging(int M)
{
	int ProcNum, ProcRank, RecvRank, SendRank, currentLoop = M + 1, message, totalSum = 0;
	MPI_Status Status;

	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	while (M != 0)
	{
		if (ProcRank == 0)
		{
			cout << "Loop: " << currentLoop - M << endl;
			SendRank = ProcRank + 1;
			RecvRank = ProcNum - 1;
			MPI_Bcast(&SendRank, 1, MPI_INT, ProcRank, MPI_COMM_WORLD);
			cout << "Senior process: " << ProcRank << " recieve message from: " << RecvRank << endl;
		}
		else
		{
			if (ProcRank == ProcNum - 1)
			{
				SendRank = 0;
				RecvRank = ProcRank - 1;
				MPI_Bcast(&SendRank, 1, MPI_INT, ProcRank, MPI_COMM_WORLD);
				cout << "Junior process: " << ProcRank << " recieve message from: " << RecvRank << endl;
			}
			else
			{
				SendRank = ProcRank - 1;
				RecvRank = ProcRank + 1;
				MPI_Bcast(&SendRank, 1, MPI_INT, ProcRank, MPI_COMM_WORLD);
				cout << "Junior process: " << ProcRank << " recieve message from: " << SendRank << endl;
			}
		}
		M--;
		MPI_Barrier(MPI_COMM_WORLD);
	}
	if (ProcRank == 0)
	{
		cout << "End " << currentLoop << " loop of " << endl;
	}
}

// lab 3
void GeneralizedMessaging(int M)
{
	int ProcNum, ProcRank, RecvRank, SendRank, currentLoop = M + 1, message, totalSum = 0;
	MPI_Status Status;

	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	int mas[4];
	int mas2[4];

	while (M != 0)
	{
		if (ProcRank == 0)
		{
			for (int i = 0; i < 4; i++)
			{
				mas[i] = i;
			}
			cout << "Loop: " << currentLoop - M << endl;
			SendRank = ProcRank + 1;
			RecvRank = ProcNum - 1;
			MPI_Scatter(&mas, 1, MPI_INT, &mas2, 1, MPI_INT, ProcRank, MPI_COMM_WORLD);
			//cout << "Senior process: " << ProcRank << " recieve message from: " << RecvRank << endl;
		}

		else
		{
			if (ProcRank == ProcNum - 1)
			{
				SendRank = 0;
				RecvRank = ProcRank - 1;
				MPI_Scatter(&SendRank, 1, MPI_INT, &RecvRank, 1, MPI_INT, ProcRank, MPI_COMM_WORLD);
				cout << "Junior process: " << ProcRank << " recieve message from: " << RecvRank << endl;
			}
			else
			{
				SendRank = ProcRank - 1;
				RecvRank = ProcRank + 1;
				MPI_Scatter(&SendRank, 1, MPI_INT, &RecvRank, 1, MPI_INT, ProcRank, MPI_COMM_WORLD);
				cout << "Junior process: " << ProcRank << " recieve message from: " << SendRank << endl;
			}
			M--;
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}

	if (ProcRank == 0) 
	{
		cout << "End " << currentLoop << " loop of " << endl;
	}
}

// lab 3 dop matrix
void GeneralizedMessagingMatrix(int M)
{
	int ProcNum, ProcRank, RecvRank, SendRank, currentLoop = M + 1, message, totalSum = 0;
	MPI_Status Status;

	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	int mas[4];
	int mas2[4];

	while (M != 0) 
	{
		if (ProcRank == 0) 
		{
			for (int i = 0; i < 4; i++)
			{
				mas[i] = i;
			}

			cout << "Loop: " << currentLoop - M << endl;
			SendRank = ProcRank + 1;
			RecvRank = ProcNum - 1;
			MPI_Scatter(&mas, 1, MPI_INT, &mas2, 1, MPI_INT, ProcRank, MPI_COMM_WORLD);
			//cout << "Senior process: " << ProcRank << " recieve message from: " << RecvRank << endl;
		}
		else
		{

		}

		if (ProcRank == 0)
		{
			cout << "Jun:" << "mas 1: ";
			for (int i = 0; i < 4; i++)
			{
				cout << mas[i] << " ";
			}
			cout << " mas 2: ";
			for (int i = 0; i < 4; i++)
			{
				cout << mas2[i] << " ";
			}
		}
		M--;
		MPI_Barrier(MPI_COMM_WORLD);
	}

	if (ProcRank == 0) 
	{
		cout << "End " << currentLoop << " loop of " << endl;
	}
}

int main(int* argc, char** argv) {
	int M = 4;

	MPI_Init(argc, &argv);

	//Messaging(M);
	//CollectiveMessaging(M);
	//GeneralizedMessaging(M);
	GeneralizedMessagingMatrix(M);

	MPI_Finalize();

	return 0;
}
