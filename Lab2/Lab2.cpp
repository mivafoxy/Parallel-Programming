#define _CRT_SECURE_NO_WARNINGS
// Lab2.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

static const int ITERATIONS_COUNT[] = { 10, 50, 100, 500, 1000 };

static const int SIMPLE_LU = 1;
static const int ROW_LU = 2;
static const int COLUMN_LU = 3;
static const int GLOBAL_LU = 4;

//  LU-���������� ��� ������ �������� ��������.
void luDecompositionSimple(double** a, int* map, int myrank, int nprocs, int n);

// LU-���������� � �������� �������� �������� �� ������.
void luDecompositionRow(double** a, int* map, int myrank, int nprocs, int n);

// LU-���������� � ������� �������� �������� �� �������.
void luDecompositionColumn(int n);

// ���������� LU-����������.
void luDecompositionGlobal(int n);

int main(int argc, char* argv[])
{
	const int steps = sizeof(ITERATIONS_COUNT) / sizeof(int);

	int myrank, nprocs;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	int chosen_one = -1;
	if (myrank == 0)
	{
		printf("Choose your destiny.\r\n");

		printf("1 - simple LU decomposition.\r\n");
		printf("2 - LU decomposition by rows.\r\n");
		printf("3 - LU decomposition by columns.\r\n");
		printf("4 - global LU decomposition.\r\n");

		scanf("%d", &chosen_one);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	double start;
	double end;
	double finish;

	FILE *f = fopen("res.txt", "w");

	for (int iterationCounts = 0; iterationCounts < steps; iterationCounts++)
	{
		int innerSteps = ITERATIONS_COUNT[iterationCounts];

		printf("Iterations count: %d \n", innerSteps);

		int *map = (int*)malloc(sizeof(int) * innerSteps);
		double **a = new double*[innerSteps];

		// ������������� �������.
		for (int row = 0; row < innerSteps; row++)
		{
			a[row] = new double[innerSteps];

			for (int column = 0; column < innerSteps; column++)
				a[row][column] = 1 / (1.0 * (row + column + 1));
		}

		// ������������� �����.
		for (int i = 0; i < innerSteps; i++)
			map[i] = i % nprocs;

		start = MPI_Wtime();

		if (chosen_one == SIMPLE_LU)
			luDecompositionSimple(a, map, myrank, nprocs, innerSteps);
		else if (chosen_one == ROW_LU)
			luDecompositionRow(a, map, myrank, nprocs, innerSteps);
		else if (chosen_one = COLUMN_LU)
			luDecompositionColumn(innerSteps);
		else if (chosen_one == GLOBAL_LU)
			luDecompositionGlobal(innerSteps);

		end = MPI_Wtime() - start;

		MPI_Reduce(&end, &finish, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		if (myrank == 0)
		{
			// ������� ������ ������� � ����.
			fprintf(f, "%d\t %lg\t \r\n", innerSteps, finish / (1.0*nprocs));

			free(map);
			delete a;
		}
	}

	MPI_Finalize();

	system("pause");
	return 0;
}

//  LU-���������� ��� ������ �������� ��������.
void luDecompositionSimple(double** a, int* map, int myrank, int nprocs, int n)
{
	for (int i = 0; i < (n - 1); i++)
	{
		if (map[i] == myrank)
		{
			for (int j = (i + 1); j < n; j++)
				a[i][j] /= a[i][i];
		}

		MPI_Bcast(&a[i][i + 1], (n - i - 1), MPI_DOUBLE, map[i], MPI_COMM_WORLD);

		for (int j = (i + 1); j < n; j++)
		{
			if (map[j] == myrank)
			{
				for (int k = i + 1; k < n; k++)
					a[j][k] -= a[j][i] * a[i][k];
			}
		}
	}
}

// LU-���������� � �������� �������� �������� �� ������.
void luDecompositionRow(double** a, int* map, int myrank, int nprocs, int n)
{
	for (int k = 0; k < (n - 1); k++)
	{
		int maxElementIndex = k;

		if (map[k] == myrank)
		{
			double max = a[k][maxElementIndex];

			for (int column = (k + 1); column < n; column++) // ����� ������������� �������� � ������� �� ������� � ��� �������.
			{
				if (map[column] == myrank && a[k][column] > max)
				{
					max = a[k][column];
					maxElementIndex = column;
				}
			}

			if (k != maxElementIndex) // ������������ ���������� ������������� �������� �����.
			{
				for (int row = 0; row < n; row++)
				{
					double temp = a[row][k];
					a[row][k] = a[row][maxElementIndex];
					a[row][maxElementIndex] = temp;
				}
			}

			for (int column = (k + 1); column < n; column++)
				a[k][column] /= a[k][k];
		}

		MPI_Status status;

		for (int row = 0; row < n; row++) // ��������� ���� ��������� ������� � ���, ��� � ��� ��� ������ ��� - �� ����������.
		{
			if (map[k] == myrank)
			{
				if (map[row] != myrank)
				{
					MPI_Send(&a[row][k], 1, MPI_DOUBLE, map[row], 1, MPI_COMM_WORLD);
					MPI_Send(&maxElementIndex, 1, MPI_INT, map[row], 1, MPI_COMM_WORLD);
					MPI_Send(&a[row][maxElementIndex], 1, MPI_DOUBLE, map[row], 1, MPI_COMM_WORLD);
				}
			}
			else
			{
				if (map[row] == myrank)
				{
					int recvIndex;

					MPI_Recv(&a[row][k], 1, MPI_DOUBLE, map[k], 1, MPI_COMM_WORLD, &status);
					MPI_Recv(&recvIndex, 1, MPI_INT, map[k], 1, MPI_COMM_WORLD, &status);
					MPI_Recv(&a[row][recvIndex], 1, MPI_DOUBLE, map[k], 1, MPI_COMM_WORLD, &status);
				}
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);

		for (int rowTemp = (k + 1); rowTemp < n; rowTemp++)
		{
			if (map[rowTemp] == myrank)
			{
				for (int column = (k + 1); column < n; column++)
					a[k][column] -= a[k][rowTemp] * a[k][column];
			}
		}
	}
}

// LU-���������� � ������� �������� �������� �� �������.
void luDecompositionColumn(int n)
{

}

// ���������� LU-����������.
void luDecompositionGlobal(int n)
{

}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
