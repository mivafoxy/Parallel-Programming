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

//  LU-разложение без выбора ведущего элемента.
void luDecompositionSimple(double** a, int* map, int myrank, int nprocs, int n);

// LU-разложение с выобором ведущего элемента по строке.
void luDecompositionRow(int n);

// LU-разложение с выбором ведущего элемента по столбцу.
void luDecompositionColumn(int n);

// Глобальное LU-разложение.
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
		int *map = (int*)malloc(sizeof(int) * innerSteps);
		double **a = new double*[innerSteps];

		// Инициализация матрицы.
		for (int row = 0; row < innerSteps; row++)
		{
			a[row] = new double[innerSteps];
			
			for (int column = 0; column < innerSteps; column++)
				a[row][column] = 1 / (1.0 * (row + column + 1));
		}

		// Инициализация карты.
		for (int i = 0; i < innerSteps; i++)
			map[i] = i % nprocs;

		start = MPI_Wtime();

		if (chosen_one == SIMPLE_LU)
			luDecompositionSimple(a, map, myrank, nprocs, innerSteps);
		else if (chosen_one == ROW_LU)
			luDecompositionRow(innerSteps);
		else if (chosen_one = COLUMN_LU)
			luDecompositionColumn(innerSteps);
		else if (chosen_one == GLOBAL_LU)
			luDecompositionGlobal(innerSteps);

		end = MPI_Wtime() - start;

		MPI_Reduce(&end, &finish, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		if (myrank == 0)
		{
			// Сделать запись времени в файл.
			fprintf(f, "%d\t %lg\t \r\n", innerSteps, finish / (1.0*nprocs));

			free(map);
			delete a;
		}
	}

	MPI_Finalize();

	system("pause");
	return 0;
}

//  LU-разложение без выбора ведущего элемента.
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

// LU-разложение с выобором ведущего элемента по строке.
void luDecompositionRow(int n)
{

}

// LU-разложение с выбором ведущего элемента по столбцу.
void luDecompositionColumn(int n)
{

}

// Глобальное LU-разложение.
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
