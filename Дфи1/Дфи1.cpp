#define _CRT_SECURE_NO_WARNINGS

#include "pch.h"
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"
#include "math.h"

static const int ITERATIONS_COUNT[] = { 10, 50, 100, 500, 1000 };

static const int RIEMAN_METHOD = 1;
static const int TRAPEZOIDAL_RULE = 2;
static const int SIMPSON_RULE = 3;

static const int PI = 3.141592653589793;

// Получить количество "столбцов интегрирования". - но это неточно. Однако это 100% связано со столбцами.
double get_h(int n)
{
	return (1.0 / n);
}

// arctan(1.0) = Pi/4 => Pi = 4 * integral(arctan).
double calculate_pi(int n, double sum)
{
	return (4 * get_h(n)*sum);
}

// Дифференциал функции, по которой вообще эту PI все ищут.
double dx_arctan(double x)
{
	return (1/(1 + x*x));
}

// Квадратурная формула прямоугольника.
double rieman_sum(int myrank, int nprocs, int n)
{
	double sum = 0;
	double h = get_h(n);

	for (int rank = myrank + 1; rank <= n; rank += nprocs) // Считаем и суммируем каждый столбец интеграла по формуле прямоугольника.
		sum += dx_arctan(h*(rank - 0.5));

	return sum;
}

// Квадратурная формула трапеций.
double trapeziodal_rule(int myrank, int nprocs, int n)
{
	double sum = 0;
	double h = get_h(n);
	
	for (int rank = myrank + 1; rank <= n - 1; rank += nprocs) // Считаем и суммируем каждый столбец интеграла по формуле трапеций.
		sum += dx_arctan(rank*h);

	int a = 0;
	int b = 1;

	if (myrank == 0)
		sum += (dx_arctan(a) + dx_arctan(b)) / 2.0;

	return sum;
}

// Формула Симпсона.
double simpsons_rule(int myrank, int nprocs, int n)
{
	double sum = 0;
	double h = get_h(n);

	sum = 2 * rieman_sum(myrank, nprocs, n) + trapeziodal_rule(myrank, nprocs, n); // Формула Симсона - это тупо сумма формул трапеции и прямоугольника. См. Теорию.

	return sum;
}

// Всё описание применений функций MPI смотри в теории. Там точно расписано.
int main(int argc, char* argv[])
{
	double pi, sum = 0;
	int myrank, nprocs;

	int chosen_one = 0;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	if (myrank == 0)
	{
		printf("Choose your method: \r\n");
		printf("1 = reiman sum\r\n");
		printf("2 = Trapezoidal rule \r\n");
		printf("3 = Simpson's rule \r\n");


		scanf("%d", &chosen_one);
	}

	if (chosen_one == TRAPEZOIDAL_RULE)
		printf("Sides are 0 and 1 because there is only true way to calculate PI!!! \n");

	double start = 0;
	double end = 0;
	double finish = -1;

	FILE *f = fopen("results.txt", "w");

	for (int iterationsCount = 0; iterationsCount < (sizeof(ITERATIONS_COUNT) / sizeof(int)); iterationsCount++)
	{
		int n = ITERATIONS_COUNT[iterationsCount];

		MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

		printf("Iterations count: %d\n", n);

		start = MPI_Wtime();
		if (chosen_one == RIEMAN_METHOD)
		{
			double riemanSum = rieman_sum(myrank, nprocs, n);
			sum = calculate_pi(n, riemanSum);
		}
		else if (chosen_one == TRAPEZOIDAL_RULE)
		{
			double trapeziodalSum = trapeziodal_rule(myrank, nprocs, n);
			sum = calculate_pi(n, trapeziodalSum);
		}
		else if (chosen_one == SIMPSON_RULE)
		{
			double simpsonSum = simpsons_rule(myrank, nprocs, n);
			sum = calculate_pi(n, simpsonSum) / 3.0;
		}
		end = MPI_Wtime() - start;

		MPI_Reduce(&sum, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&end, &finish, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		if (myrank == 0)
		{
			printf("Computed value of pi = %lg\n", pi);
			fprintf(f, "%d\t %.18lg\t %lg\t %lg\n", n, pi, finish, fabs(pi - PI));
		}
	}

	MPI_Finalize();

	system("pause");
	return 0;
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
