#include <iostream>
#include <math.h>
#include <conio.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <cmath>
#include <locale.h>

using namespace std;

struct Interval
{
	double start;
	double finish;
};

struct Approximation
{
	double Sm;
	double Sm_next;

	int m;
	int n;

	double number;

	void aprPow(double number, int m)
	{
		n *= number;
	}
};

struct NewtonCotesCoef
{
	int k;
	double c0;
	int omega[];

	void setOrder(int k)
	{
		switch (k)
		{
		case 3:
		{
				  this->k = 3;
				  c0 = 3.0 / 8.0;
				  omega[0] = 1;
				  omega[1] = 3;
				  omega[2] = 3;
				  omega[3] = 1;
				  break;
		}
		case 5:
		{
				  this->k = 5;
				  c0 = 5.0 / 288.0;
				  omega[0] = 19;
				  omega[1] = 75;
				  omega[2] = 50;
				  omega[3] = 50;
				  omega[4] = 75;
				  omega[5] = 19*2;
				  break;
		}
		default:
		{
				   this->k = 5;
				   break;
		}
		}
	}
};

NewtonCotesCoef coef;

double f(double x)
{
	return pow(x+1, 2) / (x*x*x*sqrt(2+x));
}

double IntegralMethodOfMediumRectangles(double a, double b, double n)
{
	double integral = 0;
	double h = (b - a) / n;

	/*
	for (double i = a; i < b; i += h)
	{
	integral += h * f(i + h/2);
	}
	*/

	double x = a;
	for (int i = 0; i < n; i++)
	{
		integral += h * f((2 * x + h) / 2);
		x += h;
	}


	return integral;
}

double IntegralMethodOfTrapezoid(double a, double b, double n)
{
	double integral = 0;
	double h = (b - a) / n;

	/*
	for (double i = a + h; i < b; i += h)
	{
	integral += f(i);
	}
	*/

	double x = a + h;
	for (int i = 1; i < n; i++)
	{
		integral += f(x);
		x += h;
	}

	integral *= 2;
	integral += f(a) + f(b);
	integral *= h / 2;

	return integral;
}

double IntegralMethodOfNewtonCotes(double a, double b, double n)
{
	double integral = 0;
	double h = (b - a) / n;
	double x = a;
	
	integral += coef.omega[0]*f(x);
	x += h;

	for (int i = h; i < n; i++)
	{
		integral += coef.omega[(i % coef.k) + 1] * f(x);
		x += h;
	}
	
	integral -= coef.omega[0] * f(b);
	integral *= h * coef.c0;
	
	return integral;
}


void solveIntegral(int number, int m, Interval i, double epsilon, double(*funcMethod)(double, double, double))
{
	Approximation apr;

	apr.number = number;
	apr.m = m;
	apr.n = pow(apr.number, apr.m);

	do
	{
		apr.m++;
		apr.aprPow(apr.number, apr.m);
		apr.Sm = funcMethod(i.start, i.finish, apr.n);

		apr.m++;
		apr.aprPow(apr.number, apr.m);
		apr.Sm_next = funcMethod(i.start, i.finish, apr.n);

	} while (abs(apr.Sm - apr.Sm_next) >= epsilon);

	cout
		<<
		"Точность: " << epsilon << endl <<
		"Число разбиений: " << apr.n << endl <<
		"Значение интеграла: " << apr.Sm << endl
		<<
	endl;
}

void main()
{
	setlocale(LC_ALL, "Russian");

	Interval i{ 2, 5 };
	double epsilon = 0.01;
	coef.setOrder(5);

	cout << "Метод средних прямоугольников. " << endl; solveIntegral(2, 0, i, epsilon, &IntegralMethodOfMediumRectangles);
	cout << "Метод трапеций. " << endl; solveIntegral(2, 0, i, epsilon, &IntegralMethodOfTrapezoid);
	cout << "Метод Ньютона-Котеса 5-го порядка. " << endl; solveIntegral(coef.k, 1, i, epsilon, &IntegralMethodOfNewtonCotes);

	_getch();
}