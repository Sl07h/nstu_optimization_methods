#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
using namespace std;
typedef double real;



/* Заданная функция
 * x - аргумент функции
 */
real f(real x) {
	return (x - 2)*(x - 2);
}


/* Вычисление n-го числа Фибоначчи
 * n - номер числа
 */
int F(real n) {
	int tmp, prev = 1, current = 1;
	if (n == 1 || n == 2)
		return 1;

	for (int i = 0; i < n - 2; ++i) {
		tmp = current;
		current += prev;
		prev = tmp;
	}

	return current;
}


/* Метод дихотомии
 * a,b - границы отрезка
 * E - точность
 * fileName - файл, в который будет сохранён результат
 * return количество вычислений целевой функции
 */
int calcByDichotomy(real a, real b, real E, const char *fileName) {

	ofstream fout(fileName);
	real aPrev, bPrev, x1, x2, f1, f2, delta = E / 2;
	int k = 1;

	fout << setprecision(numeric_limits<real>::digits10 + 1) << fixed;

	fout << "k \t a \t b \t (b - a) \t (bPrev - aPrev) / (b - a) \t x1 \t x2 \t f1 \t f2 " << endl;
	while (abs(b - a) >= E) {

		fout << k << "\t" << a << "\t" << b << "\t" << b - a << "\t";

		aPrev = a;
		bPrev = b;

		x1 = (a + b - delta) / 2;
		x2 = (a + b + delta) / 2;

		f1 = f(x1);
		f2 = f(x2);

		if (f1 > f2) { // [x1, b]
			a = x1;
		}
		else { // [a, x2]
			b = x2;
		}
		k++;
		fout << (bPrev - aPrev) / (b - a) << "\t" << x1 << "\t" << x2 << "\t" << f1 << "\t" << f2 << endl;
	}
	fout.close();
	return k * 2;
}


/* Метод золотого сечения
 * a,b - границы отрезка
 * E - точность
 * fileName - файл, в который будет сохранён результат
 * return количество вычислений целевой функции
 */
int calcByGoldenSectionSearch(real a, real b, real E, const char *fileName) {

	ofstream fout(fileName);
	real aPrev, bPrev, x1, x2, f1, f2, delta = E / 2;
	real alpha = (3 - sqrt(5)) / 2;
	real beta = (sqrt(5) - 1) / 2;

	int k = 1;

	fout << setprecision(numeric_limits<real>::digits10 + 1) << fixed;

	x1 = a + alpha * (b - a);
	x2 = a + beta * (b - a);
	f1 = f(x1);
	f2 = f(x2);

	fout << "k \t a \t b \t (b - a) \t (bPrev - aPrev) / (b - a) \t x1 \t x2 \t f1 \t f2 " << endl;
	while (abs(b - a) >= E) {

		fout << k << "\t" << a << "\t" << b << "\t" << b - a << "\t";
		aPrev = a;
		bPrev = b;

		if (f1 > f2) { // [x1, b]
			a = x1;

			x1 = x2;
			f1 = f2;

			x2 = a + beta * (b - a);
			f2 = f(x2);
		}
		else { // [a, x2]
			b = x2;

			x2 = x1;
			f2 = f1;

			x1 = a + alpha * (b - a);
			f1 = f(x1);
		}
		k++;
		fout << (bPrev - aPrev) / (b - a) << "\t" << x1 << "\t" << x2 << "\t" << f1 << "\t" << f2 << endl;
	}
	fout.close();
	return k + 1;
}



/* Метод Фибоначчи
 * a,b - границы отрезка
 * E - точность
 * fileName - файл, в который будет сохранён результат
 * return количество вычислений целевой функции
 */
int calcByFibonacci(real a, real b, real E, const char *fileName) {

	ofstream fout(fileName);
	real aPrev, bPrev, length = b - a, x1, x2, f1, f2, delta = E / 2;
	int k = 1, n = 0;

	while (F(n + 2) < (b - a) / E) {
		++n;
	}

	fout << setprecision(numeric_limits<real>::digits10 + 1) << fixed;

	x1 = a + F(n) * (b - a) / F(n + 2);
	x2 = a + F(n + 1) * (b - a) / F(n + 2);
	f1 = f(x1);
	f2 = f(x2);
	real Fn2 = F(n + 2);

	fout << "k \t a \t b \t (b - a) \t (bPrev - aPrev) / (b - a) \t x1 \t x2 \t f1 \t f2 " << endl;
	while (abs(b - a) >= E) {

		fout << k << "\t" << a << "\t" << b << "\t" << b - a << "\t";
		aPrev = a;
		bPrev = b;

		if (f1 > f2) { // [x1, b]
			a = x1;

			x1 = x2;
			f1 = f2;

			x2 = a + F(n - k + 2) * length / Fn2;
			f2 = f(x2);
		}
		else { // [a, x2]
			b = x2;

			x2 = x1;
			f2 = f1;

			x1 = a + F(n - k + 1) * length / Fn2;
			f1 = f(x1);
		}
		k++;
		fout << (bPrev - aPrev) / (b - a) << "\t" << x1 << "\t" << x2 << "\t" << f1 << "\t" << f2 << endl;
	}

	fout << "x = " << x1 << "\tx2 - x1 = " << x2 - x1 << endl;
	fout.close();
	return k + 1;
}


/* Поиск интервала, содержащего минимум
 * a,b - границы отрезка
 * x0 - начальная точка
 * delta - точность
 * fileName - файл, в который будет сохранён результат
 */
void findRange(real a, real b, real x0, real delta, const char *fileName) {

	ofstream fout(fileName);
	real x = x0, xNext = x0, xPrev = x0, h, fx, fxNext;
	int k = 0;
	fout << setprecision(7) << fixed;


	// Шаг 1
	if (f(x0) > f(x0 + delta)) {
		x = x0;
		xNext = x + delta;
		h = delta;
	}
	else if (f(x0) > f(x0 - delta)) {
		x = x0;
		xNext = x - delta;
		h = -delta;
	}
	else {
		fout << ++k << "\t" << x0 << "\t" << f(x0)
			<< "\t[ " << x0 - delta << " ; " << x0 + delta << " ]" << endl;
		return;
	}

	fx = f(x);
	fxNext = f(xNext);
	do {
		// Сохраняем текущий шаг, как предыдущий
		xPrev = x;
		fx = fxNext;
		x = xNext;

		// Выполняем шаг 2
		h *= 2;
		xNext = x + h;
		fxNext = f(xNext);

		fout << ++k << "\t" << x << "\t" << fx
			<< "\t[ " << xPrev << " ; " << xNext << " ]" << endl;
	} while (fx > fxNext);

	fout.close();
}



/* Построение таблицы для графика зависисмости количества вычислений
 * целевой функции от логарифма задаваемой точности E
 * a,b - границы отрезка
 * E0,E1 - границы точности E
 * step - во сколько раз повышается точность на каждой итерации
 * fileName - файл, в который будет сохранён результат
 */
void createTable(real a, real b, real E0, real E1, real step) {

	ofstream fout("Results/Table.txt");
	fout << "ln(E)\tD\tGS\tF" << endl;
	for (real E = E0; E >= E1; E /= step) {
		fout << log(E) << "\t"
			<< calcByDichotomy(a, b, E, "Results/Dichotomy.txt") << "\t"
			<< calcByGoldenSectionSearch(a, b, E, "Results/GoldenSection.txt") << "\t"
			<< calcByFibonacci(a, b, E, "Results/Fibonacci.txt") << endl;
	}
	fout.close();
}



void main() {

	real a = -2, b = 20, E = 1e-2, x0 = -20;
	calcByDichotomy(a, b, E, "Results/Dichotomy.txt");
	calcByGoldenSectionSearch(a, b, E, "Results/GoldenSection.txt");
	calcByFibonacci(a, b, E, "Results/Fibonacci.txt");

	findRange(a, b, x0, 1e-1, "Results/Range.txt");
	createTable(a, b, 1e-1, 1e-7, 2);
}
