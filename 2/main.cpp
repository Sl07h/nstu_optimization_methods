#include "head.h"
int fCalcCount;

struct methodResult
{
	real E;
	int iterationsCount;
	int fCalcCount;
	vector1D x0, x, xPrev, S, gradf;
	matrix2D A;
	real fx, fxPrev, lambda;
	void printFirstResult(std::ofstream &fout);
};


// Вывод результатов в поток для 1 таблицы
void methodResult::printFirstResult(std::ofstream &fout) {
	fout << E << "\t"
		<< iterationsCount << "\t"
		<< fCalcCount << "\t"
		<< x0 << "\t"
		<< x << "\t"
		<< fx << endl;
}


// Заданная целевая функция
// x - вектор аргументов функции
real f1(const vector1D &x) {
	// - (2 / (1 + (((x - 1) / 2)^ 2) + (((y - 1) / 1)^ 2)) +3 / (1 + (((x - 2) / 3) ^ 2) + (((y - 3) / 2) ^ 2)))

	///minimize -(2/(1 + ((x - 1)/2)^2 + ((y - 1) /1)^2) + 3/(1 + ((x - 2)/3)^2 + ((y - 3)/2)^2))
	// 1 вариант -3.5
	int A1 = 2, A2 = 3, a1 = 1, a2 = 2, b1 = 2, b2 = 3, c1 = 1, c2 = 3, d1 = 1, d2 = 2;


	// 2 вариант
	//int A1 = 1, A2 = 3, a1 = 2, a2 = 1, b1 = 3, b2 = 1, c1 = 2, c2 = 1, d1 = 3, d2 = 2;

	fCalcCount++;
	return  -(A1 / (1 + pow((x[0] - a1) / b1, 2) + pow((x[1] - c1) / d1, 2))
		+ A2 / (1 + pow((x[0] - a2) / b2, 2) + pow((x[1] - c2) / d2, 2)));
}



// Функция для исследования min=0 в точке (1,1)
// x - вектор аргументов функции
real f2(const vector1D &x) {
	fCalcCount++;
	return 100 * pow((x[1] - x[0]), 2) + pow((1 - x[0]), 2);
}


// Функция Розенброка min=0 в точке (1,1)
// x - вектор аргументов функции
real f3(const vector1D &x) {
	fCalcCount++;
	return 100 * pow((x[1] - x[0] * x[0]), 2) + pow((1 - x[0]), 2);
}





// Поиск интервала, содержащего минимум
// f - целевая одномерная функция
// a,b - искомые границы отрезка
// x - начальное приближение
// S - орты (направления)
void interval(const function<real(const vector1D &x)> &f, real &a, real &b, vector1D &x, vector1D &S)
{
	real lambda0 = 0.0;
	real delta = 1.0e-8;
	real lambda_k_minus_1 = lambda0;
	real f_k_minus_1 = f(x + S * lambda_k_minus_1);
	real lambda_k;
	real f_k;
	real lambda_k_plus_1;
	real f_k_plus_1;
	real h;
	if (f(x + S * lambda0) > f(x + S * (lambda0 + delta)))
	{
		lambda_k = lambda0 + delta;
		h = delta;
	}
	else
	{
		lambda_k = lambda0 - delta;
		h = -delta;
	}
	f_k = f(x + S * lambda_k);
	while (true)
	{
		h *= 2.0;
		lambda_k_plus_1 = lambda_k + h;
		f_k_plus_1 = f(x + S * lambda_k_plus_1);
		if (f_k > f_k_plus_1)
		{
			lambda_k_minus_1 = lambda_k;
			f_k_minus_1 = f_k;
			lambda_k = lambda_k_plus_1;
			f_k = f_k_plus_1;
		}
		else
		{
			a = lambda_k_minus_1;
			b = lambda_k_plus_1;
			if (b < a)
				swap(a, b);
			return;
		}
	}
}


// Вычисление n-го числа Фибоначии
inline real fib(int n)
{
	real sqrt5 = sqrt(5.0), pow2n = pow(2.0, n);
	return (pow(1.0 + sqrt5, n) / pow2n - pow(1.0 - sqrt5, n) / pow2n) / sqrt5;
}


/* Определение коэффициента лямбда методом Фибоначчи
 * f - минимизируемая функция
 * x - начальное значение
 * S - базис
 * E - точность
 */
real fibonacci(const function<real(const vector1D &x)> &f, vector1D &x, vector1D &S, real E)
{
	real a, b;
	interval(f, a, b, x, S);
	int iter;
	real len = fabs(a - b);
	int n = 0;
	while (fib(n) < (b - a) / E) n++;
	iter = n - 3;
	real lambda1 = a + (fib(n - 2) / fib(n)) * (b - a);
	real f1 = f(x + S * lambda1);
	real lambda2 = a + (fib(n - 1) / fib(n)) * (b - a);
	real f2 = f(x + S * lambda2);
	for (int k = 0; k < n - 3; k++)
	{
		if (f1 <= f2)
		{
			b = lambda2;
			lambda2 = lambda1;
			f2 = f1;
			lambda1 = a + (fib(n - k - 3) / fib(n - k - 1)) * (b - a);
			f1 = f(x + S * lambda1);
		}
		else
		{
			a = lambda1;
			lambda1 = lambda2;
			f1 = f2;
			lambda2 = a + (fib(n - k - 2) / fib(n - k - 1)) * (b - a);
			f2 = f(x + S * lambda2);
		}
		len = b - a;
	}
	lambda2 = lambda1 + E;
	f2 = f(x + S * lambda2);
	if (f1 <= f2)
		b = lambda1;
	else
		a = lambda1;
	return (a + b) / 2.0;
}


// Метод Розенброка
// f - оптимизиуремая функция
// x0 - начальное приближение
// E - точность
// funcname - название функции
methodResult calcByRosenbrock(const function<real(const vector1D &x)> &f, const vector1D &x0, real E, const string &funcname) {

	cout << endl << "Rosenbrock " << funcname << endl;
	ofstream fout("report/SecondTableRosenbrock_" + funcname + ".txt");
	ofstream steps("steps/Rosenbrock_" + funcname + ".txt");
	fout << fixed << setprecision(4);
	steps << fixed << setprecision(4);
	steps << x0 << endl;

	methodResult result;
	fCalcCount = 0;
	vector1D xPrev, B, x = x0;
	int maxiter = 100;
	int iterationsCount = 0, count = 0;
	real lambda1, lambda2;
	matrix2D A(2);
	A[0] = { 1.0, 0.0 };
	A[1] = { 0.0, 1.0 };


	// Начальные ортогональные направления
	matrix2D S(2);
	S[0] = { 1.0, 0.0 };
	S[1] = { 0.0, 1.0 };

	do {
		xPrev = x;

		// Минимализируем функцию в направлениях S^k_1...S^k_n
		lambda1 = fibonacci(f, x, S[0], E);
		x = x + S[0] * lambda1;
		lambda2 = fibonacci(f, x, S[1], E);
		x = x + S[1] * lambda2;

		// Построение новых ортогональных направлений при
		// сортировке лямбд в порядке убывания по абсолютным значениям
		A[0] = S[0] * lambda1 + S[1] * lambda2;
		if (fabs(lambda1 >= lambda2))
			A[1] = S[1] * lambda2;
		else
			A[1] = S[0] * lambda1;


		// Ортогонализация Грамма-Шмидта
		S[0] = A[0] / calcNormE(A[0]);
		B = A[1] - S[1] * A[1] * S[1];
		if (calcNormE(B) > E)
			S[1] = B / calcNormE(B);
		iterationsCount++;
		fout << iterationsCount << "\t"
			<< x << "\t"
			<< f(x) << "\t"
			<< S << "\t"
			<< lambda1 << "\t"
			<< x - xPrev << "\t"
			<< f(x) - f(xPrev) << "\t"
			<< A << endl;

		steps << x << endl;
	} while (abs(f(x) - f(xPrev)) > E && abs(calcNormE(x) - calcNormE(xPrev)) > E && iterationsCount < maxiter);

	fout.close();
	steps.close();

	result.E = E;
	result.iterationsCount = iterationsCount;
	result.fCalcCount = fCalcCount;
	result.x0 = x0;
	result.x = x;
	result.fx = f(x);

	return result;
}



//-----------------------------------------------------------------------------

// Градиент
vector1D grad(const function<real(const vector1D &x)> &f, const vector1D &x1) {
	const real eps = 1e-7;
	const real fx1 = f(x1);

	vector1D result(x1.size());
	vector1D x = x1;

	for (int i = 0; i < x.size(); ++i) {
		x[i] += eps;
		result[i] = (f(x) - fx1) / eps;
		x[i] -= eps;
	}
	return result;
}


// Умножение вектора a' на вектор b и получение матрицы
matrix2D multTwoVectorsToMatrix(const vector1D &a, const vector1D &b) {

	matrix2D tmp;
	tmp.resize(a.size());
	for (size_t i = 0; i < tmp.size(); i++)
	{
		tmp[i] = a[i] * b;
	}

	return tmp;
}


// Метод Бройдена
// f - оптимизиуремая функция
// x0 - начальное приближение
// E - точность
// funcname - название функции
methodResult calcByBroyden(const function<real(const vector1D &x)> &f, const vector1D &x0, real E, const string &funcname) {

	cout << endl << "Broyden " << funcname << endl;
	ofstream fout("report/SecondTableBroyden_" + funcname + ".txt");
	ofstream steps("steps/Broyden_" + funcname + ".txt");
	fout << fixed << setprecision(4);
	steps << fixed << setprecision(4);
	steps << x0 << endl;

	methodResult result;
	fCalcCount = 0;
	matrix2D A(2);
	A[0] = { 1.0, 0.0 };
	A[1] = { 0.0, 1.0 };

	vector1D x = x0, xPrev = x;
	vector1D gradf = grad(f, x);
	int iterationsCount = 0, maxiter = 100;

	do
	{
		if (iterationsCount % 2 == 0) {
			A[0] = { 1.0, 0.0 };
			A[1] = { 0.0, 1.0 };
			cout << "Method renewal" << endl;
		}
		
		if (det(A) <= 0 || !isfinite(det(A))) {
			A[0] = { 1.0, 0.0 };
			A[1] = { 0.0, 1.0 };
			cout << "Matrix isn't positive-defined" << endl;
		}

		

		if (calcNormE(gradf) <= E || iterationsCount >= maxiter) {

			result.E = E;
			result.iterationsCount = iterationsCount;
			result.fCalcCount = fCalcCount;
			result.x0 = x0;
			result.x = x;
			result.fx = f(x);

			return result;
		}

		vector1D ngradf = A * gradf;
		// получаем лямбду
		real lambda = fibonacci(f, x, ngradf, E);
		//cout << " lambda = "<< lambda << "\t";

		vector1D dx = lambda * ngradf;
		xPrev = x;
		x = x + dx;

		vector1D old_gradf = gradf;
		gradf = grad(f, x);

		vector1D dg = gradf - old_gradf;
		vector1D temp = dx - A * dg;

		// Находится очередное приближение матрицы H^(-1)
		matrix2D dA = multTwoVectorsToMatrix(temp, temp) / (temp * dg);

		if (temp * dg == 0) {
			A[0] = { 1.0, 0.0 };
			A[1] = { 0.0, 1.0 };
			cout << "Zero division" << endl;
		}
		else
			A = A + dA;

		iterationsCount++;

		fout << iterationsCount << "\t"
			<< x << "\t"
			<< f(x) << "\t"
			<< lambda << "\t"
			<< x - xPrev << "\t"
			<< f(x) - f(xPrev) << "\t"
			<< gradf << "\t"
			<< A << endl;

		steps << x << endl;
	} while (abs(f(x) - f(xPrev)) > E && abs(calcNormE(x) - calcNormE(xPrev)) > E && iterationsCount < maxiter);

	fout.close();
	steps.close();

	result.E = E;
	result.iterationsCount = iterationsCount;
	result.fCalcCount = fCalcCount;
	result.x0 = x0;
	result.x = x;
	result.fx = f(x);

	return result;
}


// Генерация первых таблиц
void makeFirstTables(const function<real(const vector1D &x)> &f, const vector1D &x0, const string &funcname) {

	ofstream foutR("report/FirstTableRosenbrock_" + funcname + ".txt");
	ofstream foutB("report/FirstTableBroyden_" + funcname + ".txt");
	foutR << fixed << setprecision(7);
	foutB << fixed << setprecision(7);

	methodResult result;

	for (double E = 1e-3; E >= 1e-7; E /= 10)
	{
		result = calcByRosenbrock(f, x0, E, funcname);
		result.printFirstResult(foutR);

		result = calcByBroyden(f, x0, E, funcname);
		result.printFirstResult(foutB);
	}


	foutR.close();
	foutB.close();
}


// Генерация первых таблиц
void exploreConvergence(const function<real(const vector1D &x)> &f, const vector1D &x0, const string &funcname) {

	ofstream foutR("report/FirstTableRosenbrock_" + funcname + ".txt");
	ofstream foutB("report/FirstTableBroyden_" + funcname + ".txt");
	foutR << fixed << setprecision(7);
	foutB << fixed << setprecision(7);

	methodResult result;

	for (double E = 1e-3; E >= 1e-7; E /= 10)
	{
		result = calcByRosenbrock(f, x0, E, funcname);
		result.printFirstResult(foutR);

		result = calcByBroyden(f, x0, E, funcname);
		result.printFirstResult(foutB);
	}


	foutR.close();
	foutB.close();
}


void main() {

	vector1D x0 = { 1.5, -0.5 };
	real E = 1e-7;

	
	/*makeFirstTables(f1, x0, "f1");
	makeFirstTables(f2, x0, "f2");
	makeFirstTables(f3, x0, "f3");*/
	

	calcByRosenbrock(f1, x0, E, "f1");
	calcByRosenbrock(f2, x0, E, "f2");
	calcByRosenbrock(f3, x0, E, "f3");

	calcByBroyden(f1, x0, E, "f1");
	calcByBroyden(f2, x0, E, "f2");
	calcByBroyden(f3, x0, E, "f3");

	string runVisualisation = "python plot.py " + to_string(x0[0]) + " " + to_string(x0[1]);
	system(runVisualisation.c_str());
}