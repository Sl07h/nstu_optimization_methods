#include "head.h"
int fCalcCount, gCalcCount, alpha, rMult;
real r;

//------------------------------------------------------------------------------
struct methodResult
{
	real E;
	int iterationsCount;
	int fCalcCount;
	int r0;
	vector1D x0, x, xPrev, S, gradf;
	matrix2D A;
	real fx, fxPrev, lambda;
	void printResultE(std::ofstream &fout);
	void printResultR(std::ofstream &fout);
	void printResultRFirst(std::ofstream &fout);
};


// Вывод результатов в поток
void methodResult::printResultE(std::ofstream &fout) {
	fout << E << "\t"
		<< iterationsCount << "\t"
		<< fCalcCount << "\t"
		<< r << "\t"
		<< x << "\t"
		<< fx << endl;
}



// Вывод результатов в поток
void methodResult::printResultR(std::ofstream &fout) {
	fout << rMult << "\t"
		<< iterationsCount << "\t"
		<< fCalcCount << "\t"
		<< r << "\t"
		<< x << "\t"
		<< fx << endl;
}


// Вывод результатов в поток
void methodResult::printResultRFirst(std::ofstream &fout) {
	fout << r0 << "\t"
		<< iterationsCount << "\t"
		<< fCalcCount << "\t"
		<< r << "\t"
		<< x << "\t"
		<< fx << endl;
}


//------------------------------------------------------------------------------
// Функция для исследования min=0 в точке (1,1)
// x - вектор аргументов функции
inline real f(const vector1D &x) {
	fCalcCount++;
	return 4 * pow((x[1] - x[0]), 2) + 3 * pow((x[0] - 1), 2);
}


// Ограничение области
inline real g(const vector1D &x) {
	gCalcCount++;
	return x[0] + x[1] + 1;
}


//------------------------------------------------------------------------------
// Штрафная функция G
inline real G1(const vector1D &x) {
	if (g(x) > 0)
		return 0.5*(g(x) + fabs(g(x)));
	else
		return 0;
}
inline real G2(const vector1D &x) {
	if (g(x) > 0)
		return pow(0.5*(g(x) + fabs(g(x))), 2);
	else
		return 0;
}
inline real G3(const vector1D &x) {
	if (g(x) > 0)
		return pow(0.5*(g(x) + fabs(g(x))), alpha);
	else
		return 0;
}
inline real G4(const vector1D &x) {
	if (g(x) <= 0)
		return -log(-g(x));
	else
		return DBL_MAX;
}
inline real G5(const vector1D &x) {
	if (g(x) <= 0)
		return -1.0 / g(x);
	else
		return DBL_MAX;
}



//------------------------------------------------------------------------------
// Минимизируемая функция Q(x,y)
real Q1(const vector1D &x) { return f(x) + r * G1(x); }
real Q2(const vector1D &x) { return f(x) + r * G2(x); }
real Q3(const vector1D &x) { return f(x) + r * G3(x); }
real Q4(const vector1D &x) { return f(x) + r * G4(x); }
real Q5(const vector1D &x) { return f(x) + r * G5(x); }



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



// Определение коэффициента лямбда методом Фибоначчи
// f - минимизируемая функция
// x - начальное значение
// S - базис
// E - точность
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
methodResult calcByRosenbrock(const function<real(const vector1D &x)> &f, const function<real(const vector1D &x)> &G, const vector1D &x0, real r0, real E, const string &funcname, bool isFineNotBarier) {

	cout << endl << "Rosenbrock " << funcname << endl;
	ofstream fout("report/tableRosenbrock_" + funcname + ".txt");
	ofstream steps("steps/Rosenbrock_" + funcname + ".txt");
	fout << scientific;
	steps << fixed << setprecision(12);
	steps << x0 << endl;

	methodResult result;
	fCalcCount = 0;
	gCalcCount = 0;
	vector1D xPrev, B, x = x0;
	int maxiter = 1000;
	int iterationsCount = 0, count = 0;
	r = r0;
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
			<< fCalcCount << "\t"
			<< r << "\t"
			<< x << "\t"
			<< f(x) << endl;

		steps << x << endl;

		if (isFineNotBarier)
			r *= rMult;
		else
			r /= rMult;

	} while (abs(f(x) - f(xPrev)) > E && abs(calcNormE(x) - calcNormE(xPrev)) > E && iterationsCount < maxiter && r*G(x) > E);

	fout.close();
	steps.close();

	result.E = E;
	result.iterationsCount = iterationsCount;
	result.fCalcCount = fCalcCount;
	result.x0 = x0;
	result.x = x;
	result.fx = f(x);
	result.r0 = r0;

	return result;
}



//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Таблицы
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


// Одиночный тест
void singleTest()
{
	real E = 1e-12;
	real r0 = 10;
	rMult = 2;
	alpha = 4;

	vector1D x0Fine = { -2, 0 };
	vector1D x0Barrier = { 0, -2 };
	
	calcByRosenbrock(Q1, G1, x0Fine, r0, E, "Q1", true);
	calcByRosenbrock(Q2, G2, x0Fine, r0, E, "Q2", true);
	calcByRosenbrock(Q3, G3, x0Fine, r0, E, "Q3", true);
	calcByRosenbrock(Q4, G4, x0Barrier, r0, E, "Q4", false);
	calcByRosenbrock(Q5, G5, x0Barrier, r0, E, "Q5", false);

	// Визуализация
	string runVisualisation = "python plot.py "
		+ to_string(x0Fine[0]) + " " + to_string(x0Fine[1]) + " "
		+ to_string(x0Barrier[0]) + " " + to_string(x0Barrier[1]);
	system(runVisualisation.c_str());
}


// Зависимость скорости сходимости метода от стратегии изменения коэффициента штрафа
void researchRMult()
{
	ofstream foutF("steps/tableRMultFines.txt");
	ofstream foutB("steps/tableRMultBarriers.txt");
	ofstream fout1("report/tableRMult1.txt");
	ofstream fout2("report/tableRMult2.txt");
	ofstream fout3("report/tableRMult3.txt");
	ofstream fout4("report/tableRMult4.txt");
	ofstream fout5("report/tableRMult5.txt");
	fout1 << scientific;
	fout2 << scientific;
	fout3 << scientific;
	fout4 << scientific;
	fout5 << scientific;

	methodResult result;
	real E = 1e-12;
	real r0 = 10;
	alpha = 8;
	vector1D x0Fine = { -2, 0 };
	vector1D x0Barrier = { 0, -2 };

	for (size_t i = 2; i <= 10; i++)
	{
		rMult = i;
		result = calcByRosenbrock(Q1, G1, x0Fine, r0, E, "Q1", true);
		result.printResultR(fout1);
		foutF << result.x << endl;
		result = calcByRosenbrock(Q2, G2, x0Fine, r0, E, "Q2", true);
		result.printResultR(fout2);
		foutF << result.x << endl;
		result = calcByRosenbrock(Q3, G3, x0Fine, r0, E, "Q3", true);
		result.printResultR(fout3);
		foutF << result.x << endl;

		result = calcByRosenbrock(Q4, G4, x0Barrier, r0, E, "Q4", false);
		result.printResultR(fout4);
		foutB << result.x << endl;
		result = calcByRosenbrock(Q5, G5, x0Barrier, r0, E, "Q5", false);
		result.printResultR(fout5);
		foutB << result.x << endl;
	}

	// Визуализация
	string runVisualisation = "python plot.py "
		+ to_string(x0Fine[0]) + " " + to_string(x0Fine[1]) + " "
		+ to_string(x0Barrier[0]) + " " + to_string(x0Barrier[1]);
	system(runVisualisation.c_str());

	foutF.close();
	foutB.close();
	fout1.close();
	fout2.close();
	fout3.close();
	fout4.close();
	fout5.close();
}



// Зависимость скорости сходимости метода от стратегии изменения коэффициента штрафа
void researchRFirst()
{
	ofstream foutF("steps/tableRFirstFines.txt");
	ofstream foutB("steps/tableRFirstBarriers.txt"); 
	ofstream fout1("report/tableRFirst1.txt");
	ofstream fout2("report/tableRFirst2.txt");
	ofstream fout3("report/tableRFirst3.txt");
	ofstream fout4("report/tableRFirst4.txt");
	ofstream fout5("report/tableRFirst5.txt");
	fout1 << scientific;
	fout2 << scientific;
	fout3 << scientific;
	fout4 << scientific;
	fout5 << scientific;

	methodResult result;
	real E = 1e-12;
	rMult = 2;
	alpha = 8;
	vector1D x0Fine = { -2, 0 };
	vector1D x0Barrier = { 0, -2 };

	for (real r0 = 2; r0 <= 64; r0 *= 2)
	{
		result = calcByRosenbrock(Q1, G1, x0Fine, r0, E, "Q1", true);
		result.printResultRFirst(fout1);
		foutF << result.x << endl;
		result = calcByRosenbrock(Q2, G2, x0Fine, r0, E, "Q2", true);
		result.printResultRFirst(fout2);
		foutF << result.x << endl;
		result = calcByRosenbrock(Q3, G3, x0Fine, r0, E, "Q3", true);
		result.printResultRFirst(fout3);
		foutF << result.x << endl;

		result = calcByRosenbrock(Q4, G4, x0Barrier, r0, E, "Q4", false);
		result.printResultRFirst(fout4);
		foutB << result.x << endl;
		result = calcByRosenbrock(Q5, G5, x0Barrier, r0, E, "Q5", false);
		result.printResultRFirst(fout5);
		foutB << result.x << endl;
	}

	// Визуализация
	string runVisualisation = "python plot.py "
		+ to_string(x0Fine[0]) + " " + to_string(x0Fine[1]) + " "
		+ to_string(x0Barrier[0]) + " " + to_string(x0Barrier[1]);
	system(runVisualisation.c_str());

	foutF.close();
	foutB.close();
	fout1.close();
	fout2.close();
	fout3.close();
	fout4.close();
	fout5.close();
}



// Зависимость скорости сходимости метода от заданной точности
void researchE()
{
	ofstream foutF("steps/tableEFines.txt");
	ofstream foutB("steps/tableEBarriers.txt");
	ofstream fout1("report/tableE1.txt");
	ofstream fout2("report/tableE2.txt");
	ofstream fout3("report/tableE3.txt");
	ofstream fout4("report/tableE4.txt");
	ofstream fout5("report/tableE5.txt");
	fout1 << scientific;
	fout2 << scientific;
	fout3 << scientific;
	fout4 << scientific;
	fout5 << scientific;

	methodResult result;
	real r0 = 10;
	rMult = 2;
	alpha = 8;
	vector1D x0Fine = { -2, 0 };
	vector1D x0Barrier = { 0, -2 };

	for (double E = 1e-3; E >= 1e-7; E /= 10)
	{
		result = calcByRosenbrock(Q1, G1, x0Fine, r0, E, "Q1", true);
		result.printResultE(fout1);
		foutF << result.x << endl;
		result = calcByRosenbrock(Q2, G2, x0Fine, r0, E, "Q2", true);
		result.printResultE(fout2);
		foutF << result.x << endl;
		result = calcByRosenbrock(Q3, G3, x0Fine, r0, E, "Q3", true);
		result.printResultE(fout3);
		foutF << result.x << endl;

		result = calcByRosenbrock(Q4, G4, x0Barrier, r0, E, "Q4", false);
		result.printResultE(fout4);
		foutB << result.x << endl;
		result = calcByRosenbrock(Q5, G5, x0Barrier, r0, E, "Q5", false);
		result.printResultE(fout5);
		foutB << result.x << endl;
	}

	// Визуализация
	string runVisualisation = "python plot.py "
		+ to_string(x0Fine[0]) + " " + to_string(x0Fine[1]) + " "
		+ to_string(x0Barrier[0]) + " " + to_string(x0Barrier[1]);
	system(runVisualisation.c_str());

	foutF.close();
	foutB.close();
	fout1.close();
	fout2.close();
	fout3.close();
	fout4.close();
	fout5.close();
}


void main()
{
	researchE();
	researchRFirst();
	researchRMult();
	singleTest();


}