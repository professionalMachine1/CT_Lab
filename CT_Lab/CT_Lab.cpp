#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include <random>
#include <armadillo>
#include "ODE.h"
#include "ORE.h"

using namespace std;
using namespace std::complex_literals;

struct TaskParametrs
{
	double m, M, I, l, Bp, Beq, g, Kf, Ks, tmax;
};

TaskParametrs P;

const string setting = "set terminal pngcairo size 800,600 enhanced font 'Times New Roman,12'\n";
const string way = "C:\\Users\\HOME\\source\\repos\\CT_Lab\\CT_Lab\\putHere\\";
void PLOT(const string& commands)
{
	FILE* gpipe = _popen("C:\\Users\\gnuplot\\bin\\gnuplot.exe -persist", "w");
	fprintf(gpipe, commands.c_str());
	_pclose(gpipe);
}

/* Инициализация и решение */

void ContinuosSystem();
void DiscreteSystem();
void LinearQuadraticControl();
void KalmanFilter();
Matrix<double> MatrixExp(Matrix<double> A);
Matrix<double> IntegrateMatrixExp(Matrix<double>& A, double h);
Vector<double> GetControl(Matrix<double>& A, const Vector<double>& B, Vector<complex<double>>& eigV); 
Matrix<double> RikkatiSolver(const Matrix<double>& A_, const Matrix<double>& B_, const Matrix<double>& Q_, const Matrix<double>& R_);
tuple<Matrix<double>, Vector<double>> InicialSyst();
Vector<double>& NonLinSystCoef(Vector<double>& p);

/* Картинки */

void GetImageDiscrSystem(DiscreteSystemParametrs& DS, Vector<double>& x_start);
void GetImageDiscrSystemObs(DiscreteSystemParametrs& DS, Vector<double>& x_start, Vector<double>& ksi_start);
void GetImageOriginalSystem(ContinuosSystemParametrs& SP, Vector<double>& pNon, Vector<double>& x_start);
void GetImageLinObserver(ContinuosSystemParametrs& SP, Vector<double>& x_start);
void GetImageNonLinObserver(ContinuosSystemParametrs& SP, Vector<double>& pNon, Vector<double>& x_start);

/* Налаживание совместимости */

Vector<double> operator+(const Vector<double>& v1, const arma::vec& v2);
Vector<double> operator+(const arma::vec& v1, const Vector<double>& v2);
Matrix<double> mat_to_Matrix(const arma::mat& mt);
arma::mat Matrix_to_mat(const Matrix<double>& mt);

#define PRINT(name) cout << #name << " = " << name << endl
#define PRINTF(name) cout << #name << " = \n" << name << endl
#define GETNAME(name) cout << #name << endl

int main()
{
	system("chcp 1251");
	system("cls");

	srand(time(0));

	P.tmax = 16;
	P.m = 0.127, P.M = 1.206, P.I = 1.2 * pow(10, -3), P.l = 0.1778, P.Kf = 1.726, P.Ks = 4.487;
	P.Beq = 5.4, P.Bp = 2.4 * pow(10, -3), P.g = 9.81;

	/* Решение для непрерывной системы */
	//GETNAME(Непрерывная система);
	//ContinuosSystem();

	/* Решение для дискретной системы */
	//GETNAME(Дискретная система);
	//DiscreteSystem();

	/* Решение для непрерывной системы */
	//GETNAME(Линейно-квадратичное управление);
	//LinearQuadraticControl();

	/* Решение для дискретной системы */
	GETNAME(Фильтр Калмана);
	KalmanFilter();



	return 0;
}

void ContinuosSystem()
{
	//Инициализируем необходимые структуры
	ContinuosSystemParametrs CSP;

	int n = 4;

	Matrix<double> AT(n);
	CSP.A = Matrix<double>(n), CSP.L = Vector<double>(n);
	CSP.B = Vector<double>(n), CSP.C = Vector<double>(n);
	Vector<double> pNon, x_start(n);
	Vector<complex<double>> eigV(n);
	CSP.C[0] = 1, CSP.C[1] = 1;

	/*

		Без наблюдателя

	*/

	//Инициализируем A, B
	tie(CSP.A, CSP.B) = InicialSyst();

	//Получаем theta с заданными собственными числами и выводим его на экран
	eigV[0] = -11.1856;
	eigV[1] = -6.33331;
	eigV[2] = -1.0;
	eigV[3] = -1.0;
	CSP.theta = GetControl(CSP.A, CSP.B, eigV);
	PRINT(CSP.theta);

	//Заполняем структуру, которая потом будет применятся в решении ДУ
	CSP.BTheta = CellRowMultiply(CSP.B, CSP.theta), CSP.Ac = CSP.A + CSP.BTheta;
	//Заполняем коэффициенты для нелинейной системы
	NonLinSystCoef(pNon);

	//Начальные условия
	x_start[0] = 0.5, x_start[1] = 0.5, x_start[2] = 0.3, x_start[3] = 0.5;

	//Строим графики исходной системы и сохраняем их в формате картинок
	GetImageOriginalSystem(CSP, pNon, x_start);

	/*

		Наблюдатель

	*/

	//Ищем L = -theta^T для наблюдателя и выводим его
	eigV[0] = -11.1856;
	eigV[1] = -6.33331;
	eigV[2] = -10.0;
	eigV[3] = 0;
	AT = TransposedMatrix(CSP.A);
	CSP.L = GetControl(AT, CSP.C, eigV);
	CSP.L = -1.0 * CSP.L;
	PRINT(CSP.L);

	//Заполняем ещё один элемент структуры
	CSP.LC = CellRowMultiply(CSP.L, CSP.C);

	x_start = Vector<double>(2 * n);
	//Начальные условия для наблюдателя
	x_start[0] = 0.3, x_start[1] = 0.5, x_start[2] = 0.3, x_start[3] = 0.5;
	//Начальные условия для исходной системы
	x_start[4] = 0.3, x_start[5] = 0.5, x_start[6] = 0.1, x_start[7] = 0.2;

	//Строим графики линейной и нелинейной систем с наблюдателем и сохраняем их в формате картинок
	GetImageLinObserver(CSP, x_start);
	GetImageNonLinObserver(CSP, pNon, x_start);
}

void DiscreteSystem()
{
	DiscreteSystemParametrs DS;
	int n = 4;
	Matrix<double> A(n), AT(n), I = Single<double>(n);
	Vector<double> B(n), C(n), x_start(n), ksi_start(n);
	Vector<complex<double>> eigV(n);
	C[0] = 1, C[1] = 1;
	DS.h = pow(10, -1);

	/*

		Без наблюдателя

	*/

	//Инициализируем Ad, Bd
	tie(A, B) = InicialSyst();
	DS.Ad = MatrixExp(A * DS.h);

	//Считаем Bd в зависимости от того, каким является управление - импульсным или кусочно-постоянным
	DS.Bd_impulse = DS.Ad * B;
	DS.Bd_const = IntegrateMatrixExp(A, DS.h) * B;

	//Получаем theta с заданными собственными числами
	eigV[0] = exp(-11.1856 * DS.h);
	eigV[1] = exp(-6.33331 * DS.h);
	eigV[2] = exp(-1.0 * DS.h);
	eigV[3] = exp(-1.0 * DS.h);
	DS.thetad_const = GetControl(DS.Ad, DS.Bd_const, eigV);
	PRINT(DS.thetad_const);
	DS.thetad_impulse = GetControl(DS.Ad, DS.Bd_impulse, eigV);
	PRINT(DS.thetad_impulse);

	DS.BTd_const = CellRowMultiply(DS.Bd_const, DS.thetad_const);
	DS.ABTd_const = DS.Ad + DS.BTd_const;

	DS.BTd_impulse = CellRowMultiply(DS.Bd_impulse, DS.thetad_impulse);
	DS.ABTd_impulse = DS.Ad + DS.BTd_impulse;

	//Начальные условия
	x_start[0] = 0.5, x_start[1] = 0.5, x_start[2] = 0.4, x_start[3] = -0.2;
	GetImageDiscrSystem(DS, x_start);

	/*

		Наблюдатель

	*/

	//Получаем L = -theta^T с заданными собственными числами
	eigV[0] = exp(-11.1856 * DS.h);
	eigV[1] = exp(-6.33331 * DS.h);
	eigV[2] = exp(-10.0 * DS.h);
	eigV[3] = 1;
	AT = TransposedMatrix(DS.Ad);
	DS.Ld = GetControl(AT, C, eigV);
	DS.Ld = -1.0 * DS.Ld;
	PRINT(DS.Ld);

	DS.LCd = CellRowMultiply(DS.Ld, C);

	//Начальные условия
	ksi_start[0] = 0.6, ksi_start[1] = 0.5, ksi_start[2] = -0.1, ksi_start[3] = 0.1;
	x_start[0] = 0.6, x_start[1] = 0.5, x_start[2] = -0.4, x_start[3] = 0.5;

	GetImageDiscrSystemObs(DS, x_start, ksi_start);

}

void LinearQuadraticControl()
{
	Matrix<double> A(4), R(1), Q(4), X(4), B(4, 1);
	Q.diag({ 5, 5, 5, 5 });
	R.diag({ 1.5 });

	tie(A, B(0)) = InicialSyst();
	X = RikkatiSolver(A, B, Q, R);

	ContinuosSystemParametrs CSP;
	CSP.RBX = -1.0 * R.i() * B.t() * X;
	CSP.Ac = A + B * CSP.RBX;

	/*

		Решаем, пишем данные в файл и делаем картинки

	*/

	Vector<double> x_start(4), p;
	x_start = { 0.57, 0.21, 0.42, -0.35 };

	double eps = pow(10, -10), h = pow(10, -3);
	int MaxIt = 20000;

	//Заполняем файлы линейной системы
	ODE4 solve(0, x_start, h, eps, 1);
	solve.SetStructParam(CSP);
	solve.WriteToFileLQR(P.tmax, "L_LQR", MaxIt);

	//Заполняем файлы нелинейной системы
	solve = ODE4(0, x_start, h, eps, 5);
	NonLinSystCoef(p);
	solve.SetVecParam(p);
	solve.SetStructParam(CSP);
	solve.WriteToFileLQR(P.tmax, "NL_LQR", MaxIt);

	string s;
	s = setting;
	s += "set grid\n";

	s += "set output 'LQR\\1R_fi.png'\n";
	s += "set multiplot layout 2,2\n";
	s += "plot '" + way + "L_LQR0.txt' with lines lc 'red' lw 2 title 'lin-{/Symbol f}', ";
	s += "'" + way + "NL_LQR0.txt' with lines lc 'blue' lw 2 title 'nonlin-{/Symbol f}'\n";

	//s += "set output 'LQR\\1R_x.png'\n";
	s += "plot '" + way + "L_LQR1.txt' with lines lc 'red' lw 2 title 'lin-x', ";
	s += "'" + way + "NL_LQR1.txt' with lines lc 'blue' lw 2 title 'nonlin-x'\n";

	//s += "set output 'LQR\\1R_dfi.png'\n";
	s += "plot '" + way + "L_LQR2.txt' with lines lc 'red' lw 2 title 'lin-d{/Symbol f}/dt', ";
	s += "'" + way + "NL_LQR2.txt' with lines lc 'blue' lw 2 title 'nonlin-d{/Symbol f}/dt'\n";

	//s += "set output 'LQR\\1R_dx.png'\n";
	s += "plot '" + way + "L_LQR3.txt' with lines lc 'red' lw 2 title 'lin-dx/dt', ";
	s += "'" + way + "NL_LQR3.txt' with lines lc 'blue' lw 2 title 'nonlin-dx/dt'\n";
	s += "unset multiplot\n";

	s += "set output 'LQR\\1R_u.png'\n";
	s += "plot '" + way + "L_LQRcontrol0.txt' with lines lc 'red' lw 2 title 'lin-u', ";
	s += "'" + way + "NL_LQRcontrol0.txt' with lines lc 'blue' lw 2 title 'nonlin-u'\n";

	PLOT(s);

}

void KalmanFilter()
{
	double h = pow(10, -1), t;
	int MaxIt;

	/* ---------------------- */

	Matrix<double> A, Ad, Bth, Ac;
	Vector<double> B, Bd, theta;
	Vector<complex<double>> eigv(4);
	tie(A, B) = InicialSyst();
	Ad = MatrixExp(A * h);
	//Кусочно-постоянное управление
	Bd = IntegrateMatrixExp(A, h) * B;
	//Импульсное управление
	//Bd = Ad * B;

	//Получаем theta с заданными собственными числами
	eigv[0] = exp(-11.1856 * h);
	eigv[1] = exp(-6.33331 * h);
	eigv[2] = exp((-1.0 + 1i) * h);
	eigv[3] = exp((-1.0 - 1i) * h);
	theta = GetControl(Ad, Bd, eigv);

	//Записали матрицу замкнутой системы
	Bth = CellRowMultiply(Bd, theta);
	Ac = Ad + Bth;

	/* ---------------------- */

	Matrix<double> W, E(4), C(2, 4), Adt, Ct;
	Vector<double> y(2), x(4), x_cor;
	C[0][0] = 1, C[1][1] = 1;

	W.diag({ 0.0004, 0.0004 });

	x = { 0.21, 0.42, -0.001, 0.2 };

	//Изначальные среднее и ковариационная матрица
	x_cor = { 0.2, 0.4, 0.0, 0.15 };
	E[0] = { 0.04, 0 ,0, 0 };
	E[1] = { 0, 0.04, 0, 0 };
	E[2] = { 0, 0, 0.04, 0 };
	E[3] = { 0, 0, 0, 0.04 };

	Adt = Ad.t(), Ct = C.t();
	
	Vector<ofstream> datax(4), datax_cor(4);
	for (int i = 0; i < 4; i++)
	{
		datax[i].open("putHere\\Kalman" + to_string(i) + ".txt", ios_base::out);
		datax_cor[i].open("putHere\\corKalman" + to_string(i) + ".txt", ios_base::out);
	}
	ofstream udata("putHere\\Kalmancontrol.txt", ios_base::out);

	t = 0;
	MaxIt = 100;

	//Support
	arma::mat aW = Matrix_to_mat(W);
	arma::vec w(2);

	for (int i = 0; i < MaxIt; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			datax[j] << t << " " << x[j] << endl;
			datax_cor[j] << t << " " << x_cor[j] << endl;
		}
		udata << t << " " << theta * x_cor << endl;

		//Делаем шаг
		x = Ad * x + Bth * x_cor;
		//x = Ac * x_cor;
		//Наблюдаемый выход
		y = C * x + arma::mvnrnd(w, aW);
		t += h;

		//Априорная оценка условной ковариационной матрицы
		E = Ad * E * Adt;
		//E = Ac * E * Ac.t();
		//Апостериорная оценка условной ковариацинной матрицы
		E = E - E * Ct * (W + C * E * Ct).i() * C * E;
		//Априорная оценка состояния объекта
		x_cor = Ac * x_cor;
		//Апостериорная оценка состояния объекта
		x_cor = x_cor + E * Ct * W.i() * (y - C * x_cor);

	}

	for (int i = 0; i < 4; i++)
	{
		datax[i].close();
		datax_cor[i].close();
	}
	udata.close();

	string s;
	s = setting;
	s += "set grid\n";

	s += "set output 'KalmanFilter\\1R_fi.png'\n";
	s += "set multiplot layout 2,2\n";
	s += "plot '" + way + "Kalman0.txt' with lines lc 'red' lw 2 title '{/Symbol f}', ";
	s += "'" + way + "corKalman0.txt' with lines lc 'blue' lw 2 title '{/Symbol f}_{cor}'\n";

	//s += "set output 'KalmanFilter\\1R_x.png'\n";
	s += "plot '" + way + "Kalman1.txt' with lines lc 'red' lw 2 title 'x', ";
	s += "'" + way + "corKalman1.txt' with lines lc 'blue' lw 2 title 'x_{cor}'\n";

	//s += "set output 'KalmanFilter\\1R_dfi.png'\n";
	s += "plot '" + way + "Kalman2.txt' with lines lc 'red' lw 2 title 'd{/Symbol f}/dt', ";
	s += "'" + way + "corKalman2.txt' with lines lc 'blue' lw 2 title '{d{/Symbol f}/dt}_{cor}'\n";

	//s += "set output 'KalmanFilter\\1R_dx.png'\n";
	s += "plot '" + way + "Kalman3.txt' with lines lc 'red' lw 2 title 'dx/dt', ";
	s += "'" + way + "corKalman3.txt' with lines lc 'blue' lw 2 title '{dx/dt}_{cor}'\n";
	s += "unset multiplot\n";

	s += "set output 'KalmanFilter\\1R_u.png'\n";
	s += "plot '" + way + "Kalmancontrol.txt' with lines lc 'green' lw 2 title 'u'\n";

	PLOT(s);

}

Matrix<double> RikkatiSolver(const Matrix<double>& A_, const Matrix<double>& B_, const Matrix<double>& Q_, const Matrix<double>& R_)
{
	Matrix<double> X(4, 4);

	ofstream file("ABQR.txt");

	for (int i = 0; i < A_.size(); i++)
	{
		for (int j = 0; j < A_[0].size(); j++)
			file << A_[i][j] << " ";
		file << endl;
	}
	file << endl;
	for (int i = 0; i < B_.size(); i++)
	{
		for (int j = 0; j < B_[0].size(); j++)
			file << B_[i][j] << " ";
		file << endl;
	}
	file << endl;
	for (int i = 0; i < Q_.size(); i++)
	{
		for (int j = 0; j < Q_[0].size(); j++)
			file << Q_[i][j] << " ";
		file << endl;
	}
	file << endl;
	for (int i = 0; i < R_.size(); i++)
	{
		for (int j = 0; j < R_[0].size(); j++)
			file << R_[i][j] << " ";
		file << endl;
	}

	system("python riccatisolve.py");

	ifstream fileX("X.txt");
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			fileX >> X[i][j];

	return X;
}

Vector<double> GetControl(Matrix<double>& A, const Vector<double>& B, Vector<complex<double>>& eigV)
{
	int n = A.size();
	//Матрица управляемости
	Matrix<double> ControlMatrix(n);
	//Единичная матрица и вектор e_n = (0, 0, ..., 1)^T
	Matrix<double> I = Single<double>(n);
	Vector<double> en(n), theta(n);
	en[3] = 1;

	//Заполняем матрицу управляемости
	ControlMatrix[0] = B;
	for (int i = 1; i < n; i++)
		ControlMatrix[i] = A * ControlMatrix[i - 1];
	ControlMatrix = TransposedMatrix(ControlMatrix);

	//Определитель матрицы управляемости
	PRINT(Determinant(ControlMatrix));

	theta = -1.0 * en * InverseMatrix(ControlMatrix);

	//Ищем коэффициенты полинома
	Vector<complex<double>> a(n + 1);

	a[1] = 1, a[0] = -eigV[0];
	for (int i = 2; i <= n; i++)
	{
		a[i] = a[i - 1];
		for (int j = i - 1; j > 0; j--)
			a[j] = a[j - 1] - a[j] * eigV[i - 1];
		a[0] = -a[0] * eigV[i - 1];
	}

	Matrix<double> Temp(n);
	//Ищем управление с помощью формулы Аккермана
	for (int i = 0; i <= n; i++)
	{
		Temp += a[i].real() * I;
		I *= A;
	}

	return theta * Temp;
}

Matrix<double> MatrixExp(Matrix<double> A)
{
	Matrix<double> GapMt = Single<double>(A.size()), res(A.size());
	double GapCoef = 1;
	for (int i = 0; i < 14; i++)
	{
		res += GapMt / GapCoef;
		GapMt *= A;
		GapCoef *= (i + 1.0);
	}
	return res;
}

Matrix<double> IntegrateMatrixExp(Matrix<double>& A, double h)
{
	Matrix<double> GapMt = Single<double>(A.size()), res(A.size());
	double GapCoef = h;
	for (int i = 1; i < 14; i++)
	{
		res += GapMt * GapCoef;
		GapMt *= A;
		GapCoef *= h / (i + 1.0);
	}
	return res;
}

tuple<Matrix<double>, Vector<double>> InicialSyst()
{
	Matrix<double> A(4);
	Vector<double> B(4);
	//Заполнили матрицу A и матрицу-столбец B
	double temp = (P.m + P.M) * (P.I + P.m * pow(P.l, 2)) - pow(P.m, 2) * pow(P.l, 2);
	A[2][0] = P.m * P.g * P.l * (P.m + P.M) / temp;
	A[2][1] = 0;
	A[2][2] = -P.Bp * (P.m + P.M) / temp;
	A[2][3] = -P.m * P.l * (P.Kf * P.Ks + P.Beq) / temp;
	B[2] = P.m * P.l * P.Kf / temp;

	A[3][0] = pow(P.m, 2) * pow(P.l, 2) * P.g / temp;
	A[3][1] = 0;
	A[3][2] = -P.m * P.l * P.Bp / temp;
	A[3][3] = -(P.I + P.m * pow(P.l, 2)) * (P.Kf * P.Ks + P.Beq) / temp;
	B[3] = P.Kf * (P.I + P.m * pow(P.l, 2)) / temp;

	A[0][2] = 1, A[1][3] = 1;
	return make_tuple(A, B);
}

Vector<double>& NonLinSystCoef(Vector<double>& p)
{
	p = Vector<double>(12);
	p[0] = P.m * P.l * P.Kf;
	p[1] = pow(P.m, 2) * pow(P.l, 2);
	p[2] = P.m * P.l * (P.Kf * P.Ks + P.Beq);
	p[3] = P.m * P.g * P.l * (P.m + P.M);
	p[4] = P.Bp * (P.m + P.M);
	p[5] = (P.I + P.m * pow(P.l, 2)) * P.Kf;
	p[6] = P.m * P.l * (P.I + P.m * pow(P.l, 2));
	p[7] = (P.Kf * P.Ks + P.Beq) * (P.I + P.m * pow(P.l, 2));
	p[8] = pow(P.m, 2) * P.g * pow(P.l, 2);
	p[9] = P.m * P.l * P.Bp;
	p[10] = (P.m + P.M) * (P.I + P.m * pow(P.l, 2));
	p[11] = pow(P.m, 2) * pow(P.l, 2);
	return p;
}

void GetImageDiscrSystem(DiscreteSystemParametrs& DS, Vector<double>& x_start)
{
	int MaxIt = 20000;
	double t;

	/*

		С кусочно-постоянным управлением

	*/

	t = 0;
	ORE solve(t, x_start, DS.h, 2);
	solve.SetSystParametrs(DS);
	solve.WriteToFileWithControl(P.tmax, "const", MaxIt);

	/*

		С импульсным управлением

	*/

	t = 0;
	solve = ORE(t, x_start, DS.h, 1);
	solve.SetSystParametrs(DS);
	solve.WriteToFileWithControl(P.tmax, "impulse", MaxIt);

	string s;
	s = setting;
	s += "set grid\n";

	s += "set output 'Discrete\\1D_R_fi.png'\n";
	s += "set multiplot layout 2,2\n";
	s += "plot '" + way + "const0.txt' with lines lc 'red' lw 2 title 'CONST-{/Symbol f}', ";
	s += "'" + way + "impulse0.txt' with lines lc 'blue' lw 2 title 'IMPULSE-{/Symbol f}'\n";

	//s += "set output 'Discrete\\1D_R_x.png'\n";
	s += "plot '" + way + "const1.txt' with lines lc 'red' lw 2 title 'CONST-x', ";
	s += "'" + way + "impulse1.txt' with lines lc 'blue' lw 2 title 'IMPUSLE-x'\n";

	//s += "set output 'Discrete\\1D_R_dfi.png'\n";
	s += "plot '" + way + "const2.txt' with lines lc 'red' lw 2 title 'CONST-d{/Symbol f}/dt', ";
	s += "'" + way + "impulse2.txt' with lines lc 'blue' lw 2 title 'IMPUSLE-d{/Symbol f}/dt'\n";

	//s += "set output 'Discrete\\1D_R_dx.png'\n";
	s += "plot '" + way + "const3.txt' with lines lc 'red' lw 2 title 'CONST-dx/dt', ";
	s += "'" + way + "impulse3.txt' with lines lc 'blue' lw 2 title 'IMPUSLE-dx/dt'\n";
	s += "unset multiplot\n";

	s += "set output 'Discrete\\1D_R_u.png'\n";
	s += "plot '" + way + "constcontrol.txt' with lines lc 'red' lw 2 title 'CONST-u', ";
	s += "'" + way + "impulsecontrol.txt' with points pt 7 lc 'blue' lw 2 title 'IMPUSLE-u'\n";
	PLOT(s);
}

void GetImageDiscrSystemObs(DiscreteSystemParametrs& DS, Vector<double>& x_start, Vector<double>& ksi_start)
{
	Vector<double> x(2 * x_start.size());
	int MaxIt = 20000;
	double t;

	/*

		Для кусочно-постоянного управления

	*/


	t = 0;
	for (int i = 0; i < x_start.size(); i++)
	{
		x[i] = ksi_start[i];
		x[i + 4] = x_start[i];
	}
	ORE solve(t, x, DS.h, 4);
	solve.SetSystParametrs(DS);
	solve.WriteToFileWithControl(P.tmax, "constOBS", MaxIt);

	string s;
	s = setting;
	s += "set grid\n";

	s += "set output 'Discrete\\1D_OBS_CONST_fi.png'\n";
	s += "set multiplot layout 2,2\n";
	s += "plot '" + way + "constOBS4.txt' with lines lc 'red' lw 2 title 'CONST-{/Symbol f}', ";
	s += "'" + way + "constOBS0.txt' with lines lc 'blue' lw 2 title 'obs-CONST-{/Symbol f}\n";

	//s += "set output 'Discrete\\1D_OBS_CONST_x.png'\n";
	s += "plot '" + way + "constOBS5.txt' with lines lc 'red' lw 2 title 'CONST-x', ";
	s += "'" + way + "constOBS1.txt' with lines lc 'blue' lw 2 title 'obs-CONST-x\n";

	//s += "set output 'Discrete\\1D_OBS_CONST_dfi.png'\n";
	s += "plot '" + way + "constOBS6.txt' with lines lc 'red' lw 2 title 'CONST-d{/Symbol f}/dt', ";
	s += "'" + way + "constOBS2.txt' with lines lc 'blue' lw 2 title 'obs-CONST-d{/Symbol f}/dt\n";

	//s += "set output 'Discrete\\1D_OBS_CONST_dx.png'\n";
	s += "plot '" + way + "constOBS7.txt' with lines lc 'red' lw 2 title 'CONST-dx/dt', ";
	s += "'" + way + "constOBS3.txt' with lines lc 'blue' lw 2 title 'obs-CONST-dx/dt\n";
	s += "unset multiplot\n";

	s += "set output 'Discrete\\1D_OBS_CONST_u.png'\n";
	s += "plot '" + way + "constOBScontrol.txt' with lines lc 'green' lw 2 title 'CONST-u'\n";
	PLOT(s);

	/*

		Для импульсного управления

	*/

	t = 0;
	for (int i = 0; i < x_start.size(); i++)
	{
		x[i] = ksi_start[i];
		x[i + 4] = x_start[i];
	}
	solve = ORE(t, x, DS.h, 3);
	solve.SetSystParametrs(DS);
	solve.WriteToFileWithControl(P.tmax, "impulseOBS", MaxIt);

	s = setting;
	s += "set grid\n";

	s += "set output 'Discrete\\1D_OBS_IMPUSLE_fi.png'\n";
	s += "set multiplot layout 2,2\n";
	s += "plot '" + way + "impulseOBS4.txt' with lines lc 'red' lw 2 title 'IMPUSLE-{/Symbol f}', ";
	s += "'" + way + "impulseOBS0.txt' with lines lc 'blue' lw 2 title 'obs-IMPUSLE-{/Symbol f}\n";

	//s += "set output 'Discrete\\1D_OBS_IMPUSLE_x.png'\n";
	s += "plot '" + way + "impulseOBS5.txt' with lines lc 'red' lw 2 title 'IMPUSLE-x', ";
	s += "'" + way + "impulseOBS1.txt' with lines lc 'blue' lw 2 title 'obs-IMPUSLE-x\n";

	//s += "set output 'Discrete\\1D_OBS_IMPUSLE_dfi.png'\n";
	s += "plot '" + way + "impulseOBS6.txt' with lines lc 'red' lw 2 title 'IMPUSLE-d{/Symbol f}/dt', ";
	s += "'" + way + "impulseOBS2.txt' with lines lc 'blue' lw 2 title 'obs-IMPUSLE-d{/Symbol f}/dt\n";

	//s += "set output 'Discrete\\1D_OBS_IMPUSLE_dx.png'\n";
	s += "plot '" + way + "impulseOBS7.txt' with lines lc 'red' lw 2 title 'IMPUSLE-dx/dt', ";
	s += "'" + way + "impulseOBS3.txt' with lines lc 'blue' lw 2 title 'obs-IMPUSLE-dx/dt\n";
	s += "unset multiplot\n";

	s += "set output 'Discrete\\1D_OBS_IMPUSLE_u.png'\n";
	s += "plot '" + way + "impulseOBScontrol.txt' with points pt 7 lc 'green' lw 2 title 'IMPUSLE-u'\n";
	PLOT(s);

}

void GetImageOriginalSystem(ContinuosSystemParametrs& SP, Vector<double>& pNon, Vector<double>& x_start)
{
	SystemState state;
	double eps = pow(10, -10), h = pow(10, -3);
	int MaxIt = 20000;

	//Заполняем файлы линейной системы
	state.t = 0, state.x = x_start;
	ODE4 solve(state.t, x_start, h, eps, 1);
	solve.SetStructParam(SP);
	solve.WriteToFileWithControl(P.tmax, "L", MaxIt);

	//Заполняем файлы нелинейной системы 
	state.t = 0, state.x = x_start;
	solve = ODE4(state.t, x_start, h, eps, 2);
	solve.SetStructParam(SP);
	solve.SetVecParam(pNon);
	solve.WriteToFileWithControl(P.tmax, "NL", MaxIt);

	string s;
	s = setting;
	s += "set grid\n";

	s += "set output 'Continuos\\1R_fi.png'\n";
	s += "set multiplot layout 2,2\n";
	s += "plot '" + way + "L0.txt' with lines lc 'red' lw 2 title 'lin-{/Symbol f}', ";
	s += "'" + way + "NL0.txt' with lines lc 'blue' lw 2 title 'nonlin-{/Symbol f}'\n";

	//s += "set output 'Continuos\\1R_x.png'\n";
	s += "plot '" + way + "L1.txt' with lines lc 'red' lw 2 title 'lin-x', ";
	s += "'" + way + "NL1.txt' with lines lc 'blue' lw 2 title 'nonlin-x'\n";

	//s += "set output 'Continuos\\1R_dfi.png'\n";
	s += "plot '" + way + "L2.txt' with lines lc 'red' lw 2 title 'lin-d{/Symbol f}/dt', ";
	s += "'" + way + "NL2.txt' with lines lc 'blue' lw 2 title 'nonlin-d{/Symbol f}/dt'\n";

	//s += "set output 'Continuos\\1R_dx.png'\n";
	s += "plot '" + way + "L3.txt' with lines lc 'red' lw 2 title 'lin-dx/dt', ";
	s += "'" + way + "NL3.txt' with lines lc 'blue' lw 2 title 'nonlin-dx/dt'\n";
	s += "unset multiplot\n";

	s += "set output 'Continuos\\1R_u.png'\n";
	s += "plot '" + way + "Lcontrol.txt' with lines lc 'red' lw 2 title 'lin-u', ";
	s += "'" + way + "NLcontrol.txt' with lines lc 'blue' lw 2 title 'nonlin-u'\n";
	PLOT(s);

}

void GetImageLinObserver(ContinuosSystemParametrs& SP, Vector<double>& x_start)
{
	int n = 8;
	SystemState state;
	double eps = pow(10, -10), h = pow(10, -3);
	int MaxIt = 20000;

	state.t = 0, state.x = x_start;
	ODE4 solve(state.t, x_start, h, eps, 3);
	solve.SetStructParam(SP);
	solve.WriteToFileWithControl(P.tmax, "LOBS", MaxIt);

	string s;
	s = setting;
	s += "set grid\n";

	s += "set output 'Continuos\\1L_OBS_fi.png'\n";
	s += "set multiplot layout 2,2\n";
	s += "plot '" + way + "LOBS4.txt' with lines lc 'red' lw 2 title 'lin-{/Symbol f}', ";
	s += "'" + way + "LOBS0.txt' with lines lc 'blue' lw 2 title 'obs-lin-{/Symbol f}'\n";

	//s += "set output 'Continuos\\1L_OBS_x.png'\n";
	s += "plot '" + way + "LOBS5.txt' with lines lc 'red' lw 2 title 'lin-x', ";
	s += "'" + way + "LOBS1.txt' with lines lc 'blue' lw 2 title 'obs-lin-x'\n";

	//s += "set output 'Continuos\\1L_OBS_dfi.png'\n";
	s += "plot '" + way + "LOBS6.txt' with lines lc 'red' lw 2 title 'lin-d{/Symbol f}/dt', ";
	s += "'" + way + "LOBS2.txt' with lines lc 'blue' lw 2 title 'obs-lin-d{/Symbol f}/dt'\n";

	//s += "set output 'Continuos\\1L_OBS_dx.png'\n";
	s += "plot '" + way + "LOBS7.txt' with lines lc 'red' lw 2 title 'lin-dx/dt', ";
	s += "'" + way + "LOBS3.txt' with lines lc 'blue' lw 2 title 'obs-lin-dx/dt'\n";
	s += "unset multiplot\n";

	s += "set output 'Continuos\\1L_OBS_u.png'\n";
	s += "plot '" + way + "LOBScontrol.txt' with lines lc 'green' lw 2 title 'obs-lin-u'\n";
	PLOT(s);

}

void GetImageNonLinObserver(ContinuosSystemParametrs& SP, Vector<double>& pNon, Vector<double>& x_start)
{
	int n = 8;
	SystemState state;
	double eps = pow(10, -10), h = pow(10, -3);
	int MaxIt = 20000;

	state.t = 0, state.x = x_start;
	ODE4 solve(state.t, x_start, h, eps, 4);
	solve.SetStructParam(SP);
	solve.SetVecParam(pNon);
	solve.WriteToFileWithControl(P.tmax, "NLOBS", MaxIt);

	string s;
	s = setting;
	s += "set grid\n";

	s += "set output 'Continuos\\1NL_OBS_fi.png'\n";
	s += "set multiplot layout 2,2\n";
	s += "plot '" + way + "NLOBS4.txt' with lines lc 'red' lw 2 title 'nonlin-{/Symbol f}', ";
	s += "'" + way + "NLOBS0.txt' with lines lc 'blue' lw 2 title 'obs-nonlin-{/Symbol f}'\n";

	//s += "set output 'Continuos\\1NL_OBS_x.png'\n";
	s += "plot '" + way + "NLOBS5.txt' with lines lc 'red' lw 2 title 'nonlin-x', ";
	s += "'" + way + "NLOBS1.txt' with lines lc 'blue' lw 2 title 'obs-nonlin-x'\n";

	//s += "set output 'Continuos\\1NL_OBS_dfi.png'\n";
	s += "plot '" + way + "NLOBS6.txt' with lines lc 'red' lw 2 title 'nonlin-d{/Symbol f}/dt', ";
	s += "'" + way + "NLOBS2.txt' with lines lc 'blue' lw 2 title 'obs-nonlin-d{/Symbol f}/dt'\n";

	//s += "set output 'Continuos\\1NL_OBS_dx.png'\n";
	s += "plot '" + way + "NLOBS7.txt' with lines lc 'red' lw 2 title 'nonlin-dx/dt', ";
	s += "'" + way + "NLOBS3.txt' with lines lc 'blue' lw 2 title 'obs-nonlin-dx/dt'\n";
	s += "unset multiplot\n";

	s += "set output 'Continuos\\1NL_OBS_u.png'\n";
	s += "plot '" + way + "NLOBScontrol.txt' with lines lc 'green' lw 2 title 'obs-nonlin-u'\n";
	PLOT(s);

}

Vector<double> operator+(const Vector<double>& v1, const arma::vec& v2)
{
	Vector<double> res(v1.size());
	for (int i = 0; i < v1.size(); i++)
		res[i] = v1[i] + v2[i];
	return res;
}

Vector<double> operator+(const arma::vec& v1, const Vector<double>& v2)
{
	Vector<double> res(v2.size());
	for (int i = 0; i < v2.size(); i++)
		res[i] = v1[i] + v2[i];
	return res;
}

Matrix<double> mat_to_Matrix(const arma::mat& mt)
{
	int n = static_cast<int>(mt.n_rows), m = static_cast<int>(mt.n_cols);
	Matrix<double> res(n, m);

	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			res[i][j] = mt(i, j);
	return res;
}

arma::mat Matrix_to_mat(const Matrix<double>& mt)
{
	int n = mt.size(), m = mt[0].size();
	arma::mat res(mt.size(), mt[0].size());

	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			res(i, j) = mt[i][j];
	return res;
}