#include "ODE.h"

ODE4::ODE4()
{
	Number = 0;
	eps = 0, h = 0, tprev = 0;
}

ODE4::ODE4(double t_, Vector<double> x_, double h_, double eps_, int number)
{
	state.t = t_;
	state.x = x_, xprev = state.x;
	X = Vector<double>(state.x.size()), Xprev = Vector<double>(state.x.size());
	h = h_, Number = number, eps = eps_;
	p = Vector<double>(0);
}

Vector<double> ODE4::RunKut4(double t_, Vector<double>& v_, double h_)
{
	int k = state.x.size();
	Vector<double> k1(k), k2(k), k3(k), k4(k), res(k);
	k1 = systValue(t_, v_);
	k2 = v_ + k1 * h_ / 2, k2 = systValue(t_ + h_ / 2, k2);
	k3 = v_ + k2 * h_ / 2, k3 = systValue(t_ + h_ / 2, k3);
	k4 = v_ + k3 * h_, k4 = systValue(t_ + h_, k4);
	res = v_ + (k1 + k2 * 2 + k3 * 2 + k4) * h_ / 6;
	return res;
}

Vector<double> ODE4::systValue(double t_, Vector<double>& v_)
{
	Vector<double> res(v_.size());
	Vector<double> uU(SP.RBX.size());
	double u, denum;
	switch (Number)
	{
	case 1:
		res = SP.Ac * v_;
		break;
	case 2:
		u = SP.theta * v_;
		denum = p[10] - p[11] * pow(cos(v_[0]), 2);
		res[0] = v_[2];
		res[1] = v_[3];
		res[2] = (p[0] * cos(v_[0]) * u - p[1] * cos(v_[0]) * sin(v_[0]) * pow(v_[2], 2)
			- p[2] * cos(v_[0]) * v_[3] + p[3] * sin(v_[0]) - p[4] * v_[2]) / denum;
		res[3] = (p[5] * u - p[6] * sin(v_[0]) * pow(v_[2], 2) - p[7] * v_[3]
			+ p[8] * cos(v_[0]) * sin(v_[0]) - p[9] * cos(v_[0]) * v_[2]) / denum;
		break;
	case 3:
		// 0, ..., 3 - ksi; 4, ..., 7 - x
		for (int i = 0; i < 4; i++)
		{
			res[i] = SP.Ac[i][0] * v_[0] + SP.Ac[i][1] * v_[1]
				+ SP.Ac[i][2] * v_[2] + SP.Ac[i][3] * v_[3];
			res[i] += SP.LC[i][0] * (v_[4] - v_[0]) + SP.LC[i][1] * (v_[5] - v_[1])
				+ SP.LC[i][2] * (v_[6] - v_[2]) + SP.LC[i][3] * (v_[7] - v_[3]);
		}
		for (int i = 0; i < 4; i++)
		{
			res[4 + i] = SP.A[i][0] * v_[4] + SP.A[i][1] * v_[5]
				+ SP.A[i][2] * v_[6] + SP.A[i][3] * v_[7];
			res[4 + i] += SP.BTheta[i][0] * v_[0] + SP.BTheta[i][1] * v_[1]
				+ SP.BTheta[i][2] * v_[2] + SP.BTheta[i][3] * v_[3];
		}
		break;
	case 4:
		u = (SP.theta[0] * v_[0] + SP.theta[1] * v_[1] + SP.theta[2] * v_[2] + SP.theta[3] * v_[3]);
		denum = p[10] - p[11] * pow(cos(v_[0]), 2);
		res[0] = v_[2];
		res[1] = v_[3];
		res[2] = (p[0] * cos(v_[0]) * u - p[1] * cos(v_[0]) * sin(v_[0]) * pow(v_[2], 2)
			- p[2] * cos(v_[0]) * v_[3] + p[3] * sin(v_[0]) - p[4] * v_[2]) / denum;
		res[3] = (p[5] * u - p[6] * sin(v_[0]) * pow(v_[2], 2) - p[7] * v_[3]
			+ p[8] * cos(v_[0]) * sin(v_[0]) - p[9] * cos(v_[0]) * v_[2]) / denum;
		for (int i = 0; i < 4; i++)
		{
			res[i] += SP.LC[i][0] * (v_[4] - v_[0]) + SP.LC[i][1] * (v_[5] - v_[1])
				+ SP.LC[i][2] * (v_[6] - v_[2]) + SP.LC[i][3] * (v_[7] - v_[3]);
		}
		denum = p[10] - p[11] * pow(cos(v_[4]), 2);
		res[4] = v_[6];
		res[5] = v_[7];
		res[6] = (p[0] * cos(v_[4]) * u - p[1] * cos(v_[4]) * sin(v_[4]) * pow(v_[6], 2)
			- p[2] * cos(v_[4]) * v_[7] + p[3] * sin(v_[4]) - p[4] * v_[6]) / denum;
		res[7] = (p[5] * u - p[6] * sin(v_[4]) * pow(v_[6], 2) - p[7] * v_[7]
			+ p[8] * cos(v_[4]) * sin(v_[4]) - p[9] * cos(v_[4]) * v_[6]) / denum;
		break;
	case 5:
		uU = SP.RBX * v_, u = uU[0];
		denum = p[10] - p[11] * pow(cos(v_[0]), 2);
		res[0] = v_[2];
		res[1] = v_[3];
		res[2] = (p[0] * cos(v_[0]) * u - p[1] * cos(v_[0]) * sin(v_[0]) * pow(v_[2], 2)
			- p[2] * cos(v_[0]) * v_[3] + p[3] * sin(v_[0]) - p[4] * v_[2]) / denum;
		res[3] = (p[5] * u - p[6] * sin(v_[0]) * pow(v_[2], 2) - p[7] * v_[3]
			+ p[8] * cos(v_[0]) * sin(v_[0]) - p[9] * cos(v_[0]) * v_[2]) / denum;
		break;
	default:
		return 0;
	}
	return res;
}

SystemState& ODE4::ComplitWithoutControl()
{
	state.x = RunKut4(state.t, state.x, h);
	state.t += h;
	return state;
}

SystemState& ODE4::ComplitWithControl(bool boolf)
{
	//count[0] - число делений, count[1] - число умножений
	Vector<double> S(state.x.size());
	double MaxS;

	xprev = state.x, Xprev = X;
	// V - Vn+1 вычисление с шагом h/2
	X = RunKut4(state.t, state.x, h / 2.0), X = RunKut4(state.t + h / 2.0, X, h / 2.0);
	//v - Vn+1 вычисленное с шагом h
	state.x = RunKut4(state.t, state.x, h);

	S = (X - state.x) / 15.0;
	//Берем модуль погрешности
	S = abs(S);
	//Находим max S
	MaxS = S[MaxValue(S)];

	tprev = state.t, state.t += h;
	if ((boolf) && (MaxS < eps / 32.0))
		h *= 2.0;
	else
		if (MaxS > eps)
		{
			BackStep();
			ComplitWithControl();
		}

	return state;
}

void ODE4::BackStep(bool f)
{
	state.t = tprev;
	state.x = xprev, X = Xprev;
	if (f)
		h *= 0.5;
}

void ODE4::WriteToFile(double max_t, std::string name, int MaxIt)
{
	Vector<std::ofstream> data(state.x.size());
	for (int i = 0; i < state.x.size(); i++)
		data[i].open("putHere\\" + name + std::to_string(i) + ".txt", std::ios_base::out);
	for (int i = 0; (i < MaxIt) && (state.t < max_t); i++)
	{
		for (int j = 0; j < state.x.size(); j++)
			data[j] << state.t << " " << state.x[j] << std::endl;
		ComplitWithControl();
	}
	for (int i = 0; i < state.x.size(); i++)
		data[i].close();
}

void ODE4::WriteToFileWithControl(double max_t, std::string name, int MaxIt)
{
	Vector<std::ofstream> data(state.x.size());
	for (int i = 0; i < state.x.size(); i++)
		data[i].open("putHere\\" + name + std::to_string(i) + ".txt", std::ios_base::out);
	std::ofstream udata("putHere\\" + name + "control.txt", std::ios_base::out);

	double temp = 0;

	for (int i = 0; (i < MaxIt) && (state.t < max_t); i++)
	{
		for (int j = 0; j < state.x.size(); j++)
			data[j] << state.t << " " << state.x[j] << std::endl;
		temp = (SP.theta[0] * state.x[0] + SP.theta[1] * state.x[1]
			+ SP.theta[2] * state.x[2] + SP.theta[3] * state.x[3]);
		udata << state.t << " " << temp << std::endl;
		ComplitWithControl();
	}
	for (int i = 0; i < state.x.size(); i++)
		data[i].close();
	udata.close();
}

void ODE4::WriteToFileLQR(double max_t, std::string name, int MaxIt)
{
	Vector<std::ofstream> data(state.x.size()), udata(SP.RBX.size());
	for (int i = 0; i < state.x.size(); i++)
		data[i].open("putHere\\" + name + std::to_string(i) + ".txt", std::ios_base::out);
	for (int i = 0; i < SP.RBX.size(); i++)
		udata[i].open("putHere\\" + name + "control" + std::to_string(i) + ".txt", std::ios_base::out);

	Vector<double> temp = 0;

	for (int i = 0; (i < MaxIt) && (state.t < max_t); i++)
	{
		for (int j = 0; j < state.x.size(); j++)
			data[j] << state.t << " " << state.x[j] << std::endl;
		temp = SP.RBX * state.x;
		for (int j = 0; j < SP.RBX.size(); j++)
			udata[j] << state.t << " " << temp[j] << std::endl;
		ComplitWithControl();
	}
	for (int i = 0; i < state.x.size(); i++)
		data[i].close();
	for (int i = 0; i < SP.RBX.size(); i++)
		udata[i].close();
}
