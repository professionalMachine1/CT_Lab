#pragma once
class ODE
{
};

#pragma once

#ifndef ODE_H
#define ODE_H

#include "Structures.h"

class ODE4
{
private:
	int Number;
	double eps, tprev, h;
	SystemState state;
	Vector<double> X, xprev, Xprev, p;
	ContinuosSystemParametrs SP;
public:
	ODE4();
	ODE4(double t_, Vector<double> x_, double h_, double eps_, int number);

	SystemState& ComplitWithoutControl();
	SystemState& ComplitWithControl(bool boolf = true);

	Vector<double> RunKut4(double t_, Vector<double>& v_, double h_);
	Vector<double> systValue(double t_, Vector<double>& v_);

	void SetVecParam(Vector<double>& par) { p = par; }
	void SetStructParam(ContinuosSystemParametrs& par) { SP = par; }

	void SetCoordinate(double t_) { state.t = t_;; }
	void SetStep(double h_) { h = h_; }
	void SetSystemCondition(Vector<double>& x_) { state.x = x_; }

	SystemState& GetData() { return state; }
	void BackStep(bool f = true);

	void WriteToFile(double max_t, std::string name, int MaxIt = 100000);
	void WriteToFileWithControl(double max_t, std::string name, int MaxIt = 100000);
	void WriteToFileLQR(double max_t, std::string name, int MaxIt = 100000);

};

#endif