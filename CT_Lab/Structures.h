#pragma once

#ifndef STRUCTURES_H_
#define STRUCTURES_H_

#include <algorithm>
#include <iostream>
#include <fstream>
#include <complex>
#include <math.h>
#include <string>
#include <tuple>
#include <list>
#include "NextGenMatrix.h"

struct DiscreteSystemParametrs
{
	Matrix<double> Ad, ABTd_const, BTd_const, ABTd_impulse, BTd_impulse, LCd;
	Vector<double> Bd_const, Bd_impulse, Ld, thetad_const, thetad_impulse;
	double h = 0;
};

struct ContinuosSystemParametrs
{
	Matrix<double> A, BTheta, LC, Ac, RBX;
	Vector<double> theta, L, B, C;
};

struct SystemState
{
	Vector<double> x;
	double t = 0;
};

#endif