#pragma once

#ifndef ORE_H
#define ORE_H

#include "Structures.h"

class ORE
{
private:
	int Number;
	double h;
	SystemState cur, prev;
	DiscreteSystemParametrs DS;
public:
	ORE();
	ORE(double t_, Vector<double>& x_, double h_, int Number_);

	SystemState& nextState();
	SystemState& GetCurState() { return cur; }

	void SetSystParametrs(DiscreteSystemParametrs& DS_) { DS = DS_; }

	void WriteToFile(double max_t, std::string name, int MaxIt = 100000);
	void WriteToFileWithControl(double max_t, std::string name, int MaxIt = 100000);

};

#endif