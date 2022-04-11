#include "ORE.h"

ORE::ORE()
{
	h = 0;
	Number = 0;
	cur = SystemState();
	DS = DiscreteSystemParametrs();
}

ORE::ORE(double t_, Vector<double>& x_, double h_, int Number_)
{
	h = h_, Number = Number_;
	cur.x = x_, cur.t = t_;
}

SystemState& ORE::nextState()
{
	prev = cur;
	switch (Number)
	{
	case 1:
		cur.x = DS.ABTd_impulse * cur.x;
		cur.t += h;
		break;
	case 2:
		cur.x = DS.ABTd_const * cur.x;
		cur.t += h;
		break;
	case 3:
		// 0, ..., 3 - ksi; 4, ..., 7 - x
		for (int i = 0; i < 4; i++)
		{
			cur.x[i] = DS.ABTd_impulse[i][0] * prev.x[0] + DS.ABTd_impulse[i][1] * prev.x[1]
				+ DS.ABTd_impulse[i][2] * prev.x[2] + DS.ABTd_impulse[i][3] * prev.x[3];
			cur.x[i] += DS.LCd[i][0] * (prev.x[4] - prev.x[0]) + DS.LCd[i][1] * (prev.x[5] - prev.x[1])
				+ DS.LCd[i][2] * (prev.x[6] - prev.x[2]) + DS.LCd[i][3] * (prev.x[7] - prev.x[3]);
		}
		for (int i = 0; i < 4; i++)
		{
			cur.x[4 + i] = DS.Ad[i][0] * prev.x[4] + DS.Ad[i][1] * prev.x[5]
				+ DS.Ad[i][2] * prev.x[6] + DS.Ad[i][3] * prev.x[7];
			cur.x[4 + i] += DS.BTd_impulse[i][0] * prev.x[0] + DS.BTd_impulse[i][1] * prev.x[1]
				+ DS.BTd_impulse[i][2] * prev.x[2] + DS.BTd_impulse[i][3] * prev.x[3];
		}
		cur.t += h;
		break;
	case 4:
		// 0, ..., 3 - ksi; 4, ..., 7 - x
		for (int i = 0; i < 4; i++)
		{
			cur.x[i] = DS.ABTd_const[i][0] * prev.x[0] + DS.ABTd_const[i][1] * prev.x[1]
				+ DS.ABTd_const[i][2] * prev.x[2] + DS.ABTd_const[i][3] * prev.x[3];
			cur.x[i] += DS.LCd[i][0] * (prev.x[4] - prev.x[0]) + DS.LCd[i][1] * (prev.x[5] - prev.x[1])
				+ DS.LCd[i][2] * (prev.x[6] - prev.x[2]) + DS.LCd[i][3] * (prev.x[7] - prev.x[3]);
		}
		for (int i = 0; i < 4; i++)
		{
			cur.x[4 + i] = DS.Ad[i][0] * prev.x[4] + DS.Ad[i][1] * prev.x[5]
				+ DS.Ad[i][2] * prev.x[6] + DS.Ad[i][3] * prev.x[7];
			cur.x[4 + i] += DS.BTd_const[i][0] * prev.x[0] + DS.BTd_const[i][1] * prev.x[1]
				+ DS.BTd_const[i][2] * prev.x[2] + DS.BTd_const[i][3] * prev.x[3];
		}
		cur.t += h;
		break;
	default:
		break;
	}
	return cur;
}

void ORE::WriteToFile(double max_t, std::string name, int MaxIt)
{
	Vector<std::ofstream> data(cur.x.size());
	for (int i = 0; i < cur.x.size(); i++)
		data[i].open("putHere\\" + name + std::to_string(i) + ".txt", std::ios_base::out);
	for (int i = 0; (i < MaxIt) && (cur.t < max_t); i++)
	{
		for (int j = 0; j < cur.x.size(); j++)
			data[j] << cur.t << " " << cur.x[j] << std::endl;
		nextState();
	}
	for (int i = 0; i < cur.x.size(); i++)
		data[i].close();
}

void ORE::WriteToFileWithControl(double max_t, std::string name, int MaxIt)
{
	Vector<std::ofstream> data(cur.x.size());
	for (int i = 0; i < cur.x.size(); i++)
		data[i].open("putHere\\" + name + std::to_string(i) + ".txt", std::ios_base::out);
	std::ofstream udata("putHere\\" + name + "control.txt", std::ios_base::out);

	double temp = 0;

	if (Number % 2 == 1)
	{
		for (int i = 0; (i < MaxIt) && (cur.t < max_t); i++)
		{
			for (int j = 0; j < cur.x.size(); j++)
				data[j] << cur.t << " " << cur.x[j] << std::endl;
			temp = cur.x[0] * DS.thetad_impulse[0] + cur.x[1] * DS.thetad_impulse[1]
				+ cur.x[2] * DS.thetad_impulse[2] + cur.x[3] * DS.thetad_impulse[3];
			udata << cur.t << " " << temp << std::endl;
			nextState();
		}
	}
	if (Number % 2 == 0)
	{
		for (int i = 0; (i < MaxIt) && (cur.t < max_t); i++)
		{
			for (int j = 0; j < cur.x.size(); j++)
				data[j] << cur.t << " " << cur.x[j] << std::endl;
			temp = cur.x[0] * DS.thetad_const[0] + cur.x[1] * DS.thetad_const[1]
				+ cur.x[2] * DS.thetad_const[2] + cur.x[3] * DS.thetad_const[3];
			udata << cur.t << " " << temp << std::endl;
			nextState();
		}
	}
	for (int i = 0; i < cur.x.size(); i++)
		data[i].close();
	udata.close();
}
