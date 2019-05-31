#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "PSI1D_Sputtering.h"

using namespace std;

int main (int argc, char* argv[])
{
	double rn[90][500];
	double an1;
	double an2;
	double am1;
	double am2;
	double es;
	double density;

	double energy, theta, ee, rnion, reion;
	ofstream ofile_rn;

/*
	// Ar -> W
	an1 = 18.0;
	am1 = 39.95;
	an2 = 74.0;
	am2 = 183.85;
	es = 11.75;
	density = 19.35;
*/

/*
	// W -> W
	an1 = 74.0;
	am1 = 183.85;
	an2 = 74.0;
	am2 = 183.85;
	es = 11.75;
	density = 19.35;
*/


	// D -> W
	an1 = 1.0;
	am1 = 2.0;
	an2 = 74.0;
	am2 = 183.85;
	es = 11.75;
	density = 19.35;


	PSI1D_Sputtering st = PSI1D_Sputtering( an1, am1, an2, am2, es, density );
	for(int i = 1; i <= 90; i++)
	{
		for(int j = 1; j <= 500; j++)
		{
			theta = 1.0 * i;
			energy = 1.0 * j;
			rn[i-1][j-1] = st.phy_sput_yield(energy, theta);
		}
	}

	ofile_rn.open("rn.txt");
	for(int i = 1; i <= 90; i++)
	{
		for(int j = 1; j <= 500; j++)
		{
			ofile_rn<<setw(20)<<rn[i-1][j-1];
		}
		ofile_rn<<endl;
	}
	ofile_rn.close();
}
