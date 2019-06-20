#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "PSI1D_Backscattering.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
//                                                   MAIN CODE
// ---------------------------------------------------------------------------------------------------------------------
int main (int argc, char* argv[])
{
	double rn[90][500], re[90][500];
	int an1;         // atomic number of incident atomic
	int am1;          // atomic mass of incident atomic (amu)
	int ne;             // number of constituent elements in the target.
	vector<int> an2;	// array for atomic numbers of the constituents.
	vector<int> nw;     // array for relative numbers of the constituents.
	double energy, theta, ee, rnion, reion;
	ofstream ofile_rn, ofile_re;

/*
	// D -> W
	an1 = 1;
	am1 = 2;
	ne = 1;
	an2.push_back(74);
	nw.push_back(1);
*/

/*
	// C -> W
	an1 = 6;
	am1 = 12;
	ne = 1;
	an2.push_back(74);
	nw.push_back(1);
*/

/*
	// Fe -> C
	an1 = 26;
	am1 = 56;
	ne = 1;
	an2.push_back(6);
	nw.push_back(1);
*/

/*
	// D -> C
	an1 = 1;
	am1 = 2;
	ne = 1;
	an2.push_back(6);
	nw.push_back(1);
*/


	// C -> C
	an1 = 6;
	am1 = 12;
	ne = 1;
	an2.push_back(6);
	nw.push_back(1);


	PSI1D_Backscattering bs = PSI1D_Backscattering( an1, am1, ne, an2, nw );
	for(int i = 1; i <= 90; i++)
	{
		for(int j = 1; j <= 500; j++)
		{
			theta = 1.0 * i;
			energy = 1.0 * j;
			bs.scatter(rnion, reion, theta, energy);
			rn[i-1][j-1] = rnion;
			re[i-1][j-1] = reion;
		}
	}

	ofile_rn.open("rn.txt");
	ofile_re.open("re.txt");
	for(int i = 1; i <= 90; i++)
	{
		for(int j = 1; j <= 500; j++)
		{
			ofile_rn<<setw(15)<<rn[i-1][j-1];
			ofile_re<<setw(15)<<re[i-1][j-1];
		}
		ofile_rn<<endl;
		ofile_re<<endl;
	}
	ofile_rn.close();
	ofile_re.close();

	theta = 0.0;
	energy = 5.0;
	bs.scatter(rnion, reion, theta, energy);
	cout<<rnion<<" "<<reion<<"  "<<theta<<"  "<<energy<<endl;

	return 0;
	
	int i= 0;
	while(rnion < 1.0)
	{
		bs.scatter(rnion, reion, theta, energy);
		i++;
	}
	cout<<rnion<<" "<<reion<<"  "<<theta<<"  "<<energy<<"  "<<i<<endl;


}
