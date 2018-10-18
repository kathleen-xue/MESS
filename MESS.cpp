#include <iostream>
#include <fstream>
#include <string>
#include <iostream>
#include <cmath> 
#include <vector>
#include <sstream> 
#include <queue> 
#include <cstdlib>  
#include <algorithm>
#include <cstdlib>
#include <iomanip>
using namespace std;

//define costs properly
//transportation costs: coordinate/location: generate random MxM symmetric matrix, and each cell will be distance from microgrid i to microgrid j

int costToCharge(int start, int end) {
	return (end - start)*3;
}

int costToStay(int M, int T, int start, int end) {
	return -1*(start - end)*20;
}

vector<vector<vector<long double> > > cheapestPathDP(long double& currMinCost, int T, int M, vector<pair<pair<int, int>, int> >& path, vector<vector<double> > transportCost, vector<vector<vector<long double> > >& dpCUBE) {
//This function finds the maximum reduction of cost for a battery to
//end up at any microgrid at any charge state on day T. If we iterate 
//this function through x batteries, the final dpCUBE will output 
//maximum reduction of cost from all the batteries for each
//microgrid after T days (this will be seen on the last level of the cube).
	
	for(int i = 1; i < T; i++) { //cycling through each day
		pair<int, int> micro;
		for(int j = 1; j < M; j++) { //cycling through the microgrids of current day
			for(int k = 0; k <= 100; k+=25) { //cycling through each of the charge states of each microgrid of current day
				for(int l = 0; l < M; l++) { //cycling through each microgrid of previous day
					for(int m = 0; m <= 100; m+=25) { //cycling through each of the charge states of previous day
						for(int n = m; n <= 100; n+=25) { //cycling through each of the *possible* charges of current day
							if(n != m) {
								long double curr = dpCUBE[i-1][l][m/25] + costToCharge(m, n) + transportCost[j][M] + transportCost[M][l] + costToStay(j, i, k, n);
								if(curr < 0) curr = 0;
								if(dpCUBE[i][j][k/25] >= curr) {
									currMinCost = dpCUBE[i-1][l][m/25];
									micro.first = l;
									micro.second = n;
								}
								dpCUBE[i][j][k/25] = min(dpCUBE[i][j][k/25], curr);
							}
							else { 
								long double curr = dpCUBE[i-1][l][m/25] + transportCost[j][l];
								if(curr < 0) curr = 0;
								if(dpCUBE[i][j][k/25] >= curr) {
									//cout << costToStay(j, i, n, n) << endl;
									currMinCost = dpCUBE[i-1][l][m/25];
									micro.first = l;
									micro.second = m;
								}
								dpCUBE[i][j][k/25] = min(dpCUBE[i][j][k/25], curr);
							}
						}
					}
				}
			}
		}
		path.push_back(make_pair(make_pair(micro.first, micro.second), currMinCost));
	}

	return dpCUBE;
}

int main () {
	int T = 0;
	int M = 0;
	int numBatteries = 0;
	cout << "Enter # days: " << endl;
	cin >> T;
	cout << "Enter # microgrids: " << endl;
	cin >> M;
	cout << "Enter # batteries: " << endl;
	cin >> numBatteries;

	vector<vector<double> > transportCost(M+1, vector<double>(M+1));

	for(int i = 0; i < M+1; i++) {
		for(int j = i; j < M+1; j++) {
			if(i == j) transportCost[i][j] = 0;
			else {
				double currRand = (double)((rand() % 5)) + 1;
				transportCost[i][j] = currRand;
				transportCost[j][i] = currRand;
			}
		}
	}

	cout << "Transportation Costs: (between microgrid i and microgrid j)" << endl;
	for(int i = 0; i < M+1; i++) {
		for(int j = 0; j < M+1; j++) {
			cout << setw(2) << transportCost[i][j] << setw(2);
		}
		cout << endl;
	}
	cout << endl;

	vector<vector<vector<long double> > > dpCUBE;

	for(int i = 0; i < T; i++) {
		vector<vector<long double> > curr1; 
		for(int j = 0; j < M; j++) {
			vector<long double> curr2;
			for(int k = 0; k <= 100; k+=25) {
				curr2.push_back(10000);
				//cout << curr2[k/25] << " ";
			}
			//cout << endl;
			curr1.push_back(curr2);
		}
		dpCUBE.push_back(curr1);
	}

	for(int i = 0; i < M; i++) {
		for(int j = 0; j <= 100; j+=25) {
			dpCUBE[0][i][j/25] = 100000;
		}
	}

	long double currMinCost = 100000;

for(int k = 0; k < numBatteries; k++) {
	vector<pair<pair<int, int>, int> > path;
	vector<vector<vector<long double> > > dpCUBE1 = cheapestPathDP(currMinCost, T, M, path, transportCost, dpCUBE);
	cout << "Minimum cost path for battery " << k+1 << ": " << endl;
	for(int i = 0; i < path.size(); i++) {
		cout << setw(2) << i+1 << setw(2) << " | microgrid " << setw(2) << path[i].first.first << setw(2) << 
		" | state: " << setw(2) << path[i].first.second << setw(2) << " | cost: " << setw(2) << path[i].second << setw(2) << endl;
	}
}


	return 0;
}










