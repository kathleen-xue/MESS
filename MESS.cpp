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
#include <map>
#include <cstdlib>
#include <iomanip>
using namespace std;

map<pair<int ,int>, int > demand;

//define costs properly
//transportation costs: coordinate/location: generate random MxM symmetric matrix, and each cell will be distance from microgrid i to microgrid j

int costToCharge(int start, int end) {
	return (end - start)*3;
}

int costToStay(int M, int T, int end, int start) {

	int d = demand[make_pair(T, M)];

	return -1 * sqrt(d) * ((d - max(d - (start - end), 0))*3);
}

void cheapestPathDP(int T, int M, vector<pair<pair<int, int>, int> >& path, vector<vector<long long> > transportCost, vector<vector<vector<long long> > >& dpCUBE) {
//This function finds the maximum reduction of cost for a battery to
//end up at any microgrid at any charge state on day T. If we iterate 
//this function through x batteries, the final dpCUBE will output 
//maximum reduction of cost from all the batteries for each
//microgrid after T days (this will be seen on the last level of the cube).

//FINDING MINIMUM COSTS	
	for(int i = 1; i <= T; i++) { //cycling through each day
		pair<int, int> micro;
		for(int j = 0; j < M; j++) { //cycling through the microgrids of current day
			for(int l = 0; l < M; l++) { //cycling through each microgrid of previous day
				for(int m = 0; m <= 100; m+=25) { //cycling through each of the charge states of previous day
					for(int n = m; n <= 100; n+=25) { //cycling through each of the *possible* charges of current day
						for(int k = n; k >= 0; k-=25) { //cycling through each of the charge states of each microgrid of current day
							if(n != m) {
								long long curr = dpCUBE[i-1][l][m/25] + costToCharge(m, n) + transportCost[j][M] + transportCost[M][l] + costToStay(j, i, k, n);
								//cout << i << " " << j << " "  << k << " " << l << " " << m << " " << curr << endl;
								// if(curr < 0) curr = 0;
								//if(dpCUBE[i][j][k/25] >= curr) {
									//currMinCost = dpCUBE[i-1][l][m/25];
									//micro.first = l;
									//micro.second = n;
								//}
								//cout << i << " " << j << " "<<  curr << endl;
								dpCUBE[i][j][k/25] = min(dpCUBE[i][j][k/25], curr);
								//cout << i << " " << j << " " << k<< " " << dpCUBE[i][j][k/25] << endl;
							}
							else { 
								long long curr = dpCUBE[i-1][l][m/25] + transportCost[j][l] + costToStay(j, i, k, n);
								//cout << curr << endl;
								// if(curr < 0) curr = 0;
								//if(dpCUBE[i][j][k/25] >= curr) {
									//cout << costToStay(j, i, n, n) << endl;
									//currMinCost = dpCUBE[i-1][l][m/25];
									//micro.first = l;
									//micro.second = m;
								//}
								dpCUBE[i][j][k/25] = min(dpCUBE[i][j][k/25], curr);
							}
						}
					}
				}
			}
		}
		//path.push_back(make_pair(make_pair(micro.first, micro.second), currMinCost));
	}

//RECONSTRUCTING PATH
	long long minimum = 1e18; //constructing path of largest negative cost (backtracking), minimum denotes current largest negative cost
	pair<int, int> micro;
	for (int j = 0; j < M; j++) {
		for (int k = 0; k <= 100; k+=25) {
			if (dpCUBE[T][j][k / 25] < minimum) { //finding first path node (from day T)
				micro.first = j; // microgrid causing largest negative cost
				micro.second = k; // battery state causing largest negative cost
				minimum = dpCUBE[T][j][k / 25];
			}
		}
	}
	path.push_back(make_pair(micro, minimum));

	for (int i = T - 1; i >= 0; i--) {
		bool found = false;
		int j = micro.first;
		int k = micro.second;
		for(int l = 0; l < M; l++) { //cycling through each microgrid of previous day
			for(int m = 0; m <= 100; m+=25) { //cycling through each of the charge states of previous day
				for(int n = k; n <= 100; n+=25) { // cycling from charge state = end state, until charge state = 100
					if (found) continue;
					long long curr = 0;
					if(n != m) {
						curr = dpCUBE[i][l][m/25] + costToCharge(m, n) + transportCost[j][M] + transportCost[M][l] + costToStay(j, i+1, k, n);
					}
					else { 
						curr = dpCUBE[i][l][m/25] + transportCost[j][l] + costToStay(j, i+1, k, n);
					}
					if (curr == minimum) { //we found microgrid and battery state causing largest negative cost
						minimum = dpCUBE[i][l][m/25];
						found = true;
						micro.first = l;
						micro.second = m;
						demand[make_pair(i+1, j)] = max(demand[make_pair(i+1, j)] - (n-k), 0); //we update the demand at microgrid j at day i+1
					}
				}
			}
		}
		path.push_back(make_pair(micro, minimum));
	}

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
	for (int i = 0; i <= T; i++) {
		for (int j = 0; j < M; j++) {
			demand[make_pair(i, j)] = rand() % 100 + 100;
		}
	}

//INITIALIZE TRANSPORT_COST MATRIX
	vector<vector<long long> > transportCost(M+1, vector<long long>(M+1));

	for(int i = 0; i < M+1; i++) {
		for(int j = i; j < M+1; j++) {
			if(i == j) transportCost[i][j] = 0;
			else {
				long long currRand = 0;
				if(i == M || j == M) currRand = (long long)((rand() % 3)) + 1;
				else currRand = (long long)((rand() % 9)) + 1;
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

//INITIALIZE DPCUBE
	vector<vector<vector<long long> > > dpCUBE;

	for(int i = 0; i <= T; i++) {
		vector<vector<long long> > curr1; 
		for(int j = 0; j < M; j++) {
			vector<long long> curr2;
			for(int k = 0; k <= 100; k+=25) {
				curr2.push_back(1000000);
				//cout << curr2[k/25] << " ";
			}
			//cout << endl;
			curr1.push_back(curr2);
		}
		dpCUBE.push_back(curr1);
	}

	for(int i = 0; i < M; i++) {
		for(int j = 0; j <= 100; j+=25) {
			if(j-100 == 0) dpCUBE[0][i][j/25] = 0; //at day 0, battery can be at 100% at any microgrid at 0 cost
			else dpCUBE[0][i][j/25] = 1000000;
		}
	}

for(int k = 0; k < numBatteries; k++) { //loops through the batteries
	for (int i = 1; i <= T; i++) {
		for (int j = 0; j < M; j++) {
			for (int l = 0; l <= 100; l += 25) {
				dpCUBE[i][j][l / 25] = 1000000;
			}
		}
	}
	for(int i = 0; i < M; i++) {
		for(int j = 0; j <= 100; j+=25) {
			if(j-100 == 0) dpCUBE[0][i][j/25] = 0;
			else dpCUBE[0][i][j/25] = 1000000;
		}
	}

	vector<pair<pair<int, int>, int> > path;
	cheapestPathDP(T, M, path, transportCost, dpCUBE);
	cout << "Minimum cost path for battery " << k+1 << ": " << endl; //outputs min cost path for each battery in standard output
	for(int i = 0; i < path.size(); i++) {
		cout << setw(2) << i+1 << setw(2) << " | microgrid " << setw(2) << path[i].first.first << setw(2) << 
		" | state: " << setw(2) << path[i].first.second << setw(2) << " | cost: " << setw(2) << path[i].second << setw(2) << endl;
	}
}


	return 0;
}










