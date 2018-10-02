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
#include <iomanip>
using namespace std;

int costToCharge(int start, int end) {
	if(end <= start) return INT_MAX;
	return (end - start)*20;
}

int costToStay(int M, int T, int start, int end) {
	if(start <= end) return INT_MAX;
	return -1*(start - end)*200;
}

vector<vector<vector<long long> > > cheapestPathDP(int T, int M, vector<pair<pair<int, int>, int> >& path) {
//This function finds the maximum reduction of cost for a battery to
//end up at any microgrid at any charge state on day T. If we iterate 
//this function through x batteries, the final dpCUBE will output 
//maximum reduction of cost from all the batteries for each
//microgrid after T days (this will be seen on the last level of the cube).
	vector<vector<int> > transportCost;
	for(int i = 0; i < M+1; i++) {
		vector<int> curr;
		for(int j = 0; j < M+1; j++) {
			if(i == j) curr.push_back(0);
			else {
				curr.push_back(abs(i-j)*3);
			}
		}
		transportCost.push_back(curr);
	}

	vector<vector<vector<long long> > > dpCUBE;

	for(int i = 0; i < T; i++) {
		vector<vector<long long> > curr1; 
		for(int j = 0; j < M; j++) {
			vector<long long> curr2;
			for(int k = 0; k <= 100; k+=25) {
				curr2.push_back(INT_MAX);
				//cout << curr2[k/25] << " ";
			}
			//cout << endl;
			curr1.push_back(curr2);
		}
		dpCUBE.push_back(curr1);
	}


	for(int i = 0; i < M; i++) {
		for(int j = 0; j <= 100; j+=25) {
			dpCUBE[0][i][j/25] = costToStay(i, 0, 100, j);
		}
	}
	
	for(int i = 1; i < T; i++) { //cycling through each day
		int currMinCost = INT_MAX;
		pair<int, int> micro;
		for(int j = 1; j < M; j++) { //cycling through the microgrids of current day
			for(int k = 0; k <= 100; k+=25) { //cycling through each of the charge states of each microgrid of current day
				for(int l = 0; l < M; l++) { //cycling through each microgrid of previous day
					for(int m = 0; m <= 100; m+=25) { //cycling through each of the charge states of previous day
						for(int n = m; n <= 100; n+=25) { //cycling through each of the *possible* charges of current day
							if(n != m) {
								if(dpCUBE[i][j][k/25] > dpCUBE[i-1][l][m/25] + costToCharge(m, n) + transportCost[j][M] + transportCost[M][l] + costToStay(j, i, n, k)) {
									currMinCost = dpCUBE[i-1][l][m/25];
									micro.first = l;
									micro.second = m;
								}
								dpCUBE[i][j][k/25] = min(dpCUBE[i][j][k/25], dpCUBE[i-1][l][m/25] + costToCharge(m, n) + 
			                						 transportCost[j][M] + transportCost[M][l] + costToStay(j, i, n, k));
							}
							else { 
								if(dpCUBE[i][j][k/25] > dpCUBE[i-1][l][m/25] + transportCost[j][l]) {
									currMinCost = dpCUBE[i-1][l][m/25];
									micro.first = l;
									micro.second = m;
								}
								dpCUBE[i][j][k/25] = min(dpCUBE[i][j][k/25], (dpCUBE[i-1][l][m/25] + transportCost[j][l]));
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
	//testing to see if segfault; will make other tests for functionality
	int T = 0;
	int M = 0;
	cout << "Enter # days: " << endl;
	cin >> T;
	cout << "Enter # microgrids: " << endl;
	cin >> M;

	vector<pair<pair<int, int>, int> > path;
	vector<vector<vector<long long> > > dpCUBE = cheapestPathDP(T, M, path);
	for(int i = 0; i < path.size(); i++) {
		cout << i+1 << " microgrid: " << path[i].first.first << " state: " << path[i].first.second << " cost: " << path[i].second << endl;
	}

	return 0;
}










