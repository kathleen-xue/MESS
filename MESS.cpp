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

struct Microgrid {
	int currState;
	int endState;
	int time;
	int ID;
	Microgrid(int i, int t, int curr) {
		ID = i;
		time = t;
		currState = curr;
	}
};

struct Network {
	int numDays;
	int numMicros;
	vector<vector<int> > edgeCosts;
	int transportCost;
	Microgrid* B;
	Microgrid* D;
	Network(int t, int m) { //t=#days, m=#microgrids
		B = new Microgrid(0,0,100);
		B->endState = 100;
		numDays = t;
		numMicros = m;
		transportCost = 50;
		for(int i = 0; i <= numMicros; i++) {
			vector<int> currEdges;
			for(int j = 0; j <= numMicros; j++) {
				if(j != i) currEdges.push_back(transportCost);
				else currEdges.push_back(0);
			}
			edgeCosts.push_back(currEdges);
		}
	}
};

vector<vector<vector<int> > > cheapestPath(int T, int M) {

	int transportCost = 50;

	vector<vector<vector<int> > > dpCUBE;
	vector<vector<int> > first;

	for(int i = 0; i < M; i++) {
		vector<int> currMicro;
		for(int j = 0; j <= 100; j += 25) {
			currMicro.push_back(costToStay(100,j));
		}
		first.push_back(currMicro);
	}

	dpCUBE.push_back(first);

	vector<vector<int> > curr;
	for(int i = 1; i < T; i++) {
		for(int j = 1; j < M; j++) {
			vector<int> currMicro;
			int currMin = INT_MAX;
			for(int k = 0; k <= 100; k += 25) {
				for(int l = 0; l < M; l++) {
					for(int m = 0; m <= 100; m += 25) {
						int bestCharge = 0;
						for(int n = m; n <= 100; n += 25) {
							//if((m == n) && dpCUBE[i-1][l][m] + )
							currMin = min(currMin, dpCUBE[i-1][l][m] + costToCharge(m, n) + 2*transportCost);
						}
					}
				}
			}
		}
	}
}

int costToCharge(int start, int end) {
	return (end - start)*200;
}

int costToStay(int start, int end) {
	return  -1*(start - end)*200;
}


int main () {
	int T = 0;
	int M = 0;
	cout << "Enter # days: " << endl;
	cin >> T;
	cout << "Enter # microgrids: " << endl;
	cin >> M;

	for(int i = 0; i <= T; i++) {

	}
	return 0;
}










