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

/*struct Microgrid {
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
};*/

int costToCharge(int start, int end) {
	return (end - start)*200;
}

int costToStay(int M, int T, int start, int end) {
	return  -1*(start - end)*200;
}

vector<vector<vector<int> > > cheapestPath(int T, int M) {

	vector<vector<int> > transportCost;
	for(int i = 0; i < M+1; i++) {
		vector<int> curr;
		for(int j = 0; j < M+1; j++) {
			if(i == j) curr.push_back(0);
			else {
				curr.push_back(abs(i-j)*30);
			}
		}
		transportCost.push_back(curr);
	}

	vector<vector<vector<int> > > dpCUBE;

	for(int i = 0; i < T; i++) {
		vector<vector<int> > curr1; 
		for(int j = 0; j < M; j++) {
			vector<int> curr2;
			for(int k = 0; k < 100; k += 25) {
				curr2.push_back(1000000);
				//cout << curr2[k/25] << " ";
			}
			//cout << endl;
			curr1.push_back(curr2);
		}
		dpCUBE.push_back(curr1);
	}

cout << "hi" << endl;

	for(int i = 0; i < M; i++) {
		for(int j = 0; j <= 100; j += 25) {
			dpCUBE[0][i][j] = costToStay(i, 0, 100, j);
		}
	}
	
	for(int i = 1; i < T; i++) {
	  for(int j = 1; j < M; j++) {
	    for(int k = 0; k <= 100; k += 25) {
	      for(int l = 0; l < M; l++) {
	        for(int m = 0; m <= 100; m += 25) {
	          for(int n = m; n <= 100; n += 25) {
	             dpCUBE[i][j][k] = min(dpCUBE[i][j][k], dpCUBE[i-1][l][m] + costToCharge(m, n) + 
	                transportCost[j][M] + transportCost[M][l] + costToStay(j, i, n, k));
	             if(n == m) {
	                dpCUBE[i][j][k] = min(dpCUBE[i][j][k], (dpCUBE[i-1][l][m] + transportCost[j][l]));
                   }
                }
              }
            }
          }
        }
      }

	return dpCUBE;
}

int main () {
	int T = 0;
	int M = 0;
	cout << "Enter # days: " << endl;
	cin >> T;
	cout << "Enter # microgrids: " << endl;
	cin >> M;

	vector<vector<vector<int> > > dpCUBE = cheapestPath(T, M);

	return 0;
}










