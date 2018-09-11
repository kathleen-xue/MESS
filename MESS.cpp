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

struct Node {
	int microgrid;
	int time;
	int chargeState;
	bool isUsed;
	vector<Node> linkedTo;

	Node(int m_i, int t, int j, bool u) { //m_i = Microgrid i, t = time, j = charge state, u = used or not
		microgrid = m_i;
		time = t;
		chargeState = j;
		isUsed = u;
	}
	Node(int k) { //k = charging level
		chargeState = k;
	}
};

struct Microgrid {
	Node n[4][2];
	vector<vector<int> > edges;
	int ID;
	int time;
	Microgrid(int i, int t) {
		ID = i;
		time = t;
		for(int i = 0; i < 4; i++) {
			vector<int> currEdges;
			for(int j = 0; j < 2; j++) {
				n[i][j] = new Node(ID, t, 1 - (i*0.25), j == 1);
			}
			for(int k = i; k < 4; k++) {
					n[i][0].linkedTo.push_back(n[k][1]);
					currEdges.push_back(0);
			}
			edges.push_back(currEdges);
		}
	}
};

struct Network {
	Node source;
	int numTimes; // number of days
	int numMicros; // number of Microgrids
	Network(int t, int m) {
		numTimes = t;
		numMicros = m;
		source = new Node(-1, 0, 100, false);
		for(int i = 0; i < numTimes; i++) {
			
		}
	}
};











