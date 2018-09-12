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
	int costToStay[100];
	int time;
	int ID;
	Microgrid(int i, int t, int curr) {
		ID = i;
		time = t;
		currState = curr;
		for(int i = 0; i <= currState; i++) {
			costToStay[i] = (currState - i)*200; 
		}
		for(int i = currState + 1; i < 100; i++) costToStay[i] = INT_MAX;
	}
	Microgrid(int charge) {
		ID = -1;
		time = -1;
		currState = charge;
		endState = 100;
		costToStay[currState] = (100-currState)*200;
	}
};

struct Network {
	int numDays;
	int numMicros;
	map<Microgrid, map<Microgrid, int> > edgeCosts;
	int transportCost;
	Network(int t, int m) {
		Microgrid B = new Microgrid(0, 0, 100);
		B.endState = 100;
		numDays = t;
		numMicros = m;
		transportCost = 0;
		map<Microgrid, int> bEdges;
		for(int i = 1; i <= numMicros; i++) {
			Microgrid newMicro = new Microgrid(i, 1, B.endState);
			bEdges.insert(pair<Microgrid, int>(newMicro, transportCost));
		}
		edgeCosts.insert(pair<Microgrid, map<Microgrid, int> >(B, bEdges));
		Microgrid currMicro = B;
		map<Microgrid, int> currMicroEdges;
		for(int i = 2; i <= numDays; i++) {
			for(int j = 1; j <= numMicros; j++) {
				Microgrid newMicro = new Microgrid(j, i, currMicro.endState);
				currMicroEdges.insert(pair<Microgrid, int>(newMicro, transportCost));
			}
			edgeCosts.insert(pair<Microgrid, map<Microgrid, int> >(currMicro, currMicroEdges));
		}
	}
};











