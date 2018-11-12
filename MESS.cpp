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

//TO DO:
//benefit per microgrid per day (use same function, just total charge of all batteries)
//total benefit - total cost
//make it generate different random values each time
//plot MESS vs Greedy for 1 random setting
//plot MESS vs Greedy for 10 random settings //BY FRIDAY
//2 plots each for [10 days,5 microgrids], [10,10], for batteries from 0-50 in increments of 2
//also do square root benefit, fourth root benefit
//write a script so it's easy to see the results for changed parameters

//DONE:
//edit path visualisation (add start state, whether battery arrived at charging station)
//each day, generate a new set of demands for each microgrid
//demand of microgrid i on day2 can only be affected by batteries on microgrid i on day2
//change demands to be not affected by previous day's demand
//compare to simpler heuristic (simple greedy algorithm)

struct path {
	vector<pair<int, int> > states; //battery states: path.states[i].first = start state, path.states[i].second = end state
	vector<int> micros; //microgrids
	vector<int> costs;
	vector<int> negs;
};

vector<vector<int> > createDemands(int M, int T) {
	vector<vector<int> > demands;
	for(int i = 0; i < T; i++) {
		vector<int> demDay_i;
		for(int j = 0; j < M; j++) {
			demDay_i.push_back((rand() % 300) + 0);
		}
		demands.push_back(demDay_i);
	}
	return demands;
}


int costToCharge(int start, int end) {
	return (end - start);
}


int costToStay(int M, int T, int end, int start, vector<vector<int> >& demands) { //at each microgrid compute total number of batteries that arrived and the amount of benefit they created total
	//return -1 * sqrt(start - end) * 10;
	int d = demands[T-1][M];
	double diff = d - abs(start - end)*0.5;
	//cout << "costToStay: " << -1 * sqrt(d) * ((d - max(diff, 0.0))*3) << endl;
	//cout << diff << endl;
	return -1 * pow(d, 0.25) * abs((d - max(diff, 0.0))*3); //created a map of microgrid demands because previous costToStay function assumed 
															 //demands didn't exist and caused batteries to stay at one microgrid throughout 
															 //duration of days
}


path greedy(int T, int M, vector<vector<long long> > transportCost, vector<vector<int> >& demands) { //simple greedy algo
	path P;
	for(int i = 0; i <= T; i++) { //initialize path
		P.states.push_back(make_pair(0,0));
		P.micros.push_back(0);
		P.costs.push_back(0);
		P.negs.push_back(0);
	}

	P.states[0] = make_pair(0, 100); //at day 0, can charge from 0 to 100 for free

	int maxDemand = 0;
	for(int i = 0; i < M; i++) {
		if(demands[0][i] > maxDemand) {
			maxDemand = demands[0][i];
			P.micros[0] = i;
		}
	}

	for(int i = 1; i <= T; i++) { //meat of the greedy
		int minCost = INT_MAX;
		for(int j = 0; j < M; j++) {
			for(int k = (P.states[i-1].second)/25; k <= 4; k++) {
				for(int l = k; l >= 0; l--) {
					int curr = costToCharge((P.states[i-1].second), k*25) + transportCost[P.micros[i-1]][M] + transportCost[M][j] + costToStay(j, i, k*25, l*25, demands) + P.costs[i-1];
					int currneg = costToCharge((P.states[i-1].second), k*25) + transportCost[P.micros[i-1]][M] + transportCost[M][j];
					if(curr < minCost) {
						minCost = curr;
						P.micros[i] = j;
						P.costs[i] = curr;
						P.states[i] = make_pair(k*25, l*25);
						P.negs[i] = currneg;
						//cout << "microgrid: " << P.micros[i] << " cost: " << P.costs[i] << " start %: " << P.states[i].first << " end %: " << P.states[i].second << endl;
						//cout << "minCost: " << minCost << endl;
					}
				}
			}
		}
	}

	for(int i = 1; i <= T; i++) {
		int margDecrease = demands[i-1][P.micros[i]] - 10*sqrt(demands[i-1][P.micros[i]] - abs(P.states[i].first - P.states[i].second));
		//cout << "demands[" << i-1 << "][" << P.micros[i] << "] " << demands[i-1][P.micros[i]] << " decreased by " << 3*sqrt(demands[i-1][P.micros[i]] - abs(P.states[i].first - P.states[i].second)) << endl;
		demands[i-1][P.micros[i]] = max(margDecrease, 0);
	}

	return P;
}


void cheapestPathDP(int T, int M, path& P, vector<vector<long long> > transportCost, vector<vector<vector<long long> > >& dpCUBE, vector<vector<int> >& demands) {
//This function finds the maximum reduction of cost for a battery to
//end up at any microgrid at any charge state on day T. If we iterate 
//this function through x batteries, the final dpCUBE will output 
//maximum reduction of cost from all the batteries for each
//microgrid after T days (this will be seen on the last level of the cube).

//FINDING MINIMUM COSTS	
	int startState = 0;
	int endState = 0;
	for(int i = 1; i <= T; i++) { //cycling through each day
		for(int j = 0; j < M; j++) { //cycling through the microgrids of current day
			for(int l = 0; l < M; l++) { //cycling through each microgrid of previous day
				for(int m = 0; m <= 4; m++) { //cycling through each of the charge states of previous day
					for(int n = m; n <= 4; n++) { //cycling through each of the *possible* charges of current day
						for(int k = n; k >= 0; k--) { //cycling through each of the charge states of each microgrid of current day
							if(n != m) {
								long long curr = dpCUBE[i-1][l][m] + costToCharge(m*25, n*25) + transportCost[j][M] + transportCost[M][l] + costToStay(j, i, k*25, n*25, demands);
								dpCUBE[i][j][k] = min(dpCUBE[i][j][k], curr);
								if(dpCUBE[i][j][k] == curr) {
									startState = n*25;
									endState = k*25;
								}
							}
							else { 
								long long curr = dpCUBE[i-1][l][m] + transportCost[j][l] + costToStay(j, i, k*25, n*25, demands);
								dpCUBE[i][j][k] = min(dpCUBE[i][j][k], curr);
								if(dpCUBE[i][j][k] == curr) {
									startState = n*25;
									endState = k*25;
								}
							}
						}
					}
				}
			}
			//int margDecrease = demands[i-1][j] - sqrt(demands[i-1][j] - abs(startState - endState));
			//cout << "margDecrease: " << margDecrease << endl;
			//demands[i-1][j] = max(margDecrease, 0);
		}
	}

//RECONSTRUCTING PATH
	long long minimum = 1e18; //constructing path of largest negative cost (backtracking), minimum denotes current largest negative cost

	for(int i = 0; i <= T; i++) {
		P.states.push_back(make_pair(0,0));
		P.micros.push_back(0);
		P.costs.push_back(0);
		P.negs.push_back(0);
	}

	//P.states[0].first = 0;

	for (int j = 0; j < M; j++) {
		for (int k = 0; k <= 4; k++) {
			if (dpCUBE[T][j][k] < minimum) { //finding first path node (from day T)
				P.micros[0] = j; // microgrid causing largest negative cost
				P.states[0].first = k*25; // battery state causing largest negative cost
				minimum = dpCUBE[T][j][k];
			}
		}
	}
	P.costs[0] = minimum;
	//path.push_back(make_pair(micro, minimum));

	for (int i = T - 1; i >= 0; i--) {
		int i_ = T - i;
		bool found = false;
		int j = P.micros[i_ - 1];
		int k = P.states[i_ - 1].first;
		for(int l = 0; l < M; l++) { //cycling through each microgrid of previous day
			for(int m = 0; m <= 4; m++) { //cycling through each of the charge states of previous day
				for(int n = k; n <= 4; n++) { // cycling from charge state = end state, until charge state = 100
					if (found) continue;
					int currneg = 0;
					long long curr = 0;
					if(n != m) {
						curr = dpCUBE[i][l][m] + costToCharge(m*25, n*25) + transportCost[j][M] + transportCost[M][l] + costToStay(j, i+1, k, n*25, demands);
						currneg = costToCharge(m*25, n*25) + transportCost[j][M] + transportCost[M][l];
					}
					else { 
						curr = dpCUBE[i][l][m] + transportCost[j][l] + costToStay(j, i+1, k, n*25, demands);
						currneg = transportCost[j][l];
					}
					if (curr == minimum) { //we found microgrid and battery state causing largest negative cost
						minimum = dpCUBE[i][l][m];
						found = true;
						P.micros[i_] = l;
						P.states[i_].second = m*25;
						P.states[i_ - 1].first = n*25;
						P.costs[i_] = minimum;
						P.negs[i_] = currneg;
						//micro.first = l;
						//micro.second = m;
						int margDecrease = demands[i][j] - 10*sqrt(demands[i][j] - abs(P.states[i_-1].first - P.states[i_-1].second));
						//cout << sqrt(demands[i][j] - abs(P.states[i_-1].first - P.states[i_-1].second)) <<" " << margDecrease << endl;
						demands[i][j] = max(margDecrease, 0);
						//demand[make_pair(i+1, j)] = max(demand[make_pair(i+1, j)] - (n-k), 0); //we update the demand at microgrid j at day i+1
					}
				}
			}
		}
		//path.push_back(make_pair(micro, minimum));
	}

}


int main () {
	int T = 0;
	int M = 0;
	int numBatteries = 0;
	int lowerBound = 0;
	int upperBound = 0;
	int increments = 0;
	cout << "MENU: " << endl;
	cout << "1 = Clear output file" << endl;
	cout << "2 = Add input" << endl;
	cout << "3 = Add bulk input (multiple numbers of batteries)" << endl;
	int userInput;
	cin >> userInput;
	if(userInput == 1) {
		ofstream f("netBenefits.txt", ios::out | ios::trunc);
		f.close();
		return 0;
	}
	else if(userInput == 2) {
		cout << "Enter # days: " << endl;
		cin >> T;
		cout << "Enter # microgrids: " << endl;
		cin >> M;
		cout << "Enter # batteries: " << endl;
		cin >> numBatteries;
		/*for (int i = 0; i <= T; i++) {
			for (int j = 0; j < M; j++) {
				demand[make_pair(i, j)] = rand() % 200;
			}
		}*/
	}
	else {
		cout << "Enter lower bound of # batteries: " << endl;
		cin >> lowerBound;
		cout << "Enter upper bound of # batteries: " << endl;
		cin >> upperBound;
		cout << "Enter increment amount: " << endl;
		cin >> increments;
		cout << "Enter # days: " << endl;
		cin >> T;
		cout << "Enter # microgrids: " << endl;
		cin >> M;
	}

	if(userInput == 2) {
		upperBound = 0;
		increments = 1;
	}
	

//INITIALIZE DEMAND MATRIX
for(int z = lowerBound; z <= upperBound; z += increments) {
	if(userInput == 3) numBatteries = z;
	vector<vector<int> > demands = createDemands(M, T); //create demands matrix
	vector<vector<int> > demandsGreed = demands;
	vector<vector<int> > demandsMicro = demands;
	cout << "Demands Matrix: (starting demand of microgrid j at day i)" << endl;
	cout << "mGrid ";
	for(int i = 0; i < M; i++) {
		cout << "\t" << i;
	}
	cout << endl;
	for(int i = 0; i < T; i++) {
		for(int j = 0; j < M; j++) {
			cout << "\t" << demands[i][j];
		}
		cout << endl;
	}

	//INITIALIZE TRANSPORT_COST MATRIX
		vector<vector<long long> > transportCost(M+1, vector<long long>(M+1));

		for(int i = 0; i < M+1; i++) {
			for(int j = i; j < M+1; j++) {
				if(i == j) transportCost[i][j] = 0;
				else {
					long long currRand = 0;
					//if(i == M || j == M) currRand = (long long)((rand() % 3)) + 1;
					currRand = (long long)((rand() % 600)); //experiment more with transportation cost
					transportCost[i][j] = currRand;
					transportCost[j][i] = currRand;
				}
			}
		}

		cout << "Transportation Costs: (between microgrid i and microgrid j)" << endl; //printing out transport costs
		for(int i = 0; i < M+1; i++) {
			for(int j = 0; j < M+1; j++) {
				cout << "\t" << transportCost[i][j];
			}
			cout << endl;
		}
		cout << endl;

	//INITIALIZE DP_CUBE
		vector<vector<vector<long long> > > dpCUBE;

		for(int i = 0; i <= T; i++) {
			vector<vector<long long> > curr1; 
			for(int j = 0; j < M; j++) {
				vector<long long> curr2;
				for(int k = 0; k <= 100; k+=25) {
					curr2.push_back(10000000);
					//cout << curr2[k/25] << " ";
				}
				//cout << endl;
				curr1.push_back(curr2);
			}
			dpCUBE.push_back(curr1);
		}

	//RUN ALGORITHM
	vector<vector<pair<int, int> > > microgridsDP;
	for(int y = 0; y <= T; y++) {
		vector<pair<int, int> > curr;
		for(int x = 0; x < M; x++) {
			curr.push_back(make_pair(0,0));
		}
		microgridsDP.push_back(curr);
	}

	vector<vector<pair<int, int> > > microgridsGreed;
	for(int y = 0; y <= T; y++) {
		vector<pair<int, int> > curr;
		for(int x = 0; x < M; x++) {
			curr.push_back(make_pair(0,0));
		}
		microgridsGreed.push_back(curr);
	}

	long long DPtotalCost = 0;
	for(int k = 0; k < numBatteries; k++) { //loops through the batteries
		for (int i = 1; i <= T; i++) { //reinitialize dpCUBE
			for (int j = 0; j < M; j++) {
				for (int l = 0; l <= 4; l++) {
					dpCUBE[i][j][l] = 10000000;
				}
			}
		}
		for(int i = 0; i < M; i++) {
			for(int j = 0; j <= 4; j++) {
				if(j == 4) dpCUBE[0][i][j] = 0;
				else dpCUBE[0][i][j] = 10000000; //at day 0, battery can be at 100% at any microgrid at 0 cost
			}
		}
		
		path P;
		cheapestPathDP(T, M, P, transportCost, dpCUBE, demands); //actually run algorithm

		//PRINT PATH
		long long totalnegs = 0;
		for(int i = 0; i < P.negs.size(); i++) {
			totalnegs += P.negs[i];
		}
		cout << "MESS ALGORITHM minimum cost path for battery " << k+1 << ": " << endl; //outputs min cost path for each battery in standard output
		for(int i = T; i >= 0; i--) {
			microgridsDP[i][P.micros[i]].first += P.states[i].first;
			//cout << "DP start charge: " << microgridsDP[P.micros.size() - i - 1][P.micros[i]].first << endl;
			microgridsDP[i][P.micros[i]].second += P.states[i].second;
			cout << "day " << T - i << '\t' << "| microgrid " << P.micros[i] << '\t' << //state is battery state after usage on day i
			"| start %: " << P.states[i].first << '\t' << "| end %: " << P.states[i].second << '\t' << 
			"| cost: " << P.costs[i] << '\t'; //cost is the maximum negative cost accrued for smart grid per battery
			if(i < T && P.states[i].first > P.states[i+1].second) cout << "| charged ";
			cout << endl;
		}
			//DPtotalCost += P.costs[0];
			DPtotalCost += totalnegs;

	}

		long long totalCost = 0;
		for(int k = 0; k < numBatteries; k++) {
			path greed = greedy(T, M, transportCost, demandsGreed);

			long long totalnegs_g = 0;
			for(int i = 0; i < greed.negs.size(); i++) {
				totalnegs_g += greed.negs[i];
			}

			cout << "SIMPLE GREEDY ALGORITHM minimum cost path for battery " << k+1 << ": " << endl;
			for(int l = 0; l <= T; l++) {
				microgridsGreed[l][greed.micros[l]].first += greed.states[l].first;
				microgridsGreed[l][greed.micros[l]].second += greed.states[l].second;
				cout << "day " << l << '\t' << "| microgrid " << greed.micros[l] << '\t' << "| start %: " << greed.states[l].first <<
				'\t' << "| end %: " << greed.states[l].second << '\t' << "| cost: " << greed.costs[l] << '\t';
				if(l > 0 && greed.states[l].first > greed.states[l-1].second) cout << "| charged ";
				cout << endl;
			}
			//totalCost += greed.costs[T];
			totalCost += totalnegs_g;
		}

		long long MESStotalSum = 0;
		long long GREEDYtotalSum = 0;

		cout << "MESS BENEFIT PER MICROGRID: " << endl;
		for(int t = 0; t < T; t++) {
			cout << "Day " << t << endl;
			for(int x = 0; x < M; x++) {
				int ans = costToStay(x, t+1, microgridsDP[t][x].second, microgridsDP[t][x].first, demandsMicro);
				//cout << ans << endl;
				cout << "Microgrid ";
				cout << x << ": " << ans << endl;
				MESStotalSum += ans;
			}
		}

		cout << endl;

		cout << "GREEDY BENEFIT PER MICROGRID: " << endl;
		for(int t = 0; t < T; t++) {
			cout << "Day " << t << endl;
			for(int x = 0; x < M; x++) {
				int ans = costToStay(x, t+1, microgridsGreed[t][x].second, microgridsGreed[t][x].first, demandsMicro);
				cout << "Microgrid " << x << ": " << ans << endl;
				GREEDYtotalSum += ans;
			}
		}
		cout << endl;

		cout << "MESS TOTAL NET BENEFIT: " << MESStotalSum + DPtotalCost << endl;
		cout << "GREEDY TOTAL NET BENEFIT: " << GREEDYtotalSum + totalCost << endl;
		cout << endl;
		cout << "MESS ALGORITHM TOTAL COST: \t" << DPtotalCost << endl;
		cout << "GREEDY ALGORITHM TOTAL COST: \t" << totalCost << endl;
		
		ofstream myfile("netBenefits.txt", ios_base::app | ios_base::out);
		myfile << (MESStotalSum + DPtotalCost) << " " << (GREEDYtotalSum + totalCost) << "\n";
		myfile.close();

	}
	
	return 0;
}










