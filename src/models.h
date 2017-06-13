#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <stack>
#include <bitset>
#include <algorithm>
#include <map>
#include <set>
#include <ctime>
#include <cmath>

#include "./mt19937ar.h"
#include "./tools.h"

using namespace std;

typedef struct {
	int u, v;
	double a, b, c;
} edge;

double sample_delay(string &delay, double a, double b);

double prob_TVIC(string &prob, double c, int ideg, double t_u, double d_uv);
double sample_latest_TVIC(string &prob, double c, int ideg, double d_uv);

double weight_TVLT(string &prob, double c, int ideg, double t_u, double d_uv);
double sample_latest_TVLT(string &prob, double c, int ideg, double d_uv);

int MonteCarlo_TVIC(int V, vector<vector<edge> > &es, vector<int> &ideg,
		string &delay, string &prob, vector<int> &S);
int MonteCarlo_TVLT(int V, vector<vector<edge> > &es, vector<int> &ideg,
		string &delay, string &prob, vector<int> &S);
