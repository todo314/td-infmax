#include "./models.h"

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
#include <assert.h>

#include "./mt19937ar.h"
#include "./tools.h"

using namespace std;

double sample_delay(string &delay, double a, double b) {
	if (delay == "weibull") {
		double x = genrand_res53();
		return b * pow(-log1p(-x), 1.0 / a);
	}
	if (delay == "weibull-ic" || delay == "weibull-lt") {
		double z = 1.0 / a;
		double pi = acos(-1);
		double e = exp(1);
		return b * sqrt(2 * pi * z) * pow(z / e, z);
	}
	assert(false);
	return -1;
}

double prob_TVIC(string &prob, double c, int ideg, double t_u, double d_uv) {
	if (prob == "exponential" || prob == "exponential-ic") {
		return 1.0 / ideg * exp(-c * (t_u + d_uv));
	}
	if (prob == "deadline") {
		return t_u + d_uv < c ? 1 : 0;
	}
	if (prob == "inverse" || prob == "inverse-ic") {
		return 1.0 / (ideg * c * (t_u + d_uv));
	}
	assert(false);
	return false;
}

double sample_latest_TVIC(string &prob, double c, int ideg, double d_uv) {
	const double x = genrand_res53();
	if (prob == "exponential") {
		return -1.0 / c * log(x * ideg) - d_uv;
	}
	if (prob == "deadline") {
		return c - d_uv;
	}
	if (prob == "inverse") {
		return 1.0 / (x * ideg * c) - d_uv;
	}
	if (prob == "exponential-ic" || prob == "inverse-ic") {
		if (x < prob_TVIC(prob, c, ideg, 0, d_uv)) {
			return numeric_limits<double>::infinity();
		} else {
			return -numeric_limits<double>::infinity();
		}
	}
	assert(false);
	return -1;
}

double weight_TVLT(string &prob, double c, int ideg, double t_u, double d_uv) {
	if (prob == "exponential" || prob == "exponential-lt") {
		return 1.0 / ideg * exp(-c * (t_u + d_uv));
	}
	if (prob == "deadline") {
		return t_u + d_uv < c ? 1.0 / ideg : 0;
	}
	if (prob == "inverse" || prob == "inverse-lt") {
		return 1.0 / (ideg * c * (t_u + d_uv + 1));
	}
	assert(false);
	return -1;
}

double sample_latest_TVLT(string &prob, double c, int ideg, double d_uv) {
	const double x = genrand_res53();
	if (prob == "exponential") {
		return -1.0 / c * log(ideg * weight_TVLT(prob, c, ideg, 0, d_uv) * x)
				- d_uv;
	}
	if (prob == "deadline") {
		if (weight_TVLT(prob, c, ideg, 0, d_uv) * x > 1.0 / ideg) {
			return -1; // IMPOSSIBLE
		} else {
			return c - d_uv;
		}
	}
	if (prob == "inverse") {
		return 1.0 / (ideg * c * weight_TVLT(prob, c, ideg, 0, d_uv) * x) - d_uv
				- 1;
	}
	if (prob == "exponential-lt" || prob == "inverse-lt") {
		if (x < weight_TVLT(prob, c, ideg, 0, d_uv)) {
			return numeric_limits<double>::infinity();
		} else {
			return -numeric_limits<double>::infinity();
		}
	}
	assert(false);
	return -1;
}

int MonteCarlo_TVIC(int V, vector<vector<edge> > &es, vector<int> &ideg,
		string &delay, string &prob, vector<int> &S) {
	priority_queue<pair<double, int>, vector<pair<double, int> >,
			greater<pair<double, int> > > Q;
	map<int, double> tau;
	for (int s : S) {
		if (tau.count(s) == 0) {
			tau[s] = 0;
			Q.push(make_pair(tau[s], s));
		}
	}

	vector<int> FR;
	for (; !Q.empty();) {
		int u = Q.top().second;
		double t = Q.top().first;
		Q.pop();
		if (t > tau[u]) {
			continue;
		}
		FR.push_back(u);
		for (auto e : es[u]) {
			double d_uv = sample_delay(delay, e.a, e.b);
			if ((tau.count(e.v) == 0 || tau[u] + d_uv < tau[e.v])
					&& genrand_res53()
							< prob_TVIC(prob, e.c, ideg[e.v], tau[u], d_uv)) {
				tau[e.v] = tau[u] + d_uv;
				Q.push(make_pair(tau[e.v], e.v));
			}
		}
	}
	return FR.size();
}

int MonteCarlo_TVLT(int V, vector<vector<edge> > &es, vector<int> &ideg,
		string &delay, string &prob, vector<int> &S) {
	priority_queue<pair<double, int>, vector<pair<double, int> >,
			greater<pair<double, int> > > Q;
	map<int, double> tau;
	for (int s : S) {
		if (tau.count(s) == 0) {
			tau[s] = 0;
			Q.push(make_pair(tau[s], s));
		}
	}
	map<int, double> ts, I;

	vector<int> FR;
	for (; !Q.empty();) {
		int u = Q.top().second;
		double t = Q.top().first;
		Q.pop();
		if (t > tau[u]) {
			continue;
		}
		FR.push_back(u);
		for (auto e : es[u]) {
			double d_uv = sample_delay(delay, e.a, e.b);
			I[e.v] += weight_TVLT(prob, e.c, ideg[e.v], tau[u], d_uv);
			if (ts.count(e.v) == 0) {
				ts[e.v] = genrand_res53();
			}

			if ((tau.count(e.v) == 0 || tau[u] + d_uv < tau[e.v])
					&& I[e.v] >= ts[e.v]) {
				tau[e.v] = tau[u] + d_uv;
				Q.push(make_pair(tau[e.v], e.v));
			}
		}
	}
	return FR.size();
}
