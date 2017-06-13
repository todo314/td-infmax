#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <limits>

#include "./models.h"
#include "./mt19937ar.h"
#include "./tools.h"

using namespace std;

vector<int> imm(int V, vector<vector<edge> > &rs, vector<int> &ideg,
		string model, string delay, string prob, int k, double eps, double ell);
int greedy(int V, vector<vector<int> > &h2v, vector<vector<int> > &v2h, int k,
		vector<int> &S);

double log_fact(int n) {
	double val = 1;
	for (int i = 1; i <= n; i++) {
		val += log(i);
	}
	return val;
}

// {n \choose k} = n!/k!(n-k)!
double log_nCk(int n, int k) {
	return log_fact(n) - log_fact(k) - log_fact(n - k);
}

vector<int> gen_RR_TVIC(int V, vector<vector<edge> > &rs, vector<int> &ideg,
		string &delay, string &prob) {

	int z = genrand_int(V);
	map<int, double> tau;
	tau[z] = numeric_limits<double>::infinity();
	priority_queue<pair<double, int> > Q;
	Q.push(make_pair(tau[z], z));

	vector<int> RR;
	for (; !Q.empty();) {
		int v = Q.top().second;
		double t = Q.top().first;
		Q.pop();
		if (t < tau[v]) {
			continue;
		}
		RR.push_back(v);

		for (auto e : rs[v]) {
			double d_uv = sample_delay(delay, e.a, e.b);
			double t_uv = sample_latest_TVIC(prob, e.c, ideg[v], d_uv);
			double t_u = min(tau[v] - d_uv, t_uv);
			if ((tau.count(e.u) == 0 || t_u > tau[e.u]) && t_u >= 0) {
				tau[e.u] = t_u;
				Q.push(make_pair(tau[e.u], e.u));
			}
		}
	}
	return RR;
}

vector<int> gen_RR_TVLT(int V, vector<vector<edge> > &rs, vector<int> &ideg,
		string &delay, string &prob) {

	int z = genrand_int(V);
	map<int, double> tau;
	tau[z] = numeric_limits<double>::infinity();
	priority_queue<pair<double, int> > Q;
	Q.push(make_pair(tau[z], z));

	vector<int> RR;
	for (; !Q.empty();) {
		const int v = Q.top().second;
		double t = Q.top().first;
		Q.pop();
		if (t < tau[v]) {
			continue;
		}
		RR.push_back(v);

		double roulette = genrand_res53();
		double sum = 0;
		double d_uv;
		int hit = -1;
		for (int i = 0; i < rs[v].size(); i++) {
			edge e = rs[v][i];
			d_uv = sample_delay(delay, e.a, e.b);
			sum += weight_TVLT(prob, e.c, ideg[v], 0, d_uv);
			if (sum >= roulette) {
				hit = i;
				break;
			}
		}
		if (hit >= 0) {
			edge e = rs[v][hit];
			double t_uv = sample_latest_TVLT(prob, e.c, ideg[v], d_uv);
			double t_u = min(tau[v] - d_uv, t_uv);
			if ((tau.count(e.u) == 0 || t_u > tau[e.u]) && t_u >= 0) {
				tau[e.u] = t_u;
				Q.push(make_pair(tau[e.u], e.u));
			}
		} else {
		}
	}
	return RR;
}

// @return deg
int greedy(int V, vector<vector<int> > &h2v, vector<vector<int> > &v2h, int k,
		vector<int> &S) {
	int H = h2v.size();
	vector<bool> dead(H, false);
	vector<int> deg(V);
	for (int v = 0; v < V; v++) {
		deg[v] = v2h[v].size();
	}

	int degS = 0;
	int last = 0;
	for (; S.size() < k;) {
		int s = 0;
		for (int v = 0; v < V; v++) {
			if (deg[v] > deg[s]) {
				s = v;
			}
		}

		degS += deg[s];
		last = degS;

		S.push_back(s);

		for (int e : v2h[s]) {
			if (!dead[e]) {
				dead[e] = true;
				for (int v : h2v[e]) {
					deg[v]--;
				}
			}
		}
	}
	return degS;
}

vector<int> imm(int V, vector<vector<edge> > &rs, vector<int> &ideg,
		string model, string delay, string prob, int k, double eps,
		double ell) {
	const double e = exp(1);
	double log_VCk = log_nCk(V, k);

	ell = ell * (1 + log(2) / log(V));
	double eps_p = sqrt(2) * eps;

	printf("ell  = %f\n", ell);
	printf("eps' = %f\n", eps_p);
	printf("log{V c k} = %f\n", log_VCk);

	double OPT_lb = 1;

	int H = 0;
	vector<vector<int> > h2v; // h2v[i][j]: j-th vertex in the i-th hyperedge
	vector<vector<int> > v2h(V); // v2h[v][h]: a hyperedge h has a vertex v

	long long int totW = 0;

	for (int i = 1; i <= log2(V) - 1; i++) {
		double x = V / pow(2, i);
		double lambda_prime = (2 + 2.0 / 3.0 * eps_p)
				* (log_VCk + ell * log(V) + log(log2(V))) * V / (eps_p * eps_p);
		double theta_i = lambda_prime / x;

		printf("i = %d\n", i);
		printf("x  = %.0f\n", x);
		printf("lambda_prime = %.0f\n", lambda_prime);
		printf("theta_i = %.0f\n", theta_i);

		for (; H <= theta_i; H++) {
			vector<int> RR;
			if (model == "tvic") {
				RR = gen_RR_TVIC(V, rs, ideg, delay, prob);
			}
			if (model == "tvlt") {
				RR = gen_RR_TVLT(V, rs, ideg, delay, prob);
			}
			h2v.push_back(vector<int>());
			for (int v : RR) {
				h2v[H].push_back(v);
				v2h[v].push_back(H);
				totW++;
			}
			vector<int>().swap(RR);
		}
		printf("H  = %d\n", H);
		printf("totW = %lld\n", totW);

		vector<int> S;
		int degS = greedy(V, h2v, v2h, k, S);
		printf("deg(S) = %d\n", degS);
		printf("Inf(S) = %f\n", 1.0 * V * degS / H);
		printf("\n");
		if (1.0 * V * degS / theta_i >= (1 + eps_p) * x) {
			OPT_lb = (1.0 * V * degS) / ((1 + eps_p) * theta_i);
			break;
		}
	}
	double lambda_star;
	{
		double alpha = sqrt(ell * log(V) + log(2));
		double beta = sqrt((1 - 1 / e) * (log_VCk + ell * log(V) + log(2)));
		double c = (1 - 1 / e) * alpha + beta;
		lambda_star = 2 * V * c * c / (eps * eps);
	}
	double theta = lambda_star / OPT_lb;
	printf("OPT_ = %.0f\n", OPT_lb);
	printf("lambda* = %.0f\n", lambda_star);
	printf("theta = %.0f\n", theta);

	for (; H <= theta; H++) {
		vector<int> RR;
		if (model == "tvic") {
			RR = gen_RR_TVIC(V, rs, ideg, delay, prob);
		}
		if (model == "tvlt") {
			RR = gen_RR_TVLT(V, rs, ideg, delay, prob);
		}
		h2v.push_back(vector<int>());
		for (int v : RR) {
			h2v[H].push_back(v);
			v2h[v].push_back(H);
		}
	}
	printf("H  = %d\n", H);

	vector<int> S;
	int degS = greedy(V, h2v, v2h, k, S);
	printf("deg(S) = %d\n", degS);
	printf("Inf(S) = %f\n", 1.0 * V * degS / H);

	return S;
}

void run(map<string, string> args) {
	string input = get_or_die(args, "graph");
	int k = atoi(get_or_die(args, "k").c_str());
	string alg = get_or_die(args, "alg");
	double eps = atof(get_or_die(args, "eps").c_str());
	double ell = atof(get_or_die(args, "ell").c_str());
	string model = get_or_die(args, "model");
	string delay = get_or_die(args, "delay");
	string prob = get_or_die(args, "prob");
	int numMC = atoi(get_or_die(args, "numMC").c_str());

	ifstream is(input.c_str());
	vector<edge> ps;
	int V = 0;
	int u, v;
	double a, b, c;
	for (; is >> u >> v >> a >> b >> c;) {
		if (u == v) {
			continue;
		}
		edge e = { u, v, a, b, c };
		V = max(V, max(u, v) + 1);
		ps.push_back(e);
	}

	vector<vector<edge> > rs(V), es(V);
	vector<int> ideg(V);
	for (auto e : ps) {
		rs[e.v].push_back(e);
		es[e.u].push_back(e);
		ideg[e.v]++;
	}

	vector<int> S;
	if (alg == "imm") {
		S = imm(V, rs, ideg, model, delay, prob, k, eps, ell);
	} else if (alg == "imm-ctic") {
		double deadline = atof(get_or_die(args, "deadline").c_str());
		vector<vector<edge> > rs1(V);
		for (auto e : ps) {
			edge e1 = e;
			e1.c = deadline;
			rs1[e1.v].push_back(e1);
		}
		S = imm(V, rs1, ideg, model, delay, "deadline", k, eps, ell);
	} else if (alg == "imm-ic") {
		S = imm(V, rs, ideg, model, delay + "-ic", prob + "-ic", k, eps, ell);
	} else if (alg == "imm-lt") {
		S = imm(V, rs, ideg, model, delay + "-lt", prob + "-lt", k, eps, ell);
	} else {
		printf("Illegal algorithm name\n");
	}

	for (int i = 0; i < S.size(); i++) {
		printf("%3d-th seed = %d\n", i, S[i]);
	}

	double inf = 0;
	for (int sim = 0; sim < numMC; sim++) {
		if (model == "tvic") {
			inf += MonteCarlo_TVIC(V, es, ideg, delay, prob, S);
		}
		if (model == "tvlt") {
			inf += MonteCarlo_TVLT(V, es, ideg, delay, prob, S);
		}
	}
	inf /= numMC;
	printf("average influence spread = %.f\n", inf);
}

int main(int argc, char *argv[]) {
	init_genrand(1089267);

	map<string, string> args;
	init_args(argc, argv, args);

	run(args);

	return 0;
}
