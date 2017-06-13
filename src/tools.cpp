#include "./tools.h"

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

using namespace std;

void init_args(int argc, char *argv[], map<string, string> &args) {
	for (int i = 1; i < argc; i++) {
		string a(argv[i]);
		if (a[0] != '-') {
			continue;
		}
		int at = a.find("=");
		string key = a.substr(1, at - 1);
		string val = a.substr(at + 1);
		args[key] = val;
	}
}

string get_or_die(map<string, string> &argv, string key) {
	if (argv.count(key) == 0) {
		printf("Illegal key %s (get_or_die)\n", key.c_str());
		exit(1);
	}
	return argv[key];
}

int genrand_int(int n) {
	return (int) (genrand_res53() * n);
}
