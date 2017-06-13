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

#include <sys/time.h>
#include <unistd.h> // getpid
#include <memory> // auto_ptr

#include "./mt19937ar.h"

using namespace std;

void init_args(int argc, char *argv[], map<string, string> &args);
string get_or_die(map<string, string> &argv, string key);
int genrand_int(int n);
