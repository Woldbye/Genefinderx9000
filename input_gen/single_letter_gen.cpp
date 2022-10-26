/**
 * @file in_gen.cpp
 * @brief Generates a dictionary of only A's and query of A's such that there will be exit links
 * @warning can take quite a long time to generate for large numbers
 */
#include <bits/stdc++.h>
#include <functional>
#include <experimental/iterator>
#include <chrono>
using namespace std;
using namespace chrono;

int main() {
	const char C = 'T';
	int N = 1000000, M = 1000000;
	cout << 1 << endl;
	cout << string(N, C) << endl;
	int K = 0, l = 1;

	set<string> out;
	set<int> seen;
	ostringstream os;
	int m;
	for (m = 0; m + l < M; ++K, l = l * 2 + 1)
	{
		seen.insert(l);
		os << string(l, C) << "\n";
	}
	if (m != M-1) {
		while (seen.count(M-m) && m > 0) {
			m--;
		}
		++K;
		seen.insert(M - m);
		os << string(M-m, C) << "\n";
	}
	cout << K << "\n"
			 << os.str();
}