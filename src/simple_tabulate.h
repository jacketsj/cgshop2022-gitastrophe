#include <bits/stdc++.h>

#include "encodings.h"

using namespace std;

void simple_tabulate(vector<string> files) {
	vector<tuple<double, int, int, string>> scores;
	for (string f : files) {
		instance ins(f, false);
		solution s(ins);
		s.from_json_default();
		int lb = ceil(1.L * ins.m / (3*ins.n - 6));
		int score = s.num_cols;
		scores.emplace_back(1.L * lb / score, lb, score, f);
	}
	sort(begin(scores), end(scores));
	for (auto [a, b, c, d] : scores) {
		cout << setw(40) << setfill(' ') << d << setw(10) << fixed << setprecision(6) << a << setw(8) << b << setw(8) << c << endl;
	}
}
