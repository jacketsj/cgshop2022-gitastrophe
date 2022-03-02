#pragma once

#include <bits/stdc++.h>

#include "encodings.h"

using namespace std;

void save_folder(string filename, string folder_to_save) {
	instance i;
	i.load_json(filename);
	// i.compute_crossings();
	solution s(i);
	s.from_json_folder(folder_to_save);
	/*
	for (int j = 0; j < 7; ++j) {
		s.compute_greedy_update_random_shuffle();
		merge_random_cols(20, 4, i, s);
		// memory hog and not useful:
		// vector<int> perm2 = random_planar_decomp_permutation(i);
		// s.compute_greedy_update(perm2);
	}
	*/
	cout << filename << ": " << (s.verify() ? "verified" : "bad") << endl;
	s.save_if_better_than_default();
}
