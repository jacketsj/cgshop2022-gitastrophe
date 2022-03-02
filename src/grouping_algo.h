#pragma once

#include <bits/stdc++.h>

#include "encodings.h"

using namespace std;

template <typename T> vector<size_t> k_means(size_t k, vector<T> x) {
	vector<size_t> ret(x.size());
	random_shuffle(x.begin(), x.end());
	vector<pt> u(k);
	for (size_t i = 0; i < k; ++i)
		u[i] = x[i].first;

	auto length = [](pt a, pt b) { return (a - b).length2(); };
	bool changed = true;
	int iters = 20000;
	cout << "performing k means" << endl;
	while (changed && iters-- > 0) {
		changed = false;
		for (int i = 0; i < x.size(); ++i)
			for (int j = 0; j < k; ++j)
				if (length(x[i].first, u[ret[i]]) > length(x[i].first, u[j])) {
					changed = true;
					ret[i] = j;
				}
		vector<vector<T>> partitions(k);
		for (size_t i = 0; i < x.size(); ++i)
			partitions[ret[i]].push_back(x[i]);
		for (size_t i = 0; i < k; ++i) {
			pt combined;
			for (auto& [p, _] : partitions[i])
				combined = combined + p;
			combined = combined / partitions[i].size();
			u[i] = combined;
		}
	}
	cout << "done performing k means" << endl;

	return ret;
}

// try to take a large merger of two colours
// returns an order to be applied greedily
vector<int> try_merge_cols(size_t k, instance& ins, solution& s,
													 vector<size_t> to_combine) {
	// for each edge in the two colours, assign a random point along the edge as
	// representative
	vector<pt> rep(ins.m);
	for (size_t i = 0; i < ins.m; ++i) {
		// get a random double from 0 to 1
		double alpha = double(rand()) / double(RAND_MAX - 1);
		rep[i] = ins.pts[ins.edges[i].first] * alpha +
						 ins.pts[ins.edges[i].second] * (1 - alpha);
	}
	// do k-means
	vector<size_t> non_included;
	set<size_t> to_combine_set;
	for (auto& c : to_combine)
		to_combine_set.insert(c);
	vector<pair<pt, size_t>> x;
	for (size_t i = 0; i < ins.m; ++i) {
		if (to_combine_set.count(s.cols[i]) > 0)
			x.emplace_back(rep[i], i);
		else
			non_included.push_back(i);
	}
	vector<size_t> mean = k_means(k, x);
	vector<vector<size_t>> partitions(k);
	for (int i = 0; i < x.size(); ++i) {
		partitions[mean[i]].push_back(x[i].second);
	}

	vector<size_t> majorities;
	vector<size_t> non_majorities;
	for (int i = 0; i < k; ++i) {
		// get the most frequent color of cluster k
		unordered_map<int, int> freq;
		int best = 0;
		for (int e : partitions[i]) {
			best = max(best, ++freq[s.cols[e]]);
		}
		int best_col = 0;
		for (auto& [col, f] : freq) {
			if (f == best) {
				best_col = col;
				break;
			}
		}
		for (int e : partitions[i]) {
			if (s.cols[e] == best_col)
				majorities.push_back(e);
			else
				non_majorities.push_back(e);
		}
	}
	random_shuffle(majorities.begin(), majorities.end());
	random_shuffle(non_majorities.begin(), non_majorities.end());
	random_shuffle(non_included.begin(), non_included.end());
	vector<int> ret;
	for (int i : majorities)
		ret.push_back(i);
	for (int i : non_majorities)
		ret.push_back(i);
	for (int i : non_included)
		ret.push_back(i);
	return ret;
}

void merge_random_cols(size_t k, size_t num_to_combine, instance& ins,
											 solution& s) {
	vector<size_t> col_perm(s.num_cols);
	for (int i = 0; i < s.num_cols; ++i)
		col_perm[i] = i;
	random_shuffle(col_perm.begin(), col_perm.end());
	col_perm.resize(min(num_to_combine, col_perm.size()));
	vector<int> recolor_order = try_merge_cols(k, ins, s, col_perm);
	solution s2(s);
	int uncol = s2.num_cols;
	set<int> recols;
	for (int i : col_perm)
		recols.insert(i);
	s2.num_cols++;
	for (int i = 0; i < ins.m; ++i) {
		if (recols.count(s.cols[i]) > 0)
			s2.cols[i] = uncol;
	}
	s2.compute_greedy_update(recolor_order);
	if (s2.is_better(s)) {
		cout << "improvement found via k-means thingy: " << s.num_cols << " -> "
				 << s2.num_cols << endl;
		s.cols = s2.cols;
		s.num_cols = s2.num_cols;
		s.save_if_better_than_default();
	} else {
		cout << "improvement not found via k-means thingy: " << s.num_cols << " -> "
				 << s2.num_cols << endl;
	}
}
