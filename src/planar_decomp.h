#pragma once

#include <bits/stdc++.h>

#include "encodings.h"

using std::unordered_set;
using std::vector;

vector<int> random_planar_decomp_permutation(instance& ins) {
	// compute a set of planar graphs in ins
	vector<unordered_set<int>> planar_graphs(1);
	// get a random permutation of the edges
	vector<size_t> p;
	for (size_t i = 0; i < ins.m; ++i)
		p.push_back(i);
	std::random_shuffle(p.begin(), p.end());
	for (size_t e : p) {
		bool next = false;
		for (size_t e2 = ins.edge_crossings[e]._Find_first(); e2 < M; e2 = ins.edge_crossings[e]._Find_next(e2))
			if (planar_graphs.back().count(e2) > 0) {
				next = true;
				break;
			}
		if (next)
			planar_graphs.emplace_back();
		planar_graphs.back().insert(e);
	}

	vector<vector<vector<int>>> adj_arr(planar_graphs.size(),
																			vector<vector<int>>(ins.n));
	for (size_t layer = 0; layer < planar_graphs.size(); ++layer) {
		auto& pg = planar_graphs[layer];
		for (auto& e : pg) {
			size_t u, v;
			tie(u, v) = ins.edges[e];
			adj_arr[layer][u].push_back(v);
			adj_arr[layer][v].push_back(u);
		}
	}
	// dumb planar graph colouring strategy: sort by 'weight': sum of degrees
	// TODO: at least remove the min weight edge to get new weights...
	vector<int> ret;
	for (size_t layer = 0; layer < planar_graphs.size(); ++layer) {
		auto& pg = planar_graphs[layer];
		vector<pair<size_t, size_t>> w;
		for (auto& e : pg) {
			size_t u, v;
			tie(u, v) = ins.edges[e];
			w.emplace_back(adj_arr[layer][u].size() + adj_arr[layer][v].size(), e);
		}
		sort(w.begin(), w.end());
		for (auto& [_, e] : w)
			ret.push_back(e);
	}
	return ret;
}
