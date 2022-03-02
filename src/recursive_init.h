#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "conflict_optimizer.h"
#include "encodings.h"

std::mt19937 rng(time(0));

instance subinstance(instance& ins, const vector<size_t>& es) {
	instance sub;
	sub.n = ins.n;
	sub.m = es.size();
	sub.pts = ins.pts;
	for (size_t i = 0; i < sub.m; i++) {
		sub.edges.push_back(ins.edges[es[i]]);
	}
	sub.edge_crossings.resize(sub.m);
	for (size_t i = 0; i < es.size(); i++) {
		for (size_t j = i + 1; j < es.size(); j++) {
			if (ins.edge_crossings[es[i]][es[j]]) {
				sub.edge_crossings[i][j] = 1;
				sub.edge_crossings[j][i] = 1;
			}
		}
	}
	return sub;
}

solution recur(instance& ins, int depth = 0) {
	solution ret(ins);
	ret.compute_greedy_update_sorted();
	if (ins.m < 50) {
		conflict_improver_2 conf(ret, ins.m < 5 ? 1 : 5, rng);
		return ret;
	}
	auto score = [](int lcnt, int rcnt, int crosscnt) {
		return std::max(lcnt, rcnt) + crosscnt;
	};
	ls best;
	int bestscore = 1e9;
	const int C = 1e8;
	for (int tt = 0; tt < 2000; tt++) {
		pt s(rng() % C, rng() % C), t(rng() % C, rng() % C);
		ls l = {s, t};
		size_t lcnt = 0, rcnt = 0, crosscnt = 0;
		for (edge e : ins.edges) {
			ls seg = ins.edge_to_ls(e);
			if (l.cross_line(seg)) {
				crosscnt++;
			} else {
				// assert(winding(s, t, seg.a) != 0);
				if (winding(s, t, seg.a) < 0)
					lcnt++;
				else
					rcnt++;
			}
		}
		if (lcnt < ins.m && rcnt < ins.m && crosscnt < ins.m &&
				score(lcnt, rcnt, crosscnt) < bestscore) {
			bestscore = score(lcnt, rcnt, crosscnt);
			best = l;
		}
	}
	// split using best
	vector<size_t> l, r, c;
	for (size_t i = 0; i < ins.m; i++) {
		edge e = ins.edges[i];
		ls seg = ins.edge_to_ls(e);
		if (best.cross_line(seg)) {
			c.push_back(i);
		} else {
			if (winding(best.a, best.b, seg.a) < 0)
				l.push_back(i);
			else
				r.push_back(i);
		}
	}
	if (l.size() == ins.m || r.size() == ins.m || c.size() == ins.m)
		return ret;
	instance subl = subinstance(ins, l), subr = subinstance(ins, r),
					 subc = subinstance(ins, c);
	solution soll = recur(subl, depth + 1), solr = recur(subr, depth + 1),
					 solc = recur(subc, depth + 1);
	cerr << ins.m << " @ depth " << depth << ": " << subl.m << " left, " << subr.m
			 << " right, " << subc.m << " crossing, "
			 << "recursive=" << std::max(soll.num_cols, solr.num_cols) << "+"
			 << solc.num_cols << "="
			 << std::max(soll.num_cols, solr.num_cols) + solc.num_cols
			 << ", greedy=" << ret.num_cols << endl;
	// for (int i : l) for (int j : r) assert(!ins.edge_crossings[i][j]);
	vector<size_t> cols(ins.m);
	// merge
	int off = std::max(soll.num_cols, solr.num_cols);
	for (size_t i = 0; i < subl.m; i++) {
		assert(soll.cols[i] < off);
		cols[l[i]] = soll.cols[i];
	}
	for (size_t i = 0; i < subr.m; i++) {
		assert(solr.cols[i] < off);
		cols[r[i]] = solr.cols[i];
	}
	for (size_t i = 0; i < subc.m; i++) {
		cols[c[i]] = off + solc.cols[i];
	}
	// try greedily incorporate c
	int cur = off;
	vector<size_t> perm(c.size());
	iota(begin(perm), end(perm), 0);
	shuffle(begin(perm), end(perm), rng);
	for (size_t col = 0; col < solc.num_cols; col++) {
		for (size_t i : perm) {
			if (solc.cols[i] != col)
				continue;
			vector<int> avail(cur, 1);
			for (int j = ins.edge_crossings[c[i]]._Find_first(); j < M;
					 j = ins.edge_crossings[c[i]]._Find_next(j)) {
				if (cols[j] < cur) {
					avail[cols[j]] = 0;
				}
			}
			for (size_t d = 0; d < cur; d++) {
				if (avail[d]) {
					cols[c[i]] = d;
					break;
				}
			}
		}
	}
	// relabel colours
	vector<int> mapping(*max_element(begin(cols), end(cols)) + 1, -1);
	int id = 0;
	for (size_t i = 0; i < ins.m; i++) {
		if (mapping[cols[i]] == -1) {
			mapping[cols[i]] = id++;
		}
		cols[i] = mapping[cols[i]];
	}
	int prv = ret.num_cols;
	if (id <= prv) {
		ret.cols = cols;
		ret.recompute_num_cols();
	}
	prv = ret.num_cols;
	conflict_improver_2 conf(ret, ins.m < 100 ? 1 : 10, rng);
	cerr << "conflicted " << prv << " -> " << ret.num_cols << endl;
	return ret;
}
