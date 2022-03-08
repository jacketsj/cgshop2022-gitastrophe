#pragma once

#define SOLUTIONS_FOLDER "temp/"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "encodings.h"

#include <nlohmann/json.hpp>

#include "verify/bentley_ottmann_any_intersection.hpp"
#include "verify/geom.hpp"
#include "verify/verify_coloring.hpp"

struct conflict_improver {
	instance ins;
	solution s;
	mt19937 rng;
	vector<int> bad_edges;
	vector<char> bad;
	vector<int64_t> qcnt;
	std::chrono::time_point<std::chrono::steady_clock> limit;
	conflict_improver(string filename, int time_limit = 1e9)
			: ins(filename), s(ins) {
		s.from_json_default();
		rng = mt19937(time(0));
		bad.resize(ins.m);
		qcnt.resize(ins.m);
		limit = std::chrono::steady_clock::now() + std::chrono::seconds(time_limit);
		init_bad_colour();
		run();
	}
	conflict_improver(solution& _s, int time_limit = 1e9) : ins(_s.ins), s(_s) {
		// s.from_json_default();
		rng = mt19937(time(0));
		bad.resize(ins.m);
		qcnt.resize(ins.m);
		limit = std::chrono::steady_clock::now() + std::chrono::seconds(time_limit);
		init_bad_colour();
		run();
	}
	void init_bad_colour() {
		fill(begin(bad), end(bad), 0);
		fill(begin(qcnt), end(qcnt), 0);
		for (size_t e = 0; e < ins.m; e++) {
			qcnt[e] = ins.deg[e];
		}
		while (!bad_edges.empty())
			bad_edges.pop_back();
		vector<int> cnt(s.num_cols);
		for (size_t i = 0; i < ins.m; i++) {
			cnt[s.cols[i]]++;
		}
		int mn = *min_element(begin(cnt), end(cnt));
		cerr << "class size: " << mn << endl;
		size_t bad_colour = rng() % s.num_cols;
		while (cnt[bad_colour] > mn)
			bad_colour = rng() % s.num_cols;
		// bad_colour = s.num_cols-1;
		// relabel bad to s.num_cols - 1
		for (size_t i = 0; i < ins.m; i++) {
			if (s.cols[i] == bad_colour) {
				bad_edges.push_back(i);
				bad[i] = 1;
				s.cols[i] = s.num_cols - 1;
			} else if (s.cols[i] == s.num_cols - 1) {
				s.cols[i] = bad_colour;
			}
		}
	}
	void run() {
		while (1) {
			int iter = 0;
			int64_t mx_qcnt = 0, mx_qsz = 0, restarts = 0;
			bool mode = 0;
			int switchiter = 1e4;
			auto bphase = s.cols;
			vector<int> bestq;
			while (!bad_edges.empty()) {
				if (iter % 1000000 == 0) {
					cerr << setw(18) << setfill(' ') << ins.id << " Running iteration "
							 << iter << ": " << bad_edges.size() << " edges remaining.\n";
				}
				/*
				if (iter % 1000000 == 0) {
					for (auto& i : qcnt) i = (i+1)/2;
				}*/
				if (iter++ % 10000 == 0 && std::chrono::steady_clock::now() >= limit) {
					return;
					// cerr << "Running iteration " << iter << ": " << bad_edges.size() <<
					// " remaining.\n";
				}
				/*
				if (iter >= switchiter - 1e4 && (bestq.empty() || bad_edges.size() <=
				bestq.size())) { bphase = s.cols; bestq  = bad_edges;
				}*/
				mx_qsz = max(mx_qsz, (int64_t)bad_edges.size());
				swap(bad_edges[rng() % bad_edges.size()], bad_edges.back());
				int e = bad_edges.back();
				bad_edges.pop_back();
				if (mode)
					qcnt[e]++;
				mx_qcnt = max(mx_qcnt, qcnt[e]);
				// if stuck, restart from a different colour
				if (qcnt[e] > 1e5) {
					// cerr << "restarting\n";
					s.from_json_default();
					init_bad_colour();
					mx_qcnt = 0;
					mx_qsz = 0;
					restarts++;
					continue;
				}
				// cerr << "try " << e << " " << ins.m << endl;
				size_t c;
				if (mode) {
					vector<int64_t> conflict_score(s.num_cols - 1);
					for (int f = ins.edge_crossings[e]._Find_first(); f < M;
							 f = ins.edge_crossings[e]._Find_next(f)) {
						if (!bad[f] && s.cols[f] < s.num_cols - 1) {
							conflict_score[s.cols[f]] += 1000 + 1LL * qcnt[f] * qcnt[f];
						}
					}
					int64_t mnscore =
							*min_element(begin(conflict_score), end(conflict_score));
					c = rng() % (s.num_cols - 1);
					while (conflict_score[c] > 1.05 * mnscore) {
						c = rng() % (s.num_cols - 1);
					}
				} else {
					vector<int> num_confs(s.num_cols - 1);
					for (int f = ins.edge_crossings[e]._Find_first(); f < M;
							 f = ins.edge_crossings[e]._Find_next(f)) {
						if (!bad[f] && s.cols[f] < s.num_cols - 1) {
							num_confs[s.cols[f]]++;
						}
					}
					int mnscore = *min_element(begin(num_confs), end(num_confs));
					c = rng() % (s.num_cols - 1);
					while (num_confs[c] > mnscore) {
						c = rng() % (s.num_cols - 1);
					}
				}
				/*
				int64_t mnscore = *min_element(begin(conflict_score),
				end(conflict_score)); size_t c = rng() % (s.num_cols - 1); while
				(conflict_score[c] > 1.1*mnscore) { c = rng() % (s.num_cols - 1);
				}*/
				bad[e] = 0;
				s.cols[e] = c;
				for (int f = ins.edge_crossings[e]._Find_first(); f < M;
						 f = ins.edge_crossings[e]._Find_next(f)) {
					if (!bad[f] && s.cols[f] == c) {
						bad[f] = 1;
						bad_edges.push_back(f);
					}
				}
				if (iter == switchiter) {
					// s.cols = bphase;
					// bad_edges = bestq;
					// bestq.clear();
					// bad = vector(ins.m, (char)0);
					// for (int i : bad_edges) bad[i] = 1;
					mode ^= 1;
					if (mode == 1)
						switchiter += 2e4;
					else
						switchiter += 1e4;
				}
			}
			cerr << "Successful improvement! " << s.num_cols << " -> ";
			s.recompute_num_cols();
			cerr << s.num_cols << " mx_qcnt=" << mx_qcnt << " mx_qsz=" << mx_qsz
					 << " restarts=" << restarts << " num_iters=" << iter << endl;
			assert(s.verify());
			s.save_if_better_than_default();
			init_bad_colour();
		}
	}
};

struct conflict_improver_2 {
	instance ins;
	solution& s;
	mt19937& rng;
	queue<size_t> bad_edges;
	vector<char> bad;
	vector<int64_t> qcnt;
	vector<size_t> cols;
	std::chrono::time_point<std::chrono::steady_clock> limit;
	conflict_improver_2(solution& _s, int time_limit, mt19937& _rng)
			: ins(_s.ins), s(_s), rng(_rng) {
		bad.resize(ins.m);
		qcnt.resize(ins.m);
		limit = std::chrono::steady_clock::now() +
						std::chrono::milliseconds(time_limit);
		cols = s.cols;
		init_bad_colour();
		run_sub();
	}
	void init_bad_colour() {
		fill(begin(bad), end(bad), 0);
		fill(begin(qcnt), end(qcnt), 0);
		while (!bad_edges.empty())
			bad_edges.pop();
		size_t bad_colour = s.num_cols - 1;
		// relabel bad to s.num_cols - 1
		for (size_t i = 0; i < ins.m; i++) {
			if (s.cols[i] == bad_colour) {
				bad_edges.push(i);
				bad[i] = 1;
			} else if (s.cols[i] == s.num_cols - 1) {
				s.cols[i] = bad_colour;
			}
		}
	}
	void run_sub() {
		// cerr << "run_sub\n";
		if (s.num_cols == 1)
			return;
		while (1) {
			int iter = 0;
			int64_t mx_qcnt = 0, mx_qsz = 0, restarts = 0;
			int mn_qsz = bad_edges.size();
			queue<size_t> bestq;
			while (!bad_edges.empty()) {
				if (iter % 10000 == 0) {
					cerr << "Running iteration " << iter << ": " << bad_edges.size()
							 << " remaining.\n";
				}
				if (iter++ % 200 == 0 && std::chrono::steady_clock::now() >= limit) {
					s.cols = cols;
					// bad_edges = bestq;
					cerr << "returning after " << iter << " iters, " << bad_edges.size()
							 << " remaining\n";
					/*
					int c = s.num_cols-1;
					vector<int> bades;
					while (!bad_edges.empty()) {
						int e = bad_edges.front(); bad_edges.pop();
						bades.push_back(e);
						s.cols[e] = c++;
					}
					for (int i : bades) bad_edges.push(i);
					assert(s.verify());*/
					return;
				}
				mx_qsz = max(mx_qsz, (int64_t)bad_edges.size());
				if (mn_qsz > bad_edges.size()) {
					mn_qsz = bad_edges.size();
					bestq = bad_edges;
				}
				int e = bad_edges.front();
				bad_edges.pop();
				qcnt[e]++;
				mx_qcnt = max(mx_qcnt, qcnt[e]);
				// if stuck, restart from a different colour
				/*
				if (qcnt[e] > 5e4) {
					//cerr << "restarting\n";
					s.from_json_default();
					init_bad_colour();
					mx_qcnt = 0;
					mx_qsz = 0;
					restarts++;
					continue;
				}*/
				// cerr << "try " << e << " " << ins.m << endl;
				vector<int64_t> conflict_score(s.num_cols - 1, 0);
				for (int f = ins.edge_crossings[e]._Find_first(); f < M;
						 f = ins.edge_crossings[e]._Find_next(f)) {
					if (!bad[f] && s.cols[f] < s.num_cols - 1) {
						conflict_score[s.cols[f]] += 10 + 1LL * qcnt[f] * qcnt[f];
					}
				}
				int64_t mnscore =
						*min_element(begin(conflict_score), end(conflict_score));
				size_t c = rng() % (s.num_cols - 1);
				while (conflict_score[c] > 1.1 * mnscore) {
					c = rng() % (s.num_cols - 1);
				}
				bad[e] = 0;
				s.cols[e] = c;
				vector<int> edges;
				for (int f = ins.edge_crossings[e]._Find_first(); f < M;
						 f = ins.edge_crossings[e]._Find_next(f)) {
					if (!bad[f] && s.cols[f] == c) {
						edges.push_back(f);
					}
				}
				shuffle(begin(edges), end(edges), rng);
				for (int f : edges) {
					bad[f] = 1;
					bad_edges.push(f);
				}
			}
			s.recompute_num_cols();
			cerr << "improved -> " << s.num_cols << endl;
			;
			if (s.num_cols == 1)
				return;
			cols = s.cols;
			return;
			init_bad_colour();
		}
	}
};
