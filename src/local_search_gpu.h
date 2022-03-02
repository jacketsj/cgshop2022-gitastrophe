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

#if __has_include(<CL/sycl.hpp>)
#include <CL/sycl.hpp>
struct local_search_improver_gpu {
	instance ins;
	solution s;
	mt19937 rng;
	queue<size_t> bad_edges;
	vector<char> bad;
	int num_badify;
	local_search_improver_gpu(string filename, int badify = 100)
			: ins(filename), s(ins), num_badify(badify) {
		s.from_json_default();
		rng = mt19937(time(0));
		bad.resize(ins.m);
		init_bad_colour();
		run();
	}
	void init_bad_colour() {
		size_t bad_colour = s.num_cols - 1;
		for (size_t i = 0; i < ins.m; i++) {
			if (s.cols[i] == bad_colour) {
				bad_edges.push(i);
				bad[i] = 1;
			}
		}
	}
	void run() {
		/*
		cl::sycl::queue q;
		cl::sycl::buffer<decltype(ins.edge_crossings[0])> edge_crossings_buf(
				ins.edge_crossings.data(), ins.m);

		*/
		int iter = 0;
		while (1) {
			if (iter++ % 1000 == 0) {
				cerr << "Running iteration " << iter << ": " << bad_edges.size()
						 << " remaining.\n";
			}
			// decolour some random stuff
			auto ncols = s.cols;
			auto nbad_edges = bad_edges;
			auto nbad = bad;
			for (int i = 0; i < num_badify; i++) {
				int k = rng() % ins.m;
				while (nbad[k]) {
					k = rng() % ins.m;
				}
				nbad[k] = 1;
				nbad_edges.push(k);
			}
			int sz = nbad_edges.size();
			// try to colour stuff
			// FIRST STEP: DO this in parallel (why not)
			vector<int> available_cols_all(sz * (s.num_cols - 1));
			vector<int> available_cols_lists(sz * (s.num_cols - 1));
			vector<mt19937> rng_arr(sz);
			vector<int> ncol_arr(sz, -1);
			{
				/*
				cl::sycl::buffer<bool, 2> available_cols_all_buf(
						available_cols_all.data(), range<2>{sz, s.num_cols - 1});
				cl::sycl::buffer<bool, 2> available_cols_lists_buf(
						available_cols_all.data(), range<2>{sz, s.num_cols - 1});
				cl::sycl::buffer<mt19937> rng_arr_buf(rng_arr.data(), sz);
				cl::sycl::buffer<int> ncol_arr_buf(ncol_arr.data(), sz);
				cl::sycl::buffer<int> nbad_edges_buf(nbad_edges.data(), sz);
				q.submit([&](cl::sycl::handler& cgh) {
					auto available_cols_all_acc =
							available_cols_all_buf.get_access<cl::sycl::access::mode::read>(
									cgh);
					auto available_cols_lists_acc =
							available_cols_lists_buf.get_access<cl::sycl::access::mode::read>(
									cgh);
					auto rng_arr_acc =
							rng_arr_buf.get_access<cl::sycl::access::mode::read>(cgh);
					auto nbad_edges_acc =
							nbad_edges_buf.get_access<cl::sycl::access::mode::read>(cgh);
					auto ncol_arr_acc =
							ncol_arr_buf.get_access<cl::sycl::access::mode::read_write>(cgh);
					auto edge_crossings_acc =
							edge_crossings_buf.get_access<cl::sycl::access::mode::read>(cgh);
					cgh.parallel_for<class do_stuff>(
							cl::sycl::range<1>{sz}, [=](cl::sycl::item<1> item_id) {
								int e = nbad_edges[item_id];
								vector<bool> available_cols(s.num_cols - 1, true);
								for (int f = ins.edge_crossings[e]._Find_first(); f < M;
										 f = ins.edge_crossings[e]._Find_next(f)) {
									if (!nbad[f]) {
										available_cols_acc[item_id][ncols[f]] = false;
									}
								}
								vector<int> avail;
								for (size_t c = 0; c < s.num_cols - 1; c++) {
									if (available_cols[c]) {
										avail.push_back(c);
									}
								}
								if (!avail.empty()) {
									nbad[e] = 0;
									ncols[e] = avail[rng() % avail.size()];
								} else {
									nbad[e] = 1;
									nbad_edges.push(e);
								}
							});
				});
				*/
			}
			// good if number of bad edges decreases
			if (nbad_edges.size() <= bad_edges.size()) {
				bad_edges = nbad_edges;
				bad = nbad;
				s.cols = ncols;
				if (bad_edges.size() == 0) {
					cerr << "Successful improvement! " << s.num_cols << " -> ";
					s.recompute_num_cols();
					cerr << s.num_cols << endl;
					assert(s.verify());
					s.save_if_better_than_default();
					init_bad_colour();
				}
			}
		}
	}
};
#endif
