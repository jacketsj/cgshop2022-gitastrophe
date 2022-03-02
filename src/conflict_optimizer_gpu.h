#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "encodings.h"

#if __has_include(<CL/sycl.hpp>)
#include <CL/sycl.hpp>
using namespace cl;
#include <syclflow.hpp>

// better seeding apparently
// https://codereview.stackexchange.com/questions/109260/seed-stdmt19937-from-stdrandom-device
template <class T = std::mt19937,
					std::size_t N = T::state_size * sizeof(typename T::result_type)>
T ProperlySeededRandomEngine() {
	std::random_device source;
	std::random_device::result_type random_data[(N - 1) / sizeof(source()) + 1];
	std::generate(std::begin(random_data), std::end(random_data),
								std::ref(source));
	std::seed_seq seeds(std::begin(random_data), std::end(random_data));
	return T(seeds);
}

template <typename T>
// using gpu_alloc = sycl::usm_allocator<T, sycl::usm::alloc::shared>;
using gpu_alloc = std::allocator<T>;
// using instance_gpu = instance_templated<gpu_alloc, vector>;
// using solution_gpu = solution_templated<gpu_alloc, vector>;
using instance_gpu = instance;
using solution_gpu = solution;

struct conflict_improver_gpu {
	instance_gpu ins;
	solution_gpu s_main;
	struct batched_jobs_vals {
		size_t m;
		size_t P;
		size_t num_cols;
		size_t* acc_main_cols;
		bitset<M>* acc_edge_crossings;
		size_t* acc_bad_edges;
		int* acc_bad;
		int* acc_qcnt;
		int* acc_conflict_score;
		size_t* acc_cols;
		mt19937* acc_rng;
		int* acc_success;
		template <typename T> T& operator()(size_t t, size_t off, T* acc) {
			return acc[off * P + t];
		}
	};
	struct batched_jobs {
		size_t m;
		size_t P;
		size_t num_cols;
		vector<size_t, gpu_alloc<size_t>> main_cols;
		// size_t* acc_main_cols;
		vector<bitset<M>, gpu_alloc<bitset<M>>> edge_crossings;
		// bitset<M>* acc_edge_crossings;
		// all batch_ vecs are size P*m or P*num_cols
		vector<size_t, gpu_alloc<size_t>> batch_bad_edges;
		vector<int, gpu_alloc<int>> batch_bad;
		vector<int, gpu_alloc<int>> batch_qcnt;
		vector<int, gpu_alloc<int>> batch_conflict_score; // size=(s.num_cols - 1)
		vector<size_t, gpu_alloc<size_t>> batch_cols;
		vector<mt19937, gpu_alloc<mt19937>> batch_rng;
		vector<int, gpu_alloc<int>> batch_success;
		// size_t* acc_bad_edges;
		// int* acc_bad;
		// int* acc_qcnt;
		// int* acc_conflict_score;
		// size_t* acc_cols;
		// mt19937* acc_rng;
		// int* acc_success;
		batched_jobs_vals v;
		template <typename T> T& operator()(size_t t, size_t off, T* acc) {
			return v(t, off, acc);
		}
		//        m(s.ins.m),
		//				P(num_threads),
		//				num_cols(s.num_cols),
		//				acc_main_cols(main_cols.data()),
		//				acc_edge_crossings(edge_crossings.data()),
		//				acc_bad_edges(batch_bad_edges.data()),
		//				acc_bad(batch_bad.data()),
		//				acc_qcnt(batch_qcnt.data()),
		//				acc_conflict_score(batch_conflict_score.data()),
		//				acc_cols(batch_cols.data()),
		//				acc_rng(batch_rng.data()),
		//				acc_success(batch_success.data())
		batched_jobs(size_t num_threads, solution_gpu& s)
				: m(s.ins.m), P(num_threads), num_cols(s.num_cols),
					main_cols(m * P, 0, get_alloc<gpu_alloc, size_t>()),
					edge_crossings(m, bitset<M>(), get_alloc<gpu_alloc, bitset<M>>()),
					batch_bad_edges(m * P, 0, get_alloc<gpu_alloc, size_t>()),
					batch_bad(m * P, 0, get_alloc<gpu_alloc, int>()),
					batch_qcnt(m * P, 0, get_alloc<gpu_alloc, int>()),
					batch_conflict_score((num_cols - 1) * P, 0,
															 get_alloc<gpu_alloc, int>()),
					batch_cols(m * P, 0, get_alloc<gpu_alloc, size_t>()),
					batch_rng(P, mt19937(time(0)), get_alloc<gpu_alloc, mt19937>()),
					batch_success(P, 0, get_alloc<gpu_alloc, int>()),
					v{s.ins.m,
						num_threads,
						s.num_cols,
						main_cols.data(),
						edge_crossings.data(),
						batch_bad_edges.data(),
						batch_bad.data(),
						batch_qcnt.data(),
						batch_conflict_score.data(),
						batch_cols.data(),
						batch_rng.data(),
						batch_success.data()} {
			for (size_t e = 0; e < m; ++e)
				edge_crossings[e] = s.ins.edge_crossings[e];
			for (size_t t = 0; t < P; ++t) {
				for (size_t e = 0; e < m; ++e) {
					batch_cols[e * P + t] = s.cols[e];
					main_cols[e * P + t] = s.cols[e];
				}
				batch_rng[t] = ProperlySeededRandomEngine();
			}
		}
	};
	template <typename T> struct acc_with_t {
		size_t P;
		T* acc;
		size_t t;
		size_t sz;
		T& operator[](size_t off) { return acc[off * P + t]; }
		acc_with_t(size_t P, T* acc, size_t t, size_t sz)
				: P(P), acc(acc), t(t), sz(sz) {}
		size_t size() const { return sz; }
	};
	conflict_improver_gpu(string filename, size_t max_iters, bool from_scratch,
												bool continuous, size_t num_threads = 0)
			: ins(filename), s_main(ins) { //, jobs(get_alloc<gpu_alloc, job>()) {
		// getting the total number of compute units
		auto num_groups =
				q.get_device().get_info<cl::sycl::info::device::max_compute_units>();
		// getting the maximum work group size per thread
		auto work_group_size =
				q.get_device().get_info<cl::sycl::info::device::max_work_group_size>();
		cout << "num_groups=" << num_groups
				 << ",work_group_size=" << work_group_size << endl;
		// building the best number of global thread

		auto total_threads = num_groups;
		size_t group_size = work_group_size;
		// auto total_threads = num_groups * work_group_size;
		if (num_threads > 0) {
			total_threads = num_threads;
			work_group_size = 1; // TODO: is this realy the best idea?
		}
		// should i make it smaller/larger?

		if (from_scratch) {
			s_main.compute_greedy_update_sorted();
			cout << (s_main.verify() ? "verified" : "bad") << endl;
		} else
			s_main.from_json_default();

		do {
			run_gpu(total_threads, max_iters, continuous, group_size);
		} while (continuous);
	}
	// scuffed bad accessor wrapper
	template <typename T> struct bad_acc_thingy {
		const size_t P;
		T a;
		const size_t offset;
		const size_t sz;
		auto operator[](size_t i) { return a[P * i + offset]; }
		bad_acc_thingy(size_t P, T a, size_t thread_id, size_t sz)
				: P(P), a(a), offset(thread_id), sz(sz) {}
		size_t size() { return sz; }
	};
	void run_gpu(size_t total_threads, size_t max_iters_before_reset,
							 bool continuous, size_t group_size) {
		std::cout << "gpu optimizer starting up... (" << total_threads
							<< " threads)" << std::endl;

		/* TODO: Fix this feature
		size_t max_iters = std::min(
				size_t(50),
				max_iters_before_reset); // frequency with which gpu comes back to cpu
		*/
		size_t max_iters = max_iters_before_reset;
		do {
			batched_jobs bj(total_threads, s_main);
			for (size_t unreset_repeat_no = 0;
					 unreset_repeat_no < max_iters_before_reset;
					 unreset_repeat_no += max_iters) {
				{
					sycl::buffer<bitset<M>> buf_edge_crossings(bj.edge_crossings.data(),
																										 bj.edge_crossings.size());
					sycl::buffer<mt19937> buf_rng(bj.batch_rng.data(),
																				bj.batch_rng.size());
					sycl::buffer<size_t> buf_bad_edges(bj.batch_bad_edges.data(),
																						 bj.batch_bad_edges.size());
					sycl::buffer<int> buf_bad(bj.batch_bad.data(), bj.batch_bad.size());
					sycl::buffer<int> buf_qcnt(bj.batch_qcnt.data(),
																		 bj.batch_qcnt.size());
					sycl::buffer<int> buf_conflict_score(bj.batch_conflict_score.data(),
																							 bj.batch_conflict_score.size());
					sycl::buffer<size_t> buf_cols(bj.batch_cols.data(),
																				bj.batch_cols.size());
					sycl::buffer<size_t> buf_main_cols(bj.main_cols.data(),
																						 bj.main_cols.size());
					sycl::buffer<int> buf_success(bj.batch_success.data(),
																				bj.batch_success.size());
					int any_success = 0;
					sycl::buffer<int> buf_any_success(&any_success, 1);

					q.submit([&](sycl::handler& cgh) {
						auto acc_main_cols =
								buf_main_cols.get_access<sycl::access::mode::read_write>(cgh);
						auto acc_edge_crossings =
								buf_edge_crossings.get_access<sycl::access::mode::read_write>(
										cgh);
						auto acc_bad_edges =
								buf_bad_edges.get_access<sycl::access::mode::read_write>(cgh);
						auto acc_bad =
								buf_bad.get_access<sycl::access::mode::read_write>(cgh);
						auto acc_qcnt =
								buf_qcnt.get_access<sycl::access::mode::read_write>(cgh);
						auto acc_conflict_score =
								buf_conflict_score.get_access<sycl::access::mode::read_write>(
										cgh);
						auto acc_cols =
								buf_cols.get_access<sycl::access::mode::read_write>(cgh);
						auto acc_rng =
								buf_rng.get_access<sycl::access::mode::read_write>(cgh);
						auto acc_success =
								buf_success.get_access<sycl::access::mode::read_write>(cgh);
						auto acc_any_success =
								buf_any_success.get_access<sycl::access::mode::read_write>(cgh);
						// cgh.parallel(
						//		sycl::range<1>{total_threads}, sycl::range<1>{group_size},
						//[=](auto group) {
						//[=](sycl::group<1> group, sycl::physical_item<1> phys_idx) {
						cgh.parallel_for(
								// sycl::range<1>{total_threads*group_size},
								sycl::nd_range<1>{total_threads, group_size},
								//[=](sycl::id<1> id) {
								[=](sycl::nd_item<1> item) {
									//, sycl::range<1>{group_size},
									//[=](auto group) {
									// size_t t = group.get_group_id(0);
									// printf("groupid=%d\n", int(t));
									// size_t t = phys_idx.get(0);
									size_t t = item.get_global_id(0) / group_size;
									size_t local_id = item.get_global_id(0) % group_size;
									// size_t t = id.get(0);
									size_t thread_id = t;

									auto m = bj.v.m;
									size_t num_cols = bj.v.num_cols;

									// TODO: Change this for repeats and stuff
									size_t num_cols_main = bj.v.num_cols;

									size_t P = bj.P;

									// use scuffed accessor wrapper
									auto main_cols =
											bad_acc_thingy(bj.P, acc_main_cols, thread_id, bj.m);
									// only one copy of edge_crossings for all threads, pretend
									// bj.P=1
									auto edge_crossings =
											bad_acc_thingy(1, acc_edge_crossings, 0, bj.m);
									auto bad_edges =
											bad_acc_thingy(bj.P, acc_bad_edges, thread_id, bj.m);
									auto bad = bad_acc_thingy(bj.P, acc_bad, thread_id, bj.m);
									auto qcnt = bad_acc_thingy(bj.P, acc_qcnt, thread_id, bj.m);
									auto conflict_score = bad_acc_thingy(
											bj.P, acc_conflict_score, thread_id, bj.num_cols - 1);
									auto cols = bad_acc_thingy(bj.P, acc_cols, thread_id, bj.m);
									auto rng = bad_acc_thingy(bj.P, acc_rng, thread_id, 1);
									auto success =
											bad_acc_thingy(bj.P, acc_success, thread_id, 1);

									// scuffed queue
									size_t qs = 0, qt = 0;
									auto push = [&](size_t x) {
										acc_bad_edges[P * (qt++) + t] = x;
										if (qt == m)
											qt = 0;
									};
									auto pop = [&]() {
										size_t ret = acc_bad_edges[P * (qs++) + t];
										if (qs == m)
											qs = 0;
										return ret;
									};
									auto empty = [&]() { return qs == qt; };
									auto size = [&]() { return ((qt - qs + m) % m); };

									// only reset the search if unreset_max_iters is hit
									// (it's a reset loop)
									if (unreset_repeat_no == 0) {
										/*
										sycl::distribute_items_and_wait(
												group, [&](sycl::s_item<1> idx) {
													for (size_t i = idx.get_local_id(group, 0); i < m;
															 i += group_size)
														acc_cols[P * i + t] = main_cols[i];
												});
												*/
										for (size_t i = 0; i < m; ++i)
											acc_cols[P * i + t] = main_cols[i];
									}

									// init bad colour (was previously a function)
									auto init_bad_colour = [&]() {
										for (size_t i = 0; i < m; ++i) {
											acc_bad[P * i + t] = 0;
											acc_qcnt[P * i + t] = 0;
										}
										while (!empty())
											pop();
										size_t bad_colour = acc_rng[t]() % num_cols;
										// relabel bad to num_cols - 1 (i.e. swap colours)
										int num_relabeled = 0;
										for (size_t i = 0; i < m; i++) {
											if (cols[i] == bad_colour) {
												push(i);
												acc_bad[P * i + t] = 1;
											} else if (cols[i] == num_cols - 1) {
												acc_cols[P * i + t] = bad_colour;
												++num_relabeled;
											}
										}
									};
									int iter = 0;
									// run for a certain number of iterations before returning to
									// cpu
									while (iter < max_iters && !acc_any_success[0]) {
										// only reset the search if unreset_max_iters is hit
										// (or if things failed)
										if (iter > 0 || unreset_repeat_no == 0)
											init_bad_colour();
										iter++;
										// terminate if anything finished (including self), or if
										// max iters reached
										while (!empty() && iter < max_iters &&
													 !acc_any_success[0]) {
											++iter;
											int e = pop();
											acc_qcnt[P * e + t]++;
											// if stuck (one edge repeated lots),
											// restart from a different colour
											if (qcnt[e] > 1000) {
												// set back to most recent best
												for (size_t i = 0; i < m; ++i)
													acc_cols[P * i + t] = main_cols[i];
												num_cols = num_cols_main;
												init_bad_colour();
												continue;
											}
											for (size_t i = 0; i < conflict_score.size(); ++i)
												acc_conflict_score[P * i + t] = 0;
											for (int f = edge_crossings[e]._Find_first(); f < M;
													 f = edge_crossings[e]._Find_next(f)) {
												if (bad[f] == 0) {
													acc_conflict_score[P * cols[f] + t] +=
															1 + qcnt[f] * qcnt[f];
												}
											}
											// get the minimum conflicting colour
											int c = 0;
											for (int i = 0; i < conflict_score.size(); ++i)
												if (conflict_score[i] < conflict_score[c])
													c = i;
											acc_bad[P * e + t] = 0;
											acc_cols[P * e + t] = c;
											for (int f = edge_crossings[e]._Find_first(); f < M;
													 f = edge_crossings[e]._Find_next(f)) {
												if (bad[f] == 0 && cols[f] == c) {
													acc_bad[P * f + t] = 1;
													push(f);
												}
											}
										}
										if (empty()) {
											// successful improvement!
											// save new best of this thread
											for (size_t e = 0; e < m; ++e)
												acc_main_cols[P * e + t] = acc_cols[P * e + t];
											// increment success counter to later take argmax
											++acc_success[t];
											// recompute num cols
											num_cols = 0;
											for (size_t e = 0; e < m; ++e)
												if (cols[e] > num_cols)
													num_cols = cols[e];
											++num_cols;
											num_cols_main = num_cols;
											// // DEBUG TODO
											// // compute frequency stats using acc_cols temporarily
											// for (size_t c = 0; c < num_cols; ++c)
											// 	acc_cols[P * c + t] = 0;
											// for (size_t e = 0; e < m; ++e)
											// 	acc_cols[P * cols[e] + t]++;
											// printf("col freqs: ");
											// for (size_t c = 0; c < num_cols; ++c)
											// 	printf("%d,", int(acc_cols[P * c + t]));
											// printf("\n");
											// // return data to  acc_cols
											// for (size_t e = 0; e < m; ++e)
											// 	acc_cols[P * e + t] = acc_main_cols[P * e + t];
											// printf("success, cols=%d\n", int(num_cols));
											acc_any_success[0] = 1;
											// TODO: Things break if you go again for some reason
											break;
										}
									}
								});
					});
					q.wait_and_throw();
				}
				q.wait_and_throw();

				{
					// get only the best performing thread
					size_t t =
							max_element(begin(bj.batch_success), end(bj.batch_success)) -
							begin(bj.batch_success);
					cout << "(pre-reset counter: " << unreset_repeat_no << "/"
							 << max_iters_before_reset << " ) success level (thread " << t
							 << "): " << bj.batch_success[t] << endl;
					if (bj.batch_success[t] > 0) {
						cout << "thread: " << t << endl;
						cout << "success!" << endl;
						for (size_t e = 0; e < bj.m; ++e)
							s_main.cols[e] = bj.main_cols[e * bj.P + t];
						size_t old_num_cols = s_main.num_cols;
						s_main.recompute_num_cols();
						cerr << "colours: " << old_num_cols << "->" << s_main.num_cols
								 << endl;
						assert(s_main.verify());
						s_main.save_if_better_than_default();
					}
				}
			}
		} while (continuous);
		std::cout << "gpu optimizer finished" << std::endl;
	}
};
#endif
