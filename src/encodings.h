#pragma once

#define SOLUTIONS_FOLDER "temp/"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include "verify/bentley_ottmann_any_intersection.hpp"
#include "verify/geom.hpp"
#include "verify/verify_coloring.hpp"

using std::allocator;
using std::bitset;
using std::ifstream;
using std::pair;
using std::string;
using std::to_string;
// using std::vector;
using json = nlohmann::json;

typedef pair<size_t, size_t> edge;

//#include "altvec.h"

struct pt {
	int64_t x, y;
	pt() : x(0), y(0) {}
	pt(int64_t _x, int64_t _y) : x(_x), y(_y) {}
  bool operator==(const pt& o) const {
    return x == o.x && y == o.y;
  }
	pt operator+(const pt& oth) const { return pt(x + oth.x, y + oth.y); }
	pt operator-(const pt& oth) const { return pt(x - oth.x, y - oth.y); }
	pt operator*(const double& a) const { return pt(x * a, y * a); }
	pt operator/(const double& a) const { return pt(x / a, y / a); }
	int64_t length2() const { return x * x + y * y; }
	friend int64_t dot(const pt& a, const pt& b) { return a.x * b.x + a.y * b.y; }
	friend int64_t cross(const pt& a, const pt& b) { return a.x * b.y - a.y * b.x; }
	friend int winding(const pt& a, const pt& b, const pt& c) {
		int64_t res = cross(c - a, b - a);
		if (res > 0)
			res = 1;
		if (res < 0)
			res = -1;
		return res;
	}
};

namespace std {
template <>
struct hash<pt> {
  size_t operator()(const pt& p) const {
    return (p.x << 32) + p.y;
  }
};
}

struct ls {
	pt a, b;
	ls() {}
	ls(const pt& _a, const pt& _b) : a(_a), b(_b) {}
	// TODO: line intersections currently do not include endpoints or parallel,
	// is this ok?
	bool cross_line(const ls& oth) const {
		// assumes this is a line, oth is a line segment
		return winding(a, b, oth.a) * winding(a, b, oth.b) == -1;
	}
	bool cross(const ls& oth) const {
		// assumes both are line segments
		return cross_line(oth) && oth.cross_line(*this);
	}
};

constexpr int M = 75001;

#if __has_include(<CL/sycl.hpp>)
#include <CL/sycl.hpp>
using namespace cl;
template <typename T>
using gpu_alloc_shared = sycl::usm_allocator<T, sycl::usm::alloc::shared>;
// TODO: Revert to gpu selector
// sycl::queue q{sycl::gpu_selector{}};
sycl::queue q{};
template <template <typename> class Alloc, typename T> Alloc<T> get_alloc() {
	if constexpr (std::is_same<Alloc<T>, gpu_alloc_shared<T>>::value) {
		static gpu_alloc_shared<T>* alloc = NULL;
		if (alloc == NULL)
			alloc = new gpu_alloc_shared<T>(q);
		return *alloc;
	} else {
		static allocator<T>* alloc = NULL;
		if (alloc == NULL)
			alloc = new allocator<T>();
		return *alloc;
	}
}
#else
template <template <typename> class Alloc, typename T> Alloc<T> get_alloc() {
	static Alloc<T>* alloc = NULL;
	if (alloc == NULL)
		alloc = new Alloc<T>();
	return *alloc;
}
#endif

/*
template <typename T> gpu_alloc<T> get_alloc<gpu_alloc, T>() {
	static gpu_alloc<T>* alloc = NULL;
	if (alloc == NULL)
		alloc = new gpu_alloc<T>(q);
	return *alloc;
}

template <typename T> allocator<T> get_alloc<allocator, T>() {
	static allocator<T>* alloc = NULL;
	if (alloc == NULL)
		alloc = new allocator<T>();
	return *alloc;
}
*/

template <template <typename> class Alloc = allocator,
					template <class T, class Allocator> class vector = std::vector>
struct instance_templated {
	string id;
	size_t n, m;
	/*
	Alloc<pt> alloc_pt;
	Alloc<edge> alloc_edge;
	Alloc<uint32_t> alloc_uint32_t;
	Alloc<vector<uint32_t, Alloc<uint32_t>>> alloc_vector_uint32_t;
	Alloc<bitset<M>> alloc_bitset;
	*/
	vector<pt, Alloc<pt>> pts;
	vector<edge, Alloc<edge>> edges;
	vector<vector<uint32_t, Alloc<uint32_t>>,
				 Alloc<vector<uint32_t, Alloc<uint32_t>>>>
			adj;
	vector<bitset<M>, Alloc<bitset<M>>> edge_crossings;
  vector<int, Alloc<int>> deg;
	/*
#if __has_include(<CL/sycl.hpp>)
	template <>
	instance_templated<gpu_alloc>(sycl::queue q)
			: alloc_pt(q), alloc_edge(q), alloc_uint32_t(q), alloc_vector_uint32_t(q),
				alloc_bitset(q), pts(alloc_pt), edges(alloc_edge),
				adj(alloc_vector_uint32_t), edge_crossings(alloc_bitset) {}
	template <>
	instance_templated<gpu_alloc>(sycl::queue q, string filename)
			: alloc_pt(q), alloc_edge(q), alloc_uint32_t(q), alloc_vector_uint32_t(q),
				alloc_bitset(q), pts(alloc_pt), edges(alloc_edge),
				adj(alloc_vector_uint32_t), edge_crossings(alloc_bitset) {
		load_json(filename);
	}
#endif
*/
	// instance_templated() = default;
	instance_templated()
			: pts(get_alloc<Alloc, pt>()), edges(get_alloc<Alloc, edge>()),
				adj(get_alloc<Alloc, vector<uint32_t, Alloc<uint32_t>>>()),
				edge_crossings(get_alloc<Alloc, bitset<M>>()),
        deg(get_alloc<Alloc, int>()) {}
	// template <> instance_templated<allocator>() = default;
	//
	//: pts(alloc_pt), edges(alloc_edge), adj(alloc_vector_uint32_t),
	//	edge_crossings(alloc_bitset) = default;
	//
	// template <>
	// instance_templated<allocator>(string filename)
	//	: pts(alloc_pt), edges(alloc_edge), adj(alloc_vector_uint32_t),
	//		edge_crossings(alloc_bitset) {
	instance_templated(string filename, bool compute_edges=1)
			: pts(get_alloc<Alloc, pt>()), edges(get_alloc<Alloc, edge>()),
				adj(get_alloc<Alloc, vector<uint32_t, Alloc<uint32_t>>>()),
				edge_crossings(get_alloc<Alloc, bitset<M>>()),
        deg(get_alloc<Alloc, int>()) {
		load_json(filename, compute_edges);
	}
	ls edge_to_ls(const edge& e) { return ls{pts[e.first], pts[e.second]}; }
	void compute_adj() {
		// vector<vector<uint32_t>, Alloc<uint32_t>>& v_adj = adj;
		adj = vector<vector<uint32_t, Alloc<uint32_t>>,
								 Alloc<vector<uint32_t, Alloc<uint32_t>>>>(
				n, vector<uint32_t, Alloc<uint32_t>>(get_alloc<Alloc, uint32_t>()),
				get_alloc<Alloc, vector<uint32_t, Alloc<uint32_t>>>());
		for (auto& e : edges) {
			adj[e.first].push_back(e.second);
			adj[e.second].push_back(e.first);
		}
	}
	void load_json(const string& fileloc, bool compute_edges=1) {
		ifstream ifs(fileloc);
		json j = json::parse(ifs);
		id = j["id"];
		n = j["n"];
		m = j["m"];
		pts = vector<pt, Alloc<pt>>(n, pt(), get_alloc<Alloc, pt>());
		for (size_t i = 0; i < n; ++i) {
			pt p;
			pts[i].x = j["x"][i];
			pts[i].y = p.y = j["y"][i];
		}
		edges = vector<edge, Alloc<edge>>(m, edge(), get_alloc<Alloc, edge>());
		for (size_t i = 0; i < m; ++i) {
			edges[i].first = j["edge_i"][i];
			edges[i].second = j["edge_j"][i];
		}
		compute_adj();
		if (compute_edges) {
			compute_crossings();
		}
	}
	// TODO: store these on disk or use the provided bentley ottmann code
	void compute_crossings() {
		edge_crossings = vector<bitset<M>, Alloc<bitset<M>>>(
				m, bitset<M>(), get_alloc<Alloc, bitset<M>>());
		std::cout << "Computing edge crossings" << std::endl;
		// Note: std::allocator used here:
		std::vector<boi::Segment> segments;
		for (size_t i = 0; i < m; ++i) {
			ls ls_i = edge_to_ls(edges[i]);
			segments.emplace_back(boi::Point(ls_i.a.x, ls_i.a.y),
														boi::Point(ls_i.b.x, ls_i.b.y));
		}
		boi::prepare_segments(segments);
		size_t total_crossings = 0;
		size_t total_checked = 0;
		for (size_t i = 0; i < m; ++i) {
			for (size_t j = i + 1; j < m; ++j) {
				if (do_intersect(segments[i], segments[j])) {
					edge_crossings[i][j] = 1;
					edge_crossings[j][i] = 1;
					++total_crossings;
				}
				++total_checked;
			}
		}
		// edge_crossings = altvec2<uint32_t>(std::move(v_edge_crossings));
		// edge_crossings = v_edge_crossings;
		std::cout << "Done computing edge crossings: " << total_crossings << "/"
							<< total_checked << " pairs cross" << std::endl;
    deg = vector<int, Alloc<int>>(m, 0, get_alloc<Alloc, int>());
    for (size_t i = 0; i < m; i++) {
      deg[i] = edge_crossings[i].count();
    }
	}
};
using instance = instance_templated<allocator, std::vector>;

template <template <typename> class Alloc = allocator,
					template <class T, class Allocator> class vector = std::vector>
struct solution_templated {
	instance_templated<Alloc, vector>& ins;
	vector<size_t, Alloc<size_t>> cols;
	size_t num_cols;
	solution_templated(const solution_templated<Alloc, vector>& _oth)
			: ins(_oth.ins), cols(_oth.cols), num_cols(_oth.num_cols) {}
	//= default;
  void operator=(const solution_templated<Alloc, vector>& o) {
    ins = o.ins;
    cols = o.cols;
    num_cols = o.num_cols;
  }
	string to_json() {
		string s;
		s += "{\n";
		s += "\"type\": \"Solution_CGSHOP2022\",\n";
		s += "\"instance\": \"" + ins.id + "\",\n";
		s += "\"num_colors\": " + to_string(num_cols) + ",\n";
		s += "\"colors\": [";
		for (size_t i = 0; i < cols.size(); ++i)
			s += to_string(cols[i]) + (i < cols.size() - 1 ? "," : "");
		s += "]\n";
		s += "}\n";
		return s;
	}
	void from_json(string filename) {
		std::ifstream ifs(filename);
		if (!ifs.is_open()) {
			std::cout << "file " << filename
								<< " does not exist (no recorded solution?)\n";
			return;
		}
		json j = json::parse(ifs);

		num_cols = j["num_colors"];
		for (size_t i = 0; i < cols.size(); ++i)
			cols[i] = j["colors"][i];
	}
	void from_json_folder(string folder) {
		if (folder[folder.size() - 1] != '/')
			folder += "/";
		string filename = string(folder) + ins.id + ".json";
		from_json(filename);
	}
	void from_json_default() { from_json_folder(SOLUTIONS_FOLDER); }
	bool is_better(const solution_templated<Alloc, vector>& oth) const {
		if (verify())
			if (!oth.verify() || num_cols < oth.num_cols)
				return true;
		return false;
	}
	bool is_better_than_default() const {
		solution_templated<Alloc, vector> oth(ins);
		oth.from_json_default();
		return is_better(oth);
	}
	void save_debug() {
		std::cout << "Doing debug save for " << ins.id
							<< ", saving to debug.json...\n";
		string filename = "debug.json";
		std::ofstream of(filename);
		of << to_json();
	}
	void save_if_better_than_default() {
		if (is_better_than_default()) {
			std::cout << "Improvement found for " << ins.id << ", saving...\n";
			string filename = string(SOLUTIONS_FOLDER) + ins.id + ".json";
			std::ofstream of(filename);
			of << to_json();
		}
	}
	void compute_trivial_sln() {
		for (size_t i = 0; i < cols.size(); ++i)
			cols[i] = i;
		num_cols = cols.size();
		;
	}
	void recompute_num_cols() {
		num_cols = 0;
		for (auto& c : cols)
			num_cols = std::max(c, num_cols);
		++num_cols;
	}
	void compute_greedy_update(const vector<int, Alloc<int>>& permutation) {
		// assumes crossings are computed for ins
		// TODO; allow permuting this order
		for (size_t i : permutation) {
			// look for smallest available col
			vector<int, Alloc<int>> available_cols(num_cols, 1,
																						 get_alloc<Alloc, int>());
			for (size_t j = ins.edge_crossings[i]._Find_first(); j < M;
					 j = ins.edge_crossings[i]._Find_next(j))
				available_cols[cols[j]] = 0;
			for (size_t c = 0; c < num_cols; ++c)
				if (available_cols[c]) {
					cols[i] = c;
					break;
				}
		}
		recompute_num_cols();
	}
	void compute_greedy_update() {
		vector<int, Alloc<int>> p(get_alloc<Alloc, int>());
		for (size_t i = 0; i < ins.m; ++i)
			p.push_back(i);
		compute_greedy_update(p);
	}
	void compute_greedy_update_random_shuffle() {
		vector<int, Alloc<int>> p(get_alloc<Alloc, int>());
		for (size_t i = 0; i < ins.m; ++i)
			p.push_back(i);
		std::random_shuffle(p.begin(), p.end());
		compute_greedy_update(p);
	}
  void compute_greedy_update_sorted() {
    vector<int, Alloc<int>> p(get_alloc<Alloc, int>());
    for (size_t i = 0; i < ins.m; ++i)
      p.push_back(i);
    std::sort(p.begin(), p.end(), [&](int i, int j) {
      return ins.edge_crossings[i].count() > ins.edge_crossings[j].count();
    });
    compute_greedy_update(p);
  }
	bool verify() const {
    //return 1;
		// std::vector<size_t> cols_default(cols.begin(), cols.end());
		std::vector<size_t> cols_default;
		for (auto& c : cols)
			cols_default.push_back(c);
		std::vector<boi::Segment> segments;
		for (size_t i = 0; i < ins.m; ++i) {
			ls ls_i = ins.edge_to_ls(ins.edges[i]);
			segments.emplace_back(boi::Point(ls_i.a.x, ls_i.a.y),
														boi::Point(ls_i.b.x, ls_i.b.y));
		}
		boi::prepare_segments(segments);
		boi::ColoringVerifier cv(segments, cols_default);
		auto ret = cv.verify();
		if (ret.has_value()) {
			std::cerr << "bad " << ret.value().message << std::endl;
		}
		return !ret.has_value();
		/*
		// fast code - should work but doesn't
		// because the edge crossings are wrong
		for (size_t i = 0; i < ins.m; ++i)
			for (size_t j : ins.edge_crossings[i])
				if (cols[i] == cols[j])
					return false;
		return true;
		*/
	}
	/*
	void compute_greedy_update_random() {
		vector<int> p;
		for (size_t i = 0; i < ins.m; ++i)
			p.push_back(i);
		random_shuffle(p.begin(), p.end());
		compute_greedy_update(p);
	}
	*/
	solution_templated(instance_templated<Alloc, vector>& _ins)
			: ins(_ins), cols(ins.m, 0, get_alloc<Alloc, size_t>()) {
		compute_trivial_sln();
		//compute_greedy_update_random_shuffle();
	}
};
using solution = solution_templated<allocator, std::vector>;
