#pragma once

#define SOLUTIONS_FOLDER "temp/"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using std::allocator;
using std::bitset;
using std::ifstream;
using std::pair;
using std::string;
using std::to_string;
using std::vector;

#include <nlohmann/json.hpp>
using json = nlohmann::json;

constexpr int M = 4001;
struct instance {
	string id;
	size_t n, m;
	vector<bitset<M>> edge_crossings;
  vector<int> deg;
  instance() = default;
	instance(string filename) {
		load_dimacs(filename);
	}
	void load_dimacs(const string& fileloc) {
    id = fileloc.substr(fileloc.find("/")+1, fileloc.size());
		ifstream ifs(fileloc);
    string s;
    n = 0;
    while (getline(ifs, s)) {
      std::istringstream in(s);
      in >> s;
      if (s == "e") {
        int u, v; in >> u >> v; u--; v--;
        //m++;
        while (edge_crossings.size() <= std::max(u, v)) edge_crossings.emplace_back();
        edge_crossings[u][v] = edge_crossings[v][u] = 1;
      }
    }
    m = edge_crossings.size();
    deg.resize(m);
    for (int u = 0; u < m; u++) {
      deg[u] = edge_crossings[u].count();
    }
	}
};

struct solution {
	instance& ins;
	vector<size_t> cols;
	size_t num_cols;
	solution(const solution& _oth)
			: ins(_oth.ins), cols(_oth.cols), num_cols(_oth.num_cols) {}
	//= default;
  void operator=(const solution& o) {
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
	bool is_better(const solution& oth) const {
		if (verify())
			if (!oth.verify() || num_cols < oth.num_cols)
				return true;
		return false;
	}
	bool is_better_than_default() const {
		solution oth(ins);
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
	}
	void recompute_num_cols() {
		num_cols = 0;
		for (auto& c : cols)
			num_cols = std::max(c, num_cols);
		++num_cols;
	}
	void compute_greedy_update(const vector<int>& permutation) {
		// assumes crossings are computed for ins
		// TODO; allow permuting this order
		for (size_t i : permutation) {
			// look for smallest available col
			vector<int> available_cols(num_cols, 1);
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
		vector<int> p;
		for (size_t i = 0; i < ins.m; ++i)
			p.push_back(i);
		compute_greedy_update(p);
	}
	void compute_greedy_update_random_shuffle() {
		vector<int> p;
		for (size_t i = 0; i < ins.m; ++i)
			p.push_back(i);
		std::random_shuffle(p.begin(), p.end());
		compute_greedy_update(p);
	}
  void compute_greedy_update_sorted() {
    vector<int> p;
    for (size_t i = 0; i < ins.m; ++i)
      p.push_back(i);
    std::sort(p.begin(), p.end(), [&](int i, int j) {
      return ins.deg[i] > ins.deg[j];
    });
    compute_greedy_update(p);
  }
	bool verify() const {
    for (int u = 0; u < ins.m; u++) {
      for (int v = u+1; v < ins.m; v++) {
        if (ins.edge_crossings[u][v] && cols[u] == cols[v]) {
          return 0;
        }
      }
    }
    return 1;
	}
	solution(instance& _ins)
			: ins(_ins), cols(ins.m, 0) {
		compute_trivial_sln();
		//compute_greedy_update_random_shuffle();
	}
};
