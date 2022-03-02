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

struct local_search_improver {
  instance ins;
  solution s;
  mt19937 rng;
  queue<size_t> bad_edges;
  vector<char> bad;
  int num_badify;
  local_search_improver(string filename, int badify = 100): ins(filename), s(ins), num_badify(badify) {
    s.from_json_default();
    rng = mt19937(time(0));
    bad.resize(ins.m);
    init_bad_colour();
    run();
  }
  void init_bad_colour() {
    fill(begin(bad), end(bad), 0);
    while (!bad_edges.empty()) bad_edges.pop();
    size_t bad_colour = rng() % s.num_cols;
    // relabel bad to s.num_cols - 1
    for (size_t i = 0; i < ins.m; i++) {
      if (s.cols[i] == bad_colour) {
        bad_edges.push(i);
        bad[i] = 1;
      } else if (s.cols[i] == s.num_cols-1) {
        s.cols[i] = bad_colour;
      }
    }
  }
  void run() {
    int iter = 0;
    int lastimprove = 0;
    while (1) {
      if (iter++ % 1000 == 0) {
        cerr << "Running iteration " << iter << ": " << bad_edges.size() << " remaining.\n";
      }
      // if stuck, restart from a different colour
      if (iter - lastimprove > 5e4) {
        s.from_json_default();
        init_bad_colour();
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
      for (int i = 0; i < sz; i++) {
        assert(!nbad_edges.empty());
        int e = nbad_edges.front(); nbad_edges.pop();
        //cerr << "try " << e << " " << ins.m << endl;
        vector<bool> available_cols(s.num_cols-1, true);
        for (int f = ins.edge_crossings[e]._Find_first(); f < M;
             f = ins.edge_crossings[e]._Find_next(f)) {
          if (!nbad[f]) {
            available_cols[ncols[f]] = false;
          }
        }
        vector<int> avail;
        for (size_t c = 0; c < s.num_cols-1; c++) {
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
      }
      // good if number of bad edges decreases
      if (nbad_edges.size() <= bad_edges.size()) {
        if (nbad_edges.size() < bad_edges.size()) {
          lastimprove = iter;
        }
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

