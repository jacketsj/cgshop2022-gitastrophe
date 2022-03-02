#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "encodings.h"

#include "verify/bentley_ottmann_any_intersection.hpp"
#include "verify/geom.hpp"
#include "verify/verify_coloring.hpp"

#include "recursive_init.h"
#include "conflict_optimizer.h"

void compute_matching_init(string file) {
  instance ins(file);
  solution ret(ins);
  // compute the matching graph
  // mapping edge -> pt pair
  vector<vector<pt>> ps(ins.m);
  for (size_t e = 0; e < ins.m; e++) {
    if (ps[e].size() >= 2) continue;
    vector<double> ls;
    for (size_t f = ins.edge_crossings[e]._Find_first(); f < ins.m; f = ins.edge_crossings[e]._Find_next(f)) {
      auto [a, b] = ins.edge_to_ls(ins.edges[e]);
      auto [c, d] = ins.edge_to_ls(ins.edges[f]);
      double lambda = 1.0 * cross(c-a, d-c) / cross(b-a, d-c);
      ls.push_back(lambda);
    }
    sort(begin(ls), end(ls));
    // find best 2 "endpoints"
    int best = -1, besti = -1, best2 = -1, besti2 = -1;
    // want [besti, besti+best), [besti2, besti2+best2) to be disjoint
    int rptr = 0;
    for (int i = 0; i < ls.size(); i++) {
      while (rptr < ls.size() && ls[rptr]-ls[i] < 0.0001) rptr++;
      if (rptr-i > best) {
        if (besti+best <= i) {
          best2 = best;
          besti2 = besti;
        }
        best = rptr-i;
        besti = i;
      } else if (rptr-i > best2 && besti+best <= i) {
        best2 = rptr-i;
        besti2 = i;
      }
    }
    if (best2 < 5) continue;
    //cerr << best << " " << best2 << endl;
    {
      auto [a, b] = ins.edge_to_ls(ins.edges[e]);
      bool bad = 0;
      pt isect = a + (b-a) * (ls[besti] + 0.00005);
      for (pt p : ps[e]) {
        if ((p-isect).length2() < 200) bad = 1;
      }
      if (!bad) ps[e].push_back(isect);
    }
    {
      auto [a, b] = ins.edge_to_ls(ins.edges[e]);
      bool bad = 0;
      pt isect = a + (b-a) * (ls[besti2] + 0.00005);
      for (pt p : ps[e]) {
        if ((p-isect).length2() < 200) bad = 1;
      }
      if (!bad) ps[e].push_back(isect);
    }
  }
  vector<pt> mp;
  vector<vector<pair<int, size_t>>> adj;
  vector<pair<int, int>> es(ins.m, pair(-1, -1));
  int ocnt = 0, bigcnt = 0, gcnt = 0;
  for (size_t e = 0; e < ins.m; e++) {
    if (ps[e].size() == 1) ocnt++;
    if (ps[e].size() == 2) gcnt++;
    if (ps[e].size() > 2) bigcnt++;
    if (ps[e].size() != 2) continue;
    int u, v;
    for (u = 0; u < mp.size(); u++) {
      if ((mp[u]-ps[e][0]).length2() <= 2500) break;
    }
    if (u >= mp.size()) mp.push_back(ps[e][0]);
    for (v = 0; v < mp.size(); v++) {
      if ((mp[v]-ps[e][1]).length2() <= 2500) break;
    }
    if (v >= mp.size()) mp.push_back(ps[e][1]);
    while (adj.size() <= max(u, v)) adj.emplace_back();
    adj[u].emplace_back(v, e);
    adj[v].emplace_back(u, e);
    es[e] = {u, v};
  }
  cerr << "1-endpoint edges: " << ocnt << endl;
  cerr << "2-endpoint edges: " << gcnt << endl;
  cerr << "3+-endpoint edges: " << bigcnt << endl;
  cerr << "graph size: " << adj.size() << endl;
  int lb = 0;
  int ub = ins.m;
  for (auto& v : adj) {
    lb = max(lb, (int)v.size());
    ub = min(ub, (int)v.size());
  }
  cerr << "lower bound: " << lb << endl;
  cerr << "min degree: " << ub << endl;
  // make some stupid 2Delta edge colouring
  vector<int> ord;
  for (size_t e = 0; e < ins.m; e++) {
    if (es[e].first != -1) ord.push_back(e);
  }
  /*
  vector<size_t> aaa(begin(ord), end(ord));
  auto sub = subinstance(ins, aaa);
  solution subs(sub);
  subs.compute_greedy_update_sorted();
  conflict_improver_2 conf(subs, 1000000, rng);*/
  sort(begin(ord), end(ord), [&](int i, int j) {
    return adj[es[i].first].size() + adj[es[i].second].size() > adj[es[j].first].size() + adj[es[j].second].size();
  });
  fill(begin(ret.cols), end(ret.cols), ins.m);
  int mxc = 0;
  for (int e : ord) {
    vector<char> bad(1000);
    auto [u, v] = es[e];
    for (auto [_, f] : adj[u]) {
      if (ret.cols[f] < ins.m) bad[ret.cols[f]] = 1;
    }
    for (auto [_, f] : adj[v]) {
      if (ret.cols[f] < ins.m) bad[ret.cols[f]] = 1;
    }
    /*
    for (size_t f = ins.edge_crossings[e]._Find_first(); f < ins.m; f = ins.edge_crossings[e]._Find_next(f)) {
      if (ret.cols[f] < ins.m) bad[ret.cols[f]] = 1;
    }*/
    for (int c = 0; c < 1000; c++) {
      if (!bad[c]) {
        ret.cols[e] = c;
        mxc = max(mxc, c);
        break;
      }
    }
  }
  cerr << "edge coloured w/ " << mxc+1 << " colours\n";
  // find a colouring with not many conflicts
  fill(begin(ret.cols), end(ret.cols), ins.m);
  int goodcnt = 0;
  vector<int> nord;
  for (int e : ord) {
    vector<int> bad(mxc+1);
    //auto [u, v] = es[e];
    /*
    for (auto [_, f] : adj[u]) {
      if (ret.cols[f] < ins.m) bad[ret.cols[f]]++;
    }
    for (auto [_, f] : adj[v]) {
      if (ret.cols[f] < ins.m) bad[ret.cols[f]]++;
    }*/
    for (size_t f = ins.edge_crossings[e]._Find_first(); f < ins.m; f = ins.edge_crossings[e]._Find_next(f)) {
      if (ret.cols[f] < ins.m) bad[ret.cols[f]] = 1;
    }
    int mn = *min_element(begin(bad), end(bad));
    if (0 && mn > 10) {
      //ret.recompute_num_cols();
      /*
      ret.num_cols = mxc+21;
      ret.cols[e] = ret.num_cols-1;
      size_t cur = ret.num_cols;
      conflict_improver_2 c(ret, 20, rng);
      if (ret.num_cols == cur) {
        int c = rng() % (mxc+20);
        while (bad[c] > mn) c = rng() % (mxc+5);
        ret.cols[e] = c;
        //nord.push_back(e);
      } else {
        goodcnt++;
      }*/
    } else {
      int c = rng() % (mxc+1);
      while (bad[c] > mn) c = rng() % (mxc+1);
      ret.cols[e] = c;
      goodcnt++;
    }
  }
  cerr << "succesfully edge coloured " << goodcnt << " edges\n";
  ord = nord;
  for (size_t e = 0; e < ins.m; e++) {
    if (es[e].first == -1) ord.push_back(e);
  }
  sort(begin(ord), end(ord), [&](int i, int j) {
    return ins.deg[i] > ins.deg[j];
  });
  for (int e : ord) {
    assert(ret.cols[e] == ins.m);
    vector<char> bad(1000);
    for (size_t f = ins.edge_crossings[e]._Find_first(); f < ins.m; f = ins.edge_crossings[e]._Find_next(f)) {
      if (ret.cols[f] < ins.m) bad[ret.cols[f]] = 1;
    }
    for (int c = 0; c < 1000; c++) {
      if (!bad[c]) {
        ret.cols[e] = c;
        break;
      }
    }
  }
  ret.recompute_num_cols();
  cerr << "found conflicting solution w/ " << ret.num_cols << " colours\n";
  for (; ret.num_cols < 170; ret.num_cols++) {
    tabucol_improver t(ret, 2e3, rng);
  }
  //tabucol_improver t(ret, 1e5, rng);
  ord.clear();
  ord.resize(ins.m);
  iota(begin(ord), end(ord), 0);
  shuffle(begin(ord), end(ord), rng);
  int badcnt = 0;
  for (int e = 0; e < ins.m; e++) {
    for (size_t f = ins.edge_crossings[e]._Find_first(); f < ins.m; f = ins.edge_crossings[e]._Find_next(f)) {
      if (ret.cols[e] == ret.cols[f]) {
        ret.cols[e] = ret.num_cols;
        badcnt++;
        break;
      }
    }
  }
  cerr << "num conflicting: " << badcnt << endl;
  ret.recompute_num_cols();
  conflict_improver c(ret);
  //ret.save_if_better_than_default();
}
