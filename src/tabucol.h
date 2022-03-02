#pragma once

#include "encodings.h"

struct tabucol_improver {
  instance& ins;
  solution& s;
  mt19937& rng;
  //const int TABU_SIZE = 10;
  //const int NUM_NEIGHS = 2e3;
  vector<vector<int>> tabu; // first iter this move becomes non-tabu
  vector<vector<double>> confs; // # conflicts between i and colour class j
  vector<size_t> cols;
  int max_iters;
  double num_confs;
  double best_confs;
  double score(int e, int f) {
    return 1;// - 1. / (1LL*ins.m*ins.m*ins.deg[e]) - 1. / (1LL*ins.m*ins.m*ins.deg[f]);
  }
  tabucol_improver(solution& _s, int _max_iters, mt19937& _rng): ins(_s.ins), s(_s), rng(_rng), max_iters(_max_iters) {
    tabu.resize(ins.m, vector(s.num_cols, 0));
    confs.resize(ins.m, vector(s.num_cols, 0.));
    num_confs = 0;
    for (size_t e = 0; e < ins.m; e++) {
      for (size_t f = ins.edge_crossings[e]._Find_first(); f < M; f = ins.edge_crossings[e]._Find_next(f)) {
        if (s.cols[e] == s.cols[f]) {
          num_confs += score(e, f);
        }
        confs[e][s.cols[f]] += score(e, f);
      }
    }
    num_confs /= 2;
    best_confs = num_confs;
    cols = s.cols;
    run();
  }
  void run() {
    cerr << "running tabu: " << fixed << setprecision(10) << num_confs << " conflicts\n";
    int stale = 0; // # of iters with same conf score
    for (int iter = 0; num_confs > 1e-4 && iter < max_iters; iter++) {
      if (iter % 10000 == 0) {
        cerr << "running iter " << iter << " w/ " << num_confs << " conflicts\n";
      }
      // generate moves
      double best = 1e18;
      vector<pair<int, int>> bestis;
      int num_bad = 0;
      for (size_t e = 0; e < ins.m; e++) {
        if (confs[e][s.cols[e]] < 1e-4) continue;
        num_bad++;
        for (int c = 0; c < s.num_cols; c++) {
          if (c == s.cols[e]) continue;
          double cur = num_confs + confs[e][c] - confs[e][s.cols[e]];
          if (iter < tabu[e][c] && cur > best_confs) {
            // skip tabu move
            //i--;
            continue;
          }
          if (cur < best) {
            best = cur;
            bestis = {{e, c}};
            //if (best < num_confs) goto make_move;
          } else if (cur == best) {
            bestis.emplace_back(e, c);
          }
        }
      }
      make_move:
      // move
      auto [e, bestc] = bestis[rng() % bestis.size()];
      for (int f = ins.edge_crossings[e]._Find_first(); f < M; f = ins.edge_crossings[e]._Find_next(f)) {
        confs[f][s.cols[e]] -= score(e, f);
        confs[f][bestc] += score(e, f);
      }
      tabu[e][s.cols[e]] = iter + 0.6 * num_bad + (rng() % 10 + 1) + stale / 1000;
      s.cols[e] = bestc;
      if (best == num_confs) stale++;
      else stale = 0;
      num_confs = best;
      if (num_confs < best_confs) {
        best_confs = num_confs;
        cols = s.cols;
      }
      if (num_confs < 1e-4) {
        cerr << "tabu success in " << iter << " iterations\n";
        s.recompute_num_cols();
        break;
      }
    }
    s.recompute_num_cols();
    cerr << "finished tabu: " << best_confs << " remaining\n";
    s.cols = cols;
  }
};

