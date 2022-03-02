#pragma once

#include "encodings.h"
#include "conflict_optimizer.h"
#include "tabucol.h"

struct head_improver {
  instance ins;
  pair<int, solution> p1, p2, elite1, elite2, best; // pair {fitness, soln}
  mt19937 rng;
  int k;
  head_improver(string filename): ins(filename), p1(0, solution(ins)), p2(0, solution(ins)), elite1(1e9, solution(ins)), elite2(1e9, solution(ins)), best(1e9, solution(ins)) {
    rng = mt19937(time(0));
    run();
  }
  void init(pair<int, solution>& p) {
    p.second.from_json_default();
    k = p.second.num_cols;
    conflict_improver_2 c(p.second, 1'000, rng);
    p.second.num_cols = k-1;
    while (!c.bad_edges.empty()) {
      size_t e = c.bad_edges.front(); c.bad_edges.pop();
      p.second.cols[e] = rng()%p.second.num_cols;
    }
    p.first = get_num_confs(p.second);
  }
  int64_t get_num_confs(const solution& s) {
    vector<bitset<M>> cols(s.num_cols);
    for (size_t e = 0; e < ins.m; e++) {
      cols[s.cols[e]][e] = 1;
    }
    int64_t ans = 0;
    for (int c = 0; c < s.num_cols; c++) {
      for (int e = cols[c]._Find_first(); e < M; e = cols[c]._Find_next(e)) {
        ans += (ins.edge_crossings[e] & cols[c]).count();
      }
    }
    return ans / 2;
  }
  int64_t score(bitset<M>& s) {
    return -s.count();
    /*
    int64_t ans = 0;
    for (int e = s._Find_first(); e < M; e = s._Find_next(e)) {
      ans += ins.m*(ins.edge_crossings[e] & s).count();
    }
    ans -= s.count();
    return ans;*/
  }
  pair<int, solution> crossover(const solution& a, const solution& b) {
    vector colsa(k-1, bitset<M>()), colsb(k-1, bitset<M>());
    for (size_t e = 0; e < ins.m; e++) {
      colsa[a.cols[e]][e] = 1;
      colsb[b.cols[e]][e] = 1;
    }
    bitset<M> need;
    need.set();
    solution chld(ins);
    for (size_t e = 0; e < ins.m; e++) {
      chld.cols[e] = rng()%(k-1);//k-2;
    }
    for (int t = 0; t < k-1; t++) {
      auto& cols = t % 2 ? colsa : colsb;
      // add the best class
      int64_t best = LLONG_MAX;
      int bestc = -1;
      for (int c = 0; c < k-1; c++) {
        auto sub = cols[c] & need;
        int64_t sc = score(sub);
        if (sc < best) {
          best = sc;
          bestc = c;
        }
      }
      auto add = cols[bestc] & need;
      need &= ~add;
      for (int e = add._Find_first(); e < M; e = add._Find_next(e)) {
        //assert((ins.edge_crossings[e] & add).count() == 0);
        chld.cols[e] = t;
      }
    }
    cerr << ins.m - (M - need.count()) << " unassigned in child\n";
    chld.num_cols = k-1;
    tabucol_improver c(chld, 1e5, rng);
    return pair((int) get_num_confs(chld), chld);
  }
  void run() {
    while (1) {
      init(p1);
      init(p2);
      init(elite2);
      elite1.first = 1e9;
      best.first = 1e9;
      run_single();
    }
  }
  const int ITER_CYCLE = 10;
  void run_single() {
    // try to eliminate a colour
    cerr << "run_single " << k << endl;
    int gen = 0, cyc = 0;
    while (1) {
      gen++;
      auto c1 = crossover(p1.second, p2.second);
      auto c2 = crossover(p2.second, p1.second);
      p1 = c1;
      p2 = c2;
      if (p1.first < elite1.first) {
        elite1 = p1;
      }
      if (p2.first < elite1.first) {
        elite1 = p1;
      }
      if (elite1.first < best.first) {
        best = elite1;
      }
      if (best.first == 0) break;
      if (gen % ITER_CYCLE == 0) {
        p1 = elite2;
        elite2 = elite1;
        elite1 = pair(1e9, solution(ins));
        cyc++;
      }
    }
    best.second.recompute_num_cols();
    cerr << "success! " << k << " -> " << best.second.num_cols << endl;
    best.second.save_if_better_than_default();
  }
};
