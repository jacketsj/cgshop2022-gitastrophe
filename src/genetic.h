#pragma once

#include "encodings.h"
#include "conflict_optimizer.h"
#include "tabucol.h"

struct genetic_improver {
  const int POP_SIZE = 20;
  const int NUM_PARENTS = 4;
  instance ins;
  vector<pair<int, solution>> pop; // pair {fitness, soln}
  vector<int> deg;
  mt19937 rng;
  int k;
  genetic_improver(string filename): ins(filename) {
    pop.resize(POP_SIZE, pair(0, solution(ins)));
    for (auto& [f, s] : pop) {
      s.from_json_default();
    }
    deg.resize(ins.m);
    for (size_t e = 0; e < ins.m; e++) {
      deg[e] = ins.edge_crossings[e].count();
    }
    rng = mt19937(time(0));
    run();
  }
  genetic_improver() {
    char c;
    string file; cin >> file;
    ins = instance(file);
    /*
    ins.n = 0;
    ins.id = "asdf";
    ins.m = 500;
    ins.edge_crossings.resize(ins.m);
    while (cin >> c) {
      if (c == 'c' || c == 'p') {
        string s; getline(cin, s);
        continue;
      }
      int u, v; cin >> u >> v; u--; v--;
      ins.edge_crossings[u][v] = ins.edge_crossings[v][u] = 1;
    }*/
    solution s(ins);
    //s.compute_greedy_update_sorted();
    s.from_json_default();
    while (1) {
      s.num_cols--;
      cerr << "try " << s.num_cols << endl;
      for (size_t e = 0; e < ins.m; e++) {
        if (s.cols[e] == s.num_cols) {
          s.cols[e] = rng() % s.num_cols;
        }
      }
      while (1) {
        tabucol_improver t(s, 1e5, rng);
        if (t.best_confs == 0) break;
      }
    }
    /*
    pop.resize(POP_SIZE, solution(ins));
    for (auto& s : pop) {
      s.compute_greedy_update_sorted();
      //s.from_json_default();
    }
    deg.resize(ins.m);
    for (size_t e = 0; e < ins.m; e++) {
      deg[e] = ins.edge_crossings[e].count();
    }
    rng = mt19937(time(0));
    run();*/
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
    int64_t ans = 0;
    for (int e = s._Find_first(); e < M; e = s._Find_next(e)) {
      ans += ins.m*(ins.edge_crossings[e] & s).count();
    }
    ans -= s.count();
    return ans;
    /*
    int64_t ans = -2LL*ins.m*ins.m*s.count();
    for (int e = s._Find_first(); e < M; e = s._Find_next(e)) {
      ans += 2LL*ins.m*ins.m*ins.m*(ins.edge_crossings[e] & s).count();
      ans -= deg[e];
    }
    return ans;*/
  }
  pair<int, vector<int>> hungarian(const vector<vector<int>> &a) {
    if (a.empty()) return {0, {}};
    int n = size(a) + 1, m = size(a[0]) + 1;
    vector<int> u(n), v(m), p(m), ans(n - 1);
    for(int i=1; i<n; i++) {
      p[0] = i;
      int j0 = 0; // add "dummy" worker 0
      vector<int> dist(m, INT_MAX), pre(m, -1);
      vector<bool> done(m + 1);
      do { // dijkstra
        done[j0] = true;
        int i0 = p[j0], j1, delta = INT_MAX;
        for(int j=1;j < m; j++) if (!done[j]) {
          auto cur = a[i0 - 1][j - 1] - u[i0] - v[j];
          if (cur < dist[j]) dist[j] = cur, pre[j] = j0;
          if (dist[j] < delta) delta = dist[j], j1 = j;
        }
        for(int j=0;j<m;j++) {
          if (done[j]) u[p[j]] += delta, v[j] -= delta;
          else dist[j] -= delta;
        }
        j0 = j1;
      } while (p[j0]);
      while (j0) { // update alternating path
        int j1 = pre[j0];
        p[j0] = p[j1], j0 = j1;
      }
    }
    for(int j=1;j<m;j++) if (p[j]) ans[p[j] - 1] = j - 1;
    return {-v[0], ans}; // min cost
  }
  int64_t dst(const solution& a, const solution& b) {
    vector c(k, vector(k, M));
    vector<bitset<M>> ac(k), bc(k);
    for (size_t e = 0; e < ins.m; e++) {
      ac[a.cols[e]][e] = 1;
      bc[b.cols[e]][e] = 1;
    }
    for (int i = 0; i < k; i++) {
      for (int j = 0; j < k; j++) {
        c[i][j] = M-(ac[i] & bc[j]).count();
      }
    }
    int64_t ret = ins.m + hungarian(c).first - k*M;
    //cerr << ret << " " << ins.m << endl;
    return ret;
  }
  bool similar(solution& a, solution& b) {
    return dst(a, b) < ins.m / 10;
  }
  void run() {
    while (1) {
      run_single();
    }
  }
  void run_single() {
    // try to eliminate a colour
    k = pop[0].second.num_cols;
    cerr << "run_single " << k << endl;
    for (int i = 0; i < POP_SIZE; i++) {
      conflict_improver_2 c(pop[i].second, 1'000, rng);
      pop[i].second.num_cols = k-1;
      while (!c.bad_edges.empty()) {
        size_t e = c.bad_edges.front(); c.bad_edges.pop();
        pop[i].second.cols[e] = rng()%pop[i].second.num_cols;
      }
      pop[i].first = get_num_confs(pop[i].second);
    }
    int gen = 0;
    while (1) {
      gen++;
      // crossover
      vector<int> ps;
      {
        set<int> pars;
        for (int i = 0; i < NUM_PARENTS; i++) {
          int p = rng() % POP_SIZE;
          while (pars.count(p)) p = rng() % POP_SIZE;
          pars.insert(p);
        }
        ps = vector<int>(begin(pars), end(pars));
      }
      vector cols(NUM_PARENTS, vector(k-1, bitset<M>()));
      for (int i = 0; i < NUM_PARENTS; i++) {
        for (size_t e = 0; e < ins.m; e++) {
          cols[i][pop[ps[i]].second.cols[e]][e] = 1;
        }
      }
      bitset<M> need;
      need.set();
      solution chld(ins);
      for (size_t e = 0; e < ins.m; e++) {
        chld.cols[e] = k-2;
      }
      for (int t = 0; t < k-1; t++) {
        // add the best class
        int64_t best = LLONG_MAX;
        int besti = -1, bestc = -1;
        for (int i = 0; i < NUM_PARENTS; i++) {
          for (int c = 0; c < k-1; c++) {
            auto sub = cols[i][c] & need;
            int64_t sc = score(sub);
            if (sc < best) {
              best = sc;
              besti = i;
              bestc = c;
            }
          }
        }
        auto add = cols[besti][bestc] & need;
        need &= ~add;
        for (int e = add._Find_first(); e < M; e = add._Find_next(e)) {
          //assert((ins.edge_crossings[e] & add).count() == 0);
          chld.cols[e] = t;
        }
      }
      cerr << ins.m - (M - need.count()) << " unassigned in child\n";
      chld.num_cols = k-1;
      tabucol_improver c(chld, 1e5, rng);
      //conflict_improver_2 c(chld, 5'000, rng);
      if (c.best_confs == 0) {
        // hooray
        cerr << "success " << k << " -> " << chld.num_cols << endl;
        chld.save_if_better_than_default();
        for (auto& [f, s] : pop) {
          s = solution(chld);
        }
        return;
      }
      int chld_confs = get_num_confs(chld);
      // diversity check
      bool bad = 0;
      for (auto& [f, s] : pop) {
        if (similar(s, chld)) {
          bad = 1;
          if (chld_confs <= f) {
            swap(chld, s);
            swap(chld_confs, f);
          }
          break;
        }
      }
      cerr << "generation " << gen << " avg fitness ";
      cerr << 1.L * accumulate(
        begin(pop),
        end(pop),
        0LL,
        [](int64_t x, const auto& a) { return x + a.first; })
      / pop.size() << endl;
      if (bad) {
        cerr << "failed similarity check\n";
        continue;
      }
      // add the child to population and kick out a thing
      pop.emplace_back(chld_confs, chld);
      sort(begin(pop), end(pop), [](const auto& a, const auto& b) { return a.first < b.first; });
      int c1 = rng() % pop.size();
      while (c1 < pop.size()/2 && rng()%2) {
        c1 = rng() % pop.size();
      }
      int minDist = 1e9;
      int c2 = 0;
      for (size_t i = 0; i < pop.size(); i++) {
        if (i < pop.size()/2 && rng()%2) continue;
        if (i == c1) continue;
        int d = dst(pop[i].second, pop[c1].second);
        if (d < minDist) {
          minDist = d;
          c2 = i;
        }
      }
      pop.erase(begin(pop)+max(c1, c2));
    }
  }
};
