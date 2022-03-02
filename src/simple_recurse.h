#include <bits/stdc++.h>

#include "recursive_init.h"

using namespace std;

void simple_recurse(string file) {
  instance ins;
  ins.load_json(file);
  solution s = recur(ins);
  cout << "cols: " << s.num_cols << endl;
  cout << (s.verify() ? "verified" : "bad") << endl;
  s.save_if_better_than_default();
}
