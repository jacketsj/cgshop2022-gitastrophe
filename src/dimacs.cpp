#include <bits/stdc++.h>

#include "dimacs_optimizer.h"

using namespace std;

void simple_conflict(string file) {
  instance ins(file);
  solution s(ins);
  s.compute_greedy_update_sorted();
  s.from_json_default();
  conflict_improver c(s);
}

void use_threads(vector<string> todos, size_t num_threads,
								 function<void(string)> process_file) {
	vector<thread> pool;
	mutex mtx;
	size_t i = 0;
	for (size_t t = 0; t < num_threads; ++t) {
		pool.emplace_back([&]() {
			while (true) {
				string filename;
				{
					lock_guard<mutex> l(mtx);
					if (i >= todos.size())
						break;
					filename = todos[i++];
					cout << "thread " << t << " processing file " << (i + 1) << "/"
							 << todos.size() << ": " << filename << endl;
				}
				process_file(filename);
			}
		});
	}
	for (size_t t = 0; t < num_threads; ++t) {
		pool[t].join();
	}
}

int main(int argc, char* argv[]) {
	string filename;
	vector<string> files;
	while (cin >> filename) {
		files.push_back(filename);
	}

	bool threaded = false;
	size_t num_threads = 1;

	for (int i = 0; i < argc; ++i) {
		if (string(argv[i]) == string("-t")) {
			threaded = true;
			num_threads = stoi(string(argv[++i]));
    }
	}

	if (threaded) {
    use_threads(files, num_threads,
                [&](string filename) { simple_conflict(filename); });
	} else {
    for (auto& filename : files) {
      simple_conflict(filename);
    }
	}
}
