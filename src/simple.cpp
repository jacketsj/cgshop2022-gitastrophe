#include <bits/stdc++.h>

#include "save_folder.h"
#include "simple.h"
#include "simple_conflict.h"
#include "simple_conflict_gpu.h"
#include "simple_recurse.h"
#include "simple_search.h"
#include "simple_search_gpu.h"
#include "simple_tabulate.h"
#include "simple_genetic.h"
#include "head_runner.h"
#include "matching_init.h"
//#include "simple_threaded.h"

using namespace std;

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
	bool optimize = false;
	bool conflict = false;
	int badify = 100;
	bool use_gpu = false;
	size_t max_iters = 100;
	bool tabulate = false;
	bool recurse = false;
	bool from_scratch = false;
	bool forever = false;
  bool genetic = false;
  bool test = false;
  bool head = false;
  bool matching = false;
	optional<string> folder_to_save;

	for (int i = 0; i < argc; ++i) {
		if (string(argv[i]) == string("-t")) {
			threaded = true;
			num_threads = stoi(string(argv[++i]));
		} else if (string(argv[i]) == string("-o")) {
			optimize = true;
		} else if (string(argv[i]) == string("-b")) {
			badify = stoi(string(argv[++i]));
		} else if (string(argv[i]) == string("--gpu")) {
			use_gpu = true;
		} else if (string(argv[i]) == string("--forever")) {
			forever = true;
		} else if (string(argv[i]) == string("--scratch")) {
			from_scratch = true;
		} else if (string(argv[i]) == string("-i")) {
			max_iters = stoi(string(argv[++i]));
		} else if (string(argv[i]) == string("-c")) {
			conflict = true;
		} else if (string(argv[i]) == string("-sf")) {
			folder_to_save = string(argv[++i]);
		} else if (string(argv[i]) == string("--table")) {
			tabulate = true;
		} else if (string(argv[i]) == string("-r")) {
			recurse = true;
		} else if (string(argv[i]) == string("-g")) {
      genetic = true;
    } else if (string(argv[i]) == string("--test")) {
      test = true;
    } else if (string(argv[i]) == string("-h")) {
      head = true;
    } else if (string(argv[i]) == string("-m")) {
      matching = true;
    }
	}

	if (use_gpu) {
		// optimize only for now (still not yet working)
		if (conflict) {
			for (auto& filename : files) {
				simple_conflict_gpu(filename, max_iters, from_scratch, forever,
														(threaded ? num_threads : 0));
			}
		} else {
			for (auto& filename : files) {
				simple_search_gpu(filename, badify);
			}
		}
	} else if (threaded) {
		if (folder_to_save.has_value()) {
			use_threads(files, num_threads, [&](string filename) {
				save_folder(filename, folder_to_save.value());
			});
		} else if (conflict) {
			use_threads(files, num_threads,
									[&](string filename) { simple_conflict(filename); });
		} else if (optimize) {
			use_threads(files, num_threads,
									[&](string filename) { simple_search(filename, badify); });
		} else if (recurse) {
			use_threads(files, num_threads,
									[&](string filename) { simple_recurse(filename); });
		} else if (genetic) {
			use_threads(files, num_threads,
									[&](string filename) { simple_genetic(filename); });
    } else if (head) {
			use_threads(files, num_threads,
									[&](string filename) { run_head(filename); });
    } else if (matching) {
			use_threads(files, num_threads,
									[&](string filename) { compute_matching_init(filename); });
    } else {
			use_threads(files, num_threads,
									[&](string filename) { simple(filename); });
		}
	} else {
		if (folder_to_save.has_value()) {
			for (auto& filename : files)
				save_folder(filename, folder_to_save.value());
		} else if (conflict) {
			for (auto& filename : files) {
				simple_conflict(filename);
			}
		} else if (optimize) {
			for (auto& filename : files) {
				simple_search(filename, badify);
			}
		} else if (tabulate) {
			simple_tabulate(files);
		} else if (recurse) {
			for (auto& filename : files) {
				simple_recurse(filename);
			}
		} else if (genetic) {
			for (auto& filename : files) {
				simple_genetic(filename);
			}
    } else if (head) {
      for (auto& filename : files) {
        run_head(filename);
      }
    } else if (test) {
      auto g = genetic_improver();
    } else if (matching) {
      for (auto& filename : files) {
        compute_matching_init(filename);
      }
    } else {
			for (auto& filename : files) {
				simple(filename);
			}
		}
	}
}
