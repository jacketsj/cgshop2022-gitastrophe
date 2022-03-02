#include <bits/stdc++.h>

#include "conflict_optimizer_gpu.h"
//#include "conflict_optimizer_gpu_comments.h"

using namespace std;

void simple_conflict_gpu(string filename, size_t max_iters, bool from_scratch,
												 bool continuous, size_t num_threads) {
#if __has_include(<CL/sycl.hpp>)
	cout << "running..." << endl;
	conflict_improver_gpu c(filename, max_iters, from_scratch, continuous,
													num_threads);
#else
	cout << "cannot start gpu conflict improver: gpu not available (do you have "
					"(hip)SYCL installed?)"
			 << endl;
#endif
}
