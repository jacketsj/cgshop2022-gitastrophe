#include <bits/stdc++.h>

#include "local_search_gpu.h"

using namespace std;

void simple_search_gpu(string filename, int badify) {
#if __has_include(<CL/sycl.hpp>)
	local_search_improver_gpu l(filename, badify);
#else
	cout << "gpu not available (do you have (hip)SYCL installed?)" << endl;
#endif
}
