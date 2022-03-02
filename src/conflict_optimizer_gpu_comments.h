#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

//#include "encodings.h"

#if __has_include(<CL/sycl.hpp>)
#include <CL/sycl.hpp>
using namespace cl;
using namespace std;

sycl::queue q{sycl::gpu_selector{}};

template <typename T>
using gpu_alloc = sycl::usm_allocator<T, sycl::usm::alloc::shared>;

struct conflict_improver_gpu {
	conflict_improver_gpu(string) { run_gpu_bad(); }
	void run_gpu_bad() {
		q = sycl::queue{sycl::gpu_selector{}};
		q.submit([](sycl::handler& cgh) {
			 cgh.parallel_for(sycl::range<1>{10}, [](sycl::id<1> id) {
				 // int i = 0;
			 });
		 }).wait();
	}
};
#endif
