#include <bits/stdc++.h>

#include "genetic.h"

using namespace std;

void simple_genetic(string file) {
	{
		instance ins(file, false);
		if (ins.m > M) return;
	}
  genetic_improver g(file);
}
