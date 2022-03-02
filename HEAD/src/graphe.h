#ifndef GRAPHE_H
#define GRAPHE_H


#include<vector>
#include<string>
#include <bitset>

#include "util/gfile.h"

using namespace std;


class Graph {
public:
  static constexpr int M = 25000;
	void buildVoisins(); // cree le vecteur des voisins de chaque sommet

	Graph(){filename=""; nbSommets=0; nbArretes=0; tVoisins=NULL;}
	Graph(string fname){filename=fname; loadGraph();}
	~Graph();
	void loadGraph();
	void loadMatrixGraph(GInputFile& infile);

	int nbSommets;
	int nbArretes;
  vector<bitset<M>> tConnect; // tableau carre des connections noeud à noeud
	vector<int>* tVoisins; // tableau qui pour chaque noeud contient la liste des voisins
	string filename;
};


#include "graphe.cpp"
#endif
