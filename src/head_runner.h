#include "encodings.h"
#include "conflict_optimizer.h"
#include <ctime>
#include <sys/time.h>
#include "../HEAD/src/graphe.h"
#include "../HEAD/src/head.h"

void run_head(string file) {
	{
		instance ins(file, false);
		if (ins.m > 25000) return;
	}
  char graphName[80] = "";
  int nbColor;
  long long nbLocalSearch = 200'000;
  int nbGeneration = -1;
  int nbRun = 1;
  int randSeed = 0;
  int nbRandCross = 5;
  int swapIter = 10;
  double swapingRate = 99;
  bool quite = false;
  bool tabucol = false;
  double tauxAcceptWorst = 1.0;
  double weightParent = 0.5;
  char outputFile[255] = "";
  int maxseconds = -1;


  instance ins(file);
  solution s(ins);
  mt19937 rng(time(0));
  Graph g;
  g.nbSommets = ins.m;
  g.tConnect.resize(ins.m);
  for (int i = 0; i < g.nbSommets; i++) {
    for (int j = 0; j < g.nbSommets; j++) {
      g.tConnect[i][j] = ins.edge_crossings[i][j];
      if (g.tConnect[i][j]) {
        g.nbArretes++;
      }
    }
  }
  g.buildVoisins();
  while (1) {
    s.from_json_default();
    nbColor = s.num_cols-1;
    double totalCpuTime=0;
    double totalHumanTime=0;
    unsigned long long totalIterations=0;
    unsigned long long totalIterationsCross=0;
    unsigned long long totalIterationsCrossWithWrongRun=0;
    int nbFound=0;
    Head* solver=new Head();

    // >> Solver Parameters
    solver->graph=&g;
    solver->nbColors=nbColor;
    solver->nbLocalSearch=nbLocalSearch;
    solver->nbGeneration=nbGeneration;
    solver->userRandSeed=randSeed;
    solver->nbRandCross=nbRandCross;
    solver->swapIter=swapIter;
    solver->swapingRate=swapingRate/100.0;
    solver->debug=!quite;
    solver->tauxAcceptWorst=tauxAcceptWorst;
    solver->weightParent=weightParent/100.0;
    solver->initSol = new Solution(&g, nbColor);
    {
      auto s2 = s;
      /*
      conflict_improver_2 c(s2, 5'000, rng);
      while (!c.bad_edges.empty()) {
        size_t e = c.bad_edges.front(); c.bad_edges.pop();
        s2.cols[e] = rng()%(s.num_cols-1);
      }*/
      vector<int> cnt(s.num_cols);
      for (int i = 0; i < g.nbSommets; i++) {
        cnt[s.cols[i]]++;
      }
      int bad = min_element(begin(cnt), end(cnt)) - begin(cnt);
      for (int i = 0; i < g.nbSommets; i++) {
        if (s2.cols[i] == bad) s2.cols[i] = rng()%(s.num_cols-1);
        else if (s2.cols[i] == s.num_cols-1) s2.cols[i] = bad;
        solver->initSol->tColor[i] = s2.cols[i];
      }
    }
    solver->initSol2 = new Solution(&g, nbColor);
    {
      auto s2 = s;
      /*
      conflict_improver_2 c(s2, 5'000, rng);
      while (!c.bad_edges.empty()) {
        size_t e = c.bad_edges.front(); c.bad_edges.pop();
        s2.cols[e] = rng()%(s.num_cols-1);
      }*/
      vector<int> cnt(s.num_cols);
      for (int i = 0; i < g.nbSommets; i++) {
        cnt[s.cols[i]]++;
      }
      int bad1 = min_element(begin(cnt), end(cnt)) - begin(cnt);
      int bad = rng() % s.num_cols;
      while (bad == bad1) bad = rng() % s.num_cols;
      /*
      int best = 1e9, bad = -1;
      for (int i = 0; i < g.nbSommets; i++) {
        if (i == bad1) continue;
        if (cnt[i] < best) {
          best = cnt[i];
          bad = i;
        }
      }*/
      for (int i = 0; i < g.nbSommets; i++) {
        if (s2.cols[i] == bad) s2.cols[i] = rng()%(s.num_cols-1);
        else if (s2.cols[i] == s.num_cols-1) s2.cols[i] = bad;
        solver->initSol2->tColor[i] = s2.cols[i];
      }
    }
    solver->maxsecondes=maxseconds;
    solver->tabucol = tabucol;
    // << solver parameters

    //// >>>>>>>>>  affichage de l'heure au debut
    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    if(!quite)
        printf ( "Start at :  %s", asctime (timeinfo) );
    //// <<<<<<<<<  affichage de l'heure de debut

    int runId;
    for (runId=0; (runId<nbRun) || ((nbRun==-1) && (nbFound==0)) ; runId++) {
      clock_t startTime=clock();
      struct timeval tv;
      gettimeofday(&tv, NULL);
      double humanTime = (double) tv.tv_sec + (double) tv.tv_usec / 1000000.0;

      //// LANCEMENT DU CALCUL /////////////////
      solver->compute();

      double elapsedTime=(clock()-startTime) / (double)CLOCKS_PER_SEC/60.0;
      gettimeofday(&tv, NULL);
      humanTime = ((double) tv.tv_sec + (double) tv.tv_usec / 1000000.0 - humanTime)/60.0;
      
      printf("nb conflicted edges: %d\n", solver->bestSol.nbEdgesConflict);
      if(!quite){
        printf("best coloring: ");
        solver->bestSol.breakSymmetry();
        solver->bestSol.print();
      }
      
      // Save coloring(s) in outputfile
      if (strlen(outputFile) != 0)
        solver->saveBestColoring(outputFile);

      /////// Affichage de l'heure
      time ( &rawtime );
      timeinfo = localtime ( &rawtime );
      if(!quite)
          printf ( "\tFinished :  %s  (cpu time: %fmin, humain time : %fmin)\n", asctime (timeinfo), elapsedTime, humanTime );
      fflush(stdout);
      ///////

      totalIterationsCrossWithWrongRun+=solver->nbIterationsCross;
      if (solver->bestSol.nbEdgesConflict==0) {
        //solver->save(elapsedTime);
        nbFound++;
        totalCpuTime+=elapsedTime;
        totalHumanTime+=humanTime;
        totalIterations+=solver->nbIterations;
        totalIterationsCross+=solver->nbIterationsCross;
      }

      fflush(stdout);
      
    }
    printf("\n#Success / #Runs : %d / %d\n", nbFound,runId);

    if (nbFound>0) {
      printf("Mean CPU time   : %0.2f min\n", totalCpuTime/nbFound);
      printf("Mean humain time: %0.2f min\n", totalHumanTime/nbFound);
      printf("Mean iterations : %0.3f (x10.6)\n", totalIterations/nbFound/1000000.0);
      printf("Mean crossover  : %llu \n", totalIterationsCross/nbFound);
    }

    printf("End\n");
    fflush(stdout);
    if (nbFound) {
      for (int i = 0; i < g.nbSommets; i++) {
        s.cols[i] = solver->bestSol.tColor[i];
      }
    }
    s.recompute_num_cols();
    //assert(s.verify());
    s.save_if_better_than_default();
  }
}
