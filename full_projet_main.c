/*
 *  main.c
 *  Projet 2022-2023
 *  Elasticite lineaire plane
 *
 *  Code de calcul
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#include "fem.h"
#include <time.h>

int main(void) {
  
  
  femGeo *theGeometry = geoGetGeometry();
  geoMeshRead("../data/mesh.txt");
  femProblem *theProblem = femElasticityRead(theGeometry, "../data/problem.txt");
  femElasticityPrint(theProblem);
  double start = clock();
  double *theSoluce = femElasticitySolve(theProblem);
  double stop = clock();
  printf("Temps de calcul : %f secondes\n", (stop - start) / CLOCKS_PER_SEC);
  fflush(stdout);
  int nNodes = theGeometry->theNodes->nNodes;
  femSolutionWrite(nNodes, 2, theSoluce, "../data/UV.txt");
  femElasticityFree(theProblem);
  geoFree();
  return 0;
}
