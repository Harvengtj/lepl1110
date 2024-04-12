#include "fem.h"

// Il faut un fifrelin generaliser ce code.....
//  (1) Ajouter l'axisymétrique !    (mandatory)
//  (2) Ajouter les conditions de Neumann !   (mandatory)
//  (3) Ajouter les conditions en normal et tangentiel !   (strongly advised)
//  (4) Et remplacer le solveur plein par un truc un fifrelin plus subtil  (mandatory)


/////////////////////////////////////////////////////////////////////////////////////
/* Zone de renumérotation */
double *theGlobalCoord;

int compare(const void *nodeOne, const void *nodeTwo) 
{
    int *iOne = (int *)nodeOne;
    int *iTwo = (int *)nodeTwo;
    double diff = theGlobalCoord[*iOne] - theGlobalCoord[*iTwo];
    if (diff < 0)    return  1;
    if (diff > 0)    return -1;
    return  0;  
}

void femMeshRenumber(femMesh *theMesh, femRenumType renumType)
{
    int i, *inverse;
    femNodes *theNodes = theMesh->nodes;
    int nNodes = theNodes->nNodes;
    int *number = theNodes->number;
    
    switch (renumType) {
        case FEM_NO :
            // Les noeuds sont numérotés dans l'ordre d'origine
            for (i = 0; i < nNodes; i++) 
                number[i] = i;
            break;
        case FEM_XNUM : 
            inverse = malloc(sizeof(int)*nNodes);
            for (i = 0; i < nNodes; i++) 
                inverse[i] = i; 
            // Les coordonnées globales des noeuds seront basées sur les coordonnées X
            theGlobalCoord = theNodes->X;
            // Cela réorganise les indices des noeuds en fonction de leurs coordonnées globales
            qsort(inverse, nNodes, sizeof(int), compare);
            // Cela renumérote les nœuds en fonction de leurs coordonnées globales triées
            for (i = 0; i < nNodes; i++)
                number[inverse[i]] = i;
            free(inverse);  
            break;
        case FEM_YNUM : 
            inverse = malloc(sizeof(int)*nNodes);
            for (i = 0; i < nNodes; i++) 
                inverse[i] = i; 
            // Les coordonnées globales des noeuds seront basées sur les coordonnées Y
            theGlobalCoord = theNodes->Y;
            // Cela réorganise les indices des noeuds en fonction de leurs coordonnées globales
            qsort(inverse, nNodes, sizeof(int), compare);
            // Cela renumérote les nœuds en fonction de leurs coordonnées globales triées
            for (i = 0; i < nNodes; i++)
                number[inverse[i]] = i;
            free(inverse);  
            break;
        default : Error("Unexpected renumbering option"); }
}

int femMeshComputeBand(femMesh *theMesh)
{
    int iElem,j,myMax,myMin,myBand,map[4];
    int nLocal = theMesh->nLocalNode;
    femNodes *theNodes = theMesh->nodes;
    int *number = theNodes->number;

    myBand = 0;
    // On parcourt chaque élément (triangle ou quad)
    for(iElem = 0; iElem < theMesh->nElem; iElem++) {
        // On parcourt chaque noeud de l'élément (triangle ou quad)
        for (j=0; j < nLocal; ++j) 
            // On récupère les (3 ou 4) indices des noeuds correspondant à l'élément
            map[j] = number[theMesh->elem[iElem*nLocal+j]];
        myMin = map[0];
        myMax = map[0];
        // On détermine l'indice le plus petit et l'indice le plus grand
        for (j=1; j < nLocal; j++) {
            myMax = fmax(map[j],myMax);
            myMin = fmin(map[j],myMin); }
        if (myBand < (myMax - myMin)) myBand = myMax - myMin; }         
    return(++myBand);
}

/////////////////////////////////////////////////////////////////////////////////////


// On assemble les éléments de la matrice A et du vecteur B (A MON AVIS A MODIFIER POUR L'ADAPTER AU SOLVEUR BANDE)
void femElasticityAssembleElements(femProblem *theProblem) {
  femBandSystem *theSystem = theProblem->bandSystem;
  femIntegration *theRule = theProblem->rule;
  femDiscrete *theSpace = theProblem->space;
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;
  femMesh *theMesh = theGeometry->theElements;
  double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
  int iElem, iInteg, iEdge, i, j, d, map[4], mapX[4], mapY[4];
  int nLocal = theMesh->nLocalNode;
  double a = theProblem->A;
  double b = theProblem->B;
  double c = theProblem->C;
  double rho = theProblem->rho;
  double gx = theProblem->gx;
  double gy = theProblem->gy;
  double **A = theSystem->A;
  double *B = theSystem->B;

  // On parcourt TOUS les éléments
  for (iElem = 0; iElem < theMesh->nElem; iElem++) {
    for (j = 0; j < nLocal; j++) {
      map[j] = theMesh->elem[iElem * nLocal + j];
      mapX[j] = 2 * map[j];
      mapY[j] = 2 * map[j] + 1;
      x[j] = theNodes->X[map[j]];
      y[j] = theNodes->Y[map[j]];
    }

    // On effectue la règle de Hammer à 3 poids pour l'élément
    for (iInteg = 0; iInteg < theRule->n; iInteg++) {
      double xsi = theRule->xsi[iInteg];
      double eta = theRule->eta[iInteg];
      double weight = theRule->weight[iInteg];
      femDiscretePhi2(theSpace, xsi, eta, phi);
      femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

      double dxdxsi = 0.0;
      double dxdeta = 0.0;
      double dydxsi = 0.0;
      double dydeta = 0.0;
      for (i = 0; i < theSpace->n; i++) {
        dxdxsi += x[i] * dphidxsi[i];
        dxdeta += x[i] * dphideta[i];
        dydxsi += y[i] * dphidxsi[i];
        dydeta += y[i] * dphideta[i];
      }
      double jac = dxdxsi * dydeta - dxdeta * dydxsi;
      if (jac < 0.0)
        printf("Negative jacobian! Your mesh is oriented in reverse. The normals will be wrong\n");
      jac = fabs(jac);

      for (i = 0; i < theSpace->n; i++) {
        dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
        dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
      }

      // On remplit la matrice A globale par les valeurs contenues dans la matrice A locale
      for (i = 0; i < theSpace->n; i++) {
        for (j = 0; j < theSpace->n; j++) {
          A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * jac * weight;
          A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * jac * weight;
          A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * jac * weight;
          A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * jac * weight;
        }
      }

      // On remplit le vecteur B
      for (i = 0; i < theSpace->n; i++) {
        B[mapX[i]] += phi[i] * gx * rho * jac * weight;
        B[mapY[i]] += phi[i] * gy * rho * jac * weight;
      }
    }
  }
}

// Calcule les intégrales de ligne associées aux conditions de Neumann
void femElasticityAssembleNeumann(femProblem *theProblem) {
  femBandSystem *theSystem = theProblem->bandSystem;
  femIntegration *theRule = theProblem->ruleEdge;
  femDiscrete *theSpace = theProblem->spaceEdge;
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;
  femMesh *theEdges = theGeometry->theEdges;

  // Ici, on a affaire à des arêtes (2D) -> tableaux à deux dimensions
  double x[2], y[2], phi[2];
  int iBnd, iElem, iInteg, iEdge, i, j, d, map[2];
  int nLocal = 2;
  double *B = theSystem->B;

  // On parcourt les conditions (une condition est associée à un domaine)
  for (iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++) {
    femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
    femBoundaryType type = theCondition->type;
    double value = theCondition->value1;

    // S'il n'y a pas de condition frontière de type Neumann, on passe à l'itération suivante
    if(type != NEUMANN_X && type != NEUMANN_Y && type != NEUMANN_N && type != NEUMANN_T){
      continue;
    }

    // On parcourt les arêtes du domaine contraint
    for (iEdge = 0; iEdge < theCondition->domain->nElem; iEdge++) {
      iElem = theCondition->domain->elem[iEdge];
      for (j = 0; j < nLocal; j++) {
        map[j] = theEdges->elem[iElem * nLocal + j];
        x[j] = theNodes->X[map[j]];
        y[j] = theNodes->Y[map[j]];
      }

      // Calcul du Jacobien pour l'intégrale de ligne
      double tx = x[1] - x[0];
      double ty = y[1] - y[0];
      double length = hypot(tx, ty);
      double jac = length / 2.0;
      
      double f_x = 0.0;
      double f_y = 0.0;
      if (type == NEUMANN_X) {
        f_x = value;
      }
      if (type == NEUMANN_Y) {
        f_y = value;
      }

      //
      // A completer :-)
      // Attention, pour le normal tangent on calcule la normale (sortante) au SEGMENT, surtout PAS celle de constrainedNodes
      // Une petite aide pour le calcul de la normale :-)
      double nx =  ty / length;
      double ny = -tx / length;
      if (type == NEUMANN_N) {
        f_x = value * nx;
        f_y = value * ny;
      }
      if (type == NEUMANN_T) {
        f_x = value * (-ny);
        f_y = value * nx;
      }

      for (iInteg = 0; iInteg < theRule->n; iInteg++) {
        double xsi = theRule->xsi[iInteg];
        double weight = theRule->weight[iInteg];
        femDiscretePhi(theSpace, xsi, phi);
        for (i = 0; i < theSpace->n; i++) {
          B[2*map[i] + 0] += jac * weight * phi[i] * f_x;
          B[2*map[i] + 1] += jac * weight * phi[i] * f_y;
        }
      }
    }
  }
}

// On applique les conditions de Dirichlet
void femElasticityApplyDirichlet(femProblem *theProblem) {
  femBandSystem *theSystem = theProblem->bandSystem;
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;

  for (int node = 0; node < theNodes->nNodes; node++) {
    femConstrainedNode *theConstrainedNode = &theProblem->constrainedNodes[node];
    // Si le noeud n'est pas contraint, on passe à l'itération suivante
    if (theConstrainedNode->type == UNDEFINED)
      continue;
    femBoundaryType type = theConstrainedNode->type;

    if (type == DIRICHLET_X) {
      double value = theConstrainedNode->value1;
      femBandSystemConstrain(theSystem, 2 * node + 0, value);
    }
    if (type == DIRICHLET_Y) {
      double value = theConstrainedNode->value1;
      femBandSystemConstrain(theSystem, 2 * node + 1, value);
    }
    if (type == DIRICHLET_XY) {
      double value_x = theConstrainedNode->value1;
      double value_y = theConstrainedNode->value2;
      femBandSystemConstrain(theSystem, 2 * node + 0, value_x);
      femBandSystemConstrain(theSystem, 2 * node + 1, value_y);
    }

    if (type == DIRICHLET_N) {
      double value = theConstrainedNode->value1;
      double nx = theConstrainedNode->nx;
      double ny = theConstrainedNode->ny;
      // A completer :-)
      femBandSystemConstrain(theSystem, 2 * node + 0, value * nx);
      femBandSystemConstrain(theSystem, 2 * node + 1, value * ny);
    }
    if (type == DIRICHLET_T) {
      double value = theConstrainedNode->value1;
      double nx = theConstrainedNode->nx;
      double ny = theConstrainedNode->ny;
      // A completer :-)
      femBandSystemConstrain(theSystem, 2 * node + 0, value * ny);
      femBandSystemConstrain(theSystem, 2 * node + 1, value * nx);
    }
    if (type == DIRICHLET_NT) {
      double value_n = theConstrainedNode->value1;
      double value_t = theConstrainedNode->value2;
      double nx = theConstrainedNode->nx;
      double ny = theConstrainedNode->ny;
      // A completer :-)
      femBandSystemConstrain(theSystem, 2 * node + 0, value_n * nx + value_t * ny);
      femBandSystemConstrain(theSystem, 2 * node + 1, value_n * ny + value_t * nx);
    }
  }
}

double *femElasticitySolve(femProblem *theProblem) {
  femElasticityAssembleElements(theProblem);
  femElasticityAssembleNeumann(theProblem);
  femElasticityApplyDirichlet(theProblem);

  double *soluce = femBandSystemEliminate(theProblem->bandSystem);
  memcpy(theProblem->soluce, soluce, theProblem->system->size * sizeof(double));
  return theProblem->soluce;
}
