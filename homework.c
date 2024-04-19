#include "fem.h"

// Il faut un fifrelin generaliser ce code.....
//  (1) Ajouter l'axisymétrique !    (mandatory)
//  (2) Ajouter les conditions de Neumann !   (mandatory)
//  (3) Ajouter les conditions en normal et tangentiel !   (strongly advised)
//  (4) Et remplacer le solveur plein par un truc un fifrelin plus subtil  (mandatory)

// Trouve la coordonnée X du barycentre de l'élément
double findBaryCenterX(int iElem, femMesh *theElements) {
  double X = 0.0;
  int nLocalNode = theElements->nLocalNode;
  for (int i = 0; i < nLocalNode; i++) {
    double coordX = theElements->nodes->X[iElem*nLocalNode + i];
    X += coordX;
  }
  return (X/nLocalNode);
}

// Trouve la coordonnée Y du barycentre de l'élement
double findBaryCenterY(int iElem, femMesh *theElements) {
  double Y = 0.0;
  int nLocalNode = theElements->nLocalNode;
  for (int i = 0; i < nLocalNode; i++) {
    double coordY = theElements->nodes->Y[iElem*nLocalNode + i];
    Y += coordY;
  }
  return (Y/nLocalNode);
}


void femElasticityAssembleElements(femProblem *theProblem) {
  femFullSystem *theSystem = theProblem->system;
  femIntegration *theRule = theProblem->rule;
  femDiscrete *theSpace = theProblem->space;
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;
  femMesh *theMesh = theGeometry->theElements;
  double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
  int iElem, iInteg, iEdge, i, j, d, map[4], mapX[4], mapY[4];
  int nLocal = theMesh->nLocalNode;
  // MODIFIE
  double a;
  double b;
  double c;
  double rho = theProblem->rho;
  double gx = theProblem->gx;
  double gy = theProblem->gy;
  double **A = theSystem->A;
  double *B = theSystem->B;

  for (iElem = 0; iElem < theMesh->nElem; iElem++) {
    // Trouve les coordonnées des barycentres
    double X = findBaryCenterX(iElem, theMesh);
    double Y = findBaryCenterY(iElem, theMesh);

    
    // Distinction des matériaux
    if (Y < 0.85 && Y > (-0.85)) {
      a = theProblem->A1;
      b = theProblem->B1;
      c = theProblem->C1;
    } else {
      a = theProblem->A2;
      b = theProblem->B2;
      c = theProblem->C2;
    }

    for (j = 0; j < nLocal; j++) {
      map[j] = theMesh->elem[iElem * nLocal + j];
      mapX[j] = 2 * map[j];
      mapY[j] = 2 * map[j] + 1;
      x[j] = theNodes->X[map[j]];
      y[j] = theNodes->Y[map[j]];
    }

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
      for (i = 0; i < theSpace->n; i++) {
        for (j = 0; j < theSpace->n; j++) {
          A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * jac * weight;
          A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * jac * weight;
          A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * jac * weight;
          A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * jac * weight;
        }
      }
      for (i = 0; i < theSpace->n; i++) {
        B[mapX[i]] += phi[i] * gx * rho * jac * weight;
        B[mapY[i]] += phi[i] * gy * rho * jac * weight;
      }
    }
  }
}

void femElasticityAssembleNeumann(femProblem *theProblem) {
  femFullSystem *theSystem = theProblem->system;
  femIntegration *theRule = theProblem->ruleEdge;
  femDiscrete *theSpace = theProblem->spaceEdge;
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;
  femMesh *theEdges = theGeometry->theEdges;
  double x[2], y[2], phi[2];
  int iBnd, iElem, iInteg, iEdge, i, j, d, map[2];
  int nLocal = 2;
  double *B = theSystem->B;

  for (iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++) {
    femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
    femBoundaryType type = theCondition->type;
    double value = theCondition->value1;

    if(type != NEUMANN_X && type != NEUMANN_Y && type != NEUMANN_N && type != NEUMANN_T){
      continue;
    }

    for (iEdge = 0; iEdge < theCondition->domain->nElem; iEdge++) {
      iElem = theCondition->domain->elem[iEdge];
      for (j = 0; j < nLocal; j++) {
        map[j] = theEdges->elem[iElem * nLocal + j];
        x[j] = theNodes->X[map[j]];
        y[j] = theNodes->Y[map[j]];
      }

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
      // double nx =  ty / length;
      // double ny = -tx / length;

      if (type == NEUMANN_N) {
        double nx =  ty / length;
        double ny = -tx / length;
        f_x += value *nx;
        f_y += value * ny;
      }
      if (type == NEUMANN_T) {
        f_x += value * tx/ length;
        f_y += value * ty/ length;
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
  femFullSystem *theSystem = theProblem->system;
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;

  for (int node = 0; node < theNodes->nNodes; node++) {
    femConstrainedNode *theConstrainedNode = &theProblem->constrainedNodes[node];
    if (theConstrainedNode->type == UNDEFINED)
      continue;
    femBoundaryType type = theConstrainedNode->type;

    if (type == DIRICHLET_X) {
      double value = theConstrainedNode->value1;
      femFullSystemConstrain(theSystem, 2 * node + 0, value);
    }
    if (type == DIRICHLET_Y) {
      double value = theConstrainedNode->value1;
      femFullSystemConstrain(theSystem, 2 * node + 1, value);
    }
    if (type == DIRICHLET_XY) {
      double value_x = theConstrainedNode->value1;
      double value_y = theConstrainedNode->value2;
      femFullSystemConstrain(theSystem, 2 * node + 0, value_x);
      femFullSystemConstrain(theSystem, 2 * node + 1, value_y);
    }

    if (type == DIRICHLET_N) {
      double value = theConstrainedNode->value1;
      double nx = theConstrainedNode->nx;
      double ny = theConstrainedNode->ny;
      // A completer 🙂
      double tx = -ny;
      double ty = nx;

      for (int i = 0; i < theSystem->size; i++) {
        
        theSystem->A[i][2 * node + 0] = tx*(tx*theSystem->A[i][2 * node + 0] + ty*theSystem->A[i][2 * node + 1]);
        theSystem->A[i][2 * node + 1] = ty*(tx*theSystem->A[i][2 * node + 0] + ty*theSystem->A[i][2 * node + 1]);

        theSystem->B[i] -= theSystem->A[i][2 * node + 0] * value * nx + theSystem->A[i][2 * node + 1] * value * ny;
      }
    }
    if (type == DIRICHLET_T) {
      double value = theConstrainedNode->value1;
      double nx = theConstrainedNode->nx;
      double ny = theConstrainedNode->ny;
      // A completer 🙂
      double tx = -ny;
      double ty = nx;

      for (int i = 0; i < theSystem->size; i++) {
        theSystem->A[i][2 * node + 0] = tx;
        theSystem->A[i][2 * node + 1] = 0.0;
        theSystem->A[2 * node + 0][i] = 0.0;
        theSystem->A[2 * node + 1][i] = 0.0;
        theSystem->B[i] -= theSystem->A[i][2 * node + 0] * value * nx + theSystem->A[i][2 * node + 1] * value * ny;
        theSystem->A[i][2 * node + 0] = 0.0;
        theSystem->A[i][2 * node + 1] = 0.0;
      }

      femFullSystemConstrain(theSystem, 2 * node + 0, -value * ny);

      for (int i = 0; i < theSystem->size; i++) {
        theSystem->B[i] -= theSystem->A[i][2 * node + 1] * value * nx;
        theSystem->A[i][2 * node + 1] = 0.0;
      }
    }
    if (type == DIRICHLET_NT) {
      double value_n = theConstrainedNode->value1;
      double value_t = theConstrainedNode->value2;
      double nx = theConstrainedNode->nx;
      double ny = theConstrainedNode->ny;
      // A completer 🙂
      femFullSystemConstrain(theSystem, 2 * node + 0, value_n * nx - value_t * ny);
      femFullSystemConstrain(theSystem, 2 * node + 1, value_n * ny + value_t * nx);
    }
  }
}




//////////////////////////////////////////////////////////////////////////////////////////////////////
/* SOLVEUR BANDE */

typedef enum {FEM_NO,FEM_XNUM,FEM_YNUM} femRenumType;

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


// Calcule la largeur de bande de la matrice A de dimensions size x size
int calculateBandSize(double **A, int size) {
  int upper_band = 0;
  int lower_band = 0;

  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      if (A[i][j] != 0) {
        int diff = i - j;
        if (diff < 0) {
          diff = -diff; // Prendre la valeur absolue si diff est négatif
          if (diff > upper_band) {
            upper_band = diff;
          }
        } else {
          if (diff > lower_band) {
            lower_band = diff;
          }
        }
      }
    }
  }

  return (upper_band > lower_band) ? upper_band : lower_band;
}

/*
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
*/

// Fonction qui renumérote les noeuds
int *femMeshRenumber(femMesh *theMesh, femRenumType renumType)
{
    int i, *inverse;
    femNodes *theNodes = theMesh->nodes;
    int nNodes = theNodes->nNodes;
    int *number = malloc(sizeof(int)*nNodes);
    
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
        default : Error("Unexpected renumbering option"); 
        }
    return number;
}

double *femBandSystemEliminate(femProblem *theProblem) {
  double **A, *B, factor;
  int i, j, k, jend, size, band;

  A = theProblem->system->A;
  B = theProblem->system->B;
  size = theProblem->system->size;

  femRenumType renumType = FEM_NO;
  // femRenumType renumType = FEM_XNUM;
  // femRenumType renumType = FEM_YNUM;

  printf("La taille du systeme est %d\n", size);

  // Calcul de la largeur de bande
  band = calculateBandSize(A, size);

  printf("La largeur de bande sans renumerotation est %d\n", band);

  // Rénumérotation des noeuds
  int *renumber = femMeshRenumber(theProblem->geometry->theElements, renumType);
  printf("Renumerotation des noeuds effectuee\n");

  // Création des matrices renumérotées
  double **Arenumbered = malloc(sizeof(double*)*size);
  for (int i = 0; i < size; i++) {
    Arenumbered[i] = malloc(sizeof(double)*size);
  }

  double *Brenumbered = malloc(sizeof(double)*size);
  for (int i = 0; i < size/2; i++) {
    Brenumbered[2*renumber[i]] = B[2*i];
    Brenumbered[2*renumber[i]+1] = B[(2*i)+1];

    for (int j = 0; j < size/2; j++) {
      Arenumbered[2*renumber[i]][2*renumber[j]] = A[2*i][2*j];
      Arenumbered[2*renumber[i]][2*renumber[j]+1] = A[2*i][(2*j)+1];
      Arenumbered[2*renumber[i]+1][2*renumber[j]] = A[(2*i)+1][2*j];
      Arenumbered[2*renumber[i]+1][2*renumber[j]+1] = A[(2*i)+1][(2*j)+1];
    }
  }
  printf("Matrices renumerotees creees\n");
  if (renumType == FEM_XNUM || renumType == FEM_YNUM) {
    band = calculateBandSize(Arenumbered, size) * 2;
  } else {
    band = calculateBandSize(Arenumbered, size);
  }
  
  printf("La largeur de bande apres renumerotation est %d\n", band);


  /* Gauss elimination */
  for (k=0; k < size; k++) {
        // Erreur si le pivot est à 0
        if ( fabs(Arenumbered[k][k]) <= 1e-4 ) {
            Error("Cannot eleminate with such a pivot"); }
        // On doit faire l'élimination gaussienne dans la largeur de bande et pas en dehors
        jend = fmin(k + band,size);
        for (i = k+1 ; i <  jend; i++) {
            factor = Arenumbered[k][i] / Arenumbered[k][k];
            for (j = i ; j < jend; j++) 
                Arenumbered[i][j] = Arenumbered[i][j] - Arenumbered[k][j] * factor;
                Brenumbered[i] = Brenumbered[i] - Brenumbered[k] * factor; }}
        
    /* Back-substitution */

    for (i = (size-1); i >= 0 ; i--) {
        factor = 0;
        jend = fmin(i + band,size);
        for (j = i+1 ; j < jend; j++) 
            factor += Arenumbered[i][j] * Brenumbered[j];
            Brenumbered[i] = ( Brenumbered[i] - factor)/Arenumbered[i][i];}
  
  printf("Echelonnement termine\n");

  double *soluce = malloc(sizeof(double)*size);

  for (i = 0; i < size/2; i++) {
    soluce[2*i] = Brenumbered[2*renumber[i]];
    soluce[(2*i)+1] = Brenumbered[2*renumber[i]+1];
  }

  printf("Derenumerotation effectuee\n");

  // Libération de la mémoire
  for (i = 0; i < size; i++) {
    free(Arenumbered[i]);
  }
  free(Arenumbered);
  free(Brenumbered);
  free(renumber);
  printf("Memoire liberee\n");

  return soluce;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////


double *femElasticitySolve(femProblem *theProblem) {
  femElasticityAssembleElements(theProblem);
  femElasticityAssembleNeumann(theProblem);
  femElasticityApplyDirichlet(theProblem);

  double *soluce = femBandSystemEliminate(theProblem);
  memcpy(theProblem->soluce, soluce, theProblem->system->size * sizeof(double));
  return theProblem->soluce;
}


