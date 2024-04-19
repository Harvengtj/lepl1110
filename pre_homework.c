#include "fem.h"

//
// Ici, vous pouvez définir votre géométrie :-)
//  (1) Raffiner intelligemment.... (yes )
//  (2) Construire la geometrie avec OpenCascade
//  (3) Construire la geometrie avec les outils de GMSH
//  (4) Obtenir la geometrie en lisant un fichier .geo de GMSH

double geoSize(double x, double y) {

  femGeo *theGeometry = geoGetGeometry();
  return theGeometry->h;
}



// Le but de cette fonction est de générer un octogone
void geoMeshGenerate(void) {

  femGeo *theGeometry = geoGetGeometry();
  double Lx = 1.0;
  double Ly = 2.0;
  theGeometry->LxPlate = Lx;
  theGeometry->LyPlate = Ly;
  theGeometry->h = Lx * 0.05;
  theGeometry->elementType = FEM_QUAD;

  geoSetSizeCallback(geoSize);

  double w = theGeometry->LxPlate;
  double h = theGeometry->LyPlate;

  double w1 = w - w / 8.0;
  double h1 = h - w / 8.0;

  int ierr; 

  // Creation de la plaque principale 
  int idPlate = gmshModelOccAddRectangle(-w/2.0,-h/2.0,0.0,w,h,-1,0.0,&ierr); 
  ErrorGmsh(ierr);
  // Creation de la plaque à decouper
  int idPlatecut = gmshModelOccAddRectangle(-w1/2.0,-h1/2.0,0.0,w1,h1,-1,0.0,&ierr); 
  ErrorGmsh(ierr);

  // Ajout des points pour les triangles à cut sur la plaque principale
  // En bas à gauche
  int idFirstPoint = gmshModelOccAddPoint(-w/2.0, -h/2.0, 0.0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  int idSecondPoint = gmshModelOccAddPoint(-w/2.0 + w/10.0, -h/2.0, 0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  int idThirdPoint = gmshModelOccAddPoint(-w/2.0, -h/2.0 + w/10.0, 0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  // En haut à gauche
  int idFourthPoint = gmshModelOccAddPoint(-w/2.0, h/2.0, 0.0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  int idFifthPoint = gmshModelOccAddPoint(-w/2.0 + w/10.0, h/2.0, 0.0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  int idSixthPoint = gmshModelOccAddPoint(-w/2.0, h/2.0 - w/10.0, 0.0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  // En haut à droite
  int idSeventhPoint = gmshModelOccAddPoint(w/2.0, h/2.0, 0.0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  int idEighthPoint = gmshModelOccAddPoint(w/2.0 - w/10.0, h/2.0, 0.0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  int idNinthPoint = gmshModelOccAddPoint(w/2.0, h/2.0 - w/10.0, 0.0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  // En bas à droite
  int idTenthPoint = gmshModelOccAddPoint(w/2.0, -h/2.0, 0.0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  int idEleventhPoint = gmshModelOccAddPoint(w/2.0 - w/10.0, -h/2.0, 0.0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  int idTwelfthPoint = gmshModelOccAddPoint(w/2.0, -h/2.0 + w/10.0, 0.0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  
  
  

  
  // Ajout des lignes pour les triangles à cut sur la plaque principale
  int idFirstLine = gmshModelOccAddLine(idFirstPoint, idSecondPoint, -1, &ierr);
  ErrorGmsh(ierr);
  int idSecondLine = gmshModelOccAddLine(idSecondPoint, idThirdPoint, -1, &ierr);
  ErrorGmsh(ierr);
  int idThirdLine = gmshModelOccAddLine(idThirdPoint, idFirstPoint, -1, &ierr);
  ErrorGmsh(ierr);
  int idFourthLine = gmshModelOccAddLine(idFourthPoint, idFifthPoint, -1, &ierr);
  ErrorGmsh(ierr);
  int idFifthLine = gmshModelOccAddLine(idFifthPoint, idSixthPoint, -1, &ierr);
  ErrorGmsh(ierr);
  int idSixthLine = gmshModelOccAddLine(idSixthPoint, idFourthPoint, -1, &ierr);
  ErrorGmsh(ierr);
  int idSeventhLine = gmshModelOccAddLine(idSeventhPoint, idEighthPoint, -1, &ierr);
  ErrorGmsh(ierr);
  int idEighthLine = gmshModelOccAddLine(idEighthPoint, idNinthPoint, -1, &ierr);
  ErrorGmsh(ierr);
  int idNinthLine = gmshModelOccAddLine(idNinthPoint, idSeventhPoint, -1, &ierr);
  ErrorGmsh(ierr);
  int idTenthLine = gmshModelOccAddLine(idTenthPoint, idEleventhPoint, -1, &ierr);
  ErrorGmsh(ierr);
  int idEleventhLine = gmshModelOccAddLine(idEleventhPoint, idTwelfthPoint, -1, &ierr);
  ErrorGmsh(ierr);
  int idTwelfthLine = gmshModelOccAddLine(idTwelfthPoint, idTenthPoint, -1, &ierr);
  ErrorGmsh(ierr);
  

  // Ajout de la boucle de courbe pour les triangles à cut sur la plaque principale
  int lines[] = {idFirstLine, idSecondLine, idThirdLine};
  int linesloop = gmshModelOccAddCurveLoop(lines, 3, -1, &ierr);
  int lines_bis[] = {idFourthLine, idFifthLine, idSixthLine};
  int linesloop_bis = gmshModelOccAddCurveLoop(lines_bis, 3, -1, &ierr);
  int lines_ter[] = {idSeventhLine, idEighthLine, idNinthLine};
  int linesloop_ter = gmshModelOccAddCurveLoop(lines_ter, 3, -1, &ierr);
  int lines_quater[] = {idTenthLine, idEleventhLine, idTwelfthLine};
  int linesloop_quater = gmshModelOccAddCurveLoop(lines_quater, 3, -1, &ierr);
  // Ajout de la surface plane pour les triangles à cut sur la plaque principale
  int surfaces[] = {linesloop};
  int idFirstSurface = gmshModelOccAddPlaneSurface(surfaces, 1, -1, &ierr);
  ErrorGmsh(ierr);
  int surfaces_bis[] = {linesloop_bis};
  int idSecondSurface = gmshModelOccAddPlaneSurface(surfaces_bis, 1, -1, &ierr);
  ErrorGmsh(ierr);
  int surfaces_ter[] = {linesloop_ter};
  int idThirdSurface = gmshModelOccAddPlaneSurface(surfaces_ter, 1, -1, &ierr);
  ErrorGmsh(ierr);
  int surfaces_quater[] = {linesloop_quater};
  int idFourthSurface = gmshModelOccAddPlaneSurface(surfaces_quater, 1, -1, &ierr);
  ErrorGmsh(ierr);

  // Ajout des points pour les triangles à cut sur la plaque à decouper
  int idFirstPointcut = gmshModelOccAddPoint(-w1 / 2.0, -h1 / 2.0, 0.0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  int idSecondPointcut = gmshModelOccAddPoint(-w1 / 2.0 + w1 / 10.0, -h1 / 2.0, 0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  int idThirdPointcut = gmshModelOccAddPoint(-w1 / 2.0, -h1 / 2.0 + w1 / 10.0, 0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  int idFourthPointcut = gmshModelOccAddPoint(-w1/2.0, h1/2.0, 0.0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  int idFifthPointcut = gmshModelOccAddPoint(-w1/2.0 + w1/10.0, h1/2.0, 0.0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  int idSixthPointcut = gmshModelOccAddPoint(-w1/2.0, h1/2.0 - w1/10.0, 0.0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  int idSeventhPointcut = gmshModelOccAddPoint(w1/2.0, h1/2.0, 0.0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  int idEighthPointcut = gmshModelOccAddPoint(w1/2.0 - w1/10.0, h1/2.0, 0.0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  int idNinthPointcut = gmshModelOccAddPoint(w1/2.0, h1/2.0 - w1/10.0, 0.0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  int idTenthPointcut = gmshModelOccAddPoint(w1/2.0, -h1/2.0, 0.0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  int idEleventhPointcut = gmshModelOccAddPoint(w1/2.0 - w1/10.0, -h1/2.0, 0.0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  int idTwelfthPointcut = gmshModelOccAddPoint(w1/2.0, -h1/2.0 + w1/10.0, 0.0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  // Ajout des lignes pour les triangles à cut sur la plaque à decouper
  int idFirstLinecut = gmshModelOccAddLine(idFirstPointcut, idSecondPointcut, -1, &ierr);
  ErrorGmsh(ierr);
  int idSecondLinecut = gmshModelOccAddLine(idSecondPointcut, idThirdPointcut, -1, &ierr);
  ErrorGmsh(ierr);
  int idThirdLinecut = gmshModelOccAddLine(idThirdPointcut, idFirstPointcut, -1, &ierr);
  ErrorGmsh(ierr);
  int idFourthLinecut = gmshModelOccAddLine(idFourthPointcut, idFifthPointcut, -1, &ierr);
  ErrorGmsh(ierr);
  int idFifthLinecut = gmshModelOccAddLine(idFifthPointcut, idSixthPointcut, -1, &ierr);
  ErrorGmsh(ierr);
  int idSixthLinecut = gmshModelOccAddLine(idSixthPointcut, idFourthPointcut, -1, &ierr);
  ErrorGmsh(ierr);
  int idSeventhLinecut = gmshModelOccAddLine(idSeventhPointcut, idEighthPointcut, -1, &ierr);
  ErrorGmsh(ierr);
  int idEighthLinecut = gmshModelOccAddLine(idEighthPointcut, idNinthPointcut, -1, &ierr);
  ErrorGmsh(ierr);
  int idNinthLinecut = gmshModelOccAddLine(idNinthPointcut, idSeventhPointcut, -1, &ierr);
  ErrorGmsh(ierr);
  int idTenthLinecut = gmshModelOccAddLine(idTenthPointcut, idEleventhPointcut, -1, &ierr);
  ErrorGmsh(ierr);
  int idEleventhLinecut = gmshModelOccAddLine(idEleventhPointcut, idTwelfthPointcut, -1, &ierr);
  ErrorGmsh(ierr);
  int idTwelfthLinecut = gmshModelOccAddLine(idTwelfthPointcut, idTenthPointcut, -1, &ierr);
  ErrorGmsh(ierr);
  // Ajout de la boucle de courbe pour les triangles à cut sur la plaque à decouper
  int linescut[] = {idFirstLinecut, idSecondLinecut, idThirdLinecut};
  int linesloopcut = gmshModelOccAddCurveLoop(linescut, 3, -1, &ierr);
  ErrorGmsh(ierr);
  int lines_biscut[] = {idFourthLinecut, idFifthLinecut, idSixthLinecut};
  int linesloop_biscut = gmshModelOccAddCurveLoop(lines_biscut, 3, -1, &ierr);
  ErrorGmsh(ierr);
  int lines_tercut[] = {idSeventhLinecut, idEighthLinecut, idNinthLinecut};
  int linesloop_tercut = gmshModelOccAddCurveLoop(lines_tercut, 3, -1, &ierr);
  ErrorGmsh(ierr);
  int lines_quatercut[] = {idTenthLinecut, idEleventhLinecut, idTwelfthLinecut};
  int linesloop_quatercut = gmshModelOccAddCurveLoop(lines_quatercut, 3, -1, &ierr);
  ErrorGmsh(ierr);
  // Ajout de la surface plane pour les triangles à cut sur la plaque à decouper
  int surfacescut[] = {linesloopcut};
  int idFirstSurfacecut = gmshModelOccAddPlaneSurface(surfacescut, 1, -1, &ierr);
  ErrorGmsh(ierr);
  int surfaces_biscut[] = {linesloop_biscut};
  int idSecondSurfacecut = gmshModelOccAddPlaneSurface(surfaces_biscut, 1, -1, &ierr);
  ErrorGmsh(ierr);
  int surfaces_tercut[] = {linesloop_tercut};
  int idThirdSurfacecut = gmshModelOccAddPlaneSurface(surfaces_tercut, 1, -1, &ierr);
  ErrorGmsh(ierr);
  int surfaces_quatercut[] = {linesloop_quatercut};
  int idFourthSurfacecut = gmshModelOccAddPlaneSurface(surfaces_quatercut, 1, -1, &ierr);
  ErrorGmsh(ierr);

  // On découpe la plaque principale avec les triangles
  int base[] = {2,idPlate};
  int firstSurface[] = {2,idFirstSurface};
  int secondSurface[] = {2,idSecondSurface};
  int thirdSurface[] = {2,idThirdSurface};
  int fourthSurface[] = {2,idFourthSurface};
  gmshModelOccCut(base,2,firstSurface,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
  ErrorGmsh(ierr);
  gmshModelOccCut(base,2,secondSurface,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
  ErrorGmsh(ierr);
  gmshModelOccCut(base,2,thirdSurface,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
  ErrorGmsh(ierr);
  gmshModelOccCut(base,2,fourthSurface,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
  ErrorGmsh(ierr);

  // On découpe la plaque à decouper avec les triangles
  int basecut[] = {2,idPlatecut};
  int firstSurfacecut[] = {2,idFirstSurfacecut};
  int secondSurfacecut[] = {2,idSecondSurfacecut};
  int thirdSurfacecut[] = {2,idThirdSurfacecut};
  int fourthSurfacecut[] = {2,idFourthSurfacecut};
  gmshModelOccCut(basecut,2,firstSurfacecut,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
  ErrorGmsh(ierr);
  gmshModelOccCut(basecut,2,secondSurfacecut,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
  ErrorGmsh(ierr);
  gmshModelOccCut(basecut,2,thirdSurfacecut,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
  ErrorGmsh(ierr);
  gmshModelOccCut(basecut,2,fourthSurfacecut,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
  ErrorGmsh(ierr);
 
  //On délimite les surfaces pour simuler 2 matériaux différents
  int idFirstPointLimit = gmshModelOccAddPoint(-w/2.0, -h/2.0 + w/4.0, 0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  int idFirstPointLimitBis = gmshModelOccAddPoint(-w1 / 2.0, -h1 / 2.0 + w1 / 4.0, 0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  int idSecondPointLimit = gmshModelOccAddPoint(-w/2.0, h/2.0 - w/4.0, 0.0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  int idSecondPointLimitBis = gmshModelOccAddPoint(-w1/2.0, h1/2.0 - w1/4.0, 0.0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  int idThirdPointLimit = gmshModelOccAddPoint(w/2.0, h/2.0 - w/4.0, 0.0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  int idThirdPointLimitBis = gmshModelOccAddPoint(w1/2.0, h1/2.0 - w1/4.0, 0.0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  int idFourthPointLimit = gmshModelOccAddPoint(w/2.0, -h/2.0 + w/4.0, 0.0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  int idFourthPointLimitBis = gmshModelOccAddPoint(w1/2.0, -h1/2.0 + w1/4.0, 0.0, 5, -1, &ierr);
  ErrorGmsh(ierr);
  
  
  // On supprime les surfaces de la plaque principale
  gmshModelOccCut(base,2,basecut,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
  ErrorGmsh(ierr);


  // On synchronise le modèle
  geoSetSizeCallback(geoSize);            
  gmshModelOccSynchronize(&ierr);       
  gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
  gmshModelMeshGenerate(2, &ierr);

  /*
  int ierr;
  double r = w / 4;
  int idRect = gmshModelOccAddRectangle(0.0, 0.0, 0.0, w, h, -1, 0.0, &ierr);
  int idDisk = gmshModelOccAddDisk(w / 2.0, h / 2.0, 0.0, r, r, -1, NULL, 0, NULL, 0, &ierr);
  int idSlit = gmshModelOccAddRectangle(w / 2.0, h / 2.0 - r, 0.0, w, 2.0 * r, -1, 0.0, &ierr);
  int rect[] = {2, idRect};
  int disk[] = {2, idDisk};
  int slit[] = {2, idSlit};

  gmshModelOccCut(rect, 2, disk, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
  gmshModelOccCut(rect, 2, slit, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
  gmshModelOccSynchronize(&ierr);

  if (theGeometry->elementType == FEM_QUAD) {
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
    gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);
    gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);
    gmshModelGeoMeshSetRecombine(2, 1, 45, &ierr);
    gmshModelMeshGenerate(2, &ierr);
  }

  if (theGeometry->elementType == FEM_TRIANGLE) {
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);
  }*/
}
void geoMeshGenerateGeo(void) {
  femGeo *theGeometry = geoGetGeometry();
  double Lx = 1.0;
  double Ly = 1.0;
  theGeometry->LxPlate = Lx;
  theGeometry->LyPlate = Ly;
  theGeometry->h = Lx * 0.05;
  theGeometry->elementType = FEM_QUAD;

  geoSetSizeCallback(geoSize);

  /*
  4 ------------------ 3
  |                    |
  |                    |
  5 ------- 6          |
             \         |
              )        |
             /         |
  8 ------- 7          |
  |                    |
  |                    |
  1 ------------------ 2
  */

  int ierr;
  double w = theGeometry->LxPlate;
  double h = theGeometry->LyPlate;
  double r = w / 4;
  double lc = theGeometry->h;

  int p1 = gmshModelGeoAddPoint(-w / 2, -h / 2, 0., lc, 1, &ierr);
  int p2 = gmshModelGeoAddPoint(w / 2, -h / 2, 0., lc, 2, &ierr);
  int p3 = gmshModelGeoAddPoint(w / 2, h / 2, 0., lc, 3, &ierr);
  int p4 = gmshModelGeoAddPoint(-w / 2, h / 2, 0., lc, 4, &ierr);
  int p5 = gmshModelGeoAddPoint(-w / 2, r, 0., lc, 5, &ierr);
  int p6 = gmshModelGeoAddPoint(0., r, 0., lc, 6, &ierr);
  int p7 = gmshModelGeoAddPoint(0., -r, 0., lc, 7, &ierr);
  int p8 = gmshModelGeoAddPoint(-w / 2, -r, 0., lc, 8, &ierr);
  int p9 = gmshModelGeoAddPoint(0., 0., 0., lc, 9, &ierr); // center of circle

  int l1 = gmshModelGeoAddLine(p1, p2, 1, &ierr);
  int l2 = gmshModelGeoAddLine(p2, p3, 2, &ierr);
  int l3 = gmshModelGeoAddLine(p3, p4, 3, &ierr);
  int l4 = gmshModelGeoAddLine(p4, p5, 4, &ierr);
  int l5 = gmshModelGeoAddLine(p5, p6, 5, &ierr);
  int l6 = gmshModelGeoAddCircleArc(p7, p9, p6, 6, 0., 0., 0., &ierr); // NB : the direction of the curve is reversed
  int l7 = gmshModelGeoAddLine(p7, p8, 7, &ierr);
  int l8 = gmshModelGeoAddLine(p8, p1, 8, &ierr);

  int lTags[] = {l1, l2, l3, l4, l5, -l6, l7, l8}; // NB : "-l6" because the curve is reversed
  int c1[] = {1};
  c1[0] = gmshModelGeoAddCurveLoop(lTags, 8, 1, 0, &ierr);
  int s1 = gmshModelGeoAddPlaneSurface(c1, 1, 1, &ierr);
  gmshModelGeoSynchronize(&ierr);

  if (theGeometry->elementType == FEM_QUAD) {
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
    gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);
    gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);
    gmshModelGeoMeshSetRecombine(2, 1, 45, &ierr);
    gmshModelMeshGenerate(2, &ierr);
  }

  if (theGeometry->elementType == FEM_TRIANGLE) {
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);
  }

  //   gmshFltkRun(&ierr);
}

void geoMeshGenerateGeoFile(const char *filename) {
  femGeo *theGeometry = geoGetGeometry();
  int ierr;
  gmshOpen(filename, &ierr);
  ErrorGmsh(ierr);
  if (theGeometry->elementType == FEM_QUAD) {
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
    gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);
    gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);
    gmshModelGeoMeshSetRecombine(2, 1, 45, &ierr);
    gmshModelMeshGenerate(2, &ierr);
  }

  if (theGeometry->elementType == FEM_TRIANGLE) {
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);
  }
  return;
}

void geoMeshGenerateMshFile(const char *filename) {
  int ierr;
  gmshOpen(filename, &ierr);
  ErrorGmsh(ierr);
  return;
}
