/*
* Tsunami simulation - file with the algorithms for the simulation
*
* We use finite elements methods with Navier-Stokes equations for the model
*
* By Arbai Faycal and Sofyan Mimouni.
*
* With the contribution and work of Vincent Legat, Université Catholique de Louvain (UCL)
*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "GL/glfw.h"
#include <unistd.h>
#include <time.h>
#define Error(a)   femError(a,__LINE__,__FILE__)
#define Warning(a) femWarning(a,  __LINE__, __FILE__)

typedef enum {FEM_TRIANGLE,FEM_QUAD,FEM_EDGE} femElementType;
typedef enum {FEM_FULL,FEM_BAND,FEM_ITER} femSolverType;
typedef enum {FEM_NO,FEM_XNUM,FEM_YNUM} femRenumType;


typedef struct {
    int *elem;
    double *X;
    double *Y;
    double *H;
    int nElem;
    int nNode;
    int nLocalNode;
} femMesh;

typedef struct {
    int elem[2];
    int node[2];
} femEdge;

typedef struct {
    femMesh *mesh;
    femEdge *edges;
    int nEdge;
    int nBoundary;
} femEdges;

typedef struct {
    int n;
    int order;
    void (*x2)(double *xsi, double *eta);
    void (*phi2)(double xsi, double eta, double *phi);
    void (*dphi2dx)(double xsi, double eta, double *dphidxsi, double *dphideta);
    void (*x1)(double *xsi);
    void (*phi1)(double xsi, double *phi);
    void (*dphi1dx)(double xsi, double *dphidxsi);

} femDiscrete;

typedef struct {
    int n;
    const double *xsi;
    const double *eta;
    const double *weight;
} femIntegration;

typedef struct {
    femSolverType type;
    void *solver;
} femSolver;

typedef struct {
    double *B;
    double **A;
    int size;
} femFullSystem;

typedef struct
{
    double *B;
    double **A;
    int size;
    int band;
} femBandSystem;

typedef struct
{
    double *R;
    double *R0;
    double *D;
    double *D0;
    double error;
    int size;
    int iter;
} femIterativeSolver;

typedef struct {
    femMesh *mesh;
    femEdges *edges;
    femDiscrete *space;
    femIntegration *rule;
    femSolver *solver;
    int size;
    int sizeLocal[2];
    int *number;
    int *map;
    double *soluce;
    double *X;
    double *Y;
} femDiffusionProblem;

typedef struct {
    femMesh *mesh;
    femEdges *edges;
    femDiscrete *space;
    femIntegration *rule1d;
    femIntegration *rule2d;
    int size;
    double *C;
    double *U;
    double *V;
    double *F;
} femAdvProblem;

typedef struct {
    femMesh *mesh;
    femEdges *edges;
    femDiscrete *space;
    femIntegration *rule1d;
    femIntegration *rule2d;
    femSolver *solver;
    int size;


    //int *map;
    //double *X;
   // double *Y;

    double *E;
    double *U;
    double *V;
    double *FE;
    double *FU;
    double *FV;
    double f0;
    double gravity;
    double omega ;
    double R;
    double beta;
    double height;
    double rho;
    double gamma;
    double L;
    double tau0;
    double timeStep;
} femShallowProblem;


typedef struct{
    femMesh *mesh;
    femEdges *edges;
    femDiscrete *space;
    femIntegration *rule1d;
    femIntegration *rule2d;
    femSolver *solver;
    int size;
    //int *map;
    //double *X;
    //double *Y;
    double *H;
    double *U;
    double *V;
    double *FH;
    double *FU;
    double *FV;
    double f0;
    double gravity;
    double omega ;
    double R;
    double beta;
    double *height;
    double rho;
    double gamma;
    double L;
    double tau0;
    double timeStep;
    double *hEdge;
    double dt;
} femTsunamiProblem;

femTsunamiProblem   *femTsunamiCreate(const char *meshFileName);
void                 femTsunamiFree(femTsunamiProblem *myProblem);
void                 femTsunamiAddIntegralsElements(femTsunamiProblem *myProblem);
void                 femTsunamiAddIntegralsEdges(femTsunamiProblem *myProblem);
void                 femTsunamiMultiplyInverseMatrix(femTsunamiProblem *myProblem);
void                 femTsunamiCompute(femTsunamiProblem *myProblem);
void 		     tsunamiWriteFile(const char *baseResultName, int iter, double *U, double *V, double *E, int nelem, int nsub);
int 		     tsunamiReadFile(const char *baseResultName, int iter, double *U, double *V, double *E, int nelem);
double 		     tsunamiInitialConditionOkada(double x, double y);

femIntegration      *femIntegrationCreate(int n, femElementType type);
void                 femIntegrationFree(femIntegration *theRule);

femMesh             *femMeshRead(const char *filename);
void                 femMeshWrite(const femMesh* myMesh, const char *filename);
void                 femMeshFree(femMesh *theMesh);

femEdges*            femEdgesCreate(femMesh *theMesh);
void                 femEdgesFree(femEdges *theEdges);
void                 femEdgesPrint(femEdges *theEdges);
int                  femEdgesCompare(const void *edgeOne, const void *edgeTwo);
void                 femEdgesMap(femEdges *theEdges, int index, int map[2][2]);

femDiscrete*         femDiscreteCreate(int n, femElementType type);
void                 femDiscreteFree(femDiscrete* mySpace);
void                 femDiscretePrint(femDiscrete* mySpace);
void                 femDiscreteXsi2(femDiscrete* mySpace, double *xsi, double *eta);
void                 femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi);
void                 femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta);
void                 femDiscreteXsi1(femDiscrete* mySpace, double *xsi);
void                 femDiscretePhi1(femDiscrete* mySpace, double xsi, double *phi);
void                 femDiscreteDphi1(femDiscrete* mySpace, double xsi, double *dphidxsi);
double               femDiscreteInterpolate(double *phi, double *U, int *map, int n);

femSolver*           femSolverFullCreate(int size);
femSolver*           femSolverBandCreate(int size,int band);
femSolver*           femSolverIterativeCreate(int size);
void                 femSolverFree(femSolver* mySolver);
void                 femSolverInit(femSolver* mySolver);
void                 femSolverPrint(femSolver* mySolver);
void                 femSolverPrintInfos(femSolver* mySolver);
double*              femSolverEliminate(femSolver* mySolver);
void                 femSolverConstrain(femSolver* mySolver, int myNode, double value);
void                 femSolverAssemble(femSolver* mySolver, double *Aloc, double *Bloc, double *Uloc, int *map, int nLoc);
double               femSolverGet(femSolver* mySolver, int i, int j);
int                  femSolverConverged(femSolver *mySolver);

femFullSystem*       femFullSystemCreate(int size);
void                 femFullSystemFree(femFullSystem* mySystem);
void                 femFullSystemInit(femFullSystem* mySystem);
void                 femFullSystemPrint(femFullSystem* mySystem);
void                 femFullSystemPrintInfos(femFullSystem* mySystem);
double*              femFullSystemEliminate(femFullSystem* mySystem);
void                 femFullSystemConstrain(femFullSystem* mySystem, int myNode, double value);
void                 femFullSystemAssemble(femFullSystem* mySystem, double *Aloc, double *Bloc, int *map, int nLoc);
double               femFullSystemGet(femFullSystem* mySystem, int i, int j);

femBandSystem*       femBandSystemCreate(int size, int band);
void                 femBandSystemFree(femBandSystem* myBandSystem);
void                 femBandSystemInit(femBandSystem *myBand);
void                 femBandSystemPrint(femBandSystem *myBand);
void                 femBandSystemPrintInfos(femBandSystem *myBand);
double*              femBandSystemEliminate(femBandSystem *myBand);
void                 femBandSystemConstrain(femBandSystem *myBand, int myNode, double myValue);
void                 femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc);
double               femBandSystemGet(femBandSystem* myBandSystem, int i, int j);

femIterativeSolver*  femIterativeSolverCreate(int size);
void                 femIterativeSolverFree(femIterativeSolver* mySolver);
void                 femIterativeSolverInit(femIterativeSolver* mySolver);
void                 femIterativeSolverPrint(femIterativeSolver* mySolver);
void                 femIterativeSolverPrintInfos(femIterativeSolver* mySolver);
double*              femIterativeSolverEliminate(femIterativeSolver* mySolver);
void                 femIterativeSolverConstrain(femIterativeSolver* mySolver, int myNode, double myValue);
void                 femIterativeSolverAssemble(femIterativeSolver* mySolver, double *Aloc, double *Bloc, double *Uloc, int *map, int nLoc);
double               femIterativeSolverGet(femIterativeSolver* mySolver, int i, int j);
int                  femIterativeSolverConverged(femIterativeSolver *mySolver);

femDiffusionProblem *femDiffusionCreate(const char *filename, femSolverType solverType, femRenumType renumType, int orderType);
void                 femDiffusionFree(femDiffusionProblem *myProblem);
void                 femDiffusionMeshLocal(const femDiffusionProblem *myProblem, const int i, int *map, double *x, double *y, double *u);
void                 femDiffusionCompute(femDiffusionProblem *myProblem);
void                 femDiffusionRenumber(femDiffusionProblem *myProblem, femRenumType renumType);
int                  femDiffusionComputeBand(femDiffusionProblem *myProblem);
void                 femDiffusionComputeMap(femDiffusionProblem *myProblem, int orderType);
void                 femDiffusionPrintMap(femDiffusionProblem *myProblem);
void                 femDiffusionPrintInfos(femDiffusionProblem *myProblem);

femAdvProblem       *femAdvCreate(const char *meshFileName);
void                 femAdvFree(femAdvProblem *myProblem);
void                 femAdvInitialCondition(femAdvProblem *myProblem);
void                 femAdvAddIntegralsTriangles(femAdvProblem *myProblem);
void                 femAdvAddIntegralsEdges(femAdvProblem *myProblem);
void                 femAdvMultiplyInverseMatrix(femAdvProblem *myProblem);
void                 femAdvTriangleMap(femAdvProblem *myProblem, int i, int map[3]);
void                 femAdvEdgeMap(femAdvProblem *myProblem, int i, int map[2][2]);

femShallowProblem   *femShallowCreate(const char *meshFileName);
void                 femShallowFree(femShallowProblem *myProblem);
void                 femShallowAddIntegralsElements(femShallowProblem *myProblem);
void                 femShallowAddIntegralsEdges(femShallowProblem *myProblem);
void                 femShallowMultiplyInverseMatrix(femShallowProblem *myProblem);
void                 femShallowCompute(femShallowProblem *myProblem);
void                 femShallowWrite(femShallowProblem *myProblem,char *filename);
void                 femShallowRead(femShallowProblem *myProblem,char *filename);

void                 femStommel(double x, double y, double *u, double *v, double *eta);

double               femMin(double *x, int n);
double               femMax(double *x, int n);
void                 femError(char *text, int line, char *file);
void                 femWarning(char *text, int line, char *file);
double               tsunamiInitialConditionOkada(double x, double y);

static const double _gaussEdge2Xsi[2]     = { 0.577350269189626,-0.577350269189626 };
static const double _gaussEdge2Weight[2]  = { 1.000000000000000, 1.000000000000000 };

static const double _gaussQuad4Xsi[4]     = {-0.577350269189626,-0.577350269189626, 0.577350269189626, 0.577350269189626 };
static const double _gaussQuad4Eta[4]     = { 0.577350269189626,-0.577350269189626,-0.577350269189626, 0.577350269189626 };
static const double _gaussQuad4Weight[4]  = { 1.000000000000000, 1.000000000000000, 1.000000000000000, 1.000000000000000 };
static const double _gaussQuad9Xsi[9]     = { 0.774596669241483, 0.000000000000000,-0.774596669241483,
                                              0.774596669241483, 0.000000000000000,-0.774596669241483,
                                              0.774596669241483, 0.000000000000000,-0.774596669241483 };
static const double _gaussQuad9Eta[9]     = { 0.774596669241483, 0.774596669241483, 0.774596669241483,
                                              0.000000000000000, 0.000000000000000, 0.000000000000000,
                                             -0.774596669241483,-0.774596669241483,-0.774596669241483 };
static const double _gaussQuad9Weight[9]  = { 0.308641975308642, 0.493827160493827, 0.308641975308642,
                                              0.493827160493827, 0.790123456790123, 0.493827160493827,
                                              0.308641975308642, 0.493827160493827, 0.308641975308642 };

static const double _gaussTri3Xsi[3]      = { 0.166666666666667, 0.666666666666667, 0.166666666666667 };
static const double _gaussTri3Eta[3]      = { 0.166666666666667, 0.166666666666667, 0.666666666666667 };
static const double _gaussTri3Weight[3]   = { 0.166666666666667, 0.166666666666667, 0.166666666666667 };
static const double _gaussTri12Xsi[12]    = { 0.249286745170910, 0.249286745170910, 0.501426509658179,
                                              0.063089014491502, 0.063089014491502, 0.873821971016996,
                                              0.310352451033785, 0.636502499121399, 0.053145049844816,
                                              0.310352451033785, 0.636502499121399, 0.053145049844816 };
static const double _gaussTri12Eta[12]    = { 0.249286745170910, 0.501426509658179, 0.249286745170910,
                                              0.063089014491502, 0.873821971016996, 0.063089014491502,
                                              0.636502499121399, 0.053145049844816, 0.310352451033785,
                                              0.053145049844816, 0.310352451033785, 0.636502499121399 };
static const double _gaussTri12Weight[12] = { 0.058393137863189, 0.058393137863189, 0.058393137863189,
                                              0.025422453185104, 0.025422453185104, 0.025422453185104,
                                              0.041425537809187, 0.041425537809187, 0.041425537809187,
                                              0.041425537809187, 0.041425537809187, 0.041425537809187 };






femIntegration *femIntegrationCreate(int n, femElementType type)
{
    femIntegration *theRule = malloc(sizeof(femIntegration));
    if (type == FEM_QUAD && n == 4) {
        theRule->n      = 4;
        theRule->xsi    = _gaussQuad4Xsi;
        theRule->eta    = _gaussQuad4Eta;
        theRule->weight = _gaussQuad4Weight; }
    else if (type == FEM_QUAD && n == 9) {
        theRule->n      = 9;
        theRule->xsi    = _gaussQuad9Xsi;
        theRule->eta    = _gaussQuad9Eta;
        theRule->weight = _gaussQuad9Weight; }
    else if (type == FEM_TRIANGLE && n == 3) {
        theRule->n      = 3;
        theRule->xsi    = _gaussTri3Xsi;
        theRule->eta    = _gaussTri3Eta;
        theRule->weight = _gaussTri3Weight; }
    else if (type == FEM_TRIANGLE && n == 12) {
        theRule->n      = 12;
        theRule->xsi    = _gaussTri12Xsi;
        theRule->eta    = _gaussTri12Eta;
        theRule->weight = _gaussTri12Weight; }
    else if (type == FEM_EDGE && n == 2) {
        theRule->n      = 2;
        theRule->xsi    = _gaussEdge2Xsi;
        theRule->eta    = NULL;
        theRule->weight = _gaussEdge2Weight; }

    else Error("Cannot create such an integration rule !");
    return theRule;
}

void femIntegrationFree(femIntegration *theRule)
{
    free(theRule);
}


void _1c0_x(double *xsi)
{
    xsi[0] = -1.0;
    xsi[1] =  1.0;
}


void _1c0_phi(double xsi, double *phi)
{
    phi[0] = (1.0 - xsi)/2.0;
    phi[1] = (1.0 + xsi)/2.0;

}

void _1c0_dphidx(double xsi, double *dphidxsi)
{
    dphidxsi[0] =  -1.0/2.0;
    dphidxsi[1] =   1.0/2.0;
}

void _3c0_x(double *xsi)
{
    xsi[0] = -1.0;
    xsi[1] = -1.0/3.0;
    xsi[1] =  1.0/3.0;
    xsi[1] =  1.0;
}

void _3c0_phi(double xsi, double *phi)
{
    phi[0] =   9./16 * (-1./3 - xsi) * ( 1./3 - xsi) * (1.   - xsi);
    phi[1] = -27./16 * (-1.   - xsi) * ( 1./3 - xsi) * (1.   - xsi);
    phi[2] =  27./16 * (-1.   - xsi) * (-1./3 - xsi) * (1.   - xsi);
    phi[3] =  -9./16 * (-1.   - xsi) * (-1./3 - xsi) * (1./3 - xsi);
}

void _3c0_dphidx(double xsi, double *dphidxsi)
{
    dphidxsi[0] =   9./16 * ( - ( 1./3 - xsi) * (1.   - xsi) - (-1./3 - xsi) * (1.   - xsi) - (-1./3 - xsi) * ( 1./3 - xsi) ) ;
    dphidxsi[1] = -27./16 * ( - ( 1./3 - xsi) * (1.   - xsi) - (-1.   - xsi) * (1.   - xsi) - (-1.   - xsi) * ( 1./3 - xsi) ) ;
    dphidxsi[2] =  27./16 * ( - (-1./3 - xsi) * (1.   - xsi) - (-1.   - xsi) * (1.   - xsi) - (-1.   - xsi) * (-1./3 - xsi) ) ;
    dphidxsi[3] =  -9./16 * ( - (-1./3 - xsi) * (1./3 - xsi) - (-1.   - xsi) * (1./3 - xsi) - (-1.   - xsi) * (-1./3 - xsi) ) ;
}



void _q1c0_x(double *xsi, double *eta)
{
    xsi[0] =  1.0;  eta[0] =  1.0;
    xsi[1] = -1.0;  eta[1] =  1.0;
    xsi[2] = -1.0;  eta[2] = -1.0;
    xsi[3] =  1.0;  eta[3] = -1.0;
}

void _q1c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = (1.0 + xsi) * (1.0 + eta) / 4.0;
    phi[1] = (1.0 - xsi) * (1.0 + eta) / 4.0;
    phi[2] = (1.0 - xsi) * (1.0 - eta) / 4.0;
    phi[3] = (1.0 + xsi) * (1.0 - eta) / 4.0;
}

void _q1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] =   (1.0 + eta) / 4.0;
    dphidxsi[1] = - (1.0 + eta) / 4.0;
    dphidxsi[2] = - (1.0 - eta) / 4.0;
    dphidxsi[3] =   (1.0 - eta) / 4.0;
    dphideta[0] =   (1.0 + xsi) / 4.0;
    dphideta[1] =   (1.0 - xsi) / 4.0;
    dphideta[2] = - (1.0 - xsi) / 4.0;
    dphideta[3] = - (1.0 + xsi) / 4.0;

}

void _q2c0_x(double *xsi, double *eta)
{
    xsi[0] =  1.0;  eta[0] =  1.0;
    xsi[1] = -1.0;  eta[1] =  1.0;
    xsi[2] = -1.0;  eta[2] = -1.0;
    xsi[3] =  1.0;  eta[3] = -1.0;
    xsi[4] =  0.0;  eta[4] =  1.0;
    xsi[5] = -1.0;  eta[5] =  0.0;
    xsi[6] =  0.0;  eta[6] = -1.0;
    xsi[7] =  1.0;  eta[7] =  0.0;
    xsi[8] =  0.0;  eta[8] =  0.0;
}

void _q2c0_phi(double xsi, double eta, double *phi)
{
    phi[0] =  xsi*(1.0+xsi)*eta*(1.0+eta)/4.0;
    phi[1] = -xsi*(1.0-xsi)*eta*(1.0+eta)/4.0;
    phi[2] =  xsi*(1.0-xsi)*eta*(1.0-eta)/4.0;
    phi[3] = -xsi*(1.0+xsi)*eta*(1.0-eta)/4.0;
    phi[4] =  (1.0-xsi*xsi)*eta*(1.0+eta)/2.0;
    phi[5] = -xsi*(1.0-xsi)*(1.0-eta*eta)/2.0;
    phi[6] = -(1.0-xsi*xsi)*eta*(1.0-eta)/2.0;
    phi[7] =  xsi*(1.0+xsi)*(1.0-eta*eta)/2.0;
    phi[8] =  (1.0-xsi*xsi)*(1.0-eta*eta);
}

void _q2c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] =  (1.0+2.0*xsi)*eta*(1.0+eta)/4.0;
    dphidxsi[1] = (-1.0+2.0*xsi)*eta*(1.0+eta)/4.0;
    dphidxsi[2] =  (1.0-2.0*xsi)*eta*(1.0-eta)/4.0;
    dphidxsi[3] = -(1.0+2.0*xsi)*eta*(1.0-eta)/4.0;
    dphidxsi[4] =       -2.0*xsi*eta*(1.0+eta)/2.0;
    dphidxsi[5] = (-1.0+2.0*xsi)*(1.0-eta*eta)/2.0;
    dphidxsi[6] =        2.0*xsi*eta*(1.0-eta)/2.0;
    dphidxsi[7] =  (1.0+2.0*xsi)*(1.0-eta*eta)/2.0;
    dphidxsi[8] =       -2.0*xsi*(1.0-eta*eta);
    dphideta[0] =  xsi*(1.0+xsi)*(1.0+2.0*eta)/4.0;
    dphideta[1] = -xsi*(1.0-xsi)*(1.0+2.0*eta)/4.0;
    dphideta[2] =  xsi*(1.0-xsi)*(1.0-2.0*eta)/4.0;
    dphideta[3] = -xsi*(1.0+xsi)*(1.0-2.0*eta)/4.0;
    dphideta[4] =  (1.0-xsi*xsi)*(1.0+2.0*eta)/2.0;
    dphideta[5] =  xsi*(1.0-xsi)*2.0*eta/2.0;
    dphideta[6] = -(1.0-xsi*xsi)*(1.0-2.0*eta)/2.0;
    dphideta[7] =  xsi*(1.0+xsi)*(-2.0*eta)/2.0;
    dphideta[8] =  (1.0-xsi*xsi)*(-2.0*eta);
}

void _q3c0_x(double *xsi, double *eta)
{
    xsi[0] =  1.0;       eta[0] =  1.0;
    xsi[1] = -1.0;       eta[1] =  1.0;
    xsi[2] = -1.0;       eta[2] = -1.0;
    xsi[3] =  1.0;       eta[3] = -1.0;
    xsi[4] =  1.0/3.0;   eta[4] =  1.0;
    xsi[5] = -1.0/3.0;   eta[5] =  1.0;
    xsi[6] = -1.0;       eta[6] =  1.0/3.0;
    xsi[7] = -1.0;       eta[7] = -1.0/3.0;
    xsi[8] = -1.0/3.0;   eta[8] = -1.0;
    xsi[9] =  1.0/3.0;   eta[9] = -1.0;
    xsi[10] =  1.0;      eta[10] = -1.0/3.0;
    xsi[11] =  1.0;      eta[11] =  1.0/3.0;
    xsi[12] =  1.0/3.0;  eta[12] =  1.0/3.0;
    xsi[13] = -1.0/3.0;  eta[13] =  1.0/3.0;
    xsi[14] = -1.0/3.0;  eta[14] = -1.0/3.0;
    xsi[15] =  1.0/3.0;  eta[15] = -1.0/3.0;
}


void _q3c0_phi(double xsi, double eta, double *phi)
{
    double fXi[4], fEta[4];
    int map[16] = {2,8,9,3,7,14,15,10,6,13,12,11,1,5,4,0};
    _3c0_phi(xsi, fXi);
    _3c0_phi(eta, fEta);
    int i, j;
    for (i = 0; i < 4; ++i)
        for (j = 0; j < 4; ++j)
            phi[map[j*4+i]] = fXi[i] * fEta[j];
}

void _q3c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{

    double fXi[4], fEta[4], dfXi[4], dfEta[4];
    int map[16] = {2,8,9,3,7,14,15,10,6,13,12,11,1,5,4,0};

    _3c0_phi(xsi, fXi);
    _3c0_phi(eta, fEta);
    _3c0_dphidx(xsi, dfXi);
    _3c0_dphidx(eta, dfEta);
    int i, j;
    for (i = 0; i < 4; ++i) {
        for (j = 0; j < 4; ++j) {
            dphidxsi[map[j*4+i]] = dfXi[i] * fEta[j];
            dphideta[map[j*4+i]] = fXi[i] * dfEta[j]; }}
}


void _p1c0_x(double *xsi, double *eta)
{
    xsi[0] =  0.0;  eta[0] =  0.0;
    xsi[1] =  1.0;  eta[1] =  0.0;
    xsi[2] =  0.0;  eta[2] =  1.0;
}

void _p1c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = 1 - xsi - eta;
    phi[1] = xsi;
    phi[2] = eta;
}

void _p1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] = -1.0;
    dphidxsi[1] =  1.0;
    dphidxsi[2] =  0.0;
    dphideta[0] = -1.0;
    dphideta[1] =  0.0;
    dphideta[2] =  1.0;

}

void _p2c0_x(double *xsi, double *eta)
{
    xsi[0] =  0.0;  eta[0] =  0.0;
    xsi[1] =  1.0;  eta[1] =  0.0;
    xsi[2] =  0.0;  eta[2] =  1.0;
    xsi[3] =  0.5;  eta[3] =  0.0;
    xsi[4] =  0.5;  eta[4] =  0.5;
    xsi[5] =  0.0;  eta[5] =  0.5;
}


void _p2c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = 1.0 - 3.0*(xsi+eta) + 2.0*(xsi+eta)*(xsi+eta);
    phi[1] = xsi*(2.0*xsi-1.0);
    phi[2] = eta*(2.0*eta-1.0);
    phi[3] = 4.0*xsi*(1.0-xsi-eta);
    phi[4] = 4.0*xsi*eta;
    phi[5] = 4.0*eta*(1.0-xsi-eta);
}

void _p2c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] = - 3.0 + 4.0*xsi + 4.0*eta;
    dphidxsi[1] = - 1.0 + 4.0*xsi          ;
    dphidxsi[2] =   0.0                    ;
    dphidxsi[3] =   4.0 - 8.0*xsi - 4.0*eta;
    dphidxsi[4] =                   4.0*eta;
    dphidxsi[5] =                 - 4.0*eta;
    dphideta[0] = - 3.0 + 4.0*xsi + 4.0*eta;
    dphideta[1] =   0.0                    ;
    dphideta[2] =  -1.0           + 4.0*eta;
    dphideta[3] =       - 4.0*xsi          ;
    dphideta[4] =         4.0*xsi          ;
    dphideta[5] =   4.0 - 4.0*xsi - 8.0*eta;
}

void _p3c0_x(double *xsi, double *eta)
{
    xsi[0] = 0.0;     eta[0] = 0.0;
    xsi[1] = 1.0;     eta[1] = 0.0;
    xsi[2] = 0.0;     eta[2] = 1.0;
    xsi[3] = 1.0/3.0; eta[3] = 0.0;
    xsi[4] = 2.0/3.0; eta[4] = 0.0;
    xsi[5] = 2.0/3.0; eta[5] = 1.0/3.0;
    xsi[6] = 1.0/3.0; eta[6] = 2.0/3.0;
    xsi[7] = 0.0;     eta[7] = 2.0/3.0;
    xsi[8] = 0.0;     eta[8] = 1.0/3.0;
    xsi[9] = 1.0/3.0; eta[9] = 1.0/3.0;
}

void _p3c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = (1.0 - xsi - eta) * (1.0/3.0 - xsi - eta) * (2.0/3.0 - xsi - eta) *  9.0/2.0;
    phi[1] =               xsi * (xsi - 1.0/3.0)       * (xsi - 2.0/3.0)       *  9.0/2.0;
    phi[2] =               eta * (eta - 1.0/3.0)       * (eta - 2.0/3.0)       *  9.0/2.0;
    phi[3] = (1.0 - xsi - eta) * (2.0/3.0 - xsi - eta) * xsi                   * 27.0/2.0;
    phi[4] = (1.0 - xsi - eta) * (xsi - 1.0/3.0)       * xsi                   * 27.0/2.0;
    phi[5] =               xsi * eta                   * (xsi - 1.0/3.0)       * 27.0/2.0;
    phi[6] =               xsi * eta                   * (eta - 1.0/3.0)       * 27.0/2.0;
    phi[7] = (1.0 - xsi - eta) * (eta - 1.0/3.0)       * eta                   * 27.0/2.0;
    phi[8] = (1.0 - xsi - eta) * (2.0/3.0 - xsi - eta) * eta                   * 27.0/2.0;
    phi[9] = (1.0 - xsi - eta) * xsi                   * eta                   * 27.0;
}

void _p3c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] = (- (1.0/3.0 - xsi - eta) * (2.0/3.0 - xsi - eta) - (1.0 - xsi - eta) * (2.0/3.0 - xsi - eta) - (1.0 - xsi - eta) * (1.0/3.0 - xsi - eta)) *  9.0/2.0;
    dphidxsi[1] = (  (xsi - 1.0/3.0) * (xsi - 2.0/3.0) +  xsi * (xsi - 2.0/3.0) + xsi * (xsi - 1.0/3.0) ) *  9.0/2.0;
    dphidxsi[2] = 0.0;
    dphidxsi[3] = ( - (2.0/3.0 - xsi - eta) * xsi - (1.0 - xsi - eta) * xsi + (1.0 - xsi - eta) * (2.0/3.0 - xsi - eta) ) * 27.0/2.0;
    dphidxsi[4] = ( - (xsi - 1.0/3.0) * xsi + (1.0 - xsi - eta) * xsi + (1.0 - xsi - eta) * (xsi - 1.0/3.0) ) * 27.0/2.0;
    dphidxsi[5] = ( eta * (xsi - 1.0/3.0) + xsi * eta ) * 27.0/2.0;
    dphidxsi[6] =   eta * (eta - 1.0/3.0) * 27.0/2.0;
    dphidxsi[7] = - (eta - 1.0/3.0) * eta * 27.0/2.0;
    dphidxsi[8] = ( - (2.0/3.0 - xsi - eta) * eta - (1.0 - xsi - eta) * eta  ) * 27.0/2.0;
    dphidxsi[9] = ( - xsi * eta + (1.0 - xsi - eta) * eta) * 27.0;
    dphideta[0] = (- (1.0/3.0 - xsi - eta) * (2.0/3.0 - xsi - eta) - (1.0 - xsi - eta) * (2.0/3.0 - xsi - eta) - (1.0 - xsi - eta) * (1.0/3.0 - xsi - eta)) *  9.0/2.0;
    dphideta[1] = 0.0;
    dphideta[2] =  (  (eta - 1.0/3.0) * (eta - 2.0/3.0) +  eta * (eta - 2.0/3.0) + eta * (eta - 1.0/3.0) ) *  9.0/2.0;
    dphideta[3] = ( - (2.0/3.0 - xsi - eta) * xsi - (1.0 - xsi - eta) * xsi ) * 27.0/2.0;
    dphideta[4] = - (xsi - 1.0/3.0)  * xsi * 27.0/2.0;
    dphideta[5] =   xsi * (xsi - 1.0/3.0) * 27.0/2.0;;
    dphideta[6] = ( xsi * (eta - 1.0/3.0) +  xsi * eta ) * 27.0/2.0;
    dphideta[7] = (- (eta - 1.0/3.0) * eta + (1.0 - xsi - eta) * eta + (1.0 - xsi - eta) * (eta - 1.0/3.0) ) * 27.0/2.0;
    dphideta[8] = ( - (2.0/3.0 - xsi - eta) * eta - (1.0 - xsi - eta) * eta + (1.0 - xsi - eta) * (2.0/3.0 - xsi - eta)) * 27.0/2.0;
    dphideta[9] = ( - xsi * eta + (1.0 - xsi - eta) * xsi) * 27.0;

}



femDiscrete *femDiscreteCreate(int n, femElementType type)
{
    femDiscrete *theSpace = malloc(sizeof(femDiscrete));
    if (type == FEM_QUAD && n == 4) {
        theSpace->n       = 4;
        theSpace->order   = 1;
        theSpace->x2      = _q1c0_x;
        theSpace->phi2    = _q1c0_phi;
        theSpace->dphi2dx = _q1c0_dphidx;
        theSpace->x1      = _1c0_x;
        theSpace->phi1    = _1c0_phi;
        theSpace->dphi1dx = _1c0_dphidx;}
    else if (type == FEM_QUAD && n == 9) {
        theSpace->n       = 9;
        theSpace->order   = 2;
        theSpace->x2      = _q2c0_x;
        theSpace->phi2    = _q2c0_phi;
        theSpace->dphi2dx = _q2c0_dphidx; }
    else if (type == FEM_QUAD && n == 16) {
        theSpace->n       = 16;
        theSpace->order   = 3;
        theSpace->x2      = _q3c0_x;
        theSpace->phi2    = _q3c0_phi;
        theSpace->dphi2dx = _q3c0_dphidx;
        theSpace->x1      = _3c0_x;
        theSpace->phi1    = _3c0_phi;
        theSpace->dphi1dx = _3c0_dphidx; }
    else if (type == FEM_TRIANGLE && n == 3) {
        theSpace->n       = 3;
        theSpace->order   = 1;
        theSpace->x2      = _p1c0_x;
        theSpace->phi2    = _p1c0_phi;
        theSpace->dphi2dx = _p1c0_dphidx;
        theSpace->x1      = _1c0_x;
        theSpace->phi1    = _1c0_phi;
        theSpace->dphi1dx = _1c0_dphidx; }
    else if (type == FEM_TRIANGLE && n == 6) {
        theSpace->n       = 6;
        theSpace->order   = 2;
        theSpace->x2      = _p2c0_x;
        theSpace->phi2    = _p2c0_phi;
        theSpace->dphi2dx = _p2c0_dphidx; }
    else if (type == FEM_TRIANGLE && n == 10) {
        theSpace->n       = 10;
        theSpace->order   = 3;
        theSpace->x2      = _p3c0_x;
        theSpace->phi2    = _p3c0_phi;
        theSpace->dphi2dx = _p3c0_dphidx;
        theSpace->x1      = _3c0_x;
        theSpace->phi1    = _3c0_phi;
        theSpace->dphi1dx = _3c0_dphidx;}
    else Error("Cannot create such a discrete space !");
    return theSpace;
}

void femDiscreteFree(femDiscrete *theSpace)
{
    free(theSpace);
}

void femDiscreteXsi2(femDiscrete* mySpace, double *xsi, double *eta)
{
    mySpace->x2(xsi,eta);
}

void femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi)
{
    mySpace->phi2(xsi,eta,phi);
}

void femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta)
{
    mySpace->dphi2dx(xsi,eta,dphidxsi,dphideta);
}

void femDiscreteXsi1(femDiscrete* mySpace, double *xsi)
{
    mySpace->x1(xsi);
}

void femDiscretePhi1(femDiscrete* mySpace, double xsi, double *phi)
{
    mySpace->phi1(xsi,phi);
}

void femDiscreteDphi1(femDiscrete* mySpace, double xsi, double *dphidxsi)
{
    mySpace->dphi1dx(xsi,dphidxsi);
}

void femDiscretePrint(femDiscrete *mySpace)
{
    int i,j;
    int n = mySpace->n;
    double xsi[16], eta[16], phi[16], dphidxsi[16], dphideta[16];

    femDiscreteXsi2(mySpace,xsi,eta);
    for (i=0; i < n; i++) {

        femDiscretePhi2(mySpace,xsi[i],eta[i],phi);
        femDiscreteDphi2(mySpace,xsi[i],eta[i],dphidxsi,dphideta);

        for (j=0; j < n; j++)  {
            printf("(xsi=%+.1f,eta=%+.1f) : ",xsi[i],eta[i]);
            printf(" phi(%d)=%+.1f",j,phi[j]);
            printf("   dphidxsi(%d)=%+.1f",j,dphidxsi[j]);
            printf("   dphideta(%d)=%+.1f \n",j,dphideta[j]);  }
        printf(" \n"); }
}


double femDiscreteInterpolate(double *phi, double *U, int *map, int n)
{
    double u = 0.0; int i;
    for (i=0; i <n; i++)
        u += phi[i]*U[map[i]];
    return u;
}

femMesh *femMeshRead(const char *filename)
{
    femMesh *theMesh = malloc(sizeof(femMesh));

    int i,trash,*elem;

    FILE* file = fopen(filename,"r");
    if (file == NULL) Error("No mesh file !");

    fscanf(file, "Number of nodes %d \n", &theMesh->nNode);
    theMesh->X = malloc(sizeof(double)*theMesh->nNode);
    theMesh->Y = malloc(sizeof(double)*theMesh->nNode);
    theMesh->H = malloc(sizeof(double)*theMesh->nNode);
    for (i = 0; i < theMesh->nNode; ++i) {
        fscanf(file,"%d : %le %le %le \n",&trash,&theMesh->X[i],&theMesh->Y[i],&theMesh->H[i]); }

    char str[256]; fgets(str, sizeof(str), file);
    if (!strncmp(str,"Number of triangles",19))  {
        sscanf(str,"Number of triangles %d \n", &theMesh->nElem);
        theMesh->elem = malloc(sizeof(int)*3*theMesh->nElem);
        theMesh->nLocalNode = 3;
        for (i = 0; i < theMesh->nElem; ++i) {
            elem = &(theMesh->elem[i*3]);
            fscanf(file,"%d : %d %d %d\n", &trash,&elem[0],&elem[1],&elem[2]); }}
    else if (!strncmp(str,"Number of quads",15))  {
        sscanf(str,"Number of quads %d \n", &theMesh->nElem);
        theMesh->elem = malloc(sizeof(int)*4*theMesh->nElem);
        theMesh->nLocalNode = 4;
        for (i = 0; i < theMesh->nElem; ++i) {
            elem = &(theMesh->elem[i*4]);
        fscanf(file,"%d : %d %d %d %d\n", &trash,&elem[0],&elem[1],&elem[2],&elem[3]); }}

    fclose(file);
    femMeshWrite(theMesh,"toto.txt");
     return theMesh;
}



void femMeshFree(femMesh *theMesh)
{
    free(theMesh->X);
    free(theMesh->Y);
    free(theMesh->H);
    free(theMesh->elem);
    free(theMesh);
}

void femMeshWrite(const femMesh *theMesh, const char *filename)
{
    int i,*elem;

    FILE* file = fopen(filename,"w");

    fprintf(file, "Number of nodes %d \n", theMesh->nNode);
    for (i = 0; i < theMesh->nNode; ++i) {
        fprintf(file,"%6d : %14.7e %14.7e \n",i,(theMesh->X[i]+1.0)/2.0,(theMesh->Y[i]+1.0)/2.0); }

    if (theMesh->nLocalNode == 4) {
        fprintf(file, "Number of quads %d \n", theMesh->nElem);
        for (i = 0; i < theMesh->nElem; ++i) {
            elem = &(theMesh->elem[i*4]);
            fprintf(file,"%6d : %6d %6d %6d %6d \n", i,elem[0],elem[1],elem[2],elem[3]);   }}
    else if (theMesh->nLocalNode == 3) {
        fprintf(file, "Number of triangles %d \n", theMesh->nElem);
        for (i = 0; i < theMesh->nElem; ++i) {
            elem = &(theMesh->elem[i*3]);
            fprintf(file,"%6d : %6d %6d %6d \n", i,elem[0],elem[1],elem[2]);   }}

    fclose(file);
}

femEdges *femEdgesCreate(femMesh *theMesh)
{
    femEdges *theEdges = malloc(sizeof(femEdges));
    int nLoc = theMesh->nLocalNode;
    int i,j,n = theMesh->nElem * nLoc;
    femEdge* edges = malloc(n * sizeof(femEdge));
    theEdges->mesh  = theMesh;
    theEdges->edges = edges;
    theEdges->nEdge = n;
    theEdges->nBoundary = n;

    for (i = 0; i < theMesh->nElem; i++) {
        int *elem = &(theMesh->elem[i*nLoc]);
        for (j = 0; j < nLoc; j++) {
            int id = i * nLoc + j;
            edges[id].elem[0] = i;
            edges[id].elem[1] = -1;
            edges[id].node[0] = elem[j];
            edges[id].node[1] = elem[(j + 1) % nLoc]; }}

    qsort(theEdges->edges, theEdges->nEdge, sizeof(femEdge), femEdgesCompare);

    int index = 0;
    int nBoundary = 0;

    for (i=0; i < theEdges->nEdge; i++) {
      if (i == theEdges->nEdge - 1 || femEdgesCompare(&edges[i],&edges[i+1]) != 0) {
              edges[index] = edges[i];
              nBoundary++; }
      else {  edges[index] = edges[i];
              edges[index].elem[1] = edges[i+1].elem[0];
              i = i+1;}
      index++; }

    theEdges->edges = realloc(edges, index * sizeof(femEdge));
    theEdges->nEdge = index;
    theEdges->nBoundary = nBoundary;
    return theEdges;
}

void femEdgesPrint(femEdges *theEdges)
{
    int i;
    for (i = 0; i < theEdges->nEdge; ++i) {
        printf("%6d : %4d %4d : %4d %4d \n",i,
               theEdges->edges[i].node[0],theEdges->edges[i].node[1],
               theEdges->edges[i].elem[0],theEdges->edges[i].elem[1]); }
}

void femEdgesFree(femEdges *theEdges)
{
    free(theEdges->edges);
    free(theEdges);
}

int femEdgesCompare(const void *edgeOne, const void *edgeTwo)
{
    int *nodeOne = ((femEdge*) edgeOne)->node;
    int *nodeTwo = ((femEdge*) edgeTwo)->node;
    int  diffMin = fmin(nodeOne[0],nodeOne[1]) - fmin(nodeTwo[0],nodeTwo[1]);
    int  diffMax = fmax(nodeOne[0],nodeOne[1]) - fmax(nodeTwo[0],nodeTwo[1]);

    if (diffMin < 0)    return  1;
    if (diffMin > 0)    return -1;
    if (diffMax < 0)    return  1;
    if (diffMax > 0)    return -1;
                        return  0;
}

void femEdgesMap(femEdges *theEdges, int index, int map[2][2])
{
    int i,j,k;
    int n = theEdges->mesh->nLocalNode;

    for (j=0; j < 2; ++j) {
        int node = theEdges->edges[index].node[j];
        for (k=0; k < 2; k++) {
            int elem = theEdges->edges[index].elem[k];
            map[k][j] = (theEdges->mesh->nElem)*n;
            if (elem >= 0) {
                for (i=0; i < n; i++) {
                    if (theEdges->mesh->elem[elem*n + i] == node) {
                        map[k][j] = elem*n + i;  }}}}}
}


femSolver *femSolverFullCreate(int size)
{
    femSolver *mySolver = malloc(sizeof(femSolver));
    mySolver->type = FEM_FULL;
    mySolver->solver = (femSolver *)femFullSystemCreate(size);
    return(mySolver);
}

femSolver *femSolverBandCreate(int size, int band)
{
    femSolver *mySolver = malloc(sizeof(femSolver));
    mySolver->type = FEM_BAND;
    mySolver->solver = (femSolver *)femBandSystemCreate(size,band);
    return(mySolver);
}

femSolver *femSolverIterativeCreate(int size)
{
    femSolver *mySolver = malloc(sizeof(femSolver));
    mySolver->type = FEM_ITER;
    mySolver->solver = (femSolver *)femIterativeSolverCreate(size);
    return(mySolver);
}

void femSolverFree(femSolver *mySolver)
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemFree((femFullSystem *)mySolver->solver); break;
        case FEM_BAND : femBandSystemFree((femBandSystem *)mySolver->solver); break;
        case FEM_ITER : femIterativeSolverFree((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
    free(mySolver);
}

void femSolverInit(femSolver *mySolver)
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemInit((femFullSystem *)mySolver->solver); break;
        case FEM_BAND : femBandSystemInit((femBandSystem *)mySolver->solver); break;
        case FEM_ITER : femIterativeSolverInit((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
}


double femSolverGet(femSolver *mySolver,int i,int j)
{
    double value = 0;
    switch (mySolver->type) {
        case FEM_FULL : value = femFullSystemGet((femFullSystem *)mySolver->solver,i,j); break;
        case FEM_BAND : value = femBandSystemGet((femBandSystem *)mySolver->solver,i,j); break;
        case FEM_ITER : value = (i==j); break;
        default : Error("Unexpected solver type"); }
    return(value);
}

void femSolverPrint(femSolver *mySolver)
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemPrint((femFullSystem *)mySolver->solver); break;
        case FEM_BAND : femBandSystemPrint((femBandSystem *)mySolver->solver); break;
        case FEM_ITER : femIterativeSolverPrint((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
}

void femSolverPrintInfos(femSolver *mySolver)
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemPrintInfos((femFullSystem *)mySolver->solver); break;
        case FEM_BAND : femBandSystemPrintInfos((femBandSystem *)mySolver->solver); break;
        case FEM_ITER : femIterativeSolverPrintInfos((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
}


void femSolverAssemble(femSolver* mySolver, double *Aloc, double *Bloc, double *Uloc,int *map, int nLoc)
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemAssemble((femFullSystem *)mySolver->solver,Aloc,Bloc,map,nLoc); break;
        case FEM_BAND : femBandSystemAssemble((femBandSystem *)mySolver->solver,Aloc,Bloc,map,nLoc); break;
        case FEM_ITER : femIterativeSolverAssemble((femIterativeSolver *)mySolver->solver,Aloc,Bloc,Uloc,map,nLoc); break;
        default : Error("Unexpected solver type"); }
}

void femSolverConstrain(femSolver *mySolver, int myNode, double myValue)
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemConstrain((femFullSystem *)mySolver->solver,myNode,myValue); break;
        case FEM_BAND : femBandSystemConstrain((femBandSystem *)mySolver->solver,myNode,myValue); break;
        case FEM_ITER : femIterativeSolverConstrain((femIterativeSolver *)mySolver->solver,myNode,myValue); break;
        default : Error("Unexpected solver type"); }
}

double *femSolverEliminate(femSolver *mySolver)
{
    double *soluce;
    switch (mySolver->type) {
        case FEM_FULL : soluce = femFullSystemEliminate((femFullSystem *)mySolver->solver); break;
        case FEM_BAND : soluce = femBandSystemEliminate((femBandSystem *)mySolver->solver); break;
        case FEM_ITER : soluce = femIterativeSolverEliminate((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
    return(soluce);
}



int femSolverConverged(femSolver *mySolver)
{
    int  testConvergence;
    switch (mySolver->type) {
        case FEM_FULL : testConvergence = 1; break;
        case FEM_BAND : testConvergence = 1; break;
        case FEM_ITER : testConvergence = femIterativeSolverConverged((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
    return(testConvergence);
}


femFullSystem *femFullSystemCreate(int size)
{
    femFullSystem *theSystem = malloc(sizeof(femFullSystem));
    theSystem->A = malloc(sizeof(double*) * size);
    theSystem->B = malloc(sizeof(double) * size * (size+1));
    theSystem->A[0] = theSystem->B + size;
    theSystem->size = size;
    int i;
    for (i=1 ; i < size ; i++)
        theSystem->A[i] = theSystem->A[i-1] + size;
    femFullSystemInit(theSystem);

    return theSystem;
}

void femFullSystemFree(femFullSystem *theSystem)
{
    free(theSystem->A);
    free(theSystem->B);
    free(theSystem);
}


void femFullSystemInit(femFullSystem *mySystem)
{
    int i,size = mySystem->size;
    for (i=0 ; i < size*(size+1) ; i++)
        mySystem->B[i] = 0;}


void femFullSystemPrint(femFullSystem *mySystem)
{
    double  **A, *B;
    int     i, j, size;
    A    = mySystem->A;
    B    = mySystem->B;
    size = mySystem->size;

    for (i=0; i < size; i++) {
        for (j=0; j < size; j++)
            if (A[i][j] == 0)  printf("         ");
            else               printf(" %+.1e",A[i][j]);
        printf(" :  %+.1e \n",B[i]); }
}

void femFullSystemPrintInfos(femFullSystem *mySystem)
{
    int  size = mySystem->size;
    printf(" \n");
    printf("    Full Gaussian elimination \n");
    printf("    Storage informations \n");
    printf("    Matrix size      : %8d\n",size);
    printf("    Bytes required   : %8d\n",(int)sizeof(double)*size*(size+1));
}

void  femFullSystemConstrain(femFullSystem *mySystem, int myNode, double myValue)
{
    double  **A, *B;
    int     i, size;

    A    = mySystem->A;
    B    = mySystem->B;
    size = mySystem->size;

    for (i=0; i < size; i++) {
        B[i] -= myValue * A[i][myNode];
        A[i][myNode] = 0; }

    for (i=0; i < size; i++)
        A[myNode][i] = 0;

    A[myNode][myNode] = 1;
    B[myNode] = myValue;
}

void femFullSystemAssemble(femFullSystem* mySystem, double *Aloc, double *Bloc, int *map, int nLoc)
{
    int i,j;
    for (i = 0; i < nLoc; i++) {
        for(j = 0; j < nLoc; j++) {
            mySystem->A[map[i]][map[j]] += Aloc[i*nLoc+j]; }
    mySystem->B[map[i]] += Bloc[i]; }
}

double* femFullSystemEliminate(femFullSystem *mySystem)
{
    double  **A, *B, factor;
    int     i, j, k, size;

    A    = mySystem->A;
    B    = mySystem->B;
    size = mySystem->size;

    /* Gauss elimination */

    for (k=0; k < size; k++) {
        if ( fabs(A[k][k]) <= 1e-8 ) {
            printf("Pivot index %d  ",k);
            printf("Pivot value %e  ",A[k][k]);
            Error("Cannot eliminate with such a pivot"); }
        for (i = k+1 ; i <  size; i++) {
            factor = A[i][k] / A[k][k];
            for (j = k+1 ; j < size; j++)
                A[i][j] = A[i][j] - A[k][j] * factor;
            B[i] = B[i] - B[k] * factor; }}

    /* Back-substitution */

    for (i = size-1; i >= 0 ; i--) {
        factor = 0;
        for (j = i+1 ; j < size; j++)
            factor += A[i][j] * B[j];
        B[i] = ( B[i] - factor)/A[i][i]; }

    return(mySystem->B);
}

femBandSystem *femBandSystemCreate(int size, int band)
{
    femBandSystem *myBandSystem = malloc(sizeof(femBandSystem));
    myBandSystem->B = malloc(sizeof(double)*size*(band+1));
    myBandSystem->A = malloc(sizeof(double*)*size);
    myBandSystem->size = size;
    myBandSystem->band = band;
    myBandSystem->A[0] = myBandSystem->B + size;
    int i;
    for (i=1 ; i < size ; i++)
        myBandSystem->A[i] = myBandSystem->A[i-1] + band - 1;
    femBandSystemInit(myBandSystem);
    return(myBandSystem);
}

void femBandSystemFree(femBandSystem *myBandSystem)
{
    free(myBandSystem->B);
    free(myBandSystem->A);
    free(myBandSystem);
}

void femBandSystemInit(femBandSystem *myBandSystem)
{
    int i;
    int size = myBandSystem->size;
    int band = myBandSystem->band;
    for (i=0 ; i < size*(band+1) ; i++)
        myBandSystem->B[i] = 0;
}

void femBandSystemPrint(femBandSystem *myBand)
{
    double  **A, *B;
    int     i, j, band, size;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;

    for (i=0; i < size; i++) {
        for (j=i; j < i+band; j++)
            if (A[i][j] == 0) printf("         ");
            else              printf(" %+.1e",A[i][j]);
        printf(" :  %+.1e \n",B[i]); }
}

void femBandSystemPrintInfos(femBandSystem *myBand)
{
    int size = myBand->size;
    int band = myBand->band;
    printf(" \n");
    printf("    Banded Gaussian elimination \n");
    printf("    Storage informations \n");
    printf("    Matrix size      : %8d\n",size);
    printf("    Matrix band      : %8d\n",band);
    printf("    Bytes required   : %8d\n",(int)sizeof(double)*size*(band+1));
}


double femBandSystemGet(femBandSystem* myBandSystem, int myRow, int myCol)
{
    double value = 0;
    if (myCol >= myRow && myCol < myRow+myBandSystem->band)  value = myBandSystem->A[myRow][myCol];
    return(value);
}

double femFullSystemGet(femFullSystem* myFullSystem, int myRow, int myCol)
{
    return(myFullSystem->A[myRow][myCol]);
}

void femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc)
{
    int i,j;
    for (i = 0; i < nLoc; i++) {
        int myRow = map[i];
        for(j = 0; j < nLoc; j++) {
            int myCol = map[j];
            if (myCol >= myRow)  myBandSystem->A[myRow][myCol] += Aloc[i*nLoc+j]; }
        myBandSystem->B[myRow] += Bloc[i]; }
}



void femBandSystemConstrain(femBandSystem *myBand, int myNode, double myValue)
{
    double  **A, *B;
    int     i, size, band, ifirst, iend;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;

    ifirst = fmax(0,myNode - band + 1);
    iend   = myNode;
    for (i=ifirst; i < iend; i++) {
        B[i] -= myValue * A[i][myNode];
        A[i][myNode] = 0; }

    ifirst = myNode+1;
    iend = fmin(myNode + band,size);
    for (i=ifirst; i < iend; i++) {
        B[i] -= myValue * A[myNode][i];
        A[myNode][i] = 0; }

    A[myNode][myNode] = 1;
    B[myNode] = myValue;
}

double  *femBandSystemEliminate(femBandSystem *myBand)
{
    double  **A, *B, factor;
    int     i, j, k, jend, size, band;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;

    /* Incomplete Cholesky factorization */

    for (k=0; k < size; k++) {
        if ( fabs(A[k][k]) <= 1e-4 ) {
            Error("Cannot eleminate with such a pivot"); }
        jend = fmin(k + band,size);
        for (i = k+1 ; i <  jend; i++) {
            factor = A[k][i] / A[k][k];
            for (j = i ; j < jend; j++)
                A[i][j] = A[i][j] - A[k][j] * factor;
            B[i] = B[i] - B[k] * factor; }}

    /* Back-substitution */

    for (i = (size-1); i >= 0 ; i--) {
        factor = 0;
        jend = fmin(i + band,size);
        for (j = i+1 ; j < jend; j++)
            factor += A[i][j] * B[j];
        B[i] = ( B[i] - factor)/A[i][i]; }

    return(myBand->B);
}



femIterativeSolver *femIterativeSolverCreate(int size)
{
    femIterativeSolver *mySolver = malloc(sizeof(femIterativeSolver));
    mySolver->R = malloc(sizeof(double)*size*4);
    mySolver->R0 = mySolver->R + size;
    mySolver->D  = mySolver->R + size*2;
    mySolver->D0 = mySolver->R + size*3;
    mySolver->size = size;
    femIterativeSolverInit(mySolver);
    return(mySolver);
}

void femIterativeSolverFree(femIterativeSolver *mySolver)
{
    free(mySolver->R);
    free(mySolver);
}

void femIterativeSolverInit(femIterativeSolver *mySolver)
{
    int i;
    mySolver->iter = 0;
    mySolver->error = 10.0e+12;
    for (i=0 ; i < mySolver->size*4 ; i++)
        mySolver->R[i] = 0;
}

void femIterativeSolverPrint(femIterativeSolver *mySolver)
{
    double  *R;
    int     i, size;
    R    = mySolver->R;
    size = mySolver->size;

    for (i=0; i < size; i++) {
        printf("%d :  %+.1e \n",i,R[i]); }
}

void femIterativeSolverPrintInfos(femIterativeSolver *mySolver)
{
    if (mySolver->iter == 1)     printf("\n    Iterative solver \n");
    printf("    Iteration %4d : %14.7e\n",mySolver->iter,mySolver->error);
}

void femIterativeSolverAssemble(femIterativeSolver* mySolver, double *Aloc, double *Bloc, double *Uloc, int *map, int nLoc)
{
    int i,j,myRow;
    if (mySolver->iter==0) {
        for (i = 0; i < nLoc; i++) {
            myRow = map[i];
            mySolver->R[myRow] += Bloc[i];
            mySolver->D[myRow] += Bloc[i];
            for(j = 0; j < nLoc; j++) {
                mySolver->R[myRow] -= Aloc[i*nLoc+j]*Uloc[j];
                mySolver->D[myRow] -= Aloc[i*nLoc+j]*Uloc[j]; }}}
    else {
    for (i = 0; i < nLoc; i++) {
        myRow = map[i];
        for(j = 0; j < nLoc; j++) {
            int myCol = map[j];
            mySolver->D0[myRow] += Aloc[i*nLoc+j] * mySolver->D[myCol]; }}}
}

void femIterativeSolverConstrain(femIterativeSolver* mySolver, int myNode, double myValue)
{
    mySolver->R[myNode] = myValue;
    mySolver->D0[myNode] = myValue;
    mySolver->D[myNode] = myValue;
}

// R
// R0  = dX
// D   = P
// D0  = AP

double *femIterativeSolverEliminate(femIterativeSolver *mySolver)
{
    mySolver->iter++;
    double error = 0.0; int i;
    double denAlpha = 0.0;
    for (i=0; i < mySolver->size; i++) {
        error += (mySolver->R[i])*(mySolver->R[i]);
        denAlpha += mySolver->D[i] * mySolver->D0[i]; }
    double alpha = error/denAlpha;

    if (mySolver->iter == 1) {
        for (i=0; i < mySolver->size; i++) {
            mySolver->R0[i] = 0.0; }}
    else {
        double numBeta = 0.0;
        for (i=0; i < mySolver->size; i++) {
            mySolver->R0[i] = alpha * mySolver->D[i];
            mySolver->R[i] = mySolver->R[i] - alpha * mySolver->D0[i];
            numBeta += mySolver->R[i] * mySolver->R[i]; }
        double beta = numBeta/error;
        for (i=0; i < mySolver->size; i++) {
            mySolver->D0[i] = 0.0;
            mySolver->D[i] = mySolver->R[i] + beta * mySolver->D[i]; }}

    mySolver->error = sqrt(error);
    return(mySolver->R0);
}



int femIterativeSolverConverged(femIterativeSolver *mySolver)
{
    int  testConvergence = 0;
    if (mySolver->iter  > 3000)     testConvergence = -1;
    if (mySolver->error < 10.0e-6)  testConvergence = 1;
    return(testConvergence);
}

femDiffusionProblem *femDiffusionCreate(const char *filename, femSolverType solverType, femRenumType renumType, int orderType)
{
    int i,band;
    int nodesOrderTriangle[] = {1,3,6,10};
    int integOrderTriangle[] = {1,3,12,12};
    int nodesOrderQuad[]     = {1,4,9,16};
    int integOrderQuad[]     = {1,4,9,9};

    femDiffusionProblem *myProblem = malloc(sizeof(femDiffusionProblem));
    myProblem->mesh  = femMeshRead(filename);
    myProblem->edges = femEdgesCreate(myProblem->mesh);
    int nNode = myProblem->mesh->nNode;
    int nElem = myProblem->mesh->nElem;
    int nEdge = myProblem->edges->nEdge;
    if (myProblem->mesh->nLocalNode == 4) {
        switch (orderType) {
            case 0 : myProblem->size = nElem; break;
            case 1 : myProblem->size = nNode; break;
            case 2 : myProblem->size = nNode + nEdge + nElem; break;
            case 3 : myProblem->size = nNode + 2*nEdge + 4*nElem; break;  }
        myProblem->space = femDiscreteCreate(nodesOrderQuad[orderType],FEM_QUAD);
        myProblem->rule = femIntegrationCreate(integOrderQuad[orderType],FEM_QUAD); }
    else if (myProblem->mesh->nLocalNode == 3) {
        switch (orderType) {
            case 0 : myProblem->size = nElem; break;
            case 1 : myProblem->size = nNode; break;
            case 2 : myProblem->size = nNode + nEdge; break;
            case 3 : myProblem->size = nNode + 2*nEdge + nElem; break;  }
        myProblem->space = femDiscreteCreate(nodesOrderTriangle[orderType],FEM_TRIANGLE);
        myProblem->rule = femIntegrationCreate(integOrderTriangle[orderType],FEM_TRIANGLE); }

    myProblem->number  = malloc(sizeof(int)*myProblem->size);
    myProblem->map  = malloc((myProblem->space->n)*sizeof(int)*myProblem->mesh->nElem);
    myProblem->soluce = malloc(sizeof(double)*myProblem->size);
    myProblem->X = malloc(sizeof(double)*myProblem->size);
    myProblem->Y = malloc(sizeof(double)*myProblem->size);
    for (i = 0; i < myProblem->size; i++)
        myProblem->soluce[i] = 0.0;
    femDiffusionComputeMap(myProblem,orderType);


    femDiffusionRenumber(myProblem,renumType);

    switch (solverType) {
        case FEM_FULL :
                myProblem->solver = femSolverFullCreate(myProblem->size); break;
        case FEM_BAND :
                band = femDiffusionComputeBand(myProblem);
                myProblem->solver = femSolverBandCreate(myProblem->size,band); break;
        case FEM_ITER :
               myProblem->solver = femSolverIterativeCreate(myProblem->size); break;
        default : Error("Unexpected solver option"); }


    return myProblem;
}

void femDiffusionFree(femDiffusionProblem *myProblem)
{
    femIntegrationFree(myProblem->rule);
    femDiscreteFree(myProblem->space);
    femEdgesFree(myProblem->edges);
    femMeshFree(myProblem->mesh);
    femSolverFree(myProblem->solver);
    free(myProblem->number);
    free(myProblem->soluce);
    free(myProblem->map);
    free(myProblem->X);
    free(myProblem->Y);
    free(myProblem);
}


void femDiffusionMeshLocal(const femDiffusionProblem *myProblem, const int iElem, int *map, double *x, double *y, double *u)
{
    int j,nLocal = myProblem->space->n;

    for (j=0; j < nLocal; ++j) {
        map[j] = myProblem->map[iElem*nLocal+j];
        x[j]   = myProblem->X[map[j]];
        y[j]   = myProblem->Y[map[j]];
        u[j]   = myProblem->soluce[map[j]];
        map[j] = myProblem->number[map[j]];
        }
}

void femDiffusionCompute(femDiffusionProblem *myProblem)
{
    femMesh *theMesh = myProblem->mesh;
    femIntegration *theRule = myProblem->rule;
    femDiscrete *theSpace = myProblem->space;
    femSolver *theSolver = myProblem->solver;
    femEdges *theEdges = myProblem->edges;
    int *number = myProblem->number;

    int n = theSpace->n;

    double Xloc[n],Yloc[n],phi[n],dphidxsi[n],dphideta[n],dphidx[n],dphidy[n],Aloc[n*n],Bloc[n],Uloc[n];
    int iEdge,iElem,iInteg,i,j,map[n];

    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (i = 0; i < theSpace->n; i++)      Bloc[i] = 0;
        for (i = 0; i < (theSpace->n)*(theSpace->n); i++) Aloc[i] = 0;
        femDiffusionMeshLocal(myProblem,iElem,map,Xloc,Yloc,Uloc);
        for (iInteg=0; iInteg < theRule->n; iInteg++) {
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            double dxdxsi = 0;
            double dxdeta = 0;
            double dydxsi = 0;
            double dydeta = 0;
            for (i = 0; i < theSpace->n; i++) {
                dxdxsi += Xloc[i]*dphidxsi[i];
                dxdeta += Xloc[i]*dphideta[i];
                dydxsi += Yloc[i]*dphidxsi[i];
                dydeta += Yloc[i]*dphideta[i]; }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            for (i = 0; i < theSpace->n; i++) {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; }
            for (i = 0; i < theSpace->n; i++) {
                for(j = 0; j < theSpace->n; j++) {
                    Aloc[i*(theSpace->n)+j] += (dphidx[i] * dphidx[j]
                                            + dphidy[i] * dphidy[j]) * jac * weight; }}
            for (i = 0; i < theSpace->n; i++) {
                Bloc[i] += phi[i] * jac *weight; }}
        femSolverAssemble(theSolver,Aloc,Bloc,Uloc,map,theSpace->n); }

    for (iEdge= 0; iEdge < theEdges->nEdge; iEdge++) {
        if (theEdges->edges[iEdge].elem[1] < 0) {
            femSolverConstrain(theSolver,number[theEdges->edges[iEdge].node[0]],0.0);
            femSolverConstrain(theSolver,number[theEdges->edges[iEdge].node[1]],0.0);
            for (j=1; j < theSpace->order; j++)
                femSolverConstrain(theSolver,number[iEdge*(theSpace->order-1)+theMesh->nNode+j-1],0.0);
            }}

    double *soluce = femSolverEliminate(theSolver);
    for (i = 0; i < myProblem->size; i++)
        myProblem->soluce[i] += soluce[number[i]];

}

static double *theGlobal;

int compare(const void *nodeOne, const void *nodeTwo)
{
    int *iOne = (int *)nodeOne;
    int *iTwo = (int *)nodeTwo;
    double diff = theGlobal[*iOne] - theGlobal[*iTwo];
    if (diff < 0)    return  1;
    if (diff > 0)    return -1;
                     return  0;
}

void femDiffusionRenumber(femDiffusionProblem *myProblem, femRenumType renumType)
{
    int i, *inverse;
    switch (renumType) {
        case FEM_NO :
            for (i = 0; i < myProblem->size; i++)
                myProblem->number[i] = i;
            break;
        case FEM_XNUM :
        case FEM_YNUM :
            inverse = malloc(sizeof(int)*myProblem->size);
            for (i = 0; i < myProblem->size; i++)
                inverse[i] = i;
            if (renumType == FEM_XNUM) theGlobal = myProblem->X;
            if (renumType == FEM_YNUM) theGlobal = myProblem->Y;
            qsort(inverse, myProblem->size, sizeof(int), compare);
            for (i = 0; i < myProblem->size; i++)
                myProblem->number[inverse[i]] = i;
            free(inverse);
            break;

        default : Error("Unexpected renumbering option"); }
}

int femDiffusionComputeBand(femDiffusionProblem *myProblem)
{
    femMesh *theMesh = myProblem->mesh;
    int iElem,j,myMax,myMin,myBand;
    int nLocal = myProblem->space->n;
    int map[nLocal];
    myBand = 0;
    for(iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; ++j)
            map[j] = myProblem->number[myProblem->map[iElem*nLocal+j]];
        myMin = map[0];
        myMax = map[0];
        for (j=1; j < nLocal; j++) {
            myMax = fmax(map[j],myMax);
            myMin = fmin(map[j],myMin); }
        if (myBand < (myMax - myMin)) myBand = myMax - myMin; }
    return(++myBand);
}




void femDiffusionComputeMap(femDiffusionProblem *myProblem, int orderType)
{
    femMesh  *theMesh = myProblem->mesh;
    femEdges *theEdges = myProblem->edges;

    int i,j,k, iElem;
    int nElem = theMesh->nElem;
    int nNode = theMesh->nNode;
    int nEdge = theEdges->nEdge;
    int nLocal = myProblem->space->n;
    int nLocalNode = myProblem->mesh->nLocalNode;
    double xsi[16],eta[16],phi[4];
    femDiscreteXsi2(myProblem->space,xsi,eta);
    femDiscrete *theLinearSpace;
    if (nLocalNode == 3) theLinearSpace = femDiscreteCreate(3,FEM_TRIANGLE);
    else        theLinearSpace = femDiscreteCreate(4,FEM_QUAD);

    if (orderType >= 1) {
        for (i=0; i < theMesh->nNode; i++) {
            myProblem->X[i] = theMesh->X[i];
            myProblem->Y[i] = theMesh->Y[i]; }
        for (i=0; i < nElem; i++)
            for (j=0; j < nLocalNode; j++)
                myProblem->map[i*nLocal+j] = myProblem->mesh->elem[i*nLocalNode+j]; }

    if (orderType == 2) {
        for (i=0; i < theEdges->nEdge; i++) {
            int node0 = theEdges->edges[i].node[0];
            int node1 = theEdges->edges[i].node[1];
            myProblem->X[i+nNode] = (theMesh->X[theEdges->edges[i].node[0]]
                                                        + theMesh->X[theEdges->edges[i].node[1]])/2.0;
            myProblem->Y[i+nNode] = (theMesh->Y[theEdges->edges[i].node[0]]
                                                        + theMesh->Y[theEdges->edges[i].node[1]])/2.0;
            for (j=0; j < 2; j++) {
                iElem = theEdges->edges[i].elem[j];
                if (iElem >= 0) {
                    int *elem = &myProblem->mesh->elem[iElem*nLocalNode];
                    for (k=0; k < nLocalNode; k++) {
                        int k1 = k+1;
                        if (k1 == nLocalNode) k1 = 0;
                        if (node0 == elem[k] && node1 == elem[k1]) myProblem->map[iElem*nLocal+nLocalNode+k] = i+nNode;
                        if (node1 == elem[k] && node0 == elem[k1]) myProblem->map[iElem*nLocal+nLocalNode+k] = i+nNode; }}}}}

    if (orderType == 2 && nLocalNode ==  4) {
            for (i=0; i < nElem; i++) {
                myProblem->map[i*nLocal+8] = i+nNode+nEdge;
                for (j=0; j < nLocalNode; j++) {
                    myProblem->X[i+nNode+nEdge] = 0.0;
                    myProblem->Y[i+nNode+nEdge] = 0.0; }
                for (j=0; j < nLocalNode; j++) {
                    myProblem->X[i+nNode+nEdge] += theMesh->X[theMesh->elem[i*nLocalNode+j]]/4.0;
                    myProblem->Y[i+nNode+nEdge] += theMesh->Y[theMesh->elem[i*nLocalNode+j]]/4.0; }}}


    if (orderType == 3) {
        for (i=0; i < theEdges->nEdge; i++) {
            int node0 = theEdges->edges[i].node[0];
            int node1 = theEdges->edges[i].node[1];
            myProblem->X[i*2+nNode] = (2.0 * theMesh->X[theEdges->edges[i].node[0]]
                                                        + theMesh->X[theEdges->edges[i].node[1]])/3.0;
            myProblem->Y[i*2+nNode] = (2.0 * theMesh->Y[theEdges->edges[i].node[0]]
                                                        + theMesh->Y[theEdges->edges[i].node[1]])/3.0;
            myProblem->X[i*2+nNode+1] = (theMesh->X[theEdges->edges[i].node[0]]
                                                        + 2.0 * theMesh->X[theEdges->edges[i].node[1]])/3.0;
            myProblem->Y[i*2+nNode+1] = (theMesh->Y[theEdges->edges[i].node[0]]
                                                        + 2.0 * theMesh->Y[theEdges->edges[i].node[1]])/3.0;

            for (j=0; j < 2; j++) {
                iElem = theEdges->edges[i].elem[j];
                if (iElem >= 0) {
                    int *elem = &myProblem->mesh->elem[iElem*nLocalNode];
                    for (k=0; k < nLocalNode; k++) {
                        int k1 = k+1;
                        if (k1 == nLocalNode) k1 = 0;
                        if (node0 == elem[k] && node1 == elem[k1]) {
                                myProblem->map[iElem*nLocal+nLocalNode+k*2] = i*2+nNode;
                                myProblem->map[iElem*nLocal+nLocalNode+k*2+1] = i*2+nNode + 1; }
                        if (node1 == elem[k] && node0 == elem[k1]) {
                                myProblem->map[iElem*nLocal+nLocalNode+k*2] = i*2+nNode + 1;
                                myProblem->map[iElem*nLocal+nLocalNode+k*2+1] = i*2+nNode; }}}}}}

     if (orderType == 3 && nLocalNode ==  3) {
            for (i=0; i < nElem; i++) {
                myProblem->map[i*nLocal+9] = i+nNode+2*nEdge;
                for (j=0; j < nLocalNode; j++) {
                    myProblem->X[i+nNode+nEdge*2] = 0.0;
                    myProblem->Y[i+nNode+nEdge*2] = 0.0; }
                for (j=0; j < nLocalNode; j++) {
                    myProblem->X[i+nNode+nEdge*2] += theMesh->X[theMesh->elem[i*nLocalNode+j]]/3.0;
                    myProblem->Y[i+nNode+nEdge*2] += theMesh->Y[theMesh->elem[i*nLocalNode+j]]/3.0; }}}

    if (orderType == 3 && nLocalNode ==  4) {
            for (i=0; i < nElem; i++) {
                for (k=0; k < 4; k++) {
                    femDiscretePhi2(theLinearSpace, xsi[12+k], eta[12+k], phi);
                    myProblem->map[i*nLocal+12+k] = i*4+nNode+2*nEdge+k;
                        for (j=0; j < nLocalNode; j++) {
                            myProblem->X[i*4+nNode+nEdge*2+k] = 0.0;
                            myProblem->Y[i*4+nNode+nEdge*2+k] = 0.0; }
                        for (j=0; j < nLocalNode; j++) {
                            myProblem->X[i*4+nNode+nEdge*2+k] += phi[j] * theMesh->X[theMesh->elem[i*nLocalNode+j]];
                            myProblem->Y[i*4+nNode+nEdge*2+k] += phi[j] * theMesh->Y[theMesh->elem[i*nLocalNode+j]]; }


                    }}}


                femDiscreteFree(theLinearSpace);

    femDiffusionPrintMap(myProblem);
}


void femDiffusionPrintMap(femDiffusionProblem *myProblem)
{
    femMesh  *theMesh = myProblem->mesh;
    femEdges *theEdges = myProblem->edges;

    int i,j;
    int nElem = theMesh->nElem;
    int nNode = theMesh->nNode;
    int nEdge = theEdges->nEdge;
    int nLocal = myProblem->space->n;
    int order = myProblem->space->order;

    printf( "Number of Elements %d \n", nElem);
    for (i = 0; i < nElem; ++i) {
        int *elem = &(myProblem->map[i*nLocal]);
        printf("%6d : ", i);
        for (j=0; j < nLocal; j++)
            printf("%6d ", elem[j]);
        printf("\n");  }

    printf( "Number of segments %d \n", nEdge);
    for (i = 0; i < nEdge;i++) {
        printf("%6d : ", i);
        printf("%6d ", theEdges->edges[i].node[0]);
        for (j=1; j < order; j++)
            printf("%6d ", i*(order-1)+nNode+j-1);
        printf("%6d ", theEdges->edges[i].node[1]);
        printf("\n");  }
    fflush(stdout);
}


void femDiffusionPrintInfos(femDiffusionProblem *myProblem)
{
    int  size = myProblem->size;
    printf(" \n");
    printf("    Diffusion problem \n");
    printf("    Number of nodal unknowns      : %8d\n",size);
}


femAdvProblem *femAdvCreate(const char *meshFileName)
{
    femAdvProblem *myProblem = malloc(sizeof(femAdvProblem));

    myProblem->mesh = femMeshRead(meshFileName);
    myProblem->edges = femEdgesCreate(myProblem->mesh);
    myProblem->space = femDiscreteCreate(3,FEM_TRIANGLE);
    myProblem->rule1d = femIntegrationCreate(2,FEM_EDGE);
    myProblem->rule2d = femIntegrationCreate(3,FEM_TRIANGLE);

    int sizeLoc = myProblem->space->n;
    int sizeGlo = myProblem->mesh->nElem * sizeLoc + 1;

    myProblem->size = sizeGlo;
    myProblem->C = malloc(sizeof(double)*sizeGlo);
    myProblem->U = malloc(sizeof(double)*sizeGlo);
    myProblem->V = malloc(sizeof(double)*sizeGlo);
    myProblem->F = malloc(sizeof(double)*sizeGlo);

    return myProblem;
}



void femAdvInitialCondition(femAdvProblem *myProblem)
{
    double  u,v,elevation;
    int     j,elem,*node;
    double *X = myProblem->mesh->X;
    double *Y = myProblem->mesh->Y;
    double *C = myProblem->C;
    double *U = myProblem->U;
    double *V = myProblem->V;

    for (elem=0; elem < myProblem->mesh->nElem; elem++) {
        node = &(myProblem->mesh->elem[elem*3]);
        for (j=0; j < 3; ++j) {
            double x = X[node[j]];
            double y = Y[node[j]];
            double c = exp(- pow((x - 0.8), 2) / 0.01) * exp(- pow((y - 0.5), 2) /0.01);
            femStommel(x, y, &u, &v, &elevation);
            C[elem*3+j] = c;
            U[elem*3+j] = u;
            V[elem*3+j] = v; }}
    int last = myProblem->size - 1;
    C[last] = 0.0; U[last] = 0.0; V[last] = 0.0;
}


void femAdvAddIntegralsTriangles(femAdvProblem *myProblem)
{
    double *B = myProblem->F;
    double *C = myProblem->C;
    double *U = myProblem->U;
    double *V = myProblem->V;

    femIntegration *theRule = myProblem->rule2d;
    femDiscrete *theSpace = myProblem->space;

    double  xLoc[3],yLoc[3],dphidx[3],dphidy[3], phi[3];
    double  xsi,eta,weight,jac;
    double  c,u,v;
    int     i,j,k,elem,mapElem[3];

    for (elem=0; elem < myProblem->mesh->nElem; elem++) {
        femAdvTriangleMap(myProblem,elem,mapElem);
        int *mapCoord = &(myProblem->mesh->elem[elem*3]);
        for (j=0; j < 3; ++j) {
            xLoc[j] = myProblem->mesh->X[mapCoord[j]];
            yLoc[j] = myProblem->mesh->Y[mapCoord[j]]; }
        jac = (xLoc[1] - xLoc[0]) * (yLoc[2] - yLoc[0]) - (yLoc[1] - yLoc[0]) * (xLoc[2] - xLoc[0]);
        dphidx[0] = (yLoc[1] - yLoc[2])/jac;
        dphidx[1] = (yLoc[2] - yLoc[0])/jac;
        dphidx[2] = (yLoc[0] - yLoc[1])/jac;
        dphidy[0] = (xLoc[2] - xLoc[1])/jac;
        dphidy[1] = (xLoc[0] - xLoc[2])/jac;
        dphidy[2] = (xLoc[1] - xLoc[0])/jac;
        for (k=0; k < theRule->n; k++) {
            xsi = theRule->xsi[k];
            eta = theRule->eta[k];
            weight = theRule->weight[k];
            femDiscretePhi2(theSpace,xsi,eta,phi);
            u = femDiscreteInterpolate(phi,U,mapElem,3);
            v = femDiscreteInterpolate(phi,V,mapElem,3);
            c = femDiscreteInterpolate(phi,C,mapElem,3);

            for (i=0; i < 3; i++) {
                B[mapElem[i]] += (u*c*dphidx[i] + v*c*dphidy[i])*jac*weight; }}}
}

void femAdvAddIntegralsEdges(femAdvProblem *myProblem)
{
    double *B = myProblem->F;
    double *C = myProblem->C;
    double *U = myProblem->U;
    double *V = myProblem->V;
    femIntegration *theRule = myProblem->rule1d;
    femDiscrete *theSpace = myProblem->space;

    double  xEdge[2],yEdge[2],phiEdge[2];
    double  xsi,weight,jac;
    double  c,u,v;
    int     i,j,k,edge,mapEdge[2][2];

    for (edge=0; edge < myProblem->edges->nEdge; edge++) {
        femAdvEdgeMap(myProblem,edge,mapEdge);
        for (j=0; j < 2; ++j) {
            int node = myProblem->edges->edges[edge].node[j];
            xEdge[j] = myProblem->mesh->X[node];
            yEdge[j] = myProblem->mesh->Y[node]; }
        double dxdxsi = xEdge[1] - xEdge[0];
        double dydxsi = yEdge[1] - yEdge[0];
        double norm = sqrt(dxdxsi*dxdxsi + dydxsi*dydxsi);
        double normal[2] = {dydxsi/norm, -dxdxsi/norm};
        jac = norm / 2.0;
        for (k=0; k < theRule->n; k++) {
            xsi = theRule->xsi[k];
            weight = theRule->weight[k];
            femDiscretePhi1(theSpace,xsi,phiEdge);
            u = (femDiscreteInterpolate(phiEdge,U,mapEdge[0],2) + femDiscreteInterpolate(phiEdge,U,mapEdge[1],2))/2.0;
            v = (femDiscreteInterpolate(phiEdge,V,mapEdge[0],2) + femDiscreteInterpolate(phiEdge,V,mapEdge[1],2))/2.0;
            double beta = u*normal[0]+ v*normal[1];
            if (beta > 0) c = femDiscreteInterpolate(phiEdge,C,mapEdge[0],2);
            else          c = femDiscreteInterpolate(phiEdge,C,mapEdge[1],2);
            for (i=0; i < 2; i++) {
                B[mapEdge[0][i]] -= beta*c*phiEdge[i] *jac*weight;
                B[mapEdge[1][i]] += beta*c*phiEdge[i] *jac*weight; }}}
}

void femAdvMultiplyInverseMatrix(femAdvProblem *myProblem)
{
    double *B = myProblem->F;
    double invA[3][3] = {{18.0,-6.0,-6.0},{-6.0,18.0,-6.0},{-6.0,-6.0,18.0}};
    double Bloc[3];

    double  xLoc[3],yLoc[3],jac;
    int     i,j,elem,mapElem[3];

    for (elem=0; elem < myProblem->mesh->nElem; elem++) {
        femAdvTriangleMap(myProblem,elem,mapElem);
        int *mapCoord = &(myProblem->mesh->elem[elem*3]);
        for (j=0; j < 3; ++j) {
            xLoc[j] = myProblem->mesh->X[mapCoord[j]];
            yLoc[j] = myProblem->mesh->Y[mapCoord[j]]; }
        jac = (xLoc[1] - xLoc[0]) * (yLoc[2] - yLoc[0]) - (yLoc[1] - yLoc[0]) * (xLoc[2] - xLoc[0]);
        for (i=0; i < 3; i++) {
            Bloc[i] = B[mapElem[i]];
            B[mapElem[i]] = 0.0; }
        for (i=0; i < 3; i++) {
            for (j=0; j < 3; j++) {
                B[mapElem[i]] += invA[i][j] * Bloc[j] / jac; }}}
}

void femAdvTriangleMap(femAdvProblem *myProblem, int index, int map[3])
{
    int j;
    for (j=0; j < 3; ++j)
        map[j] = index*3 + j;
}

void femAdvEdgeMap(femAdvProblem *myProblem, int index, int map[2][2])
{
    int i,j,k;

    for (j=0; j < 2; ++j) {
        int node = myProblem->edges->edges[index].node[j];
        for (k=0; k < 2; k++) {
            int elem = myProblem->edges->edges[index].elem[k];
            map[k][j] = (myProblem->mesh->nElem)*3;
            if (elem >= 0) {
                for (i=0; i < 3; i++) {
                    if (myProblem->mesh->elem[elem*3 + i] == node) {
                        map[k][j] = elem*3 + i;  }}}}}
}

void femAdvFree(femAdvProblem *myProblem)
{
    free(myProblem->F);
    free(myProblem->U);
    free(myProblem->V);
    free(myProblem->C);
    femIntegrationFree(myProblem->rule1d);
    femIntegrationFree(myProblem->rule2d);
    femDiscreteFree(myProblem->space);
    femEdgesFree(myProblem->edges);
    femMeshFree(myProblem->mesh);
    free(myProblem);
}

femShallowProblem *femShallowCreate(const char *meshFileName)
{
    femShallowProblem *myProblem = malloc(sizeof(femShallowProblem));
    myProblem->omega   = (2*M_PI)/86400;
    myProblem->gravity  = 9.81;
    myProblem->gamma    = 1e-7;
    myProblem->R     = 6371220.0;


    myProblem->mesh = femMeshRead(meshFileName);
    myProblem->edges = femEdgesCreate(myProblem->mesh);
    myProblem->rule1d = femIntegrationCreate(2,FEM_EDGE);
    if (myProblem->mesh->nLocalNode == 4) {
        myProblem->space = femDiscreteCreate(4,FEM_QUAD);
        myProblem->rule2d = femIntegrationCreate(9,FEM_QUAD);  }
    else if (myProblem->mesh->nLocalNode == 3) {
        myProblem->space = femDiscreteCreate(3,FEM_TRIANGLE);
        myProblem->rule2d = femIntegrationCreate(3,FEM_TRIANGLE); }
    int sizeLoc = myProblem->space->n;
    int sizeGlo = myProblem->mesh->nElem * sizeLoc + 1;
    //myProblem->solver = femSolverFullCreate(3*sizeLoc);
    //permet de gagner du temps
    myProblem->solver = femSolverBandCreate(3*sizeLoc,sizeLoc);
    myProblem->size = sizeGlo;
    myProblem->E = malloc(sizeof(double)*sizeGlo);
    myProblem->U = malloc(sizeof(double)*sizeGlo);
    myProblem->V = malloc(sizeof(double)*sizeGlo);
    myProblem->FE = malloc(sizeof(double)*sizeGlo);
    myProblem->FU = malloc(sizeof(double)*sizeGlo);
    myProblem->FV = malloc(sizeof(double)*sizeGlo);
    int i,j;
    for (i=0; i < myProblem->mesh->nElem; i++)
        {
            for(j=0;j<sizeLoc;j++)
            {
                myProblem->E[i] = 0.0;
                myProblem->U[i] = 0.0;
                myProblem->V[i] = 0.0;
            }
        }

    return myProblem;
}
//cree un nouveau femshallow



void femShallowWrite(femShallowProblem *myProblem,char *filename) {
    int i;

    FILE* file = fopen(filename,"w");
    if (file == NULL) Error("Cannot create the result file !");

    fprintf(file, "Number of values %d \n", myProblem->size);
    for (i = 0; i < myProblem->size; ++i) {
        fprintf(file,"%le;%le;%le;\n",myProblem->E[i],myProblem->U[i],myProblem->V[i]);
         }
    fclose(file);

}

void femShallowRead(femShallowProblem *myProblem,char *filename) {
    int i;

    FILE* file = fopen(filename,"r");
    if (file == NULL) Error("No result file !");
    int dummy;

    fscanf(file, "Number of values %d \n", &dummy);
    for (i = 0; i < myProblem->size; ++i) {
        fscanf(file,"%le;",&myProblem->E[i]);
        fscanf(file,"%le;",&myProblem->U[i]);
        fscanf(file,"%le;",&myProblem->V[i]);
        fscanf(file,"\n"); }
    fclose(file);

}

void femShallowFree(femShallowProblem *myProblem)
{

    free(myProblem->FE);
    free(myProblem->FU);
    free(myProblem->FV);
    free(myProblem->E);
    free(myProblem->U);
    free(myProblem->V);
    femSolverFree(myProblem->solver);
    femIntegrationFree(myProblem->rule1d);
    femIntegrationFree(myProblem->rule2d);
    femDiscreteFree(myProblem->space);
    femEdgesFree(myProblem->edges);
    femMeshFree(myProblem->mesh);
    free(myProblem);
}

void femStommel(double x, double y, double *u, double *v, double *eta)
{
    //
    // Solution analytique de Stommel dans un carre [0,1]x[0,1]
    // Modelisation de l'elevation de l'ocean Atlantique dans un carre adimensionnel
    // Ce modele que l'on attribue generalement au grand oceanographe Henry M.
    // Stommel (1920-1992), est considere comme le premier modele qualitativement correct du Gulf Stream
    //

    const double tau0 = 0.1;
    const double L = 1e6;
    const double gamm = 1e-6;
    const double rho = 1000;
    const double delta = 1;
    const double g = 9.81;
    const double h = 1000;
    const double f0 = 1e-4;
    const double beta = 0.5e-11;

    double Y = y - 0.5;
    double epsilon = gamm / (L * beta);
    double Z1 = (-1 + sqrt(1 + (2 * M_PI * delta * epsilon) * (2 * M_PI * delta * epsilon))) / (2 * epsilon);
    double Z2 = (-1 - sqrt(1 + (2 * M_PI * delta * epsilon) * (2 * M_PI * delta * epsilon))) / (2 * epsilon);
    double D  = ((exp(Z2) - 1) * Z1 + (1 - exp(Z1)) * Z2) / (exp(Z1) - exp(Z2));
    double f1 = M_PI / D * (1 + ((exp(Z2) - 1) * exp(x * Z1) + (1 - exp(Z1)) * exp(x * Z2)) / (exp(Z1) - exp(Z2)));
    double f2 = 1 / D* (((exp(Z2) - 1) * Z1 * exp(x * Z1) + (1 - exp(Z1)) * Z2 * exp(x * Z2)) / (exp(Z1) - exp(Z2)));

    eta[0] = D * tau0 * f0 * L / (M_PI * gamm * rho * delta * g * h) *
               ( - gamm / (f0 * delta * M_PI) * f2 * sin(M_PI * Y)
                 + 1 / M_PI * f1 * (cos(M_PI * Y) * (1 + beta * Y)
                 - beta / M_PI * sin(M_PI * Y) ) );
    u[0] = D * tau0 / (M_PI * gamm * rho * h) * f1 * sin(M_PI * Y);
    v[0] = D * tau0 / (M_PI * gamm * rho * delta * h) * f2 * cos(M_PI * Y);
}

double femMin(double *x, int n)
{
    double myMin = x[0];
    int i;
    for (i=1 ;i < n; i++)
        myMin = fmin(myMin,x[i]);
    return myMin;
}

double femMax(double *x, int n)
{
    double myMax = x[0];
    int i;
    for (i=1 ;i < n; i++)
        myMax = fmax(myMax,x[i]);
    return myMax;
}

void femError(char *text, int line, char *file)
{
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
    exit(69);
}

void femWarning(char *text, int line, char *file)
{
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Warning in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
}

void Cond(femShallowProblem *myProblem)
{

	double *Y = myProblem->mesh->Y;
    double *X = myProblem->mesh->X;
	double *E = myProblem->E;
	int nElem = myProblem->mesh->nElem;
	int nLocalNode = myProblem->mesh->nLocalNode;
	int x=0.0,y=0.0;
	int i,j;

	for (i=0; i<nElem; i++)
	{
            for(j=0;j<nLocalNode;j++)
            {
                x	=	X[myProblem->mesh->elem[i*nLocalNode+j]];
                y	=	Y[myProblem->mesh->elem[i*nLocalNode+j]];
                E[i*nLocalNode+j] = tsunamiInitialConditionOkada(x,y);
            }
	}

}

void tsunamiCompute(double dt, int nmax, int sub, const char *meshFileName, const char *baseResultName)
{
    femTsunamiProblem *myProblem = femTsunamiCreate(meshFileName);
    myProblem->timeStep = dt;
    double discreteTime = 0.0;
    double stop = nmax * dt;
    int iteration = 0;
    double time = 0.0;
    double timeUntilNow = 0.0;
    int percent = 0;
    myProblem->dt = dt;

    int nElem = myProblem->mesh->nElem;
    int n = myProblem->space->n;

    while (stop>=discreteTime)
    {
        if (iteration % sub == 0)
        {
            double *U = myProblem->U;
            double *V = myProblem->V;
            double *H = myProblem->H;
            time = clock();

            percent = (iteration*100)/nmax;

            if (iteration==0)
            {
                printf("Iteration : %d / %d > Elapsed time : %g (s) ; %d%%\n",iteration,nmax,time/1000.0, percent);
            }
            else
            {
                timeUntilNow = (time * (((double) nmax)/((double) iteration)) - time)/1000.0;
                printf("Iteration : %d / %d > Elapsed time : %g (s) ; Time remaining : %g (m) ; Done : %d%%\n",iteration,nmax,time/1000.0,timeUntilNow/60.0, percent);
            }
            tsunamiWriteFile(baseResultName, iteration, U, V, H, nElem, n);
        }
        femTsunamiCompute(myProblem);
        iteration += 1;
        discreteTime += dt;
    }

    femTsunamiFree(myProblem);

}

/*
*
*!!! Le code de la partie graphique du devoir est très très fortement inspiré du travail de Benjamin Spitaels et Alexis Godfrin !!!
*
*/
static int gRasterH = 800;
static int gRasterV = 600;
GLuint fontOffset;

GLubyte space[] =
{0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};

GLubyte letters[][13] = {
{0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xff, 0xc3, 0xc3, 0xc3, 0x66, 0x3c, 0x18}, // A
{0x00, 0x00, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe}, // B
{0x00, 0x00, 0x7e, 0xe7, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xe7, 0x7e}, // C
{0x00, 0x00, 0xfc, 0xce, 0xc7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc7, 0xce, 0xfc}, // D
{0x00, 0x00, 0xff, 0xc0, 0xc0, 0xc0, 0xc0, 0xfc, 0xc0, 0xc0, 0xc0, 0xc0, 0xff}, // E
{0x00, 0x00, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xfc, 0xc0, 0xc0, 0xc0, 0xff}, // F
{0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xcf, 0xc0, 0xc0, 0xc0, 0xc0, 0xe7, 0x7e}, // G
{0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xff, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3}, // H
{0x00, 0x00, 0x7e, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x7e}, // I
{0x00, 0x00, 0x7c, 0xee, 0xc6, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06}, // J
{0x00, 0x00, 0xc3, 0xc6, 0xcc, 0xd8, 0xf0, 0xe0, 0xf0, 0xd8, 0xcc, 0xc6, 0xc3}, // K
{0x00, 0x00, 0xff, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0}, // L
{0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xdb, 0xff, 0xff, 0xe7, 0xc3}, // M
{0x00, 0x00, 0xc7, 0xc7, 0xcf, 0xcf, 0xdf, 0xdb, 0xfb, 0xf3, 0xf3, 0xe3, 0xe3}, // N
{0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xe7, 0x7e}, // O
{0x00, 0x00, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe}, // P
{0x00, 0x00, 0x3f, 0x6e, 0xdf, 0xdb, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0x66, 0x3c}, // Q
{0x00, 0x00, 0xc3, 0xc6, 0xcc, 0xd8, 0xf0, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe}, // R
{0x00, 0x00, 0x7e, 0xe7, 0x03, 0x03, 0x07, 0x7e, 0xe0, 0xc0, 0xc0, 0xe7, 0x7e}, // S
{0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0xff}, // T
{0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3}, // U
{0x00, 0x00, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3}, // V
{0x00, 0x00, 0xc3, 0xe7, 0xff, 0xff, 0xdb, 0xdb, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3}, // W
{0x00, 0x00, 0xc3, 0x66, 0x66, 0x3c, 0x3c, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3}, // X
{0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3}, // Y
{0x00, 0x00, 0xff, 0xc0, 0xc0, 0x60, 0x30, 0x7e, 0x0c, 0x06, 0x03, 0x03, 0xff} // Z
};

GLubyte lowletters[][13] = {
{0x00, 0x00, 0x7d, 0xc3, 0xc3, 0xc3, 0x7f, 0x03, 0x7e, 0x00, 0x00, 0x00, 0x00}, // a
{0x00, 0x00, 0xfe, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xfe, 0xc0, 0xc0, 0xc0, 0xc0}, // b
{0x00, 0x00, 0x7f, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0x7f, 0x00, 0x00, 0x00, 0x00}, // c
{0x00, 0x00, 0x7f, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0x7f, 0x03, 0x03, 0x03, 0x03}, // d
{0x00, 0x00, 0x7e, 0xc0, 0xc0, 0xfe, 0xc3, 0xc3, 0x7e, 0x00, 0x00, 0x00, 0x00}, // e
{0x00, 0x00, 0x3c, 0x18, 0x18, 0x18, 0x18, 0x18, 0x7e, 0x18, 0x18, 0x18, 0x0e}, // f
{0x7f, 0x03, 0x7f, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0x7e, 0x00, 0x00, 0x00, 0x00}, // g
{0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xfe, 0xc0, 0xc0, 0xc0, 0xc0}, // h
{0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x00, 0x18, 0x18, 0x00}, // i
{0x70, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x00, 0x18, 0x18, 0x00}, // j
{0x00, 0x00, 0xc3, 0xc7, 0xce, 0xfc, 0xfe, 0xc7, 0xc3, 0xc0, 0xc0, 0xc0, 0xc0}, // k
{0x00, 0x00, 0x0c, 0x1c, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18}, // l
{0x00, 0x00, 0xdb, 0xdb, 0xdb, 0xdb, 0xdb, 0xdb, 0xfe, 0x00, 0x00, 0x00, 0x00}, // m
{0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xfe, 0x00, 0x00, 0x00, 0x00}, // n
{0x00, 0x00, 0x7e, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0x7e, 0x00, 0x00, 0x00, 0x00}, // o
{0xc0, 0xc0, 0xfe, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xfe, 0x00, 0x00, 0x00, 0x00}, // p
{0x03, 0x03, 0x7f, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0x7f, 0x00, 0x00, 0x00, 0x00}, // q
{0x00, 0x00, 0xc0, 0xc0, 0xc0, 0xc0, 0xe0, 0xf0, 0xdf, 0x00, 0x00, 0x00, 0x00}, // r
{0x00, 0x00, 0xfe, 0x03, 0x03, 0x7e, 0xc0, 0xc0, 0x7f, 0x00, 0x00, 0x00, 0x00}, // s
{0x00, 0x00, 0x0e, 0x18, 0x18, 0x18, 0x18, 0x18, 0x7e, 0x18, 0x18, 0x18, 0x18}, // t
{0x00, 0x00, 0x7f, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0x00, 0x00, 0x00, 0x00}, // u
{0x00, 0x00, 0x18, 0x3c, 0x66, 0x66, 0xc3, 0xc3, 0xc3, 0x00, 0x00, 0x00, 0x00}, // v
{0x00, 0x00, 0x66, 0x7e, 0xdb, 0xdb, 0xdb, 0xdb, 0xdb, 0x00, 0x00, 0x00, 0x00}, // w
{0x00, 0x00, 0xc3, 0xe7, 0x3c, 0x18, 0x3c, 0xe7, 0xc3, 0x00, 0x00, 0x00, 0x00}, // x
{0x7f, 0x03, 0x7f, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0x00, 0x00, 0x00, 0x00}, // y
{0x00, 0x00, 0xff, 0xc0, 0x70, 0x1c, 0x06, 0x03, 0xff, 0x00, 0x00, 0x00, 0x00}, // z
};

GLubyte numletters[][13] = {
{0x00, 0x00, 0x3c, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x3c}, // 0
{0x00, 0x00, 0x3c, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x78, 0x38, 0x18}, // 1
{0x00, 0x00, 0x7e, 0x60, 0x60, 0x60, 0x60, 0x3c, 0x06, 0x06, 0x66, 0x66, 0x3c}, // 2
{0x00, 0x00, 0x3c, 0x66, 0x06, 0x06, 0x06, 0x1c, 0x06, 0x06, 0x06, 0x66, 0x3c}, // 3
{0x00, 0x00, 0x06, 0x06, 0x06, 0x06, 0x06, 0x7f, 0x66, 0x36, 0x1e, 0x0e, 0x06}, // 4
{0x00, 0x00, 0x3c, 0x66, 0x06, 0x06, 0x06, 0x7c, 0x60, 0x60, 0x60, 0x60, 0x7e}, // 5
{0x00, 0x00, 0x3c, 0x66, 0x66, 0x66, 0x66, 0x66, 0x7c, 0x60, 0x60, 0x66, 0x3c}, // 6
{0x00, 0x00, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x1f, 0x06, 0x06, 0x06, 0x06, 0x7e}, // 7
{0x00, 0x00, 0x3c, 0x66, 0x66, 0x66, 0x66, 0x3c, 0x66, 0x66, 0x66, 0x66, 0x3c}, // 8
{0x00, 0x00, 0x3c, 0x66, 0x06, 0x06, 0x06, 0x3e, 0x66, 0x66, 0x66, 0x66, 0x3c}, // 9
{0x00, 0x00, 0x00, 0x18, 0x18, 0x00, 0x00, 0x18, 0x18, 0x00, 0x00, 0x00, 0x00}, // :
{0x00, 0x00, 0x30, 0x18, 0x18, 0x00, 0x00, 0x18, 0x18, 0x00, 0x00, 0x00, 0x00}, // ;
{0x00, 0x00, 0x06, 0x1c, 0x30, 0x60, 0x30, 0x1c, 0x06, 0x00, 0x00, 0x00, 0x00}, // <
{0x00, 0x00, 0x00, 0x00, 0x3c, 0x00, 0x00, 0x3c, 0x00, 0x00, 0x00, 0x00, 0x00}, // =
{0x00, 0x00, 0x60, 0x38, 0x0c, 0x06, 0x0c, 0x38, 0x60, 0x00, 0x00, 0x00, 0x00}, // >
{0x00, 0x00, 0x18, 0x18, 0x00, 0x18, 0x18, 0x18, 0x0c, 0x06, 0x06, 0x66, 0x3c}, // ?
};

GLubyte specialletters[][13] = {
{0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}, // space
{0x00, 0x00, 0x18, 0x18, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18}, // !
{0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x24, 0x24, 0x00, 0x00}, // "
{0x00, 0x00, 0x24, 0x24, 0x7e, 0x7e, 0x24, 0x7e, 0x7e, 0x24, 0x24, 0x00, 0x00}, // #
{0x00, 0x00, 0x18, 0x3c, 0x5a, 0x5a, 0x1a, 0x3c, 0x58, 0x58, 0x5a, 0x3c, 0x18}, // $
{0x00, 0x00, 0x44, 0x4a, 0x6a, 0x24, 0x30, 0x18, 0x0c, 0x24, 0x56, 0x52, 0x22}, // %
{0x00, 0x00, 0x79, 0xcf, 0xc6, 0xcf, 0x79, 0x70, 0x78, 0xcc, 0xcc, 0xcc, 0x78}, // &
{0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x08, 0x08, 0x18, 0x00, 0x00}, // '
{0x00, 0x00, 0x0c, 0x18, 0x18, 0x30, 0x30, 0x30, 0x30, 0x30, 0x18, 0x18, 0x0c}, // (
{0x00, 0x00, 0x30, 0x18, 0x18, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x18, 0x18, 0x30}, // )
{0x00, 0x00, 0x00, 0x00, 0x10, 0x54, 0x38, 0x54, 0x10, 0x00, 0x00, 0x00, 0x00}, // *
{0x00, 0x00, 0x00, 0x00, 0x10, 0x10, 0x7c, 0x10, 0x10, 0x00, 0x00, 0x00, 0x00}, // +
{0x00, 0x30, 0x18, 0x18, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}, // ,
{0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x3c, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}, // -
{0x00, 0x00, 0x18, 0x18, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}, // .
{0x00, 0x00, 0x60, 0x60, 0x30, 0x30, 0x18, 0x18, 0x18, 0x0c, 0x0c, 0x06, 0x06}, // /
};

void glMakeRasterFont(void)
{
    GLuint i, j;
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glShadeModel (GL_FLAT);
    fontOffset = glGenLists (128);

    for (i = 0,j = 'A'; i < 26; i++,j++) {
        glNewList(fontOffset + j, GL_COMPILE);
        glBitmap(8, 13, 0.0, 2.0, 10.0, 0.0, letters[i]);
        glEndList(); }

    for (i = 0,j = 'a'; i < 26; i++,j++) {
        glNewList(fontOffset + j, GL_COMPILE);
        glBitmap(8, 13, 0.0, 2.0, 10.0, 0.0, lowletters[i]);
        glEndList(); }

    for (i = 0,j = '0'; i < 16; i++,j++) {
        glNewList(fontOffset + j, GL_COMPILE);
        glBitmap(8, 13, 0.0, 2.0, 10.0, 0.0, numletters[i]);
        glEndList(); }

    for (i = 0,j = ' '; i < 16; i++,j++) {
        glNewList(fontOffset + j, GL_COMPILE);
        glBitmap(8, 13, 0.0, 2.0, 10.0, 0.0, specialletters[i]);
        glEndList(); }
    glShadeModel (GL_SMOOTH);
}

void glfemSetRasterSize(int h, int v)
{
    gRasterH = h;
    gRasterV = v;
}

void glfemDrawMessage(int h, int v, char *s)
{
    int off;
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glDisable(GL_TEXTURE_2D);
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho (0.0, gRasterH, gRasterV, 0.0, -1.0, 1.0);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    if (h < 0) h = gRasterH + h - strlen(s)*10;
    if (v < 0) v = gRasterV + v;
    glRasterPos2i(h, v);
    glListBase(fontOffset);

    if (h >= 0) {
        glRasterPos2i(h, v);
        glCallLists((GLsizei) strlen(s), GL_UNSIGNED_BYTE, (GLubyte *) s); }
    else {
        off = (h-9)/10;
        glRasterPos2i(h - off*10, v);
        if (strlen(s)+off > 0) glCallLists((GLsizei) strlen(s)+off, GL_UNSIGNED_BYTE, (GLubyte *) &s[-off]); }

    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopAttrib ();
}

void glfemInit(char *theWindowName)
{
    glfwInit();
    glfwOpenWindow(480,480,0,0,0,0,0,0,GLFW_WINDOW);
    glfemSetRasterSize(480,480);
    glfwSetWindowTitle(theWindowName);
    glShadeModel(GL_SMOOTH);
    glMakeRasterFont();
}


void tsunamiAnimate(double dt, int nmax, int sub, const char *meshFileName, const char *baseResultName)
{
    int nElem, nLocalNode, nLocal, nNodeMesh, nNodeTot, nEdge, size, i, j, index, *elem;
    double *X, *Y, *E, *U, *V;
    int width,height,mouse;
    int order =1;
    int nout = sub;
    double R = 6371220.0;

    femMesh *theMesh = femMeshRead(meshFileName);
    nLocalNode = theMesh->nLocalNode;
    X = theMesh->X;
    Y = theMesh->Y;
    elem = theMesh->elem;
    nElem = theMesh->nElem;
    nNodeMesh = theMesh->nNode;

    femEdges *edgesStruc = femEdgesCreate(theMesh);
    nEdge = edgesStruc->nEdge;

    if (nLocalNode == 3)
    {
        nLocal = 3;
        nNodeTot = nNodeMesh;
    }
    else if (nLocalNode == 4)
    {
        nLocal = 4;
        nNodeTot = nNodeMesh;
    }

    size = 3 * nLocalNode;
 	GLfloat colors[size], coord[size], coordBis[size], AngleAxeVert, AngleAxeHor, y;

 	U = malloc(sizeof(double)*nLocal*nElem);
 	V = malloc(sizeof(double)*nLocal*nElem);
 	E = malloc(sizeof(double)*nLocal*nElem);

    int *elemO2;
    double *XX, *YY;

    XX = malloc(nNodeTot * sizeof(double));
    YY = malloc(nNodeTot * sizeof(double));
    elemO2 = malloc(nLocal * nElem * sizeof(int));

    int numNode = nNodeMesh;
    int nmaxBoucle;
    if (nNodeMesh>=nElem) { nmaxBoucle = nNodeMesh; }
    else { nmaxBoucle = nElem; }

    for (i=0; i < nmaxBoucle; i++)
    {
        if (i<nNodeMesh)
        {
            XX[i] = X[i];
            YY[i] = Y[i];
        }
        if (i<nElem)
        {
            for (j=0; j < nLocalNode; j++)
            {
                elemO2[i*nLocal+j] = elem[i*nLocalNode+j];
            }
        }

    }

        if (nLocalNode == 4)
        {
            for( i=0 ; i<nElem ; i++ )
            {
                int node1 = elemO2[i*nLocal+nLocalNode];
                double coord1[2] = {XX[node1] , YY[node1]};
                int node3 = elemO2[i*nLocal+nLocalNode+2];
                double coord3[2] = {XX[node3] , YY[node3]};

                XX[numNode] = (coord1[0] + coord3[0])/2.0;
                YY[numNode] = (coord1[1] + coord3[1])/2.0;
                elemO2[i*nLocal+2*nLocalNode] = numNode;
                numNode++;
            }
        }


    AngleAxeVert = 140.0;
    AngleAxeHor = -40.0;
    y = 1.0f;
    int pause = 0;
    int d3 = 1;
    int maillage = 0;
    char start = 'Q'; // START
    char restart = 'R'; //restart
    char paus = 'D'; // PAUSE
    char zoom = 'Z'; // ZOOM
    char dezoom = 'S'; //DEZOOM
    char maillageActivate = 'W'; // Activation du maillage
    char maillageDesactivate = 'X'; // Desactivation maillage
    char d3Activate = 'A'; // Activation de la vague en 3d
    char d3Desactivate = 'E'; // Desactivation de la vague en 3d
    char vitesse1 = 'C';
    char vitesse2 = 'V';
    char vitesse3 = 'B';
    char position1 = 'T';
    char position2 = 'Y';
    char position3 = 'U';
    char textOn = 'O';
    char textOff = 'P';


    printf("========== Controls ==========\n");
    printf("START/CONTINUE: %c\n",start);
    printf("RESTART: %c\n",restart);
    printf("PAUSE : press %c\n",paus);
    printf("ZOOM +/- : press %c/%c\n",zoom,dezoom);
    printf("ROTATION : use the arrows\n");
    printf("SPEED 1/2/3 : press %c/%c/%c\n",vitesse1,vitesse2,vitesse3);
    printf("POSITION INIT/LEGAT/JAPAN : press %c/%c/%c\n",position1, position2, position3);
    printf("3D WAVE ON/OFF : press %c/%c\n",d3Activate,d3Desactivate);
    printf("MESH ON/OFF : press %c/%c\n",maillageActivate,maillageDesactivate);
    printf("TEXT IN WINDOW ON/OFF : press %c/%c\n",textOn,textOff);
    printf("QUIT : press ESC\n");
    printf("==============================\n\n");

    printf("========== Informations ==========\n");
    printf("Filename of the Mesh : %s\n",meshFileName);
    printf("Number of iterations : %d\n",nmax);
    printf("Number of frames : %d\n",nmax/nout+1);
    printf("Time step : %f\n",dt);
    printf("Order : %d\n",order);
    if (nLocalNode==3) { printf("TYPE : TRIANGLES\n"); }
    else { printf("TYPE : QUADS\n"); }
    printf("==================================\n\n");

    printf("===>   PRESS %c TO BEGIN\n\n\n",start);

    glfwInit();
   	glfwOpenWindow(840,480,0,0,0,0,1,0,GLFW_WINDOW );
	glfwSetWindowTitle( "MECA1120 Tsunami" );
	glShadeModel(GL_SMOOTH);

    glfwEnable( GLFW_STICKY_KEYS );
    glfwSwapInterval( 1 );

    GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 0.0 };
    GLfloat mat_shininess[] = { 50.0 };
    GLfloat light_position[] = { 8.0, 8.0, 8.0, 0.0 };

	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    GLfloat light_radiance[] = {1., 1., 1., 1.};

    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_radiance);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_radiance);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_NORMALIZE);

    int inc = 0;
    int fichier = 0;
    int vitesse = 1;
    int textOnScreen = 1;
    double tt = 0;
    char Message[256];

    glfemInit("");
    do
    {
    	if( pause==1 )
        {
            fichier = inc * nout;
            inc += vitesse;
        }
    	if (fichier>=nmax) { fichier = nmax; }

        glfwGetMousePos( &mouse, NULL );  mouse = 389;

        const char *basename = "%s-%08d.txt";
        char filename[256];
        sprintf(filename,basename,baseResultName,fichier);

        if ((inc-1 <= nmax/nout) && (pause==1))
        printf("===  Reading local file %s %d \n",filename,inc-1);
        {
            tsunamiReadFile(baseResultName,fichier,U,V,E,nElem);
        }

        if ((glfwGetKey(textOn) == GLFW_PRESS) && (textOnScreen!= 1)) // VITESSE
        {
            textOnScreen = 1;
            printf("===  TEXT IN WINDOW ON  ===\n");
        }
        if ((glfwGetKey(textOff) == GLFW_PRESS) && (textOnScreen!= 0)) // VITESSE
        {
            textOnScreen = 0;
            printf("===  TEXT IN WINDOW OFF  ===\n");
        }
        if ((glfwGetKey(vitesse1) == GLFW_PRESS) && (vitesse != 1)) // VITESSE
        {
            vitesse = 1;
            printf("===  SPEED 1 (x1)  ===\n");
        }
        if ((glfwGetKey(vitesse2) == GLFW_PRESS) && (vitesse != 2))
        {
            vitesse = 2;
            printf("===  SPEED 2 (x2)  ===\n");

        }
        if ((glfwGetKey(vitesse3) == GLFW_PRESS) && (vitesse != 4))
        {
            vitesse = 4;
            printf("===  SPEED 3 (x4)  ===\n");

        }
        if (glfwGetKey(position1) == GLFW_PRESS)
        {
            vitesse = 1;
            AngleAxeVert = 140.0;
            AngleAxeHor = -40.0;
            y = 1.0f;
            printf("===  INITIAL POSITION  ===\n");
        }
        if (glfwGetKey(position2) == GLFW_PRESS)
        {
            vitesse = 1;
            AngleAxeVert = 0.3f*(GLfloat)mouse + (GLfloat)tt*10.0f;
            AngleAxeHor = 0.0;
            y = 1.0f;
            printf("===  LEGAT POSITION  ===\n");
        }
        if (glfwGetKey(position3) == GLFW_PRESS) // japon
        {
            vitesse = 1;
            AngleAxeVert = 140.0;
            AngleAxeHor = -50.0;
            y = 1.0f+5.0;
            printf("===  JAPAN VIEW POSITION  ===\n");
        }
        if ((glfwGetKey(restart) == GLFW_PRESS) && (inc != 0))
        {
            inc = 0;
            fichier = 0;
            vitesse = 1;
            printf("===  RESTART  ===\n\n");
        }
        if ((glfwGetKey(paus) == GLFW_PRESS) && (pause != 0))
        {
            pause = 0;
            printf("===  PAUSE  ===\n");
        }
        if ((glfwGetKey(start) == GLFW_PRESS) && (pause != 1))
        {
            pause = 1;
            printf("===  START/CONTINUE  ===\n");
        }
        if ((glfwGetKey(d3Activate) == GLFW_PRESS) && (d3 != 1))
        {
            d3 = 1;
            printf("===  ACTIVATION 3D WAVE  ===\n");
        }
        if ((glfwGetKey(d3Desactivate) == GLFW_PRESS) && (d3 != 0))
        {
            d3 = 0;
            printf("===  DEACTIVATION 3D WAVE  ===\n");
        }
        if ((glfwGetKey(maillageActivate) == GLFW_PRESS) && (maillage != 1))
        {
            maillage = 1;
            printf("===  ACTIVATION MESH  ===\n");
        }
        if ((glfwGetKey(maillageDesactivate) == GLFW_PRESS) && (maillage != 0))
        {
            maillage = 0;
            printf("===  DEACTIVATION MESH  ===\n");
        }
        if (glfwGetKey(zoom) == GLFW_PRESS)
        {
            y += 1.0;
        }
        else if (glfwGetKey(dezoom) == GLFW_PRESS)
        {
            y -= 1.0;
        }
        if (glfwGetKey(GLFW_KEY_LEFT ) == GLFW_PRESS)
        {
            AngleAxeVert += 10.0;
        }
        else if (glfwGetKey(GLFW_KEY_RIGHT ) == GLFW_PRESS)
        {
            AngleAxeVert -= 10.0;
        }
        else if (glfwGetKey(GLFW_KEY_UP ) == GLFW_PRESS)
        {
            AngleAxeHor -= 10.0;
        }
        else if (glfwGetKey(GLFW_KEY_DOWN ) == GLFW_PRESS)
        {
            AngleAxeHor += 10.0;
        }


        glfwGetWindowSize( &width, &height );
        height = height > 0 ? height : 1;
        glViewport( 0, 0, width, height );

        glClearColor( 0.9f, 0.9f, 0.8f, 0.0f );
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        gluPerspective(65.0f,(GLfloat)width/(GLfloat)height,1.0f,100.0f);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        gluLookAt(0.0f,y,0.0f,0.0f, 20.0f, 0.0f,0.0f,0.0f,1.0f);
        glTranslatef(0.0f,14.0f,0.0f);
        glRotatef(AngleAxeVert,0.0f,0.0f,1.0f);
        glRotatef(AngleAxeHor,1.0f,0.0f,0.0f);


        GLUquadricObj *quadratic = gluNewQuadric();
        gluQuadricNormals(quadratic, GLU_SMOOTH);
        glColor3f(0.9,0.5,0.2);
        gluSphere(quadratic,5.95,50,50);

        if (textOnScreen==1)
        {
            glColor3f(1.0,0.0,0.0);
            glfemDrawMessage(10,60,"== Controls ==");
            sprintf(Message,"POSITION 1/2/3 : %c/%c/%c",position1, position2, position3);
            glfemDrawMessage(10,80,Message);
            sprintf(Message,"START/CONTINUE: %c",start);
            glfemDrawMessage(10,100,Message);
            sprintf(Message,"RESTART: %c",restart);
            glfemDrawMessage(10,120,Message);
            sprintf(Message,"PAUSE : %c",paus);
            glfemDrawMessage(10,140,Message);
            sprintf(Message,"ROTATION : arrows");
            glfemDrawMessage(10,160,Message);
            sprintf(Message,"ZOOM +/- : %c/%c",zoom,dezoom);
            glfemDrawMessage(10,180,Message);
            sprintf(Message,"SPEED 1/2/3 : %c/%c/%c",vitesse1,vitesse2,vitesse3);
            glfemDrawMessage(10,200,Message);
            sprintf(Message,"3D WAVE ON/OFF : %c/%c",d3Activate,d3Desactivate);
            glfemDrawMessage(10,220,Message);
            sprintf(Message,"MESH ON/OFF : %c/%c",maillageActivate,maillageDesactivate);
            glfemDrawMessage(10,240,Message);
            sprintf(Message,"TEXT ON/OFF : %c/%c",textOn,textOff);
            glfemDrawMessage(10,260,Message);
            sprintf(Message,"QUIT : ESC");
            glfemDrawMessage(10,280,Message);

            double tsuTime = (fichier*dt)/60.0;
            if (tsuTime <=60.0)
            {
                sprintf(Message,"Tsunami Time : %0.1f minutes", tsuTime);
            }
            else
            {
                sprintf(Message,"Tsunami Time : %0.0f hours %0.0f minutes", tsuTime/60.0, tsuTime - ( (int) (tsuTime/60.0))*60.0);
            }
            glfemDrawMessage(170,20,Message);
            glfemDrawMessage(350,60,"== Informations ==");
            if (nLocalNode==3)
            {
                sprintf(Message,"Type : TRIANGLES");
            }
            else
            {
                sprintf(Message,"Type : QUADS");
            }
            glfemDrawMessage(360,80,Message);
            sprintf(Message,"Order : %d",order);
            glfemDrawMessage(360,100,Message);
            sprintf(Message,"Time step : %0.1f",dt);
            glfemDrawMessage(360,120,Message);
            sprintf(Message,"Number iter. : %d",nmax);
            glfemDrawMessage(360,140,Message);
            sprintf(Message,"Number frames : %d",nmax/nout+1);
            glfemDrawMessage(360,160,Message);
        }

        for (i=0; i < nElem; ++i)
        {
            if(order == 1)
            {
                for (j=0; j < nLocalNode; ++j)
                {

                    index = elem[nLocalNode*i+j];
                    double value = E[nLocalNode*i+j]*10;
                    if (value < 0) value = 0;
                    if (value > 1) value = 1;
                    colors[j*3+0] = 3.5*(value)*(value);
                    colors[j*3+1] = (1-value)*(value)*3.5;
                    colors[j*3+2] = (1-value)*(1-value);
                    double x = XX[index];
                    double y = YY[index];
                    double Factor = (4*R*R + x*x + y*y)*(R/6);
                    coord[j*3+0] = 4*R*R * x / Factor;
                    coord[j*3+1] = 4*R*R * y / Factor;
                    coord[j*3+2] = (4*R*R - x*x - y*y)*R / Factor;
                    if (d3 == 1 )
                    {
                        double norm = sqrt(coord[j*3+0] * coord[j*3+0] + coord[j*3+1] * coord[j*3+1] + coord[j*3+2] * coord[j*3+2]);
                        double h = value * 0.1;
                        coordBis[j*3+0] = coord[j*3+0] + coord[j*3+0] * h/norm;
                        coordBis[j*3+1] = coord[j*3+1] + coord[j*3+1] * h/norm;
                        coordBis[j*3+2] = coord[j*3+2] + coord[j*3+2] * h/norm;
                    }
                }
                glEnableClientState(GL_VERTEX_ARRAY);
                glEnableClientState(GL_COLOR_ARRAY);
                glEnableClientState(GL_NORMAL_ARRAY);
                glVertexPointer(3, GL_FLOAT, 0, coord);
                glNormalPointer(GL_FLOAT, 0, coord);
                glColorPointer(3, GL_FLOAT, 0, colors);
                glDrawArrays(GL_POLYGON, 0, nLocalNode);
                if (d3 == 1 )
                {
                    glVertexPointer(3, GL_FLOAT, 0, coordBis);
                    glNormalPointer(GL_FLOAT, 0, coordBis);
                    glDrawArrays(GL_POLYGON, 0, nLocalNode);
                }
                glDisableClientState(GL_NORMAL_ARRAY);
                glDisableClientState(GL_COLOR_ARRAY);
                glDisableClientState(GL_VERTEX_ARRAY);

                if (maillage == 1)
                {
                    glColor3f(0.0, 0.0, 0.0);
                    glEnableClientState(GL_VERTEX_ARRAY);
                    for (j=0; j < size; ++j)
                         coord[j] = coord[j] * 1.001;
                    glVertexPointer(3, GL_FLOAT, 0, coord);
                    glDrawArrays(GL_LINE_LOOP, 0, nLocalNode);
                    glDisableClientState(GL_VERTEX_ARRAY);
                }
            }
            else
            {
                int *tabNode = malloc(nLocal * sizeof(int));
                int numNode[16];
                int newTri[4][4];
                for (j=0 ; j<nLocal ; j++ )
                {
                    tabNode[j] = elemO2[nLocal*i+j];
                }
                if (nLocalNode==3)
                {
                    numNode[0] = 0; numNode[1] = 3; numNode[2] = 5;
                    numNode[3] = 3; numNode[4] = 1; numNode[5] = 4;
                    numNode[6] = 3; numNode[7] = 4; numNode[8] = 5;
                    numNode[9] = 5; numNode[10] = 4; numNode[11] = 2;
                    newTri[0][0] = tabNode[0]; newTri[0][1] = tabNode[3]; newTri[0][2] = tabNode[5];
                    newTri[1][0] = tabNode[3]; newTri[1][1] = tabNode[1]; newTri[1][2] = tabNode[4];
                    newTri[2][0] = tabNode[3]; newTri[2][1] = tabNode[4]; newTri[2][2] = tabNode[5];
                    newTri[3][0] = tabNode[5]; newTri[3][1] = tabNode[4]; newTri[3][2] = tabNode[2];
                }
                else
                {
                    numNode[0] = 0; numNode[1] = 4; numNode[2] = 8; numNode[3] = 7;
                    numNode[4] = 4; numNode[5] = 1; numNode[6] = 5; numNode[7] = 8;
                    numNode[8] = 8; numNode[9] = 5; numNode[10] = 2; numNode[11] = 6;
                    numNode[12] = 7; numNode[13] = 8; numNode[14] = 6; numNode[15] = 3;
                    newTri[0][0] = tabNode[0]; newTri[0][1] = tabNode[4]; newTri[0][2] = tabNode[8]; newTri[0][3] = tabNode[7];
                    newTri[1][0] = tabNode[4]; newTri[1][1] = tabNode[1]; newTri[1][2] = tabNode[5]; newTri[1][3] = tabNode[8];
                    newTri[2][0] = tabNode[8]; newTri[2][1] = tabNode[5]; newTri[2][2] = tabNode[2]; newTri[2][3] = tabNode[6];
                    newTri[3][0] = tabNode[7]; newTri[3][1] = tabNode[8]; newTri[3][2] = tabNode[6]; newTri[3][3] = tabNode[3];
                }
                for( j=0 ; j<4 ; j++ )
                {
                    int k;
                    for( k=0 ; k<nLocalNode ; k++ )
                    {
                        index = newTri[j][k];
                        double value = E[nLocal*i+numNode[j*nLocalNode+k]]*10;
                        if (value < 0) value = 0;
                        if (value > 1) value = 1;
                        colors[k*3+0] = 3.5*(value)*(value);
                        colors[k*3+1] = (1-value)*(value)*3.5;
                        colors[k*3+2] = (1-value)*(1-value);
                        double x = XX[index];
                        double y = YY[index];
                        double Factor = (4*R*R + x*x + y*y)*(R/6);
                        coord[k*3+0] = 4*R*R * x / Factor;
                        coord[k*3+1] = 4*R*R * y / Factor;
                        coord[k*3+2] = (4*R*R - x*x - y*y)*R / Factor;
                        if (d3 == 1)
                        {
                            double norm = sqrt(coord[k*3+0] * coord[k*3+0] + coord[k*3+1] * coord[k*3+1] + coord[k*3+2] * coord[k*3+2]);
                            double h = value * 0.1;
                            coordBis[k*3+0] = coord[k*3+0] + coord[k*3+0] * h/norm;
                            coordBis[k*3+1] = coord[k*3+1] + coord[k*3+1] * h/norm;
                            coordBis[k*3+2] = coord[k*3+2] + coord[k*3+2] * h/norm;
                        }
                    }
                    glEnableClientState(GL_VERTEX_ARRAY);
                    glEnableClientState(GL_COLOR_ARRAY);
                    glEnableClientState(GL_NORMAL_ARRAY);
                    glVertexPointer(3, GL_FLOAT, 0, coord);
                    glNormalPointer(GL_FLOAT, 0, coord);
                    glColorPointer(3, GL_FLOAT, 0, colors);
                    glDrawArrays(GL_POLYGON, 0, nLocalNode);
                    if(d3 == 1)
                    {
                        glVertexPointer(3, GL_FLOAT, 0, coordBis);
                        glNormalPointer(GL_FLOAT, 0, coordBis);
                        glDrawArrays(GL_POLYGON, 0, nLocalNode);
                    }
                    glDisableClientState(GL_NORMAL_ARRAY);
                    glDisableClientState(GL_COLOR_ARRAY);
                    glDisableClientState(GL_VERTEX_ARRAY);
                }
                if( maillage == 1 )
                {
                    for (j=0; j < nLocalNode; ++j)
                    {
                        index = elemO2[nLocal*i+j];
                        double x = XX[index];
                        double y = YY[index];
                        double Factor = (4*R*R + x*x + y*y)*(R/6);
                        coord[j*3+0] = 4*R*R * x / Factor;
                        coord[j*3+1] = 4*R*R * y / Factor;
                        coord[j*3+2] = (4*R*R - x*x - y*y)*R / Factor;
                    }
                    glColor3f(0.0, 0.0, 0.0);
                    glEnableClientState(GL_VERTEX_ARRAY);
                    for (j=0; j < size; ++j)
                        coord[j] = coord[j] * 1.001;
                    glVertexPointer(3, GL_FLOAT, 0, coord);
                    glDrawArrays(GL_LINE_LOOP, 0, nLocalNode);
                    glDisableClientState(GL_VERTEX_ARRAY);
                }
                free(tabNode);
            }
        }
        glfwSwapBuffers();

    }
    while( glfwGetKey( GLFW_KEY_ESC ) != GLFW_PRESS && glfwGetWindowParam( GLFW_OPENED ) );

    free(U); free(V); free(E);
    free(XX); free(YY); free(elemO2);
    femMeshFree(theMesh);
    femEdgesFree(edgesStruc);

    glfwTerminate();
    exit( EXIT_SUCCESS );
}



double tsunamiInitialConditionOkada(double x, double y)
{
    double R = 6371220.0;
    double x3d = 4*R*R*x / (4*R*R + x*x + y*y);
    double y3d = 4*R*R*y / (4*R*R + x*x + y*y);
    double z3d = R*(4*R*R - x*x - y*y) / (4*R*R + x*x + y*y);
    double lat = asin(z3d/R)*180/M_PI;
    double lon = atan2(y3d,x3d)*180/M_PI;
    double lonMin = 142;
    double lonMax = 143.75;
    double latMin = 35.9;
    double latMax = 39.5;
    double olon = (lonMin+lonMax)/2;
    double olat = (latMin+latMax)/2;
    double angle = -12.95*M_PI/180;
    double lon2 = olon + (lon-olon)*cos(angle) + (lat-olat)*sin(angle);
    double lat2 = olat - (lon-olon)*sin(angle) + (lat-olat)*cos(angle);
    if ( lon2 <= lonMax && lon2 >= lonMin &&
         lat2 >= latMin && lat2 <= latMax )
            return 1.0;
    else    return 0.0;
}

void tsunamiWriteFile(const char *baseResultName, int iter, double *U, double *V, double *E, int nelem, int nsub)
{
    int i,j;
    const char *basename = "%s-%08d.txt";
    char filename[256];
    sprintf(filename,basename,baseResultName,iter);
    FILE* file = fopen(filename,"w");
    fprintf(file, "Number of elem %d \n", nelem);
    fprintf(file, "Number of local values per element %d \n", nsub);
    for (i = 0; i < nelem; ++i) {
    	for (j = 0; j < nsub; ++j) {
        	int index = i*nsub+j;
        	fprintf(file,"%d;%d;%le;%le;%le;\n",i,j,U[index],V[index],E[index]); }}
    fclose(file);
}

int tsunamiReadFile(const char *baseResultName, int iter, double *U, double *V, double *E, int nelem)
{
    int i,j,trash,nelemFile,nsub;
    const char *basename = "%s-%08d.txt";
    char filename[256];
    sprintf(filename,basename,baseResultName,iter);
    FILE* file = fopen(filename,"r");
    fscanf(file, "Number of elem %d \n", &nelemFile);
    fscanf(file, "Number of local values per element %d \n", &nsub);
    if (nelem != nelemFile) {
        printf("Error : wrong data file %d %d:-) \n",nelem,nelemFile);
        exit(0); }
    for (i = 0; i < nelem; ++i) {
    	for (j = 0; j < nsub; ++j) {
        	int index = i*nsub+j;
        	fscanf(file,"%d;%d;%le;%le;%le;\n",&trash,&trash,&U[index],&V[index],&E[index]); }}

    fclose(file);
    return nsub;
}

//--------------------------------------------------------------------------------//devoir 9


femTsunamiProblem *femTsunamiCreate (const char *meshFileName)
{

    femTsunamiProblem *myProblem = malloc(sizeof(femTsunamiProblem));
    myProblem->omega   = (2*M_PI)/86400.0;
    myProblem->gravity  = 9.81;
    myProblem->gamma    = 1e-7;
    myProblem->R     = 6371220.0;


    myProblem->mesh = femMeshRead(meshFileName);
    myProblem->edges = femEdgesCreate(myProblem->mesh);


    myProblem->rule1d = femIntegrationCreate(2,FEM_EDGE);


    if (myProblem->mesh->nLocalNode == 4) {
        myProblem->space = femDiscreteCreate(4,FEM_QUAD);
        myProblem->rule2d = femIntegrationCreate(9,FEM_QUAD);  }


    else if (myProblem->mesh->nLocalNode == 3) {
        myProblem->space = femDiscreteCreate(3,FEM_TRIANGLE);
        myProblem->rule2d = femIntegrationCreate(3,FEM_TRIANGLE); }


    int sizeLoc = myProblem->space->n;
    int sizeGlo = myProblem->mesh->nElem * sizeLoc + 1;

    //myProblem->solver = femSolverFullCreate(3*sizeLoc);
    //permet de gagner du temps

    myProblem->solver = femSolverBandCreate(3*sizeLoc,sizeLoc);
    myProblem->height = myProblem->mesh->H;
    myProblem->size = sizeGlo;

    myProblem->H = malloc(sizeof(double)*sizeGlo);
    myProblem->U = malloc(sizeof(double)*sizeGlo);
    myProblem->V = malloc(sizeof(double)*sizeGlo);
    myProblem->FH = malloc(sizeof(double)*sizeGlo);
    myProblem->FU = malloc(sizeof(double)*sizeGlo);
    myProblem->FV = malloc(sizeof(double)*sizeGlo);
    int i;
    for (i=0; i < myProblem->size; i++)
    {
        if ( i == ((myProblem->size)-1) )
        {
            myProblem->H[i] = 0.0;
            myProblem->U[i] = 0.0;
            myProblem->V[i] = 0.0;
        }
        else
        {
            int node = myProblem->mesh->elem[i];
            myProblem->H[i] = tsunamiInitialConditionOkada(myProblem->mesh->X[node], myProblem->mesh->Y[node]);
            myProblem->U[i] = 0.0;
            myProblem->V[i] = 0.0;
        }
    }
    return myProblem;

}

void femTsunamiFree(femTsunamiProblem *myProblem)
{
 free(myProblem->FH);
    free(myProblem->FU);
    free(myProblem->FV);
    free(myProblem->H);
    free(myProblem->U);
    free(myProblem->V);
    femSolverFree(myProblem->solver);
    femIntegrationFree(myProblem->rule1d);
    femIntegrationFree(myProblem->rule2d);
    femDiscreteFree(myProblem->space);
    femEdgesFree(myProblem->edges);
    femMeshFree(myProblem->mesh);
    free(myProblem);

}

double interpolate(double *phi, double *U, int *map, int n)
{
    double u = 0.0; int i;
    for (i=0; i <n; i++)
        u += phi[i]*U[map[i]];
    return u;
}

double interpolatebis(double *phi, double *U, int n)
{
    double u =0.0; int i;
    for(i=0; i<n; i++)
	u+=phi[i]*U[i];

    return u;
}
void femTsunamiAddIntegralsElements(femTsunamiProblem *myProblem)
{
    double *BH = myProblem->FH;
    double *BU = myProblem->FU;
    double *BV = myProblem->FV;
    double *E = myProblem->H;
    double *U = myProblem->U;
    double *V = myProblem->V;
    double *Y = myProblem->mesh->Y;
    double *X = myProblem->mesh->X;
    double *H = myProblem->mesh->H;

    femIntegration *theRule = myProblem->rule2d;
    femDiscrete *theSpace = myProblem->space;

    int n = theSpace->n;
    double Xloc[n],Yloc[n],phi[n],dphidxsi[n],dphideta[n],dphidx[n],dphidy[n],Hloc[n];
    double  xsi,eta,weight,jac;
    double  y,e,u,v,x,h;
    int     i,j,k,elem,mapElem[n];

    double R = myProblem->R;
    double g     = myProblem->gravity;
    double gamma = myProblem->gamma;
    double omega = myProblem->omega;

    for (elem=0; elem < myProblem->mesh->nElem; elem++) {
        int *mapCoord = &(myProblem->mesh->elem[elem*n]);
        for (j=0; j < n; ++j) {
            mapElem[j] = elem*n + j;
        	Xloc[j] = X[mapCoord[j]];
        	Yloc[j] = Y[mapCoord[j]];
		Hloc[j] = H[mapCoord[j]];
	}

        for (k=0; k < theRule->n; k++) {
            xsi = theRule->xsi[k];
            eta = theRule->eta[k];
            weight = theRule->weight[k];
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0;
            double dydeta = 0.0;
            for (i = 0; i < n; i++) {
                dxdxsi += Xloc[i]*dphidxsi[i];
                dxdeta += Xloc[i]*dphideta[i];
                dydxsi += Yloc[i]*dphidxsi[i];
                dydeta += Yloc[i]*dphideta[i]; }
            jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            for (i = 0; i < n; i++) {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; }
            x=interpolatebis(phi, Xloc, n);
	    y=interpolatebis(phi,Yloc,n);
	    h = interpolatebis(phi,Hloc, n);


            e = interpolate(phi,E,mapElem,n);
            u = interpolate(phi,U,mapElem,n);
            v = interpolate(phi,V,mapElem,n);

  	  double theta = asin((4.0 * R * R - x * x - y * y)/(4.0 * R * R + x * x + y * y));
  	  double f =2.0*omega*sin(theta);

            for (i=0; i < n; i++) {

                BH[mapElem[i]] += ((h*u*dphidx[i] + h*v*dphidy[i])*jac*weight* (4.0*R*R +x*x + y*y)/(4.0*R*R))  ;

		BH[mapElem[i]] += (phi[i]*(h*(u*x + y*v)/(R*R))*jac*weight);

                BU[mapElem[i]] += ( (f*v -gamma*u)*phi[i] + g*e*dphidx[i]*(4.0*R*R +x*x + y*y)/(4.0*R*R))*jac*weight;

		BU[mapElem[i]] += (phi[i]*jac*weight*(g*x*e/(R*R*2.0)));

                BV[mapElem[i]] += ( (-f*u-gamma*v)*phi[i] + g*e*dphidy[i]*(4.0*R*R +x*x + y*y)/(4.0*R*R) )*jac*weight ;

		BV[mapElem[i]] += phi[i]*jac*weight*(g*y*e/(R*R*2.0)); }}}

}



void femTsunamiAddIntegralsEdges(femTsunamiProblem *myProblem)
{
    double  xEdge[2],yEdge[2],phiEdge[2],hEdge[2];
    double  xsi,weight,jac;
    double  eL,eR,uL,uR,vL,vR,unL,unR;
    double  qe,qu,qv;
    int     i,j,k,edge,mapEdge[2][2];

    double *BH = myProblem->FH;
    double *BU = myProblem->FU;
    double *BV = myProblem->FV;
    double *E = myProblem->H;
    double *U = myProblem->U;
    double *V = myProblem->V;
    femIntegration *theRule = myProblem->rule1d;
    femDiscrete *theSpace = myProblem->space;
    int     sizeLoc = theSpace->n;
    int     sizeGlo = myProblem->mesh->nElem * sizeLoc + 1;

    double  g = myProblem->gravity;
    double R = myProblem->R;

    for (edge=0; edge < myProblem->edges->nEdge; edge++) {
        femEdgesMap(myProblem->edges,edge,mapEdge);
        for (j=0; j < 2; ++j) {
        	int node = myProblem->edges->edges[edge].node[j];
        	xEdge[j] = myProblem->mesh->X[node];
        	yEdge[j] = myProblem->mesh->Y[node];
		hEdge[j] = myProblem->mesh->H[node]; }

        int boundary = (mapEdge[1][0] == sizeGlo-1);

        double dxdxsi = (xEdge[1] - xEdge[0]);
        double dydxsi = (yEdge[1] - yEdge[0]);
        double norm = sqrt(dxdxsi*dxdxsi + dydxsi*dydxsi);
        double nx =  dydxsi/norm;
        double ny = -dxdxsi/norm;
        jac = norm / 2.0;
        for (k=0; k < theRule->n; k++) {
            xsi = theRule->xsi[k];
            weight = theRule->weight[k];
            femDiscretePhi1(theSpace,xsi,phiEdge);

            eL = interpolate(phiEdge,E,mapEdge[0],2);
            eR = boundary ? eL : interpolate(phiEdge,E,mapEdge[1],2);
            uL = interpolate(phiEdge,U,mapEdge[0],2);
            uR = interpolate(phiEdge,U,mapEdge[1],2);
            vL = interpolate(phiEdge,V,mapEdge[0],2);
            vR = interpolate(phiEdge,V,mapEdge[1],2);
            unL = uL*nx+ vL*ny;
            unR = boundary ? -unL : uR*nx + vR*ny;

	    double h=hEdge[0]*phiEdge[0] + hEdge[1]*phiEdge[1];

	    double x=xEdge[0]*phiEdge[0] + xEdge[1]*phiEdge[1];
	    double y=yEdge[0]*phiEdge[0] + yEdge[1]*phiEdge[1];

            qe = ((4.0*R*R + x*x + y*y)/(4.0*R*R))*0.5*h*   ( (unL+unR) + sqrt(g/h)*( eL-eR ) );
            qu = ((4.0*R*R + x*x + y*y)/(4.0*R*R))*0.5*g*nx*( ( eL+eR ) + sqrt(h/g)*(unL-unR) );
            qv = ((4.0*R*R + x*x + y*y)/(4.0*R*R))*0.5*g*ny*( ( eL+eR ) + sqrt(h/g)*(unL-unR) );

            for (i=0; i < 2; i++) {
                BH[mapEdge[0][i]] -= qe*phiEdge[i]*jac*weight;
                BU[mapEdge[0][i]] -= qu*phiEdge[i]*jac*weight;
                BV[mapEdge[0][i]] -= qv*phiEdge[i]*jac*weight;
                BH[mapEdge[1][i]] += qe*phiEdge[i]*jac*weight;
                BU[mapEdge[1][i]] += qu*phiEdge[i]*jac*weight;
                BV[mapEdge[1][i]] += qv*phiEdge[i]*jac*weight; }}}

}

void femTsunamiMultiplyInverseMatrix(femTsunamiProblem *myProblem)
{
    double *BH = myProblem->FH;
    double *BU = myProblem->FU;
    double *BV = myProblem->FV;
    femMesh *theMesh = myProblem->mesh;
    femDiscrete *theSpace = myProblem->space;
    femSolver *theSolver = myProblem->solver;
    femIntegration *theRule = myProblem->rule2d;

    int n = theSpace->n;
    double Xloc[n],Yloc[n],phi[n],dphidxsi[n],dphideta[n],Aloc[n*n],jac;
    int iElem,iInteg,i,j,mapElem[n],mapE[n],mapU[n],mapV[n];

    for (i = 0; i < n; i++)   {
        mapE[i] = i;
        mapU[i] = i + n;
        mapV[i] = i + 2*n; }

    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        femSolverInit(theSolver);
        for (i = 0; i < n*n; i++)  Aloc[i] = 0;
        int *mapCoord = &(myProblem->mesh->elem[iElem*n]);
        for (j=0; j < n; ++j) {
            mapElem[j] = iElem*n + j;
        	Xloc[j] = myProblem->mesh->X[mapCoord[j]];
        	Yloc[j] = myProblem->mesh->Y[mapCoord[j]]; }

        for (iInteg=0; iInteg < theRule->n; iInteg++) {
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            double dxdxsi = 0;
            double dxdeta = 0;
            double dydxsi = 0;
            double dydeta = 0;
            for (i = 0; i < n; i++) {
                dxdxsi += Xloc[i]*dphidxsi[i];
                dxdeta += Xloc[i]*dphideta[i];
                dydxsi += Yloc[i]*dphidxsi[i];
                dydeta += Yloc[i]*dphideta[i]; }
            jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            for (i = 0; i < n; i++) {
                for(j = 0; j < n; j++) {
                    Aloc[i*(theSpace->n)+j] += phi[i] * phi[j] * jac * weight; }}}

        femSolverAssemble(theSolver,Aloc,&BH[mapElem[0]],0,mapE,theSpace->n);
        femSolverAssemble(theSolver,Aloc,&BU[mapElem[0]],0,mapU,theSpace->n);
        femSolverAssemble(theSolver,Aloc,&BV[mapElem[0]],0,mapV,theSpace->n);
        double *soluce = femSolverEliminate(theSolver);
     	for (i = 0; i < n; i++) {
            BH[mapElem[i]] = soluce[mapE[i]];//i
            BU[mapElem[i]] = soluce[mapU[i]];//i+n
            BV[mapElem[i]] = soluce[mapV[i]];//i+ 2*n
	}}

}


void femTsunamiCompute(femTsunamiProblem *myProblem)
{
    int  size = myProblem->size, i;
    double *FE = myProblem->FH;
    double *FU = myProblem->FU;
    double *FV = myProblem->FV;
    double *E = myProblem->H;
    double *U = myProblem->U;
    double *V = myProblem->V;
    double dt = myProblem->dt;
    double Eold[size],Uold[size],Vold[size];
    double Enew[size],Unew[size],Vnew[size];

//le calcul selon Runge-Kutta d'ordre 4 a été partagé par Florentin Goyens
  for(i=0; i<size; i++)
                    {
                    	Eold[i]=E[i];
                    	Uold[i]=U[i];
                    	Vold[i]=V[i];
        		FE[i] = 0.0;
        		FU[i] = 0.0;
        		FV[i] = 0.0; }
	   	        femTsunamiAddIntegralsElements(myProblem);
    			femTsunamiAddIntegralsEdges(myProblem);
    			femTsunamiMultiplyInverseMatrix(myProblem);
                    for (i=0; i < size; i++) {
        	 	Enew[i] = E[i]+dt*FE[i]/6.0;
        		Unew[i] = U[i]+dt*FU[i]/6.0;
        		Vnew[i] = V[i]+dt*FV[i]/6.0;

        		E[i]=Eold[i]+dt*0.5*FE[i];
        		U[i]=Uold[i]+dt*0.5*FU[i];
        		V[i]=Vold[i]+dt*0.5*FV[i];
        		FE[i] = 0.0;
        		FU[i] = 0.0;
        		FV[i] = 0.0;
        		}
	   	        femTsunamiAddIntegralsElements(myProblem);
    			femTsunamiAddIntegralsEdges(myProblem);
    			femTsunamiMultiplyInverseMatrix(myProblem);
                        for (i=0; i < size; i++)
                        {
        		Enew[i] += 2*dt*FE[i]/6.0;
        		Unew[i] += 2*dt*FU[i]/6.0;
        		Vnew[i] += 2*dt*FV[i]/6.0;

        		E[i]=Eold[i]+dt*0.5*FE[i];
        		U[i]=Uold[i]+dt*0.5*FU[i];
        		V[i]=Vold[i]+dt*0.5*FV[i];
        		FE[i] = 0.0;
        		FU[i] = 0.0;
        		FV[i] = 0.0;
        		}
	   	        femTsunamiAddIntegralsElements(myProblem);
    			femTsunamiAddIntegralsEdges(myProblem);
    			femTsunamiMultiplyInverseMatrix(myProblem);
                        for (i=0; i < size; i++)
                        {

        		Enew[i] += 2*dt*FE[i]/6.0;
        		Unew[i] += 2*dt*FU[i]/6.0;
        		Vnew[i] += 2*dt*FV[i]/6.0;

        		E[i]=Eold[i]+dt*FE[i];
        		U[i]=Uold[i]+dt*FU[i];
        		V[i]=Vold[i]+dt*FV[i];
        		FE[i] = 0.0;
        		FU[i] = 0.0;
        		FV[i] = 0.0;
        		}
	   	        femTsunamiAddIntegralsElements(myProblem);
    			femTsunamiAddIntegralsEdges(myProblem);
    			femTsunamiMultiplyInverseMatrix(myProblem);
                    for(i=0;i<size; i++)
                    {
                    	E[i]=Enew[i]+dt*FE[i]/6.0;
                    	U[i]=Unew[i]+dt*FU[i]/6.0;
                    	V[i]=Vnew[i]+dt*FV[i]/6.0;
                    }

}
