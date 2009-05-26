/* ******************************************************************
   Arquivo: cutting.h
   Descrição: estruturas e cabecalhos das funcoes de cutting.c
 * *************************************************************** */

#ifndef H_CUTTING
#define H_CUTTING 1

/*  Includes para  uso do  XPRESS */
#include "xprs.h"

/* definicões do usuário */
typedef enum{false,true} Boolean;
#define TRUE true;
#define FALSE false;
#define True true;
#define False false;
#define EPSILON 0.000001 /* uma tolerância usada em alguns cálculos com doubles*/

/* limitantes */
#define MAX_CPU_TIME -1800  /* limitante de tempo: colocar negativo ! (veja manual do XPRESS) */
#define MAX_POOL 50  /* total de entradas no pool de colunas  */

/* estruturas */

typedef struct {
  /* dados da instancia do problema de corte unidimensional */
  int m; /* numero de itens */
  int n; /* numero de padroes */
  int L; /* tamanho de cada barra */
  int *w; /* comprimento de cada item */
  int *d; /* demanda de cada item */
} cuttingInstance;

typedef struct{
  /* para pool de colunas */
  double *pool[MAX_POOL]; /* pool de colunas geradas pelo pricing */
  int totPool;/* total de colunas inseridas no Pool de colunas */
  /* usados na otimizacao do problema mestre e do pricing */
  double *xstar;/* armazena a melhor solucao inteira encontrada pelo  pricing */
  double zstar;  /* armazena o valor da solucao otima do pricing */
} solution;

/* rotinas auxiliares */
void errormsg(const char *sSubName,int nLineNo,int nErrCode, XPRSprob);
void showUsage();
int CargaPricing(XPRSprob *, cuttingInstance *);
int CargaProblemaMestre(XPRSprob *probMestre, cuttingInstance *cutting);
int IsIntegerSol(XPRSprob prob, double *z_PMR, double dual[]);
void adicionaCols(XPRSprob prob, double *pool[], int totPool, cuttingInstance *);
int CargaInstancia(char *instancia, cuttingInstance *);
void CriaCarregaConfigura(int formato, XPRSprob *, XPRSprob *, cuttingInstance *, solution *);
void ImprimeCols(double *pool[], int totpool, cuttingInstance *);

/* declaracoes das callbacks */
void XPRS_CC SalvaMelhorSol(XPRSprob prob, void *psol);

#endif
