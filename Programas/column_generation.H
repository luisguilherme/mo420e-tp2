#ifndef __COLUMN_GENERATION_H__
#define __COLUMN_GENERATION_H__

#define MAX_CPU_TIME -1800 //trocar
#define MAX_POOL 20 //trocar
#define EPSILON 1e-9 //trocar

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <vector>
#include <string>
#include "instance.H"
#include "integer_program.H"
#include "xprs.h"

typedef struct {
  /* para pool de colunas */
  std::vector<std::vector<double> > pool; /* pool de colunas geradas pelo pricing */
  int totPool;/* total de colunas inseridas no Pool de colunas */
  /* usados na otimizacao do problema mestre e do pricing */
  std::vector<double> xstar;/* armazena a melhor solucao inteira encontrada pelo  pricing */
  double zstar;  /* armazena o valor da solucao otima do pricing */
} solution;

class ColumnGeneration {
  double* dual;
  solution sol;
  clock_t t1, t2;
  XPRSprob probMestre;
  XPRSprob probPricing;
  int tPMR;

  //so pra compilar, depois vejo
  int xpress_ret;
  double tempo, tempoLP, tempoPricing;
  double z_PMR;
  double melhorPrimal;
  int totalInteiros;
  int itmelhorPrimal;
  int it;
  int* mindex;

protected:
  inline bool isIntegerSol(XPRSprob prob);

public:
  ColumnGeneration(IntegerProgram& ip, IntegerProgram& pricing);
  int loadModel(IntegerProgram&, XPRSprob, std::string);
  void solveRestricted();
  void solvePricing(Instance& instance);
  void configureModel(int, IntegerProgram& ip, IntegerProgram& pricing);
  //void detSolvePricing();
};
void errormsg(const char*, int, int, XPRSprob);
void XPRS_CC salvaMelhorSol(XPRSprob prob, void *psol);
#endif
