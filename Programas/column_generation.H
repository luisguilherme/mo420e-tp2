#ifndef __COLUMN_GENERATION_H__
#define __COLUMN_GENERATION_H__

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <vector>
#include "instance.H"
#include "xprs.h"

typedef struct {
  /* para pool de colunas */
  std::vector<double> pool; /* pool de colunas geradas pelo pricing */
  int totPool;/* total de colunas inseridas no Pool de colunas */
  /* usados na otimizacao do problema mestre e do pricing */
  std::vector<double> xstar;/* armazena a melhor solucao inteira encontrada pelo  pricing */
  double zstar;  /* armazena o valor da solucao otima do pricing */
} solution;

class ColumnGeneration {
  std::vector<double> dual;
  solution sol;
  clock_t t1, t2;
  XPRSprob probMestre;
  XPRSprob probPricing;

  //so pra compilar, depois vejo
  int xpress_ret;
  double tempo, tempoLP;
  double z_PMR;

public:
  void solveReduced();
  void solvePricing(Instance& instance);
  //void detSolvePricing();
};

#endif
