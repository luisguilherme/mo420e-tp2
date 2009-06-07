#include "uls.H"
#include <cstdlib>

ULS::ULS(const ULSInstance& instance) {
  int M;

  this->instance = instance;

  nrows = 2 * instance.t;
  ncols = instance.t;

  M = 0;
  for (int i = 0; i < instance.t; i++)
    M += instance.d[i];
  M *= 2;
    
  columns = std::vector<std::vector<double> >
    (ncols, std::vector<double>(nrows, 0.0));
    
  for (int i = 0; i < ncols; i++) {
    if (i < instance.t) {
      columns[i][i] = 1;
      columns[i][instance.t + i] = 1;
    } else if (instance.t <= i && i < 2*instance.t) {
      columns[i][i] = -1;
      if (i < 2*instance.t - 1)
	columns[i+1][i] = 1;
    } else if (i >= 2*instance.t)
      columns[i][i - instance.t] = -M;
  }

}

void ULS::getParam(int& ncol, int& nrow, char** rowtype,double** rhs,
		   double** obj, int** colbeg, int** rowidx,
		   double** matval, double** lb, double** ub,
		   int& nmip, char** miptype, int** mipcol, bool& relaxed) {

  relaxed = false;

  /* Aloca memória dinamicamente para o XPress */
  nrow = nrows;
  ncol = ncols;

  *rowtype = (char *) calloc(nrow, sizeof(char));
  *rhs = (double *) calloc(nrow, sizeof(double));
  *obj = (double *) calloc(ncol, sizeof(double));
  *colbeg = (int *) calloc(ncol+1, sizeof(int));
  *lb = (double *) calloc(ncol, sizeof(double));
  *ub = (double *) calloc(ncol, sizeof(double));

  // configura tipo de restrição
  for (int i = 0; i < instance.t; i++)
    *rowtype[i] = 'E';
  for (int i = instance.t; i < nrows; i++)
    *rowtype[i] = 'L';

  // configura limites das variáveis
  for (int i = 0; i < 2*instance.t; i++) {
    *lb[i] = 0.0;
    *ub[i] = XPRS_PLUSINFINITY;
  }
  for (int i = 2*instance.t; i < ncols; i++) {
    *lb[i] = 0;
    *ub[i] = 1;
  }
    
  // aloca matriz esparsa
  int matsize = 0;
  for (int i = 0; i < ncols; i++)
    for (int j = 0; j < nrows; j++) matsize += NONZERO(columns[i][j]) ? 1 : 0;

  *matval = (double *) calloc(matsize, sizeof(double));
  *rowidx = (int *) calloc(matsize, sizeof(int));

  // configura colunas
  int off = 0;
  for (int i = 0; i < ncols; i++) {
    *colbeg[i] = off;
    for (int j = 0; j < nrows; j++)
      if (NONZERO(columns[i][j])) {
	*matval[off] = columns[i][j];
	*rowidx[off++] = j;
      }
  }

  // configura variaveis binárias
  nmip = instance.t;
  *miptype = (char *) calloc(nmip, sizeof(char));
  *mipcol = (int *) calloc(nmip, sizeof(int));
  for (int i = 0; i < nmip; i++) {
    *miptype[i] = 'B';
    *mipcol[i] = 2*instance.t + 1;
  }

  // configura rhs
  for (int i = 0; i < instance.t; i++)
    *rhs[i] = instance.d[i];
  for (int i = instance.t; i < nrows; i++)
    *rhs[i] = 0.0;

  // configura funcao objetivo
  for (int i = 0; i < instance.t; i++)
    *obj[i] = instance.p[i];
  for (int i = instance.t; i < 2*instance.t; i++)
    *obj[i] = instance.h[i-instance.t];
  for (int i = 2*instance.t; i < ncols; i++)
    *obj[i] = instance.f[i-2*instance.t];
}
