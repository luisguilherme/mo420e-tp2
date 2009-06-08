#include "uls.H"
#include <cstdlib>

ULS::ULS(const ULSInstance& instance) {
  int M;

  this->instance = instance;

  nrows = 2 * instance.t;
  ncols = 3 * instance.t;

  M = 0;
  for (int i = 0; i < instance.t; i++)
    M += instance.d[i];
  M *= 2;
    
  columns = std::vector<std::vector<double> >
    (ncols, std::vector<double>(nrows, 0.0));

  cost = std::vector<double>(ncols, 0.0);
    
  for (int i = 0; i < ncols; i++) {
    if (i < instance.t) {
      columns[i][i] = 1;
      columns[i][instance.t + i] = 1;
    } 
    else if (instance.t <= i && i < 2*instance.t) {
      columns[i][i-instance.t] = -1;
      if (i < 2*instance.t - 1)
	columns[i][i-instance.t+1] = 1;
    } 
    else if (i >= 2*instance.t)
      columns[i][i - instance.t] = -M;
  }

  // configura funcao objetivo
  for (int i = 0; i < instance.t; i++)
    cost[i] = instance.p[i];
  for (int i = instance.t; i < 2*instance.t; i++)
    cost[i] = instance.h[i-instance.t];
  for (int i = 2*instance.t; i < ncols; i++)
    cost[i] = instance.f[i-2*instance.t];
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
    (*rowtype)[i] = 'E';
  for (int i = instance.t; i < nrows; i++)
    (*rowtype)[i] = 'L';

  // configura limites das variáveis
  for (int i = 0; i < 2*instance.t; i++) {
    (*lb)[i] = 0.0;
    (*ub)[i] = XPRS_PLUSINFINITY;
  }
  for (int i = 2*instance.t; i < ncols; i++) {
    (*lb)[i] = 0;
    (*ub)[i] = 1;
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
    (*colbeg)[i] = off;
    for (int j = 0; j < nrows; j++)
      if (NONZERO(columns[i][j])) {
	(*matval)[off] = columns[i][j];
	(*rowidx)[off++] = j;
      }
  }
  (*colbeg)[ncols] = off;

  // configura variaveis binárias
  nmip = instance.t;
  *miptype = (char *) calloc(nmip, sizeof(char));
  *mipcol = (int *) calloc(nmip, sizeof(int));
  for (int i = 0; i < nmip; i++) {
    (*miptype)[i] = 'B';
    (*mipcol)[i] = 2*instance.t + i;
  }

  // configura rhs
  for (int i = 0; i < instance.t; i++)
    (*rhs)[i] = instance.d[i];
  for (int i = instance.t; i < nrows; i++)
    (*rhs)[i] = 0.0;

  // configura funcao objetivo
  for (int i = 0; i < ncols; i++)
    (*obj)[i] = cost[i];
}

double ULS::solve(std::vector <double>& piAk, std::vector <double>& sol) {
  /* Implementa o algoritmo de Wagher e Within. Esse algoritmo segue
     dois princípios: 
     1) Só se produz quando não há estoque (caso contrário, haveria
     custos de estoque inúteis)
     2) Sempre se produz para atender TODA a demanda de um ou mais
     períodos
     
     Definimos como C[t1,t2] o menor custo para atender a demanda
     entre t1 e t2 exclusive. O objetivo do algoritmo é descobrir o valor de
     C[0,T]. C[t1,t2] pode ser trocado por C[t1,k] + C [k,t2] t1 < k
     < t2, se for essa opção mais barata. */
  
  std::vector<std::vector<double> > C(instance.t);
  std::vector<int> pai(instance.t);
  /* Calcula, inicialmente, o custo de produção a cada instante, para
     atender um ou mais (podendo ser todos) os instantes posteriores */
  for(int i=0;i<instance.t;i++) {
    int tdem = 0;
    int stock = 0;
    int acc = 0;
    pai[i] = i;
    C[i] = std::vector<double>(instance.t);
    for(int j=i;j<instance.t;j++) {
      tdem += instance.d[j];
      acc += instance.d[j]*stock;
      C[i][j] = (instance.f[i]-piAk[2*instance.t+i])*(tdem != 0) + (instance.p[i]-piAk[i])*tdem + acc;
      stock += instance.h[j] - piAk[instance.t + j];
    }
  }
    
  /* Verifica se há partições mais baratas para produzir */

  std::vector<double> melhor(instance.t+1);
  melhor[0] = 0;
  for(int f=1;f<=instance.t;f++) {
    melhor[f] = C[0][f-1] + melhor[0];
    for(int i=1;i<f;i++) {
      melhor[f] = std::min(melhor[f],C[i][f-1] + melhor[i]);
    }
  }

  // for(int m=1;m<instance.t;m++) {
  //   for(int i=0;i<instance.t && i < m;i++) {
  //     for(int j=m+1;j<=instance.t;j++) {
  // 	if (C[i][j] > C[i][m] + C[m][j]) {
  // 	  C[i][j] = C[i][m] + C[m][j];
  // 	  fprintf(stderr,"Produzir em %d e %d pra cumprir demanda até %d é melhor que só em %d\n",
  // 		  i,m,j,i);
  // 	}
  //     }
  //   }
  // }

  sol = std::vector<double>(3*instance.t);
  /* Constrói vetor {x,s,y} */
  int next = instance.t-1;
  for(int f=instance.t;f>0;f=next) {

    next = 0;
    for(int j=f-1;j>=0;j--) {
      if (C[j][f-1] == melhor[f]) {
	/* produziu em j a demanda de todo mundo entre j e f-1 */
	sol[2*instance.t + j] = 1;
	next = j;
	break;
      }
    }
  }
  for(int i = 0;i<instance.t;i++) {
    if (sol[2*instance.t + i] > 0.3)
      ;
      // fprintf(stderr,"Produziu em %d\n",i);
  }

  // fprintf(stderr,"Opt: %6.1lf\n",melhor[instance.t]);
  return (melhor[instance.t]);
}
