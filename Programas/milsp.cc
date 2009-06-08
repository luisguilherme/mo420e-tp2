#include "milsp_instance.H"
#include "uls_instance.H"
#include "integer_program.H"
#include "column_generation.H"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include "definitions.H"
#include "uls.H"


class MILSP : public IntegerProgram {
  MILSPInstance instance;
public:

  void getParam(int& ncol,int& nrow,char** rowtype,double** rhs,
		double** obj,int** colbeg, int** rowidx,
		double** matval,double** lb,double** ub,
		int& nmip,char** miptype,int** mipcol,bool& relaxed) {

    relaxed = true;

    /* Aloca memória dinamicamente para o XPress */
    ncol = ncols;
    nrow = nrows;
    nmip = 0; miptype = NULL; mipcol = NULL;
    *rowtype = (char *) calloc(nrow,sizeof(char));
    *rhs = (double *) calloc(nrow,sizeof(double));
    *obj = (double *) calloc(ncol,sizeof(double));
    *colbeg = (int *) calloc(ncol+1,sizeof(int));
    *lb = (double *) calloc(ncol,sizeof(double));
    *ub = (double *) calloc(ncol,sizeof(double));

    /* Gerando tipo de cada restrição e RHS */
    int nnzero = 0;
    for(int i=0;i<instance.t;i++) {
      (*rowtype)[i] = 'L';
      (*rhs)[i] = 1;
    }

    for(int pos = instance.t, i=0;i<instance.t;i++,pos++) {
      (*rowtype)[pos] = 'L';
      (*rhs)[pos] = instance.C[i];
    }

    for(int pos = 2*instance.t, i = 0;i<instance.m;pos++,i++) {
      (*rowtype)[pos] = 'E';
      (*rhs)[pos] = 1;
    }

    int matsize = 0;
    for(int col=0;col<ncols;col++)
      for(int row=0;row<sz(columns[col]);row++)
	matsize += NONZERO(columns[col][row]);

    /* aloca matriz esparsa */
    *matval = (double *) calloc(matsize,sizeof(double));
    *rowidx = (int *) calloc(matsize,sizeof(int));

    /* Transformando colunas em dados de entrada do XPress */
    for(int col=0;col<ncols;col++) {
      (*obj)[col] = cost[col];
      (*colbeg)[col] = nnzero;
      (*lb)[col] = 0; (*ub)[col] = 1;
      for(int row=0;row<sz(columns[col]);row++) {
	if (NONZERO(columns[col][row])) {
	  (*matval)[nnzero] = columns[col][row];
	  (*rowidx)[nnzero++] = row;
	}
      }
    }
    (*colbeg)[ncols] = nnzero;

#ifdef PRINTDEBUG
    for(int i=0;i<nrow;i++) {
      for(int j=0;j<ncol;j++) {
	printf("%6.1lfl%d %c",columns[j][i],j,((j==ncol-1)?'=':'+'));
      }
      printf("%.1lf\n",(*rhs)[i]);
    }
#endif

  }


  MILSP(MILSPInstance& instance) {
    this->instance = instance;
 
    nrows = 2*instance.t + instance.m;
    ncols = instance.m * instance.t;
 
    vvi x(ncols);
    vvi y(ncols);
    vvi s(ncols);
 
    vi demanda;
    for(int prod=0;prod<instance.m;prod++) {
      demanda.pb(0);
      for(int inst=0;inst<instance.t;inst++) {
	demanda[prod] += instance.d[prod][inst];
      }
    }
 
    for(int prod=0,sch=0;prod<instance.m;prod++) {
      for(int tp=0;tp<instance.t;tp++,sch++) {
	int dem = demanda[prod];
	for(int inst=0;inst<instance.t;inst++) {
	  if (inst < tp) {
	    int toprod = instance.d[prod][inst];
	    dem -= toprod;
	    x[sch].pb(toprod);
	    y[sch].pb(1);
	    s[sch].pb(0);
	  } else if (inst == tp) { //produz toda a demanda
	    x[sch].pb(dem);
	    y[sch].pb(1);
	    dem -= instance.d[prod][inst];
	    s[sch].pb(dem);
	  } else {
	    x[sch].pb(0);
	    y[sch].pb(0);
	    dem -= instance.d[prod][inst];
	    s[sch].pb(dem);
	  } // end else
	} // end for inst
      } // end for tp
    } // end for prod
 
    columns = std::vector<std::vector<double> >(ncols);
    for(int col=0;col<ncols;col++) {
      columns[col] = std::vector<double>(2*instance.t + instance.m);
    }
 
    // custo "infinito" as colunas da identidade
    double inf = 0.0;
 
    cost = std::vector<double>(ncols);
    /* Gera as colunas baseado no x */
    for(int prod=0,col=0;prod<instance.m;prod++) {
      for(int tp=0;tp<instance.t;tp++,col++) {
	int row = 0;
 
	for(int inst=0;inst<instance.t;inst++,row++) {
	  columns[col][row] = y[col][inst];
	}
 
	for(int inst=0;inst<instance.t;inst++,row++) {
	  columns[col][row] = x[col][inst];
	}
 
	for(int k=0;k<instance.m;k++,row++) {
	  columns[col][row] = ((k == prod)?1:0);
	}
 
	cost[col] = 0;
	for(int inst=0;inst<instance.t;inst++) {
	  cost[col] += instance.p[prod][inst]*x[col][inst] +
	    instance.f[prod][inst]*y[col][inst] +
	    instance.h[prod][inst]*s[col][inst];
	}
	inf += cost[col];
 
      } // end for tp
    } // end for prod
 
    /* "Cola" uma matriz identidade para garantir factibilidade
       do problema restrito */
    int mid = 2*instance.t + instance.m;
    for (int i = 0; i < mid; i++) {
      columns.push_back(std::vector<double>(mid, 0.0));
      columns[ncols++][i] = 1.0;
      cost.push_back(inf);
    }
 
  } // end MILSP

  void addcol(std::vector<double> &sol, int np, double c,int& nz,double** custo,
	      int** colbeg, int** rowidx, double** matval, double** lb,
	      double** ub) {
    *custo = (double *) malloc(sizeof(double));
    *lb = (double *) malloc(sizeof(double));
    *ub = (double *) malloc(sizeof(double));
    *colbeg = (int *) calloc(2,sizeof(int));
    
    int tt = instance.t;
    int mm = instance.m;

    std::vector<double> col(2*tt + mm);
    for(int i=0;i<tt;i++) {
      col[i] = sol[2*tt + i];
      col[i+tt] = sol[i];
    }
    col[2*tt + np] = 1;


    columns.pb(col);
    cost.pb(c);
    ncols++;

    **custo = c;
    **lb = 0.0;
    **ub = 1.0;
    (*colbeg)[0] = 0; (*colbeg)[1] = 1;

    nz = 0;
    for (int i = 0; i < nrows; i++) if (NONZERO(col[i])) nz++;
    
    *matval = (double *) calloc(nz,sizeof(double));
    *rowidx = (int *) calloc(nz,sizeof(int));

    int mpos = 0;
        
    for(int i=0;i<nrows;i++) {
      if (NONZERO(col[i])) {
	(*rowidx)[mpos] = i;
	(*matval)[mpos++] = col[i];
      }
    }
    
  }

  std::vector<double> piProductAk(double *dual) {
    std::vector<double> piAk(3*instance.t, 0.0);

    for (int i = 0; i < instance.t; i++)
      piAk[i] = dual[i+instance.t];

    for (int i = instance.t; i < 2*instance.t; i++)
      piAk[i] = 0.0;

    for (int i = 2*instance.t; i < 3*instance.t; i++)
      piAk[i] = dual[i-2*instance.t];

    return piAk;
  }

};

int main(int argc, char* argv[]) {
  FILE* fp = fopen(argv[1], "r");
  MILSPInstance mi;

  mi.loadFrom(fp);

  MILSP prob(mi);

  std::vector<IntegerProgram*> pricing;

  // gera subproblemas de pricing para MILSP (um ULS para cada produto)
  char prefixo[16];
  for (int i = 0; i < mi.m; i++) {
    sprintf(prefixo, "prod%d", i);
    std::string nome(prefixo);
    pricing.push_back(new ULS(ULSInstance(nome,
					  mi.t,
					  mi.d[i],
					  mi.f[i],
					  mi.h[i],
					  mi.p[i])));
  }

  ColumnGeneration cg(prob, pricing);

  cg.configureModel(0, prob, pricing);
  bool goon = true;
  int ncolumns = 0;
  time_t start = time(&start);
  for (int i = 0;goon && ncolumns < 2000; i++) {
    cg.solveRestricted();
    int gencolums = cg.solvePricing();
    goon = goon && (gencolums > 0);
    ncolumns += gencolums;
    time_t now = time(&now);
    if (now - start > 20*60) break;
    // getchar();
  }

  if (cg.totalInteiros > 0) {
    printf("## Melhor solução inteira: %d\n",(int) cg.melhorPrimal);
  }
  else {
    printf("## Sem soluções inteiras\n");
  }
  
  return 0;
}
