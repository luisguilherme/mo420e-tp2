#include "milsp_instance.H"
#include "integer_program.H"
#include "column_generation.H"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#define EPS 1e-9
#define sz(a) (int) (a).size()
#define pb push_back
#define NONZERO(a) ((a) > EPS || (a) < -EPS)

class MILSP : public IntegerProgram {
  std::vector<std::vector<double> > columns; //colunas transpostas;
  std::vector<double> cost; //tamanho dinâmico
  int ncols;
  MILSPInstance instance;
public:
  Instance& getInstance() { return instance; }
  void getParam(int& ncol,int& nrow,char** rowtype,double** rhs,
		double** obj,int** colbeg, int** rowidx,
		double** matval,double** lb,double** ub,
		int& nmip,char** miptype,int** mipcol,bool& relaxed) {

    relaxed = true;
    
    /* Aloca memória dinamicamente para o XPress */
    ncol = ncols;
    nrow = 2*instance.t + instance.m;
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

  /* TODO: Gerar planos de produção que sempre dão uma solução viável.
           + Produzir sempre a capacidade maxima para todo item 
	     em todo período
	   + Com esses planos é sempre possível obter uma solução
	     da Relaxação do Meste Restrito
   */  

  MILSP(MILSPInstance& instance) {
    this->instance = instance;
    
    ncols = instance.m*instance.t;
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
	  }
	  else if (inst == tp) { //produz toda a demanda
	    x[sch].pb(dem);
	    y[sch].pb(1);
	    dem -= instance.d[prod][inst];
	    s[sch].pb(dem);
	  }
	  else {
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
    cost = std::vector<double>(ncols);
    /* Gera as colunas baseado no x */
    for(int prod=0,col=0;prod<instance.m;prod++) {
      std::cerr << prod << std::endl;
      for(int tp=0;tp<instance.t;tp++,col++) {
	int row = 0;
	std::cerr << "-- " << tp << std::endl;
	// columns[col] = std::vector<double>(2*instance.t + instance.m);
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
	  cost[col] += instance.p[prod][inst]*x[col][inst] + instance.f[prod][inst]*y[col][inst] + instance.h[prod][inst]*s[col][inst];
	}
	
	
      } // end for tp
    } // end for prod
      
  } // end MILSP
    

  void OldMILSP(MILSPInstance& instance) {
    this->instance = instance;
    
    ncols = instance.m;
    vvi x(ncols);

    for(int sch=0;sch<instance.m;sch++) {
      for(int inst=0;inst<instance.t;inst++) {
	x[sch].pb(instance.d[sch][inst]);
      }
    }
    
    /* gera as primeiras colunas */
    columns = std::vector<std::vector<double> >(ncols);
    for(int col=0;col<ncols;col++) {
      /* Em todos os instantes gera todos os esquemas */
      for(int inst=0;inst<instance.t;inst++) {
	columns[col].pb(1);
      }
      /* Produz a demanda atual, sempre (JIT) */
      for(int inst=0;inst<instance.t;inst++) {
	columns[col].pb(x[col][inst]);
      }

      /* Se o esquema é do produto ou não (Identidade) */
      for(int sch=0;sch<ncols;sch++) {
	int val = ((col == sch)?1:0);
	columns[col].pb(val);
      }

      /* custo do esquema */ 
      cost.pb(0);
      for(int inst=0;inst<instance.t;inst++) {
	cost[col] += instance.p[col][inst]*x[col][inst] + instance.f[col][inst];
      }
      
    } // fim for col
  }
    
};


int main(int argc, char* argv[]) {
  FILE* fp = fopen(argv[1], "r");
  MILSPInstance mi;
  mi.loadFrom(fp);
  MILSP prob(mi);
  ColumnGeneration cg;
  cg.configureModel(0,prob,prob);
  cg.solveRestricted();
  return(0);
}
