#include "milsp_instance.H"
#include "integer_program.H"
#include "column_generation.H"
#include <cstdio>
#include <cstdlib>
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
  void getParam(int& ncol,int& nrow,char* rowtype,double* rhs,
		double* obj,int* colbeg, int* rowidx,
		double* matval,double* lb,double* ub,
		int& nmip,char* miptype,int* mipcol,bool& relaxed) {

    relaxed = true;
    
    /* Aloca memória dinamicamente para o XPress */
    ncol = ncols;
    nrow = 2*instance.t + instance.m;
    nmip = 0; miptype = NULL; mipcol = NULL;
    rowtype = (char *) calloc(nrow,sizeof(char));
    rhs = (double *) calloc(nrow,sizeof(double));
    obj = (double *) calloc(ncol,sizeof(double));
    colbeg = (int *) calloc(ncol,sizeof(int));
    lb = (double *) calloc(ncol,sizeof(double));
    ub = (double *) calloc(ncol,sizeof(double));

    /* Gerando tipo de cada restrição e RHS */
    int nnzero = 0;
    for(int i=0;i<instance.t;i++) {
      rowtype[i] = 'L';
      rhs[i] = 1;
    }

    for(int pos = instance.t, i=0;i<instance.t;i++,pos++) {
      rowtype[pos] = 'L';
      rhs[pos] = instance.C[i];
    }
    
    for(int pos = 2*instance.t, i = 0;i<instance.m;pos++,i++) {
      rowtype[pos] = 'E';
      rhs[pos] = 1;
    }
    
    int matsize = 0;
    for(int col=0;col<ncols;col++) 
      for(int row=0;row<sz(columns[col]);row++) 
	matsize += NONZERO(columns[col][row]); 
    
    /* aloca matriz esparsa */
    matval = (double *) calloc(matsize,sizeof(double));
    rowidx = (int *) calloc(matsize,sizeof(int));

    /* Transformando colunas em dados de entrada do XPress */
    for(int col=0;col<ncols;col++) {
      obj[col] = cost[col];
      colbeg[col] = nnzero;
      lb[col] = 0; ub[col] = 1;
      for(int row=0;row<sz(columns[col]);row++) {
	if (NONZERO(columns[col][row])) {
	  matval[nnzero] = columns[col][row];
	  rowidx[nnzero++] = row;
	}
      }
    }
  }
  

  MILSP(MILSPInstance& instance) {
    this->instance = instance;
    
    ncols = instance.m;
    vvi x(ncols);
    /* as primeiras colunas a serem geradas são esquemas de produção
       just-in-time, um para cada item, que é viável linearmente */
    
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
      cost[col] = 0;
      for(int inst=0;inst<instance.t;inst++) {
	cost[col] += instance.p[col][inst]*x[col][inst] + instance.f[col][inst];
      }
      
    } // fim for col
  }
    
};


int main(int argc, char* argv[]) {
  FILE* fp = fopen(argv[1],"r");
  MILSPInstance mi;
  mi.loadFrom(fp);
  MILSP prob(mi);
  ColumnGeneration cg;
  cg.configureModel(0,prob,prob);
  cg.solveRestricted();
  return(0);
}
