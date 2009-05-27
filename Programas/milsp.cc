#include "milsp_instance.H"
#include "integer_program.H"
#include <cstdio>

class MILSP : public IntegerProgram {
  MILSPInstance instance;
  Instance& getInstance() { return instance; }
  void getParam(int& ncol,int& nrow,char* rowtype,double* rhs,
			double* obj,int* colbeg, int* rowidx,
			double* matval,double* lb,double* ub,
			int& nmip,char* miptype,int* mipcol,bool& relaxed) {
    
  }

};

int main(int argc, char* argv[]) {
  FILE* fp = fopen(argv[1],"r");
  MILSPInstance mi;
  mi.loadFrom(fp);
  return(0);
}
