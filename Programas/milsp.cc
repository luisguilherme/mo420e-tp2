#include "milsp_instance.H"
#include <cstdio>

int main(int argc, char* argv[]) {
  FILE* fp = fopen(argv[1],"r");
  MILSPInstance mi(fp);
  return(0);
}
