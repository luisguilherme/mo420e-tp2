#ifndef __INTEGER_PROGRAM_H__
#define __INTEGER_PROGRAM_H__

#include "instance.H"

class IntegerProgram { 
  Instance instance;
  int ncols, nrows;
public:
  virtual Instance& getInstance() { return instance; }
  virtual int getncols() { return ncols; }
  virtual int getnrows() { return nrows; }
  virtual void getParam(int&,int&,char**,double**,double**,int**,int**,
			double**,double**,double**,int&,char**,int**,bool&) { }
  
};

#endif
