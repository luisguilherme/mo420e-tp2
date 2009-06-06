#ifndef __INTEGER_PROGRAM_H__
#define __INTEGER_PROGRAM_H__

#include "instance.H"

class IntegerProgram { 
  Instance instance;
public:
  virtual Instance& getInstance() { return instance; }
  virtual int getncols() { return 0; }
  virtual int getnrows() { return 0; }
  virtual void getParam(int&,int&,char**,double**,double**,int**,int**,
			double**,double**,double**,int&,char**,int**,bool&) { }
  
};

#endif
