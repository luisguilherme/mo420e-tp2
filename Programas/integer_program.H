#ifndef __INTEGER_PROGRAM_H__
#define __INTEGER_PROGRAM_H__
#include "instance.H"
class IntegerProgram { 
public:
  virtual Instance& getInstance();
  virtual void getParam(int&,int&,char*,double*,double*,int*,int*,
			double*,double*,double*,int&,char*,int*,bool&);

};

#endif
