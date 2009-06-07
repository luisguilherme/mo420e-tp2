#ifndef __ULS_H__
#define __ULS_H__

#include "integer_program.H"
#include "uls_instance.H"
#include <vector>
#include "xprs.h"
#include "definitions.H"

class ULS : public IntegerProgram {
  ULSInstance instance;
  
public:
  virtual void getParam(int&,int&,char**,double**,double**,int**,int**,
			double**,double**,double**,int&,char**,int**,bool&);
  virtual void addcol(std::vector<double>&, double);
  ULS(const ULSInstance& instance);
};

#endif
