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
  void getParam(int&,int&,char**,double**,double**,int**,int**,
			double**,double**,double**,int&,char**,int**,bool&);
  double solve(std::vector<double>& sol);
  ULS(const ULSInstance& instance);
};

#endif
