#ifndef __INSTANCE_H__
#define __INSTANCE_H__

#include <cstdio>

class Instance {
public:
  virtual double* getNewObj(double* dual, int* mindex);
  virtual void loadFrom(FILE* fp);
  virtual void getParam(int&,int&,char*,double*,double*,int*,int*,
			double*,double*,double*,int&,char*,int*,bool&);
  virtual int dim(void);
};

#endif
