#ifndef __INTEGER_PROGRAM_H__
#define __INTEGER_PROGRAM_H__

#include "instance.H"
#include <vector>

class IntegerProgram { 
  Instance instance;
protected:
  std::vector<std::vector<double> > columns;
  std::vector<double> cost;
  int ncols, nrows;
public:
  virtual Instance& getInstance() { return instance; }
  virtual int getncols() { return ncols; }
  virtual std::vector<std::vector<double> >& getcolumns() { return columns; }
  virtual std::vector<double>& getcost() { return cost; }
  virtual int getnrows() { return nrows; }
  virtual void getParam(int&,int&,char**,double**,double**,int**,int**,
			double**,double**,double**,int&,char**,int**,bool&) { }
  virtual void addcol(std::vector<double>&, double) {}
};

#endif
