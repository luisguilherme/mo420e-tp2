#ifndef __ULS_INSTANCE_H__
#define __ULS_INSTANCE_H__

#include <cstdio>
#include <vector>
#include <string>
#include <cassert>
#include "instance.H"

typedef std::vector<int> vi;

class ULSInstance : public Instance {

public:

  std::string nome;
  
  int t; //total de tempo

  vi d; //demanda
  vi f; //custo de setup
  vi h; //custo de armazenagem
  vi p; //custo de produção

  ULSInstance() {}
  ULSInstance(std::string& nome, int t, vi& d, vi& f, vi& h, vi& p);
  
  virtual int dim(void) { return 0; }
};

#endif
