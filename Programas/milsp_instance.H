#ifndef __MILSP_INSTANCE_H__
#define __MILSP_INSTANCE_H__

#include <cstdio>
#include <vector>
#include <string>
#include <cassert>
#include "instance.H"

typedef std::vector<int> vi;
typedef std::vector<vi> vvi;

class MILSPInstance : public Instance {
  
public:

  std::string nome;
  
  int t; //total de tempo
  int m; //numero de itens

  vvi d; //demanda
  
  vi C; //capacidade de produção
  
  vvi f; //custo de setup
  vvi h; //custo de armazenagem
  vvi p; //custo de produção

  //MILSPInstance(string filename);
  //MILSPInstance(FILE* fp);
  virtual void loadFrom(FILE* fp);
  //MILSPInstance(vvi& d, vi& C, vvi& p, vvi& h, vvi& f);
  
};

#endif
