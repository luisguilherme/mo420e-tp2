#include "milsp_instance.H"


MILSPInstance::MILSPInstance(FILE* fp) {
  char buffer[256];
  fscanf(fp," %[^\n]",buffer);
  nome = std::string(buffer);
  fscanf(fp," %d %d",&t,&m);
  C = vi(t);
  for(int k=0;k<t;k++) {
    int ct, tt;
    fscanf(fp," %d %d",&ct,&tt);
    assert(tt == k);
    C[k] = ct;
  }
  p = vvi(m); h = vvi(m); f = vvi(m); d = vvi(m);
  for(int i=0;i<m;i++) {
    p[i] = vi(t); h[i] = vi(t); f[i] = vi(t); d[i] = vi(t);
    for(int k=0;k<t;k++) {
      int tt, ii;
      int pit,hit,fit,dit;
      fscanf(fp," %d %d %d %d %d %d",&ii,&tt,&pit,&hit,&fit,&dit);
      assert(tt == k && ii == i);
      p[i][k] = pit; h[i][k] = hit; f[i][k] = fit; d[i][k] = dit;
    }
  }
}
