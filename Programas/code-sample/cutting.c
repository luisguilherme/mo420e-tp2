/* ******************************************************************
   Algoritmo de geracao de colunas para o problema de Corte Unidimensional

   Autor: Edna Ayako Hoshino (eah@dct.ufms.br)
   baseado no exemplo do prof Cid Carvalho de Souza
   Data: 10/2005 
   Arquivo: cutting.c 

Versoes:

   - ainda NAO usa limitantes duais para a relaxacao linear do problema mestre;
   - termina quando nao existem colunas com custo reduzido
   - nao faz uso de heuristicas primais;
   - nao faz pricing parcial (ou seja, resolve de forma exata o subproblema de pricing)

   Entrada: nome do arquivo de instancia

   Formato da instancia: 
   
    m L  (m= numero de itens de corte, L=comprimento da barra)
    c_i d_i (comprimento e demanda do item i)
	
    Formulacao construida (ver cutting.ps)
**************************************************************** */
 
/* includes dos arquivos .h */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "cutting.h" 

/* ============================ */
/* Inicio do programa principal */
/* ============================ */
int main(int argc, char * argv[]){
  /* variaveis auxiliares */  
  int i,
    formato;  /* formato de impressao do otimizador na saida padrao */
  clock_t t1, t2;

  /* var para impressao de resultados */
  int totalInteiros, itmelhorPrimal, it;
  double tempo, tempoLP, tempoPricing, tPMR;

  /* var para gerar coluna */
  int *mindex, colunas, totalColunas;
  double *dual, *obj_pricing, z_PMR, melhorPrimal;

  /*   usadas para xpress */
  XPRSprob probMestre, probPricing; /* estrutura do XPRESS contendo o problema */
  int xpress_ret;			/* valor de retorno das rotinas do XPRESS */
  int xpress_status;		/* status retornado por algumas rotinas do XPRESS */
	
  /* dados da instancia */
  cuttingInstance cutting;

  /* dados da solucao */
  solution sol;

  /* ========================================= */
  /* Verifica argumentos da linha de comando */
  /* ========================================= */
	
  if (argc<3){ 
    showUsage();
    exit(-1);
  }
  /* Carrega os dados da instancia */
  if (!CargaInstancia(argv[1], &cutting)){
    printf("Erro na abertura do arquivo (%s)", argv[1]);
    exit(1);
  }
  formato = atoi(argv[2]); 
  if (formato<0 || formato>3){
    printf("Formato invalido: escolha de 0 a 3");
    exit(1);
  }

  /* ========================================= */
  /* inicializacão do XPRESS                   */
  /* ========================================= */
  xpress_ret=XPRSinit("");
  if (xpress_ret) 
    errormsg("Main: Erro de inicializacao do XPRESS.\n",__LINE__,xpress_ret, probMestre);

  /* inicializa valores de variaveis globais */
  tPMR=0;

  printf("\n==========================\nOtimizacao do LP e IP");
	    
  /* Cria, carrega e configura cada modelo */
  CriaCarregaConfigura(formato, &probMestre, &probPricing, &cutting, &sol);

  /*=========================================================================*/
  /* resolve PMR, iterativamente
  /*=========================================================================*/

  /* inicializa var */
  it=0;
  itmelhorPrimal=-1;
  tempoLP=0.0;
  tempoPricing=0.0;
  melhorPrimal=XPRS_PLUSINFINITY;
  totalInteiros=0;
	    
  /* var para resolver pricing */
  dual=(double *)malloc(cutting.m*sizeof(double));
  obj_pricing=(double *)malloc(cutting.m*sizeof(double));
  mindex=(int *)malloc(cutting.m*sizeof(int));

  /* inicializa solucao */
  sol.xstar=(double *)malloc(cutting.m*sizeof(double));
  sol.zstar=XPRS_MINUSINFINITY; 
	    
  do {
    /*========================================================================= */
    /* resolve o problema mestre reduzido                                       */
    /*========================================================================= */
    printf("\n==========================\nOtimizacao do PMR \n");
	    
    /* inicializa o pool de colunas */
    sol.totPool=0;

    XPRSinitglobal(probMestre); /* reinicializa o prob mestre reduzido*/
	    
    /* resolve o PMR com o optimizer */
    t1=clock();
    xpress_ret=XPRSminim(probMestre,NULL); /* g = algoritmo de busca B&B, NULL = PL */
    t2=clock();

    tempo=((double)(t2-t1))/CLOCKS_PER_SEC;
    tempoLP+=t2-t1;

    if (xpress_ret) 
      errormsg("Main: Erro na chamada da rotina XPRSmaxim.\n",
	       __LINE__,xpress_ret, probMestre);

    /* obtem valor das var duais e valor da solucao */
    if (IsIntegerSol(probMestre, &z_PMR, dual)){ 
      /* se solucao eh inteira --> atualiza bound primal, se melhor */
      totalInteiros++;
      if (z_PMR < melhorPrimal){
	melhorPrimal = z_PMR;
	itmelhorPrimal = it;
      }
    }
	    
    printf("\nIteracao: %d\n", it );
    printf("- Valor da solucao otima do PMR =%12.6f \n",z_PMR);
    printf("- Tempo do LP = %lf\n", tempo);

    /*========================================================================= */
    /* resolve o problema de pricing                                            */
    /*========================================================================= */
    printf("\n==========================\nOtimizacao do Pricing \n");
	    
    /* calcula os novos coeficientes das variaveis na 
       funcao objetivo do subproblema de pricing */
    for(i=0;i<cutting.m;i++){
      obj_pricing[i]=dual[i];
      mindex[i]=i;
      printf("\ndual[%d]=%lf", i, dual[i]);
    }
	    
    /* resolve o sub de pricing */
    /* uma vez maximizado, perde a base inicial, 
       e o problema tem que ser carregado de novo? */
    XPRSinitglobal(probPricing); /* reinicializa o prob de pricing */

    /* seta novos coeficientes na funcao objetivo */
    XPRSchgobj(probPricing,cutting.m,mindex,obj_pricing); 

    /* reinicializa a melhor solucao */
    sol.zstar=XPRS_MINUSINFINITY; 
    xpress_ret=XPRSsetdblcontrol(probPricing,XPRS_MIPABSCUTOFF,sol.zstar);
	    
    t1=clock();
    xpress_ret=XPRSmaxim(probPricing,"g"); /* g = algoritmo de busca B&B
					      NULL = PL */
    t2=clock();
	      
    tempo=((double)(t2-t1))/CLOCKS_PER_SEC;
    tempoPricing+=t2-t1;
	    
    /* imprime resultado do pricing */
    printf("\n==========================\nResultado do Pricing");
    printf("\nZstar (custo reduzido): %lf", 1-sol.zstar);
    printf("\nTempo do pricing: %lf\n", tempo);
    fflush(stdout);
    getchar();  

    // custo reduzido negativo
    if (1-sol.zstar<-EPSILON){
      /*=========================================================== */
      /* adiciona as colunas do pool                                */
      /*=========================================================== */
      adicionaCols(probMestre, sol.pool, sol.totPool, &cutting);
      ImprimeCols(sol.pool,sol.totPool, &cutting);
      totalColunas+=sol.totPool;
    } else
      break; /* nenhuma coluna com custo reduzido negativo =>
		solucao otima do problema mestre relaxado foi encontrada */

    /* libera mem do pool */
    for(i=0;i<sol.totPool;i++){
      if(sol.pool[i])
	free(sol.pool[i]);
    }
    it++;
  } while(1); // repete ate convergir

  tPMR=(tempoLP+tempoPricing)/CLOCKS_PER_SEC;

  /*========================================================================= */
  /* impressao final                                                          */
  /*========================================================================= */
  printf("\n==========================\nResultado do PMR\n");
  printf(".total de iteracoes........................... = %d\n",it);
  printf(".total de solucoes inteiras geradas........... = %d\n",totalInteiros);
  printf(".valor da melhor solucao inteira ............. = %lf\n",melhorPrimal); 
  printf(".it da melhor solucao inteira................. = %d\n",itmelhorPrimal);
  printf(".tempo total de LP............................ = %lf\n",tempoLP/CLOCKS_PER_SEC);
  printf(".tempo total de Pricing....................... = %lf\n",tempoPricing/CLOCKS_PER_SEC);
  printf(".valor do PMR................................. = %lf\n",z_PMR); 
  fflush(stdout);
	
  /* libera toda memoria usada no programa */
  if (dual) free(dual);
  if (obj_pricing) free(obj_pricing);
  if (mindex) free(mindex);
  if (sol.xstar) free(sol.xstar);
  XPRSdestroyprob(probMestre);
  XPRSdestroyprob(probPricing);
	
  /* finaliza */
  xpress_ret=XPRSfree();
  if (xpress_ret)
    errormsg("Main: Erro na liberacao de memoria final.\n",__LINE__,xpress_ret, probMestre);
  return 0;
}

/* =====================================================================
 * Rotina (callback) para salvar a  melhor solucao.  Roda para todo nó
 * onde acha solucao inteira.
 * callback para XPRSsetcbintsol
 * Autor: Cid Carvalho de Souza
 * =====================================================================
*/
void XPRS_CC SalvaMelhorSol(XPRSprob prob, void *psol)
{
  int i, cols, peso_aux=0, node, xpress_ret;
  double objval, *x;
  Boolean viavel;
  solution* sol;

  sol=(solution*) psol;

  /* pega o numero do nó corrente */
  xpress_ret=XPRSgetintattrib(prob,XPRS_NODES,&node);
  if (xpress_ret) 
    errormsg("SalvaMelhorSol: rotina XPRSgetintattrib.\n",__LINE__,xpress_ret, prob);

  /* recupera numero de colunas do prob */
  xpress_ret=XPRSgetintattrib(prob,XPRS_COLS,&cols);
  if (xpress_ret) 
    errormsg("SalvaMelhorSol: rotina XPRSgetintattrib\n",__LINE__,xpress_ret, prob);

  /* recupera solucao */
  x=(double *)malloc(cols*sizeof(double));

  xpress_ret=XPRSgetdblattrib(prob,XPRS_LPOBJVAL,&objval);
  if (xpress_ret) 
    errormsg("SalvaMelhorSol: rotina XPRSgetdblattrib\n",__LINE__,xpress_ret, prob);

  xpress_ret=XPRSgetsol(prob,x,NULL,NULL,NULL);
  if (xpress_ret) 
    errormsg("SalvaMelhorSol: Erro na chamada da rotina XPRSgetsol\
			\n",__LINE__,xpress_ret, prob);

  /* guarda no pool de colunas se o custo reduzido da solucao inteira é negativo */
  if (1-objval<0 && sol->totPool<MAX_POOL){
    (sol->pool)[sol->totPool] = (double*) malloc(cols*sizeof(double));
    for(i=0;i<cols;i++) 
      (sol->pool)[sol->totPool][i]=x[i];
    (sol->totPool)++;
  }

  /* testa se a solução é viável */
  viavel = true;

  /* se a solucao tiver custo melhor que a melhor solucao disponivel entao salva */
  if ((objval > sol->zstar) && viavel ) { 
    for(i=0;i<cols;i++) 
      sol->xstar[i]=x[i];
    sol->zstar=objval;
     
    /* informa xpress sobre novo incumbent */
    xpress_ret=XPRSsetdblcontrol(prob,XPRS_MIPABSCUTOFF,sol->zstar);//zstar+1.0-EPSILON);
    if (xpress_ret) 
      errormsg("SalvaMelhorSol: XPRSsetdblcontrol.\n",__LINE__,xpress_ret, prob);
  }   
  free(x);
  return;
}
/* ================================= */
/* carga do subproblema de pricing   */ 
/* ================================= */
int CargaPricing(XPRSprob *probPricing, cuttingInstance *cutting)
{
  int xpress_ret, i, j, z, nz;  /* nz = contador de entradas nao-nulas */

  /* variaveis para carga do LP */
  int ncols, nrows; /* numero de linhas e colunas da matriz A de restricoes do LP*/
  double *obj; /* coeficientes das variaveis na funcao objetivo */

  /* coeficientes nao-nulos de A: */
  int 
    *mrwind, /* 1 entrada para cada a_{ij} nao-nulo, representa i, i.e., o indice da linha */
    *mstart, /* 1 entrada para cada coluna j, mstart[j]= offset do primeiro a_{ij} nao-nulo em mrwind */
    *mnel;   /* 1 entrada para cada coluna j, representa o total de coeficientes nao-nulos na coluna j*/
  double *dmatval; /* 1 entrada para cada a_{ìj} nao-nulo, representa o valor de a_{ij} */

  /* rhs e tipo de cada restricao */
  double  *rhs; /* 1 entrada para cada restricao, representa o rhs da restricao */
  char *qrtype; /* 1 "         "   "     "      , pode ser l, g, e respectivamente para <=, >= e = */

  /* tipo especial de cada variavel */
  double *dlb, *dub; /* limites inferior e superior */
  int *mgcols;
  char *qgtype; /* indice da variavel e seu tipo */

  ncols=cutting->m; /* total de variaveis */
  nrows=1; /* total de restricoes */
	
  /* aloca espaco para os arrays que serão usados pelo XPRESS na carga do problema.*/
  obj=(double *)malloc(ncols*sizeof(double));
  mrwind=(int *)malloc(cutting->m*sizeof(int));/* uma para cada entrada nao-nula em A */
  dmatval=(double *)malloc(cutting->m*sizeof(double)); /* idem */
  mstart=(int *)malloc(ncols*sizeof(int)); /* uma entrada para cada coluna*/
  mnel=(int *)malloc(ncols*sizeof(int));
  dlb=(double *)malloc(ncols*sizeof(double));
  dub=(double *)malloc(ncols*sizeof(double));
  qgtype = malloc(ncols*sizeof(char));
  mgcols=(int *)malloc(ncols*sizeof(int));

  /* arrays com informações sobre as restrições */
  qrtype = malloc(nrows*sizeof(char));
  rhs=(double *)malloc(nrows*sizeof(double));

  for(i=0;i<ncols;i++){
    obj[i]=0.0;
    dlb[i]=0.0;    /* limite inferior das variáveis */
    dub[i]=XPRS_PLUSINFINITY;    /* limite superior das variáveis */
  }
  for(j=0;j<ncols;j++){
    mgcols[j]=j;   /* todas variáveis xj são binarias  */         
    qgtype[j]='I'; /* todas variáveis são inteiras */       
  }
  /* configuracao das restricoes */
  rhs[0]=cutting->L;
  qrtype[0]='L'; 

  /* carrega os coeficientes das variaveis y  */
  nz=0;
  for(j=0;j<ncols;j++){
    mnel[j]=1; 
    mstart[j]=nz;   /* indice em dmatval da primeira entrada nao-nula da coluna j */ 
    mrwind[nz]=0;
    dmatval[nz++]=(cutting->w)[j];
  }

  /* carga do modelo*/
  xpress_ret=XPRSloadglobal(*probPricing, "pricing", ncols, nrows,  qrtype, rhs, NULL, obj, \
			    mstart, mnel, mrwind, dmatval, dlb, dub, ncols, 0, qgtype, mgcols, \
			    NULL, NULL, NULL,  NULL, NULL);

  /* libera memoria dos vetores usado na carga do LP */
  free(qrtype); free(rhs); free(obj); free(mstart); free(mrwind); free(dmatval); 
  free(dlb);    free(dub); free(qgtype); free(mgcols); 
  free(mnel);

  return xpress_ret;
}
/* ================================= */
/* carga do problema mestre reduzido */ 
/* usa: n,m,A,c                      */
/* ================================= */
int CargaProblemaMestre(XPRSprob *probMestre, cuttingInstance *cutting)
{
  int xpress_ret, i, j, nz;  /* nz = contador de entradas nao-nulas */

  /* variaveis para carga do LP */
  int ncols, nrows; /* numero de linhas e colunas da matriz A de restricoes do LP*/
  double *obj; /* coeficientes das variaveis na funcao objetivo */

  /* coeficientes nao-nulos de A: */
  int 
    *mrwind, /* 1 entrada para cada a_{ij} nao-nulo, armazena o indice i da linha */
    *mstart, /* 1 entrada para cada coluna j, mstart[j]= offset do primeiro a_{ij} nao-nulo em mrwind */
    *mnel;   /* 1 entrada para cada coluna j, representa o total de coeficientes nao-nulos em A na coluna j*/
  double *dmatval; /* 1 entrada para cada a_{ìj} nao-nulo, representa o valor de a_{ij} */

  /* rhs e tipo de cada restricao */
  double  *rhs; /* 1 entrada para cada restricao, armazena o rhs da restricao */
  char *qrtype; /* 1 "         "   "     "      , pode ser L, G, E respectivamente para <=, >= e = */

  double *dlb, *dub; /* limites inferior e superior de cada variavel */

  ncols=cutting->m; /* total de variaveis */
  nrows=cutting->m; /* total de restricoes */
	
  /* aloca espaco para os arrays que serão usados pelo XPRESS na carga do problema  */

  obj=(double *)malloc(ncols*sizeof(double));
  dlb=(double *)malloc(ncols*sizeof(double));
  dub=(double *)malloc(ncols*sizeof(double));
  mnel=(int *)malloc(ncols*sizeof(int));

  /* arrays com informações sobre as restrições (total de m restricoes)*/
  qrtype = malloc(nrows*sizeof(char));
  rhs=(double *)malloc(nrows*sizeof(double));
  mstart=(int *)malloc(ncols*sizeof(int));/* uma entrada para cada coluna*/
  mrwind=(int *)malloc(cutting->m*sizeof(int));	/* uma para cada entrada nao-nula em A */
  dmatval=(double *)malloc(cutting->m*sizeof(double)); /* idem */

  for(i=0;i<ncols;i++){
    dlb[i]=0.0;    /* limite inferior das variáveis */
    dub[i]=XPRS_PLUSINFINITY;    /* limite superior das variáveis */
  }

  /* configuracao das restricoes */
  for(j=0;j<nrows;j++){
    rhs[j]=cutting->d[j];
    qrtype[j]='G'; 
  }

  /* carrega os coeficientes da coluna x_j */
  nz=0;
  for (j=0;j<ncols;j++){
    obj[j]=1.0; /* coeficiente da var x_j na f.o.*/
    mstart[j]=nz;   /* indice em dmatval da primeira entrada nao-nula da coluna j */ 
    mnel[j]=1;
    mrwind[nz]=j;
    dmatval[nz++]=floor(cutting->L/cutting->w[j]);
  }

  /* carga do modelo*/
  xpress_ret=XPRSloadlp(*probMestre, "Cutting", ncols, nrows,  qrtype, rhs, NULL, obj, \
			mstart, mnel, mrwind, dmatval, dlb, dub);

  /* libera memoria dos vetores usado na carga do LP */
  free(qrtype); free(rhs); free(obj); free(mstart); free(mrwind); free(dmatval); 
  free(dlb);    free(dub); 
  free(mnel);

  return xpress_ret;
}
/* ========================================= */
/* IsIntegerSol                             */
/* ========================================= */
/* obtem o valor das variaveis duais e o valor da solucao. E retorna se sol eh inteira */
int IsIntegerSol(XPRSprob prob, double *z_PMR, double dual[])
{
  int i, inteiro, solfile, colunas;
  double *lambda, *obj;
	      
  /* recupera dados do problema */
  XPRSgetintattrib(prob,XPRS_COLS,&colunas); /* recupera numero de colunas  */

  lambda = (double *)malloc(colunas*sizeof(double));
  obj=(double *)malloc(colunas*sizeof(double));
	
  /* recupera dados do problema */
  XPRSgetobj(prob,obj, 0, colunas-1); /* recupera coef das var na f.o.*/

  /* recupera a solucao do problema */
  XPRSgetintcontrol(prob,XPRS_SOLUTIONFILE,&solfile); /* guarda valor */
  XPRSsetintcontrol(prob,XPRS_SOLUTIONFILE,0);		/* reseta */
  XPRSgetsol(prob,lambda, NULL, dual, NULL);
  XPRSsetintcontrol(prob,XPRS_SOLUTIONFILE,solfile);		/* volta config anterior */
	      
  /* testar se solucao lambda eh inteira. Se for, ver se solucao eh melhor e guardar  */
  inteiro=1;
  *z_PMR=0.0;
  for(i=0;i<colunas;i++){
    if ((lambda[i]-trunc(lambda[i]))>EPSILON)
      inteiro=0;
    (*z_PMR) +=obj[i]*lambda[i];
  }

  free(lambda); free(obj);
  return inteiro;
}

/* ========================================= */
/* adicionaCols                              */
/* ========================================= */
/* adiciona as colunas do pool no prob */
void adicionaCols(XPRSprob prob, double *pool[], int totPool, cuttingInstance *cutting)
{
  int *mstart, *mrwind, i, j, nz;
  double *dmatval, *dlb, *dub, *custo;

  /* var para incluir nova coluna */
  mstart=(int *)malloc((totPool+1)*sizeof(int));
  mrwind=(int *)malloc(cutting->m*totPool*sizeof(int));
  dmatval=(double *)malloc(cutting->m*totPool*sizeof(double));
  dlb=(double *)malloc(totPool*sizeof(double));
  dub=(double *)malloc(totPool*sizeof(double));
  custo=(double *)malloc(totPool*sizeof(double));
		
  nz=0;
  for(j=0;j<totPool;j++){
    mstart[j]=nz;
    dlb[j]=0.0;
    dub[j]=XPRS_PLUSINFINITY;
    for(i=0;i<cutting->m;i++){
      if(pool[j][i]>0){
	mrwind[nz]=i;
	dmatval[nz++]=pool[j][i];
      }
    }
    /* custo da nova coluna */
    custo[j]=1.0;
  }
  mstart[j]=nz;
  /* adiciona as colunas no problema mestre */
  XPRSaddcols(prob,totPool,nz,custo,mstart,mrwind,dmatval,dlb,dub);

  free(mstart); free(mrwind); free(dmatval); free(dlb); free(dub); free(custo);
}

/* ========================================= */
/* CargaInstancia                            */
/* ========================================= */
int CargaInstancia(char *instancia, cuttingInstance *cutting)
{
  FILE *fin;
  int i,j;

  if (!(fin = fopen(instancia, "r"))){
    return(0);
  }
	  
  /* ========================================= */
  /* carga da instancia                        */
  /* ========================================= */
  fscanf(fin, "%d %d",&(cutting->m), &(cutting->L));
  printf("\ninstancia: %s m=%d L=%d", instancia,cutting->m,cutting->L);
  
  /* aloca espaco para carga da instancia */
  cutting->w=(int *) malloc(cutting->m*sizeof(int));
  cutting->d=(int *) malloc(cutting->m*sizeof(int));
  /* carga dos comprimentos e demandas das demandas  */
  for(i=0;i<cutting->m;i++){ 
    fscanf(fin,"%d", &((cutting->w)[i]));
    fscanf(fin,"%d", &((cutting->d)[i]));
  }
  fclose(fin);
  return 1;
}

/* Cria, Carrega e configura cada modelo */
void CriaCarregaConfigura(int formato, XPRSprob *probMestre, XPRSprob *probPricing, cuttingInstance *cutting, solution *sol)
{
  int xpress_ret;

  /* "cria" o problema  mestre reduzido */
  xpress_ret=XPRScreateprob(probMestre);
  if (xpress_ret) 
    errormsg("Main: Erro na initializacao do problema",__LINE__,xpress_ret, *probPricing);
  
  /* "cria" o problema Pricing */
  xpress_ret=XPRScreateprob(probPricing);
  if (xpress_ret) 
    errormsg("Main: Erro na initializacao do problema",__LINE__,xpress_ret, *probPricing);
  
  /* ==================================================================================== */
  /* Carga e configuracao dos parametros de controle do XPRESS para o subprob de pricing  */
  /* ==================================================================================== */
  
  /* limita o tempo de execucao */
  xpress_ret=XPRSsetintcontrol(*probPricing,XPRS_MAXTIME,MAX_CPU_TIME);
  if (xpress_ret) 
    errormsg("Main: Erro ao tentar setar o XPRS_MAXTIME.\n",__LINE__,xpress_ret, *probPricing);
  
  /* carga do modelo de pricing */
  xpress_ret=CargaPricing(probPricing, cutting);
  if (xpress_ret) errormsg("Main: Erro na carga do modelo.",__LINE__,xpress_ret, *probPricing);
  
  /* salva um arquivo ".lp" com o LP original */
  xpress_ret=XPRSwriteprob(*probPricing,"Pricing","l");
  if (xpress_ret) 
    errormsg("Main: Erro na chamada da rotina XPRSwriteprob.\n",__LINE__,xpress_ret, *probPricing);
  
  /* Desabilita o PRESOLVE */
  xpress_ret=XPRSsetintcontrol(*probPricing,XPRS_PRESOLVE,0);
  if (xpress_ret) 
    errormsg("Main: Erro ao desabilitar o presolve.",__LINE__,xpress_ret, *probPricing);
  
  /* impressão para conferência */
  xpress_ret=XPRSsetintcontrol(*probPricing,XPRS_MIPLOG,formato);
  if (xpress_ret) 
    errormsg("Main: Erro ao setar MIPLOG.",__LINE__,xpress_ret, *probPricing);
  
  /* Desabilita heuristica */
  xpress_ret=XPRSsetintcontrol(*probPricing,XPRS_HEURSTRATEGY,0);
  if (xpress_ret) 
    errormsg("Main: Erro ao tentar setar o XPRS_CUTSTRATEGY.\n",__LINE__,xpress_ret, *probPricing);
  
  /* Desabilita a separacao de cortes do XPRESS.  */
  xpress_ret=XPRSsetintcontrol(*probPricing,XPRS_CUTSTRATEGY,0);
  if (xpress_ret) 
    errormsg("Main: Erro ao tentar setar o XPRS_CUTSTRATEGY.\n",__LINE__,xpress_ret, *probPricing);
  
  /* callback para salvar a melhor solucao inteira encontrada */
  xpress_ret=XPRSsetcbintsol(*probPricing,SalvaMelhorSol, sol);
  if (xpress_ret) 
    errormsg("Main: Erro na chamada da rotina XPRSsetcbintsol\
			.\n",__LINE__,xpress_ret, *probPricing);
  
  /* ===================================================================================== */
  /* carga e Configuracao dos parametros de controle do XPRESS para o problema mestre red  */
  /* ===================================================================================== */
  
  /* limita a qtde de colunas a serem incluidas */
  XPRSsetintcontrol(*probMestre,XPRS_EXTRACOLS,50);
  
  /* limita o tempo de execucao */
  xpress_ret=XPRSsetintcontrol(*probMestre,XPRS_MAXTIME,MAX_CPU_TIME);
  if (xpress_ret) 
    errormsg("Main: Erro ao tentar setar o XPRS_MAXTIME.\n",__LINE__,xpress_ret, *probMestre);
  
  /* carga do problema mestre */
  xpress_ret=CargaProblemaMestre(probMestre, cutting);
  if (xpress_ret) errormsg("Main: Erro na carga do modelo.",__LINE__,xpress_ret, *probMestre);
  
  /* Desabilita o PRESOLVE */
  xpress_ret=XPRSsetintcontrol(*probMestre,XPRS_PRESOLVE,0);
  if (xpress_ret) 
    errormsg("Main: Erro ao desabilitar o presolve.",__LINE__,xpress_ret, *probMestre);

  /* salva um arquivo ".lp" com o LP original */
  xpress_ret=XPRSwriteprob(*probMestre,"Mestre","l");
  if (xpress_ret) 
    errormsg("Main: Erro na chamada da rotina XPRSwriteprob.\n",__LINE__,xpress_ret, *probMestre);
}
void ImprimeCols(double *pool[], int totpool, cuttingInstance *cutting){ /* mostra os padroes de corte gerados  */
  int i,j;
  for(j=0; j<totpool;j++){
    printf("\nPadrao de corte (item,qtde): { ");
    for(i=0;i<cutting->m;i++){
      if(pool[j][i])
	printf("(%d,%.2lf) ", i, pool[j][i]);  /* (item, qtde) */
    }
    printf(" }\n");
  }
}
/**********************************************************************************\
 * Rotina copiada de programas exemplo do XPRES ... :(
 * ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 *
 * Name:         errormsg 
 * Purpose:      Display error information about failed subroutines.
 * Arguments:    const char *sSubName   Subroutine name             
 *               int nLineNo            Line number                 
 *               int nErrCode           Error code                  
 * Return Value: None                                               
\**********************************************************************************/
void errormsg(const char *sSubName,int nLineNo,int nErrCode, XPRSprob probMestre)
{
   int nErrNo;             /* Optimizer error number */

   /* Print error message */
   printf("The subroutine %s has failed on line %d\n",sSubName,nLineNo);

   /* Append the error code, if it exists */
   if (nErrCode!=-1) printf("with error code %d\n",nErrCode);

   /* Append Optimizer error number, if available */
   if (nErrCode==32) {
      XPRSgetintattrib(probMestre,XPRS_ERRORCODE,&nErrNo);
      printf("The Optimizer error number is: %d\n",nErrNo);
   }

   /* Free memory, close files and exit */
   XPRSdestroyprob(probMestre);
   XPRSfree();
   exit(nErrCode);
}

/* mensagem de uso do programa */
void showUsage() {
  printf ("Uso: cutting <instancia> <formato-saida>\n");
  printf ("- instancia.....: nome do arquivo de instância\n");
  printf ("- formato-saida.: 0 a 3 (3 = detalhado)\n"); 
}

