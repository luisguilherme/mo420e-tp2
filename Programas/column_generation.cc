#include "column_generation.H"

void ColumnGeneration::solveReduced() {
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
    errormsg("Main: Erro na chamada da rotina XPRSminim.\n",__LINE__,xpress_ret, probMestre);

  if (isIntegerSol(probMestre)){ /* obtem valor das var duais e valor da solucao */
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

}

void ColumnGeneration::solvePricing(Instance& instance) {
  /*========================================================================= */
  /* resolve o problema de pricing                                            */
  /*========================================================================= */
  printf("\n==========================\nOtimizacao do Pricing \n");
  
  /* calcula os novos coeficientes das variaveis na funcao objetivo do subproblema de pricing */
  double* obj_pricing = instance.getNewObj(dual,mindex);
    
  /* resolve o sub de pricing */
  /* uma vez maximizado, perde a base inicial, e o problema tem que ser carregado de novo? */
  XPRSinitglobal(probPricing); /* reinicializa o prob de pricing */
  
  XPRSchgobj(probPricing,instance.dim(),mindex,obj_pricing); /* seta novos coeficientes na funcao objetivo */
  /* reinicializa a melhor solucao */
  sol.zstar=XPRS_MINUSINFINITY; 
  xpress_ret=XPRSsetdblcontrol(probPricing,XPRS_MIPABSCUTOFF,sol.zstar);
  
  t1=clock();
  xpress_ret=XPRSmaxim(probPricing,"g"); /* g = algoritmo de busca B&B, NULL = PL */
  t2=clock();
  
  tempo=((double)(t2-t1))/CLOCKS_PER_SEC;
  tempoPricing+=t2-t1;
  
  /* imprime resultado do pricing */
  printf("\n==========================\nResultado do Pricing");
  printf("\nZstar (custo reduzido): %lf", 1-sol.zstar);
  printf("\nTempo do pricing: %lf\n", tempo);
}

ColumnGeneration::ColumnGeneration(Instance& instance) {
  /* ========================================= */
  /* inicializacão do XPRESS                   */
  /* ========================================= */
  xpress_ret=XPRSinit("");
  if (xpress_ret) 
    errormsg("Main: Erro de inicializacao do XPRESS.\n",__LINE__,xpress_ret, probMestre);

  /* inicializa valores de variaveis globais */
  tPMR=0;

  printf("\n==========================\nOtimizacao do LP e IP");
}

/* Cria, Carrega e configura cada modelo */
void ColumnGeneration::configureModel(int formato, Instance& instance) {
  int xpress_ret;

  /* "cria" o problema  mestre reduzido */
  xpress_ret=XPRScreateprob(&probMestre);
  if (xpress_ret) 
    errormsg("Main: Erro na initializacao do problema",__LINE__,xpress_ret, probPricing);
  
  /* "cria" o problema Pricing */
  xpress_ret=XPRScreateprob(&probPricing);
  if (xpress_ret) 
    errormsg("Main: Erro na initializacao do problema",__LINE__,xpress_ret, probPricing);
  
  /* ==================================================================================== */
  /* Carga e configuracao dos parametros de controle do XPRESS para o subprob de pricing  */
  /* ==================================================================================== */
  

  /* ALTERAR DEPOIS!!!  */
  /* limita o tempo de execucao */
  xpress_ret=XPRSsetintcontrol(probPricing,XPRS_MAXTIME,MAX_CPU_TIME);
  if (xpress_ret) 
    errormsg("Main: Erro ao tentar setar o XPRS_MAXTIME.\n",__LINE__,xpress_ret, probPricing);
  
  /* carga do modelo de pricing */
  xpress_ret=loadModel(instance,probPricing,"Pricing");
  if (xpress_ret) errormsg("Main: Erro na carga do modelo.",__LINE__,xpress_ret, probPricing);
  
  /* salva um arquivo ".lp" com o LP original */
  xpress_ret=XPRSwriteprob(probPricing,"Pricing","l");
  if (xpress_ret) 
    errormsg("Main: Erro na chamada da rotina XPRSwriteprob.\n",__LINE__,xpress_ret, probPricing);
  
  /* Desabilita o PRESOLVE */
  xpress_ret=XPRSsetintcontrol(probPricing,XPRS_PRESOLVE,0);
  if (xpress_ret) 
    errormsg("Main: Erro ao desabilitar o presolve.",__LINE__,xpress_ret, probPricing);
  
  /* impressão para conferência */
  xpress_ret=XPRSsetintcontrol(probPricing,XPRS_MIPLOG,formato);
  if (xpress_ret) 
    errormsg("Main: Erro ao setar MIPLOG.",__LINE__,xpress_ret, probPricing);
  
  /* Desabilita heuristica */
  xpress_ret=XPRSsetintcontrol(probPricing,XPRS_HEURSTRATEGY,0);
  if (xpress_ret) 
    errormsg("Main: Erro ao tentar setar o XPRS_CUTSTRATEGY.\n",__LINE__,xpress_ret, probPricing);
  
  /* Desabilita a separacao de cortes do XPRESS.  */
  xpress_ret=XPRSsetintcontrol(probPricing,XPRS_CUTSTRATEGY,0);
  if (xpress_ret) 
    errormsg("Main: Erro ao tentar setar o XPRS_CUTSTRATEGY.\n",__LINE__,xpress_ret, probPricing);
  
  /* callback para salvar a melhor solucao inteira encontrada */
  xpress_ret=XPRSsetcbintsol(probPricing,salvaMelhorSol, &sol);
  if (xpress_ret) 
    errormsg("Main: Erro na chamada da rotina XPRSsetcbintsol\
			.\n",__LINE__,xpress_ret, probPricing);
  
  /* ===================================================================================== */
  /* carga e Configuracao dos parametros de controle do XPRESS para o problema mestre red  */
  /* ===================================================================================== */
  
  /* limita a qtde de colunas a serem incluidas */
  XPRSsetintcontrol(probMestre,XPRS_EXTRACOLS,50);
  
  /* limita o tempo de execucao */
  xpress_ret=XPRSsetintcontrol(probMestre,XPRS_MAXTIME,MAX_CPU_TIME);
  if (xpress_ret) 
    errormsg("Main: Erro ao tentar setar o XPRS_MAXTIME.\n",__LINE__,xpress_ret, probMestre);
  
  /* carga do problema mestre */
  xpress_ret=loadModel(instance, probMestre, "Mestre");
  if (xpress_ret) errormsg("Main: Erro na carga do modelo.",__LINE__,xpress_ret, probMestre);
  
  /* Desabilita o PRESOLVE */
  xpress_ret=XPRSsetintcontrol(probMestre,XPRS_PRESOLVE,0);
  if (xpress_ret) 
    errormsg("Main: Erro ao desabilitar o presolve.",__LINE__,xpress_ret, probMestre);

  /* salva um arquivo ".lp" com o LP original */
  xpress_ret=XPRSwriteprob(probMestre,"Mestre","l");
  if (xpress_ret) 
    errormsg("Main: Erro na chamada da rotina XPRSwriteprob.\n",__LINE__,xpress_ret, probMestre);
}




/* =====================================================================
 * Rotina (callback) para salvar a  melhor solucao.  Roda para todo nó
 * onde acha solucao inteira.
 * callback para XPRSsetcbintsol
 * Autor: Cid Carvalho de Souza
 * =====================================================================
*/
void XPRS_CC salvaMelhorSol(XPRSprob prob, void *psol)
{
  int i, cols, node, xpress_ret;
  double objval, *x;
  bool viavel;
  solution* sol;

  sol=(solution*) psol;

  /* pega o numero do nó corrente */
  xpress_ret=XPRSgetintattrib(prob,XPRS_NODES,&node);
  if (xpress_ret) 
    errormsg("salvaMelhorSol: rotina XPRSgetintattrib.\n",__LINE__,xpress_ret, prob);

  /* recupera numero de colunas do prob */
  xpress_ret=XPRSgetintattrib(prob,XPRS_COLS,&cols);
  if (xpress_ret) 
    errormsg("salvaMelhorSol: rotina XPRSgetintattrib\n",__LINE__,xpress_ret, prob);

  /* recupera solucao */
  x=(double *)malloc(cols*sizeof(double));

  xpress_ret=XPRSgetdblattrib(prob,XPRS_LPOBJVAL,&objval);
  if (xpress_ret) 
    errormsg("salvaMelhorSol: rotina XPRSgetdblattrib\n",__LINE__,xpress_ret, prob);

  xpress_ret=XPRSgetsol(prob,x,NULL,NULL,NULL);
  if (xpress_ret) 
    errormsg("salvaMelhorSol: Erro na chamada da rotina XPRSgetsol\
			\n",__LINE__,xpress_ret, prob);

  /* guarda no pool de colunas se o custo reduzido da solucao inteira é negativo */
  if (1-objval<0 && sol->totPool<MAX_POOL){
    (sol->pool)[sol->totPool] = std::vector<double>(cols);
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
      errormsg("salvaMelhorSol: XPRSsetdblcontrol.\n",__LINE__,xpress_ret, prob);
  }   
  free(x);
  return;
}

/* ========================================= */
/* IsIntegerSol                             */
/* ========================================= */
/* obtem o valor das variaveis duais e o valor da solucao. E retorna se sol eh inteira */
inline bool ColumnGeneration::isIntegerSol(XPRSprob prob) {
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
  z_PMR=0.0;
  for(i=0;i<colunas;i++){
    if ((lambda[i]-trunc(lambda[i]))>EPSILON)
      inteiro=0;
    z_PMR +=obj[i]*lambda[i];
  }

  free(lambda); free(obj);
  return inteiro;
}


/* Esta é a função mais importante: é ela que transforma a instância
   em um problema do XPRESS. A idéia é fazê-la o mais genérica
   possível e transferir responsabilidade para métodos da
   instância. Com isso, poderemos resolver qualquer problema, desde
   que ele saiba se transformar em uma entrada do XPRESS (via método
   getParam) 

   Mais importante, podemos usar o mesmo método para carregar o
   problema mestre e o problema de pricing. Neste caso não sei se o
   mais propício é ele receber "Instance" ou "IntegerProgram" como
   argumento.

   Os antigos CargaProblemaPricing e CargaProblemaMestre serão
   transformados em loadModel   
*/
int ColumnGeneration::loadModel(Instance& instance,XPRSprob prob,std::string name) {
  int ncol, nrow,nmip; 
  char* rowtype,* miptype;
  double* rhs,* obj,* lb,* ub,* matval;
  int* colbeg,* rowidx,*mipcol;
  bool relaxed;
  instance.getParam(ncol,nrow,rowtype,rhs,obj,colbeg,rowidx,matval,lb,
		    ub,nmip,miptype,mipcol,relaxed);
  if (relaxed) {
    xpress_ret = XPRSloadlp(prob,name.c_str(),ncol,nrow,rowtype,rhs,
			    NULL,obj,colbeg,NULL,rowidx,matval,lb,ub);
  }
  else { // integer 
    xpress_ret = XPRSloadglobal(prob,name.c_str(),ncol,nrow,rowtype,rhs,
				NULL,obj,colbeg,NULL,rowidx,matval,lb,
				ub,nmip,0,miptype,mipcol,NULL,NULL,NULL,
				NULL,NULL);
  }
  return(xpress_ret);
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



