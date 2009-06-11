#include "column_generation.H"

/* Flow da geração de colunas: (Wolsey pp 188-194)
   lê todos os dados da entrada bla bla bla
   escolhe as colunas
   while (1) {
   resolve o problema restrito (limitante primal do mestre relaxado)
   verifica se é inteira (limitante primal do mestre)
   resolve o problema de pricing (verifica se restrito é dual-viável)
     verifica se a solução do restrito é ótima (pricing)
       se não, gera colunas; continue
     verifica se é inteira (ótima para o PI)
       se sim, break;
       se não, ???;
   }
*/

ColumnGeneration::ColumnGeneration(IntegerProgram &ip, std::vector<IntegerProgram*> &pricing, int test_type=1)
   : ipMestre(ip) {

  /* ========================================= */
  /* inicializacão do XPRESS                   */
  /* ========================================= */
  xpress_ret=XPRSinit("");
  if (xpress_ret)
    errormsg("Main: Erro de inicializacao do XPRESS.\n",
	     __LINE__, xpress_ret, probMestre);

  /* inicializa valores de variaveis globais */
  tPMR = 0;
  it = 0;
  totalInteiros = 0;
  melhorPrimal = +INF;
  this->test_type = test_type;

  printf("\n==========================\nOtimizacao do LP e IP");
}

void ColumnGeneration::solveRestricted() {
  /*========================================================================= */
  /* resolve o problema mestre reduzido                                       */
  /*========================================================================= */
  printf("\n==========================\nOtimizacao do PMR \n");

  /* inicializa o pool de colunas */
  sol.totPool=0;

  /* reinicializa o prob mestre reduzido*/
  XPRSinitglobal(probMestre);

  /* resolve o PMR com o optimizer */
  t1=clock();
  xpress_ret=XPRSminim(probMestre,NULL); /* g = algoritmo de busca B&B, NULL = PL */
  t2=clock();

  tempo=((double)(t2-t1))/CLOCKS_PER_SEC;
  tempoLP+=t2-t1;

  if (xpress_ret)
    errormsg("Main: Erro na chamada da rotina XPRSminim.\n",
	     __LINE__, xpress_ret, probMestre);

  /* se solucao eh inteira --> atualiza bound primal, se melhor */
  if (isIntegerSol(probMestre)){
    totalInteiros++;
    if (z_PMR < melhorPrimal){
      fprintf(stderr,"it %d: Melhor primal: %.1lf\n",it,z_PMR);
      melhorPrimal = z_PMR;
      itmelhorPrimal = it;
    }
  }
  char mestre[128];
  sprintf(mestre,"Mestre%03d",it);
  xpress_ret = XPRSwriteprob(probMestre,mestre,"l");
  if (xpress_ret)
    errormsg("Main: Erro na chamada da rotina XPRSwriteprob.\n",__LINE__,xpress_ret, probMestre);

  printf("\nIteracao: %d\n", it );
  printf("- Valor da solucao otima do PMR =%12.6f \n",z_PMR);
  printf("- Tempo do LP = %lf\n", tempo);
  it++;
}

double *ColumnGeneration::getNewObj(int pindex, int *mindex) {
  int colunas;
  double *cost;

  XPRSgetintattrib(probPricing[pindex], XPRS_COLS, &colunas);
  cost = (double *)calloc(colunas, sizeof(double));

  // std::vector<std::vector<double> > & columns = ipPricing[pindex]->getcolumns();
  std::vector<double> & obj = ipPricing[pindex]->getcost();
  piAk = ipMestre.piProductAk(dual);

  for (int i = 0; i < colunas; i++) {
    cost[i] = obj[i] - piAk[i];
    mindex[i] = i;
  }

  // for (int i = 0; i < colunas; i++) {
  //   cost[i] = obj[i];
    // for (int j = 0; j < ipPricing[pindex]->getnrows(); j++) {
    //   cost[i] -= dual[j] * columns[i][j];
    // }
    // mindex[i] = i;
    // fprintf(stderr,"Custo reduzido da coluna %d é %6.1lf\n",i,cost[i]);
    // fprintf(stderr,"Objetivo da variavel %d é %6.1lf\n",i,obj[i]);
  // }

  return cost;
}

int ColumnGeneration::solvePricing() {
  /*========================================================================= */
  /* resolve o problema de pricing                                            */
  /*========================================================================= */
  printf("\n==========================\nOtimizacao do Pricing \n");

  // podem ser vários subproblemas de pricing
  int nrows, ncols, *mindex;
  int columns_added = 0;

  // tem pelos menos um problema de pricing
  XPRSgetintattrib(probPricing[0], XPRS_ROWS, &nrows);
  XPRSgetintattrib(probPricing[0], XPRS_COLS, &ncols);
  mindex = (int *) calloc(ncols, sizeof(int));
  
  for (int k = 0; k < (int)probPricing.size(); k++) {
    /* calcula os novos coeficientes das variaveis
       na funcao objetivo do subproblema de pricing */
    double* obj_pricing = getNewObj(k, mindex);

    /* resolve o sub de pricing */
    /* uma vez maximizado, perde a base inicial, e o problema tem que
       ser carregado de novo? */

    /* reinicializa o prob de pricing */
    XPRSinitglobal(probPricing[k]);
    
    /* seta novos coeficientes na funcao objetivo */
    XPRSchgobj(probPricing[k], ncols, mindex, obj_pricing);
    
    /* reinicializa a melhor solucao */
    sol.zstar = XPRS_PLUSINFINITY;
    xpress_ret = XPRSsetdblcontrol(probPricing[k], XPRS_MIPABSCUTOFF, sol.zstar);

    if(test_type == 1) { // usa XPRESS pra resolver o PLI
      t1=clock();
      xpress_ret = XPRSminim(probPricing[k], "g"); /* g = algoritmo de busca B&B
						    NULL = PL */
      t2=clock();
    }
    else { // test_type == 2: usa W&W pra resolver
      t1=clock();
      sol.zstar = ipPricing[k]->solve(piAk,sol.xstar);
      t2=clock();
    }
    
    
    /* verifica se coluna deve ser adicionada ao PMR */
    double zstar = sol.zstar - dual[nrows+k];
    if (fabs(zstar) < EPSILON) {
      ;
    } else if (sol.zstar - dual[nrows+k] < -EPSILON) {
      int *mstart, *mrwind, nz;
      double *dmatval, *dlb, *dub, *custo;
      double c = 0.0;
      
      for (int i = 0; i < (int)sol.xstar.size(); i++)
	c += sol.xstar[i] * ipPricing[k]->getcost()[i];
      ipMestre.addcol(sol.xstar, k, c,nz,&custo,&mstart,&mrwind,&dmatval,&dlb,&dub);
      XPRSaddcols(probMestre,1,nz,custo,mstart,mrwind,dmatval,dlb,dub);

      columns_added++;
    }

    tempo=((double)(t2-t1))/CLOCKS_PER_SEC;
    tempoPricing+=t2-t1;

    /* imprime resultado do pricing */
    printf("\n==========================\nResultado do Pricing %d", k);
    printf("\nZstar (custo reduzido): %lf", sol.zstar - dual[nrows+k]);
    printf("\nTempo do pricing: %lf\n", tempo);
    printf("\nColuna Adicionada? %s\n", columns_added ? "S" : "N");
  }

  free(mindex);

  if (columns_added == 0) {
    fprintf(stderr,"it %d: Limitante dual encontrado\n",it);
  }

  return columns_added;
}

/* Cria, Carrega e configura cada modelo */
void ColumnGeneration::configureModel(int formato,
				      IntegerProgram& ip,
				      std::vector<IntegerProgram*>& pricing) {
  int xpress_ret;

  ipPricing = pricing;

  /* "cria" o problema  mestre reduzido */
  xpress_ret = XPRScreateprob(&probMestre);
  if (xpress_ret)
    errormsg("Main: Erro na initializacao do problema",
	     __LINE__,xpress_ret, probMestre);

  /* "cria" o problema Pricing */
  probPricing.resize(pricing.size());
  for (int i = 0; i < (int)pricing.size(); i++) {
    char nome[16];
    sprintf(nome, "pricing%d", i);
    xpress_ret = XPRScreateprob(&probPricing[i]);
    if (xpress_ret)
      errormsg("Main: Erro na initializacao do problema",
	       __LINE__,xpress_ret, probPricing[i]);

    /* ==================================================================================== */
    /* Carga e configuracao dos parametros de controle do XPRESS para o subprob de pricing  */
    /* ==================================================================================== */

    /* limita o tempo de execucao */
    xpress_ret=XPRSsetintcontrol(probPricing[i], XPRS_MAXTIME,MAX_CPU_TIME);
    if (xpress_ret)
      errormsg("Main: Erro ao tentar setar o XPRS_MAXTIME.\n",__LINE__,xpress_ret, probPricing[i]);

    /* carga do modelo de pricing */
    xpress_ret=loadModel(*(pricing[i]), probPricing[i], nome);
    if (xpress_ret) errormsg("Main: Erro na carga do modelo.",__LINE__,xpress_ret, probPricing[i]);

    /* salva um arquivo ".lp" com o LP original */
    xpress_ret=XPRSwriteprob(probPricing[i], nome, "l");
    if (xpress_ret)
      errormsg("Main: Erro na chamada da rotina XPRSwriteprob.\n",__LINE__,xpress_ret, probPricing[i]);

    /* Desabilita o PRESOLVE */
    xpress_ret=XPRSsetintcontrol(probPricing[i],XPRS_PRESOLVE,0);
    if (xpress_ret)
      errormsg("Main: Erro ao desabilitar o presolve.",__LINE__,xpress_ret, probPricing[i]);

    /* impressão para conferência */
    xpress_ret=XPRSsetintcontrol(probPricing[i],XPRS_MIPLOG,formato);
    if (xpress_ret)
      errormsg("Main: Erro ao setar MIPLOG.",__LINE__,xpress_ret, probPricing[i]);

    /* Desabilita heuristica */
    xpress_ret=XPRSsetintcontrol(probPricing[i],XPRS_HEURSTRATEGY,0);
    if (xpress_ret)
      errormsg("Main: Erro ao tentar setar o XPRS_CUTSTRATEGY.\n",__LINE__,xpress_ret, probPricing[i]);

    /* Desabilita a separacao de cortes do XPRESS.  */
    xpress_ret=XPRSsetintcontrol(probPricing[i],XPRS_CUTSTRATEGY,0);
    if (xpress_ret)
      errormsg("Main: Erro ao tentar setar o XPRS_CUTSTRATEGY.\n",__LINE__,xpress_ret, probPricing[i]);

    /* callback para salvar a melhor solucao inteira encontrada */
    xpress_ret=XPRSsetcbintsol(probPricing[i],salvaMelhorSol, &sol);
    if (xpress_ret)
      errormsg("Main: Erro na chamada da rotina XPRSsetcbintsol\
  			.\n",__LINE__,xpress_ret, probPricing[i]);

  }

  /* ===================================================================================== */
  /* carga e Configuracao dos parametros de controle do XPRESS para o problema mestre red  */
  /* ===================================================================================== */

  /* limita a qtde de colunas a serem incluidas */
  XPRSsetintcontrol(probMestre,XPRS_EXTRACOLS,50);

  /* limita o tempo de execucao */
  xpress_ret = XPRSsetintcontrol(probMestre,XPRS_MAXTIME,MAX_CPU_TIME);
  if (xpress_ret)
    errormsg("Main: Erro ao tentar setar o XPRS_MAXTIME.\n",__LINE__,xpress_ret, probMestre);

  /* carga do problema mestre */
  xpress_ret = loadModel(ip, probMestre, "Mestre");
  if (xpress_ret) errormsg("Main: Erro na carga do modelo.",__LINE__,xpress_ret, probMestre);

  /* aloca variaveis duais */
  dual = (double *)calloc(ip.getnrows(), sizeof(double));

  /* Desabilita o PRESOLVE */
  xpress_ret = XPRSsetintcontrol(probMestre,XPRS_PRESOLVE,0);
  if (xpress_ret)
    errormsg("Main: Erro ao desabilitar o presolve.",__LINE__,xpress_ret, probMestre);

  /* salva um arquivo ".lp" com o LP original */
  xpress_ret = XPRSwriteprob(probMestre,"Mestre","l");
  if (xpress_ret)
    errormsg("Main: Erro na chamada da rotina XPRSwriteprob.\n",__LINE__,xpress_ret, probMestre);
}

/* ========================================= */
/* IsIntegerSol                              */
/* ========================================= */
/* obtem o valor das variaveis duais e o valor da solucao. E retorna se sol eh inteira */
bool ColumnGeneration::isIntegerSol(XPRSprob prob) {
  int i, inteiro, solfile, colunas;
  double *lambda, *obj;

  /* recupera dados do problema */
  XPRSgetintattrib(prob,XPRS_COLS,&colunas); /* recupera numero de colunas  */

  lambda = (double *)malloc(colunas*sizeof(double));
  obj = (double *)malloc(colunas*sizeof(double));

  /* recupera dados do problema */
  XPRSgetobj(prob,obj, 0, colunas-1); /* recupera coef das var na f.o.*/

  /* recupera a solucao do problema */
  XPRSgetintcontrol(prob,XPRS_SOLUTIONFILE,&solfile); /* guarda valor */
  XPRSsetintcontrol(prob,XPRS_SOLUTIONFILE,0);	      /* reseta */
  XPRSgetsol(prob, lambda, NULL, dual, NULL);
  XPRSsetintcontrol(prob,XPRS_SOLUTIONFILE,solfile);  /* volta config anterior */

  /* testar se solucao lambda eh inteira. Se for, ver se solucao eh melhor e guardar  */
  inteiro=1;
  z_PMR=0.0;
  for(i=0;i<colunas;i++){
    if ((lambda[i]-trunc(lambda[i]))>EPSILON && (ceil(lambda[i]) - lambda[i])>EPSILON)
      inteiro=0;
    z_PMR +=obj[i]*lambda[i];
  }
  fprintf(stderr,"it %d: z = %.3lf\n",it,z_PMR);
  free(lambda);
  free(obj);

  return inteiro;
}

/* Esta é a função mais importante: é ela que transforma a instância
   em um problema do XPRESS. A idéia é fazê-la o mais genérica
   possível e transferir responsabilidade para métodos da
   instância. Com isso, poderemos resolver qualquer problema, desde
   que ele saiba se transformar em uma entrada do XPRESS (via método
   getParam)

   Mais importante, podemos usar o mesmo método para carregar o
   problema mestre e o problema de pricing. Os antigos
   CargaProblemaPricing e CargaProblemaMestre serão transformados em
   loadModel
*/
int ColumnGeneration::loadModel(IntegerProgram& ip, XPRSprob prob, std::string name) {
  int ncol, nrow,nmip;
  char* rowtype,* miptype;
  double* rhs,* obj,* lb,* ub,* matval;
  int* colbeg,* rowidx,*mipcol;
  bool relaxed;

  ip.getParam(ncol,nrow,&rowtype,&rhs,&obj,&colbeg,&rowidx,&matval,&lb,
	      &ub,nmip,&miptype,&mipcol,relaxed);

  if (relaxed) {
    xpress_ret = XPRSloadlp(prob,name.c_str(),ncol,nrow,rowtype,rhs,
			    NULL,obj,colbeg,NULL,rowidx,matval,lb,ub);
  } else { // integer
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

/* =====================================================================
 * Rotina (callback) para salvar a  melhor solucao. Roda para todo nó
 * onde acha solucao inteira.
 * callback para XPRSsetcbintsol
 * Autor: Cid Carvalho de Souza
 * =====================================================================*/
void XPRS_CC salvaMelhorSol(XPRSprob prob, void *psol) {
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

//   if (1-objval<0 && sol->totPool<MAX_POOL){
//     (sol->pool)[sol->totPool] = std::vector<double>(cols);
//     for(i=0;i<cols;i++)
//       (sol->pool)[sol->totPool][i]=x[i];
//     (sol->totPool)++;
//   }

  /* testa se a solução é viável */
  viavel = true;

  /* se a solucao tiver custo melhor que a melhor solucao disponivel entao salva */
  if ((objval < sol->zstar) && viavel ) {
    sol->xstar = std::vector<double>(cols, 0.0);
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
