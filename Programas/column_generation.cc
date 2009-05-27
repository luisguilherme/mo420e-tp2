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

  if (isIntegerSol(probMestre, &z_PMR, dual)){ /* obtem valor das var duais e valor da solucao */
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
  obj_pricing = instance.getNewObj(dual,mindex);
    
  /* resolve o sub de pricing */
  /* uma vez maximizado, perde a base inicial, e o problema tem que ser carregado de novo? */
  XPRSinitglobal(probPricing); /* reinicializa o prob de pricing */
  
  XPRSchgobj(probPricing,instance.m,mindex,obj_pricing); /* seta novos coeficientes na funcao objetivo */
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
