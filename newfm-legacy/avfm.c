#include "common.h"
#include "avfm.h"


int main(int argc,char *argv[])
{
  double start_cputime=clock();

  FR_DAT fr;
  CG_PAR cg;

  ReadControl(&cg);
  ReadTop(&cg);
  ReadRange(&cg);


  MAT_DAT mat;

  //cg.n_block=1;
  if(cg.be_sparse==0)
  {
    InitialDenseMatrix(&cg,&mat);
  }
  else if(cg.be_sparse==1)
  {
    InitialSparseMatrix(&cg,&mat);
  }
  else if(cg.be_sparse==2)
  {
    InitialAccumulationMatrix(&cg,&mat);
  }


  if(cg.ba_type==0)
  {
    InitialBsplineSpace(&cg);
  }

  if(cg.be_sparse==0) ReadBinDense(&mat);
  else if(cg.be_sparse==1) ReadBinSparse(&mat);
  else if(cg.be_sparse==2) ReadBinAccumulation(&mat);

  FreeAVSpace(&cg,&mat,&fr);
  //cg.con_p=0;
  if(cg.be_sparse==1) PostOperationSparse(&cg,&mat);
  else if(cg.be_sparse==0) PostOperationDense(&cg,&mat);
  else if(cg.be_sparse==2) PostOperationAccumulationReg(&cg,&mat);

  WriteOutputFM(&cg,&mat);

  //print cpu time used
  double end_cputime=clock();
  double elapsed_cputime=((double) (end_cputime - start_cputime)) / CLOCKS_PER_SEC;
  printf("\n%f seconds used.\n",elapsed_cputime);

  return 0;
}

