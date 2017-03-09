#include "common.h"
#include "range.h"

int main(int argc,char *argv[])
{
  double start_cputime=clock();

  FR_DAT fr;
  CG_PAR cg;

  ReadAuguments(argc,argv,&fr);
  ReadControl(&cg);
  cg.use_3b=0;
  ReadTop(&cg);

  int ReadNextFrameTrr(FR_DAT *const);
  int ReadNextFrameXtc(FR_DAT *const);

  if(fr.trj_type==0) {ReadFirstFrameTrr(&cg,&fr); cg.read_next=ReadNextFrameTrr;}
  else if(fr.trj_type==1)  {ReadFirstFrameXtc(&cg,&fr); cg.read_next=ReadNextFrameXtc;}
  printf("\rStart reading frames."); fflush(stdout);

  InitialRange(&cg);

  InitialLinkedList(&cg,&fr);


  MAT_DAT mat;
  mat.mv=0; mat.mp=0;
  mat.bb=(double *)MemAlloc(3*cg.tolcgn*sizeof(double));
  mat.xx=(double *)MemAlloc(sizeof(double));

  DoFM(&cg,&mat,&fr);
  FreeRangeSapce (&cg, &mat, &fr);
  WriteOutputRange(&cg,&mat);
  
  //print cpu time used
  double end_cputime=clock();
  double elapsed_cputime=((double) (end_cputime - start_cputime)) / CLOCKS_PER_SEC;
  printf("\n%f seconds used.\n",elapsed_cputime);

  return 0;
}

