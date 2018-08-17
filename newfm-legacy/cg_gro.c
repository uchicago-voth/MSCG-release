#include "common.h"
#include "cg_gro.h"


int main(int argc,char *argv[])
{
  double start_cputime=clock(); //start to compute cputime

  FR_DAT fr; //trajectory frame data
  CG_PAR cg; //CG parameters

  ReadAuguments(argc,argv,&fr); //read command flags
  ReadControl(&cg); //read control.in
  ReadTop(&cg); //read top.in


  if(fr.trj_type==0) {ReadFirstFrameTrr(&cg,&fr); cg.read_next=ReadNextFrameTrr;}
  else if(fr.trj_type==1)  {ReadFirstFrameXtc(&cg,&fr); cg.read_next=ReadNextFrameXtc;}
  printf("\rStart reading frames."); fflush(stdout);

  MakeGro(&cg,&fr);

  //print cpu time used
  double end_cputime=clock();
  double elapsed_cputime=((double) (end_cputime - start_cputime)) / CLOCKS_PER_SEC;
  printf("\n%f seconds used.\n",elapsed_cputime);

  return 0;
}

