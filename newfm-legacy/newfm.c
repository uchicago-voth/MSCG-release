//Force matching main program

#include "common.h"

int main(int argc,char *argv[])
{
  double start_cputime=clock(); //start to compute cputime

  FR_DAT fr; //trajectory frame data
  CG_PAR cg; //CG parameters

  ReadAuguments(argc,argv,&fr); //read command flags
  ReadControl(&cg); //read control.in
  ReadTop(&cg); //read top.in
  ReadRange(&cg); //read rmin.in/rmin_b.in

  if(cg.tollj2>0 || cg.tolb2>0 || cg.tola2>0 || cg.told2>0) ReadExt(&cg); //read table.in
  if(cg.pres_con==1) ReadVirial(&fr,&cg); //read p_con.in


  //determine functions to read the trajectory
  if(fr.trj_type==0) {ReadFirstFrameTrr(&cg,&fr); cg.read_next=ReadNextFrameTrr;}
  else if(fr.trj_type==1) /*{printf("Xtc not ready yet!\n"); exit (0);} */{ReadFirstFrameXtc(&cg,&fr); cg.read_next=ReadNextFrameXtc;}
  printf("\rStart reading frames."); fflush(stdout);

  //initial cell lined-lists
  InitialLinkedList(&cg,&fr);
  if(cg.use_3b>0) InitialLinkedList_3(&cg,&fr);

  MAT_DAT mat; //matrix data

  InitialFunction(&cg); //determine other functions

  //Initial matrix
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
  else {printf("Wrong value for use_sparse!\n"); exit(1);}

  //Initial b-spline work spaces
  if(cg.ba_type==0)
  {
    InitialBsplineSpace(&cg);
  }

  FILE *sol=fopen("sol_info.out","w");
  fprintf(sol,"mm:%d; nn:%d;\n",mat.mm,mat.nn);
  fclose(sol);

  //FM main function
  DoFM(&cg,&mat,&fr);

  //Free some spaces after FM
  FreeFMSpace(&cg,&mat,&fr);

  //Calculate final results
  if(cg.be_sparse==1) PostOperationSparse(&cg,&mat);
  else if(cg.be_sparse==0) PostOperationDense(&cg,&mat);
  else if(cg.be_sparse==2) PostOperationAccumulation(&cg,&mat);

  //write output files
  WriteOutputFM(&cg,&mat);

  //print cpu time used
  double end_cputime=clock();
  double elapsed_cputime=((double) (end_cputime - start_cputime)) / CLOCKS_PER_SEC;
  printf("\n%f seconds used.\n",elapsed_cputime);
//_PGOPTI_Prof_Dump_All();
  return 0;
}

