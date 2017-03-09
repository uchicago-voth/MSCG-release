#define VERYLARGE 1000.0

void InitialRmin(int tol, double *const rmin, double *const rmax, int *const num)
{
  int i;
  for(i=0;i<tol;i++)
  {
    rmin[i]=VERYLARGE;
    rmax[i]=-VERYLARGE;
    num[i]=i+1;
  }
}

void FillRHSNULL(int shift_i, int site_i, MAT_DAT *const mat, FR_DAT *const fr)
{
}


void InitialRange(CG_PAR *const cg)
{  
  
  cg->tolb=CalTol(cg->bontype,TOL_TWO(cg->tolcgtype));
  cg->bon_m=(int *)MemAlloc(cg->tolb*sizeof(int));
  cg->bonnum=(int *)MemAlloc(cg->tolb*sizeof(int));
  CalM(cg->bontype,TOL_TWO(cg->tolcgtype),cg->bon_m);

  cg->tola=CalTol(cg->angtype,TOL_THREE(cg->tolcgtype));
  cg->ang_m=(int *)MemAlloc(cg->tola*sizeof(int));
  cg->angnum=(int *)MemAlloc(cg->tola*sizeof(int));
  CalM(cg->angtype,TOL_THREE(cg->tolcgtype),cg->ang_m);

  cg->told=CalTol(cg->dihtype,TOL_FOUR(cg->tolcgtype));
  cg->dih_m=(int *)MemAlloc(cg->told*sizeof(int));
  cg->dihnum=(int *)MemAlloc(cg->told*sizeof(int));
  CalM(cg->dihtype,TOL_FOUR(cg->tolcgtype),cg->dih_m);

  cg->ljnum=(int *)MemAlloc(TOL_TWO(cg->tolcgtype)*sizeof(int));
  cg->tollj=TOL_TWO(cg->tolcgtype);

  cg->rmin=(double *)MemAlloc(cg->tollj*sizeof(double));
  cg->rmax=(double *)MemAlloc(cg->tollj*sizeof(double));
  cg->rmin_a=(double *)MemAlloc(cg->tola*sizeof(double));
  cg->rmax_a=(double *)MemAlloc(cg->tola*sizeof(double));
  cg->rmin_b=(double *)MemAlloc(cg->tolb*sizeof(double));
  cg->rmax_b=(double *)MemAlloc(cg->tolb*sizeof(double));
  cg->rmin_d=(double *)MemAlloc(cg->told*sizeof(double));
  cg->rmax_d=(double *)MemAlloc(cg->told*sizeof(double));

  InitialRmin(cg->tollj,cg->rmin,cg->rmax,cg->ljnum);
  InitialRmin(cg->tolb,cg->rmin_b,cg->rmax_b,cg->bonnum);
  InitialRmin(cg->tola,cg->rmin_a,cg->rmax_a,cg->angnum);
  InitialRmin(cg->told,cg->rmin_d,cg->rmax_d,cg->dihnum);

  cg->tollj1=cg->tollj;
  cg->tolb1=cg->tolb;
  cg->tola1=cg->tola;
  cg->told1=cg->told;
  cg->tollj2=0;
  cg->tolb2=0;
  cg->tola2=0;
  cg->told2=0;

  cg->grid_n=(int *)MemAlloc(sizeof(int));
  cg->grid_bn=(int *)MemAlloc(sizeof(int));
  cg->grid_an=(int *)MemAlloc(sizeof(int));
  cg->grid_dn=(int *)MemAlloc(sizeof(int));

  cg->be_sparse=0;
  cg->ba_type=1;
  cg->cutoff2=VERYLARGE*VERYLARGE;

  void FMTwoBodyRange(CG_PAR *const CG, TYPE_INFO *const, FR_DAT *const, MAT_DAT *const);
  void FMThreeBodyRange(CG_PAR *const CG, TYPE_INFO *const, FR_DAT *const, MAT_DAT *const);
  void FMFourBodyRange(CG_PAR *const CG, TYPE_INFO *const, FR_DAT *const, MAT_DAT *const);

  cg->fm_nb=FMTwoBodyRange;
  cg->fm_bon=FMTwoBodyRange;
  if(cg->ang_type==0) cg->fm_ang=FMTwoBodyRange;
  else cg->fm_ang=FMThreeBodyRange;
  if(cg->dih_type==0) cg->fm_dih=FMTwoBodyRange;
  else cg->fm_dih=FMFourBodyRange;


  void SetZeroNull(MAT_DAT *const);
  cg->set_zero=SetZeroNull;

  void NullOperation(CG_PAR *const, MAT_DAT *const);
  cg->mat_op=NullOperation;

  cg->fr_op=FrameOperationNull;

  void FillVirial(MAT_DAT *const, FR_DAT *const);
  cg->fill_virial=FillVirial;

  cg->fill_rhs=FillRHSNULL;  

}

void SetZeroNull(MAT_DAT *const mat)
{}

void NullOperation(CG_PAR *const cg, MAT_DAT *const mat)
{}

void FMTwoBodyRange(CG_PAR *const cg, TYPE_INFO *const info, FR_DAT *const fr, MAT_DAT *const mat)
{
  double rr2,rx[3];
  rr2=CalDis(info->k,info->l,fr,rx);
  info->rr=sqrt(rr2);

  if(info->rmin[info->mi]>info->rr) info->rmin[info->mi]=info->rr;
  if(info->rmax[info->mi]<info->rr) info->rmax[info->mi]=info->rr;
}

void FMThreeBodyRange(CG_PAR *const cg, TYPE_INFO *const info, FR_DAT *const fr, MAT_DAT *const mat)
{
  double rx1[3],rx2[3],rr1,rr2;
  double cos_theta=CalAngleCos(info->j,info->k,info->l,fr,rx1,rx2,&rr1,&rr2);
  info->rr=acos(cos_theta)*R2D;


  if(info->rmin[info->mi]>info->rr) info->rmin[info->mi]=info->rr;
  if(info->rmax[info->mi]<info->rr) info->rmax[info->mi]=info->rr;
}

void FMFourBodyRange(CG_PAR *const cg, TYPE_INFO *const info, FR_DAT *const fr, MAT_DAT *const mat)
{
  double rx1[3],rx2[3],rx3[3],pb[3],pc[3],rpb1,rpc1,pbpc,s;
  info->rr=CalDihedral(info->i,info->j,info->k,info->l,fr,rx1,rx2,rx3,pb,pc,&rpb1,&rpc1,&pbpc,&s)*R2D+180.0;

  if(info->rmin[info->mi]>info->rr) info->rmin[info->mi]=info->rr;
  if(info->rmax[info->mi]<info->rr) info->rmax[info->mi]=info->rr;
}


void OutputRange(enum INTER_TYPE type, CG_PAR *const cg, MAT_DAT *const mat, TYPE_INFO *const info, FILE *const br_out)
{
  if(type==non || type==bon) fprintf(br_out,"%s %s ",cg->name[info->i],cg->name[info->j]);
  else if(type==ang) fprintf(br_out,"%s %s %s ",cg->name[info->i],cg->name[info->j],cg->name[info->k]);
  else if(type==dih) fprintf(br_out,"%s %s %s %s ",cg->name[info->i],cg->name[info->j],cg->name[info->k],cg->name[info->l]);

  if(fabs(info->rmax[info->mi]+VERYLARGE)<VERYSMALL_F || (info->rmin[info->mi]>cg->cutoff && type==non)) {info->rmax[info->mi]=-1.0; info->rmin[info->mi]=-1.0;}
  else if(info->rmax[info->mi]>cg->cutoff && type==non) info->rmax[info->mi]=cg->cutoff;

  fprintf(br_out,"%lf %lf\n",info->rmin[info->mi],info->rmax[info->mi]);
}


void WriteOutputRange(CG_PAR *const cg, MAT_DAT *const mat)
{
  FILE *nb_out,*bon_out;
  nb_out=OpenFile("rmin.in","w");
  bon_out=OpenFile("rmin_b.in","w");

  void OutputRangeOne(enum INTER_TYPE, CG_PAR *const, MAT_DAT *const, TYPE_INFO *const, FILE *const);
  void (*output)(enum INTER_TYPE, CG_PAR *const, MAT_DAT *const, TYPE_INFO *const, FILE *const);
  output=OutputRange;

  WriteOutput(cg,mat,nb_out,bon_out,output);

  fclose(nb_out); fclose(bon_out);
}

void FreeRangeSapce (CG_PAR *const cg, MAT_DAT *const mat, FR_DAT *const fr)
{
  int i;

  free(cg->blist_n);
  free(cg->alist_n);
  free(cg->dlist_n);
  for(i=0;i<cg->tolcgn;i++)
  {
    free(cg->blist_i[i]);
    free(cg->alist_i[i]);
    free(cg->dlist_i[i]);
  }
  free(cg->blist_i);
  free(cg->alist_i);
  free(cg->dlist_i);
  free(fr->x);
  free(fr->f);
  free(fr->list);
  free(fr->head);
  free(fr->map);
  free(mat->bb);
  if(fr->trj_type==0) xdrfile_close(fr->fp);
  else if(fr->trj_type==1) {xdrfile_close(fr->fp); xdrfile_close(fr->fp1);}
}

