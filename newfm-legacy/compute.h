int CalTwoBodyM (TYPE_INFO *const info, CG_PAR *const cg)
{
  return TwoBodyNum(cg->cgtype[info->k],cg->cgtype[info->l],cg->tolcgtype);
}

int CalThreeBodyM (TYPE_INFO *const info, CG_PAR *const cg)
{
  return ThreeBodyNum(cg->cgtype[info->j],cg->cgtype[info->k],cg->cgtype[info->l],cg->tolcgtype);
}

int CalFourBodyM (TYPE_INFO *const info, CG_PAR *const cg)
{
  return FourBodyNum(cg->cgtype[info->i],cg->cgtype[info->j],cg->cgtype[info->k],cg->cgtype[info->l],cg->tolcgtype);
}

void FrameOperationNull(CG_PAR *const cg, MAT_DAT *const mat, FR_DAT *const fr)
{
}

void InitialFunction (CG_PAR *const cg)
{
  void SetZeroSparse(MAT_DAT *const);
  void SetZeroDense(MAT_DAT *const);
  void SetZeroAccumulation(MAT_DAT *const);

  void SparseMatrixInsert(const int, const int, double *const, MAT_DAT *const);
  void DenseMatrixInsert(const int,const int, double *const, MAT_DAT *const);
  void AccumulationMatrixInsert(const int,const int, double *const, MAT_DAT *const);

  void VirialInsertDense(const int, const int, const double, MAT_DAT *const);
  void VirialInsertSparse(const int, const int, const double, MAT_DAT *const);
  void VirialInsertAccumulation(const int, const int, const double, MAT_DAT *const);

  void DenseOperation(CG_PAR *const, MAT_DAT *const);
  void IterOperation(CG_PAR *const, MAT_DAT *const);
  void SparseOperation(CG_PAR *const, MAT_DAT *const);
  void AccumulationOperation(CG_PAR *const, MAT_DAT *const);

  void FMBspline(TYPE_INFO *const);
  void FMLinear(TYPE_INFO *const);

  void FMTwoBody(CG_PAR *const, TYPE_INFO *const, FR_DAT *const, MAT_DAT *const);
  void FMThreeBody(CG_PAR *const, TYPE_INFO *const, FR_DAT *const, MAT_DAT *const);
  void FMFourBody(CG_PAR *const, TYPE_INFO *const, FR_DAT *const, MAT_DAT *const);

  void FillRHS(int, int, MAT_DAT *const, FR_DAT *const);
  void FillRHSAccumulation(int, int, MAT_DAT *const, FR_DAT *const);

  void FillVirial(MAT_DAT *const, FR_DAT *const);
  void FillVirialAccumulation(MAT_DAT *const, FR_DAT *const);

  void TableTwoBody(const int, const int, const int, const double, double *const, double *const, MAT_DAT *const);
  void TableTwoBodyAccumulation(const int, const int, const int, const double, double *const, double *const, MAT_DAT *const);
  void TableThreeBody(const int, const int, const int, const int, double *const, double *const, double *const, double *const, double *const, double *const, MAT_DAT *const);
  void TableThreeBodyAccumulation(const int, const int, const int, const int, double *const, double *const, double *const, double *const, double *const, double *const, MAT_DAT *const);
  void TableFourBody(const int, const int, const int, const int, double *const, const double, double *const, double *const, double *const, double *const, MAT_DAT *const);
  void TableFourBodyAccumulation(const int, const int, const int, const int, double *const, const double, double *const, double *const, double *const, double *const, MAT_DAT *const);


  if(cg->be_sparse==0)
  {
    cg->set_zero=SetZeroDense;
    cg->mat_insert=DenseMatrixInsert;
    if(cg->iter==0) cg->mat_op=DenseOperation;
    else if(cg->iter==1) cg->mat_op=IterOperation;
    if(cg->pres_con==1) cg->vir_insert=VirialInsertDense;
    cg->fill_rhs=FillRHS;
    cg->fill_virial=FillVirial;
    cg->table_twobody=TableTwoBody;
    cg->table_threebody=TableThreeBody;
    cg->table_fourbody=TableFourBody;
  }
  else if(cg->be_sparse==1)
  {
    cg->set_zero=SetZeroSparse;
    cg->mat_insert=SparseMatrixInsert;
    cg->mat_op=SparseOperation;
    if(cg->pres_con==1) cg->vir_insert=VirialInsertSparse;
    cg->fill_rhs=FillRHS;
    cg->fill_virial=FillVirial;
    cg->table_twobody=TableTwoBody;
    cg->table_threebody=TableThreeBody;
    cg->table_fourbody=TableFourBody;
  }
  else if(cg->be_sparse==2)
  {
    cg->set_zero=SetZeroAccumulation;
    cg->mat_insert=AccumulationMatrixInsert;
    cg->mat_op=AccumulationOperation;
    if(cg->pres_con==1) cg->vir_insert=VirialInsertAccumulation;
    cg->fill_rhs=FillRHSAccumulation;
    cg->fill_virial=FillVirialAccumulation;
    cg->table_twobody=TableTwoBodyAccumulation;
    cg->table_threebody=TableThreeBodyAccumulation;
    cg->table_fourbody=TableFourBodyAccumulation;

  }


  if(cg->ba_type==0) cg->fm_basis=FMBspline;
  else if(cg->ba_type==1) cg->fm_basis=FMLinear;

  cg->fr_op=FrameOperationNull;

  cg->fm_nb=FMTwoBody;
  cg->fm_bon=FMTwoBody;
  if(cg->ang_type==0) cg->fm_ang=FMTwoBody;
  else cg->fm_ang=FMThreeBody;
  if(cg->dih_type==0) cg->fm_dih=FMTwoBody;
  else cg->fm_dih=FMFourBody;


  void FMThreeBody(CG_PAR *const, TYPE_INFO *const, FR_DAT *const, MAT_DAT *const);
  void FMThreeBodyR(CG_PAR *const, TYPE_INFO *const, FR_DAT *const, MAT_DAT *const);  
  void FMThreeBodyL(CG_PAR *const, TYPE_INFO *const, FR_DAT *const, MAT_DAT *const);

  if(cg->use_3b==1) cg->fm_threebody=FMThreeBody;
  else if(cg->use_3b==2) 
  {
    if(cg->ba_type!=0) {printf("Three body with distance term can be only used with B-splines!\n"); exit(1);}
    cg->fm_threebody=FMThreeBodyR;
  }
  else if(cg->use_3b==3) cg->fm_threebody=FMThreeBodyL;
  
}
  



void InitialBsplineOne(const int tol, const int tol1, int *const grid_n, int *const num, double *const rmin, double *const rmax, const int n_k, gsl_bspline_workspace ***const w_bs, gsl_vector **const b_bs)
{
  int i,tn=0;
  int n_br,grid;

  *w_bs=(gsl_bspline_workspace **)MemAlloc(tol1*sizeof(gsl_bspline_workspace *));
  *b_bs=gsl_vector_alloc(n_k);

  for(i=0;i<tol;i++)
  {
    if(num[i]>0)
    {
      grid=grid_n[tn+1]-grid_n[tn];
      n_br=grid-n_k+2;
      (*w_bs)[tn]=gsl_bspline_alloc(n_k,n_br);
      gsl_bspline_knots_uniform(rmin[i],rmax[i],(*w_bs)[tn]);
      tn++;
    }
  }
}
      

void InitialBsplineSpace(CG_PAR *const cg)
{
  InitialBsplineOne(cg->tollj,cg->tollj1,cg->grid_n,cg->ljnum,cg->rmin,cg->rmax,cg->n_k,&cg->w_bs,&cg->b_bs);
  InitialBsplineOne(cg->tolb,cg->tolb1,cg->grid_bn,cg->bonnum,cg->rmin_b,cg->rmax_b,cg->nb_k,&cg->wb_bs,&cg->bb_bs);
  InitialBsplineOne(cg->tola,cg->tola1,cg->grid_an,cg->angnum,cg->rmin_a,cg->rmax_a,cg->na_k,&cg->wa_bs,&cg->ba_bs);
  InitialBsplineOne(cg->told,cg->told1,cg->grid_dn,cg->dihnum,cg->rmin_d,cg->rmax_d,cg->nd_k,&cg->wd_bs,&cg->bd_bs);

  int i,n_br,grid;
  if(cg->use_3b==1 || cg->use_3b==2)
  {
    n_br=floor(180.0/cg->space_3+0.5)+1;
    cg->b3_bs=gsl_vector_alloc(cg->n3_k);
    cg->w3_bs=(gsl_bspline_workspace **)MemAlloc(cg->tol3*sizeof(gsl_bspline_workspace *));
    for(i=0;i<cg->tol3;i++)
    {
      cg->w3_bs[i]=gsl_bspline_alloc(cg->n3_k,n_br);
      gsl_bspline_knots_uniform(cg->rmin_3[i],cg->rmax_3[i],cg->w3_bs[i]);
    }
    if(cg->use_3b==2) 
    {
      cg->b3_bs1=gsl_matrix_alloc(cg->n3_k,2);
      cg->w3_bs1=gsl_bspline_deriv_alloc(cg->n3_k);
    }
  }
}

  
void SetPBC (const int l, rvec *frx, matrix box)
{ 
  int i;
  for(i=0;i<3;i++)
  {
    if(frx[l][i]<0) frx[l][i]+=box[i][i];
    else if(frx[l][i]>=box[i][i]) frx[l][i]-=box[i][i];
  }
}

void InitialType(const enum INTER_TYPE type, TYPE_INFO *const info, CG_PAR *const cg)
{
  int  CalTwoBodyM (TYPE_INFO *const, CG_PAR *const);
  int  CalThreeBodyM (TYPE_INFO *const, CG_PAR *const);
  int  CalFourBodyM (TYPE_INFO *const, CG_PAR *const);

  switch(type)
  {
    case non : 
    {
      info->grid_base=0;
      info->grid=cg->grid_n;
      info->rmin=cg->rmin;
      info->rmax=cg->rmax;
      info->num=cg->ljnum;
      info->bs=cg->b_bs;
      info->w=cg->w_bs;
      info->space=cg->space;
      info->n_k=cg->n_k;
      info->space_es=cg->space_es;
      info->es_co=cg->es_co;
      info->cal_m=CalTwoBodyM;
      info->m_size=cg->tollj;
      info->fm_nbody=cg->fm_nb;
      info->cutoff2=cg->cutoff2;
      break;
    }
    case bon :
    {
      info->grid_base=cg->grid_n[cg->tollj1];
      info->grid=cg->grid_bn;
      info->rmin=cg->rmin_b;
      info->rmax=cg->rmax_b;
      info->num=cg->bonnum;
      info->bs=cg->bb_bs;
      info->w=cg->wb_bs;
      info->space=cg->space_b;
      info->n_k=cg->nb_k;      
      info->space_es=cg->spaceb_es;
      info->es_co=cg->esb_co;
      info->cal_m=CalTwoBodyM;
      info->m_array=cg->bon_m;
      info->m_size=cg->tolb;
      info->fm_nbody=cg->fm_bon;
      info->cutoff2=1.0e6;
      break;
    }
    case ang :
    {
      info->grid_base=cg->grid_n[cg->tollj1]+cg->grid_bn[cg->tolb1];
      info->grid=cg->grid_an;
      info->rmin=cg->rmin_a;
      info->rmax=cg->rmax_a;
      info->num=cg->angnum;
      info->bs=cg->ba_bs;
      info->w=cg->wa_bs;
      info->space=cg->space_a;
      info->n_k=cg->na_k;
      info->space_es=cg->spacea_es;
      info->es_co=cg->esa_co;
      info->cal_m=CalThreeBodyM;
      info->m_array=cg->ang_m;
      info->m_size=cg->tola;
      info->fm_nbody=cg->fm_ang;
      info->cutoff2=1.0e6;
      break;
    }
    case dih :
    {
      info->grid_base=cg->grid_n[cg->tollj1]+cg->grid_bn[cg->tolb1]+cg->grid_an[cg->tola1];
      info->grid=cg->grid_dn;
      info->rmin=cg->rmin_d;
      info->rmax=cg->rmax_d;
      info->num=cg->dihnum;
      info->bs=cg->bd_bs;
      info->w=cg->wd_bs;
      info->space=cg->space_d;
      info->n_k=cg->nd_k;
      info->space_es=cg->spaced_es;
      info->es_co=cg->esd_co;
      info->cal_m=CalFourBodyM;
      info->m_array=cg->dih_m;
      info->m_size=cg->told;
      info->fm_nbody=cg->fm_dih;
      info->cutoff2=1.0e6;
      break;
    }
    case t_b :
    {
      info->grid_base=cg->grid_n[cg->tollj1]+cg->grid_bn[cg->tolb1]+cg->grid_an[cg->tola1]+cg->grid_dn[cg->told1];
      info->rmin=cg->rmin_3;
      info->rmax=cg->rmax_3;
      info->num=cg->tbnum;
      info->bs=cg->b3_bs;
      info->w=cg->w3_bs;
      info->grid=cg->grid_3n;
      info->space=cg->space_3;
      info->n_k=cg->n3_k;
      info->cal_m=CalThreeBodyM;
      info->m_array=cg->tb_m;
      info->m_size=cg->tol3;
      //info->cutoff2=cg->cutoff_3;
      break;
    }


    default : {printf("Interaction type not found!\n"); exit(1);}
  }  

  if(cg->ba_type==0) info->n_coef=info->n_k;
  else info->n_coef=2;
}
    
double CalDis (const int k, const int l, FR_DAT *const fr, double *const rx)
{
  int i;
  for(i=0;i<3;i++)
  {
    rx[i]=fr->x[l][i]-fr->x[k][i];
    if(rx[i]>fr->box_dim2[i]) rx[i]=rx[i]-fr->box[i][i];
    else if(rx[i]<-fr->box_dim2[i]) rx[i]=rx[i]+fr->box[i][i];
  }
  double rr2=0.0;
  for(i=0;i<3;i++) rr2+=rx[i]*rx[i];  
  return rr2;
}

void Sub (const int k, const int l, FR_DAT *const fr, double *const rx)
{
  int i;
  for(i=0;i<3;i++)
  {
    rx[i]=fr->x[l][i]-fr->x[k][i];
    if(rx[i]>fr->box_dim2[i]) rx[i]=rx[i]-fr->box[i][i];
    else if(rx[i]<-fr->box_dim2[i]) rx[i]=rx[i]+fr->box[i][i];
  }
}

void CrossProduct(const double *a, const double *b, double *const c)
{
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[2]*b[0]-a[0]*b[2];
  c[2]=a[0]*b[1]-a[1]*b[0];
}

double DotProduct(const double *a, const double *b)
{
  double t=0.0;
  int i;
  for(i=0;i<3;i++) t+=a[i]*b[i];
  return t;
}

double CalAngleCos  (const int j, const int k, const int l, FR_DAT *const fr, double *const rx1, double *const rx2, double *const rr1, double *const rr2)
{
  *rr1=sqrt(CalDis (j,k,fr,rx1));
  *rr2=sqrt(CalDis(j,l,fr,rx2));

  double tx=(rx1[0]*rx2[0]+rx1[1]*rx2[1]+rx1[2]*rx2[2])/((*rr1)*(*rr2));
  double max=1.0-VERYSMALL_F;
  double min=-1.0+VERYSMALL_F;
  if(tx>max) tx=max;
  else if(tx<min) tx=min;
  return tx;
}

double CalDihedral (const int i, const int j, const int k, const int l, FR_DAT *const fr, double *const rx1, double *const rx2, double *const rx3, double *const pb, double *const pc, double *const rpb1, double *const rpc1, double *const pbpc, double *const s)
{
  Sub(k,i,fr,rx1);
  Sub(i,j,fr,rx2);
  Sub(j,l,fr,rx3);

  double rrbc=1./sqrt(DotProduct(rx2,rx2));

  //dihedral vectors
  CrossProduct(rx1,rx2,pb);
  CrossProduct(rx2,rx3,pc);

  double pb2=DotProduct(pb,pb);
  *rpb1=1./sqrt(pb2);

  double pc2=DotProduct(pc,pc);
  *rpc1=1./sqrt(pc2);

  *pbpc=DotProduct(pb,pc);
  double c=(*pbpc)*(*rpb1)*(*rpc1);
  *s=(rx2[0]*(pc[1]*pb[2]-pc[2]*pb[1])+rx2[1]*(pb[0]*pc[2]-pb[2]*pc[0])+ rx2[2]*(pc[0]*pb[1]-pc[1]*pb[0]))*((*rpb1)*(*rpc1)*rrbc);
  if ((*s)<0. && (*s)>-VERYSMALL_F) (*s)=-VERYSMALL_F;
  if ((*s)>0. && (*s)<VERYSMALL_F) (*s)=VERYSMALL_F;
  return atan2((*s),c);
}

int CalMesh (TYPE_INFO *const info)
{
  double tx;

  info->rr0=info->rr-info->rmin[info->mi];
  tx=info->rmax[info->mi]-info->rmin[info->mi]-VERYSMALL;
  if(info->rr0<0.0) info->rr0=0.0;
  else if(info->rr0>tx) info->rr0=tx;
  return (int)(info->rr0/info->space);
}

void FMBspline(TYPE_INFO *const info)
{
  info->i_mesh=CalMesh(info);
  int num=info->numi-1;
  size_t junk;
  gsl_bspline_eval_nonzero(info->rr0+info->rmin[info->mi],info->bs,&junk,&junk,info->w[num]);

  int i;
  for(i=0;i<info->n_k;i++)
  {    
    info->coef[i]=gsl_vector_get(info->bs, i);
  }
}

void FMLinear(TYPE_INFO *const info)
{
  info->i_mesh=CalMesh(info);

  int i;

  info->coef[1]=(info->rr0-info->space*info->i_mesh)/info->space;
  info->coef[0]=1.0-info->coef[1];
}
 
void FMExt(TYPE_INFO *const info)
{
  int i_mesh;
  double tr;
  info->rr0=info->rr-info->rmin[info->mi];
  tr=info->rmax[info->mi]-info->rmin[info->mi]-VERYSMALL;
  if(info->rr0<0.0) info->rr0=0.0;
  else if(info->rr0>tr) info->rr0=tr;
  i_mesh=(int)(info->rr0/info->space_es);

  int i;
  info->coef[1]=(info->rr0-info->space_es*i_mesh)/info->space_es;
  info->coef[0]=1.0-info->coef[1];

  info->coef[0]*=info->es_co[-info->numi-1][i_mesh];
  info->coef[1]*=info->es_co[-info->numi-1][i_mesh+1];

}


void FMTwoBody(CG_PAR *const cg, TYPE_INFO *const info, FR_DAT *const fr, MAT_DAT *const mat)
{
  double rr2,rx[3];
  rr2=CalDis(info->k,info->l,fr,rx);
  if(rr2>info->cutoff2) return;
  info->rr=sqrt(rr2);
//test
  if(info->rr<info->rmin[info->mi]) return;
  else if (info->rr>info->rmax[info->mi]) return;

  double ft[3];
  int i,j,tn;
  double tx[3];

  for(i=0;i<3;i++) ft[i]=rx[i]/info->rr;
  int n_tmp;
  int m_tmp=info->k+info->m_shift;
  int m_tmp1=info->l+info->m_shift;


  if(info->numi>0)
  {
    (*cg->fm_basis)(info,mat);
    n_tmp=info->grid_base+info->grid[info->numi-1]+info->i_mesh;
    for(i=0;i<info->n_coef;i++)
    {
      tn=n_tmp+i;
      for(j=0;j<3;j++) tx[j]=info->coef[i]*ft[j];
      (*cg->mat_insert)(m_tmp1,tn,tx,mat);
      for(j=0;j<3;j++) tx[j]=-tx[j];
      (*cg->mat_insert)(m_tmp,tn,tx,mat);

      if(mat->mp>0) (*cg->vir_insert)(info->inner,tn,info->coef[i]*info->rr,mat);
    }  
  }
  else
  {
    FMExt(info);
    (*cg->table_twobody)(m_tmp,m_tmp1,info->inner,info->rr,info->coef,ft,mat);
  }

}

void FMThreeBody(CG_PAR *const cg, TYPE_INFO *const info, FR_DAT *const fr, MAT_DAT *const mat)
{
  double rx1[3],rx2[3],rr1,rr2;
  double cos_theta=CalAngleCos(info->j,info->k,info->l,fr,rx1,rx2,&rr1,&rr2);
  if(rr1>info->cutoff2 || rr2>info->cutoff2) return;
  double theta=acos(cos_theta);
  info->rr=theta*R2D;
  
  int i,j,tn;
  
  int n_tmp;
  int m_tmp1=info->k+info->m_shift;
  int m_tmp2=info->l+info->m_shift;
  int m_tmp=info->j+info->m_shift;
  
  double ty1[3],ty2[3],ty[3];
  double rr_12_1,rr_11c,rr_22c,vir;
  double sin_theta=sin(theta);
  rr_12_1=1.0/(rr1*rr2*sin_theta);
  rr_11c=cos_theta/(rr1*rr1*sin_theta);
  rr_22c=cos_theta/(rr2*rr2*sin_theta);
  for(i=0;i<3;i++)
  {
    ty1[i]=-rx2[i]*rr_12_1+rr_11c*rx1[i];
    ty2[i]=-rx1[i]*rr_12_1+rr_22c*rx2[i];
    ty[i]=-ty1[i]-ty2[i];
  }

  double tx[3],tx1[3],tx2[3];

  if(info->numi>0)
  {
    (*cg->fm_basis)(info);
    n_tmp=info->grid_base+info->grid[info->numi-1]+info->i_mesh;

    for(i=0;i<info->n_coef;i++)
    {
      tn=n_tmp+i;
      for(j=0;j<3;j++) tx1[j]=ty1[j]*info->coef[i];
      (*cg->mat_insert)(m_tmp1,tn,tx1,mat);  
      for(j=0;j<3;j++) tx2[j]=ty2[j]*info->coef[i];
      (*cg->mat_insert)(m_tmp2,tn,tx2,mat);
      for(j=0;j<3;j++) tx[j]=ty[j]*info->coef[i];
      (*cg->mat_insert)(m_tmp,tn,tx,mat);
/*
      if(mat->mp>0)
      {
        vir=0.0;
        for(j=0;j<3;j++) vir+=(tx1[i]*rx1[i]+tx2[i]*rx2[i]);        
        (*cg->vir_insert)(info->inner,tn,vir,mat);
      }
*/
    }
  }
  else
  {
    FMExt(info);
    (*cg->table_threebody)(m_tmp,m_tmp1,m_tmp2,info->inner,info->coef,ty,ty1,ty2,rx1,rx2,mat);

  }   
}

void FMThreeBodyR(CG_PAR *const cg, TYPE_INFO *const info, FR_DAT *const fr, MAT_DAT *const mat)
{
  double rx1[3],rx2[3],rr1,rr2;
  double cos_theta=CalAngleCos(info->j,info->k,info->l,fr,rx1,rx2,&rr1,&rr2);
  if(rr1>info->cutoff2 || rr2>info->cutoff2) return;

  double theta=acos(cos_theta);
  info->rr=theta*R2D;
//printf("%lf\n",info->rr); 


  int i,j,tn;
  int n_tmp;
  int m_tmp1=info->k+info->m_shift;
  int m_tmp2=info->l+info->m_shift;
  int m_tmp=info->j+info->m_shift;


  double ty1[3],ty2[3],ty[3];
  double rr_12_1,rr_11c,rr_22c,vir;
  double sin_theta=sin(theta);
  rr_12_1=1.0/(rr1*rr2*sin_theta);
  rr_11c=cos_theta/(rr1*rr1*sin_theta);
  rr_22c=cos_theta/(rr2*rr2*sin_theta);
  for(i=0;i<3;i++)
  {
    ty1[i]=-rx2[i]*rr_12_1+rr_11c*rx1[i];
    ty2[i]=-rx1[i]*rr_12_1+rr_22c*rx2[i];
    ty[i]=-ty1[i]-ty2[i];
  }

  double tx[3],tx1[3],tx2[3];

  info->i_mesh=CalMesh(info);
  int num=info->numi-1;
  size_t junk;

  double rr=info->rr0+info->rmin[info->mi];
  gsl_bspline_eval_nonzero(rr,info->bs,&junk,&junk,info->w[num]);

  for(i=0;i<info->n_k;i++)
  {
    info->coef[i]=gsl_vector_get(info->bs, i);
  }

  gsl_bspline_deriv_eval_nonzero(rr,(size_t)1,cg->b3_bs1,&junk,&junk,info->w[num],cg->w3_bs1);
  
  double coef1[100];
  for(i=0;i<info->n_coef;i++)
  {
    coef1[i]=-gsl_matrix_get(cg->b3_bs1,i,1);
  }


  double s1,s2,s1_1,s2_1,rt1,rt2;
  rt1=rr1-info->cutoff2;
  rt2=rr2-info->cutoff2;
  s1=exp(cg->gamma/rt1);
  s2=exp(cg->gamma/rt2);
  s1_1=cg->gamma/(rt1*rt1)*s1;
  s2_1=cg->gamma/(rt2*rt2)*s2;
  double c1,c2,c3;
  c1=s1*s2;
  c2=s2*s1_1;
  c3=s1*s2_1;


  n_tmp=info->grid_base+info->grid[info->numi-1]+info->i_mesh;

  double tt;
  c1*=R2D;
  for(i=0;i<info->n_coef;i++)
  {
    tn=n_tmp+i;

    tt=coef1[i]*c1;
    for(j=0;j<3;j++) tx1[j]=ty1[j]*tt;
    for(j=0;j<3;j++) tx2[j]=ty2[j]*tt;

    tt=info->coef[i]*c2;
    for(j=0;j<3;j++) tx1[j]+=rx1[j]/rr1*tt;
    (*cg->mat_insert)(m_tmp1,tn,tx1,mat);
    tt=info->coef[i]*c3;
    for(j=0;j<3;j++) tx2[j]+=rx2[j]/rr2*tt;
    (*cg->mat_insert)(m_tmp2,tn,tx2,mat);
    for(j=0;j<3;j++) tx[j]=-(tx1[j]+tx2[j]);
    (*cg->mat_insert)(m_tmp,tn,tx,mat);
  }

}

void FMThreeBodyL(CG_PAR *const cg, TYPE_INFO *const info, FR_DAT *const fr, MAT_DAT *const mat)
{
  double rx1[3],rx2[3],rr1,rr2;
  double cos_theta=CalAngleCos(info->j,info->k,info->l,fr,rx1,rx2,&rr1,&rr2);
  if(rr1>info->cutoff2 || rr2>info->cutoff2) return;

  double theta=acos(cos_theta);
  info->rr=theta*R2D;

  int i,j,tn;
  int n_tmp;
  int m_tmp1=info->k+info->m_shift;
  int m_tmp2=info->l+info->m_shift;
  int m_tmp=info->j+info->m_shift;


  double ty1[3],ty2[3],ty[3];
  double rr_12_1,rr_11c,rr_22c,vir;
  double sin_theta=sin(theta);
  rr_12_1=1.0/(rr1*rr2*sin_theta);
  rr_11c=cos_theta/(rr1*rr1*sin_theta);
  rr_22c=cos_theta/(rr2*rr2*sin_theta);
  for(i=0;i<3;i++)
  {
    ty1[i]=-rx2[i]*rr_12_1+rr_11c*rx1[i];
    ty2[i]=-rx1[i]*rr_12_1+rr_22c*rx2[i];
    ty[i]=-ty1[i]-ty2[i];
  }

  double tx[3],tx1[3],tx2[3];


  double s1,s2,s1_1,s2_1,rt1,rt2;
  rt1=rr1-info->cutoff2;
  rt2=rr2-info->cutoff2;
  s1=exp(cg->gamma/rt1);
  s2=exp(cg->gamma/rt2);
  s1_1=cg->gamma/(rt1*rt1)*s1;
  s2_1=cg->gamma/(rt2*rt2)*s2;
  double c1,c2,c3;
  c1=s1*s2;
  c2=s2*s1_1;
  c3=s1*s2_1;
 
  double u,u_1;
  u=(cos_theta-cg->cos_theta0)*(cos_theta-cg->cos_theta0)*4.184;
  u_1=2.0*(cos_theta-cg->cos_theta0)*sin_theta*4.184;
  
  n_tmp=info->grid_base+info->grid[info->numi-1];
  double tt;
  //c1*=R2D;

  tn=n_tmp;
//printf("%d\n",tn);
  tt=c1*u_1;
  for(j=0;j<3;j++) tx1[j]=ty1[j]*tt;
  for(j=0;j<3;j++) tx2[j]=ty2[j]*tt;

  tt=c2*u;
  for(j=0;j<3;j++) tx1[j]+=rx1[j]/rr1*tt;
  (*cg->mat_insert)(m_tmp1,tn,tx1,mat);
  tt=c3*u;
  for(j=0;j<3;j++) tx2[j]+=rx2[j]/rr2*tt;
  (*cg->mat_insert)(m_tmp2,tn,tx2,mat);
  for(j=0;j<3;j++) tx[j]=-(tx1[j]+tx2[j]);
  (*cg->mat_insert)(m_tmp,tn,tx,mat);

}



void FMFourBody(CG_PAR *const cg,TYPE_INFO *const info, FR_DAT *const fr, MAT_DAT *const mat)
{
  double rx1[3],rx2[3],rx3[3],pb[3],pc[3],rpb1,rpc1,pbpc,s;
  info->rr=CalDihedral(info->i,info->j,info->k,info->l,fr,rx1,rx2,rx3,pb,pc,&rpb1,&rpc1,&pbpc,&s)*R2D;

  double gamma= rpb1*rpc1*(1./s);

  double cb,cc,w1[3],w2[3],w3[3],w4[3];
  cb=pbpc*rpb1*rpb1; cc=pbpc*rpc1*rpc1;

  w1[0]=(-pc[1]*rx2[2]+pc[2]*rx2[1])-cb*(-pb[1]*rx2[2]+pb[2]*rx2[1]);
  w1[1]=( pc[0]*rx2[2]-pc[2]*rx2[0])-cb*( pb[0]*rx2[2]-pb[2]*rx2[0]);
  w1[2]=(-pc[0]*rx2[1]+pc[1]*rx2[0])-cb*(-pb[0]*rx2[1]+pb[1]*rx2[0]);

  w2[0]=(-pc[1]*rx1[2]+pc[2]*rx1[1])-cb*(-pb[1]*rx1[2]+pb[2]*rx1[1]);
  w2[1]=( pc[0]*rx1[2]-pc[2]*rx1[0])-cb*( pb[0]*rx1[2]-pb[2]*rx1[0]);
  w2[2]=(-pc[0]*rx1[1]+pc[1]*rx1[0])-cb*(-pb[0]*rx1[1]+pb[1]*rx1[0]);

  w3[0]=(-pb[1]*rx3[2]+pb[2]*rx3[1])-cc*(-pc[1]*rx3[2]+pc[2]*rx3[1]);
  w3[1]=( pb[0]*rx3[2]-pb[2]*rx3[0])-cc*( pc[0]*rx3[2]-pc[2]*rx3[0]);
  w3[2]=(-pb[0]*rx3[1]+pb[1]*rx3[0])-cc*(-pc[0]*rx3[1]+pc[1]*rx3[0]);

  w4[0]=(-pb[1]*rx2[2]+pb[2]*rx2[1])-cc*(-pc[1]*rx2[2]+pc[2]*rx2[1]);
  w4[1]=( pb[0]*rx2[2]-pb[2]*rx2[0])-cc*( pc[0]*rx2[2]-pc[2]*rx2[0]);
  w4[2]=(-pb[0]*rx2[1]+pb[1]*rx2[0])-cc*(-pc[0]*rx2[1]+pc[1]*rx2[0]);

  int i,j,k,tn;
  double tx[3],ty;

  int n_tmp;
  int m_tmp=info->k+info->m_shift;
  int m_tmp1=info->l+info->m_shift;
  int m_tmp2=info->i+info->m_shift;
  int m_tmp3=info->j+info->m_shift;
  
  if(info->numi>0)
  {
    (*cg->fm_basis)(info);
    n_tmp=info->grid_base+info->grid[info->numi-1]+info->i_mesh;
    for(i=0;i<info->n_coef;i++)
    {
      ty=gamma*info->coef[i];
      tn=n_tmp+i;      

      for(j=0;j<3;j++) tx[j]=w1[j]*ty;
      (*cg->mat_insert)(m_tmp,tn,tx,mat);

      for(j=0;j<3;j++) tx[j]=w4[j]*ty;
      (*cg->mat_insert)(m_tmp1,tn,tx,mat);

      for(j=0;j<3;j++) tx[j]=(-w1[j]-w2[j]+w3[j])*ty;
      (*cg->mat_insert)(m_tmp2,tn,tx,mat);
   
      for(j=0;j<3;j++) tx[j]=(w2[j]-w3[j]-w4[j])*ty;
      (*cg->mat_insert)(m_tmp3,tn,tx,mat);
    }
  }
  else
  {
    FMExt(info);
    (*cg->table_fourbody)(m_tmp,m_tmp1,m_tmp2,m_tmp3,info->coef,gamma,w1,w2,w3,w4,mat);
  }        
}  


void CalNB(TYPE_INFO *const info, CG_PAR *const cg, MAT_DAT *const mat, FR_DAT *const fr)
{
  int i;
  for(i=0;i<cg->blist_n[info->k];i++)
  {
    if(cg->blist_i[info->k][i]==info->l) return;
  }
  for(i=0;i<cg->alist_n[info->k];i++)
  {
    if(cg->alist_i[info->k][2*i+1]==info->l) return;
  }
  for(i=0;i<cg->dlist_n[info->k];i++)
  {
    if(cg->dlist_i[info->k][3*i+2]==info->l) return;
  }

  info->m=TwoBodyNum(cg->cgtype[info->k],cg->cgtype[info->l],cg->tolcgtype);
  info->mi=info->m;
  info->numi=info->num[info->mi];
  if(info->numi==0) return;

  (*info->fm_nbody)(cg,info,fr,mat);
}    

void CalBonded(TYPE_INFO *const info, CG_PAR *const cg, MAT_DAT *const mat, FR_DAT *const fr)
{
  info->m=(*info->cal_m)(info,cg);

  info->mi=SearchIntTable(info->m_array,info->m,info->m_size);
  info->numi=info->num[info->mi];
  if(info->numi==0) return;

  (*info->fm_nbody)(cg,info,fr,mat);
}

void CalNB_3(TYPE_INFO *const info, CG_PAR *const cg, MAT_DAT *const mat, FR_DAT *const fr)
{   
  void FMThreeBody(CG_PAR *const, TYPE_INFO *const, FR_DAT *const, MAT_DAT *const);
  int i;

/*
      for(i=0;i<cg->blist_n[info->k];i++)
      { 
        if(cg->blist_i[info->k][i]==info->l) return;
      } 
      for(i=0;i<cg->alist_n[info->k];i++)
      {
        if(cg->alist_i[info->k][2*i+1]==info->l) return;
      }
      for(i=0;i<cg->dlist_n[info->k];i++)
      {
        if(cg->dlist_i[info->k][3*i+2]==info->l) return;
      }
  */  

    for(i=0;i<cg->blist_n[info->l];i++)
    {
      if(cg->blist_i[info->l][i]==info->j) return;
    }
/*
    for(i=0;i<cg->alist_n[info->l];i++)
    {
      if(cg->alist_i[info->l][2*i+1]==info->j) return;
    }
    for(i=0;i<cg->dlist_n[info->l];i++)
    {
      if(cg->dlist_i[info->l][3*i+2]==info->j) return;
    }
*/
    for(i=0;i<cg->blist_n[info->j];i++)
    {
      if(cg->blist_i[info->j][i]==info->k) return;
    }
/*
    for(i=0;i<cg->alist_n[info->j];i++)
    {
      if(cg->alist_i[info->j][2*i+1]==info->k) return;
    }
    for(i=0;i<cg->dlist_n[info->j];i++)
    {
      if(cg->dlist_i[info->j][3*i+2]==info->k) return;
    }
*/

  info->m=(*info->cal_m)(info,cg);
  info->mi=SearchIntTable(info->m_array,info->m,info->m_size);
  if(info->mi==-1) return;

  info->numi=info->num[info->mi];
  if(info->numi==0) return;
  
  info->cutoff2=cg->cutoff_3i[info->mi];
  cg->cos_theta0=cg->cos_theta0i[info->mi];
  //printf("%lf %lf\n",info->cutoff2,cg->cos_theta0); 
  (*cg->fm_threebody)(cg,info,fr,mat);
}



void FreeFMSpace(CG_PAR *const cg, MAT_DAT *const mat, FR_DAT *const fr)
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
  if(cg->tollj2>0) {for(i=0;i<cg->tollj2;i++) free(cg->es_co[i]); free(cg->es_co);}
  if(cg->tolb2>0) {for(i=0;i<cg->tolb2;i++) free(cg->esb_co[i]); free(cg->esb_co);}
  if(cg->tola2>0) {for(i=0;i<cg->tola2;i++) free(cg->esa_co[i]); free(cg->esa_co);}
  if(cg->told2>0) {for(i=0;i<cg->told2;i++) free(cg->esd_co[i]); free(cg->esd_co);}
  free(fr->x);
  free(fr->f);
  free(fr->list);
  free(fr->head);
  free(fr->map);
  if(cg->use_3b>0)
  {
    free(fr->list_3);
    free(fr->head_3);
    free(fr->map_3);
  }
  if(mat->mp>0) free(fr->b_p);
  if(cg->be_sparse==0 || cg->be_sparse==1)
  {
    if(cg->be_sparse==0 && mat->mp>0) free(mat->aa);
    free(mat->bb);
  }
  if(cg->be_sparse==1)
  {
    free(mat->mat_head);
    free(mat->iw);
    free(mat->xx_b);
    free(mat->v);
    free(mat->w);
    free(mat->se);
  }
  if(fr->trj_type==0) xdrfile_close(fr->fp);
  else if(fr->trj_type==1) {xdrfile_close(fr->fp); xdrfile_close(fr->fp1);}
  if(cg->be_sparse==2) {free(mat->work); free(mat->tau);}
}
  

void DoFM(CG_PAR *const cg, MAT_DAT *const mat, FR_DAT *const fr)
{
  int i,j,k,l;

  for(i=0;i<cg->start_fr-1;i++)
  {
    if((*cg->read_next)(fr)==0) {printf("\nCan not read a frame that needs to be skipped!\n"); exit(1);}
  }
  
  int n_step;
  if(cg->be_sparse==0) {n_step=cg->fr_n; cg->n_block=1;}
  else n_step=cg->fr_n/cg->n_block;

  int m_shift;
  int tn,tn1;
  int m33=mat->mv*3;
  int kk,nei,nei_3,m,ll,ll_3;
  int read_stat=1;
  TYPE_INFO info;
  enum INTER_TYPE type;
  mat->row_shift=0; //used only for accumulation QR method

  for(mat->outer=0;mat->outer<n_step;mat->outer++)
  {
    (*cg->set_zero)(mat);
    (*cg->fill_virial)(mat,fr);

    for(info.inner=0;info.inner<cg->n_block;info.inner++)
    {
      if(read_stat==0) {printf("\nCan not read frame %d for force matching!\n",mat->outer*cg->n_block+info.inner); exit(1);}
      info.m_shift=info.inner*cg->tolcgn; //shift row number after each frame within one block

      //loop to obtain positions and forces of each CG site
      for(l=0;l<cg->tolcgn;l++)
      {
        //enforce pbc
        SetPBC(l,fr->x,fr->box);
        (*cg->fill_rhs)(info.m_shift,l,mat,fr);
      }

      LinkedList(fr);
      type=non;
      InitialType(type,&info,cg);

      //loop for nonebonded interactions, using cell linked-list
      for(kk=0;kk<fr->ncell;kk++)
      {
        info.k=fr->head[kk];
        while (info.k>=0)
        {
          //printf("%d\n",k);
          info.l=fr->list[info.k];
          while (info.l>=0)
          { 
            CalNB(&info,cg,mat,fr);                      
            info.l=fr->list[info.l];
          }
          //do the above the 2nd time for neiboring cells
          for(nei=0;nei<13;nei++)
          {
            ll=fr->map[13*kk+nei];
            info.l=fr->head[ll]; 
            while (info.l>=0)
            {
              CalNB(&info,cg,mat,fr);
              info.l=fr->list[info.l]; 
            }
          }
          info.k=fr->list[info.k];
        }
      }

      type=bon;
      InitialType(type,&info,cg);
      for(info.k=0;info.k<cg->tolcgn;info.k++)
      {
        for(kk=0;kk<cg->blist_n[info.k];kk++)
        {
          info.l=cg->blist_i[info.k][kk];
          if(info.k<info.l) CalBonded(&info,cg,mat,fr);        
        }  
      }

      type=ang;
      InitialType(type,&info,cg);
      for(info.k=0;info.k<cg->tolcgn;info.k++)
      {
        for(kk=0;kk<cg->alist_n[info.k];kk++)
        {
          info.l=cg->alist_i[info.k][2*kk+1];
          info.j=cg->alist_i[info.k][2*kk];
          if(info.k<info.l) CalBonded(&info,cg,mat,fr);
        }
      }

      type=dih;
      InitialType(type,&info,cg);
      for(info.k=0;info.k<cg->tolcgn;info.k++)
      {
        for(kk=0;kk<cg->dlist_n[info.k];kk++)
        {
          info.l=cg->dlist_i[info.k][3*kk+2];
          info.i=cg->dlist_i[info.k][3*kk];
          info.j=cg->dlist_i[info.k][3*kk+1];
          if(info.k<info.l) CalBonded(&info,cg,mat,fr);
        }
      }

      if(cg->use_3b>0)
      {
        LinkedList_3(fr);

        type=t_b;
        InitialType(type,&info,cg);


        for(kk=0;kk<fr->ncell_3;kk++)
        {
          info.j=fr->head_3[kk];
          while (info.j>=0)
          {
            info.k=fr->head_3[kk];
            while (info.k>=0)
            {
              if(info.j!=info.k)
              {
                //three body
                info.l=fr->list_3[info.k];
                while(info.l>=0)
                {
                 if(info.l!=info.j) CalNB_3(&info,cg,mat,fr);
                 info.l=fr->list_3[info.l];
                }
                for(nei_3=0;nei_3<26;nei_3++)
                {
                  ll_3=fr->map_3[26*kk+nei_3];
                  info.l=fr->head_3[ll_3];
                  while (info.l>=0)
                  {
                    CalNB_3(&info,cg,mat,fr);
                    info.l=fr->list_3[info.l];
                  }
                }
              }
              info.k=fr->list_3[info.k];
            }

            for(nei=0;nei<26;nei++)
            {
              ll=fr->map_3[26*kk+nei];
              info.k=fr->head_3[ll];
              while (info.k>=0)
              {
                //three body
                info.l=fr->list_3[info.k];
                while(info.l>=0)
                {
                  CalNB_3(&info,cg,mat,fr);
                  info.l=fr->list_3[info.l];
                }
                for(nei_3=nei+1;nei_3<26;nei_3++)
                {
                  ll_3=fr->map_3[26*kk+nei_3];
                  info.l=fr->head_3[ll_3];
                  while (info.l>=0)
                  {
                    CalNB_3(&info,cg,mat,fr);
                    info.l=fr->list_3[info.l];
                  }
                }
                info.k=fr->list_3[info.k];
              }
            }
            info.j=fr->list_3[info.j];
          }
        }

      } //end if threebody



      //(*cg->fr_op)(cg,mat,fr);
      read_stat=(*cg->read_next)(fr);

    } //end of the inner loop 


    printf("\r%d frames have been read.",(mat->outer+1)*cg->n_block); fflush(stdout);
    (*cg->mat_op)(cg,mat); 

  } //end of the outer loop
  

}



