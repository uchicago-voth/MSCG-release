void InitialDenseMatrix(CG_PAR *const cg, MAT_DAT *const mat)
{
  mat->nn=cg->grid_n[cg->tollj1]+cg->grid_bn[cg->tolb1]+cg->grid_an[cg->tola1]+cg->grid_dn[cg->told1];
  if(cg->use_3b>0) mat->nn+=cg->grid_3n[cg->tol3];
  mat->mv=cg->tolcgn;
  if(cg->pres_con==0) {mat->mm=cg->tolcgn*3; mat->mp=0;}
  else {mat->mm=cg->tolcgn*3+1; mat->mp=1;}
  if(mat->mm*cg->fr_n<mat->nn) {printf("Need more frames!\n"); exit(1);}

  mat->na=mat->nn;
  mat->ma=mat->mm;
  mat->aa=(double *)MemAlloc(mat->mm*mat->nn*sizeof(double));
  mat->bb=(double *)MemAlloc(mat->mm*sizeof(double));
  mat->gg=(double *)MemAlloc(mat->nn*mat->nn*sizeof(double));
  mat->dd=(double *)MemAlloc(mat->nn*sizeof(double));

  int i;
  for(i=0;i<mat->nn*mat->nn;i++)
  {
    mat->gg[i]=0.0;
  }

  for(i=0;i<mat->nn;i++)
  {
    mat->dd[i]=0.0;
  }

}

void InitialAccumulationMatrix(CG_PAR *const cg, MAT_DAT *const mat)
{
  mat->nn=cg->grid_n[cg->tollj1]+cg->grid_bn[cg->tolb1]+cg->grid_an[cg->tola1]+cg->grid_dn[cg->told1];
  if(cg->use_3b>0) mat->nn+=cg->grid_3n[cg->tol3];

  mat->mv=cg->tolcgn*cg->n_block;
  if(cg->pres_con==0) {mat->mm=mat->mv*3; mat->mp=0;}
  else {mat->mm=mat->mv*3+cg->n_block; mat->mp=cg->n_block;}
  if(mat->mm<mat->nn) {printf("Please increase your block size!\n"); exit(1);}

  mat->na=mat->nn+1;
  mat->ma=mat->mm+mat->na;

  mat->aa=(double *)MemAlloc(mat->ma*mat->na*sizeof(double));
  mat->tau=(double *)MemAlloc(mat->na*sizeof(double));
  mat->dd=(double *)MemAlloc(mat->na*sizeof(double));

  int i;
  for(i=0;i<mat->ma*mat->na;i++) mat->aa[i]=0.0;
}



void InitialSparseMatrix(CG_PAR *const cg, MAT_DAT *const mat)
{
  mat->nn=cg->grid_n[cg->tollj1]+cg->grid_bn[cg->tolb1]+cg->grid_an[cg->tola1]+cg->grid_dn[cg->told1];
  if(cg->use_3b>0) mat->nn+=cg->grid_3n[cg->tol3];

  mat->mv=cg->tolcgn*cg->n_block;
  if(cg->pres_con==0) {mat->mm=mat->mv*3; mat->mp=0;}
  else {mat->mm=mat->mv*3+cg->n_block; mat->mp=cg->n_block;}
  if(mat->mm<mat->nn) {printf("Please increase your block size!\n"); exit(1);}
  mat->na=mat->nn;
  mat->ma=mat->ma;

  //for each row we have one linked list
  mat->bb=(double *)MemAlloc(mat->mm*sizeof(double));
  mat->mat_head=(struct MAT_H *)MemAlloc(mat->mv*sizeof(struct MAT_H));
  if(cg->pres_con==1) mat->aa=(double *)MemAlloc(cg->n_block*mat->nn*sizeof(double));

  mat->iw=(int *)malloc((unsigned)(mat->mm+1)*sizeof(int));
  mat->iw[0]=0;

  cg->damp = 0.0;
  cg->mm1=mat->mm+1;
  mat->n_xx=(int *)MemAlloc(mat->nn*sizeof(int));
  mat->xx=(double *)MemAlloc(mat->nn*sizeof(double));
  mat->v=(double *)MemAlloc(mat->nn*sizeof(double));
  mat->w=(double *)MemAlloc(mat->nn*sizeof(double));
  mat->se=(double *)MemAlloc(mat->nn*sizeof(double));
  if(cg->itnlim==0) cg->itnlim = 4*mat->nn; //4*nn; //max number of iterations
  cg->istop = 0;
  cg->itn = 0;
  cg->anorm = 0.0;
  cg->acond = 0.0;
  cg->rnorm = 0.0;
  cg->arnorm = 0.0;
  cg->xnorm = 0.0;

  mat->xx_b=(double *)MemAlloc(mat->nn*sizeof(double));
  //number of sampled blocks and global answer for each unkown
  int i;
  for(i=0;i<mat->nn;i++)
  {
    mat->n_xx[i]=0;
    mat->xx[i]=0.0;
  }

  //precondition temp array
  //reading "Solving lest square problems" by Lawson CL and Hanson RJ, Chapt 25 for details
  mat->h=(double *)MemAlloc(mat->nn*sizeof(double));

}

void SparseMatrixInsert(const int i, const int j, double *const x, MAT_DAT *const mat)
{
  struct MAT_SP *p, *p1, *head, *pt;
  int k;

  head=mat->mat_head[i].h;

  //if the linked list is empty
  if(head==NULL) 
  {
    pt=(struct MAT_SP *)MemAlloc(sizeof(struct MAT_SP));
    mat->mat_head[i].h=pt;    
    pt->col=j;
    for(k=0;k<3;k++) pt->valx[k]=x[k];
    pt->next=NULL;
    mat->mat_head[i].n+=1;
  }
  //insert a none-empty linked list
  else
  {
    p=head;
    p1=NULL;
    while(p!=NULL)
    {     
      //if the element exists
      if(p->col==j) 
      {
        for(k=0;k<3;k++) p->valx[k]+=x[k];
        p1=p;
        break;
      }
      //if the new element should be the first in the new liked list
      else if(p->col>j) 
      {
        if(p1==NULL)
        {
          pt=(struct MAT_SP *)MemAlloc(sizeof(struct MAT_SP));
          mat->mat_head[i].h=pt;
          pt->col=j;
          for(k=0;k<3;k++) pt->valx[k]=x[k];
          pt->next=p;
        }
        //general case
        else
        {
          pt=(struct MAT_SP *)MemAlloc(sizeof(struct MAT_SP));
          p1->next=pt;
          pt->col=j;
          for(k=0;k<3;k++) pt->valx[k]=x[k]; 
          pt->next=p;
        }
        mat->mat_head[i].n+=1;
        p1=p;
        break;
      }        
      p1=p;
      p=p1->next;
    } 
    //if the new element should be the last
    if(p1->col<j)
    {
      pt=(struct MAT_SP *)MemAlloc(sizeof(struct MAT_SP));
      p1->next=pt;
      pt->col=j;
      for(k=0;k<3;k++) pt->valx[k]=x[k];
      pt->next=NULL;
      mat->mat_head[i].n+=1;
    }
  }    

}

void DenseMatrixInsert(const int i, const int j, double *const x, MAT_DAT *const mat)
{
  int m=3*i;
  int k;
  for(k=0;k<3;k++)
  {
    mat->aa[j*mat->mm+m+k]+=x[k];
  }
}

void AccumulationMatrixInsert(const int i, const int j, double *const x, MAT_DAT *const mat)
{
  int m=3*i;
  int k;
  for(k=0;k<3;k++)
  {
    mat->aa[j*mat->ma+m+mat->row_shift+k]+=x[k];
  }
}

 
void SetZeroDense(MAT_DAT *const mat)
{
  int i;
  for(i=0;i<mat->mm*mat->nn;i++)
  {
    mat->aa[i]=0.0;
  }
}

void SetZeroSparse(MAT_DAT *const mat)
{
  int k;
  for(k=0;k<mat->mv;k++) //initializing
  {
    mat->mat_head[k].h=NULL;
    mat->mat_head[k].n=0;
  }
  for(k=0;k<mat->mp*mat->nn;k++)
  {
    mat->aa[k]=0.0;
  }
}

void SetZeroAccumulation(MAT_DAT *const mat)
{
  int k,l;
  for(k=0;k<mat->na;k++)
  {
    for(l=0;l<k;l++)
    {
      mat->aa[l*mat->ma+k]=0.0;
    }
  }

  for(k=mat->na;k<mat->ma;k++)
  {
    for(l=0;l<mat->na;l++)
    {
      mat->aa[l*mat->ma+k]=0.0;
    }
  } 
}


void FillVirial(MAT_DAT *const mat, FR_DAT *const fr)
{
  int k;
  for(k=0;k<mat->mp;k++)
  {
    mat->bb[mat->mv*3+k]=fr->b_p[mat->outer*mat->mp+k];
  }
}

void FillVirialAccumulation(MAT_DAT *const mat, FR_DAT *const fr)
{
  int k;
  for(k=0;k<mat->mp;k++)
  {
    mat->aa[mat->nn*mat->ma+mat->mv*3+mat->row_shift+k]=fr->b_p[mat->outer*mat->mp+k];
  }
} 


void FillRHS(int shift_i, int site_i, MAT_DAT *const mat, FR_DAT *const fr)
{
  int i,tn;
  tn=3*(site_i+shift_i);
  for(i=0;i<3;i++) mat->bb[tn+i]=fr->f[site_i][i];
}

void FillRHSAccumulation(int shift_i, int site_i, MAT_DAT *const mat, FR_DAT *const fr)
{
  int i,tn;
  tn=mat->nn*mat->ma+3*(site_i+shift_i)+mat->row_shift;
  for(i=0;i<3;i++) mat->aa[tn+i]=fr->f[site_i][i];
}

void VirialInsertDense(const int m, const int n, const double x, MAT_DAT *const mat)
{
  mat->aa[mat->mv*3+n*mat->mm]+=x;
}
 
void VirialInsertSparse(const int m, const int n, const double x, MAT_DAT *const mat)
{
  mat->aa[m+n*mat->mp]+=x;
}

void VirialInsertAccumulation(const int m, const int n, const double x, MAT_DAT *const mat)
{
  mat->aa[mat->mv*3+mat->row_shift+n*mat->mm]+=x;
}


void TableTwoBody(const int m_tmp, const int m_tmp1, const int inner, const double rr, double *const coef, double *const ft, MAT_DAT *const mat)
{
  int i,j;
  double tx[3];
  for(i=0;i<2;i++)
  {
    for(j=0;j<3;j++)
    {
      tx[j]=coef[i]*ft[j];
      mat->bb[3*m_tmp1+j]-=tx[j];
      mat->bb[3*m_tmp+j]+=tx[j];
    }
    if(mat->mp>0) mat->bb[mat->mv*3+inner]-=coef[i]*rr;
  }
}

void TableTwoBodyAccumulation(const int m_tmp, const int m_tmp1, const int inner, const double rr, double *const coef, double *const ft, MAT_DAT *const mat)
{
  int i,j;
  double tx[3];
  int shift=mat->nn*mat->ma;
  for(i=0;i<2;i++)
  {
    for(j=0;j<3;j++)
    {
      tx[j]=coef[i]*ft[j];
      mat->aa[shift+3*m_tmp1+j]-=tx[j];
      mat->aa[shift+3*m_tmp+j]+=tx[j];
    }
    if(mat->mp>0) mat->aa[shift+mat->mv*3+inner]-=coef[i]*rr;
  }
}

void TableThreeBody(const int m_tmp, const int m_tmp1, const int m_tmp2, const int inner, double *const coef, double *const ty, double *const ty1, double *const ty2, double *const rx1, double *const rx2, MAT_DAT *const mat)
{
  int i,j;
  double tx[3],tx1[3],tx2[3];

  for(i=0;i<2;i++)
  {
    for(j=0;j<3;j++)
    {
      tx1[j]=coef[i]*ty1[j];
      mat->bb[3*m_tmp1+j]-=tx1[j];
      tx2[j]=coef[i]*ty2[j];
      mat->bb[3*m_tmp2+j]-=tx2[j];
      tx[j]=coef[i]*ty[j];
      mat->bb[3*m_tmp+j]-=tx[j];
    }
    if(mat->mp>0) mat->bb[mat->mv*3+inner]-=(tx1[i]*rx1[i]+tx2[i]*rx2[i]);
  }
}

void TableThreeBodyAccumulation(const int m_tmp, const int m_tmp1, const int m_tmp2, const int inner, double *const coef, double *const ty, double *const ty1, double *const ty2, double *const rx1, double *const rx2, MAT_DAT *const mat)
{
  int i,j;
  double tx[3],tx1[3],tx2[3];
  int shift=mat->nn*mat->ma;
  
  for(i=0;i<2;i++)
  {
    for(j=0;j<3;j++)
    {
      tx1[j]=coef[i]*ty1[j];
      mat->aa[shift+3*m_tmp1+j]-=tx1[j];
      tx2[j]=coef[i]*ty2[j];
      mat->aa[shift+3*m_tmp2+j]-=tx2[j];
      tx[j]=coef[i]*ty[j];
      mat->aa[shift+3*m_tmp+j]-=tx[j];
    }
    if(mat->mp>0) mat->aa[shift+mat->mv*3+inner]-=(tx1[i]*rx1[i]+tx2[i]*rx2[i]);
  }
}

void TableFourBody(const int m_tmp, const int m_tmp1, const int m_tmp2, const int m_tmp3, double *const coef, const double gamma, double *const w1, double *const w2, double *const w3, double *const w4,  MAT_DAT *const mat)
{
  int i,j;
  double tx[3],ty;
  for(i=0;i<2;i++)
  {
    ty=gamma*coef[i];
    for(j=0;j<3;j++)
    {
      tx[j]=w1[j]*ty;
      mat->bb[3*m_tmp+j]-=tx[j];

      tx[j]=w4[j]*ty;
      mat->bb[3*m_tmp1+j]-=tx[j];

      tx[j]=(-w1[j]-w2[j]+w3[j])*ty;
      mat->bb[3*m_tmp2+j]-=tx[j];

      tx[j]=(w2[j]-w3[j]-w4[j])*ty;
      mat->bb[3*m_tmp3+j]-=tx[j];
    }
  }
}

void TableFourBodyAccumulation(const int m_tmp, const int m_tmp1, const int m_tmp2, const int m_tmp3, double *const coef, const double gamma, double *const w1, double *const w2, double *const w3, double *const w4,  MAT_DAT *const mat)
{
  int i,j;
  double tx[3],ty;
  int shift=mat->nn*mat->ma;

  for(i=0;i<2;i++)
  {
    ty=gamma*coef[i];
    for(j=0;j<3;j++)
    {
      tx[j]=w1[j]*ty;
      mat->aa[shift+3*m_tmp+j]-=tx[j];

      tx[j]=w4[j]*ty;
      mat->aa[shift+3*m_tmp1+j]-=tx[j];

      tx[j]=(-w1[j]-w2[j]+w3[j])*ty;
      mat->aa[shift+3*m_tmp2+j]-=tx[j];

      tx[j]=(w2[j]-w3[j]-w4[j])*ty;
      mat->aa[shift+3*m_tmp3+j]-=tx[j];
    }
  }
}

void DenseOperation(CG_PAR *const cg, MAT_DAT *const mat)
{

  char uchar='U';
  char tchar='T';
  int onei=1;
  double oned=1.0;
  dsyrk_(&uchar, &tchar, &mat->nn, &mat->mm, &cg->fr_n_1, mat->aa, &mat->mm, &oned, mat->gg, &mat->nn);
  dgemv_(&tchar, &mat->mm, &mat->nn, &cg->fr_n_1, mat->aa, &mat->mm, mat->bb, &onei, &oned, mat->dd, &onei);
}

void IterOperation(CG_PAR *const cg, MAT_DAT *const mat)
{

  char uchar='U';
  char tchar='T';
  int onei=1;
  double oned=1.0;
  //dsyrk_(&uchar, &tchar, &mat->nn, &mat->mm, &oned, mat->aa, &mat->mm, &oned, mat->gg, &mat->nn);
  dgemv_(&tchar, &mat->mm, &mat->nn, &oned, mat->aa, &mat->mm, mat->bb, &onei, &oned, mat->dd, &onei);
}



void AccumulationOperation(CG_PAR *const cg, MAT_DAT *const mat)
{ 
  int info_in;
  if(mat->outer==0)
  {
    mat->work=(double *)MemAlloc(1*sizeof(double));
    mat->lwork=-1;
    dgeqrf_(&mat->mm,&mat->na,mat->aa,&mat->ma,mat->tau,mat->work,&mat->lwork,&info_in);
    mat->lwork=mat->work[0];
    free(mat->work);
    mat->work=(double *)MemAlloc(mat->lwork*sizeof(double));
    dgeqrf_(&mat->mm,&mat->na,mat->aa,&mat->ma,mat->tau,mat->work,&mat->lwork,&info_in);
    mat->row_shift=mat->na;
  }
  else
  {
    dgeqrf_(&mat->ma,&mat->na,mat->aa,&mat->ma,mat->tau,mat->work,&mat->lwork,&info_in);
  }

} 


void SparseOperation(CG_PAR *const cg, MAT_DAT *const mat)
{
  //calculate total number of none-zero elements in this block
  int n_nz=0;
  int k,l,nt;
  double tx;

  for(k=0;k<mat->mv;k++)
  {
    n_nz+=mat->mat_head[k].n;
  }
  n_nz*=3;

  if(mat->mp>0)
  {
    for(k=0;k<mat->mp*mat->nn;k++)
    {
      if(mat->aa[k]>VERYSMALL || mat->aa[k]<-VERYSMALL) n_nz++;
    }
  }

  mat->rw=(double *)MemAlloc(n_nz*sizeof(double));
  mat->jw=(int *)MemAlloc(n_nz*sizeof(int));


  //construct CSR format and free sparse matrix linked lists

  struct MAT_SP *p,*p1;
  int num,tn,tn1;
  for(k=0;k<mat->mv;k++)
  {
    p=mat->mat_head[k].h;
    
    num=0;
    while(p!=NULL)
    {
      tn=mat->iw[3*k];
      tn1=mat->mat_head[k].n;
      mat->rw[tn+num]=p->valx[0]; 
      mat->rw[tn+tn1+num]=p->valx[1];
      mat->rw[tn+tn1*2+num]=p->valx[2];
      mat->jw[tn+num]=p->col;
      mat->jw[tn+tn1+num]=p->col;
      mat->jw[tn+tn1*2+num]=p->col;
      p1=p;
      p=p1->next;
      free(p1);
      num++;
    }
    nt=3*k;
    mat->iw[nt+1]=mat->iw[nt]+num;
    mat->iw[nt+2]=mat->iw[nt+1]+num;
    mat->iw[nt+3]=mat->iw[nt+2]+num;
  }

  if(mat->mp>0)
  {
    num=mat->iw[mat->mv*3];
    for(k=0;k<mat->mp;k++)
    {
      for(l=0;l<mat->nn;l++)
      {
        tx=mat->aa[l*mat->mp+k];
        if(tx>VERYSMALL || tx<-VERYSMALL)
        {
          mat->rw[num]=tx;
          mat->jw[num]=l;
          num++;
        }
        mat->iw[mat->mv*3+k+1]=num;
      }
    }
  }


  //a simple precondition
  for(k=0;k<mat->nn;k++)
  {
    mat->h[k]=0.0;
  }

  for(k=0;k<mat->mm;k++)
  {
    for(l=mat->iw[k];l<mat->iw[k+1];l++)
    {
      mat->h[mat->jw[l]]+=mat->rw[l]*mat->rw[l];
    }
  }

  for(k=0;k<mat->nn;k++)
  {
    if(mat->h[k]>VERYSMALL) mat->h[k]=1.0/sqrt(mat->h[k]); 
    else mat->h[k]=1.0;
  }

  for(k=0;k<mat->mm;k++)
  {
    for(l=mat->iw[k];l<mat->iw[k+1];l++)
    {
      mat->rw[l]*=mat->h[mat->jw[l]];
    }
  }
  


  //solve the equation
  pda_lsqr_(&mat->mm,&mat->nn,&cg->damp,&cg->mm1,&n_nz,&n_nz,mat->iw,mat->jw,mat->rw,mat->bb,mat->v,mat->w,mat->xx_b,
                 mat->se,&cg->atol,&cg->btol,&cg->conlim,&cg->itnlim,&cg->istop,&cg->itn,&cg->anorm,
                 &cg->acond, &cg->rnorm, &cg->arnorm, &cg->xnorm);

     
  free(mat->jw); free(mat->rw);

  for(k=0;k<mat->nn;k++)
  {
    mat->xx_b[k]*=mat->h[k];
    if((mat->xx_b[k]>VERYSMALL || mat->xx_b[k]<-VERYSMALL) && (mat->xx_b[k]<FORCE_CAP && mat->xx_b[k]>-FORCE_CAP)) 
    {
      mat->n_xx[k]++;
      mat->xx[k]+=mat->xx_b[k];
    }
  }

}


void PostOperationSparse (CG_PAR *const cg, MAT_DAT *const mat)
{
   if(cg->con_p>=2)
   {
    FILE *mat_out;
    mat_out=OpenFile("result.out","wb");
    fwrite(&mat->xx[0],sizeof(double),mat->nn,mat_out);
    fwrite(&mat->n_xx[0],sizeof(int),mat->nn,mat_out);
    fclose(mat_out);
   }

  if(cg->con_p==3) exit(0);

  int i;
  for(i=0;i<mat->nn;i++)
  {
    if(mat->n_xx[i]>0)
    {
      mat->xx[i]=mat->xx[i]/mat->n_xx[i];
    }
    else mat->xx[i]=0;
  }
  free(mat->n_xx);

}

void PostOperationDense (CG_PAR *const cg, MAT_DAT *const mat)
{

  double *dd_bak;
  if(cg->out_dd==1) 
  {
    dd_bak=(double *)MemAlloc(mat->nn*sizeof(double)); 
    for(int i=0;i<mat->nn;i++) dd_bak[i]=mat->dd[i];
  }

  if(cg->iter==1)
  {
    FILE *mat_in;
    mat_in=OpenFile("result.in","rb");
    double ttx;
    double *dd1;
    dd1=(double *)MemAlloc(mat->nn*sizeof(double));

    for(int j=0;j<mat->nn;j++)
    {
      for(int k=0;k<=j;k++)
      {
        fread(&ttx,sizeof(double),1,mat_in);
        mat->gg[j*mat->nn+k]=ttx;
      }
    }
    for(int j=0;j<mat->nn;j++)
    {
      fread(&ttx,sizeof(double),1,mat_in);
      dd1[j]=ttx; 
    }
    fclose(mat_in);
    
    for(int i=0;i<mat->nn;i++) 
    mat->dd[i]=dd1[i]-mat->dd[i];
    //mat->dd[i]=2.0*dd1[i]-mat->dd[i];
    free(dd1);

    goto s;
  }


  //save the results in binary form for parallel runs
  if(cg->con_p>=2)
  {
    FILE *mat_out;
    mat_out=OpenFile("result.out","wb");
    int i;
    for(i=0;i<mat->nn;i++)
    {
      fwrite(&mat->gg[i*mat->nn],sizeof(double),i+1,mat_out);
    }
    fwrite(&mat->dd[0],sizeof(double),mat->nn,mat_out);
    fclose(mat_out);

  }


  if(cg->con_p==3) exit(0);


  int i,j;
s:  
  for(i=0;i<mat->nn;i++)
  {
    for(j=0;j<i;j++)
    {
      mat->gg[j*mat->nn+i]=mat->gg[i*mat->nn+j];
    }
  }
/*
for(int i=0;i<mat->nn;i++)
{
  for(int j=0;j<mat->nn;j++)
  {
    if(i!=j) mat->gg[j*mat->nn+i]=0.0;
    //else printf("%.15le\n",mat->gg[j*mat->nn+i]);
  }
}
*/



/*
  double *gg_bak;
  gg_bak=(double *)MemAlloc((mat->nn*mat->nn)*sizeof(double));
  for(int i=0;i<mat->nn*mat->nn;i++) gg_bak[i]=mat->gg[i];
for(int i=0;i<79;i++) {dd_bak[i]=mat->dd[i]; printf("%le\n",dd_bak[i]);}
exit(0);

  char tmp[100];
  FILE *tf;
  for(int j=0;j<mat->nn;j++)
  {
    sprintf(tmp,"%d.dat",j);
    tf=fopen(tmp,"w");
    for(int i=0;i<mat->nn;i++) fprintf(tf,"%le\n",mat->gg[i*mat->nn+j]);
    fclose(tf); 
  }
exit(0);
*/


/*  
 char tchar='N';
 double alpha;
 double done=1.0;
 int ione=1;
 alpha=ddot(&mat->nn,mat->dd,&ione,mat->dd,&ione);
 
 double *vt;
 vt=(double *)MemAlloc((unsigned)mat->nn*sizeof(double));

 dgemv(&tchar,&mat->nn,&mat->nn,&done,mat->gg,&mat->nn,mat->dd,&ione,&done,vt,&ione);
 alpha=alpha/ddot(&mat->nn,vt,&ione,mat->dd,&ione);
 
 mat->xx=(double *)MemAlloc(mat->nn*sizeof(double));
 for(int i=0;i<mat->nn;i++) mat->xx[i]=mat->dd[i]*alpha;
*/


  //a simple preconditioning
  double *h;
  h=(double *)malloc((unsigned)mat->nn*sizeof(double));

  for(i=0;i<mat->nn;i++)
  {
    h[i]=0.0;
  }

  for(i=0;i<mat->nn;i++)
  {
    for(j=0;j<mat->nn;j++)
    {
      h[j]=h[j]+mat->gg[j*mat->nn+i]*mat->gg[j*mat->nn+i];
    }
  }

  for(i=0;i<mat->nn;i++)
  {
    if(h[i]<VERYSMALL) h[i]=1.0;
    else h[i]=1.0/sqrt(h[i]);
  }

  for(i=0;i<mat->nn;i++)
  {
    for(j=0;j<mat->nn;j++)
    {
      mat->gg[j*mat->nn+i]*=h[j];
    }
  }

/*
  char tchar='N';
 double alpha;
 double done=1.0;
 int ione=1;
 alpha=ddot(&mat->nn,mat->dd,&ione,mat->dd,&ione);

 double *vt;
 vt=(double *)MemAlloc((unsigned)mat->nn*sizeof(double));

 dgemv(&tchar,&mat->nn,&mat->nn,&done,mat->gg,&mat->nn,mat->dd,&ione,&done,vt,&ione);
 alpha=alpha/ddot(&mat->nn,vt,&ione,mat->dd,&ione);

 mat->xx=(double *)MemAlloc(mat->nn*sizeof(double));
 for(int i=0;i<mat->nn;i++) mat->xx[i]=mat->dd[i]*alpha;
*/



  //Tikhonov renormalization
  double lambda2;
  if(cg->input_lambda==1)
  {
    lambda2=cg->lambda*cg->lambda;
    for(i=0;i<mat->nn;i++)
    {
      mat->gg[i*mat->nn+i]+=lambda2;
    }
  } 



  //solving the normal equation by SVD
  int space_factor=2;
  double *sing_in;
  sing_in=(double *)malloc((unsigned int)mat->nn*sizeof(double));
  int onei=1;

  //  dgelsd(nn,nn,1,gg,nn,dd,nn,sing_in,-1.0,&irank_in,&info_in);
  //double rcond=-1.0;
  int lwork=-1;
  double *work;
  work=(double *)malloc((unsigned)1*sizeof(double));
  int irank_in,info_in,liwork,*iwork,nlvl,smlsiz=25;
  double tx=mat->nn/(smlsiz+1.0);
  nlvl=(int)(log10(tx)/log10(2))+1;
  liwork=3*mat->nn*nlvl+11*mat->nn;
  liwork=space_factor*liwork;
  iwork=(int *)malloc((unsigned)liwork*sizeof(int));

  dgelsd_(&mat->nn,&mat->nn,&onei,mat->gg,&mat->nn,mat->dd,&mat->nn,sing_in,&cg->rcond,&irank_in,work,&lwork,iwork,&info_in);
  lwork=work[0];
  free(work);
  work=(double *)malloc((unsigned)lwork*sizeof(double));
  dgelsd_(&mat->nn,&mat->nn,&onei,mat->gg,&mat->nn,mat->dd,&mat->nn,sing_in,&cg->rcond,&irank_in,work,&lwork,iwork,&info_in);

  //print singular values
  FILE *sol;
  sol=OpenFile("sol_info.out","a");
  fprintf(sol,"Singular vector:\n");
  for(i=0;i<mat->nn;i++)
  {
    fprintf(sol,"%le\n",sing_in[i]);
  }

  //final results 
  mat->xx=(double *)MemAlloc(mat->nn*sizeof(double));

  for(i=0;i<mat->nn;i++)
  { 
    mat->xx[i]=mat->dd[i]*h[i];
  }

  if(cg->iter==1)
  {
    FILE *x_in;
    x_in=OpenFile("x.in","r");
    double *x0;
    x0=(double *)MemAlloc(mat->nn*sizeof(double));
    for(int i=0;i<mat->nn;i++) fscanf(x_in,"%le",x0+i);
    fclose(x_in);
//test
    for(i=0;i<mat->nn;i++) mat->xx[i]=mat->xx[i]*cg->iter_ld+x0[i];
    free(x0);
  }

/*
  //print 2-norm of solution and residual
  double xnorm=0.0;
  for(i=0;i<mat->nn;i++)
  {
    xnorm+=mat->xx[i]*mat->xx[i];
  }
  fprintf(sol,"Solution 2-norm:\n%le\n",xnorm);
*/
  if(cg->out_dd==1)
  {
    for(int i=0;i<mat->nn;i++) mat->dd[i]=dd_bak[i];
    free(dd_bak);
  }
   
  //fclose(sol);
  //free(sing_in); free(work); free(iwork);
  //free(h); 
  free(mat->gg);
}

void PostOperationAccumulation (CG_PAR *const cg, MAT_DAT *const mat)
{
  int i,j;
  for(i=0;i<mat->na;i++)
  {
    mat->dd[i]=mat->aa[mat->nn*mat->ma+i];
  }
  SetZeroAccumulation(mat);

  //save the results in binary form
  if(cg->con_p>=2)
  { 
    FILE *mat_out;
    mat_out=OpenFile("result.out","wb");
    int i;
    for(i=0;i<mat->nn;i++)
    {
      fwrite(&mat->aa[i*mat->ma],sizeof(double),i+1,mat_out);
    }
    fwrite(&mat->dd[0],sizeof(double),mat->na,mat_out);
    fclose(mat_out);
  }
  
  if(cg->con_p==3) exit(0);


  double resid;
  resid=mat->aa[(mat->na-1)*mat->ma+mat->na-1];

  //a simple preconditioning
  double *h;
  h=(double *)malloc((unsigned)mat->nn*sizeof(double));

  for(i=0;i<mat->nn;i++)
  {
    h[i]=0.0;
  }

  for(i=0;i<mat->nn;i++)
  {
    for(j=0;j<mat->nn;j++)
    {
      h[j]=h[j]+mat->aa[j*mat->ma+i]*mat->aa[j*mat->ma+i];
    }
  }

  for(i=0;i<mat->nn;i++)
  {
    if(h[i]<VERYSMALL) h[i]=1.0;
    else h[i]=1.0/sqrt(h[i]);
  }

  for(i=0;i<mat->nn;i++)
  {
    for(j=0;j<mat->nn;j++)
    {
      mat->aa[j*mat->ma+i]*=h[j];
    }
  }

  //solving the equation by SVD
  double *sing_in;
  sing_in=(double *)MemAlloc(mat->nn*sizeof(double));
  int onei=1;

  int lwork=-1;
  double *work;
  work=(double *)MemAlloc(1*sizeof(double));
  int irank_in,info_in;
  dgelss_(&mat->nn,&mat->nn,&onei,mat->aa,&mat->ma,mat->dd,&mat->nn,sing_in,&cg->rcond,&irank_in,work,&lwork,&info_in);
  lwork=work[0];
  free(work);
  work=(double *)MemAlloc(lwork*sizeof(double));
  dgelss_(&mat->nn,&mat->nn,&onei,mat->aa,&mat->ma,mat->dd,&mat->nn,sing_in,&cg->rcond,&irank_in,work,&lwork,&info_in);

  //print singular values
  FILE *sol;
  sol=OpenFile("sol_info.out","a");
  fprintf(sol,"Singular vector:\n");
  for(i=0;i<mat->nn;i++)
  {
    fprintf(sol,"%le\n",sing_in[i]);
  }

  //final results 
  mat->xx=(double *)MemAlloc(mat->nn*sizeof(double));
  for(i=0;i<mat->nn;i++)
  { 
    mat->xx[i]=mat->dd[i]*h[i];
  }

  //print 2-norm of solution and residual
  double xnorm=0.0;
  for(i=0;i<mat->nn;i++)
  {
    xnorm+=mat->xx[i]*mat->xx[i];
  }
  fprintf(sol,"Solution 2-norm:\n%le\n",xnorm);

  fprintf(sol,"Residual 2-norm:\n%le\n",resid);

  fclose(sol);
  free(sing_in); free(work); 
  free(h); free(mat->dd); free(mat->aa);
}

