void FreeAVSpace(CG_PAR *const cg, MAT_DAT *const mat, FR_DAT *const fr)
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
  if(cg->be_sparse==2) free(mat->tau);
}

void ReadBinDense(MAT_DAT *const mat)
{
  //read the matrix
  int i,j,k;
  double tx;
  int n_av; //number of blocks
  char **mat_name; //name of matrix files
  FILE *mat_av; //input file
  mat_av=OpenFile("res_av.in","r");
  fscanf(mat_av,"%d",&n_av);
  mat_name=(char **)malloc((unsigned)n_av*sizeof(char *));
  for(i=0;i<n_av;i++)
  {
    mat_name[i]=(char *)malloc((unsigned)100*sizeof(char)); //max matrix file name length is 100
  }
  for(i=0;i<n_av;i++)
  {
    fscanf(mat_av,"%s",mat_name[i]);
  }
  fclose(mat_av);

  FILE *mat_out;

  for(i=0;i<n_av;i++)
  {
    mat_out=OpenFile(mat_name[i],"rb");
    for(j=0;j<mat->nn;j++)
    {
      for(k=0;k<=j;k++)
      {
        fread(&tx,sizeof(double),1,mat_out);
        mat->gg[j*mat->nn+k]+=tx;
      }
    }

    for(j=0;j<mat->nn;j++)
    {
      fread(&tx,sizeof(double),1,mat_out);
      mat->dd[j]+=tx;
    }
    fclose(mat_out);
  }

  for(i=0;i<mat->nn;i++)
  {
    for(j=0;j<i;j++)
    {
      mat->gg[j*mat->nn+i]=mat->gg[i*mat->nn+j];
    }
  }

}

void ReadBinAccumulation(MAT_DAT *const mat)
{
  free(mat->aa);
  mat->aa=(double *)MemAlloc(mat->na*mat->na*sizeof(double));

  //read the matrix
  int i,j,k;
  double tx;
  int n_av; //number of blocks
  char **mat_name; //name of matrix files
  FILE *mat_av; //input file
  mat_av=OpenFile("res_av.in","r");
  fscanf(mat_av,"%d",&n_av);
  if(n_av>1) {printf("Can not read more than one block for the sequential accumulation algorithm!\n"); exit(1);}
  mat_name=(char **)malloc((unsigned)n_av*sizeof(char *));

  mat_name[0]=(char *)malloc((unsigned)100*sizeof(char)); //max matrix file name length is 100

  fscanf(mat_av,"%s",mat_name[0]);
  fclose(mat_av);

  FILE *mat_out;

  mat_out=OpenFile(mat_name[0],"rb");
  for(j=0;j<mat->nn;j++)
  {
    for(k=0;k<=j;k++)
    {
      fread(&tx,sizeof(double),1,mat_out);
      mat->aa[j*mat->na+k]=tx;
    }
  }

  for(j=0;j<mat->nn;j++) mat->aa[j*mat->na+mat->nn]=0.0;

  for(j=0;j<mat->na;j++)
  {
    fread(&tx,sizeof(double),1,mat_out);
    mat->dd[j]=tx;
  }
  fclose(mat_out);

}


void ReadBinSparse(MAT_DAT *const mat)
{
  //read the results in binary form for parallel runs
  int i,j;
  int n_av;
  FILE *input,*av;
  char av_in[100];
  int *nb_xx;
  nb_xx=(int *)malloc((unsigned)mat->nn*sizeof(int));
  input=OpenFile("res_av.in","r");
  fscanf(input,"%d",&n_av);
  for(i=0;i<n_av;i++)
  {
    fscanf(input,"%s",av_in);
    av=OpenFile(av_in,"rb");
    fread(mat->xx_b,sizeof(double),mat->nn,av);
    fread(nb_xx,sizeof(int),mat->nn,av);
    fclose(av);
    for(j=0;j<mat->nn;j++)
    {
      mat->xx[j]+=mat->xx_b[j];
      mat->n_xx[j]+=nb_xx[j];
    }
  }
  fclose(input);
  free(nb_xx);
}

void PostOperationAccumulationReg (CG_PAR *const cg, MAT_DAT *const mat)
{
  int i,j,k;
  double resid;
  resid=fabs(mat->dd[mat->nn]);


  //a simple preconditioning
  double *h;
  h=(double *)MemAlloc(mat->nn*sizeof(double));

/*
if(cg->input_lambda==0){
  for(i=0;i<mat->nn;i++)
  {
    h[i]=0.0;
  }

  for(i=0;i<mat->nn;i++)
  {
    for(j=0;j<mat->nn;j++)
    {
      h[j]=h[j]+mat->aa[j*mat->na+i]*mat->aa[j*mat->na+i];
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
      mat->aa[j*mat->na+i]*=h[j];
    }
  }
}
*/
  FILE *lambda_table;
  int n_l,l_mod;
  double *lam,*lam2,tx,tx1;
  double xnorm,*xnorm_l,*resid_l,*xnorm_l1,*resid_l1;
  FILE *lout;
  
  FILE *sol;
  double *sing_in;
  sing_in=(double *)MemAlloc(mat->nn*sizeof(double));
  int onei=1;
  int lwork=-1;
  double *work;
  work=(double *)MemAlloc(1*sizeof(double));
  int irank_in,info_in;
  mat->xx=(double *)MemAlloc(mat->nn*sizeof(double));

  char ochar='O';
  char nchar='N';
  char schar='S';
  
  double *uu=(double *)MemAlloc(mat->na*mat->nn*sizeof(double));
  double *vt=(double *)MemAlloc(mat->nn*mat->nn*sizeof(double));
  dgesvd_(&schar,&schar,&mat->na,&mat->nn,mat->aa,&mat->na,sing_in,uu,&mat->na,vt,&mat->nn,work,&lwork,&info_in);
  lwork=work[0];
  free(work);
  work=(double *)MemAlloc(lwork*sizeof(double));
  dgesvd_(&schar,&schar,&mat->na,&mat->nn,mat->aa,&mat->na,sing_in,uu,&mat->na,vt,&mat->nn,work,&lwork,&info_in);

  double s,s1,s2,l2;
  double *tmp=(double *)MemAlloc(mat->nn*sizeof(double));
  double eps;
  if(cg->rcond<0) eps=DBL_EPSILON;
  else eps=cg->rcond;

  if(cg->input_lambda==0 || cg->input_lambda==1)
  {
    sol=OpenFile("sol_info.out","w");
    //print singular values
    fprintf(sol,"Singular vector:\n");
    for(i=0;i<mat->nn;i++)
    {
      fprintf(sol,"%le\n",sing_in[i]);
    }
        
    if(cg->input_lambda==0)
    {

      for(j=0;j<mat->nn;j++)
      {
        s=0.0;
        if(sing_in[j]>eps)
        {
          for(i=0;i<mat->na;i++)
          {
            s+=uu[j*mat->na+i]*mat->dd[i];
          }
          s/=sing_in[j];

        }
        tmp[j]=s;
      }

      for(j=0;j<mat->nn;j++)
      {
        s=0.0;
        for(i=0;i<mat->nn;i++)
        {
          s+=vt[j*mat->nn+i]*tmp[i];
        }
        mat->xx[j]=s;
      }
   }

   else
   {
     l2=cg->lambda*cg->lambda;
     for(j=0;j<mat->nn;j++)
      {
        s=0.0;
        for(i=0;i<mat->na;i++)
        {
          s+=uu[j*mat->na+i]*mat->dd[i];
        }
        s=s*sing_in[j]/(sing_in[j]*sing_in[j]+l2);
        tmp[j]=s;
      }

      for(j=0;j<mat->nn;j++)
      {
        s=0.0;
        for(i=0;i<mat->nn;i++)
        {
          s+=vt[j*mat->nn+i]*tmp[i];
        }
        mat->xx[j]=s;
      }
   }  

/*
    //print 2-norm of solution and residual
    xnorm=0.0;
    for(i=0;i<mat->nn;i++)
    {
      xnorm+=mat->xx[i]*mat->xx[i];
    }
    fprintf(sol,"Solution 2-norm after column scaling:\n%le\n",xnorm);


    //final results 
    for(i=0;i<mat->nn;i++)
    { 
      mat->xx[i]=mat->xx[i]*h[i];
    }

    //print 2-norm of solution and residual
    xnorm=0.0;
    for(i=0;i<mat->nn;i++)
    {
      xnorm+=mat->xx[i]*mat->xx[i];
    }
    fprintf(sol,"Solution 2-norm:\n%le\n",xnorm);

    fprintf(sol,"Residual 2-norm:\n%le\n",resid);
    fclose(sol);
*/
  } //end of if(input_lambda==0)

  else
  {
    lambda_table=fopen("lambda.in","r");
    fscanf(lambda_table,"%d%d",&n_l,&l_mod);
    lout=fopen("l.out","w");

    lam=(double *)MemAlloc(n_l*sizeof(double));
    lam2=(double *)MemAlloc(n_l*sizeof(double));
    xnorm_l=(double *)MemAlloc(n_l*sizeof(double));
    resid_l=(double *)MemAlloc(n_l*sizeof(double));
    xnorm_l1=(double *)MemAlloc(n_l*sizeof(double));

    for(k=0;k<n_l;k++)
    {
      fscanf(lambda_table,"%lf",lam+k);
      if(l_mod==1) lam[k]=pow(10.0,lam[k]);
      lam2[k]=lam[k]*lam[k];
      xnorm_l[k]=0.0;
      resid_l[k]=0.0;
      xnorm_l1[k]=0.0;
    }
      

    for(j=0;j<mat->nn;j++)
    {
      s=0.0;
      for(i=0;i<mat->na;i++)
      {
        s+=uu[j*mat->na+i]*mat->dd[i];
      }
      tx1=sing_in[j]*sing_in[j];
      for(k=0;k<n_l;k++)
      {
        tx=tx1+lam2[k];
        s2=s/tx;
        s1=s2*sing_in[j];
        xnorm_l[k]+=s1*s1;
        s1=s2*lam2[k];
        resid_l[k]+=s1*s1;
        xnorm_l1[k]+=s*s*lam2[k]*tx1/(tx*tx*tx);
      }
    }

    for(i=0;i<n_l;i++)
    {
      xnorm_l1[i]*=(-4.0/lam[i]);
      tx1=xnorm_l[i]*resid_l[i];
      tx=2.0*tx1/xnorm_l1[i]*(lam2[i]*xnorm_l1[i]*resid_l[i]+2.0*lam[i]*tx1+lam2[i]*lam2[i]*xnorm_l[i]*xnorm_l1[i])/pow(lam2[i]*xnorm_l[i]*xnorm_l[i]+resid_l[i]*resid_l[i],1.5);
      fprintf(lout,"%19.14le %19.14le %19.14le %19.14le\n",resid_l[i]+resid,xnorm_l[i],lam[i],-tx);
      
    }
   
    fclose(lout);
    free(lam); free(lam2); free(xnorm_l); free(xnorm_l1); free(resid_l);
  }
  

  free(sing_in); free(work); 
  free(h); free(mat->dd); free(mat->aa);
  if(cg->input_lambda==2) exit(0);
}
