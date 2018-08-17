void WriteOutput (CG_PAR *const cg, MAT_DAT *const mat, FILE *const nb_out, FILE *const bon_out, void (*const output)(enum INTER_TYPE, CG_PAR *const, MAT_DAT *const, TYPE_INFO *const, FILE *const))
{
  TYPE_INFO info;
  enum INTER_TYPE type;
  int i,j,k,l;

  void OutputBspline_3(enum INTER_TYPE, CG_PAR *const, MAT_DAT *const, TYPE_INFO *const, FILE *const);

  //output unknowns solved and force curves
  //nonebonded
  type=non;
  InitialType(type,&info,cg);
  for(i=0;i<cg->tolcgtype;i++)
  {
    for(j=i;j<cg->tolcgtype;j++)
    {
      info.m=TwoBodyNum(i+1,j+1,cg->tolcgtype);
      info.mi=info.m;
      info.numi=cg->ljnum[info.mi];
      info.i=i; info.j=j;
      if(info.numi>0)
      {
        (*output)(type,cg,mat,&info,nb_out);
      }
    }
  }

  //bonds
  type=bon;
  InitialType(type,&info,cg);
  for(i=0;i<cg->tolcgtype;i++)
  {
    for(j=i;j<cg->tolcgtype;j++)
    {
      info.m=TwoBodyNum(i+1,j+1,cg->tolcgtype);
      info.mi=SearchIntTable(cg->bon_m,info.m,cg->tolb);
      if(info.mi==-1) continue;
      info.numi=cg->bonnum[info.mi];
      info.i=i; info.j=j;
      if(info.numi>0)
      {
        (*output)(type,cg,mat,&info,bon_out);
      }
    }
  }

  //angles
  type=ang;
  InitialType(type,&info,cg);
  for(i=0;i<cg->tolcgtype;i++)
  {
    for(j=0;j<cg->tolcgtype;j++)
    {
      for(k=j;k<cg->tolcgtype;k++)
      {
        info.m=ThreeBodyNum(i+1,j+1,k+1,cg->tolcgtype);
        info.mi=SearchIntTable(cg->ang_m,info.m,cg->tola);
        if(info.mi==-1) continue;
        info.numi=cg->angnum[info.mi];
        info.i=i; info.j=j; info.k=k;
        if(info.numi>0)
        {
          (*output)(type,cg,mat,&info,bon_out);
        }
      }
    }
  }

  //hihedrals
  type=dih;
  InitialType(type,&info,cg);
  for(i=0;i<cg->tolcgtype;i++)
  {
    j=i;
    for(k=0;k<cg->tolcgtype;k++)
    {
      for(l=k;l<cg->tolcgtype;l++)
      {
        info.m=FourBodyNum(i+1,j+1,k+1,l+1,cg->tolcgtype);
        info.mi=SearchIntTable(cg->dih_m,info.m,cg->told);
        if(info.mi==-1) continue;
        info.numi=cg->dihnum[info.mi];
        info.i=i; info.j=j; info.k=k; info.l=l;
        if(info.numi>0)
        {
          (*output)(type,cg,mat,&info,bon_out);
        }
      }
    }
 
    for(j=i+1;j<cg->tolcgtype;j++)
    {
      for(k=0;k<cg->tolcgtype;k++)
      {
        for(l=0;l<cg->tolcgtype;l++)
        {
          info.m=FourBodyNum(i+1,j+1,k+1,l+1,cg->tolcgtype);
          info.mi=SearchIntTable(cg->dih_m,info.m,cg->told);
          if(info.mi==-1) continue;
          info.numi=cg->dihnum[info.mi];
          info.i=i; info.j=j; info.k=k; info.l=l;
          if(info.numi>0)
          {
            (*output)(type,cg,mat,&info,bon_out);
          }

        }
      }
    }
  }

  FILE *tb_out;
  if(cg->use_3b==3) tb_out=fopen("3b.dat","w");
  if(cg->use_3b>0)
  {
    type=t_b;
    InitialType(type,&info,cg);
    for(i=0;i<cg->tolcgtype;i++)
    {
      for(j=0;j<cg->tolcgtype;j++)
      {
        for(k=j;k<cg->tolcgtype;k++)
        {
          info.m=ThreeBodyNum(i+1,j+1,k+1,cg->tolcgtype);
          info.mi=SearchIntTable(cg->tb_m,info.m,cg->tol3);
          if(info.mi==-1) continue;
          info.numi=cg->tbnum[info.mi];
          info.i=i; info.j=j; info.k=k;
          if(info.numi>0)
          {
            if(cg->use_3b==3)
            {
              info.numi--;
          
              fprintf(tb_out,"%.15le\n",mat->xx[info.grid_base+info.grid[info.numi]]);
            }
            else OutputBspline_3(type,cg,mat,&info,bon_out);
          }
        }
      }
    }

  }

  if(cg->use_3b==3) fclose(tb_out);


  if(cg->be_sparse==0) free(mat->dd);
  free(mat->xx);
  free(cg->ljnum);
  free(cg->angnum);
  free(cg->bonnum);
  free(cg->dihnum);
  free(cg->bon_m);
  free(cg->ang_m);
  free(cg->dih_m);
  free(cg->rmin); free(cg->rmax);
  free(cg->rmin_b); free(cg->rmax_b);
  free(cg->rmin_a); free(cg->rmax_a);
  free(cg->rmin_d); free(cg->rmax_d);
  free(cg->grid_n);
  free(cg->grid_bn);
  free(cg->grid_an);
  free(cg->grid_dn);
  for(i=0;i<cg->tolcgtype;i++) free(cg->name[i]);
  free(cg->name);
  if(cg->ba_type==0)
  {
    for(i=0;i<cg->tollj1;i++)  gsl_bspline_free(cg->w_bs[i]); 
    for(i=0;i<cg->tolb1;i++)  gsl_bspline_free(cg->wb_bs[i]); 
    for(i=0;i<cg->tola1;i++)  gsl_bspline_free(cg->wa_bs[i]);  
    for(i=0;i<cg->told1;i++)  gsl_bspline_free(cg->wd_bs[i]); 
    free(cg->w_bs); gsl_vector_free(cg->b_bs);
    free(cg->wb_bs); gsl_vector_free(cg->bb_bs);
    free(cg->wa_bs); gsl_vector_free(cg->ba_bs);
    free(cg->wd_bs); gsl_vector_free(cg->bd_bs);
  }

  if(cg->use_3b>0)
  {
    free(cg->cutoff_3i);
    free(cg->cos_theta0i);
    free(cg->tbnum);
    free(cg->tb_m);
    free(cg->rmin_3);
    free(cg->rmax_3);
    free(cg->grid_3n);
    if(cg->ba_type==0 && cg->use_3b!=3)
    {
      for(i=0;i<cg->tol3;i++)  gsl_bspline_free(cg->w3_bs[i]); 
      free(cg->w3_bs); gsl_vector_free(cg->b3_bs);
    }
    if(cg->use_3b==2 && cg->use_3b!=3) 
    {
      gsl_matrix_free(cg->b3_bs1);
      gsl_bspline_deriv_free(cg->w3_bs1);
    }
  }


}

void WriteOutputFM(CG_PAR *const cg, MAT_DAT *const mat)
{
  FILE *br_out;
  if(cg->ba_type==0) br_out=OpenFile("b-spline.out","w");
  else br_out=NULL;

  void OutputBspline(enum INTER_TYPE, CG_PAR *const, MAT_DAT *const, TYPE_INFO *const, FILE *const);
  void OutputLinear(enum INTER_TYPE, CG_PAR *const, MAT_DAT *const, TYPE_INFO *const, FILE *const);
  void (*output)(enum INTER_TYPE, CG_PAR *const, MAT_DAT *const, TYPE_INFO *const, FILE *const);

  FILE *xout;
  if(cg->x_out==1) {
  xout=fopen("x.out","wb");
  fwrite(mat->xx,sizeof(double),mat->nn,xout);
  fclose(xout); }

  if(cg->ba_type==0) output=OutputBspline;
  else if(cg->ba_type==1) output=OutputLinear;

  WriteOutput(cg,mat,br_out,br_out,output);
  if(cg->ba_type==0) fclose(br_out);
}

char InitialOutput(enum INTER_TYPE type, CG_PAR *const cg, TYPE_INFO *const info, char *name_tmp, double *spaceo)
{
  if(type==non)
  {
    sprintf(name_tmp,"%s_%s.dat",cg->name[info->i],cg->name[info->j]);
    //if(info->rmax[info->mi]>cg->cutoff+VERYSMALL_F) info->rmax[info->mi]=cg->cutoff+VERYSMALL_F;
    *spaceo=cg->spaceo; return 'n';
  }
  else if(type==bon)
  {
    sprintf(name_tmp,"%s_%s_bon.dat",cg->name[info->i],cg->name[info->j]);
    *spaceo=cg->spaceo_b; return 'b';
  }
  else if(type==ang)
  {
    sprintf(name_tmp,"%s_%s_%s_ang.dat",cg->name[info->i],cg->name[info->j],cg->name[info->k]);
    *spaceo=cg->spaceo_a; return 'a';
  }
  else if(type==dih)   
  {     
    sprintf(name_tmp,"%s_%s_%s_%s_dih.dat",cg->name[info->i],cg->name[info->j],cg->name[info->k],cg->name[info->l]);
    *spaceo=cg->spaceo_d; return 'd';
  }
  else
  {
    sprintf(name_tmp,"%s_%s_%s.dat",cg->name[info->i],cg->name[info->j],cg->name[info->k]);
    *spaceo=cg->spaceo_3; return '3';
  }


}

void OutputBspline(enum INTER_TYPE type, CG_PAR *const cg, MAT_DAT *const mat, TYPE_INFO *const info, FILE *const br_out)
{
  int grid,n_br,k,tn;
  char name_tmp[100];
  double rr,tx;
  FILE *ff;
  char i_type;
  double spaceo;
  size_t istart,iend;
  i_type=InitialOutput(type,cg,info,name_tmp,&spaceo);

  info->numi--;
  grid=info->grid[info->numi+1]-info->grid[info->numi];
  n_br=grid-info->n_k+2;
  fprintf(br_out,"%c: %s %s %d %d %.15le %.15le\n",i_type,cg->name[info->i],cg->name[info->j],info->n_k,n_br,info->rmin[info->mi],info->rmax[info->mi]);
  for(k=0;k<grid;k++)
  {
    fprintf(br_out,"%.15le ",mat->xx[info->grid[info->numi]+k]);
  }
  fprintf(br_out,"\n");
  ff=OpenFile(name_tmp,"w");

  double max;
  max=info->rmax[info->mi];
  for(rr=((int)(info->rmin[info->mi]/spaceo)+1)*spaceo;rr<max;rr+=spaceo)
  {
    gsl_bspline_eval_nonzero(rr,info->bs,&istart,&iend,info->w[info->numi]);
    tx=0.0;
    //printf("%d %d %d\n",istart,iend,info->n_k);
    for(tn=istart;tn<=iend;tn++)
    {
      tx+=gsl_vector_get(info->bs,tn-istart)*mat->xx[info->grid_base+info->grid[info->numi]+tn];
    }
    fprintf(ff,"%lf %.15le\n",rr,tx);
  }
  fclose(ff);
}


void OutputBspline_3(enum INTER_TYPE type, CG_PAR *const cg, MAT_DAT *const mat, TYPE_INFO *const info, FILE *const br_out)
{
  int grid,n_br,k,tn;
  char name_tmp[100];
  double rr,tx,tx1;
  FILE *ff;
  char i_type;
  double spaceo;
  size_t istart,iend;
  i_type=InitialOutput(type,cg,info,name_tmp,&spaceo);

  info->numi--;
  grid=info->grid[info->numi+1]-info->grid[info->numi];
  n_br=grid-info->n_k+2;
  fprintf(br_out,"%c: %s %s %d %d %.15le %.15le\n",i_type,cg->name[info->i],cg->name[info->j],info->n_k,n_br,info->rmin[info->mi],info->rmax[info->mi]);
  for(k=0;k<grid;k++)
  {
    fprintf(br_out,"%.15le ",mat->xx[info->grid[info->numi]+k]);
  }
  fprintf(br_out,"\n");
  ff=OpenFile(name_tmp,"w");

  double max;
  max=info->rmax[info->mi];
  for(rr=((int)(info->rmin[info->mi]/spaceo)+1)*spaceo;rr<max;rr+=spaceo)
  {
    gsl_bspline_eval_nonzero(rr,info->bs,&istart,&iend,info->w[info->numi]);
    tx=0.0;
    for(tn=istart;tn<=iend;tn++)
    {
      tx+=gsl_vector_get(info->bs,tn-istart)*mat->xx[info->grid_base+info->grid[info->numi]+tn];
    }
    gsl_bspline_deriv_eval_nonzero(rr,(size_t)1,cg->b3_bs1,&istart,&iend,info->w[info->numi],cg->w3_bs1);
    tx1=0.0;
    for(tn=istart;tn<=iend;tn++)
    {
      tx1+=gsl_matrix_get(cg->b3_bs1,tn-istart,1)*mat->xx[info->grid_base+info->grid[info->numi]+tn];
    }

    fprintf(ff,"%lf %.15le %.15le\n",rr,tx,tx1);
  }
  fclose(ff);
}




void OutputLinear(enum INTER_TYPE type, CG_PAR *const cg, MAT_DAT *const mat, TYPE_INFO *const info, FILE *const br_out)
{
  int i_mesh,k;
  char name_tmp[100],name_tmp1[100],name_tmp2[100];
  double rr,tx;
  FILE *ff,*ff_b,*ff_d;
  double spaceo,ep;
  InitialOutput(type,cg,info,name_tmp,&spaceo);

  info->numi--;
  ff=OpenFile(name_tmp,"w");
  double max;
  if(type==non && fabs(info->rmax[info->mi]-cg->cutoff)<VERYSMALL_F) max=info->rmax[info->mi]+VERYSMALL_F;
  else max=info->rmax[info->mi];
  for(rr=((int)(info->rmin[info->mi]/spaceo)+1)*spaceo;rr<max;rr+=spaceo)
  {
    i_mesh=(int)((rr-info->rmin[info->mi])/info->space);
    ep=(rr-i_mesh*info->space-info->rmin[info->mi])/info->space;
    tx=mat->xx[info->grid_base+info->grid[info->numi]+i_mesh]*(1.0-ep)+mat->xx[info->grid_base+info->grid[info->numi]+i_mesh+1]*ep;
    fprintf(ff,"%lf %.15le\n",rr,tx);
  }
  fclose(ff);

  if(cg->out_sp==1)
  {
    sprintf(name_tmp1,"%s.b",name_tmp);
    ff_b=OpenFile(name_tmp1,"w");
    for(int i=info->grid_base+info->grid[info->numi];i<info->grid_base+info->grid[info->numi+1];i++)
    {
      fprintf(ff_b,"%lf %.15le\n",info->rmin[info->mi]+info->space*(i-info->grid_base-info->grid[info->numi]),mat->xx[i]);
    }
    fclose(ff_b);

    if(cg->out_dd==1)
    {
      sprintf(name_tmp2,"%s.dd",name_tmp);
      ff_d=OpenFile(name_tmp2,"w");
      for(int i=info->grid_base+info->grid[info->numi];i<info->grid_base+info->grid[info->numi+1];i++)
      {
        fprintf(ff_d,"%lf %.15le\n",info->rmin[info->mi]+info->space*(i-info->grid_base-info->grid[info->numi]),mat->dd[i]);
      }
      fclose(ff_d);
    }

  }

}

