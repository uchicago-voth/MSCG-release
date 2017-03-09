void AllocTop(CG_PAR *const cg)
{
  int tn;
  cg->cgtype=(int *)MemAlloc(cg->tolcgn*sizeof(int));

  tn=TOL_TWO(cg->tolcgtype);
  cg->bontype=(int *)MemAlloc(tn*sizeof(int));
  cg->blist_n=(int *)MemAlloc(cg->tolcgn*sizeof(int));
  cg->blist_i=(int **)MemAlloc(cg->tolcgn*sizeof(int *));
  int i;
  for(i=0;i<cg->tolcgn;i++)
  {
    cg->blist_i[i]=(int *)MemAlloc(cg->max_bon*sizeof(int));
  }

  tn=TOL_THREE(cg->tolcgtype);
  cg->angtype=(int *)MemAlloc(tn*sizeof(int));
  cg->alist_n=(int *)MemAlloc(cg->tolcgn*sizeof(int));
  cg->alist_i=(int **)MemAlloc(cg->tolcgn*sizeof(int *));
  for(i=0;i<cg->tolcgn;i++)
  {
    cg->alist_i[i]=(int *)MemAlloc(2*cg->max_ang*sizeof(int));
  }

  tn=TOL_FOUR(cg->tolcgtype);
  cg->dihtype=(int *)MemAlloc(tn*sizeof(int));
  cg->dlist_n=(int *)MemAlloc(cg->tolcgn*sizeof(int));
  cg->dlist_i=(int **)MemAlloc(cg->tolcgn*sizeof(int *));
  for(i=0;i<cg->tolcgn;i++)
  {
    cg->dlist_i[i]=(int *)MemAlloc(3*cg->max_dih*sizeof(int));
  }

  cg->name=(char **)MemAlloc(cg->tolcgtype*sizeof(char *));
  for(i=0;i<cg->tolcgtype;i++)
  {
    cg->name[i]=(char *)MemAlloc((MAX_TYPENAME_LEN+1)*sizeof(char));
  }
}

void FreeTop (CG_PAR *const cg)
{
  free(cg->cgtype);

  int i;
  for(i=0;i<cg->tolcgn;i++)
  {
   free(cg->blist_i[i]);
  }
  free(cg->bontype);
  free(cg->blist_n);
  free(cg->blist_i);

  for(i=0;i<cg->tolcgn;i++)
  {
    free(cg->alist_i[i]);
  }
  free(cg->angtype);
  free(cg->alist_n);
  free(cg->alist_i);

  for(i=0;i<cg->tolcgn;i++)
  {
    free(cg->dlist_i[i]);
  }
  free(cg->dihtype);
  free(cg->dlist_n);
  free(cg->dlist_i);

  for(i=0;i<cg->tolcgtype;i++)
  {
    free(cg->name[i]);
  }
  free(cg->name);

}

void ReadErrorTop(const int line)
{
  printf("Wrong format in top.in:line %d!\n",line);
  exit(1);
}


void ReadMol(CG_PAR *const mol, CG_PAR *const cg, FILE *const top_in, int *line)
{
  int i,j;
  int exc;
  char buff[100],par[50];
  fgets(buff,100,top_in); (*line)++;
  sscanf(buff,"%s%d%d",par,&mol->tolcgn,&exc);
  if(strcmp(par,"mol")!=0) ReadErrorTop(*line);
  mol->tolcgtype=0;
  AllocTop(mol);

  fgets(buff,100,top_in); (*line)++;
  sscanf(buff,"%s",par);
  if(strcmp(par,"sitetypes")!=0) ReadErrorTop(*line);
  for(i=0;i<mol->tolcgn;i++)
  {
    fgets(buff,100,top_in); (*line)++;
    sscanf(buff,"%d",&mol->cgtype[i]);
  }

  for(i=0;i<mol->tolcgn;i++)
  {
    mol->blist_n[i]=0;
    mol->alist_n[i]=0;
    mol->dlist_n[i]=0;
  }

  int tolb;
  int tn,tn1,tn2,tn3;
  int cg_t,cg_t1,cg_t2,cg_t3;
  int m;
  fgets(buff,100,top_in); (*line)++;
  sscanf(buff,"%s%d",par,&tolb);
  if(strcmp(par,"bonds")!=0) ReadErrorTop(*line);

  for(i=0;i<tolb;i++)
  {
    fgets(buff,100,top_in); (*line)++;
    sscanf(buff,"%d%d",&tn,&tn1);
    tn--; tn1--;
    cg_t=mol->cgtype[tn];
    cg_t1=mol->cgtype[tn1];
    m=TwoBodyNum(cg_t,cg_t1,cg->tolcgtype);
    cg->bontype[m]=1;
    mol->blist_i[tn][mol->blist_n[tn]]=tn1;
    mol->blist_n[tn]++;
    mol->blist_i[tn1][mol->blist_n[tn1]]=tn;
    mol->blist_n[tn1]++;
  }
  
  if(exc==-1)
  {
    int tola;
    fgets(buff,100,top_in); (*line)++;
    sscanf(buff,"%s%d",par,&tola);
    if(strcmp(par,"angles")!=0) ReadErrorTop(*line);
     
    for(i=0;i<tola;i++)
    {
      fgets(buff,100,top_in); (*line)++;
      sscanf(buff,"%d%d%d",&tn,&tn1,&tn2);
      tn--;tn1--;tn2--;
      cg_t=mol->cgtype[tn];
      cg_t1=mol->cgtype[tn1];
      cg_t2=mol->cgtype[tn2];
      m=ThreeBodyNum(cg_t,cg_t1,cg_t2,cg->tolcgtype);
      cg->angtype[m]=1;
      mol->alist_i[tn1][2*mol->alist_n[tn1]]=tn;
      mol->alist_i[tn1][2*mol->alist_n[tn1]+1]=tn2;
      mol->alist_n[tn1]++;
      mol->alist_i[tn2][2*mol->alist_n[tn2]]=tn;
      mol->alist_i[tn2][2*mol->alist_n[tn2]+1]=tn1;
      mol->alist_n[tn2]++;
    }

    int told;
    fgets(buff,100,top_in); (*line)++;
    sscanf(buff,"%s%d",par,&told);
    if(strcmp(par,"dihedrals")!=0) ReadErrorTop(*line);
    for(i=0;i<told;i++)
    {
      fgets(buff,100,top_in); (*line)++;
      sscanf(buff,"%d%d%d%d",&tn,&tn1,&tn2,&tn3);
      tn--;tn1--;tn2--;tn3--;
      cg_t=mol->cgtype[tn];
      cg_t1=mol->cgtype[tn1];
      cg_t2=mol->cgtype[tn2];
      cg_t3=mol->cgtype[tn3];
      m=FourBodyNum(cg_t,cg_t1,cg_t2,cg_t3,cg->tolcgtype);
      cg->dihtype[m]=1;
      mol->dlist_i[tn2][3*mol->dlist_n[tn2]]=tn;
      mol->dlist_i[tn2][3*mol->dlist_n[tn2]+1]=tn1;
      mol->dlist_i[tn2][3*mol->dlist_n[tn2]+2]=tn3;
      mol->dlist_n[tn2]++;
      mol->dlist_i[tn3][3*mol->dlist_n[tn3]]=tn1;
      mol->dlist_i[tn3][3*mol->dlist_n[tn3]+1]=tn;
      mol->dlist_i[tn3][3*mol->dlist_n[tn3]+2]=tn2;
      mol->dlist_n[tn3]++;
    }
    
  }


  else
  {
    int l,n;
    for(i=0;i<mol->tolcgn;i++)
    {
      n=0;
      for(j=0;j<mol->blist_n[i];j++)
      {
        tn=mol->blist_i[i][j];
        if(mol->cgtype[i]<=mol->cgtype[tn])
        {
          m=TwoBodyNum(mol->cgtype[i],mol->cgtype[tn],cg->tolcgtype);
          cg->bontype[m]=1;
        }
        if(exc==1) continue;
        for(l=0;l<mol->blist_n[tn];l++)
        {
          tn1=mol->blist_i[tn][l];
          if(tn1==i) continue;
          mol->alist_n[i]++;
          mol->alist_i[i][n]=tn;
          n++;
          mol->alist_i[i][n]=tn1;
          n++;
          cg_t=mol->cgtype[tn];
          cg_t1=mol->cgtype[i];
          cg_t2=mol->cgtype[tn1];
          m=ThreeBodyNum(cg_t,cg_t1,cg_t2,cg->tolcgtype);
          cg->angtype[m]=1;
        }
      }
    }

    int k;
    for(i=0;i<mol->tolcgn;i++)
    {
      if(exc<=2) continue;
      tn=0;
      for(j=0;j<mol->alist_n[i];j++)
      {
        for(k=0;k<mol->blist_n[mol->alist_i[i][2*j+1]];k++)
        {
          if(mol->blist_i[mol->alist_i[i][2*j+1]][k]==mol->alist_i[i][2*j]) continue;
          if(mol->blist_i[mol->alist_i[i][2*j+1]][k]==i) continue;
          mol->dlist_i[i][3*tn]=mol->alist_i[i][2*j];
          mol->dlist_i[i][3*tn+1]=mol->alist_i[i][2*j+1];
          mol->dlist_i[i][3*tn+2]=mol->blist_i[mol->alist_i[i][2*j+1]][k];
          tn1=FourBodyNum(mol->cgtype[mol->dlist_i[i][3*tn]],mol->cgtype[mol->dlist_i[i][3*tn+1]],mol->cgtype[i],mol->cgtype[mol->dlist_i[i][3*tn+2]],cg->tolcgtype);
          cg->dihtype[tn1]=1;
          tn++;
        }
      }
      mol->dlist_n[i]=tn;
      
    }

  }

}


void ReadTop(CG_PAR *const cg)
{
  char buff[100],par[50];
  FILE *top_in;
  top_in=OpenFile("top.in","r");
  int line=0;

  fgets(buff,100,top_in); line++;
  sscanf(buff,"%s%d",par,&cg->tolcgn);
  if(strcmp(par,"cgsites")!=0) ReadErrorTop(line);
  fgets(buff,100,top_in); line++;
  sscanf(buff,"%s%d",par,&cg->tolcgtype);
  if(strcmp(par,"cgtypes")!=0) ReadErrorTop(line);

  AllocTop(cg);

  int i;
  for(i=0;i<cg->tolcgtype;i++)
  {
    fgets(buff,100,top_in); line++;
    sscanf(buff,"%s",cg->name[i]);
  }

  int tbtype;
  int *tb_i,*tb_j,*tb_k;
  if(cg->use_3b>0)
  {
    fgets(buff,100,top_in); line++;
    sscanf(buff,"%s%d",par,&tbtype);
    if(strcmp(par,"threebody")!=0) ReadErrorTop(line);
    cg->cos_theta0i=(double *)MemAlloc(tbtype*sizeof(double));
    cg->cutoff_3i=(double *)MemAlloc(tbtype*sizeof(double));


    tb_i=(int *)MemAlloc(tbtype*sizeof(int));
    tb_j=(int *)MemAlloc(tbtype*sizeof(int));
    tb_k=(int *)MemAlloc(tbtype*sizeof(int));

    for(i=0;i<tbtype;i++)
    {
      fgets(buff,100,top_in); line++;
      sscanf(buff,"%d%d%d%lf%lf",tb_i+i,tb_j+i,tb_k+i,cg->cos_theta0i+i,cg->cutoff_3i+i);
      tb_i[i]--; tb_j[i]--; tb_k[i]--;
    }

    cg->tb_n=(int *)MemAlloc(cg->tolcgtype*sizeof(int));
    cg->tb_list=(int **)MemAlloc(cg->tolcgtype*sizeof(int *));
    for(i=0;i<cg->tolcgtype;i++) cg->tb_n[i]=0;

    for(i=0;i<tbtype;i++)
    {
      cg->tb_n[tb_i[i]]++;
    }

    for(i=0;i<cg->tolcgtype;i++)
    {
      cg->tb_list[i]=(int *)MemAlloc(cg->tb_n[i]*2*sizeof(int));
      cg->tb_n[i]=0;
    }

    for(i=0;i<tbtype;i++)
    {
      cg->tb_list[tb_i[i]][cg->tb_n[tb_i[i]]*2]=tb_j[i]+1;
      cg->tb_list[tb_i[i]][cg->tb_n[tb_i[i]]*2+1]=tb_k[i]+1;
      cg->tb_n[tb_i[i]]++;
    }
    cg->tol3=tbtype;
    free(tb_i); free(tb_j); free(tb_k);

  }

  fgets(buff,100,top_in); line++;
  int moltype;
  sscanf(buff,"%s%d",par,&moltype);
  if(strcmp(par,"moltypes")!=0) ReadErrorTop(line);

  CG_PAR *mol;
  mol=(CG_PAR *)MemAlloc(moltype*sizeof(CG_PAR));

  for(i=0;i<moltype;i++)
  {
    mol[i].max_bon=cg->max_bon;
    mol[i].max_ang=cg->max_ang;
    mol[i].max_dih=cg->max_dih;
    AllocTop(&mol[i]);
  }

  for(i=0;i<TOL_TWO(cg->tolcgtype);i++)
  {
    cg->bontype[i]=0;
  }
  for(i=0;i<TOL_THREE(cg->tolcgtype);i++)
  {
    cg->angtype[i]=0;
  }
  for(i=0;i<TOL_FOUR(cg->tolcgtype);i++)
  {
    cg->dihtype[i]=0;
  }

  int j,k,l;
  for(i=0;i<moltype;i++)
  {
    ReadMol(&mol[i],cg,top_in,&line);
  }


  int block,mol_ndx,mol_num,offset=0;
  fgets(buff,100,top_in); line++;
  sscanf(buff,"%s%d",par,&block);
  if(strcmp(par,"system")!=0) ReadErrorTop(line);
  for(i=0;i<block;i++)
  {
    fgets(buff,100,top_in); line++;
    sscanf(buff,"%d%d",&mol_ndx,&mol_num);
    mol_ndx--;
    for(j=0;j<mol_num;j++)
    {
      for(k=0;k<mol[mol_ndx].tolcgn;k++)
      {
        cg->cgtype[k+offset]=mol[mol_ndx].cgtype[k];
        cg->blist_n[k+offset]=mol[mol_ndx].blist_n[k];
        for(l=0;l<mol[mol_ndx].blist_n[k];l++)
        {
          cg->blist_i[k+offset][l]=mol[mol_ndx].blist_i[k][l]+offset;
        }
        cg->alist_n[k+offset]=mol[mol_ndx].alist_n[k];
        for(l=0;l<2*mol[mol_ndx].alist_n[k];l++)
        {
          cg->alist_i[k+offset][l]=mol[mol_ndx].alist_i[k][l]+offset;
        }
        cg->dlist_n[k+offset]=mol[mol_ndx].dlist_n[k];
        for(l=0;l<3*mol[mol_ndx].dlist_n[k];l++)
        {
          cg->dlist_i[k+offset][l]=mol[mol_ndx].dlist_i[k][l]+offset;
        }
      }
      offset+=mol[mol_ndx].tolcgn;
    }
  }

  fclose(top_in);

  for(i=0;i<moltype;i++)
  {
    FreeTop(&mol[i]);
  }
  free(mol);
}

void ReadOneRange(double *const rmin, double *const rmax, int *const num, const int tol, const int n_body, const double cutoff, const double space, FILE * range_in, int *const tol1, int *const tol2, int d_r, double spaceo, int ba_type) 
{
  int tn,tn1,i,j;
  char junk[10];
  double diff;
  tn=0; tn1=0;
  for(i=0;i<tol;i++)
  {
    for(j=0;j<n_body;j++) fscanf(range_in,"%s",junk);      
    fscanf(range_in,"%lf%lf\n",rmin+i,rmax+i);
    if(fabs(rmin[i]+1.0)<VERYSMALL_F) num[i]=0;
    else if(fabs(rmin[i]+2.0)<VERYSMALL_F) {tn--; num[i]=tn;}
    else
    {
      tn1++;
      num[i]=tn1;
/*
      diff=((int)((rmax[i]-rmin[i])/space)+1)*space-(rmax[i]-rmin[i]);
      diff*=0.5;
      //rmin[i]=rmin[i]-VERYSMALL_F;
      rmax[i]=rmax[i]+diff;
      //rmax[i]=rmax[i]+VERYSMALL_F;
      rmin[i]=rmax[i]-((int)((rmax[i]-rmin[i])/space)+1)*space;
*/
      if(ba_type==1)
      {
        rmin[i]=floor(rmin[i]/spaceo+0.5)*spaceo;
        if(rmin[i]<0.0) rmin[i]=0.0;
        rmax[i]=rmin[i]+floor((rmax[i]-rmin[i])/space+0.5)*space;
        if(d_r==-1 && fabs(rmax[i]-cutoff-space)<VERYSMALL_F) rmax[i]-=space;
      }
      else if(ba_type==0)
      {
        rmax[i]=(((int)(rmax[i]/spaceo))+1)*spaceo;
        rmin[i]=rmax[i]-((int)((rmax[i]-rmin[i])/space)+1)*space;
      }

      if (d_r==1) { rmax[i]-=180.0; rmin[i]-=180.0; }
    }
  }

  tn=0; tn1=0;
  for(i=0;i<tol;i++)
  {
    if(num[i]>0) tn++;
    else if(num[i]<0) tn1++;
  }
  *tol1=tn;
  *tol2=tn1;

}

int CalTol(int *const type, const int tol_n)
{
  int tol=0;
  int i;
  for(i=0;i<tol_n;i++)
  {
    if(type[i]==1) tol++;
  }
  return tol;
}

void CalM(int *const type, const int tol_n, int *const m)
{
  int i,tn=0;
  for(i=0;i<tol_n;i++)
  {
    if(type[i]==1) {m[tn]=i; tn++;}
  }
  free(type);
}

void CalGrid(double *const rmin, double *const rmax, int *const grid, int *const num, const int tol, const double space, const int ba_type, const int n_k)
{
  int tn=0;
  int grid_i;
  grid[0]=0;
  int i;
  for(i=0;i<tol;i++)
  {
    if(num[i]>0)
    {
      grid_i=floor((rmax[i]-rmin[i])/space+0.5)+1;
      //if(grid_i<2) {printf("Grid space too large for distance range!\n"); exit(0);};
      grid[tn+1]=grid[tn]+grid_i;
      if(ba_type==0) grid[tn+1]=grid[tn+1]+n_k-2;
      tn++;
    }
  }

}

void ReadRange(CG_PAR *const cg)
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


  FILE *range_in;
  range_in=OpenFile("rmin.in","r");

  cg->tollj1=0;
  cg->tolb1=0;
  cg->tola1=0;
  cg->told1=0;
  cg->tollj2=0;
  cg->tolb2=0;
  cg->tola2=0;
  cg->told2=0;


  ReadOneRange(cg->rmin,cg->rmax,cg->ljnum,cg->tollj,2,cg->cutoff,cg->space,range_in,&cg->tollj1,&cg->tollj2,-1,cg->spaceo,cg->ba_type);
  fclose(range_in);

  range_in=OpenFile("rmin_b.in","r");

  ReadOneRange(cg->rmin_b,cg->rmax_b,cg->bonnum,cg->tolb,2,cg->cutoff,cg->space_b,range_in,&cg->tolb1,&cg->tolb2,0,cg->spaceo_b,cg->ba_type);

  int d_r;
  ReadOneRange(cg->rmin_a,cg->rmax_a,cg->angnum,cg->tola,3,cg->cutoff,cg->space_a,range_in,&cg->tola1,&cg->tola2,0,cg->spaceo_a,cg->ba_type);
  if(cg->dih_type==1) d_r=1;
  else d_r=0;
  ReadOneRange(cg->rmin_d,cg->rmax_d,cg->dihnum,cg->told,4,cg->cutoff,cg->space_d,range_in,&cg->told1,&cg->told2,d_r,cg->spaceo_d,cg->ba_type);
  fclose(range_in);
 

  cg->grid_n=(int *)MemAlloc((cg->tollj1+1)*sizeof(int));
  cg->grid_bn=(int *)MemAlloc((cg->tolb1+1)*sizeof(int));
  cg->grid_an=(int *)MemAlloc((cg->tola1+1)*sizeof(int));
  cg->grid_dn=(int *)MemAlloc((cg->told1+1)*sizeof(int));

  CalGrid(cg->rmin,cg->rmax,cg->grid_n,cg->ljnum,cg->tollj,cg->space,cg->ba_type,cg->n_k);
  CalGrid(cg->rmin_b,cg->rmax_b,cg->grid_bn,cg->bonnum,cg->tolb,cg->space_b,cg->ba_type,cg->nb_k);
  CalGrid(cg->rmin_a,cg->rmax_a,cg->grid_an,cg->angnum,cg->tola,cg->space_a,cg->ba_type,cg->na_k);
  CalGrid(cg->rmin_d,cg->rmax_d,cg->grid_dn,cg->dihnum,cg->told,cg->space_d,cg->ba_type,cg->nd_k);   

  int i,j,tmp;
  if(cg->use_3b>0)
  {
    cg->tb_m=(int *)MemAlloc(cg->tol3*sizeof(int));
    cg->tbnum=(int *)MemAlloc(cg->tol3*sizeof(int));

    cg->rmin_3=(double *)MemAlloc(cg->tol3*sizeof(double));
    cg->rmax_3=(double *)MemAlloc(cg->tol3*sizeof(double));
    cg->grid_3n=(int *)MemAlloc((cg->tol3+1)*sizeof(int));

    for(i=0;i<cg->tol3;i++)
    {
      cg->tbnum[i]=i+1; //no input 3-body
      cg->rmin_3[i]=0.0;
      cg->rmax_3[i]=180.0;
    }

    cg->grid_3n[0]=0;
    if(cg->use_3b==3) for(i=1;i<cg->tol3+1;i++) cg->grid_3n[i]=i;
    else
    {
      if(cg->ba_type==0)
      for(i=1;i<cg->tol3+1;i++) cg->grid_3n[i]=i*(cg->n3_k-2+floor(180.0/cg->space_3+0.5)+1);
      else if(cg->ba_type==1)
      for(i=1;i<cg->tol3+1;i++) cg->grid_3n[i]=i*(floor(180.0/cg->space_3+0.5)+1);
    }

    tmp=0;
    for(i=0;i<cg->tolcgtype;i++)
    {
      for(j=0;j<cg->tb_n[i];j++)
      {
        cg->tb_m[tmp]=ThreeBodyNum(i+1,cg->tb_list[i][2*j],cg->tb_list[i][2*j+1],cg->tolcgtype);
        tmp++;
      }
    }
    free(cg->tb_n);
    for(i=0;i<cg->tolcgtype;i++) free(cg->tb_list[i]);
    free(cg->tb_list);

  }



}

void ReadErrorExt(const int line)
{
  printf("Wrong format in table.in:line %d!\n",line);
  exit(1);
}

void ReadErrorExt1(const int line)
{
  printf("Numbers of tabulated interactions from rmin.in/rmin_b.in and table.in are not consistent:line %d!\n",line);
  exit(1);
}

void ReadExt(CG_PAR *const cg)
{    
   char buff[100];
   FILE *es_in; //input for external slovent etc.
   int n_es; //number of external nonebonded interactions 
   int line=0;

   es_in=OpenFile("table.in","r");

   char par[50];

   int i,j, tn;
   int tn1;
   int m;
   int grid_es;
   int es_type;
   int offset;

   fgets(buff,100,es_in); line++;
   sscanf(buff,"%s%d%lf",par,&n_es,&cg->space_es);
   if(strcmp(par,"short_range")!=0) ReadErrorExt(line);
   if(n_es!=cg->tollj2) ReadErrorExt1(line);

   cg->es_co=(double **)MemAlloc(n_es*sizeof(double *));
   for(i=0;i<n_es;i++)
   {
     fgets(buff,100,es_in); line++;
     sscanf(buff,"%d%d",&tn,&tn1);
     m=TwoBodyNum(tn,tn1,cg->tolcgtype);
     offset=m;//SearchIntTable(cg->bon_m,m,cg->tollj);
     fgets(buff,100,es_in); line++;
     sscanf(buff,"%lf%lf",cg->rmin+offset,cg->rmax+offset);
     grid_es=floor((cg->rmax[offset]-cg->rmin[offset])/cg->space_es+0.5)+1;
     es_type=-cg->ljnum[offset]-1;
     cg->es_co[es_type]=(double *)MemAlloc(grid_es*sizeof(double));
     for(j=0;j<grid_es;j++)
     {
       fgets(buff,100,es_in); line++;
       sscanf(buff,"%lf",&cg->es_co[es_type][j]);
     }
   }

   fgets(buff,100,es_in); line++;
   sscanf(buff,"%s%d%lf",par,&n_es,&cg->spaceb_es);
   if(strcmp(par,"bond")!=0) ReadErrorExt(line);
   if(n_es!=cg->tolb2) ReadErrorExt1(line);

   cg->esb_co=(double **)MemAlloc(n_es*sizeof(double *));
   for(i=0;i<n_es;i++)
   {
     fgets(buff,100,es_in); line++;
     sscanf(buff,"%d%d",&tn,&tn1);
     m=TwoBodyNum(tn,tn1,cg->tolcgtype);
     offset=SearchIntTable(cg->bon_m,m,cg->tolb);
     fgets(buff,100,es_in); line++;
     sscanf(buff,"%lf%lf",cg->rmin_b+offset,cg->rmax_b+offset);
     grid_es=floor((cg->rmax_b[offset]-cg->rmin_b[offset])/cg->spaceb_es+0.5)+1;
     es_type=-cg->bonnum[offset]-1;
     cg->esb_co[es_type]=(double *)MemAlloc(grid_es*sizeof(double));
     for(j=0;j<grid_es;j++)
     {
       fgets(buff,100,es_in); line++;
       sscanf(buff,"%lf",&cg->esb_co[es_type][j]);
     }
   }

   n_es=0;
   fgets(buff,100,es_in); line++;
   sscanf(buff,"%s%d%lf",par,&n_es,&cg->spacea_es);
   if(strcmp(par,"angle")!=0) ReadErrorExt(line);
   if(n_es!=cg->tola2) ReadErrorExt1(line);

   cg->esa_co=(double **)MemAlloc(n_es*sizeof(double *));
   int tn2;
   for(i=0;i<n_es;i++)
   {
     fgets(buff,100,es_in); line++;
     sscanf(buff,"%d%d%d",&tn,&tn1,&tn2);
     m=ThreeBodyNum(tn,tn1,tn2,cg->tolcgtype);
     offset=SearchIntTable(cg->ang_m,m,cg->tola);
     fgets(buff,100,es_in); line++;
     sscanf(buff,"%lf%lf",cg->rmin_a+offset,cg->rmax_a+offset);
     grid_es=floor((cg->rmax_a[offset]-cg->rmin_a[offset])/cg->spacea_es+0.5)+1;
     es_type=-cg->angnum[offset]-1;
     cg->esa_co[es_type]=(double *)MemAlloc(grid_es*sizeof(double));
     for(j=0;j<grid_es;j++)
     {
       fgets(buff,100,es_in); line++;
       sscanf(buff,"%lf",&cg->esa_co[es_type][j]);
     }
   }

   n_es=0;
   fgets(buff,100,es_in); line++;
   sscanf(buff,"%s%d%lf",par,&n_es,&cg->spaced_es);
   if(strcmp(par,"dihedral")!=0) ReadErrorExt(line);
   if(n_es!=cg->told2) ReadErrorExt1(line);

   cg->esd_co=(double **)MemAlloc(n_es*sizeof(double *));
   int tn3;
   for(i=0;i<n_es;i++)
   {
     fgets(buff,100,es_in); line++;
     sscanf(buff,"%d%d%d%d",&tn,&tn1,&tn2,&tn3);
     m=FourBodyNum(tn,tn1,tn2,tn3,cg->tolcgtype);
     offset=SearchIntTable(cg->dih_m,m,cg->told);
     fgets(buff,100,es_in); line++;
     sscanf(buff,"%lf%lf",cg->rmin_d+offset,cg->rmax_d+offset);
     grid_es=floor((cg->rmax_d[offset]-cg->rmin_d[offset])/cg->spaced_es+0.5)+1;
     es_type=-cg->dihnum[offset]-1;
     cg->esd_co[es_type]=(double *)MemAlloc(grid_es*sizeof(double));
     for(j=0;j<grid_es;j++)
     {
       fgets(buff,100,es_in); line++;
       sscanf(buff,"%lf",&cg->esd_co[es_type][j]);
     }
   }

}

