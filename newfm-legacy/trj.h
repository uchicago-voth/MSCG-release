void FlagError()
{
  printf("Wrong flag!\n");
  exit (0);
}

void CheckSuffix(const char *name, const char *suffix)
{
  char temp;
  int len,pos,i;
  len=strlen(name);
  pos=-1;
  for(i=0;i<len;i++)
  {
    if(name[i]=='.') pos=i;
  }
  if(pos<0) FlagError();
  if(strcmp(&name[pos+1],suffix)!=0) FlagError();
}

void ReadAuguments(const int num_arg, char **arg, FR_DAT *const fr)
{
  if(num_arg!=3 && num_arg!=5) FlagError();
  else if(num_arg==3)
  {
    if(strcmp(arg[1],"-f")!=0) FlagError();   
    sscanf(arg[2],"%s",fr->trj_name);
    CheckSuffix(arg[2],"trr");
    fr->trj_type=0;
  }
  else if(num_arg==5)
  {
    if(strcmp(arg[1],"-f")!=0 || strcmp(arg[3],"-f1")!=0) FlagError();
    sscanf(arg[2],"%s",fr->trj_name);
    sscanf(arg[4],"%s",fr->trj_name1);
    CheckSuffix(arg[2],"xtc");
    CheckSuffix(arg[4],"xtc");
    fr->trj_type=1;
    fr->trj_type=1;
  }
}   

void ReadFirstFrameTrr(CG_PAR *const cg, FR_DAT *const fr)
{
  int junk_i;
  real junk_f;

  read_trr_natoms(fr->trj_name,&fr->natoms);
  fr->x=(rvec *)MemAlloc((fr->natoms+1)*sizeof(rvec));
  fr->f=(rvec *)MemAlloc((fr->natoms+1)*sizeof(rvec));

  if(cg->tolcgn!=fr->natoms) {printf("Warning: Number of CG sites defined in top.in is not consistent with trajectory!\n");}


  fr->fp=xdrfile_open(fr->trj_name,"r");
  if(read_trr(fr->fp,fr->natoms,&fr->step,&fr->time,&junk_f,fr->box,fr->x,NULL,fr->f)!=exdrOK)
  {printf("Can not read the first frame!\n"); exit(1);}
  int i;
  for(i=0;i<3;i++) fr->box_dim2[i]=fr->box[i][i]*0.5;

}


void ReadFirstFrameXtc(CG_PAR *const cg, FR_DAT *const fr)
{
  int junk_i;
  real junk_f;
  int n_atoms;

  read_xtc_natoms(fr->trj_name,&fr->natoms);
  read_xtc_natoms(fr->trj_name1,&n_atoms);
  if(fr->natoms!=n_atoms) {printf("Atom numbers are not consistent between two xtc files!\n"); exit(1);}

  fr->x=(rvec *)MemAlloc((fr->natoms+1)*sizeof(rvec));
  fr->f=(rvec *)MemAlloc((fr->natoms+1)*sizeof(rvec));

  if(cg->tolcgn!=fr->natoms) {printf("Warning: Number of CG sites defined in top.in is not consistent with trajectory!\n");}


  fr->fp=xdrfile_open(fr->trj_name,"r");
  fr->fp1=xdrfile_open(fr->trj_name1,"r");

  if((read_xtc(fr->fp1,fr->natoms,&fr->step,&fr->time,fr->box,fr->f,&junk_f)
     !=exdrOK) ||
  (read_xtc(fr->fp,fr->natoms,&junk_i,&junk_f,fr->box,fr->x,&junk_f)
  !=exdrOK))
  {printf("Can not read the first frame!\n"); exit(1);}
  int i;
  for(i=0;i<3;i++) fr->box_dim2[i]=fr->box[i][i]*0.5;

}


int ReadNextFrameTrr(FR_DAT *const fr)
{
  int junk_i;
  real junk_f;
  int ret;

  if(read_trr(fr->fp,fr->natoms,&fr->step,&fr->time,&junk_f,fr->box,fr->x,NULL,fr->f)==exdrOK) ret=1;
  else ret=0;
  int i;

  for(i=0;i<3;i++) fr->box_dim2[i]=fr->box[i][i]*0.5;
  return ret;
}
  

int ReadNextFrameXtc(FR_DAT *const fr)
{
  int junk_i;
  real junk_f;
  int ret;

  if((read_xtc(fr->fp1,fr->natoms,&fr->step,&fr->time,fr->box,fr->f,&junk_f)
     ==exdrOK) &&
  (read_xtc(fr->fp,fr->natoms,&junk_i,&junk_f,fr->box,fr->x,&junk_f)
  ==exdrOK)) ret=1;
  else ret=0;

  int i;
  for(i=0;i<3;i++) fr->box_dim2[i]=fr->box[i][i]*0.5;
  return ret;
}


void ReadVirial(FR_DAT *const fr, CG_PAR *const cg)
{
  //reading virial information
  //b_p is the RHS of eq. (12) in JCP,123,134105,2005
  double junk;
  fr->b_p=(double *)MemAlloc(cg->fr_n*sizeof(double));
  FILE *pp=OpenFile("p_con.in","r");
  int i;

  for(i=0;i<cg->start_fr-1;i++)
  {
    fscanf(pp,"%lf",&junk);
  }
  for(i=0;i<cg->fr_n;i++)
  {
    fscanf(pp,"%lf",fr->b_p+i);
  }
  fclose(pp);

}

//two cell link-list subroutines make_map and link_list
//refer to "Computer Simulation of Liquids" or other text books
void MakeMap(int *const map, int m0, int m1, int m2)
{
  int ix,iy,iz,imap;
  int m01=m0*m1;
  for(iz=0;iz<m2;iz++)
  {
    for(iy=0;iy<m1;iy++)
    {
      for(ix=0;ix<m0;ix++)
      {
        imap=(ix+iy*m0+iz*m01)*13;
        map[imap]=(ix+1)%m0+iy*m0+iz*m01;
        map[imap+1]=(ix+1)%m0+(iy+1)%m1*m0+iz*m01;
        map[imap+2]=(ix+1)%m0+iy*m0+(iz+1)%m2*m01;
        map[imap+3]=(ix+1)%m0+(iy-1+m1)%m1*m0+iz*m01;
        map[imap+4]=(ix+1)%m0+iy*m0+(iz-1+m2)%m2*m01;
        map[imap+5]=(ix+1)%m0+(iy+1)%m1*m0+(iz-1+m2)%m2*m01;
        map[imap+6]=(ix+1)%m0+(iy-1+m1)%m1*m0+(iz+1)%m2*m01;
        map[imap+7]=(ix+1)%m0+(iy-1+m1)%m1*m0+(iz-1+m2)%m2*m01;
        map[imap+8]=(ix+1)%m0+(iy+1)%m1*m0+(iz+1)%m2*m01;
        map[imap+9]=ix+(iy+1)%m1*m0+(iz+1)%m2*m01;
        map[imap+10]=ix+(iy+1)%m1*m0+iz*m01;
        map[imap+11]=ix+(iy+1)%m1*m0+(iz-1+m2)%m2*m01;
        map[imap+12]=ix+iy*m0+(iz+1)%m2*m01;
      }
    }
  }
}

void MakeMap_3(int *const map, int m0, int m1, int m2)
{
  int ix,iy,iz,i,j,k,imap;
  int m01=m0*m1;
  for(iz=0;iz<m2;iz++)
  {
    for(iy=0;iy<m1;iy++)
    {
      for(ix=0;ix<m0;ix++)
      {

        imap=(ix+iy*m0+iz*m01)*26;          
        for(i=-1;i<=1;i++)
        {
          for(j=-1;j<=1;j++)
          {
            for(k=-1;k<=1;k++)
            {
              if(i!=0 || j!=0 || k!=0)
              {
                map[imap]=(ix+i+m0)%m0+(iy+j+m1)%m1*m0+(iz+k+m2)%m2*m01;
                imap++;
              }
            }
          }
        }

      }
    }
  }
}


void CalCell(int *m0, int *m1, int *m2, real *cellx, real *celly, real *cellz, int *ncell, int **map, int **head, int **list, const double cutoff, FR_DAT *const fr)
{
  if(cutoff>0.5*fr->box[0][0] || cutoff>0.5*fr->box[1][1] || cutoff>0.5*fr->box[2][2])
  { printf("Cutoff is larger than the half box length!\n"); exit(1);}
  //set up cell linked-list
  *m0=(int)(fr->box[0][0]/cutoff);
  *m1=(int)(fr->box[1][1]/cutoff);
  *m2=(int)(fr->box[2][2]/cutoff);
  if((*m0)<3 || (*m1)<3 || (*m2)<3)
  {
    (*m0)=(*m1)=(*m2)=3;
    (*cellx)=(*celly)=(*cellz)=-1.0;
  }
  else
  {
    *cellx=fr->box[0][0]/(*m0);
    *celly=fr->box[1][1]/(*m1);
    *cellz=fr->box[2][2]/(*m2);
  }
  *ncell=(*m0)*(*m1)*(*m2);
  *head=(int *)MemAlloc((*ncell)*sizeof(int));
  *list=(int *)MemAlloc(fr->natoms*sizeof(int));
}


void InitialLinkedList(CG_PAR *const cg, FR_DAT *const fr)
{
  CalCell(&fr->m0,&fr->m1,&fr->m2,&fr->cellx,&fr->celly,&fr->cellz,&fr->ncell,&fr->map,&fr->head,&fr->list,cg->cutoff,fr);
  fr->map=(int *)MemAlloc((fr->ncell*13)*sizeof(int));
  MakeMap(fr->map,fr->m0,fr->m1,fr->m2);
}

void InitialLinkedList_3(CG_PAR *const cg, FR_DAT *const fr)
{
  double cut_off;
  int i;
  cut_off=cg->cutoff_3i[0];
  for(i=1;i<cg->tol3;i++) {if(cg->cutoff_3i[i]>cut_off) cut_off=cg->cutoff_3i[i];}
  CalCell(&fr->m0_3,&fr->m1_3,&fr->m2_3,&fr->cellx_3,&fr->celly_3,&fr->cellz_3,&fr->ncell_3,&fr->map_3,&fr->head_3,&fr->list_3,cut_off,fr);
  fr->map_3=(int *)MemAlloc((fr->ncell_3*26)*sizeof(int));
  MakeMap_3(fr->map_3,fr->m0_3,fr->m1_3,fr->m2_3);
}


void CalList(const real cellx, const real celly, const real cellz, const int m0, const int m1, const int ncell, int *const head, int *const list,FR_DAT *const fr)
{
  int i,icell;
  int m01=m0*m1;
  double cellx1,celly1,cellz1;
  cellx1=1.0/cellx;
  celly1=1.0/celly;
  cellz1=1.0/cellz; 


  if(cellx>0)
  {
    for(i=0;i<ncell;i++)
    {
     head[i]=-1;
    }

    for(i=0;i<fr->natoms;i++)
    {
      icell=(int)(fr->x[i][0]*cellx1)+(int)(fr->x[i][1]*celly1)*m0+(int)(fr->x[i][2]*cellz1)*m01;
      list[i]=head[icell];
      head[icell]=i;
    }
  }
  else
  {
    for(i=1;i<ncell;i++)
    {
      head[i]=-1;
    }

    head[0]=0;
    for(i=0;i<fr->natoms-1;i++)
    {
      list[i]=i+1;
    }
    list[fr->natoms-1]=-1;
  }

}

void LinkedList(FR_DAT *const fr)
{
  CalList(fr->cellx,fr->celly,fr->cellz,fr->m0,fr->m1,fr->ncell,fr->head,fr->list,fr);
}

void LinkedList_3(FR_DAT *const fr)
{
  CalList(fr->cellx_3,fr->celly_3,fr->cellz_3,fr->m0_3,fr->m1_3,fr->ncell_3,fr->head_3,fr->list_3,fr);
}




