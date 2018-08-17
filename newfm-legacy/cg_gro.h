#define MAX_MOL_NAME 5

void MakeGro(CG_PAR *const cg, FR_DAT *const fr)
{

  int i,j,k,l,n,m,ki,kk,kj;
  FILE *ff,*pp;    //input and output files
  char name_tmp[20];

  char buff[1000];
  int tmp=1;
  do {
       printf("\r%d frames have been read.",tmp); fflush(stdout);
       tmp++;
     } while((*cg->read_next)(fr));

  FILE *mol;
  int seg,*mol_num,*mol_size,n_atom,n_mol;
  char **mol_name;
  mol=OpenFile("mol.in","r");  
  fscanf(mol,"%d",&seg);
  mol_num=(int *)malloc((unsigned)seg*sizeof(int));
  mol_name=(char **)malloc((unsigned)seg*sizeof(char *));
  mol_size=(int *)malloc((unsigned)seg*sizeof(int));
  for(i=0;i<seg;i++)
  {
    mol_name[i]=(char *)malloc((unsigned)MAX_MOL_NAME*sizeof(char));
  }
  for(i=0;i<seg;i++)
  {
    fscanf(mol,"%s%d%d",mol_name[i],mol_num+i,mol_size+i);
  }
  fclose(mol);

    FILE *gro;
    gro=OpenFile("cg.gro","w");
    fprintf(gro,"CG_GRO\n%d\n",cg->tolcgn);
    n_atom=0; n_mol=0;
    for(k=0;k<seg;k++)
    {
      for(l=0;l<mol_num[k];l++)
      {
        for(m=0;m<mol_size[k];m++)
        {
          fprintf(gro,"%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",n_mol+1,mol_name[k],cg->name[cg->cgtype[n_atom]-1],n_atom+1,fr->x[n_atom][0],fr->x[n_atom][1],fr->x[n_atom][2]);
          n_atom++;
        }
        n_mol++;
      }
    }

    fprintf(gro,"%10.5f %10.5f %10.5f\n",fr->box[0][0],fr->box[1][1],fr->box[2][2]);
    fclose(gro);

  free(fr->x);
  free(fr->f);
  for(i=0;i<cg->tolcgn;i++)
  {
    free(cg->blist_i[i]);
    free(cg->alist_i[i]);
    free(cg->dlist_i[i]);
  }

  free(cg->cgtype);
  
  
  free(cg->blist_i);
  free(cg->alist_i);
  free(cg->blist_n);
  free(cg->alist_n);
  free(cg->dlist_n);
  free(cg->dlist_i);

  //free memory
  for(i=0;i<cg->tolcgtype;i++) free(cg->name[i]);
  free(cg->name);

}

