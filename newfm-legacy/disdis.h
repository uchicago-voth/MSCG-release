//some vector functions from Gromacs 4, "vec.h"
#define XX 0
#define YY 1
#define ZZ 2
static inline void rvec_sub(const rvec a,const rvec b,rvec c)
{
  real x,y,z;
  
  x=a[XX]-b[XX];
  y=a[YY]-b[YY];
  z=a[ZZ]-b[ZZ];
  
  c[XX]=x;
  c[YY]=y;
  c[ZZ]=z;
}

static inline void cprod(const rvec a,const rvec b,rvec c)
{
  c[XX]=a[YY]*b[ZZ]-a[ZZ]*b[YY];
  c[YY]=a[ZZ]*b[XX]-a[XX]*b[ZZ];
  c[ZZ]=a[XX]*b[YY]-a[YY]*b[XX];
}

static inline real cos_angle(const rvec a,const rvec b)
{
  /* 
   *                  ax*bx + ay*by + az*bz
   * cos-vec (a,b) =  ---------------------
   *                      ||a|| * ||b||
   */
  real   cos;
  int    m;
  double aa,bb,ip,ipa,ipb,ipab; /* For accuracy these must be double! */

  ip=ipa=ipb=0.0;
  for(m=0; (m<DIM); m++) {              /* 18           */
    aa   = a[m];
    bb   = b[m];
    ip  += aa*bb;
    ipa += aa*aa;
    ipb += bb*bb;
  }
  ipab = ipa*ipb;
  if (ipab > 0)
    //cos = ip*invsqrt(ipab);             /*  7           */
    cos=ip/sqrt(ipab);
  else
    cos = 1;
                                        /* 25 TOTAL     */
  if (cos > 1.0)
    return  1.0;
  if (cos <-1.0)
    return -1.0;

  return cos;
}

static inline real iprod(const rvec a,const rvec b)
{
  return (a[XX]*b[XX]+a[YY]*b[YY]+a[ZZ]*b[ZZ]);
}
//end of Gromacs functions

typedef struct
{
  int grid_b;
  int grid_a;
  int grid_d;
  double *bond;
  double *angle;
  double *dihedral;
} DIS_DAT;

void AssignDisArray (CG_PAR *const cg, DIS_DAT *const dis)
{
  dis->grid_b=floor(cg->cutoff/cg->space_b+0.5)+1; 
  dis->grid_a=floor(180.0/cg->space_a+0.5)+1;
  cg->space_a=cg->space_a/180.0*M_PI;
  dis->grid_d=floor(360.0/cg->space_d+0.5)+1;
  cg->space_d=cg->space_d/180.0*M_PI;

  cg->tolb=CalTol(cg->bontype,TOL_TWO(cg->tolcgtype));
  cg->bon_m=(int *)MemAlloc(cg->tolb*sizeof(int));
  CalM(cg->bontype,TOL_TWO(cg->tolcgtype),cg->bon_m);

  cg->tola=CalTol(cg->angtype,TOL_THREE(cg->tolcgtype));
  cg->ang_m=(int *)MemAlloc(cg->tola*sizeof(int));
  CalM(cg->angtype,TOL_THREE(cg->tolcgtype),cg->ang_m);
  
  cg->told=CalTol(cg->dihtype,TOL_FOUR(cg->tolcgtype));
  cg->dih_m=(int *)MemAlloc(cg->told*sizeof(int));
  CalM(cg->dihtype,TOL_FOUR(cg->tolcgtype),cg->dih_m);

  dis->bond=(double *)MemAlloc(cg->tolb*dis->grid_b*sizeof(double));
  dis->angle=(double *)MemAlloc(cg->tola*dis->grid_a*sizeof(double));
  dis->dihedral=(double *)MemAlloc(cg->told*dis->grid_d*sizeof(double));
  int i;
  for(i=0;i<cg->tolb*dis->grid_b;i++) dis->bond[i]=0.0;
  for(i=0;i<cg->tola*dis->grid_a;i++) dis->angle[i]=0.0;
  for(i=0;i<cg->told*dis->grid_d;i++) dis->dihedral[i]=0.0;
}

void CalDistribution(CG_PAR *const cg, FR_DAT *const fr, DIS_DAT *const dis)
{

  /* for dihedral calculation */
  int *t1, *t2, *t3;
  rvec xi, xj, xk, xl, r_ki, r_ji, r_jl, mdih, ndih;
  double cos_phi, phi, iprd, sign;

  int i,j,k,l,kk,ki,kj,kl,i_mesh,ntt,n_ang,n_tmp,m;
  double rr,rx[3],rx1[3],rr1,tx,rr1_2,rr2_2,rr_2,theta;
  //skip some frames
  for(i=0;i<cg->start_fr-1;i++)
  {
    if((ReadNextFrameTrr)(fr)==0) {printf("\nCan not read a frame that needs to be skipped!\n"); exit(1);}
  }

  int read_stat=1;
  //start the main loop
  for(j=0;j<cg->fr_n;j++)
  {
    if(read_stat==0) {printf("\nCan not read frame %d!\n",j+1); exit(1);}
    //loop to obtain positions and forces of each CG site
    for(l=0;l<cg->tolcgn;l++)
    {
      //enforce pbc 
      SetPBC(l,fr->x,fr->box);
    }

    //loop for bonded interactions, everything is like what we've done for nonebonded
    for(k=0;k<cg->tolcgn;k++)
    {
      for(kk=0;kk<cg->blist_n[k];kk++)
      {
        l=cg->blist_i[k][kk];
        if(k>l) continue;
        rr=sqrt(CalDis (k,l,fr,rx));
        if(rr<cg->cutoff)
        {
          n_ang=TwoBodyNum(cg->cgtype[k],cg->cgtype[l],cg->tolcgtype);
          ntt=SearchIntTable(cg->bon_m,n_ang,cg->tolb);

          i_mesh=(int)((rr+0.5*cg->space_b)/cg->space_b);
          n_tmp=ntt*dis->grid_b+i_mesh;
          dis->bond[n_tmp]+=1.0;        
        } 
      } 
    }

    //loop for angles, again like what we've done
    for(k=0;k<cg->tolcgn;k++)
    {
      
      for(m=0;m<cg->alist_n[k];m++)
      {

        kk=cg->alist_i[k][2*m+1];
        if(kk<k) continue;
        ki=cg->alist_i[k][2*m];
        
        rr=sqrt(CalDis (ki,k,fr,rx));
        rr1=sqrt(CalDis(ki,kk,fr,rx1));

        rr_2=rr*rr; rr1_2=rr*rr1; rr2_2=rr1*rr1;
        tx=(rx[0]*rx1[0]+rx[1]*rx1[1]+rx[2]*rx1[2])/rr1_2;
        if(tx>1.0) tx=1.0;
        else if(tx<-1.0) tx=-1.0;

        theta=acos(tx); //printf("%lf\n",theta);

        n_ang=ThreeBodyNum(cg->cgtype[ki],cg->cgtype[k],cg->cgtype[kk],cg->tolcgtype);
        ntt=SearchIntTable(cg->ang_m,n_ang,cg->tola);
         
        i_mesh=(int)((theta+0.5*cg->space_a)/cg->space_a);
        n_tmp=dis->grid_a*ntt+i_mesh;
        dis->angle[n_tmp]+=1.0;
      }
    }                     


    //loop for dihedral
    for(k=0;k<cg->tolcgn;k++)
    {

      for(m=0;m<cg->dlist_n[k];m++)
      {
        
        kl=cg->dlist_i[k][3*m+2];
        if(kl<k) continue;
        ki=cg->dlist_i[k][3*m];
        kj=cg->dlist_i[k][3*m+1];

        rr=sqrt(CalDis (k,kl,fr,rx));
        n_ang=FourBodyNum(cg->cgtype[ki],cg->cgtype[kj],cg->cgtype[k],cg->cgtype[kl],cg->tolcgtype);
        ntt=SearchIntTable(cg->dih_m,n_ang,cg->told);

        xi[0]=fr->x[ki][0]; xi[1]=fr->x[ki][1]; xi[2]=fr->x[ki][2];
        xj[0]=fr->x[kj][0]; xj[1]=fr->x[kj][1]; xj[2]=fr->x[kj][2];
        xk[0]=fr->x[k][0]; xk[1]=fr->x[k][1]; xk[2]=fr->x[k][2];
        xl[0]=fr->x[kl][0]; xl[1]=fr->x[kl][1]; xl[2]=fr->x[kl][2];

        rvec_sub(xk,xi,r_ki);
        rvec_sub(xj,xi,r_ji);
        rvec_sub(xj,xl,r_jl);

        if(r_ki[0]>fr->box_dim2[0]) r_ki[0]-=fr->box[0][0];
        else if(r_ki[0]<-0.5*fr->box_dim2[0]) r_ki[0]+=fr->box[0][0];
        if(r_ki[1]>fr->box_dim2[1]) r_ki[1]-=fr->box[1][1];
        else if(r_ki[1]<-0.5*fr->box_dim2[1]) r_ki[1]+=fr->box[1][1];
        if(r_ki[2]>fr->box_dim2[2]) r_ki[2]-=fr->box[2][2];
        else if(r_ki[2]<-0.5*fr->box_dim2[2]) r_ki[2]+=fr->box[2][2];

        if(r_ji[0]>fr->box_dim2[0]) r_ji[0]-=fr->box[0][0];
        else if(r_ji[0]<-0.5*fr->box_dim2[0]) r_ji[0]+=fr->box[0][0];
        if(r_ji[1]>fr->box_dim2[1]) r_ji[1]-=fr->box[1][1];
        else if(r_ji[1]<-0.5*fr->box_dim2[1]) r_ji[1]+=fr->box[1][1];
        if(r_ji[2]>fr->box_dim2[2]) r_ji[2]-=fr->box[2][2];
        else if(r_ji[2]<-0.5*fr->box_dim2[2]) r_ji[2]+=fr->box[2][2];

        if(r_jl[0]>fr->box_dim2[0]) r_jl[0]-=fr->box[0][0];
        else if(r_jl[0]<-0.5*fr->box_dim2[0]) r_jl[0]+=fr->box[0][0];
        if(r_jl[1]>fr->box_dim2[1]) r_jl[1]-=fr->box[1][1];
        else if(r_jl[1]<-0.5*fr->box_dim2[1]) r_jl[1]+=fr->box[1][1];
        if(r_jl[2]>fr->box_dim2[2]) r_jl[2]-=fr->box[2][2];
        else if(r_jl[2]<-0.5*fr->box_dim2[2]) r_jl[2]+=fr->box[2][2];




        cprod(r_ki,r_ji,mdih);
        cprod(r_ji,r_jl,ndih);
        cos_phi=cos_angle(mdih,ndih);
        phi=acos(cos_phi);
        iprd=iprod(r_ki,ndih);
        sign=(iprd<0.0)?-1.0:1.0;
        phi=sign*phi;
        //printf("%lf\n",phi+M_PI);


        i_mesh=(int)((M_PI+phi+0.5*cg->space_d)/cg->space_d);
        n_tmp=dis->grid_d*ntt+i_mesh;
        dis->dihedral[n_tmp]+=1.0;
      }
    }

    //read next frame
    read_stat=(*cg->read_next)(fr);
    //ReadNextFrameTrr(fr);  
    printf("\r%d frames have been read.",j+2); fflush(stdout);
  }
  //end of the main loop

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



  double sum;
  double kb=(1.3806503e-23)*cg->temper*(6.0221415e23)/1000.0;
  char name_tmp[100],name_tmp1[100];
  FILE *ff,*ff1;
  double ty;

  for(i=0;i<cg->tolcgtype;i++)
  {
    for(j=i;j<cg->tolcgtype;j++)
    {

      m=TwoBodyNum(i+1,j+1,cg->tolcgtype);
      m=SearchIntTable(cg->bon_m,m,cg->tolb);
      if(m>=0)
      {
        sprintf(name_tmp,"%s_%s_bon.dis",cg->name[i],cg->name[j]);
        sprintf(name_tmp1,"%s_%s_bon.u",cg->name[i],cg->name[j]);
        
        ff=OpenFile(name_tmp,"w");
        ff1=OpenFile(name_tmp1,"w");
        sum=0.0;
        for(k=0;k<dis->grid_b;k++)
        {
          sum+=dis->bond[m*dis->grid_b+k];
        }

        for(k=0;k<dis->grid_b;k++)
        {
          tx=k*cg->space_b;
          if(k==0 || k==dis->grid_b-1) ty=dis->bond[m*dis->grid_b+k]/sum/cg->space_b*2.0;
          else ty=dis->bond[m*dis->grid_b+k]/sum/cg->space_b/(4.0*M_PI*tx*tx);
          fprintf(ff,"%le %le\n",tx,ty*(4.0*M_PI*tx*tx));
          if(ty>1.0e-4)
          fprintf(ff1,"%le %le\n",tx,-1.0*log(ty)*kb);
        }
        fclose(ff);
        fclose(ff1);
      }
    }
  }

  for(i=0;i<cg->tolcgtype;i++)
  {
    for(j=0;j<cg->tolcgtype;j++)
    {
      for(k=j;k<cg->tolcgtype;k++)
      {
        m=ThreeBodyNum(i+1,j+1,k+1,cg->tolcgtype);
        m=SearchIntTable(cg->ang_m,m,cg->tola);

        if(m>=0)
        {
          sprintf(name_tmp,"%s_%s_%s_ang.dis",cg->name[i],cg->name[j],cg->name[k]);
          sprintf(name_tmp1,"%s_%s_%s_ang.u",cg->name[i],cg->name[j],cg->name[k]);
          ff=OpenFile(name_tmp,"w");
          ff1=OpenFile(name_tmp1,"w");
          sum=0.0;
          for(l=0;l<dis->grid_a;l++)
          {
            sum+=dis->angle[m*dis->grid_a+l];
          }

          for(l=0;l<dis->grid_a;l++)
          {
            tx=l*cg->space_a;
            if(l==0 || l==dis->grid_a-1) ty=dis->angle[m*dis->grid_a+l]/sum/cg->space_a*2.0;
            else ty=dis->angle[m*dis->grid_a+l]/sum/cg->space_a/sin(tx);
            if(cg->con_p==1) tx=tx/M_PI*180.0;
            fprintf(ff,"%le %le\n",tx,ty*sin(tx));
            if(ty>1.0e-4)
            fprintf(ff1,"%le %le\n",tx,-1.0*log(ty)*kb);
          }
          fclose(ff);
          fclose(ff1);
        }
      }
    }
  }


  int l2;
  double dih_shift;
  for (i=0;i<cg->tolcgtype;i++) {
    j=i;
    for (k=0;k<cg->tolcgtype;k++)
       for (l=k;l<cg->tolcgtype;l++) {
          m=FourBodyNum(i+1,j+1,k+1,l+1,cg->tolcgtype);
          m=SearchIntTable(cg->dih_m,m,cg->told);

          if (m>=0) {
            sprintf(name_tmp,"%s_%s_%s_%s_dih.dis",cg->name[i],cg->name[j],cg->name[k],cg->name[l]);
            sprintf(name_tmp1,"%s_%s_%s_%s_dih.u",cg->name[i],cg->name[j],cg->name[k],cg->name[l]);
            ff=OpenFile(name_tmp,"w");
            ff1=OpenFile(name_tmp1,"w");
            sum=0.0;
            for (l2=0;l2<dis->grid_d;l2++)
               sum+=dis->dihedral[m*dis->grid_d+l2];

            for (l2=0;l2<dis->grid_d;l2++) {
               tx=l2*cg->space_d;
               if(l2==0 || l2==dis->grid_d-1) ty=dis->dihedral[m*dis->grid_d+l2]/sum/cg->space_d*2.0;
               else ty=dis->dihedral[m*dis->grid_d+l2]/sum/cg->space_d;
               if(cg->con_p==1) { tx=tx/M_PI*180.0; dih_shift=180.0; }
               else dih_shift=M_PI;
               fprintf(ff,"%le %le\n",tx,ty);
               if (ty>1.0e-4)
                 fprintf(ff1,"%le %le\n",tx-dih_shift,-1.0*log(ty)*kb);
            }
            fclose(ff);
            fclose(ff1);
          }
       }

    for (j=i+1;j<cg->tolcgtype;j++)
       for (k=0;k<cg->tolcgtype;k++)
          for (l=0;l<cg->tolcgtype;l++) {
             m=FourBodyNum(i+1,j+1,k+1,l+1,cg->tolcgtype);
             m=SearchIntTable(cg->dih_m,m,cg->told); 
          
             if(m>=0) {
              sprintf(name_tmp,"%s_%s_%s_%s_dih.dis",cg->name[i],cg->name[j],cg->name[k],cg->name[l]);
              sprintf(name_tmp1,"%s_%s_%s_%s_dih.u",cg->name[i],cg->name[j],cg->name[k],cg->name[l]);
              ff=OpenFile(name_tmp,"w");
              ff1=OpenFile(name_tmp1,"w");
              sum=0.0;
              for (l2=0;l2<dis->grid_d;l2++)
                 sum+=dis->dihedral[m*dis->grid_d+l2];

              for (l2=0;l2<dis->grid_d;l2++) {
                 tx=l2*cg->space_d;
                 if(l2==0 || l2==dis->grid_d-1) ty=dis->dihedral[m*dis->grid_d+l2]/sum/cg->space_d*2.0;
                 else ty=dis->dihedral[m*dis->grid_d+l2]/sum/cg->space_d;
                 if(cg->con_p==1) { tx=tx/M_PI*180.0; dih_shift=180.0; }
                 else dih_shift=M_PI;
                 fprintf(ff,"%le %le\n",tx,ty);
                 if (ty>1.0e-4)
                   fprintf(ff1,"%le %le\n",tx-dih_shift,-1.0*log(ty)*kb);
              }  
              fclose(ff);
              fclose(ff1);
            }
          }
 }

  //free memory
  free(cg->bon_m);
  free(cg->ang_m);
  free(cg->dih_m);
  free(dis->bond);
  free(dis->angle);
  free(dis->dihedral);
  for(i=0;i<cg->tolcgtype;i++) free(cg->name[i]);
  free(cg->name);

}

