void SetTopPar(const char *par, const char *val, CG_PAR *const cg, const int line)
{
    
  if(strcmp("block_size",par)==0) sscanf(val,"%d",&cg->n_block);
  else if(strcmp("pressure_constrain",par)==0) sscanf(val,"%d",&cg->pres_con);
  else if(strcmp("start_frame",par)==0) sscanf(val,"%d",&cg->start_fr);  
  else if(strcmp("total_frame",par)==0) sscanf(val,"%d",&cg->fr_n);
  else if(strcmp("cut_off", par)==0) sscanf(val,"%lf",&cg->cutoff); 
  else if(strcmp("nonebonded_bin",par)==0) sscanf(val,"%lf",&cg->space);
  else if(strcmp("bond_bin",par)==0) sscanf(val,"%lf",&cg->space_b);
  else if(strcmp("angle_bin",par)==0) sscanf(val,"%lf",&cg->space_a);
  else if(strcmp("dihedral_bin",par)==0) sscanf(val,"%lf",&cg->space_d);
  else if(strcmp("k_nonebonded",par)==0) sscanf(val,"%d",&cg->n_k);
  else if(strcmp("k_bond",par)==0) sscanf(val,"%d",&cg->nb_k);
  else if(strcmp("k_angle",par)==0) sscanf(val,"%d",&cg->na_k);
  else if(strcmp("k_dihedral",par)==0) sscanf(val,"%d",&cg->nd_k);
  else if(strcmp("basis_type",par)==0) sscanf(val,"%d", &cg->ba_type);
  else if(strcmp("use_sparse", par)==0) sscanf(val,"%d", &cg->be_sparse);
  else if(strcmp("nonebonded_outbin", par)==0) sscanf(val,"%lf", &cg->spaceo);
  else if(strcmp("bond_outbin", par)==0) sscanf(val,"%lf", &cg->spaceo_b);
  else if(strcmp("angle_outbin", par)==0) sscanf(val,"%lf", &cg->spaceo_a);
  else if(strcmp("dihedral_outbin", par)==0) sscanf(val,"%lf", &cg->spaceo_d);
  else if(strcmp("output_style",par)==0) sscanf(val,"%d", &cg->con_p);
  else if(strcmp("atol", par)==0) sscanf(val,"%lf", &cg->atol);
  else if(strcmp("btol", par)==0) sscanf(val,"%lf", &cg->btol);
  else if(strcmp("conlim", par)==0) sscanf(val,"%lf", &cg->conlim);
  else if(strcmp("itnlim",par)==0) sscanf(val,"%d", &cg->itnlim);
  else if(strcmp("rcond",par)==0) sscanf(val,"%lf", &cg->rcond);
  else if(strcmp("max_bond",par)==0) sscanf(val,"%d", &cg->max_bon);
  else if(strcmp("max_angle",par)==0) sscanf(val,"%d", &cg->max_ang);
  else if(strcmp("max_dihedral",par)==0) sscanf(val,"%d", &cg->max_dih);
  else if(strcmp("x_out",par)==0) sscanf(val,"%d", &cg->x_out);
  else if(strcmp("iterative",par)==0) sscanf(val,"%d", &cg->iter);
  else if(strcmp("lambda",par)==0) sscanf(val,"%lf", &cg->lambda);
  else if(strcmp("input_lambda",par)==0) sscanf(val,"%d", &cg->input_lambda);
  else if(strcmp("angle_type",par)==0) sscanf(val,"%d", &cg->ang_type);
  else if(strcmp("dihedral_type",par)==0) sscanf(val,"%d", &cg->dih_type);
  else if(strcmp("threebody",par)==0) sscanf(val,"%d", &cg->use_3b);
  else if(strcmp("threebody_bin",par)==0) sscanf(val,"%lf", &cg->space_3);
  else if(strcmp("threebody_outbin",par)==0) sscanf(val,"%lf", &cg->spaceo_3);
  else if(strcmp("k_threebody",par)==0) sscanf(val,"%d", &cg->n3_k);
  else if(strcmp("output_residual",par)==0) sscanf(val,"%d", &cg->output_residual);
  else if(strcmp("threebody_gamma",par)==0) sscanf(val,"%lf", &cg->gamma);
  else if(strcmp("threebody_exclude",par)==0) sscanf(val,"%d", &cg->tb_exc);
  else if(strcmp("output_spline",par)==0) sscanf(val,"%d", &cg->out_sp);
  else if(strcmp("output_d",par)==0) sscanf(val,"%d", &cg->out_dd);
  //else if(strcmp("threebody_cutoff",par)==0) sscanf(val,"%lf", &cg->cutoff_3);
  else if(strcmp("iterative_lambda",par)==0) sscanf(val,"%lf", &cg->iter_ld);
  else if(strcmp("temperature",par)==0) sscanf(val,"%lf", &cg->temper);
  else printf("Warning: Wrong parameter name in control.in: line %d!\n",line);  
  
}

void ReadControl (CG_PAR *const cg)
{
  char buff[100];

  cg->n_block=10;
  cg->pres_con=0;
  cg->start_fr=1;
  cg->fr_n=10;
  cg->cutoff=1.0;
  cg->space=0.05;
  cg->space_b=0.05;
  cg->space_a=0.05;
  cg->space_d=0.05;
  cg->n_k=4;
  cg->nb_k=4;
  cg->na_k=4;
  cg->nd_k=4;
  cg->ba_type=0;
  cg->be_sparse=1;
  cg->spaceo=0.002;
  cg->spaceo_b=0.002;
  cg->spaceo_a=0.002;
  cg->spaceo_d=0.002;
  cg->con_p=0;
  cg->atol=1.0e-6;
  cg->btol=1.0e-6;
  cg->conlim=1.0e14;
  cg->itnlim=0;
  cg->rcond=-1.0;
  cg->max_bon=4;
  cg->max_ang=12;
  cg->max_dih=36;
  cg->x_out=0;
  cg->iter=0;
  cg->iter=0;
  cg->lambda=0.0;
  cg->input_lambda=0;
  cg->ang_type=0;
  cg->dih_type=0;
  cg->use_3b=0;
  cg->space_3=1.0;
  cg->spaceo_3=0.2;
  cg->n3_k=4;
  cg->output_residual=0;
  cg->gamma=0.12;
  //cg->cutoff_3=0.37;
  cg->tb_exc=0;
  cg->out_sp=0;
  cg->out_dd=0;
  cg->iter_ld=1.0;
  cg->temper=310;

  FILE *control_in;
  int i;
  char *left,*right;
  control_in=OpenFile("control.in","r");

  int line=1;
  while(fgets(buff,100,control_in)!=NULL)
  {
    char left[50],right[50];
    sscanf(buff,"%s%s",left,right);
    SetTopPar(left,right,cg,line);
    line++;
  }
  fclose(control_in);
  cg->cutoff2=cg->cutoff*cg->cutoff;
  if(cg->be_sparse!=0) cg->out_dd=0;
  cg->fr_n_1=1.0/cg->fr_n;
}

