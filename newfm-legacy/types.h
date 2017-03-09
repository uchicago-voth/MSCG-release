#define VERYSMALL    1.0e-14 
#define VERYSMALL_F  1.0e-6 //small number for single precision
//#define MAX_BON 4    //max bonded interactions for each cg site
//#define MAX_ANG 12   //max angle
//#define MAX_DIH 36  //max dihedral
#define FORCE_CAP 1000.0 //filter some noisy data
#define MAX_TYPENAME_LEN 4 //max length for CG type names
#define TOL_TWO(n) (n*(n+1)/2) //total possible interaction types from CG types
#define TOL_THREE(n) (n*n*(n+1)/2) //three body
#define TOL_FOUR(n) (n*n*(n*n+1)/2) //four body
const double R2D=180.0/M_PI;

typedef float real;
typedef struct
{ 
  //parameters in control.in
  int n_block; //block size 
  int pres_con; //1 to use pressure constrain
  int start_fr; //frame number to start
  int fr_n; //total frame number
  double fr_n_1; //1.0/fr_n
  double cutoff; //nonebonded cutoff
  double space; //nonebonded bin width for FM
  double space_b; //bond
  double space_a; //angle
  double space_d; //dihedral
  double spaceo; //nonebonded bin width for output
  double spaceo_b;
  double spaceo_a;
  double spaceo_d;
  int n_k; //b-spline k-value for nonebonded
  int nb_k; 
  int na_k; 
  int nd_k; 
  int ba_type; //basis set type; 0 for b-spline; 1 for linear
  int be_sparse; //1 to use sparse format
  int con_p; // output style; 0 for tables; 2 for tables+binary; 3 for binary 
  int max_bon; //max bonds connected to one CG site
  int max_ang;
  int max_dih;
  int x_out; //0: do nothing; 1: output x.out that stores the solution
  int iter;
  double lambda;
  int input_lambda; //0: do nothing; 1: calculate results from lambda; 2: calculate results from file lambda.in
  int ang_type; //0: Sergey; 1: real
  int dih_type;
  double atol; //LSQR parameter; accuracy of matrix A
  double btol; //LSQR parameter; accuracy of vector B
  double conlim; //LSQR parameter; condition number limit
  int itnlim; //LSQR parameter; max number of iterations
  double rcond; //SVD condtional number thredhold
  int use_3b; //0: no three body; 1: angular only; 2: SW with spline angular; 3: SW with fixed theta0. option 2 only for b-spline
  double space_3; //theree body bin width
  double spaceo_3;
  int n3_k; //there body b-spline k
  double gamma; //SW parameter
  int tb_exc;

  double *cutoff_3i;  
  double *cos_theta0i;//SW parameter
  double cos_theta0;

  double damp; //LSQR parameter damping factor
  int mm1; //LSQR parameter; mm-1
  int istop; //LSQR parameter 
  int itn; //LSQR parameter
  double anorm; //LSQR parameter
  double acond; //LSQR parameter
  double rnorm; //LSQR parameter
  double arnorm; //LSQR parameter
  double xnorm; //LSQR parameter
  double cutoff2; //cutoff^2
  int tolcgn; //total number of CG sites
  int tolcgtype; //total number of CG types
  int *cgtype; //CG type for each site
  int *bontype; //if one bondtype exists; 1 for existing
  int *blist_n; //number of bonds connected to one site
  int **blist_i; //bond table for one site
  int *angtype; 
  int *alist_n; 
  int **alist_i; 
  int *dihtype; 
  int *dlist_n; 
  int **dlist_i; 
  char **name; //type names
  int tollj; //toltal number of nonebonded interactions
  int tolb;
  int tola;
  int told;
  int tollj1; //total number of nonebonded interactions from FM
  int tolb1;
  int tola1;
  int told1;
  int tollj2; //from table.in
  int tolb2;
  int tola2;
  int told2;
  int *bon_m; //hash number for each bond type
  int *ang_m;
  int *dih_m;
  int *ljnum; //type index for each nonebonded; 0 for not existing; posive N for Nth interaction from FM; negtive N for Nth interaction from table.in
  int *bonnum;
  int *angnum;
  int *dihnum;
  int *grid_n; //accumulated numbers of unknowns for each nonebonded interaction
  int *grid_bn;
  int *grid_an;
  int *grid_dn;
  double *rmin; //min distance for each interaction
  double *rmax;
  double *rmin_b;
  double *rmax_b;
  double *rmin_a;
  double *rmax_a;
  double *rmin_d;
  double *rmax_d;
  double **es_co; //external force table linear spline coefficients for nonebonded
  double **esb_co;
  double **esa_co;
  double **esd_co;
  double space_es;//bin width for external table splines
  double spaceb_es;
  double spacea_es;
  double spaced_es;
  gsl_bspline_workspace **w_bs; //b-spline work spaces for nonebonded 
  gsl_bspline_workspace **wb_bs;
  gsl_bspline_workspace **wa_bs;
  gsl_bspline_workspace **wd_bs;
  gsl_vector *b_bs; //b-spline vectors for nonebonded
  gsl_vector *bb_bs;
  gsl_vector *ba_bs;
  gsl_vector *bd_bs;
  FILE *diff;
  FILE *diff_vir;
  int (*read_next)(); //read next frame 
  void (*set_zero)(); //initialize matrix for each block
  void (*mat_insert)(); //insert an element into matrix A
  void (*vir_insert)(); //insert an element for pressure constrain
  void (*mat_op)(); //matrix computation
  void (*fr_op)();
  void (*fm_basis)();
  void (*fm_nb)();
  void (*fm_bon)();
  void (*fm_ang)();
  void (*fm_dih)();
  double cutoff2_3;
  int *grid_3n;
  int *tb_m;
  gsl_bspline_workspace **w3_bs;
  gsl_bspline_deriv_workspace *w3_bs1;
  gsl_vector *b3_bs; 
  gsl_matrix *b3_bs1;
  double *rmin_3;
  double *rmax_3;
  int *tbnum;
  int *tb_n;
  int **tb_list;
  int tol3;
  int output_residual;
  int out_sp; //output splines coefficients
  int out_dd;
  double iter_ld;
  void (*fill_rhs)();
  void (*fill_virial)();
  void (*table_twobody)();
  void (*table_threebody)();
  void (*table_fourbody)();
  void (*fm_threebody)();

  double temper; //temperature
} CG_PAR; //CG parameters

typedef struct 
{ 
  int trj_type; //trajectory type
  int natoms; //total number of sites
  char trj_name[1000]; //file name 
  char trj_name1[1000]; //file name #2 (for xtc)
  XDRFILE *fp; //file number
  XDRFILE *fp1; //file number #2
  matrix box; //box dimension
  real box_dim2[3]; //half box length
  rvec *x; //positions
  rvec *f; //forces
  int ncell; //number of cells
  int *list; //linked-lists
  int *head; //head of linked-lists
  int m0; //number of cells on x-direction
  int m1;
  int m2;
  int *map; //neighbor cell map
  real cellx; //cell dimension on x-direction
  real celly;
  real cellz;
  double *b_p; //pressure constain RHS 
  int step;
  real time;
  rvec *df;
  int read_fr;
  int ncell_3; //number of cells
  int *list_3; //linked-lists
  int *head_3; //head of linked-lists
  int m0_3; //number of cells on x-direction
  int m1_3;
  int m2_3;
  int *map_3; //neighbor cell map
  real cellx_3; //cell dimension on x-direction
  real celly_3;
  real cellz_3;

} FR_DAT; //trajectory frame data
 
struct MAT_SP //sparse matrix structure, x,y,z components are stored together
{
  int col; //column number
  double valx[3]; //x,y,z compoments
  struct MAT_SP *next; //pointer to the next element in the linked list
};

struct MAT_H //head struct for mat_sp, one linked list for one row
{
  int n; //total number of non-zero elements in this row
  struct MAT_SP *h;
};

typedef struct
{ 
  //FM matrix equation is Ax=b and the normal equation is Gy=d
  int mm; //number of rows for FM matrix
  int nn; //number of columns
  int mv; //rows without pressure constrains
  int mp; //rows for pressure constrains
  double *aa; //dense matrix A
  double *bb; //vector B
  double *gg; //G
  double *dd; //d
  //for sparse matrix format
  int *iw; //row sizes
  int *jw;//column indexes
  double *rw; //element values
  int *n_xx; //number of sampled times for one unknown
  double *xx; //final answers
  double *xx_b; //answers from one block
  double *v; //temp space for LSQR
  double *w; //temp space for LSQR
  double *se; //temp space for LSQR
  double *h; //for preconditioning
  int n_nz; //toltal number of none zero values
  struct MAT_H *mat_head; //head structure for sparse format matrix A
  double residual;  
  int ma;
  int na;
  int row_shift;
  int outer;
  int lwork;
  double *work;
  double *tau;
} MAT_DAT; //matrix data

enum INTER_TYPE {non,bon,ang,dih,t_b};

typedef struct
{
  int inner; //inner loop number
  int m_shift; //shifted row number
  //bonds/nonebonded: k-l; angles: k-j-l; dihedrals: k-i-j-l;
  int k;
  int l;
  int i;
  int j;
  int m; //hash number for this interaction type from bon_m/ang_m/dih_m
  int mi; //posion in the array bon_m/ang_m/dih_m; for nonebonded m=mi
  int numi; //interaction type index from ljnum/bonnum/angnum/dihnum
  double rr; //pair distance
  double rr0; //pair distance from the rmin
  int grid_base; //when count number of unknowns, start from this number
  int n_k; //b-spline k
  int *grid; //grid array
  double *rmin; 
  double *rmax;
  int *num; //type index array
  double space; //bin width
  gsl_vector *bs; //b-spline vector
  gsl_bspline_workspace **w; //b-spline work space
  double **es_co; //external table coefficients
  double space_es; //bin width for external
  int (*cal_m)(); //function to calculate hash number
  int *m_array; //hash number array
  int m_size; //size of m_array
  double coef[100];
  int n_coef;
  double cutoff2;
  void (*fm_nbody)();
  int i_mesh;
} TYPE_INFO; //info needed for FM calculation of each pair


//calculate two-body hash number
int TwoBodyNum(int i, int j, const int tolcgtype)
{
    int tn;
    if(i>j){tn=i;i=j;j=tn;}
    return ((i-1)*tolcgtype+j-1)-i*(i-1)/2; 
}

//three-body hash number
int ThreeBodyNum(int i, int j, int k, const int tolcgtype)
{
  int tn;
  if(j>k) {tn=j;j=k;k=tn;}
  return (i-1)*tolcgtype*(tolcgtype+1)/2+(2*tolcgtype-j+2)*(j-1)/2+k-j;
}

//four-body hash number
int FourBodyNum(int a, int b, int c, int d, const int n)
{
  int tmp,n_ab,n1_ab,n1,n2;
  if(a>b)
  {
    {tmp=a; a=b; b=tmp; tmp=c; c=d; d=tmp;}
    n_ab=((a-1)*n+b-1)-a*(a-1)/2;
    n1_ab=n_ab-a;
    n1=n1_ab*n*n+a*(n*(n+1)/2);
    n2=(c-1)*n+d-1;
  }
  else if(a<b)
  {
    n_ab=((a-1)*n+b-1)-a*(a-1)/2;
    n1_ab=n_ab-a;
    n1=n1_ab*n*n+a*(n*(n+1)/2);
    n2=(c-1)*n+d-1;
  }
  else
  {
    n_ab=((a-1)*n+b-1)-a*(a-1)/2;
    n1_ab=n_ab-(a-1);
    n1=n1_ab*n*n+(a-1)*(n*(n+1)/2);
    if(c>d) {tmp=c; c=d; d=tmp;}
    n2=((c-1)*n+d-1)-c*(c-1)/2;
  }
  return (n1+n2);
}
 
//search in an integer table
int SearchIntTable(const int *a, const int m, const int size)
{
  int left, right, mid;
  left=0; right=size-1;
  int ret=-1;
  while(left<=right)
  {
    mid=(left+right)/2;
    if(m==a[mid]) {ret=mid; break;}
    else if(m<a[mid]) right=mid-1;
    else left=mid+1;
  }
  //if(ret==-1) {printf("Sth. wrong with the bonded interaction type searching!\n"); exit(0);
  return ret;
}

void * MemAlloc(const int size)
{
  void *p;
  p=malloc((unsigned)size);
  if(p==NULL) {printf("Not enough memory!\n"); exit(1);}
  return p;
}

void * MemRealloc(void * p, const int size)
{
  void *p1;
  p1=realloc(p,(unsigned)size);
  if(p1==NULL) {printf("Not enough memory!\n"); exit(1);}
  return p1;
}

FILE * OpenFile(char *file, char *mode)
{
  FILE *fp;
  fp=fopen(file,mode);
  if(fp==NULL) {printf("Can not open file: %s!\n",file); exit(1);}
  return fp;
}


