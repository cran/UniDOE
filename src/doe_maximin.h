void create_maximin(double **x,int nnew1, int np1,int nv1,CRITOPT *critopt);
void free_maximin();
double maximin_set(double **x0);
double maximin_cp_set(int ncol,int ncp,int *idx1,int *idx2);
double maximin_cp(int ncol,int ncp,int *idx1,int *idx2);
double maximin_cp1(int ncol,int i1,int i2);
double maximin_pm(int ncol, int npm,int *idx1,int *idx2);
double maximin_pm_set(int ncol, int npm,int *idx1,int *idx2);
double maximin();
double maximin_eval(double **x0);
double **maximin_x(double **xnw);
void maximin_snap(int ncol);
void maximin_reset(int ncol);
void maximin_full_snap();
void maximin_full_reset();
void maximin_global_x(double **xnew);