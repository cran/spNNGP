#ifndef R_NO_REMAP
#  define R_NO_REMAP
#endif

double exprnd(double mu);
double aterm(int n, double x, double t);
double randinvg(double mu);
double truncgamma();
double tinvgauss(double z, double t);
double samplepg(double z);
double rpg(int n, double z);
