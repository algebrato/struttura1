#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <lapacke.h>
#define ABS(a)     (((a) < 0) ? -(a) : (a))

int main()
{
  /*
   *	hbar^2/2m = 1
  */
  
  static double pi = M_PI;
  int n, npw;
  double v0, a, b, ecut, k;
  double *g, *e, *h, *wrk;
  int i, j, ij, m, lwork, info, p;
  char *V = "V";
  char *U = "U";
  FILE *out;
  
  /*  Input data
      Potential: V(x)=-V_0 for |x|<b/2, V(x)=0 for |x|>b/2
      Periodicity:    V(x+a)=V(x)
  */
  fprintf(stdout, "Parameters for potential well: V_0, a, b > ");
  scanf("%lf %lf %lf",&v0, &a, &b);
  if ( v0 <= 0 ||  a <= 0 || b <= 0 || a <= b )  { 
    fprintf(stderr, "wrong input parameters\n");
    exit(1);
  }
  fprintf (stdout,"   V_0=%f, a0=%f, b=%f\n",v0,a,b);
  /*
    Plane waves basis set: G_n=n*2pi/a, \hbar^2/2m*G^2 < Ecut
  */
  fprintf(stdout, "Cutoff for plane waves: ecut > ");
  scanf("%lf",&ecut);
  if ( ecut <= 0) {
    fprintf(stderr, "wrong input parameter\n");
    exit(2);
  }
  /*
    Number of plane waves
  */
  npw = (int) ( sqrt ( ecut/pow( 2.0*pi/a, 2) ) + 0.5 );
  npw = 2*npw+1;
  fprintf (stdout,"   ecut = %f,  # of PWs=%d\n",ecut, npw);
  /*
    Assign values of  G_n: n=0,+1,-1,+2,-2, etc
  */
  n=100;
  g  = (double *) malloc ( npw * sizeof (double) );
  e  = (double *) malloc ( npw * sizeof (double) );
  h  = (double *) malloc ( n*n * sizeof (double) );
  wrk= (double *) malloc (3*npw* sizeof (double) );

  g[0] = 0.0;
  for ( i = 1; i < npw; i+=2 ) {
    g[i  ] = (i+1)*pi/a;
    g[i+1] =-(i+1)*pi/a;
  }
  for ( ij = 0; ij < n*n; ++ij ) {
    h[ij] = 0.0;
  }
  /*
    Loop on k-vectors: k runs from -pi/a to pi/a
  */
  out = fopen("bands.out", "w");
  for ( m =-n; m <= n; m++ ) {
    k = m*pi/n/a;
    /*
      Assign values of the matrix elements of the hamiltonian 
      on the plane wave basis
    */
    ij = 0;
    for ( i = 0; i < npw; ++i ) {
      for ( j = 0; j < npw; ++j ) {
      /* NOTA BENE: the matrix h is a vector in the calling program,
         while dsyev expects a (pointer to a) fortran matrix.
         A fortran matrix is "simulated" in the following way:
         if h[ij] is the C array and  h(i,j) is a N*M Fortran matrix,
         h[ij] == h(i,j), where  ij = (j-1)*N + (i-1)  (i=1,N, j=1,M)
      */

	if ( i == j ) {
	  h[ij++] = (k+g[i])*(k+g[i]) - v0/a*b; 
	} else {
	  h[ij++] = -v0/a * sin( (g[j]-g[i])*b/2.0 ) / (g[j]-g[i])*2.0;
	}
      }
    }
    /*
      Solution [expansion coefficients are stored into h(i,j)
      j=basis function index, i= eigenvalue index]
      (beware fortran-C reversed index convention!)
    */
    lwork = 3*npw;
    dsyev_ ( V, U, &npw, h, &npw, e, wrk, &lwork, &info ); 

    if ( info != 0) {
      fprintf(stderr, "H-matrix diagonalization failed\n");
      exit(3);
    }
    
    /*printf("%f ",k);
    for(p=0; p<npw; p++) printf("%f ",e[p]);printf("\n");*/
   
    fprintf(out,"%12.6f ",k);
    for(p=0; p<npw; p++) fprintf(out,"%10.6f ",e[p]);
    fprintf(out,"\n");



  }
  fclose(out);

  return 0;

}
