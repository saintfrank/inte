#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
     
/* Computation of the integral,
     
   I = int (dx dy dz)/(2pi)^3  1/(1-cos(x)cos(y)cos(z))
     
   over (-pi,-pi,-pi) to (+pi, +pi, +pi).  The exact answer
   is Gamma(1/4)^4/(4 pi^3).  This example is taken from
   C.Itzykson, J.M.Drouffe, "Statistical Field Theory -
   Volume 1", Section 1.1, p21, which cites the original
   paper M.L.Glasser, I.J.Zucker, Proc.Natl.Acad.Sci.USA 74
   1800 (1977) */
     
/* For simplicity we compute the integral over the region 
   (0,0,0) -> (pi,pi,pi) and multiply by 8 */
     
/*double exact = 1.3932039296856768591842462603255;*/
double exact = 41.4;
     
double
g (double *k, size_t dim, void *params)
{
  int i, j;

  double scalar_product=0;
  for (i=0; i<dim; i++)
    for (j=0; j<dim; j++) {
      scalar_product += k[i]*k[j]*(i*dim*1.5+j+1.2)*0.01; // use (i*dim*0.3+j) as matrix element
    }
  /*  printf("scalar prod. = %f", (float) scalar_product); */
  return exp(-scalar_product);
}
     
void
display_results (char *title, double result, double error)
{
  printf ("%s ==================\n", title);
  printf ("result = % .6f\n", result);
  printf ("sigma  = % .6f\n", error);
  printf ("exact  = % .6f\n", exact);
  printf ("error  = % .6f = %.2g sigma\n", result - exact,
	  fabs (result - exact) / error);
}
     
int
main (void)
{
  double res, err;

#define DIM 10   
  double xl[DIM] = { 1, -2, 3, -2, 4, -6, 2, -1, 2, -3 }; /**/
  double xu[DIM] = { 1.2, 2.2, 3.5, 2.2, 4.5, 6.7, 2.3, 1.2, 2.6, 4}; /**/

  int loop=0;
  for (loop=0; loop<100; loop++) {
    const gsl_rng_type *T;
    gsl_rng *r;
     
    gsl_monte_function G = { &g, DIM, 0 };
     
    size_t calls = 50000;
     
    gsl_rng_env_setup ();
     
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
     
    /*    {
      gsl_monte_plain_state *s = gsl_monte_plain_alloc (DIM);
      gsl_monte_plain_integrate (&G, xl, xu, DIM, calls, r, s, 
				 &res, &err);
      gsl_monte_plain_free (s);
     
      display_results ("plain", res, err);
    }
    */
    {
      gsl_monte_miser_state *s = gsl_monte_miser_alloc (DIM);
      gsl_monte_miser_integrate (&G, xl, xu, DIM, calls, r, s,
				 &res, &err);     
      if (loop<2) {
	display_results ("miser", res, err); 
	
	printf ("converging...\n");
	int check_RMS=0, max=10;
	double tot=0, tot2=0;
	for (check_RMS=0; check_RMS<max; check_RMS++) {
	  gsl_monte_miser_integrate (&G, xl, xu, DIM, calls, r, s,
				     &res, &err);
	  tot+=res;
	  tot2+=res*res;
	  printf ("result = % .6f sigma = % .6f rel. err. = % .6f\n", res, err, err/res);
	}
	double n_inverse = 1./max;
	double mean = tot*n_inverse;
	double rms  = sqrt(fabs(tot2*n_inverse -mean*mean));
	printf ("mean = % .6f rms = % .6f rel. err. = % .6f\n", mean, rms, rms/mean);
      }
      gsl_monte_miser_free (s);
    }
    /*     
	   {
	   gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (DIM);
     
	   gsl_monte_vegas_integrate (&G, xl, xu, DIM, 10000, r, s,
	   &res, &err);
	   display_results ("vegas warm-up", res, err);
     
	   printf ("converging...\n");
     
	   do
	   {
	   gsl_monte_vegas_integrate (&G, xl, xu, DIM, calls/5, r, s,
	   &res, &err);
	   printf ("result = % .6f sigma = % .6f "
	   "chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
	   }
	   while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
     
	   display_results ("vegas final", res, err);
     
	   gsl_monte_vegas_free (s);
	   }
    */ 
    gsl_rng_free (r);
  } /* end loop */
     
  return 0;
}
