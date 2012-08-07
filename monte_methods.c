//
//
//  This file is created to experiment with integration methods. 5 August 2012
//
//

#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include "gsl_monte_ftk.h"
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_blas.h>   

double kaverage[5] = { 26.23,-618.97,24.06,	-130.18,-6.27 };	
	
double kernel[5][10] = { 
2.32e-04, -9.71e-01, 2.01e-03, 2.16e+00, -3.70e-03, -1.19e+00, -1.13e-02, 4.67e-02, -5.00e-02, 1.40e-02,
1.57e-02, -1.35e-02, 2.88e-02, 3.69e-02, -8.21e-02, -2.39e-02, 5.25e-01, -1.24e+00, 1.08e+00, -2.74e-01,
2.55e-01, -4.00e-04, -7.32e-01, 2.24e-03, 5.44e-01, -3.81e-03, -2.19e-01, 4.34e-02, 1.88e-01, -9.90e-02,
-1.23e-01, 7.70e-03, 2.59e-01, 2.72e-03, -6.36e-02, -9.17e-03, -5.60e-01, 2.63e-01, 4.87e-01, -3.04e-01,
-2.42e-01, -1.10e-02, 9.81e-02, 1.85e-03, 2.54e-01, 1.10e-02, -1.10e-01, -1.38e-01, -5.20e-02, 1.28e-01};

double myx[10] = {  182.399702,  115.236628,  215.728292,  112.966393,  64.668668,  126.782378,  
                    157.386171,  60.721202,  729.187580,  606.042635};

gsl_matrix_view A ; 
gsl_vector_view B ;

int nconstr = 5;
int ncoords = 10 ; 

void print_k( float * avg, float ** local_kernel );

double chi2f (double *x, size_t dim, void *params){

    double chi2(0);

    for(int i=0;i<nconstr;i++)
    {
        double s = kaverage[i];
        for(int j=0;j<ncoords;j++) {
            s += kernel[i][j]*x[j];
        }
        chi2 += s*s;
    }

    return chi2;
}


double expChi2f (double *x, size_t dim, void *params){


#ifdef STUDY_OUTPUT
    for( unsigned int i = 0 ; i < dim  ; i++  )     
        std::cout << x[i] << " "  ; 

    double mychi2 = chi2f(x, dim, params);


    std::cout << std::setprecision (15) << chi2f(x, dim, params) << " "  ; 
    std::cout << exp(- mychi2/2) << std::setprecision (2) << std::endl  ; 

#endif 
    return exp(-chi2f(x, dim, params)/2);

}

    
void
display_results (const char *title, double result, double error, size_t calls)
{
  printf ("%s ==================\n", title);
  printf ("result = % .6f\n", result);
  printf ("sigma  = % .6f\n", error);
  printf ("error  = % .6f \n", error/result );
  printf ("calls_used  = %d \n", calls );

}
    

int main (void)
{
	double res, err;
  
	#define DIM 10   
	
	double xl[DIM] = { 168,  108,192, 108, 48, 108, 140, 60,720, 600};
	double xu[DIM] = { 192, 144, 216, 144, 72, 144, 160, 80, 740, 620}; /**/
	
	
	gsl_vector_fprintf ( stdout , & B.vector, "%f");

	printf ( "%s %f\n ", "mychi2",  chi2f( myx, 10, 0) );

	int loop=0;
	for (loop=0; loop<100; loop++) {

		const gsl_rng_type *T;
		gsl_rng *r;

		gsl_monte_function G = { &expChi2f, DIM, 0   };

		size_t calls = 500000;
		size_t max_calls = 1500000;
		double min_err = 0.30;
        size_t required_calls = 0;

		gsl_rng_env_setup ();

		T = gsl_rng_default;
		r = gsl_rng_alloc (T);

		{
			gsl_monte_ftk_state *s = gsl_monte_ftk_alloc (DIM);
			gsl_monte_ftk_integrate_2 (&G, xl, xu, DIM, calls, max_calls, min_err , r, s, &required_calls	,&res, &err);
			gsl_monte_ftk_free (s);

            if(loop == 0)			
			    display_results ("ftk", res, err, required_calls);
		}

		/*{
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
		}*/
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
