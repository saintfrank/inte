#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_blas.h>   

double average[5] = { 26.23,-618.97,24.06,	-130.18,-6.27 };	

double kernel[5][11] = { 0.00,-0.97,	0.00,	2.16,	-0.00,	-1.19,	-0.01,	0.05,	-0.05,	0.01, 26.23,	
	0.02,	-0.01,	0.03,	0.04,	-0.08,	-0.02,	0.53,	-1.24,	1.08,	-0.27, -618.97	, 
	0.25,	-0.00,	-0.73,	0.00,	0.54,	-0.00,	-0.22,	0.04,	0.19,	-0.10, 24.06	,
	-0.12,	0.01,	0.26,	0.00,	-0.06,	-0.01,	-0.56,	0.26,	0.49,	-0.30, -130.18	,
	-0.24,	-0.01,	0.10,	0.00,	0.25,	0.01,	-0.11,	-0.14,	-0.05,	0.13, -6.27 };

double kernel_vec[50] = { 0.00,-0.97,	0.00,	2.16,	-0.00,	-1.19,	-0.01,	0.05,	-0.05,	0.01,	
	0.02,	-0.01,	0.03,	0.04,	-0.08,	-0.02,	0.53,	-1.24,	1.08,	-0.27	, 
	0.25,	-0.00,	-0.73,	0.00,	0.54,	-0.00,	-0.22,	0.04,	0.19,	-0.10	,
	-0.12,	0.01,	0.26,	0.00,	-0.06,	-0.01,	-0.56,	0.26,	0.49,	-0.30	,
	-0.24,	-0.01,	0.10,	0.00,	0.25,	0.01,	-0.11,	-0.14,	-0.05,	0.13 };

int nconstr = 5;
int ncoords = 10 ; 

double chi2f (double *x, size_t dim, void *params){

	double chi2 = 0;

	double ** mykernel = (double **) params;


	//gsl_matrix_view A = gsl_matrix_view_array(kernel_vec, 5, 10); 
	//gsl_vector_view X = gsl_vector_view_array(x, 10); 
	//gsl_vector_view B = gsl_vector_view_array(average, 5); 


	//gsl_blas_dgemv ( CblasNoTrans, 1, & A.matrix, & X.vector , 1, & B.vector );


	for( int i = 0 ;i<nconstr;++i)
	{
		double s = mykernel[i][10];
		for(int j=0;j<ncoords;++j) { 
			s += mykernel[i][j]*x[j];
		} 
		chi2 += s*s;
	}   

	return chi2;
}
/*
double chi2f (double *x, size_t dim, void *params){

	double chi2 = 0;

	float ** mykernel = (float **) params;

	for( int i = 0 ;i<nconstr;++i)
	{
		double s = mykernel[i][10];
		for(int j=0;j<ncoords;++j) { 
			s += mykernel[i][j]*x[j];
		} 
		chi2 += s*s;
	}   

	return chi2;
}*/

double expChi2f (double *x, size_t dim, void *params){

	//return chi2f(x, dim, params);
	return exp(-chi2f(x, dim, params)/2);
	//return 1 ; 
}

/*
double
g (double *k, size_t dim, void *params)
{
  int i, j;

  double scalar_product=0;
  for (i=0; i<dim; i++)
    for (j=0; j<dim; j++) {
      scalar_product += k[i]*k[j]*(i*dim*1.5+j+1.2)*0.01; // use (i*dim*0.3+j) as matrix element
    }
  //  printf("scalar prod. = %f", (float) scalar_product); 
  return exp(-scalar_product);
}
*/     







     
void
display_results (char *title, double result, double error)
{
  printf ("%s ==================\n", title);
  printf ("result = % .6f\n", result);
  printf ("sigma  = % .6f\n", error);
  printf ("error  = % .6f \n", error/result );
}
    

int main (void)
{
	double res, err;
  
	#define DIM 10   
	//double xl[DIM] = { 1, -2, 3, -2, 4, -6, 2, -1, 2, -3 }; /**/
        //double xu[DIM] = { 1.2, 2.2, 3.5, 2.2, 4.5, 6.7, 2.3, 1.2, 2.6, 4}; /**/

	double xl[DIM] = { 168, 108, 192, 108, 48, 108, 140, 60, 720, 600 }; /**/
	double xu[DIM] = { 192, 144, 216, 144, 72, 144, 160, 80, 740, 620}; /**/
		
	register float ** local_kernel   = (float **)calloc (5, sizeof(float*));
	for(int i = 0 ; i < 5 ; i++)
	{
		local_kernel[i] = (float *)calloc (11, sizeof(float));
		for (int j=0; j < 10 ; j++)
			local_kernel[i][j] = kernel[i][j];



		local_kernel[i][10] = average[i];
	}


	int loop=0;
	for (loop=0; loop<100; loop++) {

		const gsl_rng_type *T;
		gsl_rng *r;

		gsl_monte_function G = { &expChi2f, DIM, local_kernel  };

		size_t calls = 500000;

		gsl_rng_env_setup ();

		T = gsl_rng_default;
		r = gsl_rng_alloc (T);

		{
			gsl_monte_plain_state *s = gsl_monte_plain_alloc (DIM);
			gsl_monte_plain_integrate (&G, xl, xu, DIM, calls, r, s, 
					&res, &err);
			gsl_monte_plain_free (s);

			
			display_results ("plain", res, err);
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
