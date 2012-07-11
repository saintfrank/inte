#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
   
	double a[5][10] = { 0.00,-0.97,	0.00,	2.16,	-0.00,	-1.19,	-0.01,	0.05,	-0.05,	0.01,	
	0.02,	-0.01,	0.03,	0.04,	-0.08,	-0.02,	0.53,	-1.24,	1.08,	-0.27	,
	0.25,	-0.00,	-0.73,	0.00,	0.54,	-0.00,	-0.22,	0.04,	0.19,	-0.10	,
	-0.12,	0.01,	0.26,	0.00,	-0.06,	-0.01,	-0.56,	0.26,	0.49,	-0.30	,
	-0.24,	-0.01,	0.10,	0.00,	0.25,	0.01,	-0.11,	-0.14,	-0.05,	0.13};


double
myG (double *k, size_t dim, void *params)
{
  int i, j;

  for (int i = 0 ; i < 100 ; i ++)
  {	  
	double scalar_product=0;
  	for (i=0; i<dim; i++)
		for (j=0; j<dim; j++) 
	  	{
		  scalar_product += k[i]*k[j]*par[i][j]*0.01; // use (i*dim*0.3+j) as matrix element
	  	}	
  }

  return exp(-scalar_product);
}
     
int main (void)
{
  double res, err;

  #define DIM 10   
  
  double xl[DIM] = { 1, -2, 3, -2, 4, -6, 2, -1, 2, -3 }; /**/
  double xu[DIM] = { 1.2, 2.2, 3.5, 2.2, 4.5, 6.7, 2.3, 1.2, 2.6, 4}; /**/

  int i, j;

  myG (xl, 0, 0);


  printf("scalar prod. = %f", (float) scalar_product); */
  
  float myexp = exp(-scalar_product);

  return ;
}
