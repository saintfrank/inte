/* monte/plain.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2009 Michael Booth
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Plain Monte-Carlo. */

/* Author: MJB */

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include "gsl_monte_ftk.h"
#include <string.h> 

//#define FOLLOW_ALGO
//#define DEBUG_FOLLOW_ALGO 
#define DIM_TRIES 15
#define MAX_ZEROS 1000
#define MAX_FOLLOWING 2000

//#ifdef DEBUG_FOLLOW_ALGO 
//				printf("");
//#endif
int
gsl_monte_ftk_integrate (const gsl_monte_function * f,
                           const double xl[], const double xu[],
                           const size_t dim,
                           const size_t calls,
                           const size_t max_calls,
                           double min_err,
                           gsl_rng * r,
                           gsl_monte_ftk_state * state,
                           size_t *non_zeros, double *result, double *abserr)
{
  double vol, m = 0, q = 0;
  double *x = state->x;
  size_t n, i;

  size_t temp_non_zeros = 0 ;

  if (dim != state->dim)
    {
      GSL_ERROR ("number of dimensions must match allocated size", GSL_EINVAL);
    }

  for (i = 0; i < dim; i++)
    {
      if (xu[i] <= xl[i])
        {
          GSL_ERROR ("xu must be greater than xl", GSL_EINVAL);
        }

      if (xu[i] - xl[i] > GSL_DBL_MAX)
        {
          GSL_ERROR ("Range of integration is too large, please rescale",
                     GSL_EINVAL);
        }
    }

  /* Compute the volume of the region */

  vol = 1;

  for (i = 0; i < dim; i++)
    {
      vol *= xu[i] - xl[i];
    }


  for (n = 0; n < calls; n++)
    {
      /* Choose a random point in the integration region */

      for (i = 0; i < dim; i++)
        {
          x[i] = xl[i] + gsl_rng_uniform_pos (r) * (xu[i] - xl[i]);
        }

      {
        double fval = GSL_MONTE_FN_EVAL (f, x);

		if ( fval )
		{
			temp_non_zeros++; 
		}		

        /* recurrence for mean and variance */

        double d = fval - m;
        m += d / (n + 1.0);
        q += d * d * (n / (n + 1.0));
      }
    }

	
  *non_zeros = temp_non_zeros ; 
  *result = vol * m;

  if (calls < 2)
    {
      *abserr = GSL_POSINF;
    }
  else
    {
      *abserr = vol * sqrt (q / (calls * (calls - 1.0)));
    }

  return GSL_SUCCESS;
}

// Global follow variables TODO to be part of monte carlo state 
double xl_fol[10] ;
double xu_fol[10] ;
double x_fol[10] ;
int  amp_fol ; 
int next_fol ;
int tries ;

int zeros_in_a_row = 0 ;


void set_following ( double * x, double * xl, double * xu, int dim ){

	// This combination is saved
    memcpy ( xl_fol, xl, dim * sizeof(double) );	
    memcpy ( xu_fol, xu, dim * sizeof(double) );
    memcpy ( x_fol, x, dim * sizeof(double) );

	next_fol 	= 0 ; // Reset pointer to next dimension to search to	
	amp_fol 	= 0 ; //
	tries 		= 0 ; 
}

double get_following ( gsl_rng * r , int i ){


	if ( i == next_fol  )
	{
    	tries++;

		if ( tries < DIM_TRIES  )	
		{		
#ifdef DEBUG_FOLLOW_ALGO 
				printf("Try #%d for dimension %d\n", tries, i );
#endif
			return  xl_fol[i] + gsl_rng_uniform_pos (r) * (xu_fol[i] - xl_fol[i]);

		}
		else
		{
#ifdef DEBUG_FOLLOW_ALGO 
				printf("Now passing to next dimension : %d", i );
#endif
			tries = 0 ; 
			next_fol= (next_fol+1)%10 ; //10 is hardcoded TODO	
			return x_fol[i]; 
		}		
	}
	else
	{
		return x_fol[i]; 
	}		

}

int
gsl_monte_ftk_integrate_2 (const gsl_monte_function * f,
                         const double xl[], const double xu[],
                         const size_t dim,
                         const size_t calls,
                         const size_t max_calls,
                         double min_err,
                         gsl_rng * r,
                         gsl_monte_ftk_state * state,
                         size_t *steps_required, size_t *non_zeros, double *result, double *abserr)
{
    double vol, m = 0, q = 0;
    double *x = state->x;
    size_t n, i;
   
    // Follow Algorithm Variables, to be moved in the state
	//
    enum {FOLLOW, NOT_FOLLOW} mode;
	mode = NOT_FOLLOW ;
    int following_from = 0 ;  	
	size_t temp_non_zeros = 0;
	//
	/////////


    if (dim != state->dim)
    {
        GSL_ERROR ("number of dimensions must match allocated size", GSL_EINVAL);
    }
    
    for (i = 0; i < dim; i++)
    {
        if (xu[i] <= xl[i])
        {
            GSL_ERROR ("xu must be greater than xl", GSL_EINVAL);
        }
        
        if (xu[i] - xl[i] > GSL_DBL_MAX)
        {
            GSL_ERROR ("Range of integration is too large, please rescale",
                       GSL_EINVAL);
        }
    }
    
    /* Compute the volume of the region */
    
    vol = 1;
    
    for (i = 0; i < dim; i++)
    {
        vol *= xu[i] - xl[i];
    }
   

    //
    // SAMPLING PHASE
    //

    size_t tot_calls = calls;
    n = 2;

    double err_ratio = 1 ; 


    do  // Loop until the error is acceptable ( < min_err )
    {
    
		for (; n < tot_calls; n++)
        {

			/* Choose a random point in the integration region */


#ifdef FOLLOW_ALGO
			if ( mode == NOT_FOLLOW  )
			{
				for (i = 0; i < dim; i++)
				{
					x[i] = xl[i] + gsl_rng_uniform_pos (r) * (xu[i] - xl[i]);
				}
			}
			else 
			{	// If following
				for (i = 0; i < dim; i++)
				{
					x[i] = get_following( r, i ) ;
				}

			}
#else
			for (i = 0; i < dim; i++)
			{
					x[i] = xl[i] + gsl_rng_uniform_pos (r) * (xu[i] - xl[i]);
			}
#endif


            
            {
                double fval = GSL_MONTE_FN_EVAL (f, x);
         
#ifdef FOLLOW_ALGO
				if ( fval )
				{
						zeros_in_a_row = 0 ; 
						temp_non_zeros++; 

						if ( mode == NOT_FOLLOW  )
						{
							mode = FOLLOW ; 
							
							set_following( x, (double *) xl, (double *) xu, dim );  
							following_from = 0 ; 								
#ifdef DEBUG_FOLLOW_ALGO 
							printf("Following a Non Zero\n" );
#endif
						}
						else 
						{
							following_from++;
							
							if( zeros_in_a_row > MAX_ZEROS || following_from > MAX_FOLLOWING  )
							{
								// printf("Following too much ! Stop following.  \n") ;
#ifdef DEBUG_FOLLOW_ALGO 
								printf("Stop Following: %d ( > %d ) zeros in a row, following from %d ( > %d )steps \n", zeros_in_a_row, MAX_ZEROS , following_from, MAX_FOLLOWING );
#endif
								
								mode = NOT_FOLLOW ; 
							}		
							else						
							{
#ifdef DEBUG_FOLLOW_ALGO 
								printf("Continuingo Following: %d , following from %d steps \n", zeros_in_a_row, following_from );
#endif
							
							}

							// If continuing following
							//else
							//{
							//}    If I am already following, Then I might want to search between the two 		

						}	
				}
				else
				{
#ifdef DEBUG_FOLLOW_ALGO 
								printf("Zero This TIme ! \n ");
#endif
						zeros_in_a_row++ ; 
				}		
#endif 	

                /* recurrence for mean and variance */
                
                double d = fval - m;
                m += d / (n + 1.0);
                q += d * d * (n / (n + 1.0));
            }
        }

		tot_calls += 50000; // Adding 50000 samples at every iteration  

		*steps_required = n ; 
		*result = vol * m;

		err_ratio = vol * sqrt (q / (n * (n - 1.0))) / *result ; 

		// printf( "Performed %d iterations, err_ratio %f\n", tot_calls, err_ratio) ;

	} while( err_ratio > min_err && n < max_calls  ) ; 
    
	//
	// Final Phase
	//
    
	*non_zeros = temp_non_zeros ; 
	
	if (calls < 2)
    {
        *abserr = GSL_POSINF;
    }
    else
    {
        *abserr = vol * sqrt (q / (n * (n - 1.0)));
    }
    
    return GSL_SUCCESS;
}


gsl_monte_ftk_state *
gsl_monte_ftk_alloc (size_t dim)
{
  gsl_monte_ftk_state *s =
    (gsl_monte_ftk_state *) malloc (sizeof (gsl_monte_ftk_state));

  if (s == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for state struct",
                     GSL_ENOMEM, 0);
    }

  s->x = (double *) malloc (dim * sizeof (double));

  if (s->x == 0)
    {
      free (s);
      GSL_ERROR_VAL ("failed to allocate space for working vector",
                     GSL_ENOMEM, 0);
    }

  s->dim = dim;

  return s;
}

/* Set some default values and whatever */

int
gsl_monte_ftk_init (gsl_monte_ftk_state * s)
{
  size_t i;

  for (i = 0; i < s->dim; i++)
    {
      s->x[i] = 0.0;
    }

  return GSL_SUCCESS;
}

void
gsl_monte_ftk_free (gsl_monte_ftk_state * s)
{

  if (s == 0 )
	return;
  //RETURN_IF_NULL (s);
  
  free (s->x);
  free (s);
}
