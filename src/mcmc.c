/**********************************************************************
 *
 * mcmc_ms.c: MCMC program to do model selection in order to analyze QTL 
 *            data from a backcross
 *
 *   mcmc_ms
 *
 * Karl Broman, 7/14/01 [originally 8/27/96, 9/13/96, 9/14/96]
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>

#include "mcmc.h"
#include "sweep.h"
#include "analbc.h"
#include "perm.h"

/**********************************************************************
 * 
 * mcmc_ms:  model selection on data from a backcross, using MCMC
 *
 * input:
 * 
 *   n_ind        = number of progeny
 *
 *   tot_markers      = total number of markers (n_mar * n_chr)
 *
 *   genotypes        = matrix of size n_ind * tot_markers,
 *                      containing the genotypes (0's or 1's).  
 *                      Data should be stored by column
 *                      with each column corresponding to a different 
 *                      marker, and each row corresponding to a
 *                      different individual
 * 
 *   phenotypes       = vector of length n_ind containing the
 *                      simulated phenotypes
 *
 *   xpx              = empty matrix of size tot_markers + 2,
 *                      used as workspace
 *
 *   n_steps          = number of MCMC steps to take
 *
 *   n_qtl_id         = number of QTLs identified
 *
 *   qtl_id           = vector giving (on output) marker numbers 
 *                      identified as QTLs (length at least tot_markers)
 *
 *   neg_log_post     = negative log posterior prob (up to scale factor) 
 *                      of best model
 *   
 *   first_seen       = step on which the best model was first seen
 *
 *   index            = vector of length tot_markers + 1, used as a
 *                      workspace
 *
 *   indicate         = vector of indicators saying which of the
 *                      markers are in the INITIAL model (vector 
 *                      of length tot_markers + 1)
 *    
 *   delta            = parameter used in the prior on the model size;
 *                      larger value -> fewer markers in the model
 *                      (log of previously used delta)
 *
 *   n_qtl_id_list    = empty vector of length n_steps + 1
 *                      on output, contains no qtl in model at each step
 *
 *   post_list        = empty vector of length n_steps + 1
 *                      on output, contains all of the neg_log_post's
 *
 **********************************************************************/

void mcmc_ms(int n_ind, int tot_mar,
             int *genotypes, double *phenotypes, double *xpx, 
             int n_steps, int *n_qtl_id, int *qtl_id,
             double *neg_log_post, int *first_seen, int *index, 
	     int *indicate, double delta, int *n_qtl_id_list,
	     double *post_list)
{
  int size, sizesqm1, err;
  int i, j, k, flag;
  double cur_neg_log_post, cur_rss, next_rss;
  double prob;
  double t1, t2, t3;

  size = tot_mar + 2;
  sizesqm1 = size * size - 1;

  /* calculate x'x matrix */
  calc_xpx(n_ind, size, genotypes, phenotypes, xpx);

  /* set up index for the initial sweep */
  j = 0;
  for(i=1; i <= tot_mar; i++) {
    if(indicate[i]) {
      index[j] = i;
      j++;
    }
  }

  /* sweep initial set of markers */
  sweep(xpx, size, index, j, &err);
  
  /* save rss for initial model */
  *neg_log_post = log(xpx[sizesqm1]) + (double)indicate[0] * delta /
    (double)(n_ind);
  /*  Rprintf("    %5d %2d %10.5lf %9.5lf %9.5lf\n", 0, indicate[0], 
      xpx[sizesqm1], *neg_log_post, *neg_log_post);  */
  *n_qtl_id = indicate[0];
  *first_seen = 0;
  for(j=1, k=0; j<= tot_mar; j++) {
    if(indicate[j]) {
      qtl_id[k] = j;
      k++;
    }
  }
  post_list[0] = *neg_log_post;
  n_qtl_id_list[0] = indicate[0];

  /* index contains numbers 1, ..., tot_mar */
  for(j=0; j< tot_mar; j++) index[j] = j+1;

  /* now begin MCMC */
  for(i=1; i<= n_steps; i++) {
    
    /* permute list of markers */
    int_permute(index, tot_mar);

    for(j=0; j< tot_mar; j++) {
      /* calculate rss with and without this marker in the model */
      cur_rss = xpx[sizesqm1];
      sweep(xpx, size, index+j, 1, &err);
      next_rss = xpx[sizesqm1];

      /* calculate prob of keeping it or putting it into the model */
      t1 = -0.5 * (double)(n_ind) * log(next_rss);
      t2 = -0.5 * (double)(n_ind) * log(cur_rss);
      t3 = 0.5 * delta;
      if(indicate[index[j]])  
	prob = 1.0/(1.0 + exp(t1 + t3 - t2));
      else 
	prob = 1.0/(1.0 + exp(t2 + t3 - t1));

      /* take a bernoulli draw */
      if(unif_rand() < prob) { /* variable should be in model */
        /* var not yet in model */
        if(!indicate[index[j]]) { 
          indicate[index[j]] = 1;
          indicate[0]++;
        }
        /* variable already in model */
        else {
          sweep(xpx, size, index+j, 1, &err);
        }
      }
      else { /* variable should not be in model */
        /* var currently in model */
        if(indicate[index[j]]) { 
          indicate[index[j]] = 0;
          indicate[0]--;
        }
        /* variable not currently in model */
        else {
          sweep(xpx, size, index+j, 1, &err);
        }
      }
    }

    /* better than currest best model? */
    cur_neg_log_post = log(xpx[sizesqm1]) + (double)indicate[0] * delta /
      (double) n_ind;
    post_list[i] = cur_neg_log_post;
    n_qtl_id_list[i] = indicate[0];

    /*    Rprintf("    %5d %2d %10.5lf %9.5lf %9.5lf\n", i, indicate[0], 
	  xpx[sizesqm1], cur_neg_log_post, *neg_log_post);   */
    if(cur_neg_log_post < *neg_log_post) {  /* new model is an improvement */
      *neg_log_post = cur_neg_log_post;

      /* check that new model really is different */
      flag = 0;
      if(*n_qtl_id == indicate[0]) {
	for(j=1, k=0; j<= tot_mar; j++) {
	  if(k < *n_qtl_id && qtl_id[k] == j) { 
	    if(!indicate[j]) {
	      flag = 1;
	      break;
	    }
	    k++;
	  }
	  else
	    if(indicate[j]) {
	      flag = 1;
	      break;
	    }
	}
      }	
      else flag = 1;

      if(flag) {  /* new model really is different */
	*n_qtl_id = indicate[0];
	*first_seen = i;
	for(j=1, k=0; j<= tot_mar; j++) {
	  if(indicate[j]) {
	    qtl_id[k] = j;
	    k++;
	  }
	}
      }

    }
  }
  /*  Rprintf("\t%d\n", *n_qtl_id); */
}      


void R_mcmc_ms(int *n_ind, int *tot_mar,
	       int *genotypes, double *phenotypes, double *xpx, 
	       int *n_steps, int *n_qtl_id, int *qtl_id,
	       double *neg_log_post, int *first_seen, int *index, 
	       int *indicate, double *delta, int *n_qtl_id_list,
	       double *post_list)
{
  GetRNGstate();
  mcmc_ms(*n_ind, *tot_mar, genotypes, phenotypes, xpx, 
          *n_steps, n_qtl_id, qtl_id, neg_log_post, first_seen, 
	  index, indicate, *delta, n_qtl_id_list, post_list);
  PutRNGstate();
}

/* end of mcmc.c */
