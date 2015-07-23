/**********************************************************************
 *
 * perm.c:  Program to perform forward selection using a
 *          permutation test to determine the size of the model
 *
 *    forw_perm, fit1, forw1, R_forw_perm
 *    double_permute, random_int, int_permute
 *
 * Karl W Broman, 7/26/01, 7/09/01
 *                [originally 6/9/96, 6/10/96, 6/20/96, 9/15/96]
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>

#include "analbc.h"
#include "perm.h"
#include "simbc.h"
#include "sweep.h"
#define TOL 1e-10

/**********************************************************************
 *
 * forw_perm:   Forward selection with permutation tests
 *              Needs to be surrounded with calls to GetRNGstate()
 *              and PutRNGstate()
 *
 * input:
 *
 *   n_ind            = number of progeny
 *
 *   tot_mar          = total number of markers
 *
 *   genotypes        = matrix of size n_ind * tot_mar
 *                      containing the genotypes
 *
 *   phenotypes       = vector of phenotypes
 *
 *   index            = vector of length tot_mar
 *
 *   n_perm           = number of permutation test simulations to do at each
 *                      stage
 *
 *   alpha            = significance level for testing at each stage
 *                      [really signif level * rss]
 *
 *   n_chosen         = number of QTLs identified in the end
 *
 *   dwork            = empty matrix of doubles of size
 *                      n_ind *(tot_mar + 3)
 *
 **********************************************************************/

void R_forw_perm(int *n_ind, int *tot_mar, int *genotypes,
         double *phenotypes, int *index, int *n_perm,
         int *alpha, int *n_chosen, double *dwork)
{
  GetRNGstate();

  forw_perm(*n_ind, *tot_mar, genotypes, phenotypes, index,
        *n_perm, *alpha, n_chosen, dwork);

  PutRNGstate();

}

void forw_perm(int n_ind, int tot_mar, int *genotypes, double *phenotypes,
           int *index, int n_perm, int alpha, int *n_chosen,
           double *dwork)
{
  int i, j, k;
  int below, above, n_below, n_above;
  int best;
  double rss, minrss, rsstemp;
  double *res, *gen, *phe, *phet;

  *n_chosen = 0;
  below = alpha; above = n_perm - alpha;

  /* set up index vector */
  for(i=0; i<tot_mar; i++) index[i] = i+1;

  /* set up workspace */
  res = dwork;     /* length = n_ind */
  phe = res + n_ind; /* length = n_ind */
  phet = phe + n_ind; /* length = n_ind */
  gen = phet + n_ind; /* length = n_ind * tot_mar */

  /* re-center phenotypes and genotypes */
  /* here, rss is used for phenotype or genotype average */
  rss = minrss = 0.0;
  for(i=0; i<n_ind; i++) rss += phenotypes[i];
  rss /= (double)n_ind;
  for(i=0; i<n_ind; i++) {
    phe[i] = phenotypes[i] - rss;
    minrss += (phe[i]*phe[i]);
  }
  for(j=0; j<tot_mar; j++) {
    rss = 0;
    for(i=0; i<n_ind; i++) rss += (double)genotypes[i+j*n_ind];
    rss /= (double)n_ind;
    for(i=0; i<n_ind; i++)
      gen[i+j*n_ind] = (double)genotypes[i+j*n_ind]-rss;
  }


  for(i=0; i<tot_mar; i++) {

    rsstemp = minrss;

    /* one step of forward selection */
    forw1(n_ind, tot_mar-i, gen, phe, index+i, res, &minrss, &best);
    best += i;

    for(j=0; j<n_ind; j++) phet[j] = phe[j];

    n_below = n_above = 0;

    for(j=0; j<n_perm; j++) {
      /* permute residuals */
      double_permute(phet,n_ind);

      rss = rsstemp;
      forw1(n_ind, tot_mar-i, gen, phet, index+i, res, &rss, &k);

      /* compare RSS to cur_rss */
      if(rss < minrss) /* better than observed */
    n_below++;
      else /* worse than observed */
    n_above++;

      if(n_below >= below) return; /* permutation test failed */
      if(n_above > above) /* go to next stage */
    j=n_perm;
    }

    /* re-order index */
    k = index[best];
    index[best] = index[i];
    index[i] = k;
    (*n_chosen)++;

    /* re-center phenotypes and remaining genotypes columns */
    fit1(n_ind, gen+n_ind*(index[i]-1), phe, res, &rss);
    for(j=0; j<n_ind; j++) phe[j] = res[j];
    for(j=i+1; j<tot_mar; j++) {
      fit1(n_ind, gen+n_ind*(index[i]-1), gen+n_ind*(index[j]-1),
       res, &rss);
      for(k=0; k<n_ind; k++)
    gen[k+n_ind*(index[j]-1)] = res[k];
    }
  }
}


/* fit1: calculate the residuals from regressing y on x */

void fit1(int n_ind, double *x, double *y, double *res, double *rss)
{
  int i;
  double sxy, sxx;

  sxy = sxx = 0.0;
  for(i=0; i<n_ind; i++) {
    sxy += (x[i] * y[i]);
    sxx += (x[i] * x[i]);
  }
  if(sxx < TOL) { /* linear combination of markers already chosen */
    *rss = -1.0;
    return;
  }
  sxy /= sxx; *rss=0.0;
  for(i=0; i<n_ind; i++) {
    res[i] = y[i] - sxy * x[i];
    *rss += (res[i] * res[i]);
  }
}


/* forw1: one step of forward selection */
/* minrss should be set at some reasonable maximum */
void forw1(int n_ind, int n_mar, double *x, double *y,
       int *index, double *res, double *minrss,
       int *best)
{
  int i;
  double rss;

  /* the following is not necessary */
  /* *minrss = 0.0; */
  /* for(i=0; i<n_ind; i++) *minrss += (y[i]*y[i]); */

  for(i=0; i<n_mar; i++) {
    fit1(n_ind, x+n_ind*(index[i]-1), y, res, &rss);

    if(rss > 0 && rss < *minrss) {
      *minrss = rss;
      *best = i;
    }
  }
}




/* function for permuting residuals; requires R's unif_rand()
   this needs to be placed between calls to GetRNGstate()
   and putRNGstate() */

void double_permute(double *array, int len)
{
  int i, which;
  double tmp;

  for(i=0; i < len; i++) {
    which = random_int(i, len-1);
    tmp = array[which];
    array[which] = array[i];
    array[i] = tmp;
  }
}

void int_permute(int *array, int len)
{
  int i, which;
  int tmp;

  for(i=0; i < len; i++) {
    which = random_int(i, len-1);
    tmp = array[which];
    array[which] = array[i];
    array[i] = tmp;
  }
}

int random_int(int low, int high)
{
  return((int)(unif_rand()*(double)(high - low + 1)) + low);
}


/* end of perm.c */
