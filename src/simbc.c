/**********************************************************************
 * 
 * simbc.c:   Programs to simulate and analyze QTL data from a backcross,
 *            with unequal numbers of markers per chromosome, and at 
 *            unequal spacing
 *
 *   
 *
 * Karl Broman, 7/6/01 [originally 9/30/96, 10/7/96, 10/8/96]
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>

#include "simbc.h"

/**********************************************************************
 * 
 * simbc_mar
 *
 * Simulate marker data from a backcross (All chromosomes are of the
 *   same length, and with the same number of markers at constant
 *   spacing.  No crossover interference is assumed.) 
 *
 *
 * input:
 *
 *   n_progeny =      number of progeny
 *
 *   n_chromosomes =  number of chromosomes
 *
 *   n_markers =      vector (len n_chr) giving number of markers on each 
 *                    chromosome
 *
 *   recfrac =        vector [length sum(n_mar-1)] giving recombination
 *                    fractions between markers
 *
 *   genotypes =      empty matrix of size n_progeny * sum(n_markers), 
 *                    in which the genotypes will be placed (0's or 1's).  
 *                    Data will be stored by column, with each column 
 *                    corresponding to a different marker, and each row 
 *                    corresponding to a different individual
 * 
 **********************************************************************/

void simbc_mar(int n_progeny, int n_chromosomes, int *n_markers,
	       double *recfrac, int *genotypes)
{
  int i, j, k, mar;

  for(i=0, mar=0; i < n_chromosomes; i++) {
    /* first marker on chromosome */
    for(k=0; k < n_progeny; k++) {
      if(unif_rand() > 0.5) genotypes[k + mar * n_progeny] = 1;
      else genotypes[k + mar * n_progeny] = 0;
    }
    for(j=1, mar++; j < n_markers[i]; j++, mar++) {
      /* convert marker spacing from distance to r.f. */
      for(k=0; k < n_progeny; k++) {
	if(unif_rand() < recfrac[mar-i-1])
	  genotypes[k + mar * n_progeny] = 
	    1 - genotypes[k + (mar-1) * n_progeny];
	else
	  genotypes[k + mar * n_progeny] = 
	    genotypes[k + (mar-1) * n_progeny];
	  
      }
    }
  }
}
      
  
/**********************************************************************
 * 
 * simbc_qtl
 *
 * Simulate marker and qtl data from a backcross 
 *
 * input:
 *
 *   n_progeny =      number of progeny
 *
 *   n_chromosomes =  number of chromosomes
 *
 *   n_markers =      vector (len n_chr) giving number of markers on each 
 *                    chromosome
 *
 *   recfrac =        vector [length sum(n_mar-1)] giving recombination 
 *                    fractions between markers
 *
 *   genotypes =      empty matrix of size n_progeny * sum(n_markers), 
 *                    in which the genotypes will be placed (0's or 1's).  
 *                    Data will be stored by column, with each column 
 *                    corresponding to a different marker, and each row 
 *                    corresponding to a different individual
 * 
 *   phenotypes =     vector of length n_progeny, in which the phenotypes
 *                    will be placed.
 *
 *   n_qtl =          number of QTLs to simulate
 *
 *   qtl_chr =        number of the chrom on which the QTL sit
 *
 *   mar_to_left =    marker number to left of QTL (vec of len n_qtl)
 *                    (cumulative number: 0, 1, ..., sum(n_mar)-1 )
 *
 *   recfrac_to_left =recombination fraction (in M) between QTL and marker
 *                    to left (vector of length n_qtl)
 *
 *   effect =         effect of QTL (vector of length n_qtl)
 *
 *   sigma =          SD of environmental variation (noise)
 *
 **********************************************************************/

void simbc_qtl(int n_progeny, int n_chromosomes, int *n_markers,
	       double *recfrac, int *genotypes, 
	       double *phenotypes, int n_qtl, int *qtl_chr,
	       int *mar_to_left, double *recfrac_to_left, 
	       double *effect, double sigma)
{
  int i, k;
  double r, theta_left, theta_right;
  double r_r, r_nr;

  /* simulate marker data */
  simbc_mar(n_progeny, n_chromosomes, n_markers,
	    recfrac, genotypes);

  /* simulate environmental variation */
  for(k=0; k< n_progeny; k++) 
    phenotypes[k] =  norm_rand() * sigma;

  for(i=0; i < n_qtl; i++) {
    /* rec. frac to left and right of QTL */
    theta_left = recfrac_to_left[i];
    theta_right = 0.5*(1.0-(1.0-2.0*recfrac[mar_to_left[i] - qtl_chr[i] - 1])/
		       (1.0-2.0*theta_left));

    /* get cond'l prob of QTL genotype given marker genotypes */
    r_r = theta_left * (1.0-theta_right);
    r_r = r_r / (r_r + theta_right * (1.0-theta_left));
    r_nr = theta_left * theta_right;
    r_nr = r_nr / (r_nr + (1.0 - theta_left)*(1.0-theta_right));

    /* simulate QTL genotypes */
    r = unif_rand();
    for(k=0; k< n_progeny; k++) {
      if(genotypes[k + mar_to_left[i] * n_progeny]) {
	if(genotypes[k + (mar_to_left[i]-1) * n_progeny]) { 
	  /* both markers are 1 */
	  if(r > r_nr) /* non recombinant : QTL = 1 */
	    phenotypes[k] += effect[i];
	}
	else {
	  /* mar to left is 1; mar to right is 0 */
	  if(r > r_r) /* recomb in right interval: QTL = 1 */
	    phenotypes[k] += effect[i];
	}
      }
      else {
	if(genotypes[k + (mar_to_left[i]-1) * n_progeny]) { 
	  /* mar to left is 0; mar to right is 1 */
	  if(r < r_r) /* recomb in left interval: QTL = 1 */
	    phenotypes[k] += effect[i];
	}
	else {
	  /* both markers are 0 */
	  if(r < r_nr) /* double recombinant : QTL = 1 */
	    phenotypes[k] += effect[i];
	}
      }
    }
  }
}
      
  
void R_simbc_qtl(int *n_progeny, int *n_chromosomes, int *n_markers,
		 double *recfrac, int *genotypes, 
		 double *phenotypes, int *n_qtl, int *qtl_chr,
		 int *mar_to_left, double *recfrac_to_left, 
		 double *effect, double *sigma)
{

  GetRNGstate();

  simbc_qtl(*n_progeny, *n_chromosomes, n_markers,
	    recfrac, genotypes, phenotypes, *n_qtl, qtl_chr,
	    mar_to_left, recfrac_to_left, effect, *sigma);

  PutRNGstate();
}


/* end of simbc.c */
