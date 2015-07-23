/**********************************************************************
 *
 * analmulti.c: Programs to analyze QTL data from a backcross
 *              These functions either call multiple analysis
 *              methods or do big simulations
 *
 * anal_all, sim_null, anal_multi, anal_multi2
 * R_anal_all, R_sim_null, R_anal_multi, R_anal_multi2
 *
 * Karl Broman, 7/23/01, 7/10/01 [originally 5/21/96, 5/27/96, 5/30/96]
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>

#include "sweep.h"
#include "analbc.h"
#include "analmulti.h"
#include "perm.h"
#include "simbc.h"
#include "mcmc.h"

/**********************************************************************
 *
 * sim_null:  Simulate under the null hypothesis (no QTLs) and
 *            then call anal_anova and anal_zeng
 *
 * Input:
 *
 *   n_ind =         number of individuals
 *   n_chr =         number of chromosomes
 *   n_mar =         vector giving number of markers per chromosome
 *   tot_mar =       total number of markers
 *   recfrac =       vector of length sum(n_mar-1) giving rec fracs
 *   n_cim =         number of different steps for CIM
 *   cim_steps =     vector of steps for CIM
 *   n_sim =         number of simulations
 *   maxlod =        the output (maximum LOD scores) size n_sim * (1+n_cim)
 *   iwork =         workspace of integers [size tot_mar * (n_ind+1)]
 *   dwork =         workspace of doubles [size n_ind + (tot_mar+2)^2
 *                             + 2*tot_mar+1]
 *
 **********************************************************************/

void sim_null(int n_ind, int n_chr, int *n_mar, int tot_mar,
          double *recfrac, int n_cim, int *cim_steps,
          int n_sim, double *maxlod, int *iwork, double *dwork)
{
  int i, j, k, r, err;
  double *phenotypes, *xpx, *lod, *rss;
  int *index, *genotypes;

  /* set up workspace */
  genotypes = iwork; /* length = n_ind * tot_mar */
  index = genotypes + n_ind * tot_mar; /* length = tot_mar */
  phenotypes = dwork; /* length = n_ind */
  xpx = phenotypes + n_ind; /* length = (tot_mar+2)^2 */
  lod = xpx+(tot_mar+2)*(tot_mar+2); /* length = tot_mar */
  rss = lod + tot_mar; /* length = tot_mar + 1 */

  /* set up index */
  for(i=0; i<tot_mar; i++) index[i] = i+1;

  for(i=0; i<n_sim; i++) {

    /* simulate genotype data */
    simbc_mar(n_ind, n_chr, n_mar, recfrac, genotypes);

    /* simulate phenotype data */
    for(j=0; j<n_ind; j++)
      phenotypes[j] = norm_rand();

    /* calculate X'X matrix */
    calc_xpx(n_ind, tot_mar+2, genotypes, phenotypes, xpx);

    /* perform ANOVA */
    anal_anova(n_ind, tot_mar, xpx, lod);

    /* find maximum */
    maxlod[i] = lod[0];
    for(j=1; j<tot_mar; j++)
      if(maxlod[i] < lod[j])
    maxlod[i] = lod[j];

    /* perform forward selection */
    forward(tot_mar, xpx, cim_steps[0], index, rss);

    for(j=0, r=n_sim; j<n_cim; j++, r += n_sim) {
      /* unsweep columns */
      if(j>0) sweep(xpx, tot_mar+2, index+cim_steps[j],
            cim_steps[j-1]-cim_steps[j], &err);

      /* perform CIM */
      anal_cim(n_ind, tot_mar, xpx, lod, index, cim_steps[j], 1);

      maxlod[i+r] = lod[0];
      for(k=1; k<tot_mar; k++)
    if(maxlod[i+r] < lod[k])
      maxlod[i+r] = lod[k];
    }

  }
}


void R_sim_null(int *n_ind, int *n_chr, int *n_mar, int *tot_mar,
        double *recfrac, int *n_cim, int *cim_steps,
        int *n_sim, double *maxlod, int *iwork, double *dwork)
{
  GetRNGstate();

  sim_null(*n_ind, *n_chr, n_mar, *tot_mar, recfrac, *n_cim,
       cim_steps, *n_sim, maxlod, iwork, dwork);


  PutRNGstate();
}


/**********************************************************************
 *
 * anal_all:  Call anal_anova, anal_cim, and forward
 *
 * Input:
 *
 *   n_ind =         number of individuals
 *   tot_mar =       total number of markers
 *   genotypes =     genotype matrix [size = n_ind x tot_mar]
 *   phenotypes =    phenotype vector [length = n_ind]
 *   xpx =           workspace of size (tot_mar+2)^2
 *   n_steps =       number of initial steps for performing CIM
 *   max_steps =     maximum number of steps in forward selection
 *   lod =           vector of size tot_mar to contain the ANOVA results
 *   lodcim =        vector of size tot_mar to contain the CIM results
 *   index =         vector of size max_steps to contain the marker indices
 *                   from forward selection
 *   rss =           vector of size max_steps to contain the RSS
 *                   from forward selection
 *
 **********************************************************************/

void anal_all(int n_ind, int tot_mar, int *genotypes,
          double *phenotypes, double *xpx, int n_steps,
          int max_steps, double *lod, double *lodcim,
          int *index, double *rss)
{
  int err;

  /* calculate X'X matrix */
  calc_xpx(n_ind, tot_mar+2, genotypes, phenotypes, xpx);

  /* perform ANOVA */
  anal_anova(n_ind, tot_mar, xpx, lod);

  /* perform forward selection */
  forward(tot_mar, xpx, max_steps, index, rss);

  /* unsweep columns from forward selection that won't be used in CIM */
  if(max_steps > n_steps)
    sweep(xpx, tot_mar+2, index+n_steps, max_steps-n_steps, &err);

  /* perform CIM */
  anal_cim(n_ind, tot_mar, xpx, lodcim, index, n_steps, 1);
}

void R_anal_all(int *n_ind, int *tot_mar, int *genotypes,
        double *phenotypes, double *xpx, int *n_steps,
        int *max_steps, double *lod, double *lodcim,
        int *index, double *rss)
{
  anal_all(*n_ind, *tot_mar, genotypes, phenotypes, xpx, *n_steps,
       *max_steps, lod, lodcim, index, rss);
}


/**********************************************************************
 *
 * anal_multi: Function to simulate a backcross and analyze it by
 *            ANOVA, CIM, and forward selection with BIC and permutation
 *            tests.
 *
 * INPUT:
 *
 *     n_ind =         number of individuals
 *     n_chr =         number of chromosomes
 *     n_mar =         number of markers per chromosome
 *     recfrac =       recombination fractions between markers
 *     n_sim =         number of simulations to perform
 *     n_qtl =
 *     qtl_chr =
 *     mar_to_left =
 *     recfrac_to_left =
 *     effect =
 *     sigma =
 *     n_cim =         number of different cim_steps to use
 *     cim_steps =     number of steps of f.s. to perform before CIM
 *                     (must be decreasing)
 *     n_bic =         number of different BIC multipliers to use
 *     bic_mult =      the BIC multipliers
 *     thresh =        LOD thresholds (ANOVA followed by the CIM thresholds)
 *     drop =          drop in LOD between inferred peaks
 *     n_perm =        number of permutations in permutation tests
 *     alpha =         signif. level in perm'n tests (actually level * n_perm)
 *     n_qtl_id =      vector of size (1+n_cim+n_bic+1) * n_sim
 *     chr_id =        vector of size (1+n_cim+n_bic+1) * sum(n_mar) * n_sim
 *     mar_id =        vector of same size as chr_id
 *     iwork =         int workspace of size tot_mar*(n_ind+2)
 *     dwork =         double workspace of size (tot_mar+2)^2 +
 *                     n_ind*(tot_mar+4)+tot_mar+1
 *
 **********************************************************************/

void anal_multi(int n_ind, int n_chr, int *n_mar, double *recfrac,
        int n_sim, int n_qtl, int *qtl_chr, int *mar_to_left,
        double *recfrac_to_left, double *effect, double sigma,
        int n_cim, int *cim_steps, int max_steps, int n_bic,
        double *bic_mult, double *thresh, double *drop,
        int n_perm, int alpha, int *n_qtl_id, int *chr_id,
        int *mar_id, int *iwork, double *dwork)
{
  int i, j, k, first_on_chr, s, tot_mar, size, sizesq, err;
  int *genotypes, *ijunk, *index;
  double *phenotypes, *xpx, *lod, *djunk;

  /* important numbers */
  tot_mar = 0;
  for(i=0; i<n_chr; i++) tot_mar += n_mar[i];
  size = tot_mar + 2;
  sizesq = size*size;

  /* reform workspace: ints */
  genotypes = iwork; /* size tot_mar * n_ind */
  ijunk = genotypes + tot_mar*n_ind; /* size tot_mar */
  index = ijunk + tot_mar; /* size tot_mar */

  /* reform workspace: doubles */
  phenotypes = dwork; /* size n_ind */
  xpx = phenotypes + n_ind; /* size (tot_mar+2)^2 */
  lod = xpx + sizesq; /* size tot_mar+1 */
  djunk = lod + tot_mar+1; /* size n_ind*(tot_mar+3) */

  /* begin simulations */
  for(s=0; s<n_sim; s++) {
    /* simulate marker data */
    simbc_mar(n_ind, n_chr, n_mar, recfrac, genotypes);

    /* simulate phenotypes */
    simbc_qtl(n_ind, n_chr, n_mar, recfrac, genotypes, phenotypes,
          n_qtl, qtl_chr, mar_to_left, recfrac_to_left, effect,
          sigma);

    /* calc X'X matrix */
    calc_xpx(n_ind, size, genotypes, phenotypes, xpx);

    /* perform ANOVA and identify QTLs */
    anal_anova(n_ind, tot_mar, xpx, lod);
    identify_qtls(n_chr, n_mar, lod, thresh[0], drop[0],
          n_qtl_id+s, chr_id+tot_mar*s, mar_id+tot_mar*s,
          ijunk, djunk);
    /* perform forward selection (note: lod is really rss) */
    forward(tot_mar, xpx, max_steps, index, lod);
    /* minimize BIC */
    for(i=0; i<n_bic; i++)
      identify_qtls_bic(n_chr, n_mar, n_ind, index, lod, max_steps,
            n_qtl_id+s+n_sim*(1+n_cim+i),
            chr_id+s*tot_mar+tot_mar*n_sim*(1+n_cim+i),
            mar_id+s*tot_mar+tot_mar*n_sim*(1+n_cim+i),
            bic_mult[i]);
    /* perform CIM */
    /* unsweep columns if necessary */
    if(max_steps > cim_steps[0])
      sweep(xpx, size, index+cim_steps[0], max_steps-cim_steps[0], &err);
    for(i=0; i<n_cim; i++) {
      if(i>0)
    sweep(xpx, size, index+cim_steps[i], cim_steps[i-1]-cim_steps[i],
          &err);
      anal_cim(n_ind, tot_mar, xpx, lod, index, cim_steps[i], 1);
      identify_qtls(n_chr, n_mar, lod, thresh[i+1], drop[i+1],
            n_qtl_id+s+n_sim*(i+1),
            chr_id+s*tot_mar+tot_mar*n_sim*(1+i),
            mar_id+s*tot_mar+tot_mar*n_sim*(1+i),
            ijunk, djunk);
    }

    /* do forward selection with permutation tests */
    k = s+n_sim*(1+n_cim+n_bic);

    forw_perm(n_ind, tot_mar, genotypes, phenotypes, index,
          n_perm, alpha, n_qtl_id+k, djunk);

    /* turn index into chr_id and mar_id */
    for(i=0; i < n_qtl_id[k]; i++) {
      for(j=0, first_on_chr=1; j<n_chr; first_on_chr += n_mar[j], j++) {
    if(index[i] < first_on_chr + n_mar[j]) {
      chr_id[i+s*tot_mar+tot_mar*n_sim*(1+n_cim+n_bic)] = j+1;
      mar_id[i+s*tot_mar+tot_mar*n_sim*(1+n_cim+n_bic)] =
        index[i] - first_on_chr + 1;
      break;
    }
      }
    }

  } /* end of simulations */

}

void R_anal_multi(int *n_ind, int *n_chr, int *n_mar, double *recfrac,
          int *n_sim, int *n_qtl, int *qtl_chr, int *mar_to_left,
          double *recfrac_to_left, double *effect, double *sigma,
          int *n_cim, int *cim_steps, int *max_steps, int *n_bic,
          double *bic_mult, double *thresh, double *drop,
          int *n_perm, int *alpha, int *n_qtl_id, int *chr_id,
          int *mar_id, int *iwork, double *dwork)
{
  GetRNGstate();
  anal_multi(*n_ind, *n_chr, n_mar, recfrac, *n_sim, *n_qtl, qtl_chr,
         mar_to_left, recfrac_to_left, effect, *sigma, *n_cim,
         cim_steps, *max_steps, *n_bic, bic_mult, thresh, drop,
         *n_perm, *alpha, n_qtl_id, chr_id, mar_id, iwork, dwork);
  PutRNGstate();
}




/**********************************************************************
 *
 * anal_multi2: Function to simulate a backcross and analyze it by
 *              ANOVA, CIM, forward selection with BIC and permutation
 *              tests, and MCMC with BIC
 *
 * INPUT:
 *
 *     n_ind =         number of individuals
 *     n_chr =         number of chromosomes
 *     n_mar =         number of markers per chromosome
 *     recfrac =       recombination fractions between markers
 *     n_sim =         number of simulations to perform
 *     n_qtl =
 *     qtl_chr =
 *     mar_to_left =
 *     recfrac_to_left =
 *     effect =
 *     sigma =
 *     n_cim =         number of different cim_steps to use
 *     cim_steps =     number of steps of f.s. to perform before CIM
 *                     (must be decreasing)
 *     n_bic =         number of different BIC multipliers to use
 *     bic_mult =      the BIC multipliers
 *     thresh =        LOD thresholds (ANOVA followed by the CIM thresholds)
 *     drop =          drop in LOD between inferred peaks
 *     n_perm =        number of permutations in permutation tests
 *     alpha =         signif. level in perm'n tests (actually level * n_perm)
 *
 *
 *     n_mcmc =        number of steps in MCMC
 *     mcmc_delta =    delta value in MCMC
 *
 *     n_qtl_id =      vector of size (3+n_cim+n_bic) * n_sim
 *     chr_id =        vector of size (3+n_cim+n_bic) * sum(n_mar) * n_sim
 *     mar_id =        vector of same size as chr_id
 *     iwork =         int workspace of size tot_mar*(n_ind+4)+n_mcmc+4
 *     dwork =         double workspace of size n_ind*(tot_mar+4)+
 *                          (tot_mar+2)^2+tot_mar+n_mcmc+2
 *
 **********************************************************************/

void anal_multi2(int n_ind, int n_chr, int *n_mar, double *recfrac,
         int n_sim, int n_qtl, int *qtl_chr, int *mar_to_left,
         double *recfrac_to_left, double *effect, double sigma,
         int n_cim, int *cim_steps, int max_steps, int n_bic,
         double *bic_mult, double *thresh, double *drop,
         int n_perm, int alpha, int n_mcmc, double mcmc_delta,
         int *n_qtl_id, int *chr_id,
         int *mar_id, int *iwork, double *dwork)
{
  int i, j, k, first_on_chr, s, tot_mar, size, sizesq, err;
  int *genotypes, *ijunk, *index, *indicate, *ijunk2, ijunk3, *index2;
  double *phenotypes, *xpx, *lod, *djunk, *djunk2, djunk3;

  /* important numbers */
  tot_mar = 0;
  for(i=0; i<n_chr; i++) tot_mar += n_mar[i];
  size = tot_mar + 2;
  sizesq = size*size;

  /* reform workspace: ints */
  genotypes = iwork; /* size tot_mar * n_ind */
  ijunk = genotypes + tot_mar*n_ind; /* size tot_mar */
  index = ijunk + tot_mar; /* size tot_mar+1 */
  indicate = index + tot_mar + 1; /* size tot_mar + 1 */
  ijunk2 = indicate + tot_mar + 1; /* size n_mcmc + 1*/
  index2 = ijunk2 + n_mcmc + 1; /* size tot_mar + 1 */

  /* reform workspace: doubles */
  phenotypes = dwork; /* size n_ind */
  xpx = phenotypes + n_ind; /* size (tot_mar+2)^2 */
  lod = xpx + sizesq; /* size tot_mar+1 */
  djunk = lod + tot_mar+1; /* size n_ind*(tot_mar+3) */
  djunk2 = djunk + n_ind*(tot_mar+3); /* size n_mcmc + 1*/

  /* begin simulations */
  for(s=0; s<n_sim; s++) {
    /* simulate marker data */
    simbc_mar(n_ind, n_chr, n_mar, recfrac, genotypes);

    /* simulate phenotypes */
    simbc_qtl(n_ind, n_chr, n_mar, recfrac, genotypes, phenotypes,
          n_qtl, qtl_chr, mar_to_left, recfrac_to_left, effect,
          sigma);

    /* calc X'X matrix */
    calc_xpx(n_ind, size, genotypes, phenotypes, xpx);

    /* perform ANOVA and identify QTLs */
    anal_anova(n_ind, tot_mar, xpx, lod);
    identify_qtls(n_chr, n_mar, lod, thresh[0], drop[0],
          n_qtl_id+s, chr_id+tot_mar*s, mar_id+tot_mar*s,
          ijunk, djunk);
    /* perform forward selection (note: lod is really rss) */
    forward(tot_mar, xpx, max_steps, index, lod);
    /* minimize BIC */
    for(i=0; i<n_bic; i++)
      identify_qtls_bic(n_chr, n_mar, n_ind, index, lod, max_steps,
            n_qtl_id+s+n_sim*(1+n_cim+i),
            chr_id+s*tot_mar+tot_mar*n_sim*(1+n_cim+i),
            mar_id+s*tot_mar+tot_mar*n_sim*(1+n_cim+i),
            bic_mult[i]);

    /* save model with bic_mult[0] as starting point for MCMC */
    for(i=0; i<tot_mar+1; i++) indicate[i] = 0;
    /* I comment out the following, to start with the null model
    indicate[0] = n_qtl_id[s+n_sim*(1+n_cim)];
    for(i=0; i<n_qtl_id[s+n_sim*(1+n_cim)]; i++)
      indicate[index[i]] = 1; */


    /* perform CIM */
    /* unsweep columns if necessary */
    if(max_steps > cim_steps[0])
      sweep(xpx, size, index+cim_steps[0], max_steps-cim_steps[0], &err);
    for(i=0; i<n_cim; i++) {
      if(i>0)
    sweep(xpx, size, index+cim_steps[i], cim_steps[i-1]-cim_steps[i],
          &err);
      anal_cim(n_ind, tot_mar, xpx, lod, index, cim_steps[i], 1);
      identify_qtls(n_chr, n_mar, lod, thresh[i+1], drop[i+1],
            n_qtl_id+s+n_sim*(i+1),
            chr_id+s*tot_mar+tot_mar*n_sim*(1+i),
            mar_id+s*tot_mar+tot_mar*n_sim*(1+i),
            ijunk, djunk);
    }

    /* do forward selection with permutation tests */
    k = s+n_sim*(1+n_cim+n_bic);

    forw_perm(n_ind, tot_mar, genotypes, phenotypes, index,
          n_perm, alpha, n_qtl_id+k, djunk);

    /* turn index into chr_id and mar_id */
    for(i=0; i < n_qtl_id[k]; i++) {
      for(j=0, first_on_chr=1; j<n_chr; first_on_chr += n_mar[j], j++) {
    if(index[i] < first_on_chr + n_mar[j]) {
      chr_id[i+s*tot_mar+tot_mar*n_sim*(1+n_cim+n_bic)] = j+1;
      mar_id[i+s*tot_mar+tot_mar*n_sim*(1+n_cim+n_bic)] =
        index[i] - first_on_chr + 1;
      break;
    }
      }
    }

    /* do MCMC */
    k = s+n_sim*(1+n_cim+n_bic+1);

    mcmc_ms(n_ind, tot_mar, genotypes, phenotypes, xpx,
        n_mcmc, n_qtl_id + k, index,
        &djunk3, &ijunk3, index2,
        indicate, mcmc_delta, ijunk2, djunk2);

    /* turn index into chr_id and mar_id */
    for(i=0; i < n_qtl_id[k]; i++) {
      for(j=0, first_on_chr=1; j<n_chr; first_on_chr += n_mar[j], j++) {
    if(index[i] < first_on_chr + n_mar[j]) {
      chr_id[i+s*tot_mar+tot_mar*n_sim*(1+n_cim+n_bic+1)] = j+1;
      mar_id[i+s*tot_mar+tot_mar*n_sim*(1+n_cim+n_bic+1)] =
        index[i] - first_on_chr + 1;
      break;
    }
      }
    }


  } /* end of simulations */

}

void R_anal_multi2(int *n_ind, int *n_chr, int *n_mar, double *recfrac,
           int *n_sim, int *n_qtl, int *qtl_chr, int *mar_to_left,
           double *recfrac_to_left, double *effect, double *sigma,
           int *n_cim, int *cim_steps, int *max_steps, int *n_bic,
           double *bic_mult, double *thresh, double *drop,
           int *n_perm, int *alpha, int *n_mcmc,
           double *mcmc_delta, int *n_qtl_id, int *chr_id,
           int *mar_id, int *iwork, double *dwork)
{
  GetRNGstate();
  anal_multi2(*n_ind, *n_chr, n_mar, recfrac, *n_sim, *n_qtl, qtl_chr,
          mar_to_left, recfrac_to_left, effect, *sigma, *n_cim,
          cim_steps, *max_steps, *n_bic, bic_mult, thresh, drop,
          *n_perm, *alpha, *n_mcmc, *mcmc_delta,
          n_qtl_id, chr_id, mar_id, iwork, dwork);
  PutRNGstate();
}


/* end of analmulti.c */
