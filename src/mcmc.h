/**********************************************************************
 *
 * mcmc_ms.h: MCMC program to do model selection in order to analyze QTL
 *            data from a backcross
 *
 *   mcmc_ms
 *
 * Karl Broman, 7/14/01 [originally 8/27/96, 9/13/96, 9/14/96]
 *
 **********************************************************************/

void mcmc_ms(int n_ind, int tot_mar,
             int *genotypes, double *phenotypes, double *xpx,
             int n_steps, int *n_qtl_id, int *qtl_id,
             double *neg_log_post, int *first_seen, int *index,
         int *indicate, double delta, int *n_qtl_id_list,
         double *post_list);

void R_mcmc_ms(int *n_ind, int *tot_mar,
           int *genotypes, double *phenotypes, double *xpx,
           int *n_steps, int *n_qtl_id, int *qtl_id,
           double *neg_log_post, int *first_seen, int *index,
           int *indicate, double *delta, int *n_qtl_id_list,
           double *post_list);

/* end of mcmc.h */
