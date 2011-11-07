/**********************************************************************
 * 
 * analmulti.h: Programs to analyze QTL data from a backcross
 *              These functions either call multiple analysis
 *              methods or do big simulations
 *
 * anal_all, sim_null, anal_multi, anal_multi2
 * R_anal_all, R_sim_null, R_anal_multi, R_anal_multi2
 *
 * Karl Broman, 7/23/01, 7/10/01 [originally 5/21/96, 5/27/96, 5/30/96]
 *
 **********************************************************************/

void sim_null(int n_ind, int n_chr, int *n_mar, int tot_mar, 
	      double *recfrac, int n_cim, int *cim_steps, 
	      int n_sim, double *maxlod, int *iwork, double *dwork);

void R_sim_null(int *n_ind, int *n_chr, int *n_mar, int *tot_mar, 
		double *recfrac, int *n_cim, int *cim_steps, 
		int *n_sim, double *maxlod, int *iwork, double *dwork);

void anal_all(int n_ind, int tot_mar, int *genotypes,
	      double *phenotypes, double *xpx, int n_steps,
	      int max_steps, double *lod, double *lodcim, 
	      int *index, double *rss);

void R_anal_all(int *n_ind, int *tot_mar, int *genotypes,
		double *phenotypes, double *xpx, int *n_steps,
		int *max_steps, double *lod, double *lodcim, 
		int *index, double *rss);

void anal_multi(int n_ind, int n_chr, int *n_mar, double *recfrac,
		int n_sim, int n_qtl, int *qtl_chr, int *mar_to_left, 
		double *recfrac_to_left, double *effect, double sigma, 
		int n_cim, int *cim_steps, int max_steps, int n_bic, 
		double *bic_mult, double *thresh, double *drop, 
		int n_perm, int alpha, int *n_qtl_id, int *chr_id, 
		int *mar_id, int *iwork, double *dwork) ;

void R_anal_multi(int *n_ind, int *n_chr, int *n_mar, double *recfrac,
		  int *n_sim, int *n_qtl, int *qtl_chr, int *mar_to_left, 
		  double *recfrac_to_left, double *effect, double *sigma, 
		  int *n_cim, int *cim_steps, int *max_steps, int *n_bic, 
		  double *bic_mult, double *thresh, double *drop, 
		  int *n_perm, int *alpha, int *n_qtl_id, int *chr_id, 
		  int *mar_id, int *iwork, double *dwork);

void anal_multi2(int n_ind, int n_chr, int *n_mar, double *recfrac,
		 int n_sim, int n_qtl, int *qtl_chr, int *mar_to_left, 
		 double *recfrac_to_left, double *effect, double sigma, 
		 int n_cim, int *cim_steps, int max_steps, int n_bic, 
		 double *bic_mult, double *thresh, double *drop, 
		 int n_perm, int alpha, int n_mcmc, double mcmc_delta,
		 int *n_qtl_id, int *chr_id, 
		 int *mar_id, int *iwork, double *dwork);

void R_anal_multi2(int *n_ind, int *n_chr, int *n_mar, double *recfrac,
		   int *n_sim, int *n_qtl, int *qtl_chr, int *mar_to_left, 
		   double *recfrac_to_left, double *effect, double *sigma, 
		   int *n_cim, int *cim_steps, int *max_steps, int *n_bic, 
		   double *bic_mult, double *thresh, double *drop, 
		   int *n_perm, int *alpha, int *n_mcmc, 
		   double *mcmc_delta, int *n_qtl_id, int *chr_id, 
		   int *mar_id, int *iwork, double *dwork);

/* end of analmulti.h */
