/**********************************************************************
 * 
 * analbc.h: Programs to analyze QTL data from a backcross
 *
 *   calc_xpx, anal_anova, anal_cim, forward, 
 *   R_calc_xpx, 
 *   identify_qtls, piksrt, identify_qtls_bic,
 *   R_identify_qtls, R_identify_qtls_bic,     
 *   which_correct, R_which_incorrect
 *
 * Karl Broman, 7/13/01 [originally 5/21/96, 5/27/96, 5/30/96]
 *
 **********************************************************************/

/*
void R_calc_xpx(int *n_progeny, int *size, int *genotypes,
		double *phenotypes, double *xpx);
*/

void calc_xpx(int n_progeny, int size, int *genotypes, 
              double *phenotypes, double *xpx);

void anal_anova(int n_progeny, int tot_mar, double *xpx, double *lod);
    
void anal_cim(int n_progeny, int tot_mar, double *xpx, double *lod, 
	       int *index, int n_steps, int skip_forw);

void forward(int tot_mar, double *xpx, int max_steps,
             int *index, double *rss);

void identify_qtls(int n_chr, int *n_mar, double *lod,
                   double threshold, double drop, int *n_qtl_id,
                   int *chr_id, int *mar_id, int *lwork, 
                   double *dwork);

/*
void R_identify_qtls(int *n_chr, int *n_mar, double *lod,
		     double *threshold, double *drop, int *n_qtl_id,
		     int *chr_id, int *mar_id, int *lwork, 
		     double *dwork);
*/

void piksrt(int n, double *arr, int *larr);

void identify_qtls_bic(int n_chr, int *n_mar, int n_progeny,
                       int *index, double *rss, int n_models, 
                       int *n_qtl_id, int *chr_id, int *mar_id,
                       double multiplier);

/*
void R_identify_qtls_bic(int *n_chr, int *n_mar, int *n_progeny,
			 int *index, double *rss, int *n_models, 
			 int *n_qtl_id, int *chr_id, int *mar_id,
			 double *multiplier);
*/

void which_correct(int n_infer, int *chr_infer, int *mar_infer,
		   int n_true, int *chr_true, int *mar_true,
		   int within, int *correct, int *n_incorrect);

void R_which_correct(int *n_infer, int *chr_infer, int *mar_infer,
		     int *n_true, int *chr_true, int *mar_true,
		     int *within, int *correct, int *n_incorrect);

/* end of analbc.h */
