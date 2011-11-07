/**********************************************************************
 * 
 * perm.h:  Program to perform forward selection using a  
 *          permutation test to determine the size of the model
 *
 *    forw_perm, fit1, forw1, R_forw_perm
 *    double_permute, random_int, int_permute
 *
 * Karl W Broman, 7/26/01, 7/09/01 
 *                [originally 6/9/96, 6/10/96, 6/20/96, 9/15/96]
 *
 **********************************************************************/

void R_forw_perm(int *n_ind, int *tot_mar, int *genotypes, 
		 double *phenotypes, int *index, int *n_perm, 
		 int *alpha, int *n_chosen, double *dwork);

void forw_perm(int n_ind, int tot_mar, int *genotypes, double *phenotypes, 
	       int *index, int n_perm, int alpha, int *n_chosen, 
	       double *dwork);

void fit1(int n_ind, double *x, double *y, double *res, double *rss);

void forw1(int n_ind, int n_mar, double *x, double *y, 
	   int *index, double *res, double *minrss, 
	   int *best);

void double_permute(double *array, int len);

int random_int(int low, int high);

void int_permute(int *array, int len);

/* end of perm.h */
