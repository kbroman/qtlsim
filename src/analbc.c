/**********************************************************************
 * 
 * analbc.c: Programs to analyze QTL data from a backcross
 *
 *   calc_xpx, anal_anova, anal_cim, forward, 
 *   R_calc_xpx, 
 *   identify_qtls, piksrt, identify_qtls_bic,
 *   R_identify_qtls, R_identify_qtls_bic,     
 *   which_correct, R_which_correct
 *
 * Karl Broman, 7/13/01 [originally 5/21/96, 5/27/96, 5/30/96]
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

/**********************************************************************
 * 
 * calc_xpx
 *
 *   Calculates (x y)' (x y) given genotype phenotype data.
 *   We actually form the "corrected" SSCP matrix, with the first
 *   column already "swept"
 *
 * input:
 * 
 *   n_progeny 
 *
 *   size       = tot_mar + 2
 *
 *   genotypes 
 *
 *   phenotypes
 *
 *   xpx        = matrix of size (size)^2 in which (x y)' (x y) will
 *                be placed
 *
 **********************************************************************/

/*
void R_calc_xpx(int *n_progeny, int *size, int *genotypes,
		double *phenotypes, double *xpx)
{
  calc_xpx(*n_progeny, *size, genotypes, phenotypes, xpx);
} 
*/

void calc_xpx(int n_progeny, int size, int *genotypes, 
              double *phenotypes, double *xpx)
{
  int i, j, k, p, q, r, sizesq, sizem1;
  double a;

  sizesq = size * size;
  sizem1 = size-1;

  /* place zeros in matrix */
  for(i=0; i<sizesq; i++) xpx[i] = 0.0;
      
  /* first column */
  xpx[0] = -1.0/(double)n_progeny;
  for(i=1,p=0; i<sizem1; i++, p+= n_progeny) {
    for(j=0; j<n_progeny; j++)
      xpx[i] += (double)genotypes[j+p];
    xpx[i] /= (double)n_progeny;
  }
  for(j=0; j<n_progeny; j++) 
    xpx[sizem1] += phenotypes[j];
  xpx[sizem1] /= (double)n_progeny;

  for(j=0; j<n_progeny; j++) {
    for(i=1, p=0, r=size; i<sizem1; i++, p+=n_progeny, r+=size) {
      a = ((double)genotypes[j+p]-xpx[i]);
      xpx[i+r] += (a*a);
      for(k=(i+1),q=p+n_progeny; k<sizem1; k++, q+=n_progeny) 
	xpx[k+r] += (a*((double)genotypes[j+q]-xpx[k]));
      xpx[sizem1+r] += (a*(phenotypes[j]-xpx[sizem1]));
    }
    a = (phenotypes[j]-xpx[sizem1]);
    xpx[sizesq-1] += (a*a);
  }

  /* fill out upper triangle by symmetry */
  for(i=0, p=0; i<sizem1; i++, p+=size) 
    for(j=i+1, q= p+size; j<size; j++, q+=size) 
      xpx[i + q] = xpx[j + p];
}


/**********************************************************************
 * 
 * anal_anova
 *
 *   QTL analysis of backcross data by a lander/botstein type method,
 *   assuming QTLs are located at marker loci (thus, we simply regress
 *   the phenotypes on one marker at a time)
 *
 * input:
 *
 *   n_progeny      = number of progeny
 *  
 *   tot_mar        = total number of markers
 *
 *   xpx            = empty square matrix of size (tot_mar + 2)^2
 *                    used as workspace
 *
 *   lod            = vector of length tot_mar in which log lik
 *                    ratio will be placed
 *
 **********************************************************************/

void anal_anova(int n_progeny, int tot_mar, double *xpx, double *lod)
{
  int i, j;
  int size, sizesqm1, err;
  
  size = tot_mar + 2;
  sizesqm1 = size*size - 1;
  
  /* log RSS under null model */
  j = 1;
  lod[0] = log(xpx[sizesqm1]);
  for(i=1; i< size - 2; i++) lod[i] = lod[0];
    
  /* sweep each of the other marker columns, one at a time  */
  for(i=1; i < size - 1; i++) {
    sweep(xpx, size, &i, j, &err);
    lod[i-1] -= log(xpx[sizesqm1]);
    sweep(xpx, size, &i, j, &err);
  }
  
  /* get lod */
  for(i=0; i<size-2; i++)
    lod[i] *= ((double)n_progeny/(2.0*log(10.0)));
}
    
/**********************************************************************
 * 
 * anal_cim
 *
 *   QTL analysis of backcross data by a Cim-type method,
 *   assuming QTLs are located at marker loci, and using forward
 *   selection to choose markers for use as covariates (choosing 
 *   "n_steps" markers) 
 *
 * input:
 *
 *   n_progeny      = number of progeny
 *  
 *   tot_mar        = total number of markers
 * 
 *   xpx            = empty square matrix of size (tot_mar + 2)
 *                    used as workspace
 *
 *   lod            = vector of length tot_mar in which log lik
 *                    ratio will be placed
 *
 *   index          = index vector of length tot_mar
 *                    used as workspace
 *
 *   n_steps        = number of markers to use to "control background"
 *
 **********************************************************************/

void anal_cim(int n_progeny, int tot_mar, double *xpx, double *lod, 
	       int *index, int n_steps, int skip_forw)
{
  int i;
  int size, sizesqm1, err;
  
  size = tot_mar + 2;
  sizesqm1 = size*size - 1;
  
  /* forward selection up to size n_steps */
  if(!skip_forw) 
    forward(tot_mar, xpx, n_steps, index, lod);

  /* the first column + columns indicated in the first "n_steps" 
     positions of "index" (those chosen by forward selection) are 
     swept.  The rest of the positions in "index" have not been chosen */

  /* we now subtract and add each of the columns in the selected set,
     and then add and subtract each of the other columns */

  /* columns in selected set */
  for(i=0; i< n_steps; i++) {
    lod[index[i]-1] = -log(xpx[sizesqm1]);
    sweep(xpx, size, index+i, 1, &err);
    lod[index[i]-1] += log(xpx[sizesqm1]);
    sweep(xpx, size, index+i, 1, &err);
  }

  /* columns not in selected set */
  for(i=n_steps; i<size-2; i++) {
    lod[index[i]-1] = log(xpx[sizesqm1]);
    sweep(xpx, size, index+i, 1, &err);
    lod[index[i]-1] -= log(xpx[sizesqm1]);
    sweep(xpx, size, index+i, 1, &err);
  }

  /* get lod */
  for(i=0; i<size-2; i++)
    lod[i] *= ((double)n_progeny/(2.0*log(10.0))); 
}
    
/**********************************************************************
 * 
 * forward:  performs forward selection
 *
 * input:
 *
 *   tot_mar          = total number of markers
 *
 *   xpx              = empty matrix of size n_chr * n_mar + 2,
 *                      used as workspace
 *
 *   max_steps        = maximum number of steps to take
 *
 *   index            = vector of length tot_mar; on
 *                      output contains the indices of the markers
 *                      in the order they're added
 *
 *   rss              = vector of length max_steps + 1, which on
 *                      output contains the residual sum of squares
 *                      for the models as in "index"
 *
 **********************************************************************/

void forward(int tot_mar, double *xpx, int max_steps,
             int *index, double *rss)
{
  int i, j, k, size, err, sizesqm1;
  double a;

  size = tot_mar + 2;
  sizesqm1 = size * size - 1;

  /* sweep first column 
  i=0; 
  sweep(xpx, size, &i, 1, &err); */

  /* rss for model with no markers */
  rss[0] = xpx[sizesqm1];

  /* set up index */
  for(i=0; i<size - 2; i++) index[i]=i+1;

  for(i=0; i< max_steps; i++) {  
    rss[i+1] = xpx[sizesqm1];
    for(j=i; j < size-2; j++) {
      sweep(xpx, size, index+j, 1, &err);
      a = xpx[sizesqm1];
      sweep(xpx, size, index+j, 1, &err);

      if(a < rss[i+1]) {
        k = index[j]; 
        index[j] = index[i];
        index[i] = k;
        rss[i+1] = a;
      }
    }
    sweep(xpx, size, index+i, 1, &err);
  }
}





/**********************************************************************
 * 
 * identify_qtls:  take in lod curve, threshold, and min. drop between
 *                 peaks and output the chromosome and marker numbers
 *                 for the chosen QTLs
 *
 * input:
 *
 *    n_chr, n_mar               = info about genome size
 * 
 *    lod                        = lod curve
 *
 *    threshold                  = threshold to use
 *
 *    drop                       = value lod curve must drop in-between
 *                                 inferred peaks
 *
 *    n_qtl_id                   = on output, number of QTLs identified
 *
 *    chr_id                     = vector of length at least tot_mar
 *                                 on output it contains the chromosome
 *                                 numbers of the identified QTLs
 *
 *    mar_id                     = vector of length at least tot_mar
 *                                 on output it contains the marker
 *                                 numbers of the identified QTLs
 *
 *    lwork                      = workspace of ints of size tot_mar
 *
 *    dwork                      = workspace of doubles of size tot_mar
 *
 **********************************************************************/

void identify_qtls(int n_chr, int *n_mar, double *lod,
                   double threshold, double drop, int *n_qtl_id,
                   int *chr_id, int *mar_id, int *lwork, 
                   double *dwork)
{
  int i, j, k, n_found, flag, flag1;
  int closest, first_on_chr, cur;

  *n_qtl_id = 0;
  for(i=0, first_on_chr=0; i< n_chr; first_on_chr += n_mar[i], i++) {
    /* any above threshold ? */
    flag1=0;
    for(j=0, cur=first_on_chr; j < n_mar[i]; j++, cur++) {
      if(lod[cur] > threshold) {
	flag1=1;
	break;
      }
    }

    if(flag1) { /* at least one above threshold */
      n_found = 0;
      /* sort lods and indices */
      for(j=0, cur=first_on_chr; j< n_mar[i]; j++, cur++) {
	lwork[j] = j+1;
	dwork[j] = lod[cur];
      }
      piksrt(n_mar[i], dwork, lwork);
    
      for(j=0, cur=first_on_chr; j< n_mar[i] && dwork[j] > threshold; 
	  j++, cur++) {
	if(!n_found) {
	  /* first QTL identified on this chromosome */
	  chr_id[*n_qtl_id] = i+1;
	  mar_id[*n_qtl_id] = lwork[j];
	  (*n_qtl_id)++;
	  n_found++;
	}
	else {
	  closest = mar_id[*n_qtl_id - n_found];
	  if(n_found > 1) {
	    /* find closest identified QTL */
	    for(k= *n_qtl_id - n_found + 1; k < *n_qtl_id; k++) 
	      if(abs(closest - lwork[j]) > abs(mar_id[k] - lwork[j]))
		closest = mar_id[k];
	  }
	  flag = 0;
	  /* closest identified comes before the current one */
	  if(closest < lwork[j]) {
	    for(k=closest+1; k<lwork[j]; k++) 
	      if(lod[first_on_chr + k-1] < dwork[j] - drop) {
		flag = 1;
		break;
	      }
	  }
	  else { 
	    for(k=lwork[j]+1; k<closest; k++) 
	      if(lod[first_on_chr + k-1] < dwork[j] - drop) {
		flag = 1;
		break;
	      }
	  }
	  if(flag) {
	    chr_id[*n_qtl_id] = i+1;
	    mar_id[*n_qtl_id] = lwork[j];
	    (*n_qtl_id)++;
	    n_found++;
	  }
	}
      } 
    } /* at least one identified QTL on chromosome */
  } /* loop over chromosomes */
}

/*
void R_identify_qtls(int *n_chr, int *n_mar, double *lod,
		     double *threshold, double *drop, int *n_qtl_id,
		     int *chr_id, int *mar_id, int *lwork, 
		     double *dwork)
{
  identify_qtls(*n_chr, n_mar, lod, *threshold, *drop, n_qtl_id,
		chr_id, mar_id, lwork, dwork);
}
*/

/**********************************************************************
 * 
 * piksrt
 *
 * This function was take from Numerical Recipes in C (Sec 8.1:
 * straight insertion and Shell's method).  It should suffice for my 
 * very simple sorting needs.
 *
 * Input:
 *
 *     n   = length of vector to sort
 *   
 *   arr   = pointer to array to be sorted into ascending order;
 *           on output, it contains the sorted array
 * 
 *   larr  = pointer to array to be sorted aint side arr 
 *
 **********************************************************************/

void piksrt(int n, double *arr, int *larr)
{
  int i,j, b;
  double a;
  
  for(j=1; j<n; j++) {        /* pick out each element in turn */
    a=arr[j];
    b = larr[j];
    i=j-1;
    while(i>=0 && arr[i]<a) {   /* look for the place to insert it */
      arr[i+1]=arr[i];
      larr[i+1] = larr[i];
      i--;
    }
    arr[i+1]=a;                /* insert it */ 
    larr[i+1]=b;
  }
}



/**********************************************************************
 * 
 * identify_qtls_bic:  take the results of forward or backward selection
 *                     and use BIC to choose a model and pick out the 
 *                     location of QTLs
 *
 * input:
 *
 *    n_chr, n_mar               = info about genome size
 * 
 *    n_progeny
 *
 *    index                      = index giving order of markers to
 *                                 be placed in model
 *
 *    rss                        = corresponding residual sum of squares
 *                                 for the set of markers
 * 
 *    n_models                   = number of models generated (used when
 *                                 forward selection is not carried out
 *                                 completely)
 *
 *    n_qtl_id                   = on output, number of QTLs identified
 *
 *    chr_id                     = vector of length at least n_chr * n_mar,
 *                                 on output it contains the chromosome
 *                                 numbers of the identified QTLs
 *
 *    mar_id                     = vector of length at least n_chr * n_mar,
 *                                 on output it contains the marker
 *                                 numbers of the identified QTLs
 *
 *    multiplier                 = 1.0 to use BIC, 2.0 to use 2*BIC, etc
 *
 **********************************************************************/

void identify_qtls_bic(int n_chr, int *n_mar, int n_progeny,
                       int *index, double *rss, int n_models, 
                       int *n_qtl_id, int *chr_id, int *mar_id,
                       double multiplier)
{
  int i, j, first_on_chr;
  double a, b;

  /* BIC for null model */
  a = log(rss[0]); 
  *n_qtl_id = 0;

  /* calculate BIC for the rest of the models,
     and determine number of QTLs in best model */
  for(i=1; i<=n_models; i++) {
    b = log(rss[i]) + multiplier*(double)i*log((double)n_progeny)/
      (double)n_progeny;
    
    if(b < a) {
      *n_qtl_id = i;
      a = b;
    }
  }
  
  /* fill up chr_id and mar_id with chosen QTLs */
  if(*n_qtl_id) {
    for(i=0; i < *n_qtl_id; i++) {
      for(j=0, first_on_chr=1; j<n_chr; first_on_chr += n_mar[j], j++) {
	if(index[i] < first_on_chr + n_mar[j]) { 
	  chr_id[i] = j+1;
	  mar_id[i] = index[i] - first_on_chr + 1; 
	  break;
	}
      }
    }
  }
}



/*
void R_identify_qtls_bic(int *n_chr, int *n_mar, int *n_progeny,
			 int *index, double *rss, int *n_models, 
			 int *n_qtl_id, int *chr_id, int *mar_id,
			 double *multiplier)
{
     identify_qtls_bic(*n_chr, n_mar, *n_progeny, index, rss, 
		       *n_models, n_qtl_id, chr_id, mar_id,
		       *multiplier);
}
*/


/**********************************************************************
 * 
 * which_correct: figure out which QTLs were correctly identified
 *                and how many were incorrectly identified
 *
 * n_infer = inferred number of QTLs
 * chr_infer = chromosome number for each inferred QTL (length n_infer)
 * mar_infer = marker number for each inferred QTL (length n_infer)
 * n_true = true number of QTLs
 * chr_true = chromosome number for each true QTL (length n_true)
 * mar_true = marker number for each true QTL (length n_true)
 * within = how many markers away can an inferred QTL be froma true QTL?
 * correct = on output, vector of length n_true with 0/1 indicating
 *           whether a qtl was correctly identified
 * n_incorrect = on output, vector of length two giving the number of
 *               extraneous loci on the same chromosome as a QTL and
 *               the total number of the number of extraneous loci
 *
 **********************************************************************/

void which_correct(int n_infer, int *chr_infer, int *mar_infer,
		   int n_true, int *chr_true, int *mar_true,
		   int within, int *correct, int *n_incorrect)
{
  int i, j, flag;
  
  for(i=0; i<n_true;  i++) correct[i] = 0;
  n_incorrect[0] = n_incorrect[1] = 0;
  for(i=0; i<n_infer; i++) {
    flag = 0;
    for(j=0; j<n_true; j++) {
      /* correctly identified? */
      if(correct[j] != 1 && chr_infer[i] == chr_true[j] &&
	 mar_infer[i] >= mar_true[j] - within &&
	 mar_infer[i] <= mar_true[j] + within) {
	flag = 1;
	correct[j] = 1;
	break;
      }
    }
    if(flag==0) {
      for(j=0; j<n_true; j++) {
	if(chr_infer[i]==chr_true[j]) { 
	  flag = 1;
	  break;
	}
      }
      if(flag==1) (n_incorrect[0])++;
      else (n_incorrect[1])++;
    }
  }

}

void R_which_correct(int *n_infer, int *chr_infer, int *mar_infer,
		     int *n_true, int *chr_true, int *mar_true,
		     int *within, int *correct, int *n_incorrect)
{     
  which_correct(*n_infer, chr_infer, mar_infer, *n_true, 
		chr_true, mar_true, *within, correct, n_incorrect);
}

/* end of analbc.c */
