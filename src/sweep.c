/**********************************************************************
 *
 *
 *                        sweep.c
 *
 *                  Karl Broman (s243ap)
 *                     Assignment #4
 *                    3 December 1993
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>

#include "sweep.h"


/**********************************************************************
 *
 *      This function "sweeps" the columns of the matrix
 *  "matrix" which are listed in the vector "index".
 *  <size> is the # of rows in <matrix>.  <len> is the
 *  length of <index>.
 *
 *      The indices are in 0, 1, 2, ..., size - 1.
 *
 **********************************************************************/

void sweep(double *matrix, int size, int *index, int len, int *err)
{
  int i, j, k;
  double b, d;

  *err=0;
  for(k=0; k<len; k++) {
    d=matrix[index[k]+index[k]*size];
    if(fabs(d)<=ZERO) {
      warning("Problem in sweep: The pivot is very near zero:\n\t%.12lf\n",d);
      *err=1;
      return;
    }
    for(i=0; i<size; i++) matrix[index[k]+i*size] /= d;
    for(i=0; i<size; i++) {
      if(i!=index[k]) {
        b=matrix[i+index[k]*size];
        for(j=0; j<size; j++)
          matrix[i+j*size] -= (b*matrix[index[k]+j*size]);
        matrix[i+index[k]*size] = -b/d;
      }
    }
    matrix[index[k]+index[k]*size] = 1./d;
  }
}



/**********************************************************************
 *
 *     This is a little wrapper for the above function
 *  sweep() which makes it easier to call from S.
 *
 **********************************************************************/

void s_sweep(double *a, int *n, int *i, int *l, int *err)
{
  sweep(a, *n, i, *l, err);
}


/*            This is the end of <sweep.c>                 */
