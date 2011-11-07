/**********************************************************************
 * 
 * simbc.h:   Programs to simulate and analyze QTL data from a backcross,
 *            with unequal numbers of markers per chromosome, and at 
 *            unequal spacing
 *
 *   
 *
 * Karl Broman, 7/6/01 [originally 9/30/96, 10/7/96, 10/8/96]
 *
 **********************************************************************/

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
	       double *recfrac, int *genotypes);
  
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
	       double *effect, double sigma);
  
void R_simbc_qtl(int *n_progeny, int *n_chromosomes, int *n_markers,
		 double *recfrac, int *genotypes, 
		 double *phenotypes, int *n_qtl, int *qtl_chr,
		 int *mar_to_left, double *recfrac_to_left, 
		 double *effect, double *sigma);

/* end of simbc.h */
