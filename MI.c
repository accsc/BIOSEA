/*********************************************************************
 2014, 2015. Alvaro Cortes Cabrera 

 ---------------------------------------------------------------------
 This program reads two sets with Biological fingerprints in TSV format
 and applies the BioSEA method with the input parameters.

 It outputs the results per molecule in the query set

 Input TSV file:

 First line: Title. Names of the assays (short codes) separated by tabs
 Next lines: Values. Each one separated by a tab. If the value is 
             missing then a tab follows inmediatly to the other. 
	     No value.

 0.1 - Initial version. 

 ---------------------------------------------------------------------

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.


********************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


/* Show how to use the program. Name is argv[0] so I can print 
   the name of the acutal program on the screen 
*/
void show_use(char *name)
{
	fprintf(stderr,"%s FP_file1 FP_file2 cutoff a1 b1 a2 b2\n\n",name);
	fprintf(stderr,"FP_file1. TSV file with the fingerprints of the query molecule/s\n");
	fprintf(stderr,"FP_file2. TSV file with the fingerprints of the reference molecules\n");
	fprintf(stderr,"Cutoff. Threshold to consider the pairwise similarity. Depends on fingerprint\n");
	fprintf(stderr,"a1. A parameter of the power-law to calculate the mean score given a number of molecules. aN**b\n");
	fprintf(stderr,"b1. B parameter of the power-law to calculate the mean score given a number of molecules. aN**b\n");
        fprintf(stderr,"a2. A parameter of the power-law to calculate the std deviation score given a number of molecules. aN**b\n");
        fprintf(stderr,"b2. B parameter of the power-law to calculate the std deviation score given a number of molecules. aN**b\n\n");
	fflush(stderr);
}

/**
 * Translated version of the digamma (psi) function
 * from a old Fortran 77 version. 
 *
 * X is the number. Returns the value.
 */
double psi(double x){
	double xa, pi, el, s, x2, ps, a1, a2, a3, a4, a5, a6, a7, a8;
	int n = 0;
	int i = 0;
	xa = fabs(x);
        pi=3.141592653589793;
        el=0.5772156649015329;
        s=0.0;

	if( x == (int)x && x <= 0)
	 return 300;

	if( xa == (int) xa)
	{
		n = (int)xa;
		s = 0.0;

		for( i = 1; i < n; i++)
		{
			s = s + 1.0/(double)i;
		}
		ps = -el+s;

	}else if( xa+0.5 == ((int)xa)+0.5){
		n = xa - 0.5;
                for( i = 1; i <= (int) n; i++)
                        s += 1.0/((2.0*(double) i )-1.0);
		ps = -el+2.0*s-1.386294361119891;

	}else{
	 	if( xa < 10.0)
		{
			n = 10 - (int) xa;

	                for( i = 1; i < (int) n; i++)
			{
				s += 1.0 / (xa+(double)i);
			}
			xa += (double) n;

		}
	
           x2=1.0/(xa*xa);
           a1=-.8333333333333E-01;
           a2=.83333333333333333E-02;
           a3=-.39682539682539683E-02;
           a4=.41666666666666667E-02;
           a5=-.75757575757575758E-02;
           a6=.21092796092796093E-01;
           a7=-.83333333333333333E-01;
           a8=.4432598039215686;

           ps=log(xa)-.5/xa+x2*(((((((a8*x2+a7)*x2+a6)*x2+a5)*x2+a4)*x2+a3)*x2+a2)*x2+a1);
           ps=ps-s;
	}
	if ( x < 0.0)
        ps=ps-pi*cos(pi*x)/sin(pi*x)-1.0/x;

	return ps;
}

/** 
 * Digamma function by Tailor expansion. Much simpler. 
 * However I'm using the old one for compatibility reasons
 * x is the number. Returns the value.
 **/
double psi_simple(double x) {
  double result = 0, xx, xx2, xx4;
  for ( ; x < 7; ++x)
    result -= 1/x;
  x -= 1.0/2.0;
  xx = 1.0/x;
  xx2 = xx*xx;
  xx4 = xx2*xx2;
  result += log(x)+(1./24.)*xx2-(7.0/960.0)*xx4+(31.0/8064.0)*xx4*xx2-(127.0/30720.0)*xx4*xx4;
  return result;
}


/**
 * Computes Mutual Information based on kNN for two fingerprints (vectors)
 * using the Kraskov, Stoegbauer and Grassberger method.
 *
 * fp1. First fingerprint of length len
 * fp2. Second fingerprint of length len
 * norm_flag. If 1, normalize the MI, if else, nothing (to get the entropy). 
 *            MI(X;X) = H(X)
 * incommon. Returns the number of shared results
 * len. Length of the fingerprints
 * k_param. Use k neighbors for MI estimation
 *
 * returns the actual MI value
 */
double MI_knn( double * fp1, double *fp2, int norm_flag, int *incommon, int len, int k_param)
{
	int i = 0, j = 0, n = 0, k = 0, swpi = 0, sum_x = 0, sum_y = 0;
	double res = 0.0, swp = 0.0, maxe = 0.0f, dig = 0.0f, digav = 0.0f;
	double *x = NULL, *y = NULL, *r = NULL, *dx = NULL, *dy = NULL, *rdup = NULL;
	int *rank = NULL, *order = NULL;
	int change = 0;
	double cdx = 0.0f, cdy = 0.0f, H_X = 1.0f, H_Y = 1.0f;
	int in2 = 0;
	int sum_nx = 0, sum_ny = 0;

	/* Just the parts that we have in common */
	x = (double*) calloc(sizeof(double),len);
	y = (double*) calloc(sizeof(double),len);

	for( i = 0; i < len; i++)
	{

		/* 999.9 and beyond means that the value is not available */
                if( (fp1[i] < 999.9) && (fp2[i] < 999.9))
		{
			x[j] = fp1[i];
			y[j] = fp2[i];
			j++;
		}
	}

	/* Let them know how many assays we have */
	n = j;
	*incommon = n;

	if( n < (k_param+2))
	{
		free(x); free(y);
		return 0.0f;
	}


	/* Memory allocation */
	r = (double*) calloc(sizeof(double),n);
	rdup = (double*) calloc(sizeof(double),n);
	rank = (int*) calloc(sizeof(int),n);
	order = (int*) calloc(sizeof(int),n);
	dx = (double*) calloc(sizeof(double),n);
	dy = (double*) calloc(sizeof(double),n);

	for( i = 0; i < n; i++)
	{
		for( j = 0; j < n; j++)
		{
			cdx = fabs(x[i]-x[j]);
			cdy = fabs(y[i]-y[j]);
			if (cdx > cdy)
				r[j] = cdx;
			else
	                        r[j] = cdy;

			rdup[j] = r[j];
			order[j] = j;
			dx[j] = cdx;
			dy[j] = cdy;
		}

		change = 1;
		/* Bubble sort (face palm). No comments */
		while( change != 0) 
		{
			change = 0;
			for( k = 0; k < n-1; k++)
			{
				if( rdup[k] > rdup[k+1])
				{
				 change = 1;
                                 swp = rdup[k+1];
                                 rdup[k+1] = rdup[k];
                                 rdup[k] = swp;
                                 swpi = order[k+1];
				 order[k+1] = order[k];
				 order[k] = swpi;
				}
			}
		}

		/* Get rank */
		for( j = 0; j < n; j++)
                 rank[order[j]] = j + 1;


	        maxe = -9999;
	        for ( k = 0; k < n; k++)
		{
	          if( rank[k] <= (k_param+1))
		  {
	           if( r[k] > maxe)
	            maxe = r[k];
                  }
	        }

		sum_nx = 0;
		sum_ny = 0;
		for( k = 0; k < n; k++)
		{
			if( dx[k] < maxe)
			  sum_nx++;
			if( dy[k] < maxe)
                          sum_ny++;
		}
		dig = psi(sum_nx)+psi(sum_ny);
		digav += dig;
	}
	digav /= (double) n;
        res = psi(n) + psi((double)k_param) - digav;

	/* Normalization here means we get H(X) and H(Y) by calling again 
  	   this function but this time with norm_flag off and the same fp
           each time so H(fp1) = MI(fp1;fp1) and H(fp2) = MI(fp2;fp2)
           Then normalization means we divide by SQRT(Hx*Hy)
        */
	if ( norm_flag == 1)
	{
		H_X = MI_knn(fp1,fp1,0, &in2, len, k_param);
		H_Y = MI_knn(fp2,fp2,0, &in2, len, k_param);
		res /= sqrt(H_X*H_Y);
	}

	/* free all those horrible mallocs */
        free(x); free(y); free(r); free(rank); free(order);
	free(dx); free(dy); free(rdup);

	return res;
}


int main(int argc, char *argv[])
{

	FILE *input = NULL;
	char *buffer_line = NULL;
	char *point = NULL, *token = NULL;
	double **matrix = NULL, **matrix2 = NULL; 
        char **mol_names = NULL, **mol_names2 = NULL;
	int lines = 0, i = 0, j = 0, k = 0, z = 0, lines2 = 0, a = 0, b = 0;
	size_t l = 0;
	int n_column = 0, n_column2 = 0;
	char number[50];
	char *buffer2 = NULL;
	int incommon = 0;
	double r = 0.0f, tc = 0.0;

	int exact = 0, total = 0, real_total = 0;
	double sum = 0.0f;
	double cutoff = 0.0, sigma = 1.0;
	double zscore = 0.0f, av = 0.0f, sd = 0.0f, prob = 0.0f;
	double shape = 0.0f, scale = 0.0, location = 0.0f;

	float powerlaw1_a = 0.0f, powerlaw1_b = 0.0f, powerlaw2_a = 0.0f, powerlaw2_b = 0.0f;

	if( argc < 8)
	{
		fprintf(stderr,"Error. Missing arguments\n");
		fflush(stderr);
		show_use(argv[0]);
		return -1;
	}

        /* The actual zscore threshold for taking a contribution into account */
        cutoff = atof(argv[3]);

	/* Power law parameters for average and standard deviation in SEA */
	powerlaw1_a = atof(argv[4]);
	powerlaw1_b = atof(argv[5]);
	powerlaw2_a = atof(argv[6]);
	powerlaw2_b = atof(argv[7]);

	/* Check parameters */
	fprintf(stderr,"Parameters: Cutoff: %8.3f, Average: %E*N(**%E), Std Dev: %E*N(**%E)\n",cutoff,
		powerlaw1_a,powerlaw1_b,powerlaw2_a,powerlaw2_b);
	fflush(stderr);


	if( (input = fopen(argv[1],"r")) == NULL)
	{
		fprintf(stderr,"Cannot open the query set.\n");
		fflush(stderr);
		return -1;
	}

	/* Buffer to load lines
           Keep in mind that we cannot read lines of more than 8kb 
           Modify is FP is huge or something
        */
	buffer_line = (char *) calloc(sizeof(char),8096);

	/* Count number of lines = molecules+1 */
	while( (fgets(buffer_line,8095,input)) != NULL)
	{
		++lines;
	}

	rewind(input);
	fgets(buffer_line,8095,input); /* Title */
	
	/* strtok does not work with empty fields */
	token = strtok_r(buffer_line,"\t",&point);
	while( token != NULL)
	{
		++n_column;
		token = strtok_r(NULL,"\t",&point);
	}
	/* We expect the file to have a title line **with**
           all the columns. We use this number 
           as the final number of columns.
        */

	rewind(input);
        fgets(buffer_line,8095,input); /* Title */

	fprintf(stderr,"Molecules on query set: %i\n",lines-1);
	fprintf(stderr,"Columns on query set: %i\n",n_column);

	/* For the actual values in the fingerprint */
	matrix = (double **)calloc( sizeof( double *), lines+1);
	/* Name of the molecules. Should be the first column */
        mol_names = (char **) calloc(sizeof(char *), lines);

	/* Allocate the thing */
	for( i = 0; i < (lines-1); i++)
	{
		matrix[i] = (double *) calloc(sizeof(double), n_column);
                mol_names[i] = (char *) calloc(sizeof(char), 100); /* More than enough */
	}
	
	/* Read matrix */
	i = 0;
        while( (fgets(buffer_line,8095,input)) != NULL)
	{
		/* strtok_r is not working with missing data! */
		/* This works! */
		buffer2 = buffer_line;

		k = strcspn(buffer_line,"\t");
		if( k > 99) /* Cut molecular name */
                {
			strncpy(mol_names[i],buffer2,99);
		}else{
	                strncpy(mol_names[i],buffer2,k);
		}

		/* Move the pointer to the next field */
		for( z = 0; z <= k; z++) 
			buffer2++;

		/* Read the actual values */
		for( j = 0; j < n_column; j++)
		{
			k = strcspn(buffer2,"\t\n");
			if( k > 48)
                        {
				fprintf(stderr,"Huge number. Something is wrong. Please check input format\n");
				fflush(stderr);
		                for( i = 0; i < (lines-1); i++)
		                {
		                        free(matrix[i]);
                                        free(mol_names[i]);
		                }
		
		                free(matrix);
		                free(mol_names);
		                free(buffer_line);
				return -1;
			}
			strncpy(number,buffer2,k);
			if (k > 1)
			{
			 number[k] = '\0';
			 matrix[i][j] = atof(number);
	                 for( z = 0; z <= k; z++)
	                         buffer2++;
			}else{ /* There is a gap here. We fill it with 999.9f. So be advised. */
			 matrix[i][j] = 999.9f;
			 buffer2++;
			}
		}
		i++;
	}

	/* If it is huge you want to know what this is doing */
	fprintf(stderr,"Query set in memory\n");
	fflush(stderr);

	/* We are done with the first set */
	fclose(input);


	if( (input = fopen(argv[2],"r")) == NULL)
	{
		fprintf(stderr,"Cannot open the reference set.\n");
		fflush(stderr);
	        for( i = 0; i < (lines-1); i++)
	        {
	                free(matrix[i]);
                        free(mol_names[i]);
	        }

	        free(matrix);
	        free(mol_names);
	        free(buffer_line);
		return -1;
	}

	/* Count the fingerprints in the second set */
	while( (fgets(buffer_line,8095,input)) != NULL)
	{
		++lines2;
	}

	rewind(input);
	fgets(buffer_line,8095,input); /* Title */

        /* strtok does not work with empty fields */
        token = strtok_r(buffer_line,"\t",&point);
        while( token != NULL)
        {
                ++n_column2;
                token = strtok_r(NULL,"\t",&point);
        }

	if( n_column != n_column2)
	{
		fprintf(stderr,"Both sets have a different number of columns (%i and %i). Please check your input file.\n",n_column,n_column2);
		fflush(stderr);
                for( i = 0; i < (lines-1); i++)
                {
                        free(matrix[i]);
                        free(mol_names[i]);
                }

                free(matrix);
                free(mol_names);
                free(buffer_line);
		return -1;
	}

	if( (matrix2 = (double **)calloc( sizeof( double *), lines2+1)) == NULL)
	{
		fprintf(stderr,"Error allocating memory. No matrix for molecules: %i\n",lines2-1);
		fflush(stderr);
                for( i = 0; i < (lines-1); i++)
                {
                        free(matrix[i]);
                        free(mol_names[i]);
                }

                free(matrix);
                free(mol_names);
                free(buffer_line);
		return -1;
	}

	if( (mol_names2 = (char **) calloc(sizeof(char *), lines2)) == NULL)
	{
		fprintf(stderr,"Error allocating memory.\n");
		fflush(stderr);
                for( i = 0; i < (lines-1); i++)
                {
                        free(matrix[i]);
                        free(mol_names[i]);
                }

                free(matrix);
                free(mol_names);
                free(matrix2);
                free(buffer_line);
		return -1;
	}

	for( i = 0; i < (lines2-1); i++)
	{
		matrix2[i] = (double *) calloc(sizeof(double), n_column);
		mol_names2[i] = (char *) calloc(sizeof(char), 100); /* More than enough */
	}

	i = 0;
	while( (fgets(buffer_line,8095,input)) != NULL)
	{
		/* strtok_r is not working with missing data! */
		/* This works! */
		buffer2 = buffer_line;
		k = strcspn(buffer_line,"\t");
                if( k > 99) /* Cut molecular name */
                {
                        strncpy(mol_names2[i],buffer2,98);
                }else{
                        strncpy(mol_names2[i],buffer2,k);
                }

		for( z = 0; z <= k; z++)
			buffer2++;

		for( j = 0; j < n_column; j++)
		{
			k = strcspn(buffer2,"\t\n");
                        if( k > 48)
                        {
                                fprintf(stderr,"Huge number. Something is wrong. Please check input format\n");
				fflush(stderr);
                                for( i = 0; i < (lines-1); i++)
                                {
                                        free(matrix[i]);
					free(mol_names[i]);
                                }
                                free(matrix);
				free(mol_names);
                                for( i = 0; i < (lines2-1); i++)
                                {
                                        free(matrix2[i]);
                                        free(mol_names2[i]);
                                }
                                free(mol_names2);
                                free(matrix2);
                                free(buffer_line);
                                return -1;
                        }
			strncpy(number,buffer2,k);
			if (k > 1)
			{
			 number[k] = '\0';
			 matrix2[i][j] = atof(number);
			 for( z = 0; z <= k; z++)
				 buffer2++;
			}else{ /* Again. 999.9 means no value */
			 matrix2[i][j] = 999.9f;
			 buffer2++;
			}
		}
		i++;
	}

	/* This should be smaller */
	fprintf(stderr,"Reference set in memory\n");
	fflush(stderr);

	fclose(input);

	/* Easier to remember */
        total = lines2-1;

	/* Print info */
        fprintf(stdout,"CID ZSCORE N EXACT SUM N_REAL\n");
	fflush(stdout);

        for( i = 0; i< lines-1; i++)
        {
		exact = 0; /* Number of times the molecules is in the reference set */
		sum = 0.0; /* Raw sum of scores */
		real_total = 0; /* Number of possible comparisons (enough assays) */
	/* Run in parallel for each query molecule */
        #pragma omp parallel for \
                    private(j,incommon,tc,sigma,av,sd,zscore,prob) \
                    shared(matrix, matrix2, exact, i,n_column, total) \
                    reduction(+ : sum, real_total) \
                    schedule(static)
                for( j = 0; j< lines2-1; j++)
                {
				incommon = 0;

				/* Compare molecular names. Skip if same */
                                if( strcmp(mol_names[i],mol_names2[j]) == 0){
					exact += 1;
					continue;
				}

				/* Calculate normalized MI for the two fingerprints */
				tc = MI_knn(matrix[i],matrix2[j],1,&incommon,n_column-1,10);
			
				/* Minimum number of assays in common for considering the info */
				if( incommon <= 12)
				 continue;

				/* Normalize again based on the number of assays for the calculation */
				/* These numbers are based on random comparions using the HTSfps from PubChem */
                                sigma = 5.130027E-01 * powf( (float)incommon, -0.787372);

				real_total += 1;

				/* Higher than threshold? */
				if( tc/sigma > cutoff)
				{
					sum += tc/sigma;
                                        #ifdef DEBUG
                                           fprintf(stdout,"%s %s %f %f %i\n",mol_names[i],mol_names2[j],tc/sigma,tc,incommon);
                                           fflush(stdout);
                                        #endif
				}
                }


		/* Calculate the average expected at random*/
		av = powerlaw1_a*pow(total,powerlaw1_b);
                /* Calculate the standard deviation expected at random*/
                sd = powerlaw2_a*pow(total,powerlaw2_b);

                zscore = (sum - av)/sd;

		/* Print results */
                fprintf(stdout,"%s %f %i %i %f %i\n" ,mol_names[i],zscore,total,exact,sum,real_total);
		fflush(stdout); /* Slower but you have the complete information till it fails */
        }

	/* Clean up */
        for( i = 0; i < (lines-1); i++)
        {
                free(matrix[i]);
		free(mol_names[i]);
        }

        for( i = 0; i < (lines2-1); i++)
        {
                free(matrix2[i]);
                free(mol_names2[i]);
        }
	free(matrix);
	free(mol_names);
	free(mol_names2);
	free(matrix2);
	free(buffer_line);
	exit(0);
}
