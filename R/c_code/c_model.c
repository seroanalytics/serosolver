 /* Set up likelihood for vector of infection history */

#include <math.h>
#include <stdio.h>

/* Inputs needed:
	- test strain data vector
	- infection history vector
	- d.ij vector
 */

void c_model2_sr(int *nin, int *itot, int *nsin, double *x, double *x1, double *titre, 
                  double *titrepred, double *dd, double *dd2, int *ntheta, 
                  double *theta, int *inputtestyr)
{
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/* Calculate lambda */
	
	int n = nin[0]; // length of infection history
	// int total_inf = itot[0]; // total infections
	int nsamp = nsin[0]; // number of strains tested against
	// double T_1 = theta[1];
	double T_2 = theta[2];
	double wane = theta[3];
	double mu = theta[0];
	double mu2 = theta[5]; // as sigma = theta[4]
	
	// This to be made an argument of the function -- gives test year 
	int t_sample = inputtestyr[0]; 
  	double yrTitre[nsamp];
  	int maskedInfectionHistory[n];
  	double distanceFromTest[n];
  	double cumInfectionHistory[n];
	
	#define MAX(a,b) ((a) < (b) ? (b) : (a)) // define MAX function for use later
  
	/* Add for loop over k*/
	
	int k;
	int i;
	int j;
	int m;
	
	double xx2; 
	
	for (k=0; k<nsamp; k++){ // Iterate over samples tested against

	  xx2=0;

	  // Make a masked infection history
	  j=(t_sample-1); // fix test year. Need (t-1) as index from 0
	  for (m=0;m<n;m++) {
	    if (m <= j) {
	      maskedInfectionHistory[m] = x[m];
	    } else {
	      maskedInfectionHistory[m] = 0;
	    }
	  }
	  
	  // Make an index for waning
	  for (m=0;m<n;m++) {
		  // distanceFromTest[m] = exp(-wane * (j-m )); // Distance from test year 
		  distanceFromTest[m] = MAX(0, 1 - wane * (j-m) ); // Distance from test year
	  }

	  // Make a cumulative infection history
	  cumInfectionHistory[0] = maskedInfectionHistory[0];
	  for (m=1;m<n;m++) {
	    cumInfectionHistory[m] = cumInfectionHistory[m-1] + 
	      maskedInfectionHistory[m];
	  }
	  	    
		/* Calculate expected titre	- note k indexed from 0 */
	    /* Note that waning is currently linked with back boosting */
	    /* dd is long term cross-reaction, dd2 is short-term */

		for (i=0; i<n; i++){
			x1[i] = maskedInfectionHistory[i] *
			  // exp(-1.0 * T_2 * ( cumInfectionHistory[i]  - 1.0)) *
			  MAX(0, 1.0 - T_2 * ( cumInfectionHistory[i]  - 1.0)) * // Antigenic seniority
			  //(pow(1.0 + T_1 , (total_inf - cumInfectionHistory[i])) ) * REMOVED Tau 1
			  (mu * dd[k*n+i] + mu2 * dd2[k*n+i] * distanceFromTest[i] );
		}
	
		for (i=0; i<n; i++){
			xx2 =  xx2 + x1[i];
		}
	
	  yrTitre[k] = xx2;
	
	
	} // end sample loop (k)
	
	for (k=0;k<nsamp;k++) {
	  titrepred[k] = yrTitre[k]; 
	}
	
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/* Compare to observed titre*/

}