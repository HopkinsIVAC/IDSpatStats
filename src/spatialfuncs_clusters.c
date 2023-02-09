#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>


/*****************************************************************/
/*pi function optimized for iterating through a list of two types*/
/*finding ones that fulfill the distance requirement.            */
/*                                                               */
/*@param type_1  the  first type column                          */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_pi_typed_wts (int *type,
		   double *x,
		   double *y,
		   double *weight,
		   int *len,
		   int *typeA,
		   int *typeB,
		   double *r_low,
		   double *r,
		   int *len_r,
		   int *inds,
		   double *rc) {

  int i,j,k;
  int num_cnt, denom_cnt; /*counters for those filling conditions*/
  double dist;

  /*repeat the calculation for all r*/
  for (i=0;i<*len_r;i++) {

    /*zero out counts*/
    num_cnt = 0;
    denom_cnt = 0;

    if (*typeA != -1) {

      for (j=0;j<*len;j++) {
	if (type[j] != *typeA) continue;

	for (k=0;k<*len;k++) {
	  /*ignore pairs of the same individual*/
	  if (inds[k]==inds[j]) continue;

	  dist = sqrt(pow(x[j]-x[k],2)+pow(y[j]-y[k],2));
	  if ((dist<=r[i])  & (dist>=r_low[i])) denom_cnt = denom_cnt + weight[k]*weight[j];

	  if (type[k] != *typeB) continue;
	  if ((dist<=r[i])  & (dist>=r_low[i])) num_cnt = num_cnt + weight[k]*weight[j];
	}
      }

    } else {
      Rprintf("To be implemented\n");
      return;
    }
    //Rprintf("%d/%d\n",num_cnt,denom_cnt);//DEBUG
    rc[i] = (double)num_cnt/denom_cnt;
  }
}



/*****************************************************************/
/*pi function optimized for iterating through a list of two types*/
/*finding ones that fulfill the distance requirement.            */
/*                                                               */
/*@param type_1  the  first type column                          */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_pi_clustsurvey (double *p,
		   double *x,
		   double *y,
		   int *s,
		   int *len,
		   double *r_low,
		   double *r,
		   int *len_r,
		   int *inds,
		   double *rc) {

  int i,j,k;
  int skip_r, skip_phi;  /*counters for those filling conditions*/
  double dist;
  double si, sj, sum_sj;
  double num_typeAj, num_typeAi, sum_typeAi, sum_typeAj, phi, sum_phi;
  double mean_probsuscept;

  	//repeat the calculation for all r
	for (k=0;k<*len_r;k++) {
			
		sum_phi = 0;
		skip_r = 1;
		sum_typeAi = 0;
		si = 0;
		
		for (i=0;i<*len;i++) {
			
			// Skip cluster if no individuals with the characteristic
			if (p[i] == 0) continue; 
			skip_phi = 1;
			phi = 0;
			sum_sj = 0;		
			num_typeAi = 0;
			sum_typeAj = 0;			
			
			
			// Cycle through paired clusters with cluster j
			//  - need to add s and r to each cluster within the correct distance to the sums
			//  - need to sum the sum_typeAj; 
			for (j=0;j<*len;j++) {
				num_typeAj = 0;
				sj = 0;
                
                // Calculate distance between clusters i and j
				dist = sqrt(pow(x[i]-x[j],2)+pow(y[i]-y[j],2));
				  
				if ((dist<=r[k]) & (dist>=r_low[k])) {
					sj = s[j];  // number of individuals in the cluster j
					num_typeAj = p[j]*sj;  // number of individuals of type A in cluster i (type A = susceptible)
					
					// Subtract 1 from total number in the cluster and 1 from number susceptible if looking at contact within cluster.
					// Cannot infect self
					if (inds[i]==inds[j]) {
						sj = s[j] - 1;
						num_typeAj = num_typeAj - 1;
					}
					
					sum_typeAj = sum_typeAj + num_typeAj; // Sum of individuals in clusters j at distance r[k] from cluster i
					sum_sj = sum_sj + sj;  // sum of s[j] at distance r[k] from cluster i
					skip_phi = 0; // don't skip phi
				}
			}
			
			if (skip_phi == 1) continue;
			
			si = s[i];  // number of individuals in the cluster i
			num_typeAi = p[i]*si;  // number of individuals of type A in cluster i (type A = susceptibles)
			sum_typeAi = sum_typeAi + num_typeAi;	  // sum of number of type A in primary clusters [i]
			mean_probsuscept = sum_typeAj / sum_sj;

			phi = num_typeAi*mean_probsuscept;
			sum_phi = sum_phi + phi;
			skip_r = 0;
		}

		//Skip if no phi values for that distance of r
		if (skip_r == 1) {
			rc[k] = NA_REAL;
		} else {
			rc[k] = (double)sum_phi/sum_typeAi;
		}

	}
}





/*****************************************************************/
/*pi function optimized for iterating through a list of two types*/
/*finding ones that fulfill the distance requirement.            */
/*                                                               */
/*@param type_1  the first type column                           */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param alpha   the weighted number of individuals/cluster      */
/*@param delta   the weighted number of non-vacc indivs/cluster  */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_pi_a_clustsurvey_wts (double *p,
                    		   double *x,
                    		   double *y,
                    		   int *s,
                    		   double *weight,
                    		   double *delta,
                    		   double *alpha,
                    		   int *len,
                    		   double *r_low,
                    		   double *r,
                    		   int *len_r,
                    		   int *inds,
                    		   double *rc,
                    		   int *remove_self) {
    

    int i,j,k;
    int skip_r, skip_phi;  /*counters for those filling conditions*/
    double sum_numerator_i, sum_denom_i, sum_delta_j, sum_alpha_j;
    double dist;

  //repeat the calculation for all r
	for (k=0;k<*len_r;k++) {
	    sum_numerator_i = 0;
	    sum_denom_i = 0;
		skip_r = 1;

		for (i=0;i<*len;i++) {
		  // Skip cluster if no individuals with the characteristic
		  if (p[i] == 0) continue;
			
			sum_delta_j = 0;
			sum_alpha_j = 0;
			skip_phi = 1;

			// Cycle through paired clusters with cluster j
			//  - need to add s and r to each cluster within the correct distance to the sums
			//  - need to sum the sum_typeAj; 
			for (j=0;j<*len;j++) {
        // Check distance of clusters
				dist = sqrt(pow(x[i]-x[j],2) + pow(y[i]-y[j],2));
				  
				if ((dist<=r[k]) & (dist>=r_low[k])) {
				  
				  if ((inds[j] == inds[i]) & (*remove_self == 1)) {
				      // Weighted number of susceptibles
				      sum_delta_j = sum_delta_j + delta[j] - weight[j];
				      // Weighted number of individuals
				      sum_alpha_j = sum_alpha_j + alpha[j] - weight[j];
				  } 
				  else if ((inds[j] != inds[i]) | (*remove_self == 0)) {
				      sum_delta_j = sum_delta_j + delta[j];
	    			  sum_alpha_j = sum_alpha_j + alpha[j];
				  }
				  
				  if (sum_alpha_j > 0) {
				      skip_phi = 0;
				      skip_r = 0;
				  }
				}
			}
			
			if (skip_phi == 1) continue;
			
			sum_numerator_i = sum_numerator_i + (delta[i] * sum_delta_j);
			sum_denom_i = sum_denom_i + (delta[i] * sum_alpha_j);
		}

		//Skip if no phi values for that distance of r
		if (skip_r == 1){
		  rc[k] = NA_REAL;
		}else {
		  rc[k] = (double) sum_numerator_i / sum_denom_i;
		}
	}
}











/*****************************************************************/
/*pi function optimized for iterating through a list of two types*/
/*finding ones that fulfill the distance requirement.            */
/*                                                               */
/*@param type_1  the first type column                           */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param alpha   the weighted number of individuals/cluster      */
/*@param delta   the weighted number of non-vacc indivs/cluster  */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_pi_b_clustsurvey_wts (double *p,
                             double *x,
                             double *y,
                             int *s,
                             double *weight,
                             double *delta,
                             double *alpha,
                             int *len,
                             double *r_low,
                             double *r,
                             int *len_r,
                             int *inds,
                             double *rc,
                             int *remove_self) {
    
    int i,j,k;
    int skip_r, skip_phi;  /*counters for those filling conditions*/
    double sum_numerator_i, sum_denom_i, sum_delta_j, sum_alpha_j;
    double dist;
    
    //repeat the calculation for all r
    for (k=0;k<*len_r;k++) {
        sum_numerator_i = 0;
        sum_denom_i = 0;
        skip_r = 1;
        
        for (i=0;i<*len;i++) {

            sum_delta_j = 0;
            sum_alpha_j = 0;
            skip_phi = 1;
            
            // Cycle through paired clusters with cluster j
            //  - need to add s and r to each cluster within the correct distance to the sums
            //  - need to sum the sum_typeAj; 
            for (j=0;j<*len;j++) {
                // Check distance of clusters
                dist = sqrt(pow(x[i]-x[j],2) + pow(y[i]-y[j],2));
                
                if ((dist<=r[k]) & (dist>=r_low[k])) {
                    
                    // adjust for self
                    if ((inds[j] == inds[i]) & (*remove_self == 1)) {
                        // Weighted number of susceptibles
                        sum_delta_j = sum_delta_j + delta[j] - weight[j]; 
                        // Weighted number of individuals
                        sum_alpha_j = sum_alpha_j + alpha[j] - weight[j];
                    } 
                    else if ((inds[j] != inds[i]) | (*remove_self == 0)) {
                        sum_delta_j = sum_delta_j + delta[j];
                        sum_alpha_j = sum_alpha_j + alpha[j];
                    }
                    
                    if (sum_alpha_j > 0) {
                        skip_phi = 0;
                        skip_r = 0;
                    }
                }
            }
            
            if (skip_phi == 1) continue;
            
            sum_numerator_i = sum_numerator_i + (alpha[i] * sum_delta_j);
            sum_denom_i = sum_denom_i + (alpha[i] * sum_alpha_j);
        }
        
        //Skip if no phi values for that distance of r
        if (skip_r == 1){
            rc[k] = NA_REAL;
        }else {
            rc[k] = (double) sum_numerator_i / sum_denom_i;
        }
    }
}

















/*****************************************************************/
/*pi function optimized for iterating through a list of two types*/
/*finding ones that fulfill the distance requirement.            */
/*                                                               */
/*@param type_1  the first type column                           */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param alpha   the weighted number of individuals/cluster      */
/*@param delta   the weighted number of non-vacc indivs/cluster  */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_pi_a_typed_clustsurvey_wts (double *p,
                                     int *type_a,
                                     int *type_b,
                                     double *x,
                                     double *y,
                                     int *s,
                                     double *weight,
                                     double *delta,
                                     double *alpha,
                                     int *len,
                                     int *typeA,
                                     int *typeB,
                                     double *r_low,
                                     double *r,
                                     int *len_r,
                                     int *inds,
                                     double *rc,
                                     int *remove_self) {
        int i,j,k;
    int skip_r, skip_phi;  /*counters for those filling conditions*/
    double sum_numerator_i, sum_denom_i, sum_delta_j, sum_alpha_j;
    double dist;
    
    //repeat the calculation for all r
    for (k=0;k<*len_r;k++) {
        sum_numerator_i = 0;
        sum_denom_i = 0;
        skip_r = 1;
        
        ///Loop through all Type A points
        for (i=0;i<*len;i++) {
            
            // Skip any that are not Type A
            if (type_a[i] != *typeA) continue;
            
            // Skip cluster if no individuals with the characteristic
            if (p[i] == 0) continue;
            
            sum_delta_j = 0;
            sum_alpha_j = 0;
            skip_phi = 1;
            
            // Cycle through paired clusters with cluster j
            //  - need to add s and r to each cluster within the correct distance to the sums
            //  - need to sum the sum_typeAj; 
            for (j=0;j<*len;j++) {
                
                // Skip any pairs that are not Type B
                if (type_b[j] != *typeB) continue;
                
                // Check distance of clusters
                dist = sqrt(pow(x[i]-x[j],2) + pow(y[i]-y[j],2));
                
                if ((dist<=r[k]) & (dist>=r_low[k])) {
                    
                    if ((inds[j] == inds[i]) & (*remove_self == 1)) {
                        sum_delta_j = sum_delta_j + delta[j] - weight[j];
                        sum_alpha_j = sum_alpha_j + alpha[j] - weight[j];
                    } 
                    else if ((inds[j] != inds[i]) | (*remove_self == 0)) {
                        sum_delta_j = sum_delta_j + delta[j];
                        sum_alpha_j = sum_alpha_j + alpha[j];
                    }
                    
                    if (sum_alpha_j > 0) {
                        skip_phi = 0;
                        skip_r = 0;
                    }
                }
            }
            
            if (skip_phi == 1) continue;
            
            sum_numerator_i = sum_numerator_i + (delta[i] * sum_delta_j);
            sum_denom_i = sum_denom_i + (delta[i] * sum_alpha_j);
        }
        
        //Skip if no phi values for that distance of r
        if (skip_r == 1){
            rc[k] = NA_REAL;
        }else {
            rc[k] = (double) sum_numerator_i / sum_denom_i;
        }
    }
}





/*****************************************************************/
/*pi function optimized for iterating through a list of two types*/
/*finding ones that fulfill the distance requirement.            */
/*                                                               */
/*@param type_1  the first type column                           */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param alpha   the weighted number of individuals/cluster      */
/*@param delta   the weighted number of non-vacc indivs/cluster  */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_pi_b_typed_clustsurvey_wts (double *p,
                                   int *type_a,
                                   int *type_b,
                                   double *x,
                                   double *y,
                                   int *s,
                                   double *weight,
                                   double *delta,
                                   double *alpha,
                                   int *len,
                                   int *typeA,
                                   int *typeB,
                                   double *r_low,
                                   double *r,
                                   int *len_r,
                                   int *inds,
                                   double *rc,
                                   int *remove_self) {
    
    int i,j,k;
    int skip_r, skip_phi;  /*counters for those filling conditions*/
    double sum_numerator_i, sum_denom_i, sum_delta_j, sum_alpha_j;
    double dist;
    
    //repeat the calculation for all r
    for (k=0;k<*len_r;k++) {
        sum_numerator_i = 0;
        sum_denom_i = 0;
        skip_r = 1;
        
        ///Loop through all Type A points
        for (i=0;i<*len;i++) {
            
            // Skip any that are not Type A
            if (type_a[i] != *typeA) continue;
            
            sum_delta_j = 0;
            sum_alpha_j = 0;
            skip_phi = 1;
            
            // Cycle through paired clusters with cluster j
            //  - need to add s and r to each cluster within the correct distance to the sums
            //  - need to sum the sum_typeAj; 
            for (j=0;j<*len;j++) {
                
                // Skip any pairs that are not Type B
                if (type_b[j] != *typeB) continue;
                
                // Check distance of clusters
                dist = sqrt(pow(x[i]-x[j],2) + pow(y[i]-y[j],2));
                
                if ((dist<=r[k]) & (dist>=r_low[k])) {
                    
                    if ((inds[j] == inds[i]) & (*remove_self == 1)) {
                        sum_delta_j = sum_delta_j + delta[j] - weight[j];
                        sum_alpha_j = sum_alpha_j + alpha[j] - weight[j];
                    } 
                    else if ((inds[j] != inds[i]) | (*remove_self == 0)) {
                        sum_delta_j = sum_delta_j + delta[j];
                        sum_alpha_j = sum_alpha_j + alpha[j];
                    }
                    
                    if (sum_alpha_j > 0) {
                        skip_phi = 0;
                        skip_r = 0;
                    }
                }
            }
            
            if (skip_phi == 1) continue;
            
            sum_numerator_i = sum_numerator_i + (alpha[i] * sum_delta_j);
            sum_denom_i = sum_denom_i + (alpha[i] * sum_alpha_j);
        }
        
        //Skip if no phi values for that distance of r
        if (skip_r == 1){
            rc[k] = NA_REAL;
        }else {
            rc[k] = (double) sum_numerator_i / sum_denom_i;
        }
    }
}









/*****************************************************************/
/*pi function optimized for iterating through a list of two types*/
/*finding ones that fulfill the distance requirement.            */
/*                                                               */
/*@param type_1  the first type column                           */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param alpha   the weighted number of individuals/cluster      */
/*@param delta   the weighted number of non-vacc indivs/cluster  */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_pi_a_clustsurvey_hh_wts (double *p,
                             double *x,
                             double *y,
                             int *s,
                             double *weight,
                             double *delta,
                             double *alpha,
                             double *rho1,
                             double *rho2,
                             int *len,
                             double *r_low,
                             double *r,
                             int *len_r,
                             int *inds,
                             double *rc,
                             int *remove_self) {
    
    int i,j,k;
    int skip_r, skip_phi;  /*counters for those filling conditions*/
    double sum_numerator_i, sum_denom_i, sum_delta_j, sum_alpha_j;
    double dist;
    
    //repeat the calculation for all r
    for (k=0;k<*len_r;k++) {
        sum_numerator_i = 0;
        sum_denom_i = 0;
        skip_r = 1;
        
        for (i=0;i<*len;i++) {
            // Skip cluster if no individuals with the characteristic
            if (p[i] == 0) continue;
            
            sum_delta_j = 0;
            sum_alpha_j = 0;
            skip_phi = 1;
            
            // Cycle through paired clusters with cluster j
            //  - need to add s and r to each cluster within the correct distance to the sums
            //  - need to sum the sum_typeAj; 
            for (j=0;j<*len;j++) {
                // Check distance of clusters
                dist = sqrt(pow(x[i]-x[j],2) + pow(y[i]-y[j],2));
                
                if ((dist<=r[k]) & (dist>=r_low[k])) {
                    
                    // If same cluster remove household contact contribution
                    if (inds[j] == inds[i]) {
                        if (*remove_self == 0) {
                            sum_delta_j = sum_delta_j + (delta[j])*(1-rho1[j]);
                            sum_alpha_j = sum_alpha_j + (alpha[j])*(1-rho2[j]);
                        } 
                        if (*remove_self == 1) {
                            sum_delta_j = sum_delta_j + (delta[j]-weight[j])*(1-rho1[j]);
                            sum_alpha_j = sum_alpha_j + (alpha[j]-weight[j])*(1-rho2[j]);
                        }
                    }
                    else if (inds[j] != inds[i]) {
                        sum_delta_j = sum_delta_j + delta[j];
                        sum_alpha_j = sum_alpha_j + alpha[j];
                    }
                    
                    if (sum_alpha_j>0) {
                        skip_phi = 0;
                        skip_r = 0;
                    }
                }
            }
                
            if (skip_phi == 1) continue;
                
            sum_numerator_i = sum_numerator_i + (delta[i] * sum_delta_j);
            sum_denom_i = sum_denom_i + (delta[i] * sum_alpha_j);
        }
            
            //Skip if no phi values for that distance of r
        if (skip_r == 1){
            rc[k] = NA_REAL;
        }else {
            rc[k] = (double) sum_numerator_i / sum_denom_i;
        }
    }
}









/*****************************************************************/
/*pi function optimized for iterating through a list of two types*/
/*finding ones that fulfill the distance requirement.            */
/*                                                               */
/*@param type_1  the first type column                           */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param alpha   the weighted number of individuals/cluster      */
/*@param delta   the weighted number of non-vacc indivs/cluster  */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_pi_b_clustsurvey_hh_wts (double *p,
                                double *x,
                                double *y,
                                int *s,
                                double *weight,
                                double *delta,
                                double *alpha,
                                double *rho1,
                                double *rho2,
                                int *len,
                                double *r_low,
                                double *r,
                                int *len_r,
                                int *inds,
                                double *rc,
                                int *remove_self) {
    
    int i,j,k;
    int skip_r, skip_phi;  /*counters for those filling conditions*/
    double sum_numerator_i, sum_denom_i, sum_delta_j, sum_alpha_j;
    double dist;
    
    //repeat the calculation for all r
    for (k=0;k<*len_r;k++) {
        sum_numerator_i = 0;
        sum_denom_i = 0;
        skip_r = 1;
        
        for (i=0;i<*len;i++) {

            sum_delta_j = 0;
            sum_alpha_j = 0;
            skip_phi = 1;
            
            // Cycle through paired clusters with cluster j
            //  - need to add s and r to each cluster within the correct distance to the sums
            //  - need to sum the sum_typeAj; 
            for (j=0;j<*len;j++) {
                // Check distance of clusters
                dist = sqrt(pow(x[i]-x[j],2) + pow(y[i]-y[j],2));
                
                if ((dist<=r[k]) & (dist>=r_low[k])) {
                    
                    // If same cluster remove household contact contribution
                    if (inds[j] == inds[i]) {
                        if (*remove_self == 0) {
                            sum_delta_j = sum_delta_j + (delta[j])*(1-rho1[j]);
                            sum_alpha_j = sum_alpha_j + (alpha[j])*(1-rho2[j]);
                        } 
                        if (*remove_self == 1) {
                            sum_delta_j = sum_delta_j + (delta[j]-weight[j])*(1-rho1[j]);
                            sum_alpha_j = sum_alpha_j + (alpha[j]-weight[j])*(1-rho2[j]);
                        }
                    }
                    else if (inds[j] != inds[i]) {
                        sum_delta_j = sum_delta_j + delta[j];
                        sum_alpha_j = sum_alpha_j + alpha[j];
                    }
                    
                    if (sum_alpha_j>0) {
                        skip_phi = 0;
                        skip_r = 0;
                    }
                }
            }
            
            if (skip_phi == 1) continue;
            
            sum_numerator_i = sum_numerator_i + (alpha[i] * sum_delta_j);
            sum_denom_i = sum_denom_i + (alpha[i] * sum_alpha_j);
        }
        
        //Skip if no phi values for that distance of r
        if (skip_r == 1){
            rc[k] = NA_REAL;
        }else {
            rc[k] = (double) sum_numerator_i / sum_denom_i;
        }
    }
}











/*****************************************************************/
/*pi function optimized for iterating through a list of two types*/
/*finding ones that fulfill the distance requirement.            */
/*                                                               */
/*@param type_1  the first type column                           */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param alpha   the weighted number of individuals/cluster      */
/*@param delta   the weighted number of non-vacc indivs/cluster  */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_pi_a_typed_clustsurvey_hh_wts (double *p,
                                        int *type_a,
                                        int *type_b,
                                        double *x,
                                        double *y,
                                        int *s,
                                        double *weight,
                                        double *delta,
                                        double *alpha,
                                        double *rho1,
                                        double *rho2,
                                        int *len,
                                        int *typeA,
                                        int *typeB,
                                        double *r_low,
                                        double *r,
                                        int *len_r,
                                        int *inds,
                                        double *rc,
                                        int *remove_self) {
  
    int i,j,k;
    int skip_r, skip_phi;  /*counters for those filling conditions*/
    double sum_numerator_i, sum_denom_i, sum_delta_j, sum_alpha_j;
    double dist;
    
    //repeat the calculation for all r
    for (k=0;k<*len_r;k++) {
        sum_numerator_i = 0;
        sum_denom_i = 0;
        skip_r = 1;
        
        for (i=0;i<*len;i++) {
            
            // Skip any that are not Type A
            if (type_a[i] != *typeA) continue;
            
            // Skip cluster if no individuals with the characteristic
            if (p[i] == 0) continue;
            
            sum_delta_j = 0;
            sum_alpha_j = 0;
            skip_phi = 1;
            
            // Cycle through paired clusters with cluster j
            //  - need to add s and r to each cluster within the correct distance to the sums
            //  - need to sum the sum_typeAj; 
            for (j=0;j<*len;j++) {
                
                // Skip any pairs that are not Type B
                if (type_b[j] != *typeB) continue;
                
                // Check distance of clusters
                dist = sqrt(pow(x[i]-x[j],2) + pow(y[i]-y[j],2));
                
                if ((dist<=r[k]) & (dist>=r_low[k])) {
                    
                    // If same cluster remove household contact contribution
                    if (inds[j] == inds[i]) {
                        if (*remove_self == 0) {
                            sum_delta_j = sum_delta_j + (delta[j])*(1-rho1[j]);
                            sum_alpha_j = sum_alpha_j + (alpha[j])*(1-rho2[j]);
                        } 
                        if (*remove_self == 1) {
                            sum_delta_j = sum_delta_j + (delta[j]-weight[j])*(1-rho1[j]);
                            sum_alpha_j = sum_alpha_j + (alpha[j]-weight[j])*(1-rho2[j]);
                        }
                    }
                    else if (inds[j] != inds[i]) {
                        sum_delta_j = sum_delta_j + delta[j];
                        sum_alpha_j = sum_alpha_j + alpha[j];
                    }
                    
                    if (sum_alpha_j>0) {
                        skip_phi = 0;
                        skip_r = 0;
                    }
                }
            }
            
            if (skip_phi == 1) continue;
            
            sum_numerator_i = sum_numerator_i + (delta[i] * sum_delta_j);
            sum_denom_i = sum_denom_i + (delta[i] * sum_alpha_j);
        }
        
        //Skip if no phi values for that distance of r
        if (skip_r == 1){
            rc[k] = NA_REAL;
        }else {
            rc[k] = (double) sum_numerator_i / sum_denom_i;
        }
    }
}









/*****************************************************************/
/*pi function optimized for iterating through a list of two types*/
/*finding ones that fulfill the distance requirement.            */
/*                                                               */
/*@param type_1  the first type column                           */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param alpha   the weighted number of individuals/cluster      */
/*@param delta   the weighted number of non-vacc indivs/cluster  */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_pi_b_typed_clustsurvey_hh_wts (double *p,
                                        int *type_a,
                                        int *type_b,
                                        double *x,
                                        double *y,
                                        int *s,
                                        double *weight,
                                        double *delta,
                                        double *alpha,
                                        double *rho1,
                                        double *rho2,
                                        int *len,
                                        int *typeA,
                                        int *typeB,
                                        double *r_low,
                                        double *r,
                                        int *len_r,
                                        int *inds,
                                        double *rc,
                                        int *remove_self) {
        
    int i,j,k;
    int skip_r, skip_phi;  /*counters for those filling conditions*/
    double sum_numerator_i, sum_denom_i, sum_delta_j, sum_alpha_j;
    double dist;
    
    //repeat the calculation for all r
    for (k=0;k<*len_r;k++) {
        sum_numerator_i = 0;
        sum_denom_i = 0;
        skip_r = 1;
        
        for (i=0;i<*len;i++) {
            
            // Skip any that are not Type A
            if (type_a[i] != *typeA) continue;
            
            sum_delta_j = 0;
            sum_alpha_j = 0;
            skip_phi = 1;
            
            // Cycle through paired clusters with cluster j
            //  - need to add s and r to each cluster within the correct distance to the sums
            //  - need to sum the sum_typeAj; 
            for (j=0;j<*len;j++) {
                
                // Skip any pairs that are not Type B
                if (type_b[j] != *typeB) continue;
                
                // Check distance of clusters
                dist = sqrt(pow(x[i]-x[j],2) + pow(y[i]-y[j],2));
                
                if ((dist<=r[k]) & (dist>=r_low[k])) {
                    
                    // If same cluster remove household contact contribution
                    if (inds[j] == inds[i]) {
                        if (*remove_self == 0) {
                            sum_delta_j = sum_delta_j + (delta[j])*(1-rho1[j]);
                            sum_alpha_j = sum_alpha_j + (alpha[j])*(1-rho2[j]);
                        } 
                        if (*remove_self == 1) {
                            sum_delta_j = sum_delta_j + (delta[j]-weight[j])*(1-rho1[j]);
                            sum_alpha_j = sum_alpha_j + (alpha[j]-weight[j])*(1-rho2[j]);
                        }
                    }
                    else if (inds[j] != inds[i]) {
                        sum_delta_j = sum_delta_j + delta[j];
                        sum_alpha_j = sum_alpha_j + alpha[j];
                    }
                    
                    if (sum_alpha_j>0) {
                        skip_phi = 0;
                        skip_r = 0;
                    }
                }
            }
            
            if (skip_phi == 1) continue;
            
            sum_numerator_i = sum_numerator_i + (alpha[i] * sum_delta_j);
            sum_denom_i = sum_denom_i + (alpha[i] * sum_alpha_j);
        }
        
        //Skip if no phi values for that distance of r
        if (skip_r == 1){
            rc[k] = NA_REAL;
        }else {
            rc[k] = (double) sum_numerator_i / sum_denom_i;
        }
    }
}







/*****************************************************************/
/* tau function for typed data                                   */
/* interates through a list of two types finding ones            */
/*that fulfill the distance requirement.                         */
/* this version takes in coordinates instead of distances        */
/*                                                               */
/*@param type_1  the  first type column                          */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_tau_typed_wts (int *type,
		    double *x,
		    double *y,
			double *weight,
		    int *len,
		    int *typeA,
		    int *typeB,
		    double *r_low,
		    double *r,
		    int *len_r,
		    int *inds,
		    double *rc) {

    int i = 0;
    double divisor;
    double tmp_r_low = 0;
    double tmp_r = DBL_MAX;
    int tmp_len_r = 1;

    /*get the divisor in the pi function*/
    get_pi_typed_wts(type,x,y,weight,len,typeA,typeB,&tmp_r_low,&tmp_r,
	   &tmp_len_r,inds,&divisor);

    /*get the main pi function*/
    get_pi_typed_wts(type,x,y,weight,len,typeA,typeB,r_low,r,len_r,inds,rc);

    for (i = 0; i < *len_r; i++) {
      rc[i] = rc[i]/divisor;
    }
}



/*****************************************************************/
/* tau function for typed data                                   */
/* interates through a list of two types finding ones            */
/*that fulfill the distance requirement.                         */
/* this version takes in coordinates instead of distances        */
/*                                                               */
/*@param type_1  the  first type column                          */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_tau_clustsurvey (double *p,
		    double *x,
		    double *y,
			int *s,
		    int *len,
		    double *r_low,
		    double *r,
		    int *len_r,
		    int *inds,
		    double *rc) {

    int i = 0;
    double divisor;
    double tmp_r_low = 0;
    double tmp_r = DBL_MAX;
    int tmp_len_r = 1;

    /*get the divisor in the pi function*/
    get_pi_clustsurvey(p,x,y,s,len,&tmp_r_low,&tmp_r,&tmp_len_r,inds,&divisor);

    /*get the main pi function*/
    get_pi_clustsurvey(p,x,y,s,len,r_low,r,len_r,inds,rc);

    for (i = 0; i < *len_r; i++) {
      rc[i] = rc[i]/divisor;
    }
}



/*****************************************************************/
/* tau function for typed data                                   */
/* interates through a list of two types finding ones            */
/*that fulfill the distance requirement.                         */
/* this version takes in coordinates instead of distances        */
/*                                                               */
/*@param type_1  the  first type column                          */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_tau_clustsurvey_wts (double *p,
                      		double *x,
                      		double *y,
                      		int *s,
                      	    double *weight,
                      		double *delta,
                      		double *alpha,
                      		int *len,
                      		double *r_low,
                      		double *r,
                      		int *len_r,
                      		int *inds,
                      		double *rc,
                      		double *numerator,
                      		double *denominator,
                      		int *remove_self) {

    int i = 0;
    //double numerator[*len_r];
    //double denominator[*len_r];

    /*get the numerator pi function*/
    get_pi_a_clustsurvey_wts(p,x,y,s,weight,delta,alpha,len,r_low,r,len_r,inds,numerator,remove_self);
    
    /*get the denominator pi function*/
    get_pi_b_clustsurvey_wts(p,x,y,s,weight,delta,alpha,len,r_low,r,len_r,inds,denominator,remove_self);
    
    for (i = 0; i < *len_r; i++) {
      rc[i] = numerator[i]/denominator[i];
    }
}









/*****************************************************************/
/* tau function for typed data                                   */
/* interates through a list of two types finding ones            */
/*that fulfill the distance requirement.                         */
/* this version takes in coordinates instead of distances        */
/*                                                               */
/*@param type_1  the  first type column                          */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_tau_typed_clustsurvey_wts (double *p,
                                    int *type_a,
                                    int *type_b,
                                      double *x,
                                      double *y,
                                      int *s,
                                      double *weight,
                                      double *delta,
                                      double *alpha,
                                      int *len,
                                      int *typeA,
                                      int *typeB,
                                      double *r_low,
                                      double *r,
                                      int *len_r,
                                      int *inds,
                                      double *rc,
                                      double *numerator,
                                      double *denominator,
                                      int *remove_self) {
            
    int i = 0;

    /*get the numerator pi function*/
    get_pi_a_typed_clustsurvey_wts(p,type_a,type_b,x,y,s,weight,delta,alpha,len,typeA,typeB,r_low,r,len_r,inds,numerator,remove_self);

    /*get the denominator pi function*/
    get_pi_b_typed_clustsurvey_wts(p,type_a,type_b,x,y,s,weight,delta,alpha,len,typeA,typeB,r_low,r,len_r,inds,denominator,remove_self);

    for (i = 0; i < *len_r; i++) {
        rc[i] = numerator[i]/denominator[i];
    }
}


/*****************************************************************/
/* tau function for typed data                                   */
/* interates through a list of two types finding ones            */
/*that fulfill the distance requirement.                         */
/* this version takes in coordinates instead of distances        */
/*                                                               */
/*@param type_1  the  first type column                          */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_tau_typed_clustsurvey_hh_wts (double *p,
                                    int *type_a,
                                    int *type_b,
                                    double *x,
                                    double *y,
                                    int *s,
                                    double *weight,
                                    double *delta,
                                    double *alpha,
                                    double *rho1,
                                    double *rho2,
                                    int *len,
                                    int *typeA,
                                    int *typeB,
                                    double *r_low,
                                    double *r,
                                    int *len_r,
                                    int *inds,
                                    double *rc,
                                    double *numerator,
                                    double *denominator,
                                    int *remove_self) {
    
    int i = 0;
    
    /*get the numerator pi function*/
    get_pi_a_typed_clustsurvey_hh_wts(p,type_a,type_b,x,y,s,weight,delta,alpha,rho1,rho2,len,typeA,typeB,r_low,r,len_r,inds,numerator,remove_self);
    
    /*get the denominator pi function*/
    get_pi_b_typed_clustsurvey_hh_wts(p,type_a,type_b,x,y,s,weight,delta,alpha,rho1,rho2,len,typeA,typeB,r_low,r,len_r,inds,denominator,remove_self);
    
    for (i = 0; i < *len_r; i++) {
        rc[i] = numerator[i]/denominator[i];
    }
}







/*****************************************************************/
/* tau function for cluster survey data                                   */
/* interates through a list of two types finding ones            */
/*that fulfill the distance requirement.                         */
/* this version takes in coordinates instead of distances        */
/*                                                               */
/*@param type_1  the  first type column                          */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_tau_clustsurvey_hh_wts (double *p,
                              double *x,
                              double *y,
                              int *s,
                              double *weight,
                              double *delta,
                              double *alpha,
                              double *rho1,
                              double *rho2,
                              int *len,
                              double *r_low,
                              double *r,
                              int *len_r,
                              int *inds,
                              double *rc,
                              double *numerator,
                              double *denominator,
                              int *remove_self) {
    int i = 0;

    /*get the divisor in the pi function*/
    get_pi_a_clustsurvey_hh_wts(p,x,y,s,weight,delta,alpha,rho1,rho2,len,r_low,r,len_r,inds,numerator,remove_self);
    
    /*get the main pi function*/
    get_pi_b_clustsurvey_hh_wts(p,x,y,s,weight,delta,alpha,rho1,rho2,len,r_low,r,len_r,inds,denominator,remove_self);
    
    for (i = 0; i < *len_r; i++) {
        rc[i] = numerator[i]/denominator[i];
    }
}